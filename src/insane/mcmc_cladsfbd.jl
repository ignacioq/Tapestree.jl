#=

clads birth-death MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 28 07 2025
=#




"""
    insane_cladsbd(tree    ::sT_label;
                   λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                   α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                   σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                   μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                   niter   ::Int64                 = 1_000,
                   nthin   ::Int64                 = 10,
                   nburn   ::Int64                 = 200,
                   nflush  ::Int64                 = nthin,
                   ofile   ::String                = string(homedir(), "/cladsbd"),
                   λi      ::Float64               = NaN,
                   αi      ::Float64               = 0.0,
                   σλi     ::Float64               = 0.1,
                   μi      ::Float64               = NaN,
                   ϵi      ::Float64               = 0.2,
                   pupdp   ::NTuple{6,Float64}     = (1e-3, 1e-3, 1e-3, 1e-4, 0.1, 0.2),
                   prints  ::Int64                 = 5,
                   stn     ::Float64               = 0.5,
                   survival::Bool                  = true,
                   mxthf   ::Float64               = 0.1,
                   tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for clads birth-death.
"""
function insane_cladsbd(tree    ::sT_label;
                        λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                        μ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                        α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                        σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                        σμ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                        niter   ::Int64                 = 1_000,
                        nthin   ::Int64                 = 10,
                        nburn   ::Int64                 = 200,
                        nflush  ::Int64                 = nthin,
                        ofile   ::String                = string(homedir(), "/cladsbd"),
                        λi      ::Float64               = NaN,
                        μi      ::Float64               = NaN,
                        ϵi      ::Float64               = 0.2,
                        αi      ::Float64               = 0.0,
                        σλi     ::Float64               = 0.1,
                        σμi     ::Float64               = 0.1,
                        pupdp   ::NTuple{6,Float64}     = (1e-3, 1e-3, 1e-3, 1e-4, 0.2, 0.2),
                        prints  ::Int64                 = 5,
                        stnλ    ::Float64               = 0.5,
                        stnμ    ::Float64               = 0.5,
                        survival::Bool                  = true,
                        mxthf   ::Float64               = 0.1,
                        tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n  = ntips(tree)
  th = treeheight(tree)

  # turn to logarithmic terms
  λ0_prior = (log(λ0_prior[1]), 2.0*log(λ0_prior[2]))
  μ0_prior = (log(μ0_prior[1]), 2.0*log(μ0_prior[2]))

  surv = 0   # condition on survival of 0, 1, or 2 starting lineages
  rmλ  = 0.0 # condition on first speciation event
  if iszero(e(tree)) 
    rmλ  += 1.0
    surv += survival ? 2 : 0
  else
    surv += survival ? 1 : 0
  end

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ, th * mxthf)

  # starting parameters
  λc, μc = λi, μi
  if isnan(λi) || isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  end

  # make a decoupled tree
  Ξ = make_Ξ(idf, λc, μc, cTbd)

  # survival
  mc = m_surv_cladsbd(th, log(λc), log(μc), αi, σλi, σμi, 1_000, surv)

  # parameter updates (1: α, 2: σλ, 3: σμ, 4: scale, 5: internal, 6: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death clads"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, σμc, mc, ns, ne, ddλ, ssλ, ssμ, stnλ, stnμ =
    mcmc_burn_cladsbd(Ξ, idf, λ0_prior, μ0_prior, α_prior, σλ_prior, σμ_prior, 
      nburn, αi, σλi, σμi, mc, th, rmλ, surv, stnλ, stnμ, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cladsbd(Ξ, idf, llc, prc, αc, σλc, σμc, mc, th, rmλ, surv, ns, ne, 
      ddλ, ssλ, ssμ, stnλ, stnμ, λ0_prior, μ0_prior, α_prior, 
      σλ_prior, σμ_prior, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_cladsbd(Ξ       ::Vector{cTbd},
                      idf     ::Vector{iBffs},
                      λ0_prior::NTuple{2,Float64},
                      μ0_prior::NTuple{2,Float64},
                      α_prior ::NTuple{2,Float64},
                      σλ_prior::NTuple{2,Float64},
                      σμ_prior::NTuple{2,Float64},
                      nburn   ::Int64,
                      αc      ::Float64,
                      σλc     ::Float64,
                      σμc     ::Float64,
                      mc      ::Float64,
                      th      ::Float64,
                      rmλ     ::Float64,
                      surv    ::Int64,
                      stnλ    ::Float64,
                      stnμ    ::Float64,
                      pup     ::Array{Int64,1},
                      prints  ::Int64)

MCMC burn-in chain for `cladsbd`.
"""
function mcmc_burn_cladsbd(Ξ       ::Vector{cTbd},
                           idf     ::Vector{iBffs},
                           λ0_prior::NTuple{2,Float64},
                           μ0_prior::NTuple{2,Float64},
                           α_prior ::NTuple{2,Float64},
                           σλ_prior::NTuple{2,Float64},
                           σμ_prior::NTuple{2,Float64},
                           nburn   ::Int64,
                           αc      ::Float64,
                           σλc     ::Float64,
                           σμc     ::Float64,
                           mc      ::Float64,
                           th      ::Float64,
                           rmλ     ::Float64,
                           surv    ::Int64,
                           stnλ    ::Float64,
                           stnμ    ::Float64,
                           pup     ::Array{Int64,1},
                           prints  ::Int64)

  # starting likelihood and prior
  lλ0 = lλ(Ξ[1])
  llc = llik_clads(Ξ, idf, αc, σλc, σμc) - rmλ*lλ0 + log(mc) + prob_ρ(idf)
  prc = logdnorm(lλ0,       λ0_prior[1], λ0_prior[2])   +
        logdnorm(lμ(Ξ[1]),  μ0_prior[1], μ0_prior[2])   +
        logdnorm(αc,         α_prior[1],  α_prior[2]^2) +
        logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])   +
        logdinvgamma(σμc^2, σμ_prior[1], σμ_prior[2])

  el  = lastindex(idf)                          # number of branches
  ns  = sum(x -> Float64(d2(x) > 0), idf) - rmλ # number of speciation events in likelihood
  ne  = 0.0                                     # number of extinction events
  λfs = Float64[]
  μfs = Float64[]

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, ssμ = _dd_ss(Ξ, idf, αc)

  # for scale tuning
  ltn = zero(Int64)
  lup = lacλ = lacμ = zero(Float64)

  pbar = Progress(nburn, dt = prints, desc = "burning mcmc...", barlen = 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update drift
      if pupi === 1

        llc, prc, αc, mc = 
          update_α!(αc, lλ(Ξ[1]), lμ(Ξ[1]), σλc, σμc, 
            2.0*(ns + rmλ), ddλ, llc, prc, mc, th, surv, α_prior)

        # update ssλ with new drift `α`
        ssλ = _ss(Ξ, idf, αc)

      # update speciation diffusion rate
      elseif pupi === 2

        llc, prc, σλc, mc = 
          update_σλ!(σλc, lλ(Ξ[1]), lμ(Ξ[1]), αc, σμc, ssλ, 
            2.0*(ns + rmλ), llc, prc, mc, th, surv, σλ_prior)

      # update extinction diffusion rate
      elseif pupi === 3

        llc, prc, σμc, mc = 
          update_σμ!(σμc, lλ(Ξ[1]), lμ(Ξ[1]), αc, σλc, ssμ, 
            2.0*(ns + rmλ), llc, prc, mc, th, surv, σμ_prior)

      # update scale
      elseif pupi === 4

        llc, prc, mc, accλ, accμ = 
          update_scale!(Ξ, idf, αc, σλc, σμc, llc, prc, ns, ne, 
            stnλ, stnμ, mc, th, surv, λ0_prior, μ0_prior)

        lacλ += accλ
        lacμ += accμ
        lup  += 1.0

      # update internal
      elseif pupi === 5

        bix = fIrand(el) + 1

        llc, prc, ddλ, ssλ, ssμ, mc =
          update_internal!(bix, Ξ, idf, αc, σλc, σμc, llc, prc, 
            ddλ, ssλ, ssμ, mc, th, λ0_prior, μ0_prior, surv)

      # forward simulation
      else

        bix = fIrand(el) + 1

        llc, ddλ, ssλ, ssμ, ns, ne =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, llc, ddλ, ssλ, ssμ, ns, 
            ne, λfs, μfs)
      end
    end

    ltn += 1
    if ltn === 100
      stnλ = min(2.0, tune(stnλ, lacλ/lup))
      stnμ = min(2.0, tune(stnμ, lacμ/lup))
      ltn = zero(Int64)
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, σμc, mc, ns, ne, ddλ, ssλ, ssμ, stnλ, stnμ
end




"""
    mcmc_cladsbd(Ξ       ::Vector{cTbd},
                 idf     ::Vector{iBffs},
                 llc     ::Float64,
                 prc     ::Float64,
                 αc      ::Float64,
                 σλc     ::Float64,
                 σμc     ::Float64,
                 mc      ::Float64,
                 th      ::Float64,
                 rmλ     ::Float64,
                 surv    ::Int64,
                 ns      ::Float64,
                 ne      ::Float64,
                 ddλ     ::Float64,
                 ssλ     ::Float64,
                 ssμ     ::Float64,
                 stnλ    ::Float64,
                 stnμ    ::Float64,
                 λ0_prior::NTuple{2,Float64},
                 μ0_prior::NTuple{2,Float64},
                 α_prior ::NTuple{2,Float64},
                 σλ_prior::NTuple{2,Float64},
                 σμ_prior::NTuple{2,Float64},
                 pup     ::Vector{Int64},
                 niter   ::Int64,
                 nthin   ::Int64,
                 nflush  ::Int64,
                 ofile   ::String,
                 prints  ::Int64)

MCMC chain for pure-birth diffusion.
"""
function mcmc_cladsbd(Ξ       ::Vector{cTbd},
                      idf     ::Vector{iBffs},
                      llc     ::Float64,
                      prc     ::Float64,
                      αc      ::Float64,
                      σλc     ::Float64,
                      σμc     ::Float64,
                      mc      ::Float64,
                      th      ::Float64,
                      rmλ     ::Float64,
                      surv    ::Int64,
                      ns      ::Float64,
                      ne      ::Float64,
                      ddλ     ::Float64,
                      ssλ     ::Float64,
                      ssμ     ::Float64,
                      stnλ    ::Float64,
                      stnμ    ::Float64,
                      λ0_prior::NTuple{2,Float64},
                      μ0_prior::NTuple{2,Float64},
                      α_prior ::NTuple{2,Float64},
                      σλ_prior::NTuple{2,Float64},
                      σμ_prior::NTuple{2,Float64},
                      pup     ::Vector{Int64},
                      niter   ::Int64,
                      nthin   ::Int64,
                      nflush  ::Int64,
                      ofile   ::String,
                      prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin = lit = sthin = zero(Int64)

  # parameter results
  r   = Array{Float64,2}(undef, nlogs, 8)

  treev = cTbd[]           # make Ξ vector
  io    = IOBuffer()       # buffer 
  el    = lastindex(idf)   # number of branches
  λfs   = Float64[]
  μfs   = Float64[]

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\tmu_root\talpha\tsigma_lambda\tsigma_mu\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αc = αc, σλc = σλc, σμc = σμc, mc = mc, ns = ns, ne = ne, ddλ = ddλ, ssλ = ssλ, ssμ = ssμ, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update drift
            if pupi === 1

              llc, prc, αc, mc = 
                update_α!(αc, lλ(Ξ[1]), lμ(Ξ[1]), σλc, σμc, 
                  2.0*(ns + rmλ), ddλ, llc, prc, mc, th, surv, α_prior)

              # update ssλ with new drift `α`
              ssλ = _ss(Ξ, idf, αc)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, σμc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update speciation diffusion rate
            elseif pupi === 2

              llc, prc, σλc, mc = 
                update_σλ!(σλc, lλ(Ξ[1]), lμ(Ξ[1]), αc, σμc, ssλ, 
                  2.0*(ns + rmλ), llc, prc, mc, th, surv, σλ_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, σμc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update extinction diffusion rate
            elseif pupi === 3

              llc, prc, σμc, mc = 
                update_σμ!(σμc, lλ(Ξ[1]), lμ(Ξ[1]), αc, σλc, ssμ, 
                  2.0*(ns + rmλ), llc, prc, mc, th, surv, σμ_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, σμc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update scale
            elseif pupi === 4

              llc, prc, mc, accλ, accμ = 
                update_scale!(Ξ, idf, αc, σλc, σμc, llc, prc, ns, ne, 
                  stnλ, stnμ, mc, th, surv, λ0_prior, μ0_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, σμc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update internal λ
            elseif pupi === 5

              bix = fIrand(el) + 1

              llc, prc, ddλ, ssλ, ssμ, mc =
                update_internal!(bix, Ξ, idf, αc, σλc, σμc, llc, prc, 
                  ddλ, ssλ, ssμ, mc, th, λ0_prior, μ0_prior, surv)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, σμc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update by forward simulation
            else

              bix = fIrand(el) + 1

              llc, ddλ, ssλ, ssμ, ns, ne =
                update_fs!(bix, Ξ, idf, αc, σλc, σμc, llc, ddλ, ssλ, ssμ, ns, 
                  ne, λfs, μfs)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, σμc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end
            end
          end

          # log parameters
          lthin += 1
          if lthin === nthin
            lit += 1
            @inbounds begin
              r[lit,1] = Float64(it)
              r[lit,2] = llc
              r[lit,3] = prc
              r[lit,4] = exp(lλ(Ξ[1]))
              r[lit,5] = exp(lμ(Ξ[1]))
              r[lit,6] = αc
              r[lit,7] = σλc
              r[lit,8] = σμc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
                  exp(lλ(Ξ[1])),'\t', exp(lμ(Ξ[1])),'\t', αc, '\t', 
                  σλc, '\t', σμc, '\n')
            flush(of)
            ibuffer(io, couple(Ξ, idf, 1))
            write(io, '\n')
            write(tf, take!(io))
            flush(tf)
            sthin = zero(Int64)
          end

          next!(pbar)
        end
      end
    end
  end

  return r, treev
end




"""
    update_α!(αc     ::Float64,
              λ0     ::Float64,
              μ0     ::Float64,
              σλ     ::Float64,
              σμ     ::Float64,
              L      ::Float64,
              ddλ    ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              mc     ::Float64,
              th     ::Float64,
              surv   ::Int64,
              α_prior::NTuple{2,Float64})

Gibbs update for `α`.
"""
function update_α!(αc     ::Float64,
                   λ0     ::Float64,
                   μ0     ::Float64,
                   σλ     ::Float64,
                   σμ     ::Float64,
                   L      ::Float64,
                   ddλ    ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   surv   ::Int64,
                   α_prior::NTuple{2,Float64})

  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2
  αp  = rnorm((ddλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  mp  = m_surv_cladsbd(th, λ0, μ0, αp, σλ, σμ, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*ddλ*(αp - αc)/L) + llr
    prc += llrdnorm_x(αp, αc, ν, τ2)
    αc   = αp
    mc   = mp
  end

  return llc, prc, αc, mc
end




"""
    update_σλ!(σλc     ::Float64,
               λ0      ::Float64,
               μ0      ::Float64,
               α       ::Float64,
               σμ      ::Float64,
               ssλ     ::Float64,
               n       ::Float64,
               llc     ::Float64,
               prc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               surv   ::Int64,
               σλ_prior::NTuple{2,Float64})

Gibbs update for `σλ`.
"""
function update_σλ!(σλc     ::Float64,
                    λ0      ::Float64,
                    μ0      ::Float64,
                    α       ::Float64,
                    σμ      ::Float64,
                    ssλ     ::Float64,
                    n       ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    surv   ::Int64,
                    σλ_prior::NTuple{2,Float64})

  σλ_p1, σλ_p2 = σλ_prior

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)
  σλp  = sqrt(σλp2)

  mp  = m_surv_cladsbd(th, λ0, μ0, α, σλp, σμ, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += ssλ*(1.0/σλc^2 - 1.0/σλp2) - n*(log(σλp/σλc)) + llr
    prc += llrdinvgamma(σλp2, σλc^2, σλ_p1, σλ_p2)
    σλc  = σλp
    mc   = mp
  end

  return llc, prc, σλc, mc
end




"""
    update_σμ!(σμc     ::Float64,
               λ0      ::Float64,
               μ0      ::Float64,
               α       ::Float64,
               σλ      ::Float64,
               ssμ     ::Float64,
               n       ::Float64,
               llc     ::Float64,
               prc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               surv    ::Int64,
               σμ_prior::NTuple{2,Float64})

Gibbs update for `σμ`.
"""
function update_σμ!(σμc     ::Float64,
                    λ0      ::Float64,
                    μ0      ::Float64,
                    α       ::Float64,
                    σλ      ::Float64,
                    ssμ     ::Float64,
                    n       ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    surv    ::Int64,
                    σμ_prior::NTuple{2,Float64})

  σμ_p1, σμ_p2 = σμ_prior

  # Gibbs update for σ
  σμp2 = randinvgamma(σμ_p1 + 0.5 * n, σμ_p2 + ssμ)
  σμp  = sqrt(σμp2)

  mp  = m_surv_cladsbd(th, λ0, μ0, α, σλ, σμp, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += ssμ*(1.0/σμc^2 - 1.0/σμp2) - n*(log(σμp/σμc)) + llr
    prc += llrdinvgamma(σμp2, σμc^2, σμ_p1, σμ_p2)
    σμc  = σμp
    mc   = mp
  end

  return llc, prc, σμc, mc
end




"""
    update_scale!(Ξ       ::Vector{cTbd},
                  idf     ::Vector{iBffs},
                  α       ::Float64,
                  σλ      ::Float64,
                  σμ      ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ns      ::Float64,
                  ne      ::Float64,
                  stnλ    ::Float64,
                  stnμ    ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  surv    ::Int64,
                  λ0_prior::NTuple{2,Float64},
                  μ0_prior::NTuple{2,Float64})

Update scale for speciation and extinction.
"""
function update_scale!(Ξ       ::Vector{cTbd},
                       idf     ::Vector{iBffs},
                       α       ::Float64,
                       σλ      ::Float64,
                       σμ      ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ns      ::Float64,
                       ne      ::Float64,
                       stnλ    ::Float64,
                       stnμ    ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       surv    ::Int64,
                       λ0_prior::NTuple{2,Float64},
                       μ0_prior::NTuple{2,Float64})

  # extract integrated rates
  irλ, irμ = _irbd(Ξ)

  accλ = accμ = 0.0

  lλ0 = lλ(Ξ[1])[1]
  lμ0 = lμ(Ξ[1])[1]

  ## start with speciation
  # sample log(scaling factor)
  s = randn()*stnλ

  # likelihood ratio
  llr = ns * s + (1.0 - exp(s)) * irλ

  # prior ratio
  prr = llrdnorm_x(lλ0 + s, lλ0, λ0_prior[1], λ0_prior[2]) 

  lU = -randexp()

  if lU < llr + prr + log(1000.0/mc)

    # survival
    mp   = m_surv_cladsbd(th, lλ0 + s, lμ0, α, σλ, σμ, 1_000, surv)
    llr += log(mp/mc)

    if -randexp() < llr + prr
      accλ += 1.0
      llc  += llr
      prc  += prr
      mc   = mp
      scale_rate!(Ξ,   addlλ!, s)
      scale_rate!(idf, addlλ!, s)
      lλ0  += s
    end
  end

  ## continue with extinction
  # sample log(scaling factor)
  s = randn()*stnμ

  # likelihood ratio
  llr = ne * s + (1.0 - exp(s)) * irμ

  # prior ratio
  prr = llrdnorm_x(lμ0 + s, lμ0, μ0_prior[1], μ0_prior[2]) 

  lU = -randexp()

  if lU < llr + prr + log(1000.0/mc)

    # add survival ratio
    mp   = m_surv_cladsbd(th, lλ0, lμ0 + s, α, σλ, σμ, 1_000, surv)
    llr += log(mp/mc)

    if lU < llr + prr
      accμ += 1.0
      llc  += llr
      prc  += prr
      mc    = mp
      scale_rate!(Ξ,   addlμ!, s)
      scale_rate!(idf, addlμ!, s)
    end
  end

  return llc, prc, mc, accλ, accμ
end




"""
    update_internal!(bix     ::Int64,
                     Ξ       ::Vector{cTbd},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     σμ      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ddλ     ::Float64,
                     ssλ     ::Float64,
                     ssμ     ::Float64,
                     mc      ::Float64,
                     th      ::Float64,
                     λ0_prior::NTuple{2,Float64},
                     μ0_prior::NTuple{2,Float64},
                     surv    ::Int64)

Make an update for an internal branch and its descendants.
"""
function update_internal!(bix     ::Int64,
                          Ξ       ::Vector{cTbd},
                          idf     ::Vector{iBffs},
                          α       ::Float64,
                          σλ      ::Float64,
                          σμ      ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          ddλ     ::Float64,
                          ssλ     ::Float64,
                          ssμ     ::Float64,
                          mc      ::Float64,
                          th      ::Float64,
                          λ0_prior::NTuple{2,Float64},
                          μ0_prior::NTuple{2,Float64},
                          surv    ::Int64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  it   = iszero(i1) # is terminal
  i2   = d2(bi)
  ia   = pa(bi)
  root = iszero(ia)
  λa = μa = NaN  # ancestral speciation and extinction

  # if crown root
  if root && iszero(e(ξi))
    llc, prc, ddλ, ssλ, ssμ, mc =
      _crown_update!(ξi, Ξ[i1], Ξ[i2], α, σλ, σμ, llc, prc, ddλ, ssλ, ssμ, 
        mc, th, λ0_prior, μ0_prior, surv)
    λa, μa = lλ(ξi), lμ(ξi)
    setλt!(bi, λa)
    setμt!(bi, μa)
  else
    # if stem
    if root
      # if cladogenetic branch
      if i2 > 0
        ξ1, ξ2 = Ξ[i1], Ξ[i2]
        eds, λ1, λ2, μ1, μ2 = 0.0, lλ(ξ1), lλ(ξ2), lμ(ξ1), lμ(ξ2)
      # if mid branch
      else
        eds, λ1, λ2, μ1, μ2 = 
          downstreamλμs(bix, Ξ, idf, 0.0, NaN, NaN, NaN, NaN)
      end

      llc, prc, ddλ, ssλ, ssμ, mc, λi, μi = 
        _stem_update!(ξi, eds, λ1, λ2, μ1, μ2, α, σλ, σμ, llc, prc, 
          ddλ, ssλ, ssμ, mc, th, λ0_prior, μ0_prior, surv)

      # set new λ & μ downstream, if necessary
      setdownstreamλμ!(λi, μi, bix, Ξ, idf)

      # if there are speciation events in stem branch
      if !istip(ξi)
        eds, λ1, λ2, μ1, μ2 = 
          downstreamλμs(bix, Ξ, idf, 0.0, NaN, NaN, NaN, NaN)

        # updates within the parent branch
        llc, ddλ, ssλ, ssμ, λx, μx = 
          _update_internal!(ξi.d1, bi, eas, λa, μa, α, σλ, σμ, 
            eds, λ1, λ2, μ1, μ2, llc, ddλ, ssλ, ssμ, false)
        llc, ddλ, ssλ, ssμ, λx, μx = 
          _update_internal!(ξi.d2, bi, eas, λa, μa, α, σλ, σμ, 
            eds, λ1, λ2, μ1, μ2, llc, ddλ, ssλ, ssμ, false)
      end

    # if *not* root
    else

      # find cladogenetic ancestor
      eas, λa, μa, il = upstreamλμ(ia, Ξ, idf, 0.0, λa, μa)

      # check if mid branch does not lead to the stem branch
      if pa(idf[il]) > 0

        # it non-terminal branch
        eds, λ1, λ2, μ1, μ2 = 0.0, NaN, NaN, NaN, NaN
        if !it
          # if cladogenetic branch
          if i2 > 0
            ξ1, ξ2 = Ξ[i1], Ξ[i2]
            eds, λ1, λ2, μ1, μ2 = 0.0, lλ(ξ1), lλ(ξ2), lμ(ξ1), lμ(ξ2)
          # if mid branch
          else
            eds, λ1, λ2, μ1, μ2 = 
              downstreamλμs(i1, Ξ, idf, 0.0, NaN, NaN, NaN, NaN)
          end
        end

        ll0 = llc

        # updates within the parent branch
        llc, ddλ, ssλ, ssμ, λx, μx = 
          _update_internal!(ξi, bi, eas, λa, μa, α, σλ, σμ, eds, λ1, λ2, μ1, μ2,
            llc, ddλ, ssλ, ssμ, it)

        # if update, update up- and down-stream
        if ll0 != llc
          λi, μi = lλ(ξi), lμ(ξi)
          setupstreamλμ!(λi, μi, ia, Ξ, idf)
          if !it && iszero(i2)
            lξi = fixtip(ξi)
            setdownstreamλμ!(lλ(lξi), lμ(lξi), i1, Ξ, idf)
          end
        end
      end
    end
  end

  return llc, prc, ddλ, ssλ, ssμ, mc
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{cTbd},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               σμ   ::Float64,
               llc  ::Float64,
               ddλ  ::Float64,
               ssλ  ::Float64,
               ssμ  ::Float64,
               ns   ::Float64,
               ne   ::Float64,
               λfs  ::Vector{Float64},
               μfs  ::Vector{Float64})

Forward simulation proposal function for pure birth diffusion.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{cTbd},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    σμ   ::Float64,
                    llc  ::Float64,
                    ddλ  ::Float64,
                    ssλ  ::Float64,
                    ssμ  ::Float64,
                    ns   ::Float64,
                    ne   ::Float64,
                    λfs  ::Vector{Float64},
                    μfs  ::Vector{Float64})

  bi = idf[bix]
  ξc = Ξ[bix]
  ia = pa(bi)

  λa = μa = NaN
  # if following a speciation event
  if ia > 0 
    ba = idf[ia]
    if d2(ba) > 0
      λa, μa = λt(ba), μt(ba)
    end
  end

  ddλr = ssλr = ssμr = zero(Float64)
  llr  = NaN
  # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, ξc, λa, μa, α, σλ, σμ)

  # if mid
  elseif iszero(d2(bi))
    ξp, llr, ddλr, ssλr, ssμr = 
      fsbi_m(bi, idf, ξc, Ξ, λa, μa, α, σλ, σμ, λfs, μfs)

  # if internal
  else
    if e(bi) > 0.0
      ξ1, ξ2 = Ξ[d1(bi)], Ξ[d2(bi)]
      ξp, llr, ddλr, ssλr, ssμr = 
        fsbi_i(bi, ξc, λa, lλ(ξ1), lλ(ξ2), μa, lμ(ξ1), lμ(ξ2), 
          α, σλ, σμ, λfs, μfs)
    end
  end

  # if accepted
  if isfinite(llr)

    llc, ddλ, ssλ, ssμ, ns, ne = 
      llik_cladsbd_track!(ξc, α, σλ, σμ, llc, ddλ, ssλ, ssμ, ns, ne, -)
    llc, ddλ, ssλ, ssμ, ns, ne = 
      llik_cladsbd_track!(ξp, α, σλ, σμ, llc, ddλ, ssλ, ssμ, ns, ne, +)

    # first change from ancestor
    if !isnan(λa)
      λp, λc, = lλ(ξp), lλ(ξc)
      μp, μc, = lμ(ξp), lμ(ξc)
      llc += llrdnorm_x(λp, λc, λa + α, σλ^2) + llrdnorm_x(μp, μc, μa, σμ^2)
      ddλ += λp - λc
      ssλ += 0.5*((λp - λa - α)^2 - (λc - λa - α)^2)
      ssμ += 0.5*((μp - μa)^2 - (μc - μa)^2)
    end

    # update quantities
    ddλ += ddλr
    ssλ += ssλr
    ssμ += ssμr
    llc += llr

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ssλ, ssμ, ns, ne
end




"""
    fsbi_t(bi::iBffs,
           ξi::cTbd,
           λa::Float64,
           μa::Float64,
           α ::Float64,
           σλ::Float64,
           σμ::Float64)

Forward simulation for terminal branch `bi`.
"""
function fsbi_t(bi::iBffs,
                ξi::cTbd,
                λa::Float64,
                μa::Float64,
                α ::Float64,
                σλ::Float64,
                σμ::Float64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # if does **not** come from a cladogenetic event
  if isnan(λa)
    λi, μi = lλ(ξi), lμ(ξi)
  else
    λi = rnorm(λa + α, σλ)
    μi = rnorm(μa,     σμ)
  end

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_cladsbd_t(e(bi), λi, μi, α, σλ, σμ, lc, lU, iρi, 0, 1, 500)

  if na > 0 && isfinite(llr)

    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, NaN
  end
end




"""
    fsbi_m(bi ::iBffs,
           idf::Vector{iBffs},
           ξi ::cTbd,
           Ξ  ::Vector{cTbd},
           λa ::Float64,
           μa ::Float64,
           α  ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           λfs::Vector{Float64},
           μfs::Vector{Float64})

Forward simulation for mid branch `bi`.
"""
function fsbi_m(bi ::iBffs,
                idf::Vector{iBffs},
                ξi ::cTbd,
                Ξ  ::Vector{cTbd},
                λa ::Float64,
                μa ::Float64,
                α  ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                λfs::Vector{Float64},
                μfs::Vector{Float64})

  # if does **not** come from a cladogenetic event
  if isnan(λa)
    λi, μi = lλ(ξi), lμ(ξi)
  else
    λi = rnorm(λa + α, σλ)
    μi = rnorm(μa,     σμ)
  end

  # forward simulation during branch length
  empty!(λfs)
  empty!(μfs)
  t0, na, nn = _sim_cladsbd_i(e(bi), λi, μi, α, σλ, σμ, 0, 1, 500, λfs, μfs)

  if na < 1 || nn > 499
    return t0, NaN, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr  = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # search for next lλ1 and lλ2 if the exist
  i1 = d1(bi)
  eds, λ1, λ2, μ1, μ2 = downstreamλμs(i1, Ξ, idf, 0.0, NaN, NaN, NaN, NaN)

  ## choose most likely lineage to fix
  # if downstream is tip
  ddr = ssλr = ssμr = 0.0
  if isnan(λ1)
    wt, λp, μp, pp, λc, μc, pc, acr = wfix_m(ξi, e(bi), λfs, μfs, eds, acr)
  # if downstream is cladogenetic
  else
    wt, λp, μp, pp, λc, μc, pc, acr, ddr, ssλr, ssμr = 
      wfix_m(ξi, e(bi), λfs, μfs, eds, λ1, λ2, μ1, μ2, α, σλ, σμ, acr)
  end

  if lU < acr

    # fix the tip
    if wt <= div(na,2)
      fixtip1!(t0, wt, 0)
    else
      fixtip2!(t0, na - wt + 1, 0)
    end

    # simulated remaining tips until the present
    t0, na, nn, acr =
      tip_sims!(t0, tf(bi), α, σλ, σμ, acr, lU, iρi, na, nn)

    if lU < acr

      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc)
      if isfinite(λ1)
        llr += λp - λc
      end
      setni!(bi, na)                     # set new ni

      # downstream change
      setdownstreamλμ!(λp, μp, i1, Ξ, idf)

      return t0, llr, ddr, ssλr, ssμr
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    wfix_m(ξi ::cTbd,
           ei ::Float64,
           λfs::Vector{Float64},
           μfs::Vector{Float64},
           eds::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for middle `i` branches with downstream **tips**.
"""
function wfix_m(ξi ::cTbd,
                ei ::Float64,
                λfs::Vector{Float64},
                μfs::Vector{Float64},
                eds::Float64,
                acr::Float64)

  # select best from proposal
  sp, wt, λp, μp, pp = 0.0, 0, NaN, NaN, -Inf
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    p   = exp(- eds * (exp(λfi) + exp(μfi)))
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      μp  = μfi
      wt  = i
    end
  end

  # extract current λs and μs at time `t` and estimate ratio
  empty!(λfs)
  empty!(μfs)
  λc, μc = _λμat!(ξi, ei, λfs, μfs, 0.0, NaN, NaN)

  sc, pc = 0.0, NaN
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    p   = exp(- eds * (exp(λfi) + exp(μfs[i])))
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp/sc)

  return wt, λp, μp, pp, λc, μc, pc, acr
end




"""
    wfix_m(ξi ::cTbd,
           ei ::Float64,
           λfs::Vector{Float64},
           μfs::Vector{Float64},
           eds::Float64,
           λ1 ::Float64,
           λ2 ::Float64,
           μ1 ::Float64, 
           μ2 ::Float64,
           α  ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for middle `i` branches with downstream **cladogenetic** daughters.
"""
function wfix_m(ξi ::cTbd,
                ei ::Float64,
                λfs::Vector{Float64},
                μfs::Vector{Float64},
                eds::Float64,
                λ1 ::Float64,
                λ2 ::Float64,
                μ1 ::Float64, 
                μ2 ::Float64,
                α  ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                acr::Float64)

  # select best from proposal
  sp, wt, λp, μp, pp = 0.0, 0, NaN, NaN, -Inf
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    p   = dnorm2(λ1, λ2, λfi + α, σλ) * dnorm2(μ1, μ2, μfi, σμ) *
          exp(- eds * (exp(λfi) + exp(μfi)))
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      μp  = μfi
      wt  = i
    end
  end

  # extract current λs and μs at time `t` and estimate ratio
  empty!(λfs)
  empty!(μfs)
  λc, μc = _λμat!(ξi, ei, λfs, μfs, 0.0, NaN, NaN)

  sc, pc = 0.0, NaN
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    p   = dnorm2(λ1, λ2, λfi + α, σλ) * dnorm2(μ1, μ2, μfi, σμ) *
          exp(- eds * (exp(λfi) + exp(μfi)))
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood and acceptance ratio
  acr += log(sp/sc) + λp - λc
  ddr  = 2.0*(λc - λp)
  ssλr = 0.5*((λ1 - λp - α)^2 + (λ2 - λp - α)^2 - 
              (λ1 - λc - α)^2 - (λ2 - λc - α)^2)
  ssμr = 0.5*((μ1 - μp)^2 + (μ2 - μp)^2 - (μ1 - μc)^2 - (μ2 - μc)^2)

  return wt, λp, μp, pp, λc, μc, pc, acr, ddr, ssλr, ssμr
end




"""
    fsbi_i(bi ::iBffs,
           ξi ::cTbd,
           λa ::Float64,
           μa ::Float64,
           λ1 ::Float64,
           λ2 ::Float64,
           μ1 ::Float64,
           μ2 ::Float64,
           α  ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           λfs::Vector{Float64},
           μfs::Vector{Float64})

Forward simulation for internal branch `bi`
"""
function fsbi_i(bi ::iBffs,
                ξi ::cTbd,
                λa ::Float64,
                λ1 ::Float64,
                λ2 ::Float64,
                μa ::Float64,
                μ1 ::Float64,
                μ2 ::Float64,
                α  ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                λfs::Vector{Float64},
                μfs::Vector{Float64})

  # if does **not** come from a cladogenetic event
  if isnan(λa)
    λi, μi = lλ(ξi), lμ(ξi)
  else
    λi = rnorm(λa + α, σλ)
    μi = rnorm(μa,     σμ)
  end

  # forward simulation during branch length
  empty!(λfs)
  empty!(μfs)
  t0, na, nn = _sim_cladsbd_i(e(bi), λi, μi, α, σλ, σμ, 0, 1, 500, λfs, μfs)

  if na < 1 || nn > 499
    return t0, NaN, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr  = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # choose most likely lineage to fix
  wt, λp, μp, pp, λc, μc, pc, acr, ddr, ssλr, ssμr= 
    wfix_i(ξi, e(bi), λfs, μfs, λ1, λ2, μ1, μ2, α, σλ, σμ, acr)

  if lU < acr

    # fix the tip
    if wt <= div(na,2)
      fixtip1!(t0, wt, 0)
    else
      fixtip2!(t0, na - wt + 1, 0)
    end

    # simulated remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, acr, lU, iρi, na, nn)
    end

    if lU < acr
      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc) + λp - λc
      setni!(bi, na)                     # set new ni
      setλt!(bi, λp)                     # set new λt
      setμt!(bi, μp)                     # set new λt

      return t0, llr, ddr, ssλr, ssμr
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    wfix_i(ξi ::cTbd,
           ei ::Float64,
           λfs::Vector{Float64},
           μfs::Vector{Float64},
           λ1 ::Float64,
           λ2 ::Float64,
           μ1 ::Float64,
           μ2 ::Float64,
           α  ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for bifurcating `i` branches.
"""
function wfix_i(ξi ::cTbd,
                ei ::Float64,
                λfs::Vector{Float64},
                μfs::Vector{Float64},
                λ1 ::Float64,
                λ2 ::Float64,
                μ1 ::Float64,
                μ2 ::Float64,
                α  ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                acr::Float64)

  # select best from proposal
  sp, wt, λp, μp, pp = 0.0, 0, NaN, NaN, -Inf
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    p   = dnorm2(λ1, λ2, λfi + α, σλ) * dnorm2(μ1, μ2, μfi, σμ)
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      μp  = μfi
      wt  = i
    end
  end

  # extract current λs and μs at time `t` and estimate ratio
  empty!(λfs)
  empty!(μfs)
  λc, μc = _λμat!(ξi, ei, λfs, μfs, 0.0, NaN, NaN)

  sc, pc = 0.0, NaN
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    p   = dnorm2(λ1, λ2, λfs[i] + α, σλ) * dnorm2(μ1, μ2, μfs[i], σμ)
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood and acceptance ratio
  acr += log(sp/sc) + λp - λc
  ddr  = 2.0*(λc - λp)
  ssλr = 0.5*((λ1 - λp - α)^2 + (λ2 - λp - α)^2 - 
              (λ1 - λc - α)^2 - (λ2 - λc - α)^2)
  ssμr = 0.5*((μ1 - μp)^2 + (μ2 - μp)^2 - (μ1 - μc)^2 - (μ2 - μc)^2)

  return wt, λp, μp, pp, λc, μc, pc, acr, ddr, ssλr, ssμr
end




"""
    tip_sims!(tree::cTbd,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::cTbd,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

 if lU < lr && na < 500

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nn, lr =
          _sim_cladsbd_it(t, lλ(tree), lμ(tree), α, σλ, σμ, 
            lr, lU, iρi, na-1, nn, 500)

        if isnan(lr) || nn > 499
          return tree, na, nn, NaN
        end

        setproperty!(tree, :iμ, isextinct(stree))
        sete!(tree, e(tree) + e(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr = 
        tip_sims!(tree.d1, t, α, σλ, σμ, lr, lU, iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, α, σλ, σμ, lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end





