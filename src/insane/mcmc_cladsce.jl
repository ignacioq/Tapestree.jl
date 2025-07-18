#=

clads constant-extinction MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    insane_cladsce(tree    ::sT_label;
                   λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                   α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                   σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                   μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                   niter   ::Int64                 = 1_000,
                   nthin   ::Int64                 = 10,
                   nburn   ::Int64                 = 200,
                   nflush  ::Int64                 = nthin,
                   ofile   ::String                = string(homedir(), "/cladsce"),
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

Run insane for clads contant-extinction.
"""
function insane_cladsce(tree    ::sT_label;
                        λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                        α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                        σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                        μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                        niter   ::Int64                 = 1_000,
                        nthin   ::Int64                 = 10,
                        nburn   ::Int64                 = 200,
                        nflush  ::Int64                 = nthin,
                        ofile   ::String                = string(homedir(), "/cladsce"),
                        λi      ::Float64               = NaN,
                        αi      ::Float64               = 0.0,
                        σλi     ::Float64               = 0.1,
                        μi      ::Float64               = NaN,
                        ϵi      ::Float64               = 0.2,
                        pupdp   ::NTuple{6,Float64}     = (1e-3, 1e-3, 1e-3, 1e-4, 0.2, 0.2),
                        prints  ::Int64                 = 5,
                        stn     ::Float64               = 0.5,
                        survival::Bool                  = true,
                        mxthf   ::Float64               = 0.1,
                        tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)

  # turn to logarithmic terms
  λ0_prior = (log(λ0_prior[1]), 2.0*log(λ0_prior[2]))

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
  Ξ = make_Ξ(idf, λc, cTce)

  # survival
  mc = m_surv_cladsce(th, log(λc), αi, σλi, μc, 1_000, surv)

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σ, 3: scale, 4: gbm, 5: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running clads with constant μ"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, μc, mc, ns, ne, ddλ, ssλ, L, stn =
    mcmc_burn_cladsce(Ξ, idf, λ0_prior, α_prior, σλ_prior, μ_prior, nburn, 
      αi, σλi, μc, mc, th, rmλ, surv, stn, inodes, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cladsce(Ξ, idf, llc, prc, αc, σλc, μc, mc, th, rmλ, surv, ns, ne, 
      ddλ, ssλ, L, stn, λ0_prior, α_prior, σλ_prior, μ_prior, 
      inodes, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_cladsce(Ξ       ::Vector{cTce},
                      idf     ::Vector{iBffs},
                      λ0_prior::NTuple{2,Float64},
                      α_prior ::NTuple{2,Float64},
                      σλ_prior::NTuple{2,Float64},
                      μ_prior ::NTuple{2,Float64},
                      nburn   ::Int64,
                      αc      ::Float64,
                      σλc     ::Float64,
                      μc      ::Float64,
                      mc      ::Float64,
                      th      ::Float64,
                      rmλ     ::Float64,
                      surv    ::Int64,
                      stn     ::Float64,
                      inodes  ::Array{Int64,1},
                      pup     ::Array{Int64,1},
                      prints  ::Int64)

MCMC burn-in chain for `pbd`.
"""
function mcmc_burn_cladsce(Ξ       ::Vector{cTce},
                           idf     ::Vector{iBffs},
                           λ0_prior::NTuple{2,Float64},
                           α_prior ::NTuple{2,Float64},
                           σλ_prior::NTuple{2,Float64},
                           μ_prior ::NTuple{2,Float64},
                           nburn   ::Int64,
                           αc      ::Float64,
                           σλc     ::Float64,
                           μc      ::Float64,
                           mc      ::Float64,
                           th      ::Float64,
                           rmλ     ::Float64,
                           surv    ::Int64,
                           stn     ::Float64,
                           inodes  ::Array{Int64,1},
                           pup     ::Array{Int64,1},
                           prints  ::Int64)

  # starting likelihood and prior
  lλ0 = lλ(Ξ[1])
  llc = llik_clads(Ξ, idf, αc, σλc, μc) - rmλ*lλ0 + log(mc) + prob_ρ(idf)
  prc = logdnorm(lλ0,       λ0_prior[1], λ0_prior[2])   +
        logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])   +
        logdnorm(αc,         α_prior[1],  α_prior[2]^2) +
        logdgamma(μc,        μ_prior[1],  μ_prior[2])

  L   = treelength(Ξ)      # tree length
  nin = lastindex(inodes)                       # number of internal nodes
  el  = lastindex(idf)                          # number of branches
  ns  = sum(x -> Float64(d2(x) > 0), idf) - rmλ # number of speciation events in likelihood
  ne  = 0.0                                     # number of extinction events
  λfs = Float64[]

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ = _ss_dd(Ξ, idf, lλ, αc)

  # for scale tuning
  ltn = zero(Int64)
  lup = lac = zero(Float64)

  pbar = Progress(nburn, dt = prints, desc = "burning mcmc...", barlen = 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update drift
      if pupi === 1

        llc, prc, αc, mc = 
          update_α!(αc, lλ(Ξ[1]), σλc, μc, 2.0*(ns + rmλ), ddλ, llc, prc, 
            mc, th, surv, α_prior)

        # update ssλ with new drift `α`
        ssλ = _ss(Ξ, idf, lλ, αc)

      # update diffusion
      elseif pupi === 2

        llc, prc, σλc, mc = 
          update_σ!(σλc, lλ(Ξ[1]), σλc, μc, ssλ, 2.0*(ns + rmλ), llc, prc, 
            mc, th, surv, σλ_prior)

      # update extinction
      elseif pupi === 3

        llc, prc, μc, mc =
          update_μ!(μc, lλ(Ξ[1]), αc, σλc, llc, prc, ne, L, mc, th, surv,
             μ_prior)

      # update scale
      elseif pupi === 4

        llc, prc, mc, acc = 
          update_scale!(Ξ, idf, αc, σλc, μc, llc, prc, ns, stn, 
            mc, th, surv, λ0_prior)

        lac += acc
        lup += 1.0

      # update internal
      elseif pupi === 5

        bix = inodes[fIrand(nin) + 1]

        llc, prc, ddλ, ssλ, mc =
          update_internal!(bix, Ξ, idf, αc, σλc, μc, llc, prc, ddλ, ssλ, 
            mc, th, λ0_prior, surv)

      # forward simulation
      else

        bix = fIrand(el) + 1

        llc, ddλ, ssλ, ns, ne, L =
          update_fs!(bix, Ξ, idf, αc, σλc, μc, llc, ddλ, ssλ, ns, 
            ne, L, λfs)
      end
    end

    ltn += 1.0
    if ltn === 100
      stn = tune(stn, lac/lup)
      ltn = zero(Int64)
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, μc, mc, ns, ne, ddλ, ssλ, L, stn
end




"""
    mcmc_cladsce(Ξ       ::Vector{cTce},
                 idf     ::Vector{iBffs},
                 llc     ::Float64,
                 prc     ::Float64,
                 αc      ::Float64,
                 σλc     ::Float64,
                 μc     ::Float64,
                 mc      ::Float64,
                 th      ::Float64,
                 rmλ     ::Float64,
                 surv    ::Int64,
                 ns      ::Float64,
                 ne      ::Float64,
                 L       ::Float64,
                 stn     ::Float64,
                 λ0_prior::NTuple{2,Float64},
                 α_prior ::NTuple{2,Float64},
                 σλ_prior::NTuple{2,Float64},
                 μ_prior ::NTuple{2,Float64},
                 inodes  ::Array{Int64,1},
                 pup     ::Vector{Int64},
                 niter   ::Int64,
                 nthin   ::Int64,
                 nflush  ::Int64,
                 ofile   ::String,
                 prints  ::Int64)

MCMC chain for pure-birth diffusion.
"""
function mcmc_cladsce(Ξ       ::Vector{cTce},
                      idf     ::Vector{iBffs},
                      llc     ::Float64,
                      prc     ::Float64,
                      αc      ::Float64,
                      σλc     ::Float64,
                      μc     ::Float64,
                      mc      ::Float64,
                      th      ::Float64,
                      rmλ     ::Float64,
                      surv    ::Int64,
                      ns      ::Float64,
                      ne      ::Float64,
                      ddλ     ::Float64,
                      ssλ     ::Float64,
                      L       ::Float64,
                      stn     ::Float64,
                      λ0_prior::NTuple{2,Float64},
                      α_prior ::NTuple{2,Float64},
                      σλ_prior::NTuple{2,Float64},
                      μ_prior ::NTuple{2,Float64},
                      inodes  ::Array{Int64,1},
                      pup     ::Vector{Int64},
                      niter   ::Int64,
                      nthin   ::Int64,
                      nflush  ::Int64,
                      ofile   ::String,
                      prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin = lit = sthin = zero(Int64)

  nin = lastindex(inodes)  # number of internal nodes
  el  = lastindex(idf)     # number of branches

  # parameter results
  r   = Array{Float64,2}(undef, nlogs, 7)

  λfs   = Float64[]
  treev = cTce[]  # make Ξ vector
  io    = IOBuffer() # buffer 

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\talpha\tsigma_lambda\tmu\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αc = αc, σλc = σλc, μc = μc, mc = mc, ns = ns, ne = ne, L = L, ssλ = ssλ, ddλ = ddλ, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update drift
            if pupi === 1

              llc, prc, αc, mc = 
                update_α!(αc, lλ(Ξ[1]), σλc, μc, 2.0*(ns + rmλ), ddλ, llc, prc, 
                  mc, th, surv, α_prior)

              # update ssλ with new drift `α`
              ssλ = _ss(Ξ, idf, lλ, αc)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, μc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update diffusion rate
            elseif pupi === 2

              llc, prc, σλc, mc = 
                update_σ!(σλc, lλ(Ξ[1]), σλc, μc, ssλ, 2.0*(ns + rmλ), llc, prc, 
                  mc, th, surv, σλ_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, μc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update extinction
            elseif pupi === 3

              llc, prc, μc, mc =
                update_μ!(μc, lλ(Ξ[1]), αc, σλc, llc, prc, ne, L, mc, th, surv,
                   μ_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, μc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update scale
            elseif pupi === 4

              llc, prc, mc, acc = 
                update_scale!(Ξ, idf, αc, σλc, μc, llc, prc, ns, stn, 
                  mc, th, surv, λ0_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, μc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update internal λ
            elseif pupi === 5

              bix = inodes[fIrand(nin) + 1]

              llc, prc, ddλ, ssλ, mc =
                update_internal!(bix, Ξ, idf, αc, σλc, μc, llc, prc, ddλ, ssλ, 
                  mc, th, λ0_prior, surv)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, μc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update by forward simulation
            else

              bix = fIrand(el) + 1

              llc, ddλ, ssλ, ns, ne, L =
                update_fs!(bix, Ξ, idf, αc, σλc, μc, llc, ddλ, ssλ, ns, 
                  ne, L, λfs)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, μc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
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
              r[lit,5] = αc
              r[lit,6] = σλc
              r[lit,7] = μc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
                  exp(lλ(Ξ[1])),'\t', αc, '\t', σλc, '\t', μc, '\n')
            flush(of)
            ibuffer(io, couple(Ξ, idf, 1))
            write(io, '\n')
            write(tf, take!(io))
            flush(tf)
            sthin = zero(Int64)
          end

          next!(pbar)
        end

        return r, treev
      end
    end
  end
end




"""
    update_α!(αc     ::Float64,
              λ0     ::Float64,
              σλ     ::Float64,
              μ      ::Float64,
              L      ::Float64,
              ddλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              mc     ::Float64,
              th     ::Float64,
              crown  ::Int64,
              δt     ::Float64,
              srδt   ::Float64,
              α_prior::NTuple{2,Float64})

Gibbs update for `α`.
"""
function update_α!(αc     ::Float64,
                   λ0     ::Float64,
                   σλ     ::Float64,
                   μ      ::Float64,
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

  mp  = m_surv_cladsce(th, λ0, αp, σλ, μ, 1_000, surv)
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
    update_σ!(σλc     ::Float64,
              λ0      ::Float64,
              α       ::Float64,
              μ       ::Float64,
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
function update_σ!(σλc     ::Float64,
                   λ0      ::Float64,
                   α       ::Float64,
                   μ       ::Float64,
                   ssλ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   surv   ::Int64,
                   σλ_prior::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)
  σλp  = sqrt(σλp2)

  mp  = m_surv_cladsce(th, λ0, α, σλp, μ, 1_000, surv)
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
    update_μ!(μc     ::Float64,
              λ0     ::Float64,
              α      ::Float64,
              σλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              ne     ::Float64,
              L      ::Float64,
              mc     ::Float64,
              th     ::Float64,
              surv   ::Int64,
              μ_prior::NTuple{2,Float64})

Gibbs-MH update for `μ`.
"""
function update_μ!(μc     ::Float64,
                   λ0     ::Float64,
                   α      ::Float64,
                   σλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   ne     ::Float64,
                   L      ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   surv   ::Int64,
                   μ_prior::NTuple{2,Float64})

  μp  = randgamma(μ_prior[1] + ne, μ_prior[2] + L)

  mp  = m_surv_cladsce(th, λ0, α, σλ, μp, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += ne * log(μp/μc) + L * (μc - μp) + llr
    prc += llrdgamma(μp, μc, μ_prior[1], μ_prior[2])
    μc   = μp
    mc   = mp
  end

  return llc, prc, μc, mc
end




"""
    update_scale!(α       ::Float64,
                  σλ      ::Float64,
                  μ       ::Float64,
                  Ξ       ::Vector{cTce},
                  idf     ::Vector{iBffs},
                  α       ::Float64,
                  σλ      ::Float64,
                  μ       ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ns      ::Float64,
                  stn     ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  surv    ::Int64,
                  λ0_prior::NTuple{2,Float64})

Update scale for speciation.
"""
function update_scale!(Ξ       ::Vector{cTce},
                       idf     ::Vector{iBffs},
                       α       ::Float64,
                       σλ      ::Float64,
                       μ       ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ns      ::Float64,
                       stn     ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       surv    ::Int64,
                       λ0_prior::NTuple{2,Float64})

  # sample log(scaling factor)
  s = randn()*stn

  lλ0 = lλ(Ξ[1])

  # prior ratio
  prr = llrdnorm_x(lλ0 + s, lλ0, λ0_prior[1], λ0_prior[2]) 

  # survival
  mp  = m_surv_cladsce(th, lλ0 + s, α, σλ, μ, 1_000, surv)

  # likelihood ratio
  ir  = _ir(Ξ)
  llr = ns * s + (1.0 - exp(s)) * ir + log(mp/mc)

  acc = 0.0
  if -randexp() < llr + prr
    acc += 1.0
    llc += llr
    prc += prr
    mc  = mp
    scale_rateλ!(Ξ, s)
    scale_rate!(idf, s)
  end

  return llc, prc, mc, acc
end




"""
    update_internal!(bix     ::Int64,
                     Ξ       ::Vector{cTce},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ddλ     ::Float64,
                     ssλ     ::Float64,
                     mc      ::Float64,
                     th      ::Float64,
                     λ0_prior::NTuple{2,Float64},
                     surv    ::Int64)

Make an update for an internal branch and its descendants.
"""
function update_internal!(bix     ::Int64,
                          Ξ       ::Vector{cTce},
                          idf     ::Vector{iBffs},
                          α       ::Float64,
                          σλ      ::Float64,
                          μ       ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          ddλ     ::Float64,
                          ssλ     ::Float64,
                          mc      ::Float64,
                          th      ::Float64,
                          λ0_prior::NTuple{2,Float64},
                          surv    ::Int64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  i2   = d2(bi)
  ia   = pa(bi)
  ξ1   = Ξ[i1]
  ξ2   = Ξ[i2]
  root = iszero(ia)
  λa   = NaN  # ancestral speciation

  # if crown root
  if root && iszero(e(ξi))
    llc, prc, ddλ, ssλ, mc =
      _crown_update!(ξi, ξ1, ξ2, α, σλ, μ, llc, prc, ddλ, ssλ, mc, th, 
        λ0_prior, surv)
    λa = lλ(ξi)
    setλt!(bi, λa)
  else
    # if stem
    if root
      if istip(ξi)
        llc, prc, ddλ, ssλ, mc = 
          _stem_update!(ξi, lλ(ξ1), lλ(ξ2), 
            α, σλ, μ, llc, prc, ddλ, ssλ, mc, th, λ0_prior, surv)
        λa = lλ(ξi)
      else
        llc, prc, ddλ, ssλ, mc = 
          _stem_update!(ξi, lλ(ξ1.d1), lλ(ξ2.d2),
            α, σλ, μ, llc, prc, ddλ, ssλ, mc, th, λ0_prior, surv)

        # updates within the stem daughter branches
        llc, ddλ, ssλ = 
          _update_internal!(ξi.d1, lλ(ξi), α, σλ, llc, ddλ, ssλ, false)
        llc, ddλ, ssλ = 
          _update_internal!(ξi.d2, lλ(ξi), α, σλ, llc, ddλ, ssλ, false)

        # get fixed tip and ancestral rate
        lξi, λa = fixtip(ξi, λa)

        # make fixed branch update
        llc, ddλ, ssλ, λa = 
          update_triad!(lξi, ξ1, ξ2, λa, α, σλ, llc, ddλ, ssλ)
      end
    else
      λa = λt(idf[ia])

      # updates within the parent branch
      llc, ddλ, ssλ, λx = 
        _update_internal!(ξi, λa, α, σλ, llc, ddλ, ssλ, false)

      # get fixed tip and ancestral rate
      lξi, λa = fixtip(ξi, λa)

      # make fixed branch update
      llc, ddλ, ssλ, λa = 
        update_triad!(lξi, ξ1, ξ2, λa, α, σλ, llc, ddλ, ssλ)
    end

    # set fixed `λ(t)` in branch
    setλt!(bi, λa)
  end

  # # carry on updates in the daughters
  llc, ddλ, ssλ, λx = 
    _update_internal!(ξ1, λa, α, σλ, llc, ddλ, ssλ, iszero(d1(idf[i1])))
  llc, ddλ, ssλ, λx = 
    _update_internal!(ξ2, λa, α, σλ, llc, ddλ, ssλ, iszero(d1(idf[i2])))

  return llc, prc, ddλ, ssλ, mc
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{cTce},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               μ    ::Float64,
               llc  ::Float64,
               ddλ  ::Float64,
               ssλ  ::Float64,
               ns   ::Float64,
               ne   ::Float64,
               L    ::Float64,
               λfs  ::Vector{Float64})

Forward simulation proposal function for pure birth diffusion.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{cTce},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    μ    ::Float64,
                    llc  ::Float64,
                    ddλ  ::Float64,
                    ssλ  ::Float64,
                    ns   ::Float64,
                    ne   ::Float64,
                    L    ::Float64,
                    λfs  ::Vector{Float64})

  bi = idf[bix]
  ξc = Ξ[bix]
  ia = pa(bi)

  λa = NaN
  # if following a speciation event
  if ia > 0 && d2(idf[ia]) > 0
    λa = λt(idf[ia])
  end

  # if terminal
  ssλr = ddλr = zero(Float64)
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, λa, α, σλ, μ)

  # if mid
  elseif iszero(d2(bi))
    ξp, llr =fsbi_m(bi, ξc, λa, lλ(Ξ[d1(bi)]), α, σλ, μ, λfs)

  # if internal
  else
    ξp, llr, ddλr, ssλr = 
      fsbi_i(bi, ξc, λa, lλ(Ξ[d1(bi)]), lλ(Ξ[d2(bi)]), α, σλ, μ, λfs)
  end

  # if accepted
  if isfinite(llr)

    llc, ddλ, ssλ, ns, ne, L = 
      llik_cladsce_track!(ξc, α, σλ, μ, llc, ddλ, ssλ, ns, ne, L, -)
    llc, ddλ, ssλ, ns, ne, L = 
      llik_cladsce_track!(ξp, α, σλ, μ, llc, ddλ, ssλ, ns, ne, L, +)

    # first change from ancestor
    if ia > 0 && d2(idf[ia]) > 0
      λp, λc = lλ(ξp), lλ(ξc)
      llc += llrdnorm_x(λp, λc, λa + α, σλ^2)
      ddλ += λp - λc
      ssλ += 0.5*((λp - λa - α)^2 - (λc - λa - α)^2)
    end

    # update quantities
    ddλ += ddλr
    ssλ += ssλr
    llc += llr

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ssλ, ns, ne, L
end




"""
    fsbi_t(bi::iBffs,
           λa::Float64,
           α ::Float64,
           σλ::Float64,
           μ ::Float64)

Forward simulation for terminal branch `bi`.
"""
function fsbi_t(bi::iBffs,
                λa::Float64,
                α ::Float64,
                σλ::Float64,
                μ ::Float64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_cladsce_t(e(bi), rnorm(λa + α, σλ), α, σλ, μ, lc, lU, iρi, 0, 1, 500)

  if na > 0 && isfinite(llr)
    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, NaN
  end
end




"""
    fsbi_m(bi  ::iBffs,
           ξc  ::cTce,
           λa  ::Float64,
           λ1  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           μ   ::Float64,
           λfs ::Vector{Float64})

Forward simulation for internal branch `bi`
"""
function fsbi_m(bi ::iBffs,
                ξc ::cTce,
                λa ::Float64,
                λ1 ::Float64,
                α  ::Float64,
                σλ ::Float64,
                μ  ::Float64,
                λfs::Vector{Float64})

  if isnan(λa)
    λi = lλ(ξc)
  else
    λi = rnorm(λa + α, σλ)
  end

  empty!(λfs)

  # forward simulation during branch length
  t0, na, nn = _sim_cladsce_i(e(bi), λi, α, σλ, μ, 0, 1, 500, λfs)

  if na < 1 || nn > 499
    return t0, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr  = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # choose most likely lineage to fix
  wt, λp, acr = wfix_m(t0, ξc, e(bi), λfs, λ1, α, σλ, acr)

  # fix the tip
  if wt <= div(na,2)
    bo, ix, λa, ei = fixtip1!(t0, wt, 0, NaN, NaN, λ1)
  else
    bo, ix, λa, ei = fixtip2!(t0, na - wt + 1, 0, NaN, NaN, λ1)
  end

  # if there were speciation events
  if !isnan(λa)
    acr += llrdnorm_x(λ1, λp, λa + α, σλ^2) + ei*(exp(λp) - exp(λ1))
  end

  if lU < acr

    # simulated remaining tips until the present
    t0, na, nn, acr =
      tip_sims!(t0, tf(bi), α, σλ, μ, acr, lU, iρi, na, nn)

    if lU < acr
      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
      setni!(bi, na)                     # set new ni

      return t0, llr
    end
  end

  return t0, NaN
end




"""
    wfix_m(ξp ::cTce,
           ξc ::cTce,
           ei ::Float64,
           λfs::Vector{Float64},
           λ1 ::Float64,
           α  ::Float64,
           σλ ::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for bifurcating `i` branches.
"""
function wfix_m(ξp ::cTce,
                ξc ::cTce,
                ei ::Float64,
                λfs::Vector{Float64},
                λ1 ::Float64,
                α  ::Float64,
                σλ ::Float64,
                acr::Float64)

  # select best from proposal
  sp, i, wt, λp, pp = 0.0, 0, 0, NaN, -Inf
  for λfi in λfs
    i += 1
    p  = dnorm(λ1, λfi, σλ)
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      wt  = i
    end
  end

  # extract current xis and estimate ratio
  empty!(λfs)
  λc = _λat!(ξi, ei, λfs, 0.0, NaN)

  sc, pc = 0.0, NaN
  for λfi in λfs
    p   = dnorm(λ1, λfi, σλ)
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log((pc * sp)/(pp * sc))

  return wt, λp, acr
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::cTce,
           λa  ::Float64,
           λ1  ::Float64,
           λ2  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           μ   ::Float64,
           λfs ::Vector{Float64})

Forward simulation for internal branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::cTce,
                λa  ::Float64,
                λ1  ::Float64,
                λ2  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                μ   ::Float64,
                λfs ::Vector{Float64})

  if isnan(λa)
    λi = lλ(ξc)
  else
    λi = rnorm(λa + α, σλ)
  end

  empty!(λfs)

  # forward simulation during branch length
  t0, na, nn = _sim_cladsce_i(e(bi), λi, α, σλ, μ, 0, 1, 500, λfs)

  if na < 1 || nn > 499
    return t0, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr  = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # choose most likely lineage to fix
  wt, λp, pp, λc, pc, acr, ddr, ssr = 
    wfix_i(ξc, e(bi), λfs, λ1, λ2, α, σλ, acr)

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
        tip_sims!(t0, tf(bi), α, σλ, μ, acr, lU, iρi, na, nn)
    end

    if lU < acr
      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc) + λp - λc
      setni!(bi, na)                     # set new ni
      setλt!(bi, λp)                     # set new λt

      return t0, llr, ddr, ssr
    end
  end

  return t0, NaN, NaN, NaN
end




"""
    wfix_i(ξi ::cTce,
           ei ::Float64,
           λfs::Vector{Float64},
           λ1 ::Float64,
           λ2 ::Float64,
           α  ::Float64,
           σλ ::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for bifurcating `i` branches.
"""
function wfix_i(ξi ::cTce,
                ei ::Float64,
                λfs::Vector{Float64},
                λ1 ::Float64,
                λ2 ::Float64,
                α  ::Float64,
                σλ ::Float64,
                acr::Float64)

  # select best from proposal
  sp, i, wt, λp, pp = 0.0, 0, 0, NaN, -Inf
  for λfi in λfs
    i += 1
    p   = dnorm2(λ1, λ2, λfi + α, σλ)
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      wt  = i
    end
  end

  # extract current xis and estimate ratio
  empty!(λfs)
  λc = _λat!(ξi, ei, λfs, 0.0, NaN)

  sc, pc = 0.0, NaN
  for λfi in λfs
    p   = dnorm2(λ1, λ2, λfi + α, σλ)
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood and acceptance ratio
  acr += log(sp/sc) + λp - λc
  ddr  = 2.0*(λc - λp)
  ssr  = 0.5*((λ1 - λp - α)^2 + (λ2 - λp - α)^2 - 
              (λ1 - λc - α)^2 - (λ2 - λc - α)^2)

  return wt, λp, pp, λc, pc, acr, ddr, ssr
end




"""
    tip_sims!(tree::cTce,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::cTce,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
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
          _sim_cladsce_it(t, lλ(tree), α, σλ, μ, lr, lU, iρi, na-1, nn, 500)

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
        tip_sims!(tree.d1, t, α, σλ, μ, lr, lU, iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, α, σλ, μ, lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end





