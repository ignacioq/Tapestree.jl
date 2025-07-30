#=

clads constant-turnover MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 25 07 2025
=#




"""
    insane_cladsct(tree    ::sT_label;
                   λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                   α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                   σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                   ϵ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                   niter   ::Int64                 = 1_000,
                   nthin   ::Int64                 = 10,
                   nburn   ::Int64                 = 200,
                   nflush  ::Int64                 = nthin,
                   ofile   ::String                = string(homedir(), "/cladsct"),
                   λi      ::Float64               = NaN,
                   αi      ::Float64               = 0.0,
                   σλi     ::Float64               = 0.1,
                   ϵi      ::Float64               = NaN,
                   ϵi      ::Float64               = 0.2,
                   pupdp   ::NTuple{6,Float64}     = (1e-3, 1e-3, 1e-3, 1e-4, 0.1, 0.2),
                   prints  ::Int64                 = 5,
                   stn     ::Float64               = 0.5,
                   survival::Bool                  = true,
                   mxthf   ::Float64               = 0.1,
                   tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for clads contant-extinction.
"""
function insane_cladsct(tree    ::sT_label;
                        λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                        α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                        σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                        ϵ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                        niter   ::Int64                 = 1_000,
                        nthin   ::Int64                 = 10,
                        nburn   ::Int64                 = 200,
                        nflush  ::Int64                 = nthin,
                        ofile   ::String                = string(homedir(), "/cladsct"),
                        λi      ::Float64               = NaN,
                        αi      ::Float64               = 0.0,
                        σλi     ::Float64               = 0.1,
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
  λc, ϵc = λi, ϵi
  if isnan(λi) || isnan(ϵi)
    λc, ϵc = moments(Float64(n), th, ϵi)
  end

  # make a decoupled tree
  Ξ = make_Ξ(idf, λc, cTct)

  # survival
  mc = m_surv_cladsct(th, log(λc), αi, σλi, ϵc, 1_000, surv)

  # parameter updates (1: α, 2: σ, 3: ϵ, 4: scale, 5: internal, 6: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running clads with constant turnover (μ(t) = ϵλ(t))"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, ϵc, mc, ns, ne, ddλ, ssλ, seλ, stn =
    mcmc_burn_cladsct(Ξ, idf, λ0_prior, α_prior, σλ_prior, ϵ_prior, nburn, 
      αi, σλi, ϵc, mc, th, rmλ, surv, stn, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cladsct(Ξ, idf, llc, prc, αc, σλc, ϵc, mc, th, rmλ, surv, ns, ne, 
      ddλ, ssλ, seλ, stn, λ0_prior, α_prior, σλ_prior, ϵ_prior, 
      pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_cladsct(Ξ       ::Vector{cTct},
                      idf     ::Vector{iBffs},
                      λ0_prior::NTuple{2,Float64},
                      α_prior ::NTuple{2,Float64},
                      σλ_prior::NTuple{2,Float64},
                      ϵ_prior ::NTuple{2,Float64},
                      nburn   ::Int64,
                      αc      ::Float64,
                      σλc     ::Float64,
                      ϵc      ::Float64,
                      mc      ::Float64,
                      th      ::Float64,
                      rmλ     ::Float64,
                      surv    ::Int64,
                      stn     ::Float64,
                      pup     ::Array{Int64,1},
                      prints  ::Int64)

MCMC burn-in chain for `pbd`.
"""
function mcmc_burn_cladsct(Ξ       ::Vector{cTct},
                           idf     ::Vector{iBffs},
                           λ0_prior::NTuple{2,Float64},
                           α_prior ::NTuple{2,Float64},
                           σλ_prior::NTuple{2,Float64},
                           ϵ_prior ::NTuple{2,Float64},
                           nburn   ::Int64,
                           αc      ::Float64,
                           σλc     ::Float64,
                           ϵc      ::Float64,
                           mc      ::Float64,
                           th      ::Float64,
                           rmλ     ::Float64,
                           surv    ::Int64,
                           stn     ::Float64,
                           pup     ::Array{Int64,1},
                           prints  ::Int64)

  # starting likelihood and prior
  lλ0 = lλ(Ξ[1])
  llc = llik_clads(Ξ, idf, αc, σλc, ϵc) - rmλ*lλ0 + log(mc) + prob_ρ(idf)
  prc = logdnorm(lλ0,       λ0_prior[1], λ0_prior[2])   +
        logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])   +
        logdnorm(αc,         α_prior[1],  α_prior[2]^2) +
        logdgamma(ϵc,        ϵ_prior[1],  ϵ_prior[2])

  el  = lastindex(idf)                          # number of branches
  ns  = sum(x -> Float64(d2(x) > 0), idf) - rmλ # number of speciation events in likelihood
  ne  = 0.0                                     # number of extinction events
  λfs = Float64[]

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, seλ = _dd_ss_seλ(Ξ, idf, αc)

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
          updatect_α!(αc, lλ(Ξ[1]), σλc, ϵc, 2.0*(ns + rmλ), ddλ, llc, prc, 
            mc, th, surv, α_prior)

        # update ssλ with new drift `α`
        ssλ = _ss(Ξ, idf, αc)

      # update diffusion
      elseif pupi === 2

        llc, prc, σλc, mc = 
          updatect_σ!(σλc, lλ(Ξ[1]), αc, ϵc, ssλ, 2.0*(ns + rmλ), llc, prc, 
            mc, th, surv, σλ_prior)

      # update extinction
      elseif pupi === 3

        llc, prc, ϵc, mc =
          update_ϵ!(ϵc, lλ(Ξ[1]), αc, σλc, llc, prc, ne, seλ, mc, th, 
            surv, ϵ_prior)

      # update scale
      elseif pupi === 4

        llc, prc, mc, seλ, acc = 
          update_scale!(Ξ, idf, αc, σλc, ϵc, llc, prc, ns, ne, stn, 
            mc, th, surv, λ0_prior)

        lac += acc
        lup += 1.0

      # update internal
      elseif pupi === 5

        bix = fIrand(el) + 1

        llc, prc, ddλ, ssλ, seλ, mc =
          update_internal!(bix, Ξ, idf, αc, σλc, ϵc, llc, prc, ddλ, ssλ, 
            seλ, mc, th, λ0_prior, surv)

      # forward simulation
      else

        bix = fIrand(el) + 1

        llc, ddλ, ssλ, seλ, ns, ne =
          update_fs!(bix, Ξ, idf, αc, σλc, ϵc, llc, ddλ, ssλ, seλ, ns, 
            ne, λfs)
      end
    end

    ltn += 1.0
    if ltn === 100
      stn = tune(stn, lac/lup)
      ltn = zero(Int64)
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, ϵc, mc, ns, ne, ddλ, ssλ, seλ, stn
end




"""
    mcmc_cladsct(Ξ       ::Vector{cTct},
                 idf     ::Vector{iBffs},
                 llc     ::Float64,
                 prc     ::Float64,
                 αc      ::Float64,
                 σλc     ::Float64,
                 ϵc      ::Float64,
                 mc      ::Float64,
                 th      ::Float64,
                 rmλ     ::Float64,
                 surv    ::Int64,
                 ns      ::Float64,
                 ne      ::Float64,
                 ddλ     ::Float64,
                 ssλ     ::Float64,
                 seλ     ::Float64,
                 stn     ::Float64,
                 λ0_prior::NTuple{2,Float64},
                 α_prior ::NTuple{2,Float64},
                 σλ_prior::NTuple{2,Float64},
                 ϵ_prior ::NTuple{2,Float64},
                 pup     ::Vector{Int64},
                 niter   ::Int64,
                 nthin   ::Int64,
                 nflush  ::Int64,
                 ofile   ::String,
                 prints  ::Int64)

MCMC chain for pure-birth diffusion.
"""
function mcmc_cladsct(Ξ       ::Vector{cTct},
                      idf     ::Vector{iBffs},
                      llc     ::Float64,
                      prc     ::Float64,
                      αc      ::Float64,
                      σλc     ::Float64,
                      ϵc      ::Float64,
                      mc      ::Float64,
                      th      ::Float64,
                      rmλ     ::Float64,
                      surv    ::Int64,
                      ns      ::Float64,
                      ne      ::Float64,
                      ddλ     ::Float64,
                      ssλ     ::Float64,
                      seλ     ::Float64,
                      stn     ::Float64,
                      λ0_prior::NTuple{2,Float64},
                      α_prior ::NTuple{2,Float64},
                      σλ_prior::NTuple{2,Float64},
                      ϵ_prior ::NTuple{2,Float64},
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
  r   = Array{Float64,2}(undef, nlogs, 7)

  λfs   = Float64[]
  treev = cTct[]           # make Ξ vector
  io    = IOBuffer()       # buffer 
  el    = lastindex(idf)   # number of branches

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\talpha\tsigma_lambda\tepsilon\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αc = αc, σλc = σλc, ϵc = ϵc, mc = mc, ns = ns, ne = ne, ddλ = ddλ, ssλ = ssλ, seλ = seλ, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update drift
            if pupi === 1

              llc, prc, αc, mc = 
                updatect_α!(αc, lλ(Ξ[1]), σλc, ϵc, 2.0*(ns + rmλ), ddλ, llc, prc, 
                  mc, th, surv, α_prior)

              # update ssλ with new drift `α`
              ssλ = _ss(Ξ, idf, αc)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, ϵc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update diffusion rate
            elseif pupi === 2

              llc, prc, σλc, mc = 
                updatect_σ!(σλc, lλ(Ξ[1]), αc, ϵc, ssλ, 2.0*(ns + rmλ), llc, prc, 
                  mc, th, surv, σλ_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, ϵc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update extinction
            elseif pupi === 3

              llc, prc, ϵc, mc =
                update_ϵ!(ϵc, lλ(Ξ[1]), αc, σλc, llc, prc, ne, seλ, mc, th, 
                  surv, ϵ_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, ϵc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update scale
            elseif pupi === 4

              llc, prc, mc, seλ, acc = 
                update_scale!(Ξ, idf, αc, σλc, ϵc, llc, prc, ns, ne, stn, 
                  mc, th, surv, λ0_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, ϵc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update internal λ
            elseif pupi === 5

              bix = fIrand(el) + 1

              llc, prc, ddλ, ssλ, seλ, mc =
                update_internal!(bix, Ξ, idf, αc, σλc, ϵc, llc, prc, ddλ, ssλ, 
                  seλ, mc, th, λ0_prior, surv)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, ϵc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update by forward simulation
            else

              bix = fIrand(el) + 1

              llc, ddλ, ssλ, seλ, ns, ne =
                update_fs!(bix, Ξ, idf, αc, σλc, ϵc, llc, ddλ, ssλ, seλ, ns, 
                  ne, λfs)

              # ll0 = llik_clads(Ξ, idf, αc, σλc, ϵc) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
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
              r[lit,7] = ϵc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
                  exp(lλ(Ξ[1])),'\t', αc, '\t', σλc, '\t', ϵc, '\n')
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
    updatect_α!(αc     ::Float64,
              λ0     ::Float64,
              σλ     ::Float64,
              ϵ      ::Float64,
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
function updatect_α!(αc     ::Float64,
                     λ0     ::Float64,
                     σλ     ::Float64,
                     ϵ      ::Float64,
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

  mp  = m_surv_cladsct(th, λ0, αp, σλ, ϵ, 1_000, surv)
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
    updatect_σ!(σλc     ::Float64,
              λ0      ::Float64,
              α       ::Float64,
              ϵ       ::Float64,
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
function updatect_σ!(σλc     ::Float64,
                     λ0      ::Float64,
                     α       ::Float64,
                     ϵ       ::Float64,
                     ssλ     ::Float64,
                     n       ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     mc      ::Float64,
                     th      ::Float64,
                     surv    ::Int64,
                     σλ_prior::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)
  σλp  = sqrt(σλp2)

  mp  = m_surv_cladsct(th, λ0, α, σλp, ϵ, 1_000, surv)
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
    update_ϵ!(ϵc     ::Float64,
              λ0     ::Float64,
              α      ::Float64,
              σλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              ne     ::Float64,
              seλ    ::Float64,
              mc     ::Float64,
              th     ::Float64,
              surv   ::Int64,
              ϵ_prior::NTuple{2,Float64})

Gibbs-MH update for `ϵ`.
"""
function update_ϵ!(ϵc     ::Float64,
                   λ0     ::Float64,
                   α      ::Float64,
                   σλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   ne     ::Float64,
                   seλ    ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   surv   ::Int64,
                   ϵ_prior::NTuple{2,Float64})

  ϵp  = randgamma(ϵ_prior[1] + ne, ϵ_prior[2] + seλ)

  mp  = m_surv_cladsct(th, λ0, α, σλ, ϵp, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += ne * log(ϵp/ϵc) + seλ * (ϵc - ϵp) + llr
    prc += llrdgamma(ϵp, ϵc, ϵ_prior[1], ϵ_prior[2])
    ϵc   = ϵp
    mc   = mp
  end

  return llc, prc, ϵc, mc
end




"""
    update_scale!(Ξ       ::Vector{cTct},
                  idf     ::Vector{iBffs},
                  α       ::Float64,
                  σλ      ::Float64,
                  ϵ       ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ns      ::Float64,
                  ne      ::Float64,
                  stn     ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  surv    ::Int64,
                  λ0_prior::NTuple{2,Float64})

Update scale for speciation.
"""
function update_scale!(Ξ       ::Vector{cTct},
                       idf     ::Vector{iBffs},
                       α       ::Float64,
                       σλ      ::Float64,
                       ϵ       ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ns      ::Float64,
                       ne      ::Float64,
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
  mp  = m_surv_cladsct(th, lλ0 + s, α, σλ, ϵ, 1_000, surv)

  # likelihood ratio
  ir  = _ir(Ξ)
  llr = ns * s + ne * s + (1.0 - exp(s)) * (1.0 + ϵ) * ir + log(mp/mc)

  acc = 0.0
  if -randexp() < llr + prr
    acc  += 1.0
    llc  += llr
    prc  += prr
    mc   = mp
    ir  *= exp(s)
    scale_rate!(Ξ,   addlλ!, s)
    scale_rate!(idf, addlλ!, s)
  end

  return llc, prc, mc, ir, acc
end




"""
    update_internal!(bix     ::Int64,
                     Ξ       ::Vector{cTct},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     ϵ       ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ddλ     ::Float64,
                     ssλ     ::Float64,
                     seλ     ::Float64,
                     mc      ::Float64,
                     th      ::Float64,
                     λ0_prior::NTuple{2,Float64},
                     surv    ::Int64)

Make an update for an internal branch and its descendants.
"""
function update_internal!(bix     ::Int64,
                          Ξ       ::Vector{cTct},
                          idf     ::Vector{iBffs},
                          α       ::Float64,
                          σλ      ::Float64,
                          ϵ       ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          ddλ     ::Float64,
                          ssλ     ::Float64,
                          seλ     ::Float64,
                          mc      ::Float64,
                          th      ::Float64,
                          λ0_prior::NTuple{2,Float64},
                          surv    ::Int64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  it   = iszero(i1) # is terminal
  i2   = d2(bi)
  ia   = pa(bi)
  root = iszero(ia)
  λa   = NaN  # ancestral speciation

  # if crown root
  if root && iszero(e(ξi))
    llc, prc, ddλ, ssλ, mc =
      _crown_update!(ξi, Ξ[i1], Ξ[i2], α, σλ, ϵ, llc, prc, ddλ, ssλ, mc, th, 
        λ0_prior, surv)
    λa = lλ(ξi)
    setλt!(bi, λa)
  else
    # if stem
    if root
      
      eds, λ1, λ2 = 0.0, NaN, NaN
      # if cladogenetic branch
      if i2 > 0
        eds, λ1, λ2 = 0.0, lλ(Ξ[i1]), lλ(Ξ[i2])
      # if mid branch
      else
        eds, λ1, λ2 = downstreamλs(bix, Ξ, idf, 0.0, NaN, NaN)
      end

      llc, prc, ddλ, ssλ, seλ, mc, λi = 
        _stem_update!(ξi, eds, λ1, λ2, 
          α, σλ, ϵ, llc, prc, ddλ, ssλ, seλ, mc, th, λ0_prior, surv)

      # set new λ downstream, if necessary
      setdownstreamλ!(λi, bix, Ξ, idf)

      # if there are speciation events in stem branch
      if !istip(ξi)
        eds, λ1, λ2 = downstreamλs(i1, Ξ, idf, 0.0, NaN, NaN)

        # updates within the parent branch
        llc, ddλ, ssλ, seλ, λx = 
          _update_internal!(ξi.d1, bi, eas, λi, α, σλ, ϵ, eds, λ1, λ2, 
            llc, ddλ, ssλ, seλ, false)
        llc, ddλ, ssλ, seλ, λx = 
          _update_internal!(ξi.d2, bi, eas, λi, α, σλ, ϵ, eds, λ1, λ2, 
            llc, ddλ, ssλ, seλ, false)

        setdownstreamλ!(λi, i1, Ξ, idf)
      end

    # if *not* root
    else

      # find cladogenetic ancestor
      eas, λa, il = upstreamλ(ia, Ξ, idf, 0.0, λa)

      # check if mid branch does not lead to the stem branch
      if pa(idf[il]) > 0

        # it non-terminal branch
        eds, λ1, λ2 = 0.0, NaN, NaN
        if !it
          # if cladogenetic branch
          if i2 > 0
            eds, λ1, λ2 = 0.0, lλ(Ξ[i1]), lλ(Ξ[i2])
          # if mid branch
          else
            eds, λ1, λ2 = downstreamλs(i1, Ξ, idf, 0.0, NaN, NaN)
          end
        end

        ll0 = llc

        # updates within the parent branch
        llc, ddλ, ssλ, seλ, λx = 
          _update_internal!(ξi, bi, eas, λa, α, σλ, ϵ, eds, λ1, λ2, llc, 
            ddλ, ssλ, seλ, it)

        # if update, update up- and down-stream
        if ll0 != llc
          λi = lλ(ξi)
          setupstreamλ!(λi, ia, Ξ, idf)
          λi = lλ(fixtip(ξi))
          (!it && iszero(i2)) && setdownstreamλ!(λi, i1, Ξ, idf)
        end
      end
    end
  end

  return llc, prc, ddλ, ssλ, seλ, mc
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{cTct},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               ϵ    ::Float64,
               llc  ::Float64,
               ddλ  ::Float64,
               ssλ  ::Float64,
               seλ  ::Float64,
               ns   ::Float64,
               ne   ::Float64,
               λfs  ::Vector{Float64})

Forward simulation proposal function for pure birth diffusion.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{cTct},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    ϵ    ::Float64,
                    llc  ::Float64,
                    ddλ  ::Float64,
                    ssλ  ::Float64,
                    seλ  ::Float64,
                    ns   ::Float64,
                    ne   ::Float64,
                    λfs  ::Vector{Float64})

  bi = idf[bix]
  ξc = Ξ[bix]
  ia = pa(bi)

  λa = NaN
  # if following a speciation event
  if ia > 0 && d2(idf[ia]) > 0
    λa = λt(idf[ia])
  end

  ddλr = ssλr = seλr = zero(Float64)
  llr  = NaN
  # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, ξc, λa, α, σλ, ϵ)

  # if mid
  elseif iszero(d2(bi))
    ξp, llr, ddλr, ssλr, seλr = 
      fsbi_m(bi, idf, ξc, Ξ, λa, α, σλ, ϵ, λfs)

  # if internal
  else
    if e(bi) > 0.0
      ξp, llr, ddλr, ssλr = 
        fsbi_i(bi, ξc, λa, lλ(Ξ[d1(bi)]), lλ(Ξ[d2(bi)]), α, σλ, ϵ, λfs)
    end
  end

  # if accepted
  if isfinite(llr)

    llc, ddλ, ssλ, seλ, ns, ne = 
      llik_cladsct_track!(ξc, α, σλ, ϵ, llc, ddλ, ssλ, seλ, ns, ne, -)
    llc, ddλ, ssλ, seλ, ns, ne = 
      llik_cladsct_track!(ξp, α, σλ, ϵ, llc, ddλ, ssλ, seλ, ns, ne, +)

    # first change from ancestor
    if !isnan(λa)
      λp, λc = lλ(ξp), lλ(ξc)
      llc += llrdnorm_x(λp, λc, λa + α, σλ^2)
      ddλ += λp - λc
      ssλ += 0.5*((λp - λa - α)^2 - (λc - λa - α)^2)
    end

    # update quantities
    ddλ += ddλr
    ssλ += ssλr
    seλ += seλr
    llc += llr

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ssλ, seλ, ns, ne
end




"""
    fsbi_t(bi::iBffs,
           ξi ::cTct,
           λa::Float64,
           α ::Float64,
           σλ::Float64,
           ϵ ::Float64)

Forward simulation for terminal branch `bi`.
"""
function fsbi_t(bi::iBffs,
                ξi ::cTct,
                λa::Float64,
                α ::Float64,
                σλ::Float64,
                ϵ ::Float64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # if does **not** come from a cladogenetic event
  ncl = isnan(λa)
  if ncl
    λi = lλ(ξi)
  else
    λi = rnorm(λa + α, σλ)
  end

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_cladsct_t(e(bi), λi, α, σλ, ϵ, lc, lU, iρi, 0, 1, 500)

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
           ξc ::cTct,
           Ξ  ::Vector{cTct},
           λa ::Float64,
           α  ::Float64,
           σλ ::Float64,
           ϵ  ::Float64,
           λfs::Vector{Float64})

Forward simulation for internal branch `bi`
"""
function fsbi_m(bi ::iBffs,
                idf::Vector{iBffs},
                ξi ::cTct,
                Ξ  ::Vector{cTct},
                λa ::Float64,
                α  ::Float64,
                σλ ::Float64,
                ϵ  ::Float64,
                λfs::Vector{Float64})

  # if does **not** come from a cladogenetic event
  ncl = isnan(λa)

  if ncl
    λi = lλ(ξi)
  else
    λi = rnorm(λa + α, σλ)
  end

  # forward simulation during branch length
  empty!(λfs)
  t0, na, nn = _sim_cladsct_i(e(bi), λi, α, σλ, ϵ, 0, 1, 500, λfs)

  if na < 1 || nn > 499
    return t0, NaN, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # search for next lλ1 and lλ2 if the exist
  i1 = d1(bi)
  eds, λ1, λ2 = downstreamλs(i1, Ξ, idf, 0.0, NaN, NaN)

  ## choose most likely lineage to fix
  # if downstream is tip
  ddr = ssr = ser = 0.0
  if isnan(λ1)
    wt, λp, pp, λc, pc, acr, ser = wfix_m(ξi, e(bi), λfs, eds, ϵ, acr)
  # if downstream is cladogenetic
  else
    wt, λp, pp, λc, pc, acr, ddr, ssr, ser = 
      wfix_m(ξi, e(bi), λfs, eds, λ1, λ2, α, σλ, ϵ, acr)
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
      tip_sims!(t0, tf(bi), α, σλ, ϵ, acr, lU, iρi, na, nn)

    if lU < acr

      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc)
      if isfinite(λ1)
        llr += λp - λc
      end
      setni!(bi, na) # set new ni

      # downstream change
      setdownstreamλ!(λp, i1, Ξ, idf)

      return t0, llr, ddr, ssr, ser
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    wfix_m(ξi ::cTct,
           ei ::Float64,
           λfs::Vector{Float64},
           eds::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for middle `i` branches with downstream **tips**.
"""
function wfix_m(ξi ::cTct,
                ei ::Float64,
                λfs::Vector{Float64},
                eds::Float64,
                ϵ  ::Float64,
                acr::Float64)

  # select best from proposal
  sp, i, wt, λp, pp = 0.0, 0, 0, NaN, -Inf
  for λfi in λfs
    i += 1
    p  = exp(- eds * exp(λfi) * (1.0 + ϵ))
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
    p   = exp(- eds * exp(λfi) * (1.0 + ϵ))
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp/sc)
  ser  = eds * (exp(λp) - exp(λc))

  return wt, λp, pp, λc, pc, acr, ser
end




"""
    wfix_m(ξi ::cTct,
           ei ::Float64,
           λfs::Vector{Float64},
           eds::Float64,
           λ1 ::Float64,
           λ2 ::Float64,
           α  ::Float64,
           σλ ::Float64,
           ϵ  ::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for middle `i` branches with downstream **cladogenetic** daughters.
"""
function wfix_m(ξi ::cTct,
                ei ::Float64,
                λfs::Vector{Float64},
                eds::Float64,
                λ1 ::Float64,
                λ2 ::Float64,
                α  ::Float64,
                σλ ::Float64,
                ϵ  ::Float64,
                acr::Float64)

  # select best from proposal
  sp, i, wt, λp, pp = 0.0, 0, 0, NaN, -Inf
  for λfi in λfs
    i  += 1
    p   = dnorm2(λ1, λ2, λfi + α, σλ) * exp(- eds * exp(λfi) * (1.0 + ϵ))
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
    p   = dnorm2(λ1, λ2, λfi + α, σλ) * exp(- eds * exp(λfi) * (1.0 + ϵ))
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
  ser  = eds * (exp(λp) - exp(λc))

  return wt, λp, pp, λc, pc, acr, ddr, ssr, ser
end




"""
    fsbi_i(bi  ::iBffs,
           ξi  ::cTct,
           λa  ::Float64,
           λ1  ::Float64,
           λ2  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           ϵ   ::Float64,
           λfs ::Vector{Float64})

Forward simulation for internal branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξi  ::cTct,
                λa  ::Float64,
                λ1  ::Float64,
                λ2  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                ϵ   ::Float64,
                λfs ::Vector{Float64})

  # if does **not** come from a cladogenetic event
  ncl = isnan(λa)
  if ncl
    λi = lλ(ξi)
  else
    λi = rnorm(λa + α, σλ)
  end

  empty!(λfs)

  # forward simulation during branch length
  t0, na, nn = _sim_cladsct_i(e(bi), λi, α, σλ, ϵ, 0, 1, 500, λfs)

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
    wfix_i(ξi, e(bi), λfs, λ1, λ2, α, σλ, acr)

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
        tip_sims!(t0, tf(bi), α, σλ, ϵ, acr, lU, iρi, na, nn)
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
    wfix_i(ξi ::cTct,
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
function wfix_i(ξi ::cTct,
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
    tip_sims!(tree::cTct,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              ϵ   ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::cTct,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   ϵ   ::Float64,
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
          _sim_cladsct_it(t, lλ(tree), α, σλ, ϵ, lr, lU, iρi, na-1, nn, 500)

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
        tip_sims!(tree.d1, t, α, σλ, ϵ, lr, lU, iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, α, σλ, ϵ, lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end



