#=

Anagenetic occurrence birth-death MCMC using forward simulation

Jérémy Andréoletti

v(^-^v)

Created 20 09 2023
=#




"""
    insane_gbmobd(tree    ::sTf_label,
                  ωtimes  ::Vector{Float64};
                  λa_prior::NTuple{2,Float64} = (1.5, 0.5),
                  μa_prior::NTuple{2,Float64} = (1.5, 1.0),
                  α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                  σλ_prior::NTuple{2,Float64} = (3.0, 0.5),
                  σμ_prior::NTuple{2,Float64} = (3.0, 0.5),
                  ψ_prior ::NTuple{2,Float64} = (1.0, 1.0),
                  ω_prior ::NTuple{2,Float64} = (1.0, 1.0),
                  ψω_epoch::Vector{Float64}   = Float64[],
                  f_epoch ::Vector{Int64}     = Int64[],
                  niter   ::Int64             = 1_000,
                  nthin   ::Int64             = 10,
                  nburn   ::Int64             = 200,
                  nflush  ::Int64             = nthin,
                  ofile   ::String            = homedir(),
                  tune_int::Int64             = 100,
                  ϵi      ::Float64           = 0.2,
                  λi      ::Float64           = NaN,
                  μi      ::Float64           = NaN,
                  ψi      ::Float64           = NaN,
                  ωi      ::Float64           = NaN,
                  ωtni    ::Float64           = 1.0,
                  obj_ar  ::Float64           = 0.4,
                  αi      ::Float64           = 0.0,
                  σλi     ::Float64           = 0.01,
                  σμi     ::Float64           = 0.01,
                  pupdp   ::NTuple{6,Float64} = (0.01, 0.01, 0.01, 0.01, 0.1, 0.2),
                  δt      ::Float64           = 1e-3,
                  survival::Bool              = true,
                  mxthf   ::Float64           = Inf,
                  prints  ::Int64             = 5,
                  tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for occurrence birth-death diffusion `fbdd`.
"""
function insane_gbmobd(tree    ::sTf_label,
                       ωtimes  ::Vector{Float64};
                       λa_prior::NTuple{2,Float64} = (1.5, 0.5),
                       μa_prior::NTuple{2,Float64} = (1.5, 1.0),
                       α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                       σλ_prior::NTuple{2,Float64} = (3.0, 0.5),
                       σμ_prior::NTuple{2,Float64} = (3.0, 0.5),
                       ψ_prior ::NTuple{2,Float64} = (1.0, 1.0),
                       ω_prior ::NTuple{2,Float64} = (1.0, 0.2),
                       ψω_epoch::Vector{Float64}   = Float64[],
                       f_epoch ::Vector{Int64}     = Int64[0],
                       niter   ::Int64             = 1_000,
                       nthin   ::Int64             = 10,
                       nburn   ::Int64             = 200,
                       nflush  ::Int64             = nthin,
                       ofile   ::String            = homedir(),
                       tune_int::Int64             = 100,
                       ϵi      ::Float64           = 0.2,
                       λi      ::Float64           = NaN,
                       μi      ::Float64           = NaN,
                       ψi      ::Float64           = NaN,
                       ωi      ::Float64           = NaN,
                       ωtni    ::Float64           = 1.0,
                       obj_ar  ::Float64           = 0.4,
                       αi      ::Float64           = 0.0,
                       σλi     ::Float64           = 0.1,
                       σμi     ::Float64           = 0.1,
                       pupdp   ::NTuple{6,Float64} = (0.01, 0.01, 0.01, 0.01, 0.1, 0.2),
                       δt      ::Float64           = 1e-3,
                       survival::Bool              = true,
                       mxthf   ::Float64           = Inf,
                       prints  ::Int64             = 5,
                       tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1,round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)
  LTT  = ltt(tree)

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # only include epochs where the tree occurs
  sort!(ψω_epoch, rev = true)
  tix = findfirst(x -> x < th, ψω_epoch)
  if !isnothing(tix)
    ψω_epoch = ψω_epoch[tix:end]
  end
  nep  = lastindex(ψω_epoch) + 1

  # make initial fossils per epoch vector
  if lastindex(f_epoch) !== nep
    f_epoch = fill(0, nep)
  end

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # estimate branch split (multiple of δt)
  ndts = floor(th * mxthf/δt)
  maxt = δt * ndts

  # make fix tree directory
  idf = make_idf(tree, tρ, maxt)

  # starting parameters
  if isnan(λi) && isnan(μi) && isnan(ψi)
    # if only one tip
    if isone(n)
      λc = prod(λa_prior)
      μc = prod(μa_prior)
    else
      λc, μc = moments(Float64(n), th, ϵi)
    end
    # if no sampled fossil
    nf = nfossils(tree)
    if iszero(nf)
      ψc = prod(ψ_prior)
    else
      ψc = Float64(nf)/Float64(treelength(tree))
    end
    
    nω = lastindex(ωtimes)
    # if no fossil occurrences
    if iszero(nω)
      ωc = ω_prior[1]/ω_prior[2]
    else
      ωc = Float64(nω)/treelength(tree)
    end
  else
    λc, μc, ψc, ωc = λi, μi, ψi, ωi
  end

  # make ψ and ω vectors for each epoch
  ψc = fill(ψc, nep)
  ωc = fill(ωc, nep)

  # count the number of occurrences in each epoch
  sort!(ωtimes,   rev=true)
  sort!(ψω_epoch, rev=true)

  if !isempty(ψω_epoch)
    nω = zeros(Int, nep)
    ep = 1

    for t in ωtimes
      while ep < nep && t < ψω_epoch[ep]
        ep += 1
      end
      nω[ep] += 1
    end
  else
    nω = [nω]
  end

  # condition on survival of 0, 1, or 2 starting lineages
  surv = 0
  if survival 
    if iszero(e(tree)) 
      if def1(tree)
        surv += (ntipsalive(tree.d1) > 0)
        if def2(tree)
          surv += (ntipsalive(tree.d2) > 0)
        end
      end
    else
      surv += (ntipsalive(tree) > 0)
    end
  end

  # M attempts of survival
  mc = m_surv_gbmbd(th, log(λc), log(μc), αi, σλi, σμi, δt, srδt, 5_000, surv)

  # make a decoupled tree
  Ξ = make_Ξ(idf, λc, μc, αi, σλi, σμi, δt, srδt, iTfbd)

  # set end of fix branch speciation times and get vector of internal branches
  # and make epoch start vectors and indices for each `ξ`
  inodes = Int64[]
  eixi   = Int64[]
  eixf   = Int64[]
  bst    = Float64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    if d1(bi) > 0
      push!(inodes, i)
    end
    tib = ti(bi)
    ei  = findfirst(x -> x < tib, ψω_epoch)
    ei  = isnothing(ei) ? nep : ei
    ef  = findfirst(x -> x < tf(bi), ψω_epoch)
    ef  = isnothing(ef) ? nep : ef
    push!(bst, tib)
    push!(eixi, ei)
    push!(eixf, ef)
  end

  # parameter updates (1: α, 2: σλ & σμ, 3: ψ, 4: ω, 5: gbm, 6: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(6)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running occurrence birth-death diffusion"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, σμc, ψc, ωc, ωtn, mc, L, LTT  =
    mcmc_burn_gbmobd(Ξ, idf, ωtimes, LTT, tune_int, λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior,
      ψ_prior, ω_prior, ψω_epoch, f_epoch, nburn, αi, σλi, σμi, ψc, ωc, ωtni, scalef, mc, nω, 
      th, surv, δt, srδt, bst, eixi, eixf, inodes, pup, prints)

  # mcmc
  r, treev =
    mcmc_gbmobd(Ξ, idf, ωtimes, LTT, llc, prc, αc, σλc, σμc, ψc, ωc, ωtn, mc, nω, L, th, surv,
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, 
      ψ_prior, ω_prior, ψω_epoch, f_epoch, δt, srδt, bst, eixi, eixf, inodes, pup, 
      niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_gbmobd(Ξ       ::Vector{iTfbd},
                     idf     ::Vector{iBffs},
                     ωtimes  ::Vector{Float64},
                     LTT     ::Ltt,
                     tune_int::Int64,
                     λa_prior::NTuple{2,Float64},
                     μa_prior::NTuple{2,Float64},
                     α_prior ::NTuple{2,Float64},
                     σλ_prior::NTuple{2,Float64},
                     σμ_prior::NTuple{2,Float64},
                     ψ_prior ::NTuple{2,Float64},
                     ω_prior ::NTuple{2,Float64},
                     ψω_epoch::Vector{Float64},
                     f_epoch ::Vector{Int64},
                     nburn   ::Int64,
                     αc      ::Float64,
                     σλc     ::Float64,
                     σμc     ::Float64,
                     ψc      ::Vector{Float64},
                     ωc      ::Vector{Float64},
                     ωtni    ::Float64,
                     scalef  ::Function,
                     mc      ::Float64,
                     nω      ::Vector{Int64},
                     th      ::Float64,
                     surv    ::Int64,
                     δt      ::Float64,
                     srδt    ::Float64,
                     bst     ::Vector{Float64},
                     eixi    ::Vector{Int64},
                     eixf    ::Vector{Int64},
                     inodes  ::Array{Int64,1},
                     pup     ::Array{Int64,1},
                     prints  ::Int64)


MCMC burn-in chain for `fbdd`.
"""
function mcmc_burn_gbmobd(Ξ       ::Vector{iTfbd},
                          idf     ::Vector{iBffs},
                          ωtimes  ::Vector{Float64},
                          LTT     ::Ltt,
                          tune_int::Int64,
                          λa_prior::NTuple{2,Float64},
                          μa_prior::NTuple{2,Float64},
                          α_prior ::NTuple{2,Float64},
                          σλ_prior::NTuple{2,Float64},
                          σμ_prior::NTuple{2,Float64},
                          ψ_prior ::NTuple{2,Float64},
                          ω_prior ::NTuple{2,Float64},
                          ψω_epoch::Vector{Float64},
                          f_epoch ::Vector{Int64},
                          nburn   ::Int64,
                          αc      ::Float64,
                          σλc     ::Float64,
                          σμc     ::Float64,
                          ψc      ::Vector{Float64},
                          ωc      ::Vector{Float64},
                          ωtni    ::Float64,
                          scalef  ::Function,
                          mc      ::Float64,
                          nω      ::Vector{Int64},
                          th      ::Float64,
                          surv    ::Int64,
                          δt      ::Float64,
                          srδt    ::Float64,
                          bst     ::Vector{Float64},
                          eixi    ::Vector{Int64},
                          eixf    ::Vector{Int64},
                          inodes  ::Array{Int64,1},
                          pup     ::Array{Int64,1},
                          prints  ::Int64)

  L   = treelength(Ξ, ψω_epoch, bst, eixi) # tree length

  # initialize acceptance log
  ltn = 0
  lup = 0.0
  lac = 0.0
  ωtn = ωtni

  # likelihood
  lλ0  = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, ψω_epoch, bst, eixi, δt, srδt) + 
        ω_llik(ωtimes, ωc, ψω_epoch, LTT) -
        (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ0 + log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2,        σλ_prior[1], σλ_prior[2])  +
        logdinvgamma(σμc^2,        σμ_prior[1], σμ_prior[2])  +
        logdnorm(αc,               α_prior[1],  α_prior[2]^2) +
        logdgamma(exp(lλ0),          λa_prior[1], λa_prior[2]) +
        logdgamma(exp(lμ(Ξ[1])[1]), μa_prior[1], μa_prior[2]) +
        sum(logdgamma.(ψc, ψ_prior[1], ψ_prior[2]))  +
        sum(logdgamma.(ωc, ω_prior[1], ω_prior[2]))

  L            = treelength(Ξ, ψω_epoch, bst, eixi) # tree length
  nf           = nfossils(idf, ψω_epoch, f_epoch)   # number of fossilization events per epoch
  dλ           = deltaλ(Ξ)                         # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ξ, αc)                    # sum squares in λ and μ
  nin          = lastindex(inodes)                 # number of internal nodes
  el           = lastindex(idf)                    # number of branches
  nep          = lastindex(ψc)

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup
      # print("\n", ["α","σλ & σμ","ψ","ω","gbm update","forward simulation"][pupi], "\n")

      # update α
      if pupi === 1

        llc, prc, αc, mc  =
          update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, sum(L), dλ, llc, prc,
            mc, th, surv, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nλ = sss_gbm(Ξ, αc)

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nλ,
            llc, prc, mc, th, surv, δt, srδt, σλ_prior, σμ_prior)

      # psi update
      elseif pupi === 3

        llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # ω proposal
      elseif pupi === 4

        llc, prc, lac, lup = update_ω!(llc, prc, ωc, ψω_epoch, ωtimes, LTT, nω, L, lac, lup, ωtn, ω_prior)

      # gbm update
      elseif pupi === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, prc, dλ, ssλ, ssμ, mc, th,
            surv, δt, srδt, λa_prior, μa_prior)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, nλ, L, LTT =
          update_fs!(bix, Ξ, idf, ωtimes, LTT, αc, σλc, σμc, ψc, ωc, llc, dλ, ssλ, ssμ, nλ, L,
            ψω_epoch, δt, srδt, eixi, eixf)

      end

      # log tuning parameters
      ltn += 1
      if ltn == tune_int
        if lup > 0.0
          ωtn = scalef(ωtn,lac/lup)
        end
        ltn = 0
      end
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, σμc, ψc, ωc, ωtn, mc, L, LTT
end




"""
    mcmc_gbmobd(Ξ       ::Vector{iTfbd},
                idf     ::Vector{iBffs},
                ωtimes  ::Vector{Float64},
                LTT     ::Ltt,
                llc     ::Float64,
                prc     ::Float64,
                αc      ::Float64,
                σλc     ::Float64,
                σμc     ::Float64,
                ψc      ::Vector{Float64},
                ωc      ::Vector{Float64},
                ωtn     ::Float64,
                mc      ::Float64,
                nω      ::Vector{Int64},
                L       ::Vector{Float64},
                th      ::Float64,
                surv    ::Int64,
                λa_prior::NTuple{2,Float64},
                μa_prior::NTuple{2,Float64},
                α_prior ::NTuple{2,Float64},
                σλ_prior::NTuple{2,Float64},
                σμ_prior::NTuple{2,Float64},
                ψ_prior ::NTuple{2,Float64},
                ω_prior ::NTuple{2,Float64},
                ψω_epoch::Vector{Float64},
                f_epoch ::Vector{Int64},
                δt      ::Float64,
                srδt    ::Float64,
                bst     ::Vector{Float64},
                eixi    ::Vector{Int64},
                eixf    ::Vector{Int64},
                inodes  ::Array{Int64,1},
                pup     ::Vector{Int64},
                niter   ::Int64,
                nthin   ::Int64,
                nflush  ::Int64,
                ofile   ::String,
                prints  ::Int64)

MCMC chain for `gbmbd`.
"""
function mcmc_gbmobd(Ξ       ::Vector{iTfbd},
                     idf     ::Vector{iBffs},
                     ωtimes  ::Vector{Float64},
                     LTT     ::Ltt,
                     llc     ::Float64,
                     prc     ::Float64,
                     αc      ::Float64,
                     σλc     ::Float64,
                     σμc     ::Float64,
                     ψc      ::Vector{Float64},
                     ωc      ::Vector{Float64},
                     ωtn     ::Float64,
                     mc      ::Float64,
                     nω      ::Vector{Int64},
                     L       ::Vector{Float64},
                     th      ::Float64,
                     surv    ::Int64,
                     λa_prior::NTuple{2,Float64},
                     μa_prior::NTuple{2,Float64},
                     α_prior ::NTuple{2,Float64},
                     σλ_prior::NTuple{2,Float64},
                     σμ_prior::NTuple{2,Float64},
                     ψ_prior ::NTuple{2,Float64},
                     ω_prior ::NTuple{2,Float64},
                     ψω_epoch::Vector{Float64},
                     f_epoch ::Vector{Int64},
                     δt      ::Float64,
                     srδt    ::Float64,
                     bst     ::Vector{Float64},
                     eixi    ::Vector{Int64},
                     eixf    ::Vector{Int64},
                     inodes  ::Array{Int64,1},
                     pup     ::Vector{Int64},
                     niter   ::Int64,
                     nthin   ::Int64,
                     nflush  ::Int64,
                     ofile   ::String,
                     prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0
  
  nf           = nfossils(idf, ψω_epoch, f_epoch)   # number of fossilization events per epoch
  dλ           = deltaλ(Ξ)                         # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ξ, αc)                    # sum squares in λ and μ
  nin          = lastindex(inodes)                 # number of internal nodes
  el           = lastindex(idf)                    # number of branches
  nep          = lastindex(ψc)

  # parameter results
  nep = lastindex(ψc)
  r = Array{Float64,2}(undef, nlogs, 8 + 2*nep)

  # make Ξ vector
  treev = iTfbd[]

  # number of branches and of triads
  nbr  = lastindex(idf)

  # flush to file
  sthin = 0

  function check_pr(pupi::Int64, i::Int64)
    pr0 = logdinvgamma(σλc^2,        σλ_prior[1], σλ_prior[2])  +
          logdinvgamma(σμc^2,        σμ_prior[1], σμ_prior[2])  +
          logdnorm(αc,               α_prior[1],  α_prior[2]^2) +
          logdgamma(exp(lλ(Ξ[1])[1]),          λa_prior[1], λa_prior[2]) +
          logdgamma(exp(lμ(Ξ[1])[1]), μa_prior[1], μa_prior[2]) +
          sum(logdgamma.(ψc, ψ_prior[1], ψ_prior[2]))  +
          sum(logdgamma.(ωc, ω_prior[1], ω_prior[2]))
    if !isapprox(pr0, prc, atol = 1e-5)
       error(string("Wrong prior computation during the ", ["α","σλ & σμ","ψ","ω","gbm update","forward simulation"][pupi], 
                    " update, at iteration ", i, ": pr0=", pr0, " and prc-pr0=", prc-pr0))
    end
  end

  function check_ll(pupi::Int64, i::Int64)
    ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, ψω_epoch, bst, eixi, δt, srδt) +
          ω_llik(ωtimes, ωc, ψω_epoch, LTT) - 
          (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
    if !isapprox(ll0, llc, atol = 1e-5)
       error(string("Wrong likelihood computation during the ", ["α","σλ & σμ","ψ","ω","gbm update","forward simulation"][pupi], 
                    " update, at iteration ", i, ": ll0=", ll0, " and llc-ll0=", llc-ll0))
    end
  end

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\tmu_root\talpha\tsigma_lambda\tsigma_mu\t"*join(["psi"*(isone(nep) ? "" : string("_",i)) for i in 1:nep], "\t")*"\t"*join(["omega"*(isone(nep) ? "" : string("_",i)) for i in 1:nep], "\t")*"\n")
    flush(of)

    open(ofile*".txt", "w") do tf


      pbar = Progress(niter, prints, "running mcmc...", 20)

      for it in Base.OneTo(niter)

        shuffle!(pup)

        # parameter updates
        for pupi in pup

          # update α
          if pupi === 1

            llc, prc, αc, mc  =
              update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, sum(L), dλ, llc, prc,
                mc, th, surv, δt, srδt, α_prior)

            # update ssλ with new drift `α`
            ssλ, ssμ, nλ = sss_gbm(Ξ, αc)

          # σλ & σμ update
          elseif pupi === 2

            llc, prc, σλc, σμc, mc =
              update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nλ,
                llc, prc, mc, th, surv, δt, srδt, σλ_prior, σμ_prior)

          # psi update
          elseif pupi === 3

            llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

          # ω proposal
          elseif pupi === 4

            llc, prc = update_ω!(llc, prc, ωc, ψω_epoch, ωtimes, LTT, nω, L, ωtn, ω_prior)

          # gbm update
          elseif pupi === 5

            nix = ceil(Int64,rand()*nin)
            bix = inodes[nix]

            llc, prc, dλ, ssλ, ssμ, mc =
              update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, prc, dλ, ssλ, ssμ, mc, th,
                surv, δt, srδt, λa_prior, μa_prior)

          # forward simulation update
          else

            bix = ceil(Int64,rand()*el)

            llc, dλ, ssλ, ssμ, nλ, L, LTT =
              update_fs!(bix, Ξ, idf, ωtimes, LTT, αc, σλc, σμc, ψc, ωc, llc, dλ, ssλ, 
                ssμ, nλ, L, ψω_epoch, δt, srδt, eixi, eixf)

          end

          # check_pr(pupi, it)
          # check_ll(pupi, it)
        end

        # log parameters
        lthin += 1
        if lthin === nthin
          lit += 1
          @inbounds begin
            r[lit,1] = Float64(it)
            r[lit,2] = llc
            r[lit,3] = prc
            r[lit,4] = exp(lλ(Ξ[1])[1])
            r[lit,5] = exp(lμ(Ξ[1])[1])
            r[lit,6] = αc
            r[lit,7] = σλc
            r[lit,8] = σμc
            @avx for i in Base.OneTo(nep)
              r[lit,8 + i] = ψc[i]
            end
            @avx for i in Base.OneTo(nep)
              r[lit,8 + nep + i] = ωc[i]
            end
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthin += 1
        if sthin === nflush
          write(of, 
            string(Float64(it), "\t", llc, "\t", prc, "\t", 
              exp(lλ(Ξ[1])[1]),"\t", exp(lμ(Ξ[1])[1]), "\t", αc, "\t",
               σλc, "\t", σμc, "\t", join(ψc, "\t"), "\t", join(ωc, "\t"), "\n"))
          flush(of)
          write(tf, 
            string(istring(couple(Ξ, idf, 1)), "\n"))
          flush(tf)
          sthin = 0
        end
        next!(pbar)
      end
    end
  end

  return r, treev
end




"""
    uupdate_fs!(bix   ::Int64,
                Ξ     ::Vector{iTfbd},
                idf   ::Vector{iBffs},
                ωtimes::Vector{Float64},
                LTT   ::Ltt,
                α     ::Float64,
                σλ    ::Float64,
                σμ    ::Float64,
                ψ     ::Vector{Float64},
                ω     ::Vector{Float64},
                llc   ::Float64,
                dλ    ::Float64,
                ssλ   ::Float64,
                ssμ   ::Float64,
                nλ    ::Float64,
                L     ::Vector{Float64},
                ψωts  ::Vector{Float64},
                δt    ::Float64,
                srδt  ::Float64,
                eixi  ::Vector{Int64},
                eixf  ::Vector{Int64})

Forward simulation proposal function for `gbmobd`.
"""
function update_fs!(bix   ::Int64,
                    Ξ     ::Vector{iTfbd},
                    idf   ::Vector{iBffs},
                    ωtimes::Vector{Float64},
                    LTT   ::Ltt,
                    α     ::Float64,
                    σλ    ::Float64,
                    σμ    ::Float64,
                    ψ     ::Vector{Float64},
                    ω     ::Vector{Float64},
                    llc   ::Float64,
                    dλ    ::Float64,
                    ssλ   ::Float64,
                    ssμ   ::Float64,
                    nλ    ::Float64,
                    L     ::Vector{Float64},
                    ψωts  ::Vector{Float64},
                    δt    ::Float64,
                    srδt  ::Float64,
                    eixi  ::Vector{Int64},
                    eixf  ::Vector{Int64})

  bi  = idf[bix]
  ξc  = Ξ[bix]
  ixi = eixi[bix]

  # terminal branch
  if iszero(d1(bi))

    # fossil terminal branch
    if isfossil(bi)
      ixf = eixf[bix]
      # print("\n") ; @show "fsbi_f"

      ξp, LTTp, llr = fsbi_f(bi, ξc, α, σλ, σμ, ψ, ω, ψωts, ωtimes, LTT, ixi, ixf, δt, srδt)

      # if terminal but not successful proposal, update extinct
      if !isfinite(llr)
        # print("\n") ; @show "fsbi_et"
        ξp, LTTp, llr = fsbi_et(bi, ξc, iTfbd_wofe(ξc), α, σλ, σμ, ψ, ω, ψωts, ωtimes, LTT, ixf,
          δt, srδt)
      end

    # non-fossil terminal branch
    else
      # print("\n") ; @show "fsbi_t"
      ξp, LTTp, llr = fsbi_t(bi, ξc, α, σλ, σμ, ψ, ω, ψωts, ωtimes, LTT, ixi, δt, srδt)
    end
    drλ = ssrλ = ssrμ = 0.0

  # internal non-bifurcating branch
  elseif iszero(d2(bi))
    # print("\n") ; @show "fsbi_m"

    ξp, LTTp, llr, drλ, ssrλ, ssrμ =
      fsbi_m(bi, ξc, Ξ[d1(bi)], α, σλ, σμ, ψ, ω, ψωts, ωtimes, LTT, ixi, eixf[bix], δt, srδt)

  # internal bifurcating branch
  else
    # print("\n") ; @show "fsbi_i"

    ξp, LTTp, llr, drλ, ssrλ, ssrμ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, σμ, ψ, ω, ψωts, ωtimes, LTT, 
        ixi, eixf[bix], δt, srδt)
  end

  # @show ξc, ξp
  # @show round.((ti(bi) - ltt(ξp).t[1]) .+ ltt(ξp).t, digits=2)
  # @show (ti(bi) - ltt(ξp).t[1]) .+ ltt(ξp).t
  # @show idf[bix]

  if isfinite(llr)
    tii = ti(bi)

    # function show_rec(ξ)
    #   @show e(ξ)
    #   def1(ξ) && show_rec(ξ.d1)
    #   def2(ξ) && show_rec(ξ.d2)
    # end

    nep = lastindex(ψωts) + 1

    ll1, ixd, dλ1, ssλ1, ssμ1, nλ1 =
      llik_gbm_ss(ξp, α, σλ, σμ, ψ, tii, ψωts, ixi, δt, srδt, nep)
    ll0, ixd, dλ0, ssλ0, ssμ0, nλ0 =
      llik_gbm_ss(ξc, α, σλ, σμ, ψ, tii, ψωts, ixi, δt, srδt, nep)

    # update llr, ssλ, ssμ, nλ, sns
    llc += ll1  - ll0  + llr
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    ssμ += ssμ1 - ssμ0 + ssrμ
    nλ  += nλ1  - nλ0
    # update tree lengths
    Lc = zeros(Float64, nep)
    _treelength!(ξc, tii, Lc, ψωts, ixi, nep)
    _treelength!(ξp, tii, L,  ψωts, ixi, nep)
    @avx for i in Base.OneTo(nep)
      L[i] -= Lc[i]
    end

    # set new decoupled tree
    Ξ[bix] = ξp

    # set new LTT
    LTT = LTTp
 
  end
  # @assert ltt_rm_artefacts(couple(Ξ, idf, 1)).n == LTTp.n "[[ltt_rm_artefacts(couple(Ξ, idf, 1)).n, ltt_rm_artefacts(couple(Ξ, idf, 1)).t], [LTTp.n, LTTp.t]] = ", [[ltt_rm_artefacts(couple(Ξ, idf, 1)).n, ltt_rm_artefacts(couple(Ξ, idf, 1)).t], [LTTp.n, LTTp.t]], eventimes(couple(Ξ, idf, 1))

  return llc, dλ, ssλ, ssμ, nλ, L, LTT
end




"""
    fsbi_t(bi    ::iBffs,
           ξc    ::iTfbd,
           α     ::Float64,
           σλ    ::Float64,
           σμ    ::Float64,
           ψ     ::Vector{Float64},
           ω     ::Vector{Float64},
           ψωts  ::Vector{Float64},
           ωtimes::Vector{Float64},
           LTT   ::Ltt,
           ix    ::Int64,
           δt    ::Float64,
           srδt  ::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi    ::iBffs,
                ξc    ::iTfbd,
                α     ::Float64,
                σλ    ::Float64,
                σμ    ::Float64,
                ψ     ::Vector{Float64},
                ω     ::Vector{Float64},
                ψωts  ::Vector{Float64},
                ωtimes::Vector{Float64},
                LTT   ::Ltt,
                ix    ::Int64,
                δt    ::Float64,
                srδt  ::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  nep = lastindex(ψωts) + 1
  ξp, na, nn, llr =
    _sim_gbmfbd_t(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, ψωts, ix, nep,
      δt, srδt, lc, lU, Iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr) && (treelength(ξc)!=treelength(ξp))
        
    llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ix)

    if lU < llr + llrLTTp != 0.0

      _fixrtip!(ξp, na) # fix random tip
      setni!(bi, na)    # set new ni

      return ξp, LTTp, llr + llrLTTp
    end
  end
  
  return ξp, LTT, -Inf
end




"""
     fsbi_f(bi    ::iBffs,
            ξc    ::iTfbd,
            α     ::Float64,
            σλ    ::Float64,
            σμ    ::Float64,
            ψ     ::Vector{Float64},
            ω     ::Vector{Float64},
            ψωts  ::Vector{Float64},
            ωtimes::Vector{Float64},
            LTT   ::Ltt,
            ixi   ::Int64,
            ixf   ::Int64,
            δt    ::Float64,
            srδt  ::Float64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_f(bi    ::iBffs,
                ξc    ::iTfbd,
                α     ::Float64,
                σλ    ::Float64,
                σμ    ::Float64,
                ψ     ::Vector{Float64},
                ω     ::Vector{Float64},
                ψωts  ::Vector{Float64},
                ωtimes::Vector{Float64},
                LTT   ::Ltt,
                ixi   ::Int64,
                ixf   ::Int64,
                δt    ::Float64,
                srδt  ::Float64)

   # forward simulation during branch length
  nep = lastindex(ψωts) + 1
  ξp, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ,
      ψωts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return ξp, LTT, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if lU < acr

    _fixrtip!(ξp, na) # fix random tip

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(ξp, tf(bi), α, σλ, σμ, ψ, ψωts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr

      # fossilize extant tip
      fossilizefixedtip!(ξp)

      # if terminal fossil branch
      tx, na, nn, acr =
        fossiltip_sim!(ξp, tf(bi), α, σλ, σμ, ψ, ψωts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)

      if lU < acr

        llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ixi)
          
        if lU < acr + llrLTTp != 0.0

          llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi)) + llrLTTp
          setnt!(bi, ntp)      # set new nt
          setni!(bi, na)       # set new ni

          return ξp, LTTp, llr
        end
      end
    end
  end

  return ξp, LTT, NaN
end




"""
    fsbi_et(ξp    ::iTfbd,
            bi    ::iBffs,
            α     ::Float64,
            σλ    ::Float64,
            σμ    ::Float64,
            ψ     ::Vector{Float64},
            ω     ::Vector{Float64},
            ψωts  ::Vector{Float64},
            ωtimes::Vector{Float64},
            LTT   ::Ltt,
            ixf   ::Int64,
            δt    ::Float64,
            srδt  ::Float64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_et(bi    ::iBffs,
                 ξc    ::iTfbd,
                 ξp    ::iTfbd,
                 α     ::Float64,
                 σλ    ::Float64,
                 σμ    ::Float64,
                 ψ     ::Vector{Float64},
                 ω     ::Vector{Float64},
                 ψωts  ::Vector{Float64},
                 ωtimes::Vector{Float64},
                 LTT   ::Ltt,
                 ixf   ::Int64,
                 δt    ::Float64,
                 srδt  ::Float64)

  nep = lastindex(ψωts) + 1
  lU  = -randexp()            # log-probability
  nac = ni(bi)                # current ni
  Iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # if terminal fossil branch
  tx, na, nn, acr =
    fossiltip_sim!(ξp, tf(bi), α, σλ, σμ, ψ, ψωts, ixf, nep, δt, srδt,
      acr, lU, Iρi, 1, 1)

  if lU < acr && (treelength(ξc)!=treelength(ξp))

    llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ixf)

    if lU < acr + llrLTTp != 0.0
      llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi)) + llrLTTp
      setni!(bi, na)       # set new ni

      return ξp, LTTp, llr
    end
  end

  return ξp, LTT, NaN
end



"""
    fsbi_m(bi    ::iBffs,
           ξc    ::iTfbd,
           ξ1    ::iTfbd,
           α     ::Float64,
           σλ    ::Float64,
           σμ    ::Float64,
           ψ     ::Vector{Float64},
           ω     ::Vector{Float64},
           ψωts  ::Vector{Float64},
           ωtimes::Vector{Float64},
           LTT   ::Ltt,
           ixi   ::Int64,
           ixf   ::Int64,
           δt    ::Float64,
           srδt  ::Float64)

Forward simulation for fossil internal branch `bi`.
"""
function fsbi_m(bi    ::iBffs,
                ξc    ::iTfbd,
                ξ1    ::iTfbd,
                α     ::Float64,
                σλ    ::Float64,
                σμ    ::Float64,
                ψ     ::Vector{Float64},
                ω     ::Vector{Float64},
                ψωts  ::Vector{Float64},
                ωtimes::Vector{Float64},
                LTT   ::Ltt,
                ixi   ::Int64,
                ixf   ::Int64,
                δt    ::Float64,
                srδt  ::Float64)

  # forward simulation during branch length
  nep = lastindex(ψωts) + 1
  ξp, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ,
      ψωts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return ξp, LTT, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample and fix random  tip
  λf, μf = fixrtip!(ξp, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, ssrμ, λ1p, μ1p, =
    _daughter_update!(ξ1, λf, μf, α, σλ, σμ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(ξp, tf(bi), α, σλ, σμ, ψ, ψωts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      # fossilize extant tip
      isfossil(bi) && fossilizefixedtip!(ξp)

      llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ixi)

      if lU < acr + llrLTTp != 0.0

        llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi)) + llrLTTp
        setnt!(bi, ntp)                       # set new nt
        setni!(bi, na)                        # set new ni
        l1 = lastindex(λ1p)
        unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
        unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector

        return ξp, LTTp, llr, drλ, ssrλ, ssrμ
      end
    end
  end

  return ξp, LTT, NaN, NaN, NaN, NaN
end




"""
    fsbi_i(bi    ::iBffs,
           ξc    ::iTfbd,
           ξ1    ::iTfbd,
           ξ2    ::iTfbd,
           α     ::Float64,
           σλ    ::Float64,
           σμ    ::Float64,
           ψ     ::Vector{Float64},
           ω     ::Vector{Float64},
           ψωts  ::Vector{Float64},
           ωtimes::Vector{Float64},
           LTT   ::Ltt,
           ixi   ::Int64,
           ixf   ::Int64,
           δt    ::Float64,
           srδt  ::Float64)

Forward simulation for internal branch `bi`.
"""
function fsbi_i(bi    ::iBffs,
                ξc    ::iTfbd,
                ξ1    ::iTfbd,
                ξ2    ::iTfbd,
                α     ::Float64,
                σλ    ::Float64,
                σμ    ::Float64,
                ψ     ::Vector{Float64},
                ω     ::Vector{Float64},
                ψωts  ::Vector{Float64},
                ωtimes::Vector{Float64},
                LTT   ::Ltt,
                ixi   ::Int64,
                ixf   ::Int64,
                δt    ::Float64,
                srδt  ::Float64)

  # forward simulation during branch length
  nep = lastindex(ψωts) + 1
  ξp, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ,
      ψωts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return ξp, LTT, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample and fix random  tip
  λf, μf = fixrtip!(ξp, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, ssrμ, λ1p, λ2p, μ1p, μ2p =
    _daughters_update!(ξ1, ξ2, λf, μf, α, σλ, σμ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(ξp, tf(bi), α, σλ, σμ, ψ, ψωts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr
      
      llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ixi)

      if lU < acr + llrLTTp != 0.0
        na -= 1

        llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi)) + llrLTTp
        setnt!(bi, ntp)                       # set new nt
        setni!(bi,  na)                       # set new ni
        setλt!(bi,  λf)                       # set new λt
        l1 = lastindex(λ1p)
        l2 = lastindex(λ2p)
        unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
        unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 1 λ vector
        unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector
        unsafe_copyto!(lμ(ξ2), 1, μ2p, 1, l2) # set new daughter 1 μ vector

        return ξp, LTTp, llr, drλ, ssrλ, ssrμ
      end
    end
  end

  return ξp, LTT, NaN, NaN, NaN, NaN
end





