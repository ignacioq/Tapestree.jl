#=

Anagenetic `fbdd` MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmfbd(tree    ::sTf_label;
                  λa_prior::NTuple{2,Float64}     = (1.5, 1.0),
                  μa_prior::NTuple{2,Float64}     = (1.5, 1.0),
                  αλ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                  αμ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                  σλ_prior::NTuple{2,Float64}     = (3.0, 0.5),
                  σμ_prior::NTuple{2,Float64}     = (3.0, 0.5),
                  ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                  ψ_epoch ::Vector{Float64}       = Float64[],
                  f_epoch ::Vector{Int64}         = Int64[0],
                  niter   ::Int64                 = 1_000,
                  nthin   ::Int64                 = 10,
                  nburn   ::Int64                 = 200,
                  nflushθ ::Int64                 = Int64(ceil(niter/5_000)),
                  nflushΞ ::Int64                 = Int64(ceil(niter/100)),
                  ofile   ::String                = string(homedir(), "/fbdd"),
                  tune_int::Int64                 = 100,
                  ϵi      ::Float64               = 0.2,
                  λi      ::Float64               = NaN,
                  μi      ::Float64               = NaN,
                  ψi      ::Float64               = NaN,
                  αλi      ::Float64              = 0.0,
                  αμi      ::Float64              = 0.0,
                  σλi     ::Float64               = 0.1,
                  σμi     ::Float64               = 0.1,
                  pupdp   ::NTuple{7,Float64}     = (0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.2),
                  δt      ::Float64               = 1e-3,
                  survival::Bool                  = true,
                  mxthf   ::Float64               = 0.1,
                  prints  ::Int64                 = 5,
                  stnλ    ::Float64               = 0.5,
                  stnμ    ::Float64               = 0.5,
                  tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for fossilized birth-death diffusion `fbdd`.
"""
function insane_gbmfbd(tree    ::sTf_label;
                       λa_prior::NTuple{2,Float64}     = (1.5, 1.0),
                       μa_prior::NTuple{2,Float64}     = (1.5, 1.0),
                       αλ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                       αμ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                       σλ_prior::NTuple{2,Float64}     = (3.0, 0.5),
                       σμ_prior::NTuple{2,Float64}     = (3.0, 0.5),
                       ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                       ψ_epoch ::Vector{Float64}       = Float64[],
                       f_epoch ::Vector{Int64}         = Int64[0],
                       niter   ::Int64                 = 1_000,
                       nthin   ::Int64                 = 10,
                       nburn   ::Int64                 = 200,
                       nflushθ ::Int64                 = Int64(ceil(niter/5_000)),
                       nflushΞ ::Int64                 = Int64(ceil(niter/100)),
                       ofile   ::String                = string(homedir(), "/fbdd"),
                       tune_int::Int64                 = 100,
                       ϵi      ::Float64               = 0.2,
                       λi      ::Float64               = NaN,
                       μi      ::Float64               = NaN,
                       ψi      ::Float64               = NaN,
                       αλi     ::Float64               = 0.0,
                       αμi     ::Float64               = 0.0,
                       σλi     ::Float64               = 0.1,
                       σμi     ::Float64               = 0.1,
                       pupdp   ::NTuple{7,Float64}     = (0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.2),
                       δt      ::Float64               = 1e-3,
                       survival::Bool                  = true,
                       mxthf   ::Float64               = 0.1,
                       prints  ::Int64                 = 5,
                       stnλ    ::Float64               = 0.5,
                       stnμ    ::Float64               = 0.5,
                       tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1,round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # only include epochs where the tree occurs
  sort!(ψ_epoch, rev = true)
  tix = findfirst(x -> x < th, ψ_epoch)
  if !isnothing(tix)
    ψ_epoch = ψ_epoch[tix:end]
  end
  nep  = lastindex(ψ_epoch) + 1

  # make initial fossils per epoch vector
  lep = lastindex(f_epoch)
  if lep !== nep
    if sum(f_epoch) > 0
      if lep > nep
        f_epoch = f_epoch[(end-nep+1):end]
      else 
        for i in Base.OneTo(nep - lep)
          pushfirst!(f_epoch, 0)
        end
      end
    else
      f_epoch = fill(0, nep)
    end
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
  if isnan(λi) || isnan(μi) || isnan(ψi)
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
      ψc = ψ_prior[1]/ψ_prior[2]
    else
      ψc = Float64(nf)/Float64(treelength(tree))
    end
  else
    λc, μc, ψc = λi, μi, ψi
  end

  # make ψ vector
  ψc = fill(ψc, nep)

  # condition on survival of 0, 1, or 2 starting lineages
  surv = 0
  if survival 
    if iszero(e(tree)) 
      if def1(tree)
        surv += Int64(anyalive(tree.d1))
        if def2(tree)
          surv += Int64(anyalive(tree.d2))
        end
      end
    else
      surv += Int64(anyalive(tree))
    end
  end

  # M attempts of survival
  mc = m_surv_gbmfbd(th, log(λc), log(μc), αλi, αμi, σλi, σμi, 
         δt, srδt, 1_000, surv)

  # make a decoupled tree
  Ξ = make_Ξ(idf, λc, μc, αλi, αμi, σλi, σμi, δt, srδt, iTfbd)

  # set end of fix branch speciation times and get vector of internal branches
  # and make epoch start vectors and indices for each `ξ`
  inodes = Int64[]
  eixi   = Int64[]
  eixf   = Int64[]
  bst    = Float64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    d1(bi) > 0 && push!(inodes, i)
    tib = ti(bi)
    ei  = findfirst(x -> x < tib, ψ_epoch)
    ei  = isnothing(ei) ? nep : ei
    ef  = findfirst(x -> x < tf(bi), ψ_epoch)
    ef  = isnothing(ef) ? nep : ef
    push!(bst, tib)
    push!(eixi, ei)
    push!(eixf, ef)
  end

  # parameter updates (1: αλ, 2: αμ, 3: σλ & σμ, 4: ψ, 5: scale, 6: gbm, 7: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running fossilized birth-death diffusion"

  # burn-in phase
  Ξ, idf, llc, prc, αλc, αμc, σλc, σμc, ψc, mc, ns, ne, stnλ, stnμ =
    mcmc_burn_gbmfbd(Ξ, idf, λa_prior, μa_prior, αλ_prior, αμ_prior, σλ_prior, σμ_prior,
      ψ_prior, ψ_epoch, f_epoch, nburn, tune_int, αλi, αμi, σλi, σμi, ψc, mc, th, surv, 
      stnλ, stnμ, δt, srδt, bst, eixi, eixf, inodes, pup, prints)

  # mcmc
  r, treev =
    mcmc_gbmfbd(Ξ, idf, llc, prc, αλc, αμc, σλc, σμc, ψc, mc, th, surv, ns, ne, 
      stnλ, stnμ, λa_prior, μa_prior, αλ_prior, αμ_prior, σλ_prior, σμ_prior, 
      ψ_prior, ψ_epoch, f_epoch, δt, srδt, bst, eixi, eixf, inodes, pup, 
      niter, nthin, nflushθ, nflushΞ, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_gbmfbd(Ξ       ::Vector{iTfbd},
                     idf     ::Vector{iBffs},
                     λa_prior::NTuple{2,Float64},
                     μa_prior::NTuple{2,Float64},
                     αλ_prior::NTuple{2,Float64},
                     αμ_prior::NTuple{2,Float64},
                     σλ_prior::NTuple{2,Float64},
                     σμ_prior::NTuple{2,Float64},
                     ψ_prior ::NTuple{2,Float64},
                     ψ_epoch ::Vector{Float64},
                     f_epoch ::Vector{Int64},
                     nburn   ::Int64,
                     tune_int::Int64,
                     αλc     ::Float64,
                     αμc     ::Float64,
                     σλc     ::Float64,
                     σμc     ::Float64,
                     ψc      ::Vector{Float64},
                     mc      ::Float64,
                     th      ::Float64,
                     surv    ::Int64,
                     stnλ    ::Float64, 
                     stnμ    ::Float64,
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
function mcmc_burn_gbmfbd(Ξ       ::Vector{iTfbd},
                          idf     ::Vector{iBffs},
                          λa_prior::NTuple{2,Float64},
                          μa_prior::NTuple{2,Float64},
                          αλ_prior::NTuple{2,Float64},
                          αμ_prior::NTuple{2,Float64},
                          σλ_prior::NTuple{2,Float64},
                          σμ_prior::NTuple{2,Float64},
                          ψ_prior ::NTuple{2,Float64},
                          ψ_epoch ::Vector{Float64},
                          f_epoch ::Vector{Int64},
                          nburn   ::Int64,
                          tune_int::Int64,
                          αλc     ::Float64,
                          αμc     ::Float64,
                          σλc     ::Float64,
                          σμc     ::Float64,
                          ψc      ::Vector{Float64},
                          mc      ::Float64,
                          th      ::Float64,
                          surv    ::Int64,
                          stnλ    ::Float64, 
                          stnμ    ::Float64,
                          δt      ::Float64,
                          srδt    ::Float64,
                          bst     ::Vector{Float64},
                          eixi    ::Vector{Int64},
                          eixf    ::Vector{Int64},
                          inodes  ::Array{Int64,1},
                          pup     ::Array{Int64,1},
                          prints  ::Int64)

  nsi = (iszero(e(Ξ[1])) && !isfossil(idf[1]))
  llc = llik_gbm(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) -
        nsi * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2,         σλ_prior[1], σλ_prior[2])  +
        logdinvgamma(σμc^2,         σμ_prior[1], σμ_prior[2])  +
        logdnorm(αλc,               αλ_prior[1], αλ_prior[2]^2) +
        logdnorm(αμc,               αμ_prior[1], αμ_prior[2]^2) +
        logdgamma(exp(lλ(Ξ[1])[1]), λa_prior[1], λa_prior[2]) +
        logdgamma(exp(lμ(Ξ[1])[1]), μa_prior[1], μa_prior[2]) +
        sum(logdgamma.(ψc,          ψ_prior[1],  ψ_prior[2]))

  L   = treelength(Ξ, ψ_epoch, bst, eixi)        # tree length
  nf  = nfossils(idf, ψ_epoch, f_epoch)          # number of fossilization events per epoch
  nin = lastindex(inodes)                        # number of internal nodes
  el  = lastindex(idf)                           # number of branches
  nep = lastindex(ψc)                            # number of epochs
  ns  = sum(x -> Float64(d2(x) > 0), idf) - nsi  # number of speciation events in likelihood
  ne  = Float64(ntipsextinct(Ξ))                 # number of extinction events in likelihood

  ddλ, ddμ, ssλ, ssμ, nλ = _ss_dd(Ξ, αλc, αμc)

  # for scale tuning
  ltn = 0
  lup = lacλ = lacμ = 0.0

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  function check_pr(pupi::Int64, i::Int64)
    pr0 = logdinvgamma(σλc^2,            σλ_prior[1], σλ_prior[2])  +
          logdinvgamma(σμc^2,            σμ_prior[1], σμ_prior[2])  +
          logdnorm(αλc,                  αλ_prior[1],  αλ_prior[2]^2) +
          logdnorm(αμc,                  αμ_prior[1],  αμ_prior[2]^2) +
          logdgamma(exp(lλ(Ξ[1])[1]),    λa_prior[1], λa_prior[2]) +
          logdgamma(exp(lμ(Ξ[1])[1]),    μa_prior[1], μa_prior[2]) +
          sum(logdgamma.(ψc, ψ_prior[1], ψ_prior[2]))
    if !isapprox(pr0, prc, atol = 1e-4)
       error(string("Wrong prior computation during the ", ["αλ","αμ","σλ & σμ","ψ","λ0&μ0","gbm update","forward simulation"][pupi], 
                    " update, at iteration ", i, ": pr0=", pr0, " and prc-pr0=", prc-pr0))
    end
  end

  function check_ll(pupi::Int64, i::Int64)
    ll0 = llik_gbm(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) - (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
    if !isapprox(ll0, llc, atol = 1e-4)
       error(string("Wrong likelihood computation during the ", ["αλ","αμ","σλ & σμ","ψ","λ0&μ0","gbm update","forward simulation"][pupi], 
                    " update, at iteration ", i, ": ll0=", ll0, " and llc-ll0=", llc-ll0))
       # @warn string("Wrong likelihood computation during the ", ["αλ","αμ","σλ & σμ","ψ","λ0&μ0","gbm update","forward simulation"][pupi], 
       #              " update, at iteration ", i, ": ll0=", ll0, " and llc-ll0=", llc-ll0)
    end
  end

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update αλ
      if pupi === 1

        llc, prc, αλc, mc =
          update_αλ!(αλc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αμc, σλc, σμc, sum(L), 
            ddλ, llc, prc, mc, th, surv, δt, srδt, αλ_prior)

        # update ssλ with new drift `αλc`
        ssλ = _ss(Ξ, lλ, αλc)

      # update αμ
      elseif pupi === 2

        llc, prc, αμc, mc =
          update_αμ!(αμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αλc, σλc, σμc, sum(L), 
            ddμ, llc, prc, mc, th, surv, δt, srδt, αμ_prior)

        # update ssλ with new drift `αμc`
        ssμ = _ss(Ξ, lμ, αμc)

      # σλ & σμ update
      elseif pupi === 3

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αλc, αμc, ssλ, ssμ, nλ,
            llc, prc, mc, th, surv, δt, srδt, σλ_prior, σμ_prior)

      # psi update
      elseif pupi === 4

        llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # update scale
      elseif pupi === 5

        llc, prc, accλ, accμ, mc = 
          update_scale!(Ξ, idf, αλc, αμc, σλc, σμc, llc, prc, ns, ne, 
            stnλ, stnμ, mc, th, surv, δt, srδt, λa_prior, μa_prior)

        lacλ += accλ
        lacμ += accμ
        lup += 1.0

      # gbm update
      elseif pupi === 6

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ddλ, ddμ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αλc, αμc, σλc, σμc, llc, prc, ddλ, ddμ, 
            ssλ, ssμ, mc, th, surv, δt, srδt, λa_prior, μa_prior)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, ddλ, ddμ, ssλ, ssμ, nλ, ns, ne, L =
          update_fs!(bix, Ξ, idf, αλc, αμc, σλc, σμc, ψc, llc, ddλ, ddμ, 
            ssλ, ssμ, nλ, ns, ne, L, ψ_epoch, δt, srδt, eixi, eixf)

      end

      # check_pr(pupi, i)
      # check_ll(pupi, i)

    end

    ltn += 1
    if ltn === tune_int

      stnλ = min(2.0, tune(stnλ, lacλ/lup))
      stnμ = min(2.0, tune(stnμ, lacμ/lup))
      ltn = 0
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αλc, αμc, σλc, σμc, ψc, mc, ns, ne, stnλ, stnμ
end




"""
    mcmc_gbmfbd(Ξ       ::Vector{iTfbd},
                idf     ::Vector{iBffs},
                llc     ::Float64,
                prc     ::Float64,
                αλc     ::Float64,
                αμc     ::Float64,
                σλc     ::Float64,
                σμc     ::Float64,
                ψc      ::Vector{Float64},
                mc      ::Float64,
                th      ::Float64,
                surv    ::Int64,
                ns      ::Float64,
                ne      ::Float64,
                stnλ    ::Float64, 
                stnμ    ::Float64,
                λa_prior::NTuple{2,Float64},
                μa_prior::NTuple{2,Float64},
                αλ_prior::NTuple{2,Float64},
                αμ_prior::NTuple{2,Float64},
                σλ_prior::NTuple{2,Float64},
                σμ_prior::NTuple{2,Float64},
                ψ_prior ::NTuple{2,Float64},
                ψ_epoch ::Vector{Float64},
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
                nflushθ ::Int64,
                nflushΞ ::Int64,
                ofile   ::String,
                prints  ::Int64)

MCMC chain for `fbdd`.
"""
function mcmc_gbmfbd(Ξ       ::Vector{iTfbd},
                     idf     ::Vector{iBffs},
                     llc     ::Float64,
                     prc     ::Float64,
                     αλc     ::Float64,
                     αμc     ::Float64,
                     σλc     ::Float64,
                     σμc     ::Float64,
                     ψc      ::Vector{Float64},
                     mc      ::Float64,
                     th      ::Float64,
                     surv    ::Int64,
                     ns      ::Float64,
                     ne      ::Float64,
                     stnλ    ::Float64, 
                     stnμ    ::Float64,
                     λa_prior::NTuple{2,Float64},
                     μa_prior::NTuple{2,Float64},
                     αλ_prior::NTuple{2,Float64},
                     αμ_prior::NTuple{2,Float64},
                     σλ_prior::NTuple{2,Float64},
                     σμ_prior::NTuple{2,Float64},
                     ψ_prior ::NTuple{2,Float64},
                     ψ_epoch ::Vector{Float64},
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
                     nflushθ ::Int64,
                     nflushΞ ::Int64,
                     ofile   ::String,
                     prints  ::Int64)

  # logging
  nlogs = fld(niter, nthin)
  lthin = lit = sthinθ = sthinΞ =  0

  L   = treelength(Ξ, ψ_epoch, bst, eixi) # tree length
  nf  = nfossils(idf, ψ_epoch, f_epoch)   # number of fossilization events per epoch
  nin = lastindex(inodes)                 # number of internal nodes
  el  = lastindex(idf)                    # number of branches
  nep = lastindex(ψc)

  ddλ, ddμ, ssλ, ssμ, nλ = _ss_dd(Ξ, αλc, αμc)

  # parameter results
  r = Array{Float64,2}(undef, nlogs, 9 + nep)

  treev = iTfbd[]    # make tree vector
  io    = IOBuffer() # buffer 

  function check_pr(pupi::Int64, i::Int64)
    pr0 = logdinvgamma(σλc^2,            σλ_prior[1], σλ_prior[2])  +
          logdinvgamma(σμc^2,            σμ_prior[1], σμ_prior[2])  +
          logdnorm(αλc,                  αλ_prior[1], αλ_prior[2]^2) +
          logdnorm(αμc,                  αμ_prior[1], αμ_prior[2]^2) +
          logdgamma(exp(lλ(Ξ[1])[1]),    λa_prior[1], λa_prior[2]) +
          logdgamma(exp(lμ(Ξ[1])[1]),    μa_prior[1], μa_prior[2]) +
          sum(logdgamma.(ψc, ψ_prior[1], ψ_prior[2]))
    if !isapprox(pr0, prc, atol = 1e-4)
       error(string("Wrong prior computation during the ", ["αλ","αμ","σλ & σμ","ψ","λ0&μ0","gbm update","forward simulation"][pupi], 
                    " update, at iteration ", i, ": pr0=", pr0, " and prc-pr0=", prc-pr0))
    end
  end

  function check_ll(pupi::Int64, i::Int64)
    ll0 = llik_gbm(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) - (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
    if !isapprox(ll0, llc, atol = 1e-4)
       error(string("Wrong likelihood computation during the ", ["αλ","αμ","σλ & σμ","ψ","λ0&μ0","gbm update","forward simulation"][pupi], 
                    " update, at iteration ", i, ": ll0=", ll0, " and llc-ll0=", llc-ll0))
       # @warn string("Wrong likelihood computation during the ", ["αλ","αμ","σλ & σμ","ψ","λ0&μ0","gbm update","forward simulation"][pupi], 
       #              " update, at iteration ", i, ": ll0=", ll0, " and llc-ll0=", llc-ll0)
    end
  end

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\tmu_root\talpha_lambda\talpha_mu\tsigma_lambda\tsigma_mu\t"*join(["psi"*(isone(nep) ? "" : string("_",i)) for i in 1:nep], '\t')*'\n')
    flush(of)

    open(ofile*".txt", "w") do tf


      pbar = Progress(niter, prints, "running mcmc...", 20)

      for it in Base.OneTo(niter)

        shuffle!(pup)

        # parameter updates
        for pupi in pup
          # @show ["αλ","αμ","σλ & σμ","ψ","λ0&μ0","gbm update","forward simulation"][pupi]

          # update αλ
          if pupi === 1
            llc, prc, αλc, mc =
              update_αλ!(αλc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αμc, σλc, σμc, sum(L), 
                ddλ, llc, prc, mc, th, surv, δt, srδt, αλ_prior)

            # update ssλ with new drift `αλc`
            ssλ = _ss(Ξ, lλ, αλc)

          # update αμ
          elseif pupi === 2

            llc, prc, αμc, mc =
              update_αμ!(αμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αλc, σλc, σμc, sum(L), 
                ddμ, llc, prc, mc, th, surv, δt, srδt, αμ_prior)

            # update ssλ with new drift `αμc`
            ssμ = _ss(Ξ, lμ, αμc)

          # σλ & σμ update
          elseif pupi === 3

            llc, prc, σλc, σμc, mc =
              update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αλc, αμc, ssλ, ssμ, 
                nλ, llc, prc, mc, th, surv, δt, srδt, σλ_prior, σμ_prior)

          # psi update
          elseif pupi === 4

            llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

          # update scale
          elseif pupi === 5

            llc, prc, accλ, accμ, mc = 
              update_scale!(Ξ, idf, αλc, αμc, σλc, σμc, llc, prc, ns, ne, 
                stnλ, stnμ, mc, th, surv, δt, srδt, λa_prior, μa_prior)

          # gbm update
          elseif pupi === 6

            nix = ceil(Int64,rand()*nin)
            bix = inodes[nix]

            llc, prc, ddλ, ddμ, ssλ, ssμ, mc =
              update_gbm!(bix, Ξ, idf, αλc, αμc, σλc, σμc, llc, prc, ddλ, ddμ, 
                ssλ, ssμ, mc, th, surv, δt, srδt, λa_prior, μa_prior)

          # forward simulation update
          else

            bix = ceil(Int64,rand()*el)

            llc, ddλ, ddμ, ssλ, ssμ, nλ, ns, ne, L =
              update_fs!(bix, Ξ, idf, αλc, αμc, σλc, σμc, ψc, llc, ddλ, ddμ, 
                ssλ, ssμ, nλ, ns, ne, L, ψ_epoch, δt, srδt, eixi, eixf)

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
            r[lit,6] = αλc
            r[lit,7] = αμc
            r[lit,8] = σλc
            r[lit,9] = σμc
            @turbo for i in Base.OneTo(nep)
              r[lit, 9 + i] = ψc[i]
            end
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthinθ += 1
        if sthinθ === nflushθ
          write(of, 
            string(Float64(it), "\t", llc, "\t", prc, "\t", 
              exp(lλ(Ξ[1])[1]),"\t", exp(lμ(Ξ[1])[1]), "\t", αλc, "\t", 
              αμc, "\t", σλc, "\t", σμc, "\t", join(ψc, "\t"), "\n"))
          flush(of)
          sthinθ = 0
        end
        sthinΞ += 1
        if sthinΞ === nflushΞ
          ibuffer(io, couple(Ξ, idf, 1))
          write(io, '\n')
          write(tf, take!(io))
          flush(tf)
          sthinΞ = 0
        end
        next!(pbar)
      end
    end
  end

  return r, treev
end




"""
    update_αλ!(αλ      ::Float64,
               λ0      ::Float64,
               μ0      ::Float64,
               αμ      ::Float64,
               σλ      ::Float64,
               σμ      ::Float64,
               L       ::Float64,
               ddλ     ::Float64,
               llc     ::Float64,
               prc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               surv    ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               αλ_prior::NTuple{2,Float64})

Gibbs update for `αλ`.
"""
function update_αλ!(αλc     ::Float64,
                    λ0      ::Float64,
                    μ0      ::Float64,
                    αμ      ::Float64,
                    σλ      ::Float64,
                    σμ      ::Float64,
                    L       ::Float64,
                    ddλ     ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    surv    ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    αλ_prior::NTuple{2,Float64})

  ν   = αλ_prior[1]
  τ2  = αλ_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2
  αλp  = rnorm((ddλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  mp  = m_surv_gbmfbd(th, λ0, μ0, αλp, αμ, σλ, σμ, δt, srδt, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σλ2*(αλc^2 - αλp^2 + 2.0*ddλ*(αλp - αλc)/L) + llr
    prc += llrdnorm_x(αλp, αλc, ν, τ2)
    αλc  = αλp
    mc   = mp
  end

  return llc, prc, αλc, mc
end




"""
    update_αμ!(αμ      ::Float64,
               λ0      ::Float64,
               μ0      ::Float64,
               αλ      ::Float64,
               σλ      ::Float64,
               σμ      ::Float64,
               L       ::Float64,
               ddλ     ::Float64,
               llc     ::Float64,
               prc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               surv    ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               αμ_prior::NTuple{2,Float64})

Gibbs update for `αμ`.
"""
function update_αμ!(αμc     ::Float64,
                    λ0      ::Float64,
                    μ0      ::Float64,
                    αλ      ::Float64,
                    σλ      ::Float64,
                    σμ      ::Float64,
                    L       ::Float64,
                    ddμ     ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    surv    ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    αμ_prior::NTuple{2,Float64})

  ν   = αμ_prior[1]
  τ2  = αμ_prior[2]^2
  σμ2 = σμ^2
  rs  = σμ2/τ2
  αμp  = rnorm((ddμ + rs*ν)/(rs + L), sqrt(σμ2/(rs + L)))

  mp  = m_surv_gbmfbd(th, λ0, μ0, αλ, αμp, σλ, σμ, δt, srδt, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σμ2*(αμc^2 - αμp^2 + 2.0*ddμ*(αμp - αμc)/L) + llr
    prc += llrdnorm_x(αμp, αμc, ν, τ2)
    αμc  = αμp
    mc   = mp
  end

  return llc, prc, αμc, mc
end




"""
    update_σ!(σλc     ::Float64,
              σμc     ::Float64,
              λ0      ::Float64,
              μ0      ::Float64,
              α       ::Float64,
              ssλ     ::Float64,
              ssμ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              mc      ::Float64,
              th      ::Float64,
              surv    ::Bool,
              δt      ::Float64,
              srδt    ::Float64,
              σλ_prior::NTuple{2,Float64},
              σμ_prior::NTuple{2,Float64})

Gibbs update for `σλ` and `σμ`.
"""
function update_σ!(σλc     ::Float64,
                   σμc     ::Float64,
                   λ0      ::Float64,
                   μ0      ::Float64,
                   αλ      ::Float64,
                   αμ      ::Float64,
                   ssλ     ::Float64,
                   ssμ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   surv    ::Int64,
                   δt      ::Float64,
                   srδt    ::Float64,
                   σλ_prior::NTuple{2,Float64},
                   σμ_prior::NTuple{2,Float64})

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_prior[1] + 0.5 * n, σλ_prior[2] + ssλ)
  σμp2 = randinvgamma(σμ_prior[1] + 0.5 * n, σμ_prior[2] + ssμ)

  σλp = sqrt(σλp2)
  σμp = sqrt(σμp2)

  mp  = m_surv_gbmfbd(th, λ0, μ0, αλ, αμ, σλp, σμp, δt, srδt, 1_000, surv)

  llr = log(mp/mc)

  if -randexp() < llr
    llc += ssλ*(1.0/σλc^2 - 1.0/σλp2) - n*(log(σλp/σλc)) +
           ssμ*(1.0/σμc^2 - 1.0/σμp2) - n*(log(σμp/σμc)) +
           llr
    prc += llrdinvgamma(σλp2, σλc^2, σλ_prior[1], σλ_prior[2]) +
           llrdinvgamma(σμp2, σμc^2, σμ_prior[1], σμ_prior[2])
    σλc  = σλp
    σμc  = σμp
    mc   = mp
  end

  return llc, prc, σλc, σμc, mc
end




"""
    update_scale!(Ξ       ::Vector{T},
                  idf     ::Vector{iBffs},
                  αλ      ::Float64,
                  αμ      ::Float64,
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
                  δt      ::Float64,
                  srδt    ::Float64,
                  λa_prior::NTuple{2,Float64},
                  μa_prior::NTuple{2,Float64}) where {T <: iTfbd}

Update scale for speciation and extinction.
"""
function update_scale!(Ξ       ::Vector{T},
                       idf     ::Vector{iBffs},
                       αλ      ::Float64,
                       αμ      ::Float64,
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
                       δt      ::Float64,
                       srδt    ::Float64,
                       λa_prior::NTuple{2,Float64},
                       μa_prior::NTuple{2,Float64}) where {T <: iTfbd}

  irλ, irμ = _ir(Ξ)

  accλ = accμ = 0.0

  lλ0c = lλ(Ξ[1])[1]
  lμ0c = lμ(Ξ[1])[1]

  # sample log(scaling factor)
  Δlλ  = randn()*stnλ
  lλ0p = lλ0c + Δlλ

  # likelihood ratio
  llr = ns * Δlλ + (1.0 - exp(Δlλ)) * irλ

  # prior ratio
  prr = llrdgamma(exp(lλ0p), exp(lλ0c), λa_prior[1], λa_prior[2])

  lU = -randexp()

  if lU < llr + prr + log(1000.0/mc)

    # add survival ratio
    mp = m_surv_gbmfbd(th, lλ0p, lμ0c, αλ, αμ, σλ, σμ, δt, srδt, 1_000, surv)
    llr += log(mp/mc)

    if lU < llr + prr
      accλ += 1.0
      llc  += llr
      prc  += prr
      mc    = mp
      scale_rate!(Ξ, lλ, Δlλ)
      scale_rate!(idf, Δlλ)
    end
  end

  # sample log(scaling factor)
  Δlμ  = randn()*stnμ
  lμ0p = lμ0c + Δlμ

  # likelihood ratio  
  llr = ne * Δlμ + (1.0 - exp(Δlμ)) * irμ
  prr = llrdgamma(exp(lμ0p), exp(lμ0c), μa_prior[1], μa_prior[2])

  lU = -randexp()

  if lU < llr + prr + log(1000.0/mc)

    # add survival ratio
    mp = m_surv_gbmfbd(th, lλ0c, lμ0p, αλ, αμ, σλ, σμ, δt, srδt, 1_000, surv)
    llr += log(mp/mc)

    if lU < llr + prr
      accμ += 1.0
      llc  += llr
      prc  += prr
      mc    = mp
      scale_rate!(Ξ, lμ, Δlμ)
    end
  end

  return llc, prc, accλ, accμ, mc
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTfbd},
                idf  ::Vector{iBffs},
                αλ   ::Float64,
                αμ   ::Float64,
                σλ   ::Float64,
                σμ   ::Float64,
                llc  ::Float64,
                prc  ::Float64,
                ddλ  ::Float64,
                ddμ  ::Float64,
                ssλ  ::Float64,
                ssμ  ::Float64,
                mc   ::Float64,
                th   ::Float64,
                surv ::Int64,
                δt   ::Float64,
                srδt ::Float64,
                λa_prior::NTuple{2,Float64},
                μa_prior::NTuple{2,Float64})

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTfbd},
                     idf  ::Vector{iBffs},
                     αλ   ::Float64,
                     αμ   ::Float64,
                     σλ   ::Float64,
                     σμ   ::Float64,
                     llc  ::Float64,
                     prc  ::Float64,
                     ddλ  ::Float64,
                     ddμ  ::Float64,
                     ssλ  ::Float64,
                     ssμ  ::Float64,
                     mc   ::Float64,
                     th   ::Float64,
                     surv ::Int64,
                     δt   ::Float64,
                     srδt ::Float64,
                     λa_prior::NTuple{2,Float64},
                     μa_prior::NTuple{2,Float64})

  @inbounds begin
    ξi   = Ξ[bix]
    bi   = idf[bix]
    i1   = d1(bi)
    i2   = d2(bi)
    ξ1   = Ξ[i1]
    root = iszero(pa(bi))

    if root && iszero(e(bi))

      # if stem fossil
      if isfossil(bi)
        llc, prc, ddλ, ddμ, ssλ, ssμ, mc =
          _fstem_update!(ξi, ξ1, αλ, αμ, σλ, σμ, llc, prc, ddλ, ddμ, ssλ, ssμ, 
            mc, th, δt, srδt, λa_prior, μa_prior, surv)
      # if crown
      else
        llc, prc, ddλ, ddμ, ssλ, ssμ, mc =
          _crown_update!(ξi, ξ1, Ξ[i2], αλ, αμ, σλ, σμ, llc, prc, ddλ, ddμ, ssλ, ssμ, 
            mc, th, δt, srδt, λa_prior, μa_prior, surv)
        setλt!(bi, lλ(ξi)[1])
      end
    else
      # if stem
      if root
        llc, prc, ddλ, ddμ, ssλ, ssμ, mc =
          _stem_update!(ξi, αλ, αμ, σλ, σμ, llc, prc, ddλ, ddμ, ssλ, ssμ,
            mc, th, δt, srδt, λa_prior, μa_prior, surv)
      end

      # updates within the parent branch
      llc, ddλ, ddμ, ssλ, ssμ =
        _update_gbm!(ξi, αλ, αμ, σλ, σμ, llc, ddλ, ddμ, ssλ, ssμ, 
          δt, srδt, false)

      # get fixed tip
      lξi = fixtip(ξi)

      # if mid branch
      if iszero(i2)

        # make between decoupled trees duo node update
        llc, ssλ, ssμ =
          update_duo!(lλ(lξi), lλ(ξ1), lμ(lξi), lμ(ξ1), e(lξi), e(ξ1),
            fdt(lξi), fdt(ξ1), αλ, αμ, σλ, σμ, llc, ssλ, ssμ, δt, srδt)

      # if internal branch
      else
        ξ2 = Ξ[i2]
        # make between decoupled trees trio node update
        llc, ddλ, ddμ, ssλ, ssμ, λf =
          update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), lμ(lξi), lμ(ξ1), lμ(ξ2),
            e(lξi), e(ξ1), e(ξ2), fdt(lξi), fdt(ξ1), fdt(ξ2),
            αλ, αμ, σλ, σμ, llc, ddλ, ddμ, ssλ, ssμ, δt, srδt)

        # set fixed `λ(t)` in branch
        setλt!(bi, λf)
      end
    end

    # carry on updates in the daughters
    llc, ddλ, ddμ, ssλ, ssμ =
      _update_gbm!(ξ1, αλ, αμ, σλ, σμ, llc, ddλ, ddμ, ssλ, ssμ, δt, srδt,
        iszero(d1(idf[i1])))
    if i2 > 0
      llc, ddλ, ddμ, ssλ, ssμ =
        _update_gbm!(Ξ[i2], αλ, αμ, σλ, σμ, llc, ddλ, ddμ, ssλ, ssμ, 
          δt, srδt, iszero(d1(idf[i2])))
    end
  end

  return llc, prc, ddλ, ddμ, ssλ, ssμ, mc
end




"""
    update_fs!(bix ::Int64,
               Ξ   ::Vector{iTfbd},
               idf ::Vector{iBffs},
               αλ  ::Float64,
               αμ  ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               ψ   ::Vector{Float64},
               llc ::Float64,
               ddλ ::Float64,
               ddμ ::Float64,
               ssλ ::Float64,
               ssμ ::Float64,
               nλ  ::Float64,
               ns  ::Float64,
               ne  ::Float64,
               L   ::Vector{Float64},
               ψts ::Vector{Float64},
               δt  ::Float64,
               srδt::Float64,
               eixi::Vector{Int64},
               eixf::Vector{Int64})

Forward simulation proposal function for `gbmfbd`.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{iTfbd},
                    idf ::Vector{iBffs},
                    αλ  ::Float64,
                    αμ  ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Vector{Float64},
                    llc ::Float64,
                    ddλ ::Float64,
                    ddμ ::Float64,
                    ssλ ::Float64,
                    ssμ ::Float64,
                    nλ  ::Float64,
                    ns  ::Float64,
                    ne  ::Float64,
                    L   ::Vector{Float64},
                    ψts ::Vector{Float64},
                    δt  ::Float64,
                    srδt::Float64,
                    eixi::Vector{Int64},
                    eixf::Vector{Int64})

  bi  = idf[bix]
  ξc  = Ξ[bix]
  ixi = eixi[bix]

  # terminal branch
  if iszero(d1(bi))

    drλ = drμ = ssrλ = ssrμ = 0.0
    # fossil terminal branch
    if isfossil(bi)
      ixf = eixf[bix]

      ξp, llr = fsbi_t(bi, ξc, αλ, αμ, σλ, σμ, ψ, ψts, ixi, ixf, δt, srδt)

      # if terminal but not successful proposal, update extinct
      if !isfinite(llr)
        ξp, llr = fsbi_et(iTfbd_wofe(ξc), bi, αλ, αμ, σλ, σμ, ψ, ψts, ixf,
          δt, srδt)
      end

    # non-fossil terminal branch
    else
      ξp, llr = fsbi_t(bi, ξc, αλ, αμ, σλ, σμ, ψ, ψts, ixi, δt, srδt)
    end

  # internal non-bifurcating branch
  elseif iszero(d2(bi))

    ξp, llr, drλ, drμ, ssrλ, ssrμ =
      fsbi_m(bi, ξc, Ξ[d1(bi)], αλ, αμ, σλ, σμ, ψ, ψts, ixi, eixf[bix], δt, srδt)

  # internal bifurcating branch
  else

    ξp, llr, drλ, drμ, ssrλ, ssrμ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], αλ, αμ, σλ, σμ, ψ, ψts, 
        ixi, eixf[bix], δt, srδt)
  end

  if isfinite(llr)
    tii = ti(bi)

    nep = lastindex(ψts) + 1

    ll1, ixd, ddλ1, ddμ1, ssλ1, ssμ1, nλ1, ns1, ne1 =
      llik_gbm_ss(ξp, αλ, αμ, σλ, σμ, ψ, tii, ψts, ixi, δt, srδt, nep, 0.0, 0.0)
    ll0, ixd, ddλ0, ddμ0, ssλ0, ssμ0, nλ0, ns0, ne0 =
      llik_gbm_ss(ξc, αλ, αμ, σλ, σμ, ψ, tii, ψts, ixi, δt, srδt, nep, 0.0, 0.0)

    # update quantities
    llc += ll1  - ll0  + llr
    ddλ += ddλ1 - ddλ0 + drλ
    ddμ += ddμ1 - ddμ0 + drμ
    ssλ += ssλ1 - ssλ0 + ssrλ
    ssμ += ssμ1 - ssμ0 + ssrμ
    nλ  += nλ1  - nλ0
    ns  += ns1  - ns0
    ne  += ne1  - ne0

    # update tree lengths
    Lc = zeros(Float64, nep)
    _treelength!(ξc, tii, Lc, ψts, ixi, nep)
    _treelength!(ξp, tii, L,  ψts, ixi, nep)
    @turbo for i in Base.OneTo(nep)
      L[i] -= Lc[i]
    end

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ddμ, ssλ, ssμ, nλ, ns, ne, L
end




"""
    fsbi_t(bi::iBffs,
           ξc  ::iTfbd,
           αλ  ::Float64,
           αμ  ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Vector{Float64},
           ψts ::Vector{Float64},
           ix  ::Int64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs,
                ξc  ::iTfbd,
                αλ  ::Float64,
                αμ  ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ix  ::Int64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nn, llr =
    _sim_gbmfbd_t(e(bi), lλ(ξc)[1], lμ(ξc)[1], αλ, αμ, σλ, σμ, ψ, ψts, ix, nep,
      δt, srδt, lc, lU, Iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr)
    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, -Inf
  end
end




"""
     fsbi_t(bi  ::iBffs,
            ξc  ::iTfbd,
            αλ  ::Float64,
            αμ  ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            ψ   ::Vector{Float64},
            ψts ::Vector{Float64},
            ix  ::Int64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTfbd,
                αλ  ::Float64,
                αμ  ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ixi ::Int64,
                ixf ::Int64,
                δt  ::Float64,
                srδt::Float64)

   # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], αλ, αμ, σλ, σμ, ψ,
      ψts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return t0, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if lU < acr

    _fixrtip!(t0, na) # fix random tip

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr

      # fossilize extant tip
      fossilizefixedtip!(t0)

      # if terminal fossil branch
      tx, na, nn, acr =
        fossiltip_sim!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)

      if lU < acr

        llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
        setnt!(bi, ntp)      # set new nt
        setni!(bi, na)       # set new ni

        return t0, llr
      end
    end
  end

  return t0, NaN
end




"""
    fsbi_et(t0  ::iTfbd,
            bi  ::iBffs,
            αλ  ::Float64,
            αμ  ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            ψ   ::Vector{Float64},
            ψts ::Vector{Float64},
            ixf ::Int64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_et(t0  ::iTfbd,
                 bi  ::iBffs,
                 αλ  ::Float64,
                 αμ  ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 ψ   ::Vector{Float64},
                 ψts ::Vector{Float64},
                 ixf ::Int64,
                 δt  ::Float64,
                 srδt::Float64)

  nep = lastindex(ψts) + 1
  lU  = -randexp()            # log-probability
  nac = ni(bi)                # current ni
  Iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # if terminal fossil branch
  tx, na, nn, acr =
    fossiltip_sim!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
      acr, lU, Iρi, 1, 1)

  if lU < acr

    llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
    setni!(bi, na)       # set new ni

    return t0, llr
  end

  return t0, NaN
end



"""
    fsbi_m(bi  ::iBffs,
           ξc  ::iTfbd,
           ξ1  ::iTfbd,
           αλ  ::Float64,
           αμ  ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Vector{Float64},
           ψts ::Vector{Float64},
           ixi ::Int64,
           ixf ::Int64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for fossil internal branch `bi`.
"""
function fsbi_m(bi  ::iBffs,
                ξc  ::iTfbd,
                ξ1  ::iTfbd,
                αλ  ::Float64,
                αμ  ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ixi ::Int64,
                ixf ::Int64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], αλ, αμ, σλ, σμ, ψ,
      ψts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return t0, NaN, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample and fix random  tip
  λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, drμ, ssrλ, ssrμ, λ1p, μ1p =
    _daughter_update!(ξ1, λf, μf, αλ, αμ, σλ, σμ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      # fossilize extant tip
      isfossil(bi) && fossilizefixedtip!(t0)

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setnt!(bi, ntp)                       # set new nt
      setni!(bi, na)                        # set new ni
      l1 = lastindex(λ1p)
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector

      return t0, llr, drλ, drμ, ssrλ, ssrμ
    end
  end

  return t0, NaN, NaN, NaN, NaN, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTfbd,
           ξ1  ::iTfbd,
           ξ2  ::iTfbd,
           αλ  ::Float64,
           αμ  ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Vector{Float64},
           ψts ::Vector{Float64},
           ixi ::Int64,
           ixf ::Int64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for internal branch `bi`.
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTfbd,
                ξ1  ::iTfbd,
                ξ2  ::iTfbd,
                αλ  ::Float64,
                αμ  ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ixi ::Int64,
                ixf ::Int64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], αλ, αμ, σλ, σμ, ψ,
      ψts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return t0, NaN, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample and fix random  tip
  λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, drμ, ssrλ, ssrμ, λ1p, λ2p, μ1p, μ2p =
    _daughters_update!(ξ1, ξ2, λf, μf, αλ, αμ, σλ, σμ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setnt!(bi, ntp)                       # set new nt
      setni!(bi,  na)                       # set new ni
      setλt!(bi,  λf)                       # set new λt
      l1 = lastindex(λ1p)
      l2 = lastindex(λ2p)
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 1 λ vector
      unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector
      unsafe_copyto!(lμ(ξ2), 1, μ2p, 1, l2) # set new daughter 1 μ vector

      return t0, llr, drλ, drμ, ssrλ, ssrμ
    end
  end

  return t0, NaN, NaN, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTfbd,
              t   ::Float64,
              αλ  ::Float64,
              αμ  ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              ψ   ::Vector{Float64},
              ψts ::Vector{Float64},
              ix  ::Int64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTfbd,
                   t   ::Float64,
                   αλ  ::Float64,
                   αμ  ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   ψ   ::Vector{Float64},
                   ψts ::Vector{Float64},
                   ix  ::Int64,
                   nep ::Int64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

  if lU < lr && nn < 1_000

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)
        lμ0  = lμ(tree)
        l    = lastindex(lλ0)

        # simulate
        stree, na, nn, lr =
          _sim_gbmfbd_it(max(δt-fdti, 0.0), t, lλ0[l], lμ0[l], αλ, αμ, σλ, σμ, 
            ψ, ψts, ix, nep, δt, srδt, lr, lU, Iρi, na-1, nn, 1_000)

        if !isfinite(lr) || nn > 999
          return tree, na, nn, NaN
        end

        setproperty!(tree, :iμ, isextinct(stree))
        sete!(tree, e(tree) + e(stree))

        lλs = lλ(stree)
        lμs = lμ(stree)

        if lastindex(lλs) === 2
          setfdt!(tree, fdt(tree) + fdt(stree))
        else
          setfdt!(tree, fdt(stree))
        end

        pop!(lλ0)
        pop!(lμ0)
        popfirst!(lλs)
        popfirst!(lμs)
        append!(lλ0, lλs)
        append!(lμ0, lμs)

        # merge to current tip
        if def1(stree)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr =
        tip_sims!(tree.d1, t, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr =
        tip_sims!(tree.d2, t, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::iTfbd,
                   t   ::Float64,
                   αλ  ::Float64,
                   αμ  ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   ψ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

Continue simulation until time `t` for the fixed tip in `tree`.
"""
function fossiltip_sim!(tree::iTfbd,
                        t   ::Float64,
                        αλ  ::Float64,
                        αμ  ::Float64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        ψ   ::Vector{Float64},
                        ψts ::Vector{Float64},
                        ix  ::Int64,
                        nep ::Int64,
                        δt  ::Float64,
                        srδt::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        Iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64)


  if lU < lr && nn < 1_000
    if istip(tree)

      nep = lastindex(ψts) + 1
      stree, na, nn, lr =
        _sim_gbmfbd_it(t, lλ(tree)[end], lμ(tree)[end], αλ, αμ, σλ, σμ, ψ,
          ψts, ix, nep, δt, srδt, lr, lU, Iρi, na-1, nn, 1_000)

      if !isfinite(lr) || nn > 999
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr =
        fossiltip_sim!(tree.d1, t, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, Iρi, na, nn)
    else
      tree.d2, na, nn, lr =
        fossiltip_sim!(tree.d2, t, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end



