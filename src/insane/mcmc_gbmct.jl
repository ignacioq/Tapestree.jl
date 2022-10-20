#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmct(tree    ::sT_label,
                 out_file::String;
                 λa_prior::NTuple{2,Float64} = (1.0, 1.0),
                 α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                 ϵ_prior ::NTuple{2,Float64} = (0.0, 10.0),
                 niter   ::Int64             = 1_000,
                 nthin   ::Int64             = 10,
                 nburn   ::Int64             = 200,
                 tune_int::Int64             = 100,
                 αi      ::Float64           = 0.0,
                 λi      ::Float64           = NaN,
                 σλi     ::Float64           = 0.01,
                 ϵi      ::Float64           = 0.2,
                 ϵtni    ::Float64           = 1.0,
                 obj_ar  ::Float64           = 0.234,
                 pupdp   ::NTuple{5,Float64} = (0.1,0.1,0.1,0.2,0.2),
                 ntry    ::Int64             = 2,
                 nlim    ::Int64             = 500,
                 δt      ::Float64           = 5e-3,
                 prints  ::Int64             = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for GBM birth-death.
"""
function insane_gbmct(tree    ::sT_label,
                      out_file::String;
                      λa_prior::NTuple{2,Float64} = (1.0, 1.0),
                      α_prior ::NTuple{2,Float64} = (0.0, 0.5),
                      σλ_prior::NTuple{2,Float64} = (3.0, 0.5),
                      ϵ_prior ::NTuple{2,Float64} = (0.0, 10.0),
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      tune_int::Int64             = 100,
                      αi      ::Float64           = 0.0,
                      λi      ::Float64           = NaN,
                      σλi     ::Float64           = 0.01,
                      ϵi      ::Float64           = 0.2,
                      ϵtni    ::Float64           = 0.1,
                      obj_ar  ::Float64           = 0.234,
                      pupdp   ::NTuple{5,Float64} = (0.01, 0.01, 0.01, 0.1, 0.2),
                      ntry    ::Int64             = 2,
                      nlim    ::Int64             = 500,
                      δt      ::Float64           = 1e-3,
                      prints  ::Int64             = 5,
                      survival::Bool              = true,
                      mxthf   ::Float64           = Inf,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)
  crown = survival ? Int64(iszero(e(tree))) : 2

  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # estimate branch split (multiple of δt)
  ndts = floor(th * mxthf/δt)
  maxt = δt * ndts

  # make fix tree directory
  idf = make_idf(tree, tρ, maxt)

   # starting parameters (using method of moments)
  if isnan(λi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc = λi
  end
  ϵc = ϵi

  # make a decoupled tree
  Ξ = make_Ξ(idf, λc, αi, σλi, δt, srδt, iTct)

  # survival
  mc = m_surv_gbmct(th, log(λc), αi, σλi, ϵc, δt, srδt, 1_000, crown)

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf))  if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σλ, 3: ϵ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  @info "running birth-death gbm with constant ϵ"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, ϵc, ϵtn, mc =
    mcmc_burn_gbmct(Ξ, idf, λa_prior, α_prior, σλ_prior, ϵ_prior,
      nburn, tune_int, αi, σλi, ϵc, ϵtni, mc, th, crown, δt, srδt, inodes, pup,
       prints, scalef)

  # mcmc
  R, Ξv =
    mcmc_gbmct(Ξ, idf, llc, prc, αc, σλc, ϵc, ϵtn, mc, th, crown,
      λa_prior, α_prior, σλ_prior, ϵ_prior, niter, nthin, δt, srδt,
      inodes, pup, prints)

  pardic = Dict(("lambda_root"   => 1,
                 "alpha"        => 2,
                 "sigma_lambda" => 3,
                 "epsilon"      => 4))

  write_ssr(R, pardic, out_file)

  return R, Ξv
end




"""
    mcmc_burn_gbmct(Ξ       ::Vector{iTct},
                    idf     ::Vector{iBffs},
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    ϵ_prior ::NTuple{2,Float64},
                    nburn   ::Int64,
                    tune_int::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    ϵc      ::Float64,
                    ϵtn     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    crown   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Vector{Int64},
                    pup     ::Vector{Int64},
                    prints  ::Int64,
                    scalef  ::Function)

MCMC burn-in chain for `gbmct`.
"""
function mcmc_burn_gbmct(Ξ       ::Vector{iTct},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         ϵ_prior ::NTuple{2,Float64},
                         nburn   ::Int64,
                         tune_int::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         ϵc      ::Float64,
                         ϵtn     ::Float64,
                         mc      ::Float64,
                         th      ::Float64,
                         crown   ::Int64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Vector{Int64},
                         pup     ::Vector{Int64},
                         prints  ::Int64,
                         scalef  ::Function)

  ltn = 0
  lup = 0.0
  lac = 0.0

  λ0  = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) - Float64(crown > 0) * λ0 + 
        log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])  +
        logdgamma(exp(λ0),   λa_prior[1], λa_prior[2])  +
        logdnorm(αc,        α_prior[1],  α_prior[2]^2) +
        logdunif(ϵc,        ϵ_prior[1],  ϵ_prior[2])

  # maximum bounds according to uniform priors
  ϵxpr  = ϵ_prior[2]


  L       = treelength(Ξ)      # tree length
  dlλ     = deltaλ(Ξ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)     # sum squares in λ
  Σλ      = Σλ_gbm(Ξ)          # sum of λ
  ne      = 0.0                # number of extinction events
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  # number of branches
  nbr  = lastindex(idf)

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      if pupi === 1

        llc, prc, αc, mc =
          update_α_ϵ!(αc, lλ(Ξ[1])[1], σλc, ϵc, L, dlλ, llc, prc, mc, th, 
            crown, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

      elseif pupi === 2

        llc, prc, σλc, mc =
          update_σ_ϵ!(σλc, lλ(Ξ[1])[1], αc, ϵc, ssλ, nλ, llc, prc, mc, th, 
            crown, δt, srδt, σλ_prior)

      elseif pupi === 3

        llc, ϵc, mc, lac =
          update_ϵ!(ϵc, lλ(Ξ[1])[1], αc, σλc, llc, mc, th, crown, ϵtn,
            lac, ne, Σλ, δt, srδt, ϵxpr)

        lup += 1.0

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dlλ, ssλ, Σλ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, ϵc, llc, dlλ, ssλ, Σλ, mc, th,
            δt, srδt)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dlλ, ssλ, Σλ, nλ, ne, L =
          update_fs!(bix, Ξ, idf, αc, σλc, ϵc, llc, dlλ, ssλ, Σλ, nλ, ne, L,
            δt, srδt)
      end
    end

    # log tuning parameters
    ltn += 1
    if ltn === tune_int
      ϵtn = scalef(ϵtn, lac/lup)
      ltn = 0
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, ϵc, ϵtn, mc
end




"""
    mcmc_gbmct(Ξ       ::Vector{iTct},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               ϵc      ::Float64,
               ϵtn     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               crown   ::Int64,
               λa_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               ϵ_prior ::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Array{Int64,1},
               prints  ::Int64)

MCMC chain for `gbmct`.
"""
function mcmc_gbmct(Ξ       ::Vector{iTct},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    ϵc      ::Float64,
                    ϵtn     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    crown   ::Int64,
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    ϵ_prior ::NTuple{2,Float64},
                    niter   ::Int64,
                    nthin   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # maximum bounds according to uniform priors
  ϵxpr  = ϵ_prior[2]

  L       = treelength(Ξ)            # tree length
  dlλ     = deltaλ(Ξ)                # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)           # sum squares in λ
  Σλ      = Σλ_gbm(Ξ)                # sum of λ
  ne      = Float64(ntipsextinct(Ξ)) # number of extinction events
  nin     = lastindex(inodes)        # number of internal nodes
  el      = lastindex(idf)           # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 7)

  # make Ξ vector
  Ξv = iTct[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      if pupi === 1

        llc, prc, αc, mc =
          update_α_ϵ!(αc, lλ(Ξ[1])[1], σλc, ϵc, L, dlλ, llc, prc, mc, th, crown,
            δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) + log(mc) + prob_ρ(idf) - Float64(crown > 0) * lλ(Ξ[1])[1]
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      elseif pupi === 2

        llc, prc, σλc, mc =
          update_σ_ϵ!(σλc, lλ(Ξ[1])[1], αc, ϵc, ssλ, nλ, llc, prc, mc, th, crown,
            δt, srδt, σλ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) + log(mc) + prob_ρ(idf) - Float64(crown > 0) * lλ(Ξ[1])[1]
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      elseif pupi === 3

        llc, ϵc, mc =
          update_ϵ!(ϵc, lλ(Ξ[1])[1], αc, σλc, llc, mc, th, crown, ϵtn,
            ne, Σλ, δt, srδt, ϵxpr)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) + log(mc) + prob_ρ(idf) - Float64(crown > 0) * lλ(Ξ[1])[1]
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dlλ, ssλ, Σλ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, ϵc, llc, dlλ, ssλ, Σλ, mc, th,
            δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) + log(mc) + prob_ρ(idf) - Float64(crown > 0) * lλ(Ξ[1])[1]
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dlλ, ssλ, Σλ, nλ, ne, L =
          update_fs!(bix, Ξ, idf, αc, σλc, ϵc, llc, dlλ, ssλ, Σλ, nλ, ne, L,
            δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) + log(mc) + prob_ρ(idf) - Float64(crown > 0) * lλ(Ξ[1])[1]
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end
      end

    end

    # log parameters
    lthin += 1
    if lthin === nthin
      lit += 1
      @inbounds begin
        R[lit,1] = Float64(lit)
        R[lit,2] = llc
        R[lit,3] = prc
        R[lit,4] = exp(lλ(Ξ[1])[1])
        R[lit,5] = αc
        R[lit,6] = σλc
        R[lit,7] = ϵc
        push!(Ξv, couple(Ξ, idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ξv
end




"""
    update_fs!(bix    ::Int64,
               Ξ      ::Vector{iTct},
               idf    ::Vector{iBffs},
               α      ::Float64,
               σλ     ::Float64,
               ϵ      ::Float64,
               llc    ::Float64,
               dlλ     ::Float64,
               ssλ    ::Float64,
               Σλ     ::Float64,
               nλ     ::Float64,
               ne     ::Float64,
               L      ::Float64,
               δt     ::Float64,
               srδt   ::Float64)

Forward simulation proposal function for `gbmct`.
"""
function update_fs!(bix    ::Int64,
                    Ξ      ::Vector{iTct},
                    idf    ::Vector{iBffs},
                    α      ::Float64,
                    σλ     ::Float64,
                    ϵ      ::Float64,
                    llc    ::Float64,
                    dlλ     ::Float64,
                    ssλ    ::Float64,
                    Σλ     ::Float64,
                    nλ     ::Float64,
                    ne     ::Float64,
                    L      ::Float64,
                    δt     ::Float64,
                    srδt   ::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]

  # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, ξc, α, σλ, ϵ, δt, srδt)
    drλ = ssrλ = Σrλ  = 0.0

  # if mid
  elseif iszero(d2(bi))
    ξp, llr, drλ, ssrλ, Σrλ =
      fsbi_m(bi, ξc, Ξ[d1(bi)], α, σλ, ϵ, δt, srδt)

  # if internal
  else
    ξp, llr, drλ, ssrλ, Σrλ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, ϵ, δt, srδt)
  end

  if isfinite(llr)
    ll1, dlλ1, ssλ1, Σλ1, nλ1 = llik_gbm_ssλ(ξp, α, σλ, ϵ, δt, srδt)
    ll0, dlλ0, ssλ0, Σλ0, nλ0 = llik_gbm_ssλ(ξc, α, σλ, ϵ, δt, srδt)

    # update llc, ssλ, nλ, Σλ, ne, L
    llc += llr  + ll1  - ll0
    dlλ += dlλ1  - dlλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    Σλ  += Σλ1  - Σλ0  + Σrλ
    nλ  += nλ1  - nλ0
    ne  += ntipsextinct(ξp) - ntipsextinct(ξc)
    L   += treelength(ξp)   - treelength(ξc)

    Ξ[bix] = ξp          # set new tree
  end

  return llc, dlλ, ssλ, Σλ, nλ, ne, L
end




"""
    fsbi_t(bi  ::iBffs,
           ξc  ::iTct,
           α   ::Float64,
           σλ  ::Float64,
           ϵ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTct,
                α   ::Float64,
                σλ  ::Float64,
                ϵ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_gbmct_t(e(bi), lλ(ξc)[1], α, σλ, ϵ, δt, srδt, lc, lU, Iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr)
    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, -Inf
  end
end




"""
    fsbi_m(bi  ::iBffs,
           ξc  ::iTct,
           ξ1  ::iTct,
           ξ2  ::iTct,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           ϵ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_m(bi  ::iBffs,
                ξc  ::iTct,
                ξ1  ::iTct,
                α   ::Float64,
                σλ  ::Float64,
                ϵ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  t0, na, nn =
    _sim_gbmct(e(bi), lλ(ξc)[1], α, σλ, ϵ, δt, srδt, 0, 1, 1_000)

  if na < 1 || nn > 999
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # fix random tip
  λf = fixrtip!(t0, na, NaN)

  llrd, acrd, drλ, ssrλ, Σrλ, λ1p =
    _daughter_update!(ξ1, λf, α, σλ, ϵ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, ϵ, δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      l1  = lastindex(λ1p)
      setnt!(bi, ntp)                    # set new nt
      setni!(bi, na)                     # set new ni
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector

      return t0, llr, drλ, ssrλ, Σrλ
    else
      return t0, NaN, NaN, NaN, NaN
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTct,
           ξ1  ::iTct,
           ξ2  ::iTct,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           ϵ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTct,
                ξ1  ::iTct,
                ξ2  ::iTct,
                α   ::Float64,
                σλ  ::Float64,
                ϵ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  t0, na, nn =
    _sim_gbmct(e(bi), lλ(ξc)[1], α, σλ, ϵ, δt, srδt, 0, 1, 1_000)

  if na < 1 || nn > 999
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # fix random tip
  λf = fixrtip!(t0, na, NaN)

  llrd, acrd, drλ, ssrλ, Σrλ, λ1p, λ2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, ϵ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, ϵ, δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      l1  = lastindex(λ1p)
      l2  = lastindex(λ2p)
      setnt!(bi, ntp)                    # set new nt
      setni!(bi, na)                     # set new ni
      setλt!(bi, λf)                     # set new λt
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 2 λ vector

      return t0, llr, drλ, ssrλ, Σrλ
    else
      return t0, NaN, NaN, NaN, NaN
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTct,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              ϵ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTct,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   ϵ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn ::Int64)

  if lU < lr && nn < 1_000

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)

        # simulate
        stree, na, nn, lr =
          _sim_gbmct_it(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, ϵ, δt, srδt,
                     lr, lU, Iρi, na-1, nn, 1_000)

        if isnan(lr) || nn > 999
          return tree, na, nn, NaN
        end

        setproperty!(tree, :iμ, isextinct(stree))
        sete!(tree, e(tree) + e(stree))

        lλs = lλ(stree)

        if lastindex(lλs) === 2
          setfdt!(tree, fdt(tree) + fdt(stree))
        else
          setfdt!(tree, fdt(stree))
        end

        pop!(lλ0)
        popfirst!(lλs)
        append!(lλ0, lλs)

        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr =
        tip_sims!(tree.d1, t, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr =
        tip_sims!(tree.d2, t, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTct},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                ϵ    ::Float64,
                llc  ::Float64,
                dlλ   ::Float64,
                ssλ  ::Float64,
                Σλ   ::Float64,
                mc   ::Float64,
                th   ::Float64,
                δt   ::Float64,
                srδt ::Float64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTct},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     ϵ    ::Float64,
                     llc  ::Float64,
                     dlλ   ::Float64,
                     ssλ  ::Float64,
                     Σλ   ::Float64,
                     mc   ::Float64,
                     th   ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)
  @inbounds begin

    ξi   = Ξ[bix]
    bi   = idf[bix]
    i1   = d1(bi)
    i2   = d2(bi)
    ξ1   = Ξ[i1]
    root = iszero(pa(bi))

    # if crown root
    if root && iszero(e(bi))
      ξ2 = Ξ[i2]
      llc, dlλ, ssλ, Σλ, mc =
        _crown_update!(ξi, ξ1, ξ2, α, σλ, ϵ, llc, dlλ, ssλ, Σλ, mc, th,
          δt, srδt)
      setλt!(bi, lλ(ξi)[1])
    else
      # if stem branch
      if root
        llc, dlλ, ssλ, Σλ, mc =
          _stem_update!(ξi, α, σλ, ϵ, llc, dlλ, ssλ, Σλ, mc, th,
            δt, srδt)
      end

      # updates within the parent branch
      llc, dlλ, ssλ, Σλ =
        _update_gbm!(ξi, α, σλ, ϵ, llc, dlλ, ssλ, Σλ, δt, srδt, false)

      # get fixed tip
      lξi = fixtip(ξi)

      if iszero(i2)

        llc, ssλ, Σλ =
          update_duo_ϵ!(lλ(lξi), lλ(ξ1), e(lξi), e(ξ1), fdt(lξi), fdt(ξ1), 
            α, σλ, ϵ, llc, ssλ, Σλ, δt, srδt)

      # if internal branch
      else
        ξ2 = Ξ[i2]

        # make node update between decoupled trees
        llc, dlλ, ssλ, Σλ =
          update_triad_ϵ!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2),
            fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, ϵ, llc, dlλ, ssλ, Σλ, δt, srδt)

        # set fixed `λ(t)` in branch
        setλt!(bi, lλ(ξ1)[1])
      end
    end

    # carry on updates in the daughters
    llc, dlλ, ssλ, Σλ = 
      _update_gbm!(ξ1, α, σλ, ϵ, llc, dlλ, ssλ, Σλ, δt, srδt, 
        iszero(d1(idf[i1])))

    if i2 > 0
      ξ2 = Ξ[i2]
      llc, dlλ, ssλ, Σλ = 
        _update_gbm!(ξ2, α, σλ, ϵ, llc, dlλ, ssλ, Σλ, δt, srδt, 
          iszero(d1(idf[i2])))
    end
  end

  return llc, dlλ, ssλ, Σλ, mc
end




"""
    update_α_ϵ!(αc     ::Float64,
                λ0     ::Float64,
                σλ     ::Float64,
                ϵ      ::Float64,
                L      ::Float64,
                dlλ     ::Float64,
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
function update_α_ϵ!(αc     ::Float64,
                     λ0     ::Float64,
                     σλ     ::Float64,
                     ϵ      ::Float64,
                     L      ::Float64,
                     dlλ     ::Float64,
                     llc    ::Float64,
                     prc    ::Float64,
                     mc     ::Float64,
                     th     ::Float64,
                     crown  ::Int64,
                     δt     ::Float64,
                     srδt   ::Float64,
                     α_prior::NTuple{2,Float64})

  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2
  αp  = rnorm((dlλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  mp  = m_surv_gbmct(th, λ0, αp, σλ, ϵ, δt, srδt, 1_000, crown)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*dlλ*(αp - αc)/L) + llr
    prc += llrdnorm_x(αp, αc, ν, τ2)
    αc   = αp
    mc   = mp
  end

  return llc, prc, αc, mc
end




"""
    update_σ_ϵ!(σλc     ::Float64,
                λ0      ::Float64,
                α       ::Float64,
                ϵ       ::Float64,
                ssλ     ::Float64,
                n       ::Float64,
                llc     ::Float64,
                prc     ::Float64,
                mc      ::Float64,
                th      ::Float64,
                crown   ::Int64,
                δt      ::Float64,
                srδt    ::Float64,
                σλ_prior::NTuple{2,Float64})

Gibbs update for `σλ`.
"""
function update_σ_ϵ!(σλc     ::Float64,
                     λ0      ::Float64,
                     α       ::Float64,
                     ϵ       ::Float64,
                     ssλ     ::Float64,
                     n       ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     mc      ::Float64,
                     th      ::Float64,
                     crown   ::Int64,
                     δt      ::Float64,
                     srδt    ::Float64,
                     σλ_prior::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)
  σλp  = sqrt(σλp2)

  mp  = m_surv_gbmct(th, λ0, α, σλp, ϵ, δt, srδt, 1_000, crown)
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
    update_ϵ!(ϵc   ::Float64,
              λ0   ::Float64,
              α    ::Float64,
              σλ   ::Float64,
              llc  ::Float64,
              mc   ::Float64,
              th   ::Float64,
              crown::Int64,
              ϵtn  ::Float64,
              lac  ::Float64,
              ne   ::Float64,
              Σλ   ::Float64,
              δt   ::Float64,
              srδt ::Float64,
              ϵxpr ::Float64)

MCMC update for `ϵ` with acceptance log.
"""
function update_ϵ!(ϵc   ::Float64,
                   λ0   ::Float64,
                   α    ::Float64,
                   σλ   ::Float64,
                   llc  ::Float64,
                   mc   ::Float64,
                   th   ::Float64,
                   crown::Int64,
                   ϵtn  ::Float64,
                   lac  ::Float64,
                   ne   ::Float64,
                   Σλ   ::Float64,
                   δt   ::Float64,
                   srδt ::Float64,
                   ϵxpr ::Float64)

  ϵp  = abs(addupt(ϵc, ϵtn))::Float64
  llr = ne*log(ϵp/ϵc) + Σλ*(ϵc - ϵp)
  prr = ϵp > ϵxpr ? -Inf : 0.0

  # log probability
  lU = -randexp()

  # check if valid proposal before doing survival conditioning simulation
  if lU < llr + prr + log(1000.0/mc)

    mp   = m_surv_gbmct(th, λ0, α, σλ, ϵp, δt, srδt, 1_000, crown)
    llr += log(mp/mc)

    if lU < llr
      llc += llr
      ϵc   = ϵp
      mc   = mp
      lac += 1.0
    end
  end

  return llc, ϵc, mc, lac
end




"""
    update_ϵ!(ϵc   ::Float64,
              λ0   ::Float64,
              α    ::Float64,
              σλ   ::Float64,
              llc  ::Float64,
              mc   ::Float64,
              th   ::Float64,
              crown::Int64,
              ϵtn  ::Float64,
              ne   ::Float64,
              Σλ   ::Float64,
              δt   ::Float64,
              srδt ::Float64,
              ϵxpr ::Float64)

MCMC update for `ϵ`.
"""
function update_ϵ!(ϵc   ::Float64,
                   λ0   ::Float64,
                   α    ::Float64,
                   σλ   ::Float64,
                   llc  ::Float64,
                   mc   ::Float64,
                   th   ::Float64,
                   crown::Int64,
                   ϵtn  ::Float64,
                   ne   ::Float64,
                   Σλ   ::Float64,
                   δt   ::Float64,
                   srδt ::Float64,
                   ϵxpr ::Float64)

  ϵp  = abs(addupt(ϵc, ϵtn))::Float64
  llr = ne*log(ϵp/ϵc) + Σλ*(ϵc - ϵp)
  prr = ϵp > ϵxpr ? -Inf : 0.0

  # log probability
  lU = -randexp()

  # check if valid proposal before doing survival conditioning simulation
  if lU < llr + prr + log(1000.0/mc)

    mp   = m_surv_gbmct(th, λ0, α, σλ, ϵp, δt, srδt, 1_000, crown)
    llr += log(mp/mc)

    if lU < llr
      llc += llr
      ϵc   = ϵp
      mc   = mp
    end
  end

  return llc, ϵc, mc
end




