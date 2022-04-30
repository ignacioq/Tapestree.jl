#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#




"""
    insane_gbmpb(tree    ::sT_label,
                 X       ::Dict{String, Float64},
                 out_file::String;
                 δt      ::Float64               = 1e-2,
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 σλi     ::Float64               = 0.1,
                 αi      ::Float64               = 0.0,
                 βλi     ::Float64               = 0.0,
                 prints  ::Int64                 = 5,
                 pupdp   ::NTuple{4,Float64}     = (0.2, 0.2, 0.3, 0.3),
                 α_prior ::NTuple{2,Float64}     = (0.0,  10.0),
                 σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                 x0_prior::NTuple{2,Float64}     = (0.0,  10.0),
                 βλ_prior::NTuple{2,Float64}     = (0.0,  10.0),
                 σx_prior::NTuple{2,Float64}     = (0.05, 0.05),
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for GBM pure-birth.
"""
function insane_gbmpb(tree    ::sT_label,
                      X       ::Dict{String, Float64},
                      out_file::String;
                      δt      ::Float64               = 1e-2,
                      niter   ::Int64                 = 1_000,
                      nthin   ::Int64                 = 10,
                      nburn   ::Int64                 = 200,
                      σλi     ::Float64               = 0.1,
                      αi      ::Float64               = 0.0,
                      βλi     ::Float64               = 0.0,
                      prints  ::Int64                 = 5,
                      pupdp   ::NTuple{7,Float64}     = (0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                      α_prior ::NTuple{2,Float64}     = (0.0,  10.0),
                      σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                      x0_prior::NTuple{2,Float64}     = (0.0,  10.0),
                      βλ_prior::NTuple{2,Float64}     = (0.0,  10.0),
                      σx_prior::NTuple{2,Float64}     = (0.05, 0.05),
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # lλ root node
  lλa = log(λmle_cpb(tree))

  # make fix tree directory
  idf, xr, σxc = make_idf(tree, tρ, X)

  # make a decoupled tree
  Ξ = make_Ξ(idf, xr, lλa, αi, σλi, σxc, δt, srδt, iTpbX)

  # set end of fix branch speciation times and
  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    setλt!(bi, lλ(Ξ[i])[end])
    !it(bi) && push!(inodes, i)
  end

  # parameter updates (1: α, 2: σ, 3: gbm, 4: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(7)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running pure-birth gbm and trait evolution"

  # burn-in phase
  llc, prc, αc, σλc, βλc, σxc =
    mcmc_burn_gbmpb(Ξ, idf, α_prior, σλ_prior, x0_prior, βλ_prior, σx_prior, 
      nburn, αi, σλi, βλi, σxc, δt, srδt, inodes, pup, prints)

  # mcmc
  r, Ξv, αc, σλc, βλc, σxc = mcmc_gbmpb(Ξ, idf, llc, prc, αc, σλc, βλc, σxc,
    α_prior, σλ_prior, x0_prior, βλ_prior, σx_prior, niter, nthin, δt, srδt, 
    inodes, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "alpha"        => 2,
                 "sigma_lambda" => 3,
                 "x0"           => 4,
                 "beta_lambda"  => 5,
                 "sigma_x"      => 6))

  write_ssr(r, pardic, out_file)

  return r, Ξv
end





"""
    mcmc_burn_gbmpb(Ξ       ::Vector{iTpbX},
                    idf     ::Vector{iBffs},
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)


MCMC burn-in chain for GBM pure-birth.
"""
function mcmc_burn_gbmpb(Ξ       ::Vector{iTpbX},
                         idf     ::Vector{iBffs},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         x0_prior::NTuple{2,Float64},
                         βλ_prior::NTuple{2,Float64},
                         σx_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         βλc     ::Float64,
                         σxc     ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  nsi = iszero(e(Ξ[1])) ? 1.0 : 0.0

  # starting likelihood and prior
  llc = llik_gbm(Ξ, idf, αc, σλc, βλc, σxc, δt, srδt) - nsi*lλ(Ξ[1])[1] + 
        prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2]) +
        logdnorm(αc,  α_prior[1], α_prior[2]^2)       +
        logdnorm(βλc, βλ_prior[1], βλ_prior[2]^2)       +
        logdinvgamma(σxc^2, σx_prior[1], σx_prior[2])

  L            = treelength(Ξ)       # tree length
  dλ           = deltaλ(Ξ)           # delta change in λ
  ssλ, ssx, nx = sss_gbm(Ξ, αc, βλc) # sum squares in λ
  nin          = lastindex(inodes)   # number of internal nodes
  el           = lastindex(idf)      # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      ## parameter updates
      # update drift
      if p === 1

        llc, prc, αc = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, nx = sss_gbm(Ξ, αc)

      # update `λ` diffusion rate
      elseif p === 2

        llc, prc, σλc = update_σ!(σλc, ssλ, nx, llc, prc, σλ_prior)

      # update beta
      elseif p === 3

        # to do!

      # update `x` diffusion rate
      elseif p === 4

        llc, prc, σxc = update_σx!(σxc, ssx, nx, llc, prc, σx_prior)

      # update `x` bm
      elseif p === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ssx =
          update_x!(bix, Ξ, idf, σxc, llc, prc, ssx, δt, srδt, x0_prior)

      # update gbm
      elseif p === 6

        # nix = ceil(Int64,rand()*nin)
        # bix = inodes[nix]
        bix = 1

        llc, dλ, ssλ =
          update_gbm!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, δt, srδt)

      # forward simulation
      else
        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssx, nx, L =
          update_fs!(bix, Ξ, idf, αc, σλc, βλc, σxc, llc, dλ, 
            ssλ, ssx, nx, L, δt, srδt)

      end
    end

    next!(pbar)
  end

  return llc, prc, αc, σλc, βλc, σxc
end




"""
    mcmc_gbmpb(Ξ       ::Vector{iTpbX},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               λ0_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Array{Int64,1},
               prints  ::Int64)

MCMC chain for GBM pure-birth.
"""
function mcmc_gbmpb(Ξ       ::Vector{iTpbX},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    βλc     ::Float64,
                    σxc     ::Float64,
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    x0_prior::NTuple{2,Float64},
                    βλ_prior::NTuple{2,Float64},
                    σx_prior::NTuple{2,Float64},
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

  r = Array{Float64,2}(undef, nlogs, 9)

  # make Ξ vector
  Ξv = iTpbX[]

  L            = treelength(Ξ)       # tree length
  dλ           = deltaλ(Ξ)           # delta change in λ
  ssλ, ssx, nx = sss_gbm(Ξ, αc, βλc) # sum squares in λ
  nin          = lastindex(inodes)   # number of internal nodes
  el           = lastindex(idf)      # number of branches

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      ## parameter updates
      # update drift
      if p === 1

        llc, prc, αc = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, nx = sss_gbm(Ξ, αc)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, βλc, σxc, δt, srδt) - lλ(Ξ[1])[1] + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    t ll0, llc, it, p
        #    return
        # end

      # update diffusion rate
      elseif p === 2

        llc, prc, σλc = update_σ!(σλc, ssλ, nx, llc, prc, σλ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, βλc, σxc, δt, srδt) - lλ(Ξ[1])[1] + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, it, p
        #    return
        # end

      # update beta
      elseif p === 3

        # to do!

      # update `x` diffusion rate
      elseif p === 4

        llc, prc, σxc = update_σx!(σxc, ssx, nx, llc, prc, σx_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, βλc, σxc, δt, srδt) - lλ(Ξ[1])[1] + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, it, p
        #    return
        # end

      # update `x` bm
      elseif p === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ssx =
          update_x!(bix, Ξ, idf, σxc, llc, prc, ssx, δt, srδt, x0_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, βλc, σxc, δt, srδt) - lλ(Ξ[1])[1] + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, it, p
        #    return
        # end

      # update gbm
      elseif p === 6

        # nix = ceil(Int64,rand()*nin)
        # bix = inodes[nix]
        bix = 1

        llc, dλ, ssλ =
          update_gbm!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, βλc, σxc, δt, srδt) - lλ(Ξ[1])[1] + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, it, p
        #    return
        # end

      # update by forward simulation
      else
        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssx, nx, L =
          update_fs!(bix, Ξ, idf, αc, σλc, βλc, σxc, llc, dλ, 
            ssλ, ssx, nx, L, δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, βλc, σxc, δt, srδt) - lλ(Ξ[1])[1] + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, it, p
        #    return
        # end
      end
    end

    # log parameters
    lthin += 1
    if lthin === nthin
      lit += 1
      @inbounds begin
        r[lit,1] = Float64(lit)
        r[lit,2] = llc
        r[lit,3] = prc
        r[lit,4] = exp(lλ(Ξ[1])[1])
        r[lit,5] = αc
        r[lit,6] = σλc
        r[lit,7] = xv(Ξ[1])[1]
        r[lit,8] = βλc
        r[lit,9] = σxc
        push!(Ξv, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return r, Ξv, αc, σλc, βλc, σxc
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{iTpbX},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               βλ   ::Float64,
               σx   ::Float64,
               llc  ::Float64,
               dλ   ::Float64,
               ssλ  ::Float64,
               nx   ::Float64,
               L    ::Float64,
               δt   ::Float64,
               srδt ::Float64)

Forward simulation proposal function for constant birth-death.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{iTpbX},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    βλ   ::Float64,
                    σx   ::Float64,
                    llc  ::Float64,
                    dλ   ::Float64,
                    ssλ  ::Float64,
                    ssx  ::Float64,
                    nx   ::Float64,
                    L    ::Float64,
                    δt   ::Float64,
                    srδt ::Float64)

  bi = idf[bix]
  ξc = Ξ[bix]

  # if terminal
  if it(bi)
    ξp, llr = fsbi_t(bi, ξc, α, σλ, βλ, σx, δt, srδt)
    drλ  = 0.0
    ssrλ = 0.0
    ssrx = 0.0
  # if internal
  else
    ξp, llr, drλ, ssrλ, ssrx =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, βλ, σx, δt, srδt)
  end

  # if accepted
  if isfinite(llr)
    ll1, dλ1, ssλ1, ssx1, nx1 = llik_gbm_ss(ξp, α, σλ, βλ, σx, δt, srδt)
    ll0, dλ0, ssλ0, ssx0, nx0 = llik_gbm_ss(ξc, α, σλ, βλ, σx, δt, srδt)

    # update llr, ssλ, nx, L
    llc += ll1  - ll0 + llr
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    ssx += ssx1 - ssx0 + ssrx
    nx  += nx1  - nx0
    L   += treelength(ξp) - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, dλ, ssλ, ssx, nx, L
end




"""
    fsbi_t(bi  ::iBffs,
           ξc  ::iTpbX,
           α   ::Float64,
           σλ  ::Float64,
           βλ  ::Float64,
           σx  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTpbX,
                α   ::Float64,
                σλ  ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  xist = Float64[]
  xfst = Float64[]
  est  = Float64[]

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_gbmpb_t(e(bi), lλ(ξc)[1], α, σλ, xv(ξc)[1], βλ, σx, δt, srδt, lc, lU, Iρi, 
      0, 1, 500, xist, xfst, est)

  if isnan(llr)
    return t0, NaN
  end

  # if fix `x` node
  if fx(bi)
    #get fix `x` and edge
    xc  = fixed_xt(ξc)

    # sample tip according to fix `x` value
    acr = 0.0
    wp = Float64[]
    @simd for i in Base.OneTo(na)
      srt = sqrt(est[i])
      wi  = dnorm_bm(xfst[i], xist[i], srt*σx)/dnorm_bm(xc, xist[i], srt*σx)
      push!(wp, wi)
      acr += wi
    end
    acr = log(acr)

    # sample tip
    wti = sample(wp)

    if wti <= div(na,2)
      fixtip1!(t0, wti,      0, xc, σx, δt, srδt)
    else
      fixtip2!(t0, na-wti+1, 0, xc, σx, δt, srδt)
    end
  # if unfix `x` node
  else
    _fixrtip!(t0, na)
    acr = 0.0
  end

  if isfinite(acr) && lU <  acr + llr
    setni!(bi, na) # set new ni
    return t0, llr
  end

  return t0, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTpbX,
           ξ1  ::iTpbX,
           ξ2  ::iTpbX,
           α   ::Float64,
           σλ  ::Float64,
           βλ  ::Float64,
           σx  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTpbX,
                ξ1  ::iTpbX,
                ξ2  ::iTpbX,
                α   ::Float64,
                σλ  ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  λsp = Float64[]
  xsp = Float64[]

  # forward simulation during branch length
  t0, na = 
    _sim_gbmpb_i(e(bi), lλ(ξc)[1], α, σλ, xv(ξc)[1], βλ, σx, δt, srδt, 
      1, 500, λsp, xsp)

  if na >= 500
    return t0, NaN, NaN, NaN, NaN
  end

  # get current speciation rates at branch time
  λsc = λst(bi)

  e1  = e(ξ1)
  sr1 = sqrt(e1)
  e2  = e(ξ2)
  sr2 = sqrt(e2)
  λ1c = lλ(ξ1)
  λ2c = lλ(ξ2)
  l1  = lastindex(λ1c)
  l2  = lastindex(λ2c)
  λ1  = λ1c[l1]
  λ2  = λ2c[l2]

  # current acceptance ratio
  ac = 0.0
  for λi in λsc
    ac += exp(λi) * dnorm_bm(λi, λ1 - α*e1, sr1*σλ) *
                    dnorm_bm(λi, λ2 - α*e2, sr2*σλ)
  end
  ac = log(ac)

  # proposed acceptance ratio
  wp = Float64[]
  ap = 0.0
  for λi in λsp
    wi  = exp(λi) * dnorm_bm(λi, λ1 - α*e1, sr1*σλ) *
                    dnorm_bm(λi, λ2 - α*e2, sr2*σλ)
    ap += wi
    push!(wp, wi)
  end
  ap = log(ap)

  if isinf(ap)
    return t0, NaN, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = ap - ac

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample tip
  wti = sample(wp)
  λf  = λsp[wti]
  xp  = xsp[wti]

  llrd, acrd, drλ, ssrλ, ssrx, λ1p, λ2p, x1p, x2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, xp, βλ, σx, δt, srδt)

  acr += acrd

  if lU < acr

    if wti <= div(na,2)
      fixtip1!(t0, wti, 0)
    else
      fixtip2!(t0, na - wti + 1, 0)
    end

    # simulated remaining tips until the present
    if na > 1
      t0, na, acr =
        tip_sims!(t0, tf(bi), α, σλ, βλ, σx, δt, srδt, acr, lU, Iρi, na)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setni!( bi, na)                       # set new ni
      setλt!( bi, λf)                       # set new λt
      setλst!(bi, λsp)                      # set new λst
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(λ1c, 1, λ2p, 1, l2) # set new daughter 2 λ vector
      unsafe_copyto!(xv(ξ1), 1, x1p, 1, l1) # set new daughter 1 x vector
      unsafe_copyto!(xv(ξ2), 1, x2p, 1, l2) # set new daughter 2 x vector

      return t0, llr, drλ, ssrλ, ssrx
    else
      return t0, NaN, NaN, NaN, NaN
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTpbX,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              βλ  ::Float64,
              σx  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTpbX,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   βλ  ::Float64,
                   σx  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64)

 if lU < lr && na < 1_000

    if istip(tree)
      if !isfix(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)
        xv0  = xv(tree)

        # simulate
        stree, na, lr =
          _sim_gbmpb_it(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, 
            βλ, xv0[end], σx, δt, srδt, lr, lU, Iρi, na, 500)

        if isnan(lr) || na >= 500
          return tree, na, NaN
        end

        sete!( tree, e(tree) + e(stree))

        lλs = lλ(stree)
        xvs = xv(stree)

        if lastindex(lλs) === 2
          setfdt!(tree, fdt(tree) + fdt(stree))
        else
          setfdt!(tree, fdt(stree))
        end

        pop!(lλ0)
        pop!(xv0)
        popfirst!(lλs)
        popfirst!(xvs)
        append!(lλ0, lλs)
        append!(xv0, xvs)

        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, lr = 
        tip_sims!(tree.d1, t, α, σλ, βλ, σx, δt, srδt, lr, lU, Iρi, na)
      tree.d2, na, lr = 
        tip_sims!(tree.d2, t, α, σλ, βλ, σx, δt, srδt, lr, lU, Iρi, na)
    end

    return tree, na, lr
  end

  return tree, na, NaN
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTpbX},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                llc  ::Float64,
                dλ   ::Float64,
                ssλ  ::Float64,
                δt   ::Float64,
                srδt ::Float64)

Make a `gbm` update for an interna branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTpbX},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  ξ1   = Ξ[d1(bi)]
  ξ2   = Ξ[d2(bi)]
  ter1 = it(idf[d1(bi)])
  ter2 = it(idf[d2(bi)])

  # if crown root
  if iszero(pa(bi)) && iszero(e(ξi))
    llc, dλ, ssλ =
      _crown_update!(ξi, ξ1, ξ2, α, σλ, llc, dλ, ssλ, δt, srδt)
    setλt!(bi, lλ(ξi)[1])
  else
    # if stem
    if iszero(pa(bi))
      llc, dλ, ssλ = _stem_update!(ξi, α, σλ, llc, dλ, ssλ, δt, srδt)
    end

    # # updates within the parent branch
    # llc, dλ, ssλ = _update_gbm!(ξi, α, σλ, llc, dλ, ssλ, δt, srδt, false)

    # # get fixed tip
    # lξi = fixtip(ξi)

    # # make between decoupled trees node update
    # llc, dλ, ssλ = update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2),
    #   fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, llc, dλ, ssλ, δt, srδt)

    # # set fixed `λ(t)` in branch
    # setλt!(bi, lλ(lξi)[end])
  end

  # # carry on updates in the daughters
  # llc, dλ, ssλ = _update_gbm!(ξ1, α, σλ, llc, dλ, ssλ, δt, srδt, ter1)
  # llc, dλ, ssλ = _update_gbm!(ξ2, α, σλ, llc, dλ, ssλ, δt, srδt, ter2)

  return llc, dλ, ssλ
end




"""
    update_x!(bix     ::Int64,
              Ξ       ::Vector{T},
              idf     ::Vector{iBffs},
              σx      ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              ssx     ::Float64,
              δt      ::Float64,
              srδt    ::Float64,
              x0_prior::NTuple{2, Float64}) where {T <: iTX}

Make a `gbm` update for an internal branch and its descendants.
"""
function update_x!(bix     ::Int64,
                   Ξ       ::Vector{T},
                   idf     ::Vector{iBffs},
                   σx      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   ssx     ::Float64,
                   δt      ::Float64,
                   srδt    ::Float64,
                   x0_prior::NTuple{2, Float64}) where {T <: iTX}

  ξi  = Ξ[bix]
  bi  = idf[bix]
  id1 = d1(bi)
  id2 = d2(bi)
  ξ1  = Ξ[id1]
  ξ2  = Ξ[id2]

  root = iszero(pa(bi))
  # if cron root
  if root && iszero(e(bi))
    if !fx(bi)
      llc, prc, ssx =
         _crown_update_x!(ξi, ξ1, ξ2, σx, llc, prc, ssx, δt, srδt, x0_prior)
    end
  else
    # if stem
    if root
      if !fx(bi)
        llc, prc, ssx = 
          _stem_update_x!(ξi, σx, llc, prc, ssx, δt, srδt, x0_prior)
      end
    end

    # updates within the parent branch
    llc, ssx = _update_x!(ξi, σx, llc, ssx, δt, srδt, false)

    # get fixed tip
    lξi = fixtip(ξi)

    # make between decoupled trees node update
    if !fx(bi)
      llc, ssx = _update_triad_x!(lξi, ξ1, ξ2, σx, llc, ssx, δt, srδt)
    end
  end

  # carry on updates in the daughters
  llc, ssx = 
    _update_x!(ξ1, σx, llc, ssx, δt, srδt, !fx(idf[id1]) && it(idf[id1]))
  llc, ssx = 
    _update_x!(ξ2, σx, llc, ssx, δt, srδt, !fx(idf[id2]) && it(idf[id2]))

  return llc, prc, ssx
end




"""
    _stem_update_x!(ξi      ::T,
                    σx      ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    ssx     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    x0_prior::NTuple{2,Float64}) where {T <: iTX}

Make stem update for trait.
"""
function _stem_update_x!(ξi      ::T,
                         σx      ::Float64,
                         llc     ::Float64,
                         prc     ::Float64,
                         ssx     ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         x0_prior::NTuple{2,Float64}) where {T <: iTX}

  m0, σx0 = x0_prior
  σx02    = σx0^2
  xvi     = xv(ξi)
  li      = lastindex(xvi)
  xo      = xvi[1]
  m       = xvi[li]
  el      = e(ξi)
  s2      = el * σx^2

  # gibbs sampling
  xn = rnorm((m * σx02 + m0 * s2) / (s2 + σx02), sqrt(s2 * σx02 / (s2 + σx02)))

  xvp  = Vector{Float64}(undef, li)
  fdti = fdt(ξi)

  bb!(xvp, xn, m, σx, δt, fdti, srδt)

  llr, ssrx = llr_ssx(xvp, xvi, σx, δt, fdti, srδt)

  unsafe_copyto!(xvi, 1, xvp, 1, li)

  # update llc, prc and ssx
  llc += llr
  prc += llrdnorm_x(xn, xo, m0, σx02)
  ssx += ssrx

  return llc, prc, ssx
end




"""
    _crown_update_x!(ξi      ::T,
                     ξ1      ::T,
                     ξ2      ::T,
                     σx      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ssx     ::Float64,
                     δt      ::Float64,
                     srδt    ::Float64,
                     x0_prior::NTuple{2,Float64})  where {T <: iTX}

Make crown update for trait.
"""
function _crown_update_x!(ξi      ::T,
                          ξ1      ::T,
                          ξ2      ::T,
                          σx      ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          ssx     ::Float64,
                          δt      ::Float64,
                          srδt    ::Float64,
                          x0_prior::NTuple{2,Float64}) where {T <: iTX}

  m0, σx0 = x0_prior
  σx02    = σx0^2

  xc1     = xv(ξ1)
  xc2     = xv(ξ2)
  l1      = lastindex(xc1)
  l2      = lastindex(xc2)
  xo      = xc1[1]
  x1      = xc1[l1]
  x2      = xc2[l2]
  e1      = e(ξ1)
  e2      = e(ξ2)

  inve   = 1.0/(e1 + e2)
  m      = (e2 * inve * x1 + e1 * inve * x2)
  sigma2 = e1 * e2 * inve * σx^2

  # gibbs sampling
  xn = rnorm((m * σx02 + m0 * sigma2) / (sigma2 + σx02),
             sqrt(sigma2 * σx02 / (sigma2 + σx02)))

  xc1  = xv(ξ1)
  xc2  = xv(ξ2)
  xp1  = Vector{Float64}(undef, lastindex(xc1))
  xp2  = Vector{Float64}(undef, lastindex(xc2))
  fdt1 = fdt(ξ1)
  fdt2 = fdt(ξ2)

  fill!(xv(ξi), xn)
  bb!(xp1, xn, x1, σx, δt, fdt1, srδt)
  bb!(xp2, xn, x2, σx, δt, fdt2, srδt)

  llr1, ssrx1 = llr_ssx(xp1, xc1, σx, δt, fdt1, srδt)
  llr2, ssrx2 = llr_ssx(xp2, xc2, σx, δt, fdt2, srδt)

  unsafe_copyto!(xc1, 1, xp1, 1, l1)
  unsafe_copyto!(xc2, 1, xp2, 1, l2)

  # update llc, prc and ssx
  llc += llr1 + llr2
  prc += llrdnorm_x(xn, xo, m0, σx02)
  ssx += ssrx1 + ssrx2

  return llc, prc, ssx
end




"""
    _update_x!(tree::T,
               σx  ::Float64,
               llc ::Float64,
               ssx ::Float64,
               δt  ::Float64,
               srδt::Float64,
               ufx ::Bool) where {T <: iTX}

Do gbm updates on a decoupled tree recursively.
"""
function _update_x!(tree::T,
                    σx  ::Float64,
                    llc ::Float64,
                    ssx ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    ufx ::Bool) where {T <: iTX}

  if def1(tree)
    llc, ssx = _update_triad_x!(tree, tree.d1, tree.d2, σx, llc, ssx, δt, srδt)
    llc, ssx = _update_x!(tree.d1, σx, llc, ssx, δt, srδt, ufx)
    llc, ssx = _update_x!(tree.d2, σx, llc, ssx, δt, srδt, ufx)
  elseif isfix(tree)
    llc, ssx = _update_tip_x!(tree, σx, llc, ssx, δt, srδt, ufx)
  else
    llc, ssx = _update_tip_x!(tree, σx, llc, ssx, δt, srδt, true)
  end

  return llc, ssx
end




"""
    _update_triad_x!(ξi  ::T,
                     ξ1  ::T,
                     ξ2  ::T,
                     σx  ::Float64,
                     llc ::Float64,
                     ssx ::Float64
                     δt  ::Float64,
                     srδt::Float64) where {T <: iTX}

Make gibbs node update for trait.
"""
function _update_triad_x!(ξi  ::T,
                          ξ1  ::T,
                          ξ2  ::T,
                          σx  ::Float64,
                          llc ::Float64,
                          ssx ::Float64,
                          δt  ::Float64,
                          srδt::Float64) where {T <: iTX}

  xca  = xv(ξi)
  xc1  = xv(ξ1)
  xc2  = xv(ξ2)
  la   = lastindex(xca)
  l1   = lastindex(xc1)
  l2   = lastindex(xc2)
  xa   = xca[1]
  xo   = xc1[1]
  x1   = xc1[l1]
  x2   = xc2[l2]
  xpa  = Vector{Float64}(undef, la)
  xp1  = Vector{Float64}(undef, l1)
  xp2  = Vector{Float64}(undef, l2)
  fdta = fdt(ξi)
  fdt1 = fdt(ξ1)
  fdt2 = fdt(ξ2)
  ea   = e(ξi)
  e1   = e(ξ1)
  e2   = e(ξ2)

  # gibbs sampling
  if iszero(ea) || iszero(e1) || iszero(e2)
    xn = xc1[1]
  else
    xn = trioprop(xa, x1, x2, ea, e1, e2, σx)
  end

  bb!(xpa, xa, xn, σx, δt, fdta, srδt)
  bb!(xp1, xn, x1, σx, δt, fdt1, srδt)
  bb!(xp2, xn, x2, σx, δt, fdt2, srδt)

  llra, ssrxa = llr_ssx(xpa, xca, σx, δt, fdta, srδt)
  llr1, ssrx1 = llr_ssx(xp1, xc1, σx, δt, fdt1, srδt)
  llr2, ssrx2 = llr_ssx(xp2, xc2, σx, δt, fdt2, srδt)

  unsafe_copyto!(xca, 1, xpa, 1, la)
  unsafe_copyto!(xc1, 1, xp1, 1, l1)
  unsafe_copyto!(xc2, 1, xp2, 1, l2)

  # update llc, prc and ssx
  llc += llra + llr1 + llr2
  ssx += ssrxa + ssrx1 + ssrx2

  return llc, ssx
end




"""
    _update_tip_x!(ξi  ::T,
                   σx  ::Float64,
                   llc ::Float64,
                   ssx ::Float64
                   δt  ::Float64,
                   srδt::Float64,
                   ufx ::Bool) where {T <: iTX}

Make gibbs node update for trait.
"""
function _update_tip_x!(ξi  ::T,
                        σx  ::Float64,
                        llc ::Float64,
                        ssx ::Float64,
                        δt  ::Float64,
                        srδt::Float64,
                        ufx ::Bool) where {T <: iTX}

  xvi = xv(ξi)
  li  = lastindex(xvi)
  xa  = xvi[1]
  ea  = e(ξi)

  # gibbs sampling
  if ufx
    xn = rnorm(xa, sqrt(ea)*σx)
  else
    xn = xvi[li]
  end

  xvp  = Vector{Float64}(undef, li)
  fdti = fdt(ξi)

  bb!(xvp, xa, xn, σx, δt, fdti, srδt)

  llr, ssrx = llr_ssx(xvp, xvi, σx, δt, fdti, srδt)

  unsafe_copyto!(xvi, 1, xvp, 1, li)

  # update llc, prc and ssx
  llc += llr
  ssx += ssrx

  return llc, ssx
end



