#=

Anagenetic `gbmbd` MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmbd(tree    ::sT_label,
                 X       ::Dict{String, Float64},
                 out_file::String;
                 λa_prior::NTuple{2,Float64}     = (0.0, 100.0),
                 μa_prior::NTuple{2,Float64}     = (0.0, 100.0),
                 α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                 σμ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                 x0_prior::NTuple{2,Float64}     = (0.0,  10.0),
                 βλ_prior::NTuple{2,Float64}     = (0.0,  10.0),
                 σx_prior::NTuple{2,Float64}     = (0.05, 0.05),
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 ϵi      ::Float64               = 0.2,
                 λi      ::Float64               = NaN,
                 μi      ::Float64               = NaN,
                 αi      ::Float64               = 0.0,
                 σλi     ::Float64               = 0.01,
                 σμi     ::Float64               = 0.01,
                 βλi     ::Float64               = 0.0,
                 pupdp   ::NTuple{4,Float64}     = (0.0, 0.1, 0.2, 0.2),
                 δt      ::Float64               = 1e-2,
                 prints  ::Int64                 = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `gbm-bd`.
"""
function insane_gbmbd(tree    ::sT_label,
                      X       ::Dict{String, Float64},
                      out_file::String;
                      λa_prior::NTuple{2,Float64}     = (0.0, 100.0),
                      μa_prior::NTuple{2,Float64}     = (0.0, 100.0),
                      α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                      σμ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                      x0_prior::NTuple{2,Float64}     = (0.0,  10.0),
                      βλ_prior::NTuple{2,Float64}     = (0.0,  10.0),
                      σx_prior::NTuple{2,Float64}     = (0.05, 0.05),
                      niter   ::Int64                 = 1_000,
                      nthin   ::Int64                 = 10,
                      nburn   ::Int64                 = 200,
                      ϵi      ::Float64               = 0.2,
                      λi      ::Float64               = NaN,
                      μi      ::Float64               = NaN,
                      αi      ::Float64               = 0.0,
                      σλi     ::Float64               = 0.01,
                      σμi     ::Float64               = 0.01,
                      βλi     ::Float64               = 0.0,
                      pupdp   ::NTuple{7,Float64}     = (0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                      δt      ::Float64               = 1e-2,
                      prints  ::Int64                 = 5,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n     = ntips(tree)
  th    = treeheight(tree)
  δt   *= max(0.1,round(th, RoundDown, digits = 2))
  srδt  = sqrt(δt)
  crown = Int64(iszero(e(tree)))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxc = make_idf(tree, tρ, X)

   # starting parameters (using method of moments)
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end
  mc = m_surv_gbmbd(th, log(λc), log(μc), αi, σλi, σμi, δt, srδt, 500, crown)

  # make a decoupled tree
  Ξ = make_Ξ(idf, xr, log(λc), log(μc), αi, σλi, σμi, σxc, δt, srδt, iTbdX)

  # set end of fix branch speciation times and
  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    setλt!(bi, lλ(Ξ[i])[end])
    if !it(bi)
      push!(inodes, i)
    end
  end

  # parameter updates (1: α, 2: σλ, 3: σμ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(7)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death gbm and trait evolution"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, σμc, βλc, σxc, mc =
    mcmc_burn_gbmbd(Ξ, idf, λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior,
      x0_prior, βλ_prior, σx_prior, nburn, αi, σλi, σμi, βλi, σxc, mc, th,
      crown, δt, srδt, inodes, pup, prints)

  # mcmc
  R, Ξv =
    mcmc_gbmbd(Ξ, idf, llc, prc, αc, σλc, σμc, βλc, σxc, mc, th, crown, 
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, x0_prior, βλ_prior, 
      σx_prior, niter, nthin, δt, srδt, inodes, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "mu_root"      => 2,
                 "alpha"        => 3,
                 "sigma_lambda" => 4,
                 "sigma_mu"     => 5,
                 "x0"           => 6,
                 "beta_lambda"  => 7,
                 "sigma_x"      => 8))

  write_ssr(R, pardic, out_file)

  return R, Ξv
end




"""
    mcmc_burn_gbmbd(Ξ       ::Vector{iTbdX},
                    idf     ::Vector{iBffs},
                    λa_prior::NTuple{2,Float64},
                    μa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
                    x0_prior::NTuple{2,Float64},
                    βλ_prior::NTuple{2,Float64},
                    σx_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    βλc     ::Float64,
                    σxc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    crown    ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

MCMC burn-in chain for `gbmbd`.
"""
function mcmc_burn_gbmbd(Ξ       ::Vector{iTbdX},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         μa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         σμ_prior::NTuple{2,Float64},
                         x0_prior::NTuple{2,Float64},
                         βλ_prior::NTuple{2,Float64},
                         σx_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         σμc     ::Float64,
                         βλc     ::Float64,
                         σxc     ::Float64,
                         mc      ::Float64,
                         th      ::Float64,
                         crown    ::Int64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  λ0  = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, σμc, βλc, σxc, δt, srδt) - Float64(crown)*λ0 +
        log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2,        σλ_prior[1], σλ_prior[2])   +
        logdinvgamma(σμc^2,        σμ_prior[1], σμ_prior[2])   +
        logdnorm(αc,               α_prior[1], α_prior[2]^2)   +
        logdunif(exp(λ0),          λa_prior[1], λa_prior[2])   +
        logdunif(exp(lμ(Ξ[1])[1]), μa_prior[1], μa_prior[2])   +
        logdnorm(βλc,              βλ_prior[1], βλ_prior[2]^2) +
        logdinvgamma(σxc^2,        σx_prior[1], σx_prior[2])

  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L                 = treelength(Ξ)      # tree length
  dλ                = deltaλ(Ξ)         # delta change in λ
  ssλ, ssμ, ssx, nx = sss_gbm(Ξ, αc, βλc)  #sum squares in λ and μ and X
  nin               = lastindex(inodes) # number of internal nodes
  el                = lastindex(idf)    # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for p in pup

      # update α
      if p === 1

        llc, prc, αc, mc  =
          update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, L, dλ, llc, prc,
            mc, th, crown, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nx = sss_gbm(Ξ, αc)

      # σλ & σμ update
      elseif p === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nx,
            llc, prc, mc, th, crown, δt, srδt, σλ_prior, σμ_prior)

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

      # gbm update
      elseif p === 6

        # nix = ceil(Int64,rand()*nin)
        # bix = inodes[nix]
        bix = 1

        llc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
            δt, srδt, lλxpr, lμxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, ssx, nx, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, βλc, σxc, llc, dλ,
            ssλ, ssμ, ssx, nx, L, δt, srδt)
      end
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, σμc, βλc, σxc, mc
end




"""
    mcmc_gbmbd(Ξ       ::Vector{iTbdX},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               σμc     ::Float64,
               βλc     ::Float64,
               σxc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               crown   ::Int64,
               λa_prior::NTuple{2,Float64},
               μa_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               σμ_prior::NTuple{2,Float64},
               x0_prior::NTuple{2,Float64},
               βλ_prior::NTuple{2,Float64},
               σx_prior::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Vector{Int64},
               prints  ::Int64)

MCMC chain for `gbmbd`.
"""
function mcmc_gbmbd(Ξ       ::Vector{iTbdX},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    βλc     ::Float64,
                    σxc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    crown   ::Int64,
                    λa_prior::NTuple{2,Float64},
                    μa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
                    x0_prior::NTuple{2,Float64},
                    βλ_prior::NTuple{2,Float64},
                    σx_prior::NTuple{2,Float64},
                    niter   ::Int64,
                    nthin   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Vector{Int64},
                    prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # crown or crown conditioning
  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L                 = treelength(Ξ)      # tree length
  dλ                = deltaλ(Ξ)         # delta change in λ
  ssλ, ssμ, ssx, nx = sss_gbm(Ξ, αc, βλc)  #sum squares in λ and μ and X
  nin               = lastindex(inodes) # number of internal nodes
  el                = lastindex(idf)    # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 11)

  # make Ξ vector
  Ξv = iTbdX[]

  # number of branches and of triads
  nbr  = lastindex(idf)

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    # parameter updates
    for p in pup

      # update α
      if p === 1

        llc, prc, αc, mc  =
          update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, L, dλ, llc, prc,
            mc, th, crown, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nx = sss_gbm(Ξ, αc)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, βλc, σxc, δt, srδt) - Float64(crown) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # σλ & σμ update
      elseif p === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nx,
            llc, prc, mc, th, crown, δt, srδt, σλ_prior, σμ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, βλc, σxc, δt, srδt) - Float64(crown) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # update beta
      elseif p === 3

        # to do!

      # update `x` diffusion rate
      elseif p === 4

        llc, prc, σxc = update_σx!(σxc, ssx, nx, llc, prc, σx_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, βλc, σxc, δt, srδt) - Float64(crown) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # update `x` bm
      elseif p === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ssx =
          update_x!(bix, Ξ, idf, σxc, llc, prc, ssx, δt, srδt, x0_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, βλc, σxc, δt, srδt) - Float64(crown) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # gbm update
      elseif p === 6

        # nix = ceil(Int64,rand()*nin)
        # bix = inodes[nix]
        bix = 1

        llc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
            δt, srδt, lλxpr, lμxpr)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, βλc, σxc, δt, srδt) - Float64(crown) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, ssx, nx, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, βλc, σxc, llc, dλ,
            ssλ, ssμ, ssx, nx, L, δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, βλc, σxc, δt, srδt) - Float64(crown) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end
      end
    end

    # log parameters
    lthin += 1
    if lthin === nthin
      lit += 1
      @inbounds begin
        R[lit,1]  = Float64(lit)
        R[lit,2]  = llc
        R[lit,3]  = prc
        R[lit,4]  = exp(lλ(Ξ[1])[1])
        R[lit,5]  = exp(lμ(Ξ[1])[1])
        R[lit,6]  = αc
        R[lit,7]  = σλc
        R[lit,8]  = σμc
        R[lit,9]  = xv(Ξ[1])[1]
        R[lit,10] = βλc
        R[lit,11] = σxc

        push!(Ξv, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ξv
end




"""
    update_fs!(bix ::Int64,
               Ξ   ::Vector{iTbdX},
               idf ::Vector{iBffs},
               α   ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               βλ  ::Float64,
               σx  ::Float64,
               llc ::Float64,
               dλ  ::Float64,
               ssλ ::Float64,
               ssμ ::Float64,
               nx  ::Float64,
               L   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Forward simulation proposal function for `gbmbd`.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{iTbdX},
                    idf ::Vector{iBffs},
                    α   ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    βλ  ::Float64,
                    σx  ::Float64,
                    llc ::Float64,
                    dλ  ::Float64,
                    ssλ ::Float64,
                    ssμ ::Float64,
                    ssx ::Float64,
                    nx  ::Float64,
                    L   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]

  # if terminal
  if it(bi)
    ξp, llr = fsbi_t(bi, ξc, α, σλ, σμ, βλ, σx, δt, srδt)
    drλ  = 0.0
    ssrλ = 0.0
    ssrμ = 0.0
    ssrx = 0.0
  # if internal
  else
    ξp, llr, drλ, ssrλ, ssrμ, ssrx =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, σμ, βλ, σx, δt, srδt)
  end

  # if accepted
  if isfinite(llr)
    ll1, dλ1, ssλ1, ssμ1, ssx1, nx1 =
      llik_gbm_ss(ξp, α, σλ, σμ, βλ, σx, δt, srδt)
    ll0, dλ0, ssλ0, ssμ0, ssx0, nx0 =
      llik_gbm_ss(ξc, α, σλ, σμ, βλ, σx, δt, srδt)

    # update llr, ssλ, nx, sns, ne, L,
    llc += ll1  - ll0  + llr
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    ssμ += ssμ1 - ssμ0 + ssrμ
    ssx += ssx1 - ssx0 + ssrx
    nx  += nx1  - nx0
    L   += treelength(ξp)   - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, dλ, ssλ, ssμ, ssx, nx, L
end




"""
    fsbi_t(bi  ::iBffs,
           ξc  ::iTbdX,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           βλ  ::Float64,
           σx  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTbdX,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
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
    _sim_gbmbd_t(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, xv(ξc)[1], βλ, σx,
      δt, srδt, lc, lU, Iρi, 0, 1, 500, xist, xfst, est)

  if na < 0 || isnan(llr)
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

    if isfinite(acr) && lU <  acr + llr

      # sample tip
      wti = sample(wp)
      if wti <= div(na,2)
        fixtip1!(t0, wti,      0, xc, σx, δt, srδt)
      else
        fixtip2!(t0, na-wti+1, 0, xc, σx, δt, srδt)
      end

      setni!(bi, na) # set new ni
      return t0, llr
    end
  # if unfix `x` node
  else
    if lU < llr
      _fixrtip!(t0, na)
      setni!(bi, na) # set new ni
      return t0, llr
    end
  end

  return t0, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTbdX,
           ξ1  ::iTbdX,
           ξ2  ::iTbdX,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           βλ  ::Float64,
           σx  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTbdX,
                ξ1  ::iTbdX,
                ξ2  ::iTbdX,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  λsp = Float64[]
  μsp = Float64[]
  xsp = Float64[]

  t0, nn =
    _sim_gbmbd_i(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, xv(ξc)[1], βλ, σx,
      δt, srδt, 1, 500, λsp, μsp, xsp)

  na = lastindex(λsp)

  if na < 1 || nn >= 500
    return t0, NaN, NaN, NaN, NaN, NaN
  end

  # get current speciation rates at branch time
  λsc = λst(bi)
  μsc = μst(bi)

  e1  = e(ξ1)
  sr1 = sqrt(e1)
  e2  = e(ξ2)
  sr2 = sqrt(e2)
  λ1c = lλ(ξ1)
  λ2c = lλ(ξ2)
  μ1c = lμ(ξ1)
  μ2c = lμ(ξ2)
  l1  = lastindex(λ1c)
  l2  = lastindex(λ2c)
  λ1  = λ1c[l1]
  λ2  = λ2c[l2]
  μ1  = μ1c[l1]
  μ2  = μ2c[l2]

  # proposed acceptance ratio
  wp = Float64[]
  ap = 0.0
  @simd for i in Base.OneTo(lastindex(λsp))
    λi = λsp[i]
    μi = μsp[i]
    wi  = exp(λi) * dnorm_bm(λi, λ1 - α*e1, sr1*σλ) *
                    dnorm_bm(λi, λ2 - α*e2, sr2*σλ) *
                    dnorm_bm(μi, μ1, sr1*σμ)        *
                    dnorm_bm(μi, μ2, sr2*σμ)
    ap += wi
    push!(wp, wi)
  end
  ap = log(ap)

  if isinf(ap)
    return t0, NaN, NaN, NaN, NaN, NaN
  end

  # current acceptance ratio
  ac = 0.0
  @simd for i in Base.OneTo(lastindex(λsc))
    λi = λsc[i]
    μi = μsc[i]
    ac += exp(λi) * dnorm_bm(λi, λ1 - α*e1, sr1*σλ) *
                    dnorm_bm(λi, λ2 - α*e2, sr2*σλ) *
                    dnorm_bm(μi, μ1, sr1*σμ)        *
                    dnorm_bm(μi, μ2, sr2*σμ)
  end
  ac = log(ac)

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
  μf  = μsp[wti]
  xp  = xsp[wti]

  llrd, acrd, drλ, ssrλ, ssrμ, ssrx, λ1p, λ2p, μ1p, μ2p, x1p, x2p =
    _daughters_update!(ξ1, ξ2, λf, μf, α, σλ, σμ, xp, βλ, σx, δt, srδt)

  acr += acrd

  if lU < acr

     # fix sampled tip
    if wti <= div(na,2)
      fixtip1!(t0, wti, 0)
    else
      fixtip2!(t0, na-wti+1, 0)
    end

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, βλ, σx, δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setni!( bi, na)                       # set new ni
      setλt!( bi, λf)                       # set new λt
      setλst!(bi, λsp)                      # set new λst
      setμst!(bi, μsp)                      # set new μst
      unsafe_copyto!(λ1c,    1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(λ2c,    1, λ2p, 1, l2) # set new daughter 2 λ vector
      unsafe_copyto!(μ1c,    1, μ1p, 1, l1) # set new daughter 1 μ vector
      unsafe_copyto!(μ2c,    1, μ2p, 1, l2) # set new daughter 2 μ vector
      unsafe_copyto!(xv(ξ1), 1, x1p, 1, l1) # set new daughter 1 x vector
      unsafe_copyto!(xv(ξ2), 1, x2p, 1, l2) # set new daughter 2 x vector

      return t0, llr, drλ, ssrλ, ssrμ, ssrx
    end
  end

  return t0, NaN, NaN, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTbdX,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              βλ  ::Float64,
              σx  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTbdX,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   βλ  ::Float64,
                   σx  ::Float64,
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
        xv0  = xv(tree)
        l    = lastindex(lλ0)

        # simulate
        stree, na, nn, lr =
          _sim_gbmbd_it(max(δt-fdti, 0.0), t, lλ0[l], lμ0[l], α, σλ, σμ,
            xv0[l], βλ, σx, δt, srδt, lr, lU, Iρi, na-1, nn, 500)

        if isnan(lr) || nn >= 500
          return tree, na, nn, NaN
        end

        setproperty!(tree, :iμ, isextinct(stree))
        sete!(tree, e(tree) + e(stree))

        lλs = lλ(stree)
        lμs = lμ(stree)
        xvs = xv(stree)

        if lastindex(lλs) === 2
          setfdt!(tree, fdt(tree) + fdt(stree))
        else
          setfdt!(tree, fdt(stree))
        end

        pop!(lλ0)
        pop!(lμ0)
        pop!(xv0)
        popfirst!(lλs)
        popfirst!(lμs)
        popfirst!(xvs)
        append!(lλ0, lλs)
        append!(lμ0, lμs)
        append!(xv0, xvs)

        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr =
        tip_sims!(tree.d1, t, α, σλ, σμ, βλ, σx, δt, srδt, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr =
        tip_sims!(tree.d2, t, α, σλ, σμ, βλ, σx, δt, srδt, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




