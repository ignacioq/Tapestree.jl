#=

Anagenetic `gbmfbd` MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmfbd(tree    ::sTf_label,
                  out_file::String;
                  λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                  μa_prior::NTuple{2,Float64} = (0.0, 100.0),
                  α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                  σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                  σμ_prior::NTuple{2,Float64} = (0.05, 0.05),
                  ψ_prior ::NTuple{2,Float64} = (1.0, 1.0),
                  niter   ::Int64             = 1_000,
                  nthin   ::Int64             = 10,
                  nburn   ::Int64             = 200,
                  ϵi      ::Float64           = 0.2,
                  λi      ::Float64           = NaN,
                  μi      ::Float64           = NaN,
                  ψi      ::Float64           = NaN,
                  αi      ::Float64           = 0.0,
                  σλi     ::Float64           = 0.01,
                  σμi     ::Float64           = 0.01,
                  pupdp   ::NTuple{5,Float64} = (0.0, 0.1, 0.1, 0.2, 0.2),
                  δt      ::Float64           = 1e-2,
                  prints  ::Int64             = 5,
                  tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `gbm-bd`.
"""
function insane_gbmfbd(tree    ::sTf_label,
                       X       ::Dict{String, Float64},
                       out_file::String;
                       λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                       μa_prior::NTuple{2,Float64} = (0.0, 100.0),
                       α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                       σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                       σμ_prior::NTuple{2,Float64} = (0.05, 0.05),
                       ψ_prior ::NTuple{2,Float64} = (1.0, 1.0),
                       x0_prior::NTuple{2,Float64} = (0.0,  10.0),
                       βλ_prior::NTuple{2,Float64} = (0.0,  10.0),
                       σx_prior::NTuple{2,Float64} = (0.05, 0.05),
                       niter   ::Int64             = 1_000,
                       nthin   ::Int64             = 10,
                       nburn   ::Int64             = 200,
                       ϵi      ::Float64           = 0.2,
                       λi      ::Float64           = NaN,
                       μi      ::Float64           = NaN,
                       ψi      ::Float64           = NaN,
                       αi      ::Float64           = 0.0,
                       σλi     ::Float64           = 0.01,
                       σμi     ::Float64           = 0.01,
                       βλi     ::Float64           = 0.0,
                       pupdp   ::NTuple{8,Float64} = (0.1, 0.1, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1),
                       δt      ::Float64           = 1e-2,
                       prints  ::Int64             = 5,
                       tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1,round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxc = make_idf(tree, tρ, X)

  # starting parameters
  if isnan(λi) && isnan(μi) && isnan(ψi)
    # if only one tip
    if isone(n)
      λc = prod(λ_prior)
      μc = prod(μ_prior)
    else
      λc, μc = moments(Float64(n), th, ϵi)
    end
    # if no sampled fossil
    if iszero(nfossils(tree))
      ψc = prod(ψ_prior)
    else
      ψc = Float64(nfossils(tree))/Float64(treelength(tree))
    end
  else
    λc, μc, ψc = λi, μi, ψi
  end

  # define conditioning
  if ntipsalive(tree) > 0
    # if crown conditioning
    if def1(tree) && def2(tree) &&
       ntipsalive(tree.d1) > 0 && ntipsalive(tree.d2) > 0
      crown = 1
    # if crown conditioning
    else
      crown = 0
    end
  # no survival
  else
    crown = 2
  end
  mc = m_surv_gbmbd(th, log(λc), log(μc), αi, σλi, σμi, δt, srδt, 500, crown)

  # make a decoupled tree
  Ξ = make_Ξ(idf, xr, log(λc), log(μc), αi, σλi, σμi, σxc, δt, srδt, iTfbdX)

  # set end of fix branch speciation times and
  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    setλt!(bi, lλ(Ξ[i])[end])
    if !it(idf[i]) || isfossil(idf[i])
      push!(inodes, i)
    end
  end

  # parameter updates (1: α, 2: σλ & σμ, 3: ψ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(8)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running fossilized birth-death gbm"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, σμc, ψc, βλc, σxc, mc  =
    mcmc_burn_gbmbd(Ξ, idf, λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, 
      ψ_prior, x0_prior, βλ_prior, σx_prior, nburn, αi, σλi, σμi, ψc, βλi, σxc,
      mc, th, crown, δt, srδt, inodes, pup, prints)

  # mcmc
  R, Ξv =
    mcmc_gbmbd(Ξ, idf, llc, prc, αc, σλc, σμc, ψc, βλc, σxc, mc, th, crown,
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, ψ_prior, x0_prior, 
      βλ_prior, σx_prior, niter, nthin, δt, srδt, inodes, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "mu_root"      => 2,
                 "alpha"        => 3,
                 "sigma_lambda" => 4,
                 "sigma_mu"     => 5,
                 "psi"          => 6,
                 "x0"           => 7,
                 "beta_lambda"  => 8,
                 "sigma_x"      => 9))

  write_ssr(R, pardic, out_file)

  return R, Ξv
end




"""
    mcmc_burn_gbmbd(Ξ       ::Vector{iTfbdX},
                    idf     ::Vector{iBffs},
                    λa_prior::NTuple{2,Float64},
                    μa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc     ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    crown    ::Bool,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

MCMC burn-in chain for `gbmbd`.
"""
function mcmc_burn_gbmbd(Ξ       ::Vector{iTfbdX},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         μa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         σμ_prior::NTuple{2,Float64},
                         ψ_prior ::NTuple{2,Float64},
                         x0_prior::NTuple{2,Float64},
                         βλ_prior::NTuple{2,Float64},
                         σx_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         σμc     ::Float64,
                         ψc      ::Float64,
                         βλc     ::Float64,
                         σxc     ::Float64,
                         mc      ::Float64,
                         th      ::Float64,
                         crown   ::Int64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  λ0  = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, δt, srδt) - crown*λ0 +
        log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2,        σλ_prior[1], σλ_prior[2])   +
        logdinvgamma(σμc^2,        σμ_prior[1], σμ_prior[2])   +
        logdnorm(αc,               α_prior[1],  α_prior[2]^2)  +
        logdunif(exp(λ0),          λa_prior[1], λa_prior[2])   +
        logdunif(exp(lμ(Ξ[1])[1]), μa_prior[1], μa_prior[2])   +
        logdgamma(ψc,              ψ_prior[1],  ψ_prior[2])    +
        logdnorm(βλc,              βλ_prior[1], βλ_prior[2]^2) +
        logdinvgamma(σxc^2,        σx_prior[1], σx_prior[2])

  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L                 = treelength(Ξ)      # tree length
  nf                = Float64(nfossils(Ξ)) # number of fossilization events
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

      # psi update
      elseif p === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # update beta
      elseif p === 4

        # to do!

      # update `x` diffusion rate
      elseif p === 5

        llc, prc, σxc = update_σx!(σxc, ssx, nx, llc, prc, σx_prior)

      # update `x` bm
      elseif p === 6

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ssx =
          update_x!(bix, Ξ, idf, σxc, llc, prc, ssx, δt, srδt, x0_prior)

      # gbm update
      elseif p === 7

        # nix = ceil(Int64,rand()*nin)
        # bix = inodes[nix]
        bix = 1

        llc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
            crown, δt, srδt, lλxpr, lμxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, ssx, nx, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, llc, dλ, 
            ssλ, ssμ, ssx, nx, L, δt, srδt)

      end
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, σμc, ψc, βλc, σxc, mc
end




"""
    mcmc_gbmbd(Ξ       ::Vector{iTfbdX},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               σμc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               crown    ::Bool,
               λa_prior::NTuple{2,Float64},
               μa_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               σμ_prior::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Vector{Int64},
               prints  ::Int64)

MCMC chain for `gbmbd`.
"""
function mcmc_gbmbd(Ξ       ::Vector{iTfbdX},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    ψc      ::Float64,
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
                    ψ_prior ::NTuple{2,Float64},
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
  nf                = Float64(nfossils(Ξ)) # number of fossilization events
  dλ                = deltaλ(Ξ)         # delta change in λ
  ssλ, ssμ, ssx, nx = sss_gbm(Ξ, αc, βλc)  #sum squares in λ and μ and X
  nin               = lastindex(inodes) # number of internal nodes
  el                = lastindex(idf)    # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 12)

  # make Ξ vector
  Ξv = iTfbdX[]

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

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, δt, srδt) - crown*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # σλ & σμ update
      elseif p === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nx,
            llc, prc, mc, th, crown, δt, srδt, σλ_prior, σμ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, δt, srδt) - crown*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # psi update
      elseif p === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, δt, srδt) - crown*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # update beta
      elseif p === 4

        # to do!

      # update `x` diffusion rate
      elseif p === 5

        llc, prc, σxc = update_σx!(σxc, ssx, nx, llc, prc, σx_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, δt, srδt) - crown*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # update `x` bm
      elseif p === 6

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ssx =
          update_x!(bix, Ξ, idf, σxc, llc, prc, ssx, δt, srδt, x0_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, δt, srδt) - crown*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # gbm update
      elseif p === 7

        # nix = ceil(Int64,rand()*nin)
        # bix = inodes[nix]
        bix = 1

        llc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
            crown, δt, srδt, lλxpr, lμxpr)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, δt, srδt) - crown*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, p, Ξ
        #    return
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, ssx, nx, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, llc, dλ, 
            ssλ, ssμ, ssx, nx, L, δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, βλc, σxc, δt, srδt) - crown*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
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
        R[lit,9]  = ψc
        R[lit,10] = xv(Ξ[1])[1]
        R[lit,11] = βλc
        R[lit,12] = σxc

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
               Ξ   ::Vector{iTfbdX},
               idf ::Vector{iBffs},
               α   ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               ψ   ::Float64,
               llc ::Float64,
               dλ  ::Float64,
               ssλ ::Float64,
               ssμ ::Float64,
               nx  ::Float64,
               L   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Forward simulation proposal function for `gbmfbd`.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{iTfbdX},
                    idf ::Vector{iBffs},
                    α   ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Float64,
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

  bi = idf[bix]
  ξc = Ξ[bix]

  if it(bi)
    if isfossil(bi)
      ξp, llr = fsbi_ft(bi, ξc, α, σλ, σμ, ψ, βλ, σx, δt, srδt)
    else
      ξp, llr = fsbi_t(bi, ξc, α, σλ, σμ, ψ, βλ, σx, δt, srδt)
    end
    drλ  = 0.0
    ssrλ = 0.0
    ssrμ = 0.0
    ssrx = 0.0
  else
    if isfossil(bi)
      ξp, llr, drλ, ssrλ, ssrμ, ssrx =
        fsbi_i(bi, ξc, Ξ[d1(bi)], α, σλ, σμ, ψ, βλ, σx, δt, srδt)
    else
      ξp, llr, drλ, ssrλ, ssrμ, ssrx =
        fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, σμ, ψ, βλ, σx, δt, srδt)
    end
  end

  if isfinite(llr)
    ll1, dλ1, ssλ1, ssμ1, ssx1, nx1 =
      llik_gbm_ss(ξp, α, σλ, σμ, ψ, βλ, σx, δt, srδt)
    ll0, dλ0, ssλ0, ssμ0, ssx0, nx0 =
      llik_gbm_ss(ξc, α, σλ, σμ, ψ, βλ, σx, δt, srδt)

    # update llr, ssλ, ssμ, nx, sns, L
    llc += ll1  - ll0  + llr
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    ssμ += ssμ1 - ssμ0 + ssrμ
    ssx += ssx1 - ssx0 + ssrx
    nx  += nx1  - nx0
    L   += treelength(ξp) - treelength(ξc)

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, dλ, ssλ, ssμ, ssx, nx, L
end




"""
    fsbi_t(bi::iBffs,
           λ0  ::Float64,
           μ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Float64,
           βλ  ::Float64,
           σx  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs,
                ξc  ::iTfbdX,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Float64,
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
    _sim_gbmfbd_t(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ,xv(ξc)[1], βλ, σx,
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
    fsbi_ft(bi  ::iBffs,
            ξc  ::iTfbdX,
            α   ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            ψ   ::Float64,
            βλ  ::Float64,
            σx  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_ft(bi  ::iBffs,
                 ξc  ::iTfbdX,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 ψ   ::Float64,
                 βλ  ::Float64,
                 σx  ::Float64,
                 δt  ::Float64,
                 srδt::Float64)


  lU = -randexp() # log-probability

  # add sampling fraction
  nac = ni(bi)                # current ni
  Iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = - Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if fx(bi)

    xist = Float64[]
    xfst = Float64[]
    est  = Float64[]

    # forward simulation during branch length
    t0, nf, nn =
      _sim_gbmfbd_ifx(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, 
        xv(ξc)[1], βλ, σx, δt, srδt, 0, 1, 500, xist, xfst, est)

    na = lastindex(xfst)

    if na < 1 || nf > 0 || nn >= 500
      return t0, NaN
    end

    xc = fossil_xt(ξc)

    # sample tip according to fix `x` value
    wp = Float64[]
    ap = 0.0
    @simd for i in Base.OneTo(na)
      srt = sqrt(est[i])
      wi  = dnorm_bm(xfst[i], xist[i], srt*σx)/dnorm_bm(xc, xist[i], srt*σx)
      push!(wp, wi)
      ap += wi
    end
    acr += log(ap)

    if isfinite(acr) && lU < acr

      # sample tip
      wti = sample(wp)
      if wti <= div(na,2)
        fixtip1!(t0, wti,      0, xc, σx, δt, srδt)
      else
        fixtip2!(t0, na-wti+1, 0, xc, σx, δt, srδt)
      end
    else
      return t0, NaN
    end
  else

    t0, na, nf, nn =
      _sim_gbmfbd_i(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, 
        xv(ξc)[1], βλ, σx, δt, srδt, 0, 0, 1, 500)

    if na < 1 || nf > 0 || nn >= 500
      return t0, NaN
    end

    _fixrtip!(t0, na) # fix random tip
  end

  # simulate remaining tips until the present
  if na > 1
    tx, na, nn, acr =
      tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, βλ, σx, δt, srδt, acr, lU, Iρi, 
        na, nn)
  end

  if lU < acr

    # fossilize extant tip
    fossilizefixedtip!(t0)

    # if terminal fossil branch
    tx, na, nn, acr =
      fossiltip_sim!(t0, tf(bi), α, σλ, σμ, ψ, βλ, σx, δt, srδt, acr, lU, Iρi,
        na, nn)

    if lU < acr

      llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setni!(bi, na)       # set new ni

      return t0, llr
    end
  end

  return t0, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTfbdX,
           ξ1  ::iTfbdX,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Float64,
           βλ  ::Float64,
           σx  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for fossil internal branch `bi`.
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTfbdX,
                ξ1  ::iTfbdX,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # if fix `x` node
  if fx(bi)

    λsp  = Float64[]
    μsp  = Float64[]
    xist = Float64[]
    xfst = Float64[]
    est  = Float64[]

    t0, nf, nn =
      _sim_gbmfbd_i(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, 
        xv(ξc)[1], βλ, σx, δt, srδt, 0, 1, 500, λsp, μsp, xist, xfst, est)

    na = lastindex(λsp)

    if na < 1 || nf > 0 || nn >= 500
      return t0, NaN, NaN, NaN, NaN, NaN
    end

    # get current speciation rates at branch time
    λsc = λst(bi)
    μsc = μst(bi)

    e1  = e(ξ1)
    sr1 = sqrt(e1)
    λ1c = lλ(ξ1)
    μ1c = lμ(ξ1)
    l1  = lastindex(λ1c)
    λ1  = λ1c[l1]
    μ1  = μ1c[l1]

    xc = fossil_xt(ξc)

    # proposed acceptance ratio
    wp = Float64[]
    ap = 0.0
    @simd for i in Base.OneTo(lastindex(λsp))
      srt = sqrt(est[i])
      wi  = dnorm_bm(λsp[i], λ1 - α*e1, sr1*σλ) * dnorm_bm(μsp[i], μ1, sr1*σμ) *
            dnorm_bm(xfst[i], xist[i], srt*σx)/dnorm_bm(xc, xist[i], srt*σx)
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
      ac += dnorm_bm(λsc[i], λ1 - α*e1, sr1*σλ) * dnorm_bm(μsc[i], μ1, sr1*σμ)
    end
    ac = log(ac)

    acr = ap - ac
    xp  = xc

    # sample tip
    wti = sample(wp)

  else

    λsp = Float64[]
    μsp = Float64[]
    xsp = Float64[]

    t0, nf, nn =
      _sim_gbmfbd_i(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, 
        xv(ξc)[1], βλ, σx, δt, srδt, 0, 1, 500, λsp, μsp, xsp)

    na = lastindex(λsp)

    if na < 1 || nf > 0 || nn >= 500
      return t0, NaN, NaN, NaN, NaN, NaN
    end

    # get current speciation rates at branch time
    λsc = λst(bi)
    μsc = μst(bi)

    e1  = e(ξ1)
    sr1 = sqrt(e1)
    λ1c = lλ(ξ1)
    μ1c = lμ(ξ1)
    l1  = lastindex(λ1c)
    λ1  = λ1c[l1]
    μ1  = μ1c[l1]

    # proposed acceptance ratio
    wp = Float64[]
    ap = 0.0
    @simd for i in Base.OneTo(lastindex(λsp))
      wi  = dnorm_bm(λsp[i], λ1 - α*e1, sr1*σλ) * dnorm_bm(μsp[i], μ1, sr1*σμ)
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
      ac += dnorm_bm(λsc[i], λ1 - α*e1, sr1*σλ) * dnorm_bm(μsc[i], μ1, sr1*σμ)
    end
    ac = log(ac)

    acr = ap - ac

    # sample tip
    wti = sample(wp)
    xp  = xsp[wti]
  end

  lU = -randexp() # log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  λf  = λsp[wti]
  μf  = μsp[wti]

  llrd, acrd, drλ, ssrλ, ssrμ, ssrx, λ1p, μ1p, x1p =
    _daughter_update!(ξ1, λf, μf, α, σλ, σμ, xp, βλ, σx, δt, srδt)

  acr += acrd

  if lU < acr

    # fix tip
    if fx(bi)
      if wti <= div(na,2)
        fixtip1!(t0, wti,      0, xc, σx, δt, srδt)
      else
        fixtip2!(t0, na-wti+1, 0, xc, σx, δt, srδt)
      end
    else
      if wti <= div(na,2)
        fixtip1!(t0, wti, 0)
      else
        fixtip2!(t0, na-wti+1, 0)
      end
    end

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, βλ, σx, 
          δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      # fossilize extant tip
      fossilizefixedtip!(t0)
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setni!( bi, na)                       # set new ni
      setλst!(bi, λsp)                      # set new λst
      setμst!(bi, μsp)                      # set new μst
      unsafe_copyto!(λ1c,    1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(μ1c,    1, μ1p, 1, l1) # set new daughter 1 μ vector
      unsafe_copyto!(xv(ξ1), 1, x1p, 1, l1) # set new daughter 1 x vector

      return t0, llr, drλ, ssrλ, ssrμ, ssrx
    end
  end

  return t0, NaN, NaN, NaN, NaN, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTfbdX,
           ξ1  ::iTfbdX,
           ξ2  ::iTfbdX,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           βλ  ::Float64,
           σx  ::Float64,
           ψ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTfbdX,
                ξ1  ::iTfbdX,
                ξ2  ::iTfbdX,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  λsp = Float64[]
  μsp = Float64[]
  xsp = Float64[]

  t0, nf, nn =
    _sim_gbmfbd_i(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, xv(ξc)[1], βλ, σx, 
      δt, srδt, 0, 1, 500, λsp, μsp, xsp)

  na = lastindex(λsp)

  if na < 1 || nf > 0 || nn >= 500
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

    if wti <= div(na,2)
      fixtip1!(t0, wti, 0)
    else
      fixtip2!(t0, na-wti+1, 0)
    end

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, βλ, σx, δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setni!(bi,  na)                       # set new ni
      setλt!(bi,  λf)                       # set new λt
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
    tip_sims!(tree::iTfbdX,
              t   ::Float64,
              α   ::Float64,
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

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTfbdX,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   ψ   ::Float64,
                   βλ  ::Float64,
                   σx  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

  if lU < lr && nn < 500

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)
        lμ0  = lμ(tree)
        xv0  = xv(tree)
        l    = lastindex(lλ0)

        # simulate
        stree, na, nn, lr =
          _sim_gbmfbd_it(max(δt-fdti, 0.0), t, lλ0[l], lμ0[l], α, σλ, σμ, ψ,
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

        # merge to current tip
        if def1(stree)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr =
        tip_sims!(tree.d1, t, α, σλ, σμ, ψ, βλ, σx, 
          δt, srδt, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr =
        tip_sims!(tree.d2, t, α, σλ, σμ, ψ, βλ, σx, 
          δt, srδt, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::iTfbdX,
                   t   ::Float64,
                   α   ::Float64,
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
function fossiltip_sim!(tree::iTfbdX,
                        t   ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        ψ   ::Float64,
                        βλ  ::Float64,
                        σx  ::Float64,
                        δt  ::Float64,
                        srδt::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        Iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64)

  if lU < lr && nn < 500
    if istip(tree)
      stree, na, nn, lr =
        _sim_gbmfbd_it(t, lλ(tree)[end], lμ(tree)[end], α, σλ, σμ, ψ, 
          xv(tree)[end], βλ, σx, δt, srδt, lr, lU, Iρi, na-1, nn, 500)

      if isnan(lr) || nn >= 500
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr =
        fossiltip_sim!(tree.d1, t, α, σλ, σμ, ψ, βλ, σx, 
          δt, srδt, lr, lU, Iρi, na, nn)
    else
      tree.d2, na, nn, lr =
        fossiltip_sim!(tree.d2, t, α, σλ, σμ, ψ, βλ, σx, 
          δt, srδt, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTfbdX},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                σμ   ::Float64,
                llc  ::Float64,
                dλ   ::Float64,
                ssλ  ::Float64,
                ssμ  ::Float64,
                mc   ::Float64,
                th   ::Float64,
                crown ::Int64,
                δt   ::Float64,
                srδt ::Float64,
                lλxpr::Float64,
                lμxpr::Float64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTfbdX},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     σμ   ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     ssμ  ::Float64,
                     mc   ::Float64,
                     th   ::Float64,
                     crown ::Int64,
                     δt   ::Float64,
                     srδt ::Float64,
                     lλxpr::Float64,
                     lμxpr::Float64)

  ξi  = Ξ[bix]
  bi  = idf[bix]
  if !it(bi)
    id1 = d1(bi)
    ξ1  = Ξ[id1]
    if !isfossil(bi)
      id2 = d2(bi)
      ξ2  = Ξ[id2]
    end
  end

  root = iszero(pa(bi))
  # if crown
  if root && iszero(e(bi))
    llc, dλ, ssλ, ssμ, mc =
      _crown_update!(ξi, ξ1, ξ2, α, σλ, σμ, llc, dλ, ssλ, ssμ, mc, th,
        δt, srδt, lλxpr, lμxpr, crown)
    setλt!(bi, lλ(ξi)[1])
  else
    # if stem
    if root
      llc, dλ, ssλ, ssμ, mc =
        _stem_update!(ξi, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt,
          lλxpr, lμxpr)
    end

    # updates within the parent branch
    llc, dλ, ssλ, ssμ =
      _update_gbm!(ξi, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, it(bi))

    if !it(bi)
      # get fixed tip
      lξi = fixtip(ξi)

      # make between decoupled trees node update
      if isfossil(bi)
        llc, ssλ, ssμ, λf =
          update_duo!(lλ(lξi), lλ(ξ1), lμ(lξi), lμ(ξ1), e(lξi), e(ξ1),
            fdt(lξi), fdt(ξ1), α, σλ, σμ, llc, ssλ, ssμ, δt, srδt)
      else
        llc, dλ, ssλ, ssμ, λf =
          update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), lμ(lξi), lμ(ξ1), lμ(ξ2),
            e(lξi), e(ξ1), e(ξ2), fdt(lξi), fdt(ξ1), fdt(ξ2),
            α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt)
      end

      # set fixed `λ(t)` in branch
      setλt!(bi, lλ(lξi)[end])
    end
  end

  if !it(bi)
    # carry on updates in the daughters
    llc, dλ, ssλ, ssμ =
      _update_gbm!(ξ1, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, it(idf[id1]))
    if !isfossil(bi)
      llc, dλ, ssλ, ssμ =
        _update_gbm!(ξ2, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, it(idf[id2]))
    end
  end

  return llc, dλ, ssλ, ssμ, mc
end





"""
    update_x!(bix     ::Int64,
              Ξ       ::Vector{iTfbdX},
              idf     ::Vector{iBffs},
              σx      ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              ssx     ::Float64,
              δt      ::Float64,
              srδt    ::Float64,
              x0_prior::NTuple{2, Float64})

Make a `gbm` update for an internal branch and its descendants.
"""
function update_x!(bix     ::Int64,
                   Ξ       ::Vector{iTfbdX},
                   idf     ::Vector{iBffs},
                   σx      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   ssx     ::Float64,
                   δt      ::Float64,
                   srδt    ::Float64,
                   x0_prior::NTuple{2, Float64})

  ξi  = Ξ[bix]
  bi  = idf[bix]
  if !it(bi)
    id1 = d1(bi)
    id2 = d2(bi)
    ξ1  = Ξ[id1]
    if !isfossil(bi)
      ξ2  = Ξ[id2]
    end
  end

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
    llc, ssx = _update_x!(ξi, σx, llc, ssx, δt, srδt, !fx(bi))

    if !it(bi)
    # get fixed tip
      lξi = fixtip(ξi)

      # make between decoupled trees node update
      if isfossil(bi)
        llc, ssx = _update_duo_x!(lξi, ξ1, σx, llc, ssx, δt, srδt, !fx(bi))
      else
        llc, ssx = _update_triad_x!(lξi, ξ1, ξ2, σx, llc, ssx, δt, srδt)
      end
    end
  end

  if !it(bi)
    # carry on updates in the daughters
    llc, ssx = 
      _update_x!(ξ1, σx, llc, ssx, δt, srδt, !fx(idf[id1]) && it(idf[id1]))
    if !isfossil(bi)
    llc, ssx = 
      _update_x!(ξ2, σx, llc, ssx, δt, srδt, !fx(idf[id2]) && it(idf[id2]))
    end
  end

  return llc, prc, ssx
end





"""
    _update_x!(tree::iTfbdX,
               σx  ::Float64,
               llc ::Float64,
               ssx ::Float64,
               δt  ::Float64,
               srδt::Float64,
               ufx ::Bool)

Do gbm updates on a decoupled tree recursively.
"""
function _update_x!(tree::iTfbdX,
                    σx  ::Float64,
                    llc ::Float64,
                    ssx ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    ufx ::Bool)

  if def1(tree)
    if def2(tree)
      llc, ssx = 
        _update_triad_x!(tree, tree.d1, tree.d2, σx, llc, ssx, δt, srδt)
      llc, ssx = _update_x!(tree.d1, σx, llc, ssx, δt, srδt, ufx)
      llc, ssx = _update_x!(tree.d2, σx, llc, ssx, δt, srδt, ufx)
    else
      llc, ssx = _update_duo_x!(tree, tree.d1, σx, llc, ssx, δt, srδt, ufx)
      llc, ssx = _update_x!(tree.d1, σx, llc, ssx, δt, srδt, ufx)
    end
  elseif isfix(tree)
    llc, ssx = _update_tip_x!(tree, σx, llc, ssx, δt, srδt, ufx)
  else
    llc, ssx = _update_tip_x!(tree, σx, llc, ssx, δt, srδt, true)
  end

  return llc, ssx
end




"""
    _update_duo_x!(ξi  ::iTfbdX,
                   ξ1  ::iTfbdX,
                   σx  ::Float64,
                   llc ::Float64,
                   ssx ::Float64
                   δt  ::Float64,
                   srδt::Float64,
                   ufx ::Bool)

Make gibbs node update for trait.
"""
function _update_duo_x!(ξi  ::iTfbdX,
                        ξ1  ::iTfbdX,
                        σx  ::Float64,
                        llc ::Float64,
                        ssx ::Float64,
                        δt  ::Float64,
                        srδt::Float64,
                        ufx ::Bool)

  xca  = xv(ξi)
  xc1  = xv(ξ1)
  la   = lastindex(xca)
  l1   = lastindex(xc1)
  xa   = xca[1]
  x1   = xc1[l1]
  xpa  = Vector{Float64}(undef, la)
  xp1  = Vector{Float64}(undef, l1)
  fdta = fdt(ξi)
  fdt1 = fdt(ξ1)
  ea   = e(ξi)
  e1   = e(ξ1)

  # gibbs sampling
  if !ufx || iszero(ea) || iszero(e1)
    xn = xc1[1]
  else
    xn = duoprop(xa, x1, ea, e1, σx)
  end

  bb!(xpa, xa, xn, σx, δt, fdta, srδt)
  bb!(xp1, xn, x1, σx, δt, fdt1, srδt)

  llra, ssrxa = llr_ssx(xpa, xca, σx, δt, fdta, srδt)
  llr1, ssrx1 = llr_ssx(xp1, xc1, σx, δt, fdt1, srδt)

  unsafe_copyto!(xca, 1, xpa, 1, la)
  unsafe_copyto!(xc1, 1, xp1, 1, l1)

  # update llc, prc and ssx
  llc += llra + llr1
  ssx += ssrxa + ssrx1

  return llc, ssx
end



