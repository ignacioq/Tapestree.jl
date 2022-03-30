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
  idf = make_idf(tree, tρ)

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
      stem = 1
    # if stem conditioning
    else
      stem = 0
    end
  # no survival
  else
    stem = 2
  end
  mc = m_surv_gbmbd(th, log(λc), log(μc), αi, σλi, σμi, δt, srδt, 500, stem)

  # make a decoupled tree
  Ξ = make_Ξ(idf, log(λc), log(μc), αi, σλi, σμi, δt, srδt, iTfbd)

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
  for i in Base.OneTo(5)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running fossilized birth-death gbm"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, σμc, ψc, mc  =
    mcmc_burn_gbmbd(Ξ, idf,
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, ψ_prior,
      nburn, αi, σλi, σμi, ψc, mc, th, stem, δt, srδt, inodes, pup, prints)

  # mcmc
  R, Ξv =
    mcmc_gbmbd(Ξ, idf, llc, prc, αc, σλc, σμc, ψc, mc, th, stem,
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, ψ_prior,
      niter, nthin, δt, srδt, inodes, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "mu_root"      => 2,
                 "alpha"        => 3,
                 "sigma_lambda" => 4,
                 "sigma_mu"     => 5,
                 "psi"          => 6))

  write_ssr(R, pardic, out_file)

  return R, Ξv
end




"""
    mcmc_burn_gbmbd(Ξ       ::Vector{iTfbd},
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
                    stem    ::Bool,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

MCMC burn-in chain for `gbmbd`.
"""
function mcmc_burn_gbmbd(Ξ       ::Vector{iTfbd},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         μa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         σμ_prior::NTuple{2,Float64},
                         ψ_prior ::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         σμc     ::Float64,
                         ψc      ::Float64,
                         mc      ::Float64,
                         th      ::Float64,
                         stem    ::Int64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  λ0  = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, δt, srδt) - stem*λ0 +
        log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2,        σλ_prior[1], σλ_prior[2])  +
        logdinvgamma(σμc^2,        σμ_prior[1], σμ_prior[2])  +
        logdnorm(αc,               α_prior[1],  α_prior[2]^2) +
        logdunif(exp(λ0),          λa_prior[1], λa_prior[2])  +
        logdunif(exp(lμ(Ξ[1])[1]), μa_prior[1], μa_prior[2])

  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L            = treelength(Ξ)        # tree length
  nf           = Float64(nfossils(Ξ)) # number of fossilization events
  dλ           = deltaλ(Ξ)            # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ξ, αc)       # sum squares in λ and μ
  nin          = lastindex(inodes)    # number of internal nodes
  el           = lastindex(idf)       # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update α
      if pupi === 1

        llc, prc, αc, mc  =
          update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, L, dλ, llc, prc,
            mc, th, stem, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nλ = sss_gbm(Ξ, αc)

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nλ,
            llc, prc, mc, th, stem, δt, srδt, σλ_prior, σμ_prior)

      # psi update
      elseif pupi === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
            stem, δt, srδt, lλxpr, lμxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, nλ, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, ψc, llc, dλ, ssλ, ssμ, nλ, L,
            δt, srδt)

      end
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, σμc, ψc, mc
end




"""
    mcmc_gbmbd(Ξ       ::Vector{iTfbd},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               σμc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               stem    ::Bool,
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
function mcmc_gbmbd(Ξ       ::Vector{iTfbd},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    ψc      ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    stem    ::Int64,
                    λa_prior::NTuple{2,Float64},
                    μa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
                    ψ_prior ::NTuple{2,Float64},
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

  # crown or stem conditioning
  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L            = treelength(Ξ)        # tree length
  nf           = Float64(nfossils(Ξ)) # number of fossilization events
  dλ           = deltaλ(Ξ)            # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ξ, αc)       # sum squares in λ and μ
  nin          = lastindex(inodes)    # number of internal nodes
  el           = lastindex(idf)       # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 9)

  # make Ξ vector
  Ξv = iTfbd[]

  # number of branches and of triads
  nbr  = lastindex(idf)

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update α
      if pupi === 1

        llc, prc, αc, mc  =
          update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, L, dλ, llc, prc,
            mc, th, stem, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nλ = sss_gbm(Ξ, αc)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, δt, srδt) - stem*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ξ
        #    return
        # end

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nλ,
            llc, prc, mc, th, stem, δt, srδt, σλ_prior, σμ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, δt, srδt) - stem*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ξ
        #    return
        # end

      # psi update
      elseif pupi === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, δt, srδt) - stem*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ξ
        #    return
        # end

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
            stem, δt, srδt, lλxpr, lμxpr)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, δt, srδt) - stem*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ξ
        #    return
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, nλ, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, ψc, llc, dλ, ssλ, ssμ, nλ, L,
            δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, δt, srδt) - stem*lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ξ
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
        R[lit,5] = exp(lμ(Ξ[1])[1])
        R[lit,6] = αc
        R[lit,7] = σλc
        R[lit,8] = σμc
        R[lit,9] = ψc
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
               Ξ   ::Vector{iTfbd},
               idf ::Vector{iBffs},
               α   ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               ψ   ::Float64,
               llc ::Float64,
               dλ  ::Float64,
               ssλ ::Float64,
               ssμ ::Float64,
               nλ  ::Float64,
               L   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Forward simulation proposal function for `gbmfbd`.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{iTfbd},
                    idf ::Vector{iBffs},
                    α   ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Float64,
                    llc ::Float64,
                    dλ  ::Float64,
                    ssλ ::Float64,
                    ssμ ::Float64,
                    nλ  ::Float64,
                    L   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  bi = idf[bix]
  ξc = Ξ[bix]

  if it(bi) 
    if isfossil(bi)
      ξp, llr = fsbi_t(bi, ξc, α, σλ, σμ, ψ, δt, srδt)
      drλ  = 0.0
      ssrλ = 0.0
      ssrμ = 0.0
    else
      ξp, llr = fsbi_t(bi, lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, δt, srδt)
      drλ  = 0.0
      ssrλ = 0.0
      ssrμ = 0.0
    end
  else
    if isfossil(bi)
      ξp, llr, drλ, ssrλ, ssrμ =
        fsbi_i(bi, ξc, Ξ[d1(bi)], lλ(ξc)[1], lμ(ξc)[1], 
          α, σλ, σμ, ψ, δt, srδt)
    else
      ξp, llr, drλ, ssrλ, ssrμ =
        fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], lλ(ξc)[1], lμ(ξc)[1], 
          α, σλ, σμ, ψ, δt, srδt)
    end
  end

  if isfinite(llr)
    ll1, dλ1, ssλ1, ssμ1, nλ1 = llik_gbm_ss(ξp, α, σλ, σμ, ψ, δt, srδt)
    ll0, dλ0, ssλ0, ssμ0, nλ0 = llik_gbm_ss(ξc, α, σλ, σμ, ψ, δt, srδt)

    # update llr, ssλ, ssμ, nλ, sns, L
    llc += ll1  - ll0  + llr
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    ssμ += ssμ1 - ssμ0 + ssrμ
    nλ  += nλ1  - nλ0
    L   += treelength(ξp) - treelength(ξc)

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, dλ, ssλ, ssμ, nλ, L
end




"""
    fsbi_t(bi::iBffs, 
           λ0  ::Float64,
           μ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs, 
                λ0  ::Float64,
                μ0  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_gbmfbd_t(e(bi), λ0, μ0, α, σλ, σμ, ψ, δt, srδt, lc, lU, Iρi, 
      0, 1, 500)

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
           ξ1  ::iTfbd,
           λ0  ::Float64,
           μ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTfbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  t0, na, nf, nn =
    _sim_gbmfbd_i(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, δt, srδt, 
      0, 0, 1, 500)

  if na < 1 || nf > 0 || nn >= 500
    return t0, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if lU < acr

    λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr

      # fossilize extant tip
      fossilizefixedtip!(t0)

      # if terminal fossil branch
      tx, na, nn, acr = 
        fossiltip_sim!(t0, tf(bi), α, σλ, σμ, ψ, δt, srδt, acr, lU, Iρi, 
          na, nn)

      if lU < acr

        llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
        setni!(bi, na)       # set new ni
        setnt!(bi, ntp)      # set new nt
        setλt!(bi, λf)       # set new λt

        return t0, llr
      end
    end
  end

  return t0, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTfbd,
           ξ1  ::iTfbd,
           λ0  ::Float64,
           μ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTfbd,
                ξ1  ::iTfbd,
                λ0  ::Float64,
                μ0  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  t0, na, nf, nn =
    _sim_gbmfbd_i(e(bi), λ0, μ0, α, σλ, σμ, ψ, δt, srδt, 0, 0, 1, 500)

  if na < 1 || nf > 0 || nn >= 500
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, ssrμ, λ1p, μ1p, =
    _daughter_update!(ξ1, λf, μf, α, σλ, σμ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      # fossilize extant tip
      fossilizefixedtip!(t0)

      l1  = lastindex(λ1p)
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setni!(bi, na)       # set new ni
      setnt!(bi, ntp)      # set new nt
      setλt!(bi, λf)       # set new λt

      return t0, llr, drλ, ssrλ, ssrμ
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTfbd,
           ξ1  ::iTfbd,
           ξ2  ::iTfbd,
           λ0  ::Float64,
           μ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTfbd,
                ξ1  ::iTfbd,
                ξ2  ::iTfbd,
                λ0  ::Float64,
                μ0  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  t0, na, nf, nn =
    _sim_gbmfbd_i(e(bi), λ0, μ0, α, σλ, σμ, ψ, δt, srδt, 0, 0, 1, 500)

  if na < 1 || nf > 0 || nn >= 500
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, ssrμ, λ1p, λ2p, μ1p, μ2p =
      _daughters_update!(ξ1, ξ2, λf, μf, α, σλ, σμ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      l1 = lastindex(λ1p)
      l2 = lastindex(λ2p)
      setni!(bi, na)                        # set new ni
      setnt!(bi, ntp)                       # set new nt
      setλt!(bi, λf)                        # set new λt
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 2 λ vector
      unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector
      unsafe_copyto!(lμ(ξ2), 1, μ2p, 1, l2) # set new daughter 2 μ vector

      return t0, llr, drλ, ssrλ, ssrμ
    end
  end

  return t0, NaN, NaN, NaN, NaN
end





"""
    tip_sims!(tree::iTfbd,
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
function tip_sims!(tree::iTfbd,
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

  if lU < lr && nn < 500

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)
        lμ0  = lμ(tree)
        l    = lastindex(lλ0)

        # simulate
        stree, na, nn, lr = 
          _sim_gbmfbd_it(max(δt-fdti, 0.0), t, lλ0[l], lμ0[l], α, σλ, σμ, ψ, 
            δt, srδt, lr, lU, Iρi, na-1, nn, 500)

        if isnan(lr) || nn >= 500
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
        tip_sims!(tree.d1, t, α, σλ, σμ, ψ, δt, srδt, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, α, σλ, σμ, ψ, δt, srδt, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::iTfbd,
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
function fossiltip_sim!(tree::iTfbd,
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

  if lU < lr && nn < 500
    if istip(tree)
      stree, na, nn, lr = 
        _sim_gbmfbd_it(t, lλ(tree)[end], lμ(tree)[end], α, σλ, σμ, ψ, δt, srδt, 
          lr, lU, Iρi, na-1, nn, 500)

      if isnan(lr) || nn >= 500
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr = 
        fossiltip_sim!(tree.d1, t, α, σλ, σμ, ψ, δt, srδt, lr, lU, Iρi, na, nn)
    else
      tree.d2, na, nn, lr = 
        fossiltip_sim!(tree.d2, t, α, σλ, σμ, ψ, δt, srδt, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTfbd},
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
                stem ::Int64,
                δt   ::Float64,
                srδt ::Float64,
                lλxpr::Float64,
                lμxpr::Float64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTfbd},
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
                     stem ::Int64,
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
        δt, srδt, lλxpr, lμxpr, stem)
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


