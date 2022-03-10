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
  stem = !iszero(e(tree))

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
  mc = m_surv_gbmbd(th, log(λc), log(μc), αi, σλi, σμi, δt, srδt, 500, stem)

  # make a decoupled tree
  Ξ = iTfbd[]
  iTfbd!(Ξ, tree, δt, srδt, log(λc), log(μc), αi, σλi, σμi)

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
  Ξ, idf, llc, prc, αc, σλc, σμc, mc  =
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
                 "psi"     => 5))

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
                         stem    ::Bool,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)
  λ0  = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, δt, srδt) - !stem*λ0 + 
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
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, ψc, llc, dλ, ssλ, ssμ, mc, th,
            stem, δt, srδt, lλxpr, lμxpr)

      # forward simulation update
      else

        """
        here
        """


        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, nλ, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, nλ, L,
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

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # crown or stem conditioning
  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L            = treelength(Ξ)     # tree length
  dλ           = deltaλ(Ξ)         # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ξ, αc)    # sum squares in λ and μ
  nin          = lastindex(inodes) # number of internal nodes
  el           = lastindex(idf)    # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 8)

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

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ξ
        #    return
        # end

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nλ,
            llc, prc, mc, th, stem, δt, srδt, σλ_prior, σμ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ξ
        #    return
        # end

      # gbm update
      elseif pupi === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
            stem, δt, srδt, lλxpr, lμxpr)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ξ
        #    return
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, nλ, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, nλ, L,
            δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
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
        push!(Ξv, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ξv
end




"""
    update_fs!(bix    ::Int64,
               Ξ      ::Vector{iTfbd},
               idf    ::Vector{iBffs},
               α      ::Float64,
               σλ     ::Float64,
               σμ     ::Float64,
               llc    ::Float64,
               dλ     ::Float64,
               ssλ    ::Float64,
               ssμ    ::Float64,
               nλ     ::Float64,
               L      ::Float64,
               δt     ::Float64,
               srδt   ::Float64)

Forward simulation proposal function for `gbmbd`.
"""
function update_fs!(bix    ::Int64,
                    Ξ      ::Vector{iTfbd},
                    idf    ::Vector{iBffs},
                    α      ::Float64,
                    σλ     ::Float64,
                    σμ     ::Float64,
                    llc    ::Float64,
                    dλ     ::Float64,
                    ssλ    ::Float64,
                    ssμ    ::Float64,
                    nλ     ::Float64,
                    L      ::Float64,
                    δt     ::Float64,
                    srδt   ::Float64)

  bi  = idf[bix]
  itb = it(bi) # if is terminal

  ξc  = Ξ[bix]
  if !itb
    ξ1  = Ξ[d1(bi)]
    ξ2  = Ξ[d2(bi)]
  end

  # forward simulate an internal branch
  ξp, np, ntp, λf, μf = fsbi(bi, lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, δt, srδt)

  # check for survival or non-exploding simulation
  if ntp > 0

    ρbi = ρi(bi) # get branch sampling fraction
    nc  = ni(bi) # current ni
    ntc = nt(bi) # current nt

    # if terminal branch
    if itb
      llr  = log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc))
      acr  = llr
      drλ  = 0.0
      ssrλ = 0.0
      ssrμ = 0.0
    else
      np -= 1
      llr = log((1.0 - ρbi)^(np - nc))
      acr = llr + log(Float64(ntp)/Float64(ntc))
      # change daughters
      if isfinite(acr)

        llrd, acrd, drλ, ssrλ, ssrμ, λ1p, λ2p, μ1p, μ2p =
          _daughters_update!(ξ1, ξ2, λf, μf, α, σλ, σμ, δt, srδt)

        llr += llrd
        acr += acrd
      else
        return llc, dλ, ssλ, ssμ, nλ, L
      end
    end

    # MH ratio
    if -randexp() < acr

      ll1, dλ1, ssλ1, ssμ1, nλ1 = llik_gbm_ss(ξp, α, σλ, σμ, δt, srδt)
      ll0, dλ0, ssλ0, ssμ0, nλ0 = llik_gbm_ss(ξc, α, σλ, σμ, δt, srδt)

      # update llr, ssλ, nλ, sns, ne, L,
      llc += ll1  - ll0 + llr
      dλ  += dλ1  - dλ0  + drλ
      ssλ += ssλ1 - ssλ0 + ssrλ
      ssμ += ssμ1 - ssμ0 + ssrμ
      nλ  += nλ1  - nλ0
      L   += treelength(ξp)   - treelength(ξc)

      Ξ[bix] = ξp          # set new tree
      setni!(bi, np)       # set new ni
      setnt!(bi, ntp)      # set new nt
      setλt!(bi, λf)       # set new λt
      if !itb
        λ1c = lλ(ξ1)
        λ2c = lλ(ξ2)
        l1  = lastindex(λ1c)
        l2  = lastindex(λ2c)
        unsafe_copyto!(λ1c, 1, λ1p, 1, l1) # set new daughter 1 λ vector
        unsafe_copyto!(λ2c, 1, λ2p, 1, l2) # set new daughter 2 λ vector
        unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector
        unsafe_copyto!(lμ(ξ2), 1, μ2p, 1, l2) # set new daughter 2 μ vector
      end
    end
  end

  return llc, dλ, ssλ, ssμ, nλ, L
end




"""
    fsbi(bi  ::iBffs,
            λ0  ::Float64,
            μ0  ::Float64,
            α   ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi(bi  ::iBffs,
              λ0  ::Float64,
              μ0  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64)

  # times
  tfb = tf(bi)

  # forward simulation during branch length
  t0, na, nsp = _sim_gbmbd(e(bi), λ0, μ0, α, σλ, σμ, δt, srδt, 0, 1, 1_000)

  if na < 1 || nsp >= 1_000
    return iTfbd(0.0, 0.0, 0.0, false, false, Float64[], Float64[]),
      0, 0, NaN, NaN
  end

  nat = na

  if isone(na)
    f, λf, μf = fixalive!(t0, NaN, NaN)

    return t0, na, nat, λf, μf
  elseif na > 1
    # fix random tip
    λf, μf = fixrtip!(t0, na, NaN, NaN)

    if !it(bi)
      # add tips until the present
      tx, na, nsp = tip_sims!(t0, tfb, α, σλ, σμ, δt, srδt, na, nsp)

      if na < 1 || nsp >= 1_000
        return iTfbd(0.0, 0.0, 0.0, false, false, Float64[], Float64[]),
          0, 0, NaN, NaN
      end
    end

    return t0, na, nat, λf, μf
  end

  return iTfbd(0.0, 0.0, 0.0, false, false, Float64[], Float64[]),
      0, 0, NaN, NaN
end




"""
    tip_sims!(tree::iTfbd,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              na  ::Int64,
              nsp ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTfbd,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   na  ::Int64,
                   nsp ::Int64)

  if istip(tree)
    if !isfix(tree) && isalive(tree)

      fdti = fdt(tree)
      lλ0  = lλ(tree)
      lμ0  = lμ(tree)
      l    = lastindex(lλ0)

      # simulate
      stree, na, nsp =
        _sim_gbmbd(max(δt-fdti, 0.0), t, lλ0[l], lμ0[l], α, σλ, σμ, δt, srδt,
                   na - 1, nsp, 1_000)

      if na < 1 || nsp >= 1_000
        return tree, na, nsp
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

      if isdefined(stree, :d1)
        tree.d1 = stree.d1
        tree.d2 = stree.d2
      end
    end
  else
    tree.d1, na, nsp = tip_sims!(tree.d1, t, α, σλ, σμ, δt, srδt, na, nsp)
    tree.d2, na, nsp = tip_sims!(tree.d2, t, α, σλ, σμ, δt, srδt, na, nsp)
  end

  return tree, na, nsp
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
                stem ::Bool,
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
                     stem ::Bool,
                     δt   ::Float64,
                     srδt ::Float64,
                     lλxpr::Float64,
                     lμxpr::Float64)

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
  # if crown
  if root && !stem
    llc, dλ, ssλ, ssμ, mc =
      _crown_update!(ξi, ξ1, ξ2, α, σλ, σμ, llc, dλ, ssλ, ssμ, mc, th,
        δt, srδt, lλxpr, lμxpr)
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
        llc, dλ, ssλ, ssμ, λf =
          update_duo!(lλ(lξi), lλ(ξ1), lμ(lξi), lμ(ξ1), e(lξi), e(ξ1), 
            fdt(lξi), fdt(ξ1), α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt)
      else
        llc, dλ, ssλ, ssμ, λf =
          update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), lμ(lξi), lμ(ξ1), lμ(ξ2),
            e(lξi), e(ξ1), e(ξ2), fdt(lξi), fdt(ξ1), fdt(ξ2),
            α, σλ, σμ, ψ, llc, dλ, ssλ, ssμ, δt, srδt)

        # set fixed `λ(t)` in branch
        setλt!(bi, lλ(lξi)[end])
      end
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


