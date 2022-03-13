#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmct(tree    ::sT_label,
                 out_file::String;
                 λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                 α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                 ϵ_prior ::NTuple{2,Float64} = (0.0, 100.0),
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
                 δt      ::Float64           = 1e-2,
                 prints  ::Int64             = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for GBM birth-death.
"""
function insane_gbmct(tree    ::sT_label,
                      out_file::String;
                      λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                      α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                      ϵ_prior ::NTuple{2,Float64} = (0.0, 100.0),
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
                      δt      ::Float64           = 1e-2,
                      prints  ::Int64             = 5,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)
  stem = !iszero(e(tree))

  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  idf = make_idf(tree, tρ)

   # starting parameters (using method of moments)
  if isnan(λi)
    λc, μc = moments(Float64(n), ti(idf[1]), ϵi)
  else
    λc = λi
  end
  ϵc = ϵi
  mc = m_surv_gbmct(th, log(λc), αi, σλi, ϵc, δt, srδt, 500, stem)

  # make a decoupled tree
  Ξ = iTct[]
  iTct!(Ξ, tree, δt, srδt, log(λc), αi, σλi)

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
      nburn, tune_int, αi, σλi, ϵc, ϵtni, mc, th, stem, δt, srδt, inodes, pup,
       prints, scalef)

  # mcmc
  R, Ξv =
    mcmc_gbmct(Ξ, idf, llc, prc, αc, σλc, ϵc, ϵtn, mc, th, stem,
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
                    stem    ::Bool,
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
                         stem    ::Bool,
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
  nsi = stem ? 0.0 : λ0

  llc = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) - nsi + log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])        +
        logdunif(exp(λ0), λa_prior[1], λa_prior[2]) +
        logdnorm(αc, α_prior[1], α_prior[2]^2)               +
        logdunif(ϵc, ϵ_prior[1], ϵ_prior[2])

  lλxpr = log(λa_prior[2])
  ϵxpr  = ϵ_prior[2]

  L       = treelength(Ξ)      # tree length
  dλ      = deltaλ(Ξ)          # delta change in λ
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
          update_α_ϵ!(αc, lλ(Ξ[1])[1], σλc, ϵc, L, dλ, llc, prc, mc, th, stem,
            δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

      elseif pupi === 2

        llc, prc, σλc, mc =
          update_σ_ϵ!(σλc, lλ(Ξ[1])[1], αc, ϵc, ssλ, nλ, llc, prc, mc, th, stem,
            δt, srδt, σλ_prior)

      elseif pupi === 3

        llc, ϵc, mc, lac =
          update_ϵ!(ϵc, lλ(Ξ[1])[1], αc, σλc, llc, mc, th, stem, ϵtn,
            lac, ne, Σλ, δt, srδt, ϵxpr)

        lup += 1.0

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, Σλ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, ϵc, llc, dλ, ssλ, Σλ, mc, th, stem,
            δt, srδt, lλxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, Σλ, nλ, ne, L =
          update_fs!(bix, Ξ, idf, αc, σλc, ϵc, llc, dλ, ssλ, Σλ, nλ, ne, L,
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
                prints  ::Int64

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
                    stem    ::Bool,
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
  lλxpr = log(λa_prior[2])
  ϵxpr  = ϵ_prior[2]

  L       = treelength(Ξ)            # tree length
  dλ      = deltaλ(Ξ)                # delta change in λ
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
          update_α_ϵ!(αc, lλ(Ξ[1])[1], σλc, ϵc, L, dλ, llc, prc, mc, th, stem,
            δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      elseif pupi === 2

        llc, prc, σλc, mc =
          update_σ_ϵ!(σλc, lλ(Ξ[1])[1], αc, ϵc, ssλ, nλ, llc, prc, mc, th, stem,
            δt, srδt, σλ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      elseif pupi === 3

        llc, ϵc, mc =
          update_ϵ!(ϵc, lλ(Ξ[1])[1], αc, σλc, llc, mc, th, stem, ϵtn,
            ne, Σλ, δt, srδt, ϵxpr)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, Σλ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, ϵc, llc, dλ, ssλ, Σλ, mc, th, stem,
            δt, srδt, lλxpr)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, Σλ, nλ, ne, L =
          update_fs!(bix, Ξ, idf, αc, σλc, ϵc, llc, dλ, ssλ, Σλ, nλ, ne, L,
            δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, ϵc, δt, srδt) - lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
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
               Ξ      ::Vector{iTct},
               idf    ::Vector{iBffs},
               α      ::Float64,
               σλ     ::Float64,
               ϵ      ::Float64,
               llc    ::Float64,
               dλ     ::Float64,
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
                    dλ     ::Float64,
                    ssλ    ::Float64,
                    Σλ     ::Float64,
                    nλ     ::Float64,
                    ne     ::Float64,
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
  ξp, ntp, np, λf  = fsbi_ct(bi, lλ(ξc)[1], α, σλ, ϵ, δt, srδt)

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
      Σrλ  = 0.0
    else
      np  -= 1
      llr = log((1.0 - ρbi)^(np - nc))
      acr = llr + log(Float64(ntp)/Float64(ntc))
      # change daughters
      if isfinite(acr)

        llrd, acrd, drλ, ssrλ, Σrλ, λ1p, λ2p =
          _daughters_update!(ξ1, ξ2, λf, α, σλ, ϵ, δt, srδt)

        llr += llrd
        acr += acrd
      else
        return llc, dλ, ssλ, Σλ, nλ, ne, L
      end
    end

    # MH ratio
    if -randexp() < acr

      ll1, dλ1, ssλ1, Σλ1, nλ1 = llik_gbm_ssλ(ξp, α, σλ, ϵ, δt, srδt)
      ll0, dλ0, ssλ0, Σλ0, nλ0 = llik_gbm_ssλ(ξc, α, σλ, ϵ, δt, srδt)

      # update llr, ssλ, nλ, sns, ne, L,
      llc += llr  + ll1  - ll0
      dλ  += dλ1  - dλ0  + drλ
      ssλ += ssλ1 - ssλ0 + ssrλ
      Σλ  += Σλ1  - Σλ0  + Σrλ
      nλ  += nλ1  - nλ0
      ne  += ntipsextinct(ξp) - ntipsextinct(ξc)
      L   += treelength(ξp)   - treelength(ξc)

      Ξ[bix] = ξp          # set new tree
      setni!(bi, np)       # set new ni
      setnt!(bi, ntp)      # set new nt
      setλt!(bi, λf)       # set new λt
      if !itb
        copyto!(lλ(ξ1), λ1p) # set new daughter 1 λ vector
        copyto!(lλ(ξ2), λ2p) # set new daughter 2 λ vector
      end
    end
  end

  return llc, dλ, ssλ, Σλ, nλ, ne, L
end




"""
    fsbi_ct(bi  ::iBffs,
            ξc  ::iTct,
            λ0  ::Float64,
            α   ::Float64,
            σλ  ::Float64,
            ϵ   ::Float64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_ct(bi  ::iBffs,
                 λ0  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 ϵ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

  # times
  tfb = tf(bi)

  # forward simulation during branch length
  t0, na, nsp = _sim_gbmct(e(bi), λ0, α, σλ, ϵ, δt, srδt, 0, 1, 1_000)

  if na < 1 || nsp >= 1_000
    return iTct(0.0, 0.0, 0.0, false, false, Float64[]), 0, 0, NaN
  end

  nat = na

  if isone(na)
    f, λf = fixalive!(t0, NaN)

    return t0, nat, na, λf
  elseif na > 1
    # fix random tip
    λf = fixrtip!(t0, na, NaN)

    if !it(bi)
      # add tips until the present
      tx, na, nsp = tip_sims!(t0, tfb, α, σλ, ϵ, δt, srδt, na, nsp)

      if nsp >= 1_000
        return iTct(0.0, 0.0, 0.0, false, false, Float64[]), 0, 0, NaN
      end
    end

    return t0, nat, na, λf
  end

  return iTct(0.0, 0.0, 0.0, false, false, Float64[]), 0, 0, NaN
end



"""
    tip_sims!(tree::iTct,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              ϵ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTct,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   ϵ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   na  ::Int64,
                   nsp ::Int64)

  if istip(tree)
    if !isfix(tree) && isalive(tree)

      fdti = fdt(tree)
      lλ0  = lλ(tree)

      # simulate
      stree, na, nsp =
        _sim_gbmct(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, ϵ, δt, srδt,
                   na - 1, nsp, 1_000)

      if nsp >= 1_000
        return tree, na, nsp
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
    tree.d1, na, nsp = tip_sims!(tree.d1, t, α, σλ, ϵ, δt, srδt, na, nsp)
    tree.d2, na, nsp = tip_sims!(tree.d2, t, α, σλ, ϵ, δt, srδt, na, nsp)
  end

  return tree, na, nsp
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTct},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                ϵ    ::Float64,
                llc  ::Float64,
                dλ   ::Float64,
                ssλ  ::Float64,
                δt   ::Float64,
                srδt ::Float64,
                lλxpr::Float64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTct},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     ϵ    ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     Σλ   ::Float64,
                     mc   ::Float64,
                     th   ::Float64,
                     stem ::Bool,
                     δt   ::Float64,
                     srδt ::Float64,
                     lλxpr::Float64)

  @inbounds begin
    ξi   = Ξ[bix]
    bi   = idf[bix]
    ξ1   = Ξ[d1(bi)]
    ξ2   = Ξ[d2(bi)]
    ter1 = it(idf[d1(bi)])
    ter2 = it(idf[d2(bi)])

    root = iszero(pa(bi))
    # if crown root
    if root && !stem
      llc, dλ, ssλ, Σλ, mc =
        _crown_update!(ξi, ξ1, ξ2, α, σλ, ϵ, llc, dλ, ssλ, Σλ, mc, th,
          δt, srδt, lλxpr)
      setλt!(bi, lλ(ξi)[1])
    else
      # if stem branch
      if root
        llc, dλ, ssλ, Σλ, mc =
          _stem_update!(ξi, α, σλ, ϵ, llc, dλ, ssλ, Σλ, mc, th, δt, srδt, lλxpr)
      end

      # updates within the stem branch in stem conditioning
      llc, dλ, ssλ, Σλ =
        _update_gbm!(ξi, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, false)

      # get fixed tip
      lξi = fixtip(ξi)

      # make node update between decoupled trees
      llc, dλ, ssλ, Σλ =
        update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2),
          fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt)

      # set fixed `λ(t)` in branch
      setλt!(bi, lλ(lξi)[end])
    end

    # carry on updates in the daughters
    llc, dλ, ssλ, Σλ =
      _update_gbm!(ξ1, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, ter1)
    llc, dλ, ssλ, Σλ =
      _update_gbm!(ξ2, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, ter2)
  end

  return llc, dλ, ssλ, Σλ, mc
end




"""
    update_α_ϵ!(αc     ::Float64,
                λ0     ::Float64,
                σλ     ::Float64,
                ϵ      ::Float64,
                L      ::Float64,
                dλ     ::Float64,
                llc    ::Float64,
                prc    ::Float64,
                mc     ::Float64,
                th     ::Float64,
                stem   ::Bool,
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
                     dλ     ::Float64,
                     llc    ::Float64,
                     prc    ::Float64,
                     mc     ::Float64,
                     th     ::Float64,
                     stem   ::Bool,
                     δt     ::Float64,
                     srδt   ::Float64,
                     α_prior::NTuple{2,Float64})

  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2
  αp  = rnorm((dλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  mp  = m_surv_gbmct(th, λ0, αp, σλ, ϵ, δt, srδt, 500, stem)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/L) + llr
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
                stem    ::Bool,
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
                     stem    ::Bool,
                     δt      ::Float64,
                     srδt    ::Float64,
                     σλ_prior::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)
  σλp  = sqrt(σλp2)

  mp  = m_surv_gbmct(th, λ0, α, σλp, ϵ, δt, srδt, 500, stem)
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
              stem ::Bool,
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
                   stem ::Bool,
                   ϵtn  ::Float64,
                   lac  ::Float64,
                   ne   ::Float64,
                   Σλ   ::Float64,
                   δt   ::Float64,
                   srδt ::Float64,
                   ϵxpr ::Float64)

  ϵp  = mulupt(ϵc, ϵtn)::Float64
  mp  = m_surv_gbmct(th, λ0, α, σλ, ϵp, δt, srδt, 500, stem)
  ϵr  = log(ϵp/ϵc)
  llr = ne*ϵr + Σλ*(ϵc - ϵp) + log(mp/mc)

  prr = ϵp > ϵxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + ϵr)
    llc += llr
    ϵc   = ϵp
    mc   = mp
    lac += 1.0
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
              stem ::Bool,
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
                   stem ::Bool,
                   ϵtn  ::Float64,
                   ne   ::Float64,
                   Σλ   ::Float64,
                   δt   ::Float64,
                   srδt ::Float64,
                   ϵxpr ::Float64)

  ϵp  = mulupt(ϵc, ϵtn)::Float64
  mp  = m_surv_gbmct(th, λ0, α, σλ, ϵp, δt, srδt, 500, stem)
  ϵr  = log(ϵp/ϵc)
  llr = ne*ϵr + Σλ*(ϵc - ϵp) + log(mp/mc)

  prr = ϵp > ϵxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + ϵr)
    llc += llr
    ϵc   = ϵp
    mc   = mp
  end

  return llc, ϵc, mc
end

