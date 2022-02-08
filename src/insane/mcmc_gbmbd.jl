#=

Anagenetic `gbmbd` MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmbd(tree    ::sT_label, 
                 out_file::String;
                 λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                 μa_prior::NTuple{2,Float64} = (0.0, 100.0),
                 α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                 σμ_prior::NTuple{2,Float64} = (0.05, 0.05),
                 niter   ::Int64             = 1_000,
                 nthin   ::Int64             = 10,
                 nburn   ::Int64             = 200,
                 ϵi      ::Float64           = 0.2,
                 λi      ::Float64           = NaN,
                 μi      ::Float64           = NaN,
                 αi      ::Float64           = 0.0,
                 σλi     ::Float64           = 0.01, 
                 σμi     ::Float64           = 0.01,
                 pupdp   ::NTuple{4,Float64} = (0.1, 0.1, 0.2, 0.2),
                 δt      ::Float64           = 1e-2,
                 prints  ::Int64             = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `gbm-bd`.
"""
function insane_gbmbd(tree    ::sT_label, 
                      out_file::String;
                      λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                      μa_prior::NTuple{2,Float64} = (0.0, 100.0),
                      α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                      σμ_prior::NTuple{2,Float64} = (0.05, 0.05),
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      ϵi      ::Float64           = 0.2,
                      λi      ::Float64           = NaN,
                      μi      ::Float64           = NaN,
                      αi      ::Float64           = 0.0,
                      σλi     ::Float64           = 0.01, 
                      σμi     ::Float64           = 0.01,
                      pupdp   ::NTuple{4,Float64} = (0.0, 0.1, 0.2, 0.2),
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
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ)

   # starting parameters (using method of moments)
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end
  mc = m_surv_gbmbd(th, log(λc), log(μc), αi, σλi, σμi, δt, srδt, 1_000, stem)

  # make a decoupled tree
  Ψ = iTgbmbd[]
  iTgbmbd!(Ψ, tree, δt, srδt, log(λc), log(μc), αi, σλi, σμi)

 # set end of fix branch speciation times and
  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    setλt!(bi, lλ(Ψ[i])[end])
    if !it(bi)
      push!(inodes, i)
    end
  end

  # parameter updates (1: α, 2: σλ, 3: σμ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death gbm"

  # burn-in phase
  Ψ, idf, llc, prc, αc, σλc, σμc, mc  =
    mcmc_burn_gbmbd(Ψ, idf, λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, 
      nburn, αi, σλi, σμi, mc, th, stem, δt, srδt, inodes, pup, prints)

  # mcmc
  R, Ψv =
    mcmc_gbmbd(Ψ, idf, llc, prc, αc, σλc, σμc, mc, th, stem,
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, niter, nthin, δt, srδt, 
      inodes, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "mu_root"      => 2,
                 "alpha"        => 3,
                 "sigma_lambda" => 4,
                 "sigma_mu"     => 5))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmbd(Ψ       ::Vector{iTgbmbd},
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
function mcmc_burn_gbmbd(Ψ       ::Vector{iTgbmbd},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         μa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         σμ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
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
  λ0  = lλ(Ψ[1])[1]
  nsi = stem ? 0.0 : λ0

  llc = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) - nsi + log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2,        σλ_prior[1], σλ_prior[2]) +
        logdinvgamma(σμc^2,        σμ_prior[1], σμ_prior[2]) +
        logdnorm(αc,               α_prior[1], α_prior[2]^2) +
        logdunif(exp(λ0),          λa_prior[1], λa_prior[2]) +
        logdunif(exp(lμ(Ψ[1])[1]), μa_prior[1], μa_prior[2])

  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L            = treelength(Ψ)      # tree length
  dλ           = deltaλ(Ψ)         # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ψ, αc)    # sum squares in λ and μ
  nin          = lastindex(inodes) # number of internal nodes
  el           = lastindex(idf)    # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update α
      if pupi === 1

        llc, prc, αc, mc  =
          update_α!(αc, lλ(Ψ[1])[1], lμ(Ψ[1])[1], σλc, σμc, L, dλ, llc, prc, 
            mc, th, stem, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nλ = sss_gbm(Ψ, αc)

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc, mc = 
          update_σ!(σλc, σμc, lλ(Ψ[1])[1], lμ(Ψ[1])[1], αc, ssλ, ssμ, nλ, 
            llc, prc, mc, th, stem, δt, srδt, σλ_prior, σμ_prior)

      # gbm update
      elseif pupi === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, ssμ, mc = 
          update_gbm!(bix, Ψ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th, 
            stem, δt, srδt, lλxpr, lμxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, nλ, L = 
          update_fs!(bix, Ψ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, nλ, L,
            δt, srδt)
      end
    end

    next!(pbar)
  end

  return Ψ, idf, llc, prc, αc, σλc, σμc, mc
end




"""
    mcmc_gbmbd(Ψ       ::Vector{iTgbmbd},
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
function mcmc_gbmbd(Ψ       ::Vector{iTgbmbd},
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

  L            = treelength(Ψ)     # tree length
  dλ           = deltaλ(Ψ)         # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ψ, αc)    # sum squares in λ and μ
  nin          = lastindex(inodes) # number of internal nodes
  el           = lastindex(idf)    # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 8)

  # make Ψ vector
  Ψv = iTgbmbd[]

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
          update_α!(αc, lλ(Ψ[1])[1], lμ(Ψ[1])[1], σλc, σμc, L, dλ, llc, prc, 
            mc, th, stem, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nλ = sss_gbm(Ψ, αc)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) - lλ(Ψ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ψ
        #    return 
        # end

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc, mc = 
          update_σ!(σλc, σμc, lλ(Ψ[1])[1], lμ(Ψ[1])[1], αc, ssλ, ssμ, nλ, 
            llc, prc, mc, th, stem, δt, srδt, σλ_prior, σμ_prior)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) - lλ(Ψ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ψ
        #    return 
        # end

      # gbm update
      elseif pupi === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, ssμ, mc = 
          update_gbm!(bix, Ψ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th, 
            stem, δt, srδt, lλxpr, lμxpr)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) - lλ(Ψ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ψ
        #    return 
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, nλ, L = 
          update_fs!(bix, Ψ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, nλ, L,
            δt, srδt)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) - lλ(Ψ[1])[1] + log(mc) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, pupi, Ψ
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
        R[lit,4] = exp(lλ(Ψ[1])[1])
        R[lit,5] = exp(lμ(Ψ[1])[1])
        R[lit,6] = αc
        R[lit,7] = σλc
        R[lit,8] = σμc
        push!(Ψv, couple(deepcopy(Ψ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ψv
end




"""
    update_fs!(bix    ::Int64,
               Ψ      ::Vector{iTgbmbd},
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
                    Ψ      ::Vector{iTgbmbd},
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

  ψc  = Ψ[bix]
  if !itb
    ψ1  = Ψ[d1(bi)]
    ψ2  = Ψ[d2(bi)]
  end

  # forward simulate an internal branch
  ψp, np, ntp, λf, μf = fsbi(bi, lλ(ψc)[1], lμ(ψc)[1], α, σλ, σμ, δt, srδt)

  # check for survival or non-exploding simulation
  if np > 0

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
          _daughters_update!(ψ1, ψ2, λf, μf, α, σλ, σμ, δt, srδt)

        llr += llrd
        acr += acrd
      else
        return llc, dλ, ssλ, ssμ, nλ, L
      end
    end

    # MH ratio
    if -randexp() < acr

      ll1, dλ1, ssλ1, ssμ1, nλ1 = llik_gbm_ss(ψp, α, σλ, σμ, δt, srδt)
      ll0, dλ0, ssλ0, ssμ0, nλ0 = llik_gbm_ss(ψc, α, σλ, σμ, δt, srδt)

      # update llr, ssλ, nλ, sns, ne, L,
      llc += ll1  - ll0 + llr
      dλ  += dλ1  - dλ0  + drλ
      ssλ += ssλ1 - ssλ0 + ssrλ
      ssμ += ssμ1 - ssμ0 + ssrμ
      nλ  += nλ1  - nλ0
      L   += treelength(ψp)   - treelength(ψc)

      Ψ[bix] = ψp          # set new tree
      setni!(bi, np)       # set new ni
      setnt!(bi, ntp)      # set new nt
      setλt!(bi, λf)       # set new λt
      if !itb
        λ1c = lλ(ψ1)
        λ2c = lλ(ψ2)
        l1  = lastindex(λ1c)
        l2  = lastindex(λ2c)
        unsafe_copyto!(λ1c, 1, λ1p, 1, l1) # set new daughter 1 λ vector
        unsafe_copyto!(λ2c, 1, λ2p, 1, l2) # set new daughter 2 λ vector
        unsafe_copyto!(lμ(ψ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector
        unsafe_copyto!(lμ(ψ2), 1, μ2p, 1, l2) # set new daughter 2 μ vector
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
    return iTgbmbd(), 0, 0, 0.0, 0.0
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
      tx, na = tip_sims!(t0, tfb, α, σλ, σμ, δt, srδt, na)
    end

    return t0, na, nat, λf, μf
  end

  return iTgbmbd(), 0, 0, 0.0, 0.0
end




"""
    tip_sims!(tree::iTgbmbd,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::iTgbmbd,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   na  ::Int64)

  if istip(tree) 
    if !isfix(tree) && isalive(tree)

      fdti = fdt(tree)
      lλ0  = lλ(tree)
      lμ0  = lμ(tree)
      l    = lastindex(lλ0)

      # simulate
      stree, na, nsp = 
        _sim_gbmbd(max(δt-fdti, 0.0), t, lλ0[l], lμ0[l], α, σλ, σμ, δt, srδt, 
                   na - 1, 1, 1_000)

      if !isdefined(stree, :lλ)
        return tree, 1_000
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
    tree.d1, na = tip_sims!(tree.d1, t, α, σλ, σμ, δt, srδt, na)
    tree.d2, na = tip_sims!(tree.d2, t, α, σλ, σμ, δt, srδt, na)
  end

  return tree, na
end




"""
    update_gbm!(bix  ::Int64,
                Ψ    ::Vector{iTgbmbd},
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
                     Ψ    ::Vector{iTgbmbd},
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

  @inbounds begin
    ψi   = Ψ[bix]
    bi   = idf[bix]
    ψ1   = Ψ[d1(bi)]
    ψ2   = Ψ[d2(bi)]
    ter1 = it(idf[d1(bi)]) 
    ter2 = it(idf[d2(bi)])

    root = iszero(pa(bi))
    # if crown
    if root && !stem
      llc, dλ, ssλ, ssμ, mc = 
        _crown_update!(ψi, ψ1, ψ2, α, σλ, σμ, llc, dλ, ssλ, ssμ, mc, th, 
          δt, srδt, lλxpr, lμxpr)
      setλt!(bi, lλ(ψi)[1])
    else
      # if stem
      if root 
        llc, dλ, ssλ, ssμ, mc = 
          _stem_update!(ψi, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, 
            lλxpr, lμxpr)
      end

      # updates within the parent branch
      llc, dλ, ssλ, ssμ = _update_gbm!(ψi, α, σλ, σμ, llc, dλ, ssλ, ssμ, 
        δt, srδt, false)

      # get fixed tip 
      lψi = fixtip(ψi) 

      # make between decoupled trees node update
      llc, dλ, ssλ, ssμ, λf = 
        update_triad!(lλ(lψi), lλ(ψ1), lλ(ψ2), lμ(lψi), lμ(ψ1), lμ(ψ2), 
          e(lψi), e(ψ1), e(ψ2), fdt(lψi), fdt(ψ1), fdt(ψ2), 
          α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt)

      # set fixed `λ(t)` in branch
      setλt!(bi, lλ(lψi)[end])
    end

    # carry on updates in the daughters
    llc, dλ, ssλ, ssμ = 
      _update_gbm!(ψ1, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter1)
    llc, dλ, ssλ, ssμ = 
      _update_gbm!(ψ2, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter2)
  end

  return llc, dλ, ssλ, ssμ, mc
end




"""
    update_α!(αc     ::Float64,
              λ0     ::Float64,
              μ0     ::Float64,
              σλ     ::Float64,
              σμ     ::Float64,
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
function update_α!(αc     ::Float64,
                   λ0     ::Float64,
                   μ0     ::Float64,
                   σλ     ::Float64,
                   σμ     ::Float64,
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

  mp  = m_surv_gbmbd(th, λ0, μ0, αp, σλ, σμ, δt, srδt, 1_000, stem)
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
    update_σ!(σλc     ::Float64,
              σμc     ::Float64,
              α       ::Float64,
              λ0      ::Float64,
              μ0      ::Float64,
              ssλ     ::Float64,
              ssμ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              mc      ::Float64,
              th      ::Float64,
              stem    ::Bool,
              δt      ::Float64,
              srδt    ::Float64,
              σλ_prior::NTuple{2,Float64},
              σμ_prior::NTuple{2,Float64})

Gibbs update for `σλ` and `σμ`.
"""
function update_σ!(σλc     ::Float64,
                   σμc     ::Float64,
                   α       ::Float64,
                   λ0      ::Float64,
                   μ0      ::Float64,
                   ssλ     ::Float64,
                   ssμ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   stem    ::Bool,
                   δt      ::Float64,
                   srδt    ::Float64,
                   σλ_prior::NTuple{2,Float64},
                   σμ_prior::NTuple{2,Float64})

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_prior[1] + 0.5 * n, σλ_prior[2] + ssλ)
  σμp2 = randinvgamma(σμ_prior[1] + 0.5 * n, σμ_prior[2] + ssμ)

  σλp = sqrt(σλp2)
  σμp = sqrt(σμp2)

  mp  = m_surv_gbmbd(th, λ0, μ0, α, σλp, σμp, δt, srδt, 1_000, stem)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += ssλ*(1.0/σλc^2 - 1.0/σλp^2) - n*(log(σλp/σλc)) + 
           ssμ*(1.0/σμc^2 - 1.0/σμp^2) - n*(log(σμp/σμc)) + 
           llr
    prc += llrdinvgamma(σλp2, σλc^2, σλ_prior[1], σλ_prior[2]) + 
           llrdinvgamma(σμp2, σμc^2, σμ_prior[1], σμ_prior[2])
    σλc  = σλp
    σμc  = σμp
    mc   = mp
  end

  return llc, prc, σλc, σμc, mc
end



