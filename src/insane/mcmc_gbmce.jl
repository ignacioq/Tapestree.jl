#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmce(tree    ::sT_label, 
                 out_file::String;
                 λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                 α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                 μ_prior ::NTuple{2,Float64} = (0.0, 100.0),
                 niter   ::Int64             = 1_000,
                 nthin   ::Int64             = 10,
                 nburn   ::Int64             = 200,
                 tune_int::Int64             = 100,
                 λi      ::Float64           = NaN,
                 αi      ::Float64           = 0.0,
                 σλi     ::Float64           = 0.01,
                 μi      ::Float64           = NaN,
                 ϵi      ::Float64           = 0.2,
                 μtni    ::Float64           = 1.0, 
                 obj_ar  ::Float64           = 0.234,
                 pupdp   ::NTuple{5,Float64} = (0.1,0.1,0.1,0.2,0.2),
                 nlim    ::Int64             = 500,
                 δt      ::Float64           = 1e-2,
                 prints  ::Int64             = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for GBM birth-death.
"""
function insane_gbmce(tree    ::sT_label, 
                      out_file::String;
                      λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                      α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                      μ_prior ::Float64           = 0.1,
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      tune_int::Int64             = 100,
                      λi      ::Float64           = NaN,
                      αi      ::Float64           = 0.0,
                      σλi     ::Float64           = 0.01,
                      μi      ::Float64           = NaN,
                      ϵi      ::Float64           = 0.2,
                      μtni    ::Float64           = 1.0, 
                      obj_ar  ::Float64           = 0.234,
                      pupdp   ::NTuple{5,Float64} = (0.1,0.1,0.1,0.2,0.2),
                      nlim    ::Int64             = 500,
                      δt      ::Float64           = 1e-2,
                      prints  ::Int64             = 5,
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

  # make fix tree directory
  idf = make_idf(tree, tρ)

   # starting parameters (using method of moments)
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), ti(idf[1]), ϵi)
  else
    λc, μc = λi, μi
  end

  # make a decoupled tree
  Ψ = iTgbmce[]
  iTgbmce!(Ψ, tree, δt, srδt, log(λc), αi, σλi)

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

  # parameter updates (1: α, 2: σλ, 3: μ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # conditioning functions
  sns = (BitVector(), BitVector(), BitVector())
  snodes! = make_snodes(idf, !iszero(e(tree)), iTgbmce)
  snodes!(Ψ, sns)
  scond, scond0 = make_scond(idf, !iszero(e(tree)), iTgbmce)

  @info "running birth-death gbm with constant μ"

  # burn-in phase
  Ψ, idf, llc, prc, αc, σλc, μc, μtn, sns =
    mcmc_burn_gbmce(Ψ, idf, λa_prior, α_prior, σλ_prior, μ_prior, 
      nburn, tune_int, αi, σλi, μc, μtni, sns, δt, srδt, inodes, pup, 
      prints, scalef, snodes!, scond, scond0)

  # mcmc
  R, Ψv =
    mcmc_gbmce(Ψ, idf, llc, prc, αc, σλc, μc, μtn, sns,
      λa_prior, α_prior, σλ_prior, μ_prior, niter, nthin, δt, srδt, 
      inodes, pup, prints, snodes!, scond, scond0)

  pardic = Dict(("lambda_root"  => 1,
                 "alpha"        => 2,
                 "sigma_lambda" => 3,
                 "mu"           => 4))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmce(Ψ       ::Vector{iTgbmce},
                    idf     ::Vector{iBffs},
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    μ_prior ::NTuple{2,Float64},
                    nburn   ::Int64,
                    tune_int::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    μc      ::Float64,
                    μtn     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    pup     ::Array{Int64,1},
                    prints  ::Int64,
                    scalef  ::Function,
                    snodes! ::Function,
                    scond   ::Function,
                    scond0  ::Function)

MCMC burn-in chain for `gbmce`.
"""
function mcmc_burn_gbmce(Ψ       ::Vector{iTgbmce},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         μ_prior ::Float64,
                         nburn   ::Int64,
                         tune_int::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         μc      ::Float64,
                         μtn     ::Float64,
                         sns     ::NTuple{3,BitVector},
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Vector{Int64},
                         pup     ::Vector{Int64},
                         prints  ::Int64,
                         scalef  ::Function,
                         snodes! ::Function,
                         scond   ::Function,
                         scond0  ::Function)

  ltn = 0
  lup = 0.0
  lac = 0.0

  llc = llik_gbm(Ψ, idf, αc, σλc, μc, δt, srδt) + 
        scond(Ψ, μc, sns) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])        + 
        logdunif(exp(lλ(Ψ[1])[1]), λa_prior[1], λa_prior[2]) +
        logdnorm(αc, α_prior[1], α_prior[2]^2)               +
        logdexp(μc, μ_prior)

  # maximum bounds according to unfiorm priors
  lλxpr = log(λa_prior[2])

  L       = treelength(Ψ)      # tree length
  dλ      = deltaλ(Ψ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ψ, αc)     # sum squares in λ
  ne      = 0.0                # number of extinction events
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update drift
      if pupi === 1

        llc, prc, αc  = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ψ, αc)

      # update sigma
      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior)

      # update extinction
      elseif pupi === 3

        llc, prc, μc, lac  = update_μ!(Ψ, llc, prc, μc, μtn, lac, ne, L, sns, 
          μ_prior, scond)

        lup += 1.0

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ = 
          update_gbm!(bix, Ψ, idf, αc, σλc, μc, llc, dλ, ssλ, sns, δt, 
            srδt, lλxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, nλ, ne, L = 
          update_fs!(bix, Ψ, idf, αc, σλc, μc, llc, dλ, ssλ, nλ, ne, L, 
            sns, δt, srδt, snodes!, scond0)
      end
    end

    # log tuning parameters
    ltn += 1
    if ltn === tune_int
      μtn = scalef(μtn, lac/lup)
      ltn = 0
    end

    next!(pbar)
  end

  return Ψ, idf, llc, prc, αc, σλc, μc, μtn, sns
end







"""
    mcmc_gbmce(Ψ       ::Vector{iTgbmce},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               μc      ::Float64,
               λa_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               μ_prior ::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Array{Int64,1},
               prints  ::Int64,
               snodes! ::Function,
               scond   ::Function,
               scond0  ::Function)

MCMC chain for `gbmce`.
"""
function mcmc_gbmce(Ψ       ::Vector{iTgbmce},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    μc      ::Float64,
                    μtn     ::Float64,
                    sns     ::NTuple{3,BitVector},
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    μ_prior ::Float64,
                    niter   ::Int64,
                    nthin   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64,
                    snodes! ::Function,
                    scond   ::Function,
                    scond0  ::Function)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # maximum bounds according to unfiorm priors
  lλxpr = log(λa_prior[2])

  L       = treelength(Ψ)            # tree length
  dλ      = deltaλ(Ψ)                # delta change in λ
  ssλ, nλ = sss_gbm(Ψ, αc)           # sum squares in λ
  ne      = Float64(ntipsextinct(Ψ)) # number of extinction events
  nin     = lastindex(inodes)        # number of internal nodes
  el      = lastindex(idf)           # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 7)

  # make Ψ vector
  Ψv = iTgbmce[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # check for extinct

      # update σλ or σμ
      if pupi === 1

        llc, prc, αc  = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ψ, αc)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, μc, δt, srδt) + scond(Ψ, μc, sns) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
        #    return 
        # end

      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, μc, δt, srδt) + scond(Ψ, μc, sns) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
        #    return 
        # end

      elseif pupi === 3

        llc, prc, μc  = 
          update_μ!(Ψ, llc, prc, μc, μtn, ne, L, sns, μ_prior, scond)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, μc, δt, srδt) + scond(Ψ, μc, sns) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
        #    return 
        # end

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ = 
          update_gbm!(bix, Ψ, idf, αc, σλc, μc, llc, dλ, ssλ, sns, δt, 
            srδt, lλxpr)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, μc, δt, srδt) + scond(Ψ, μc, sns) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
        #    return 
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, nλ, ne, L = 
          update_fs!(bix, Ψ, idf, αc, σλc, μc, llc, dλ, ssλ, nλ, ne, L, 
            sns, δt, srδt, snodes!, scond0)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, μc, δt, srδt) + scond(Ψ, μc, sns) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
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
        R[lit,5] = αc
        R[lit,6] = σλc
        R[lit,7] = μc
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
               Ψ      ::Vector{iTgbmce},
               idf    ::Vector{iBffs},
               α      ::Float64,
               σλ     ::Float64,
               μ      ::Float64,
               llc    ::Float64,
               dλ     ::Float64,
               ssλ    ::Float64,
               nλ     ::Float64,
               ne     ::Float64,
               L      ::Float64,
               sns    ::NTuple{3,BitVector},
               δt     ::Float64,
               srδt   ::Float64,
               snodes!::Function, 
               scond0 ::Function)

Forward simulation proposal function for `gbmce`.
"""
function update_fs!(bix    ::Int64,
                    Ψ      ::Vector{iTgbmce},
                    idf    ::Vector{iBffs},
                    α      ::Float64,
                    σλ     ::Float64,
                    μ      ::Float64,
                    llc    ::Float64,
                    dλ     ::Float64,
                    ssλ    ::Float64,
                    nλ     ::Float64,
                    ne     ::Float64,
                    L      ::Float64,
                    sns    ::NTuple{3,BitVector},
                    δt     ::Float64,
                    srδt   ::Float64,
                    snodes!::Function, 
                    scond0 ::Function)

  bi  = idf[bix]
  itb = it(bi) # if is terminal
  ρbi = ρi(bi) # get branch sampling fraction
  nc  = ni(bi) # current ni
  ntc = nt(bi) # current nt

  ψc  = Ψ[bix]
  if !itb
    ψ1  = Ψ[d1(bi)]
    ψ2  = Ψ[d2(bi)]
  end

  # forward simulate an internal branch
  ψp, np, ntp, λf = fsbi_ce(bi, lλ(ψc)[1], α, σλ, μ, δt, srδt)

  # check for survival or non-exploding simulation
  if np > 0

    # if terminal branch
    if itb
      llr  = log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc))
      acr  = llr
      drλ  = 0.0
      ssrλ = 0.0
    else
      np -= 1
      llr = log((1.0 - ρbi)^(np - nc))
      acr = llr + log(Float64(ntp)/Float64(ntc))
      # change daughters
      if isfinite(acr)

        llrd, acrd, drλ, ssrλ, λ1p, λ2p = 
          _daughters_update!(ψ1, ψ2, λf, α, σλ, μ, δt, srδt)

        llr += llrd
        acr += acrd
      else
        return llc, dλ, ssλ, nλ, ne, L
      end
    end

    # MH ratio
    if -randexp() < acr

      ll1, dλ1, ssλ1, nλ1 = llik_gbm_ssλ(ψp, α, σλ, μ, δt, srδt)
      ll0, dλ0, ssλ0, nλ0 = llik_gbm_ssλ(ψc, α, σλ, μ, δt, srδt)

      # if stem or crown conditioned
      scn = (iszero(pa(bi)) && e(bi) > 0.0) || 
             (isone(pa(bi)) && iszero(e(Ψ[1])))
      if scn
        llr += scond0(ψp, μ, itb) - scond0(ψc, μ, itb)
      end

      # update llr, ssλ, nλ, sns, ne, L,
      llr += ll1  - ll0
      dλ  += dλ1  - dλ0  + drλ
      ssλ += ssλ1 - ssλ0 + ssrλ
      nλ  += nλ1  - nλ0
      ne  += ntipsextinct(ψp) - ntipsextinct(ψc)
      L   += treelength(ψp)   - treelength(ψc)

      Ψ[bix] = ψp          # set new tree
      llc += llr           # set new likelihood
      if scn
        snodes!(Ψ, sns)    # set new sns
      end
      setni!(bi, np)       # set new ni
      setnt!(bi, ntp)      # set new nt
      setλt!(bi, λf)       # set new λt
      if !itb
        copyto!(lλ(ψ1), λ1p) # set new daughter 1 λ vector
        copyto!(lλ(ψ2), λ2p) # set new daughter 2 λ vector
      end
    end
  end

  return llc, dλ, ssλ, nλ, ne, L
end




"""
    fsbi_ce(bi  ::iBffs,
            λ0  ::Float64,
            α   ::Float64,
            σλ  ::Float64,
            μ   ::Float64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_ce(bi  ::iBffs,
                 λ0  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 μ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

  # times
  tfb = tf(bi)

  # forward simulation during branch length
  t0, na, nsp = _sim_gbmce(e(bi), λ0, α, σλ, μ, δt, srδt, 0, 1, 1_000)

  if nsp >= 1_000
    return iTgbmce(), 0, 0, 0.0
  end

  nat = na

  if isone(na)
    f, λf = fixalive!(t0, NaN)

    return t0, na, nat, λf
  elseif na > 1
    # fix random tip
    λf = fixrtip!(t0, na, NaN)

    if !it(bi)
      # add tips until the present
      tx, na = tip_sims!(t0, tfb, α, σλ, μ, δt, srδt, na)
    end

    return t0, na, nat, λf
  end

  return iTgbmce(), 0, 0, 0.0
end




"""
    tip_sims!(tree::iTgbmce,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::iTgbmce,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   na  ::Int64)

  if istip(tree) 
    if !isfix(tree) && isalive(tree)

      fdti = fdt(tree)
      lλ0  = lλ(tree)

      # simulate
      stree, na, nsp = 
        _sim_gbmce(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, μ, δt, srδt, 
                   na - 1, 1, 1_000)

      if !isdefined(stree, :lλ)
        return tree, 1_000
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
    tree.d1, na = tip_sims!(tree.d1, t, α, σλ, μ, δt, srδt, na)
    tree.d2, na = tip_sims!(tree.d2, t, α, σλ, μ, δt, srδt, na)
  end

  return tree, na
end





"""
    update_gbm!(bix  ::Int64,
                Ψ    ::Vector{iTgbmce},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                μ    ::Float64,
                llc  ::Float64,
                dλ   ::Float64,
                ssλ  ::Float64,
                sns  ::NTuple{3,BitVector},
                δt   ::Float64,
                srδt ::Float64,
                lλxpr::Float64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ψ    ::Vector{iTgbmce},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     μ    ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     sns  ::NTuple{3,BitVector},
                     δt   ::Float64,
                     srδt ::Float64,
                     lλxpr::Float64)

  @inbounds begin
    ψi   = Ψ[bix]
    bi   = idf[bix]
    ψ1   = Ψ[d1(bi)]
    ψ2   = Ψ[d2(bi)]
    ter1 = it(idf[d1(bi)]) 
    ter2 = it(idf[d2(bi)])

    cn = false
    # if crown root
    if iszero(pa(bi)) && iszero(e(bi))
      llc, dλ, ssλ = 
        _crown_update!(ψi, ψ1, ψ2, α, σλ, μ, llc, dλ, ssλ, δt, srδt, lλxpr)
      setλt!(bi, lλ(ψi)[1])

      # carry on updates in the crown daughter branches
      llc, dλ, ssλ = 
        _update_gbm!(ψ1, α, σλ, μ, llc, dλ, ssλ, δt, srδt, ter1, sns[2], 1)
      llc, dλ, ssλ = 
        _update_gbm!(ψ2, α, σλ, μ, llc, dλ, ssλ, δt, srδt, ter2, sns[3], 1)
    else
      # if stem branch
      if iszero(pa(bi))
        llc, dλ, ssλ = 
          _stem_update!(ψi, α, σλ, μ, llc, dλ, ssλ, δt, srδt, lλxpr)

        # updates within the stem branch in stem conditioning
        llc, dλ, ssλ = 
          _update_gbm!(ψi, α, σλ, μ, llc, dλ, ssλ, δt, srδt, false, sns[1], 1)

        # if observed node should be conditioned
        cn = sns[1][end]

      # if crown branch
      elseif isone(pa(bi)) && iszero(e(Ψ[1]))
        wsn = bix === d1(idf[pa(bi)]) ? 2 : 3
        sni = sns[wsn]
        # updates within the crown branch with crown conditioning
        llc, dλ, ssλ = 
          _update_gbm!(ψi, α, σλ, μ, llc, dλ, ssλ, δt, srδt, false, sni, 1)

        # if observed node should be conditioned
        if lastindex(sni) > 0
          cn = sni[end]
        end
      else
        # updates within the parent branch
        llc, dλ, ssλ = _update_gbm!(ψi, α, σλ, μ, llc, dλ, ssλ, δt, srδt, false)
      end

      # get fixed tip 
      lψi = fixtip(ψi) 

      # make between decoupled trees node update
      llc, dλ, ssλ = update_triad!(lλ(lψi), lλ(ψ1), lλ(ψ2), e(lψi), e(ψ1), e(ψ2), 
        fdt(lψi), fdt(ψ1), fdt(ψ2), α, σλ, μ, llc, dλ, ssλ, δt, srδt, cn)

      # set fixed `λ(t)` in branch
      setλt!(bi, lλ(lψi)[end])

      # carry on updates in the daughters
      llc, dλ, ssλ = _update_gbm!(ψ1, α, σλ, μ, llc, dλ, ssλ, δt, srδt, ter1)
      llc, dλ, ssλ = _update_gbm!(ψ2, α, σλ, μ, llc, dλ, ssλ, δt, srδt, ter2)
    end
  end

  return llc, dλ, ssλ
end






"""
    update_μ!(psi  ::Vector{iTgbmce},
              llc  ::Float64,
              μc   ::Float64,
              μtn  ::Float64,
              lac  ::Float64,
              ne   ::Float64,
              L    ::Float64,
              μmxpr::Float64,
              scond::Function)

MCMC update for `σ` with acceptance log.
"""
function update_μ!(psi  ::Vector{iTgbmce},
                   llc  ::Float64,
                   prc  ::Float64,
                   μc   ::Float64,
                   μtn  ::Float64,
                   lac  ::Float64,
                   ne   ::Float64,
                   L    ::Float64,
                   sns  ::NTuple{3,BitVector},
                   μprior::Float64,
                   scond::Function)

  # parameter proposal
  μp = mulupt(μc, μtn)::Float64

  # log likelihood and prior ratio
  μr   = log(μp/μc)
  llr  = ne*μr + L*(μc - μp) + scond(psi, μp, sns) - scond(psi, μc, sns)

  # prior ratio
  prr  = llrdexp_x(μp, μc, μprior)

  if -randexp() < (llr + prr + μr)
    μc   = μp
    llc += llr
    prc += prr
    lac += 1.0
  end

  return llc, prc, μc, lac
end




"""
    update_μ!(psi  ::Vector{iTgbmce},
              llc  ::Float64,
              μc   ::Float64,
              μtn  ::Float64,
              ne   ::Float64,
              L    ::Float64,
              sns  ::NTuple{3,BitVector},
              μmxpr::Float64,
              scond::Function)

MCMC update for `μ`.
"""
function update_μ!(psi  ::Vector{iTgbmce},
                   llc  ::Float64,
                   prc  ::Float64,
                   μc   ::Float64,
                   μtn  ::Float64,
                   ne   ::Float64,
                   L    ::Float64,
                   sns  ::NTuple{3,BitVector},
                   μprior::Float64,
                   scond::Function)

  # parameter proposal
  μp = mulupt(μc, μtn)::Float64

  # log likelihood and prior ratio
  μr   = log(μp/μc)
  llr  = ne*μr + L*(μc - μp) + scond(psi, μp, sns) - scond(psi, μc, sns)

  # prior ratio
  prr  = llrdexp_x(μp, μc, μprior)

  if -randexp() < (llr + prr + μr)
    llc += llr
    prc += prr
    μc   = μp
  end

  return llc, prc, μc
end

