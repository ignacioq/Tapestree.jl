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

  # make a decoupled tree
  Ψ = iTgbmct[]
  iTgbmct!(Ψ, tree, δt, srδt, log(λc), αi, σλi)

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

  # parameter updates (1: α, 2: σλ, 3: ϵ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # conditioning functions
  sns = (BitVector(), BitVector(), BitVector())
  snodes! = make_snodes(idf, !iszero(e(tree)), iTgbmct)
  snodes!(Ψ, sns)
  scond, scond0 = make_scond(idf, !iszero(e(tree)), iTgbmct)

  @info "running birth-death gbm with constant ϵ"

  # burn-in phase
  Ψ, idf, llc, prc, αc, σλc, ϵc, ϵtn, sns =
    mcmc_burn_gbmct(Ψ, idf, λa_prior, α_prior, σλ_prior, ϵ_prior, 
      nburn, tune_int, αi, σλi, ϵc, ϵtni, sns, δt, srδt, inodes, pup,
       prints, scalef, snodes!, scond, scond0)

  # mcmc
  R, Ψv = 
    mcmc_gbmct(Ψ, idf, llc, prc, αc, σλc, ϵc, ϵtn, sns,
      λa_prior, α_prior, σλ_prior, ϵ_prior, niter, nthin, δt, srδt, 
      inodes, pup, prints, snodes!, scond, scond0)

  pardic = Dict(("lambda_root"   => 1,
                 "alpha"        => 2,
                 "sigma_lambda" => 3,
                 "epsilon"      => 4))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmct(Ψ       ::Vector{iTgbmct},
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

MCMC burn-in chain for `gbmct`.
"""
function mcmc_burn_gbmct(Ψ       ::Vector{iTgbmct},
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

  llc = llik_gbm(Ψ, idf, αc, σλc, ϵc, δt, srδt) + 
        scond(Ψ, ϵc, sns) + prob_ρ(idf)

  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])        + 
        logdunif(exp(lλ(Ψ[1])[1]), λa_prior[1], λa_prior[2]) +
        logdnorm(αc, α_prior[1], α_prior[2]^2)               +
        logdunif(ϵc, ϵ_prior[1], ϵ_prior[2])

  lλxpr = log(λa_prior[2])
  ϵxpr  = ϵ_prior[2]

  L       = treelength(Ψ)      # tree length
  dλ      = deltaλ(Ψ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ψ, αc)     # sum squares in λ
  Σλ      = Σλ_gbm(Ψ)          # sum of λ
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

        llc, prc, αc  = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ψ, αc)

      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior)

      elseif pupi === 3

        llc, ϵc, lac = update_ϵ!(Ψ, llc, ϵc, ϵtn, lac, ne, Σλ, sns, ϵxpr, scond)

        lup += 1.0

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, Σλ = 
          update_gbm!(bix, Ψ, idf, αc, σλc, ϵc, llc, dλ, ssλ, Σλ, δt, 
            srδt, lλxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, Σλ, nλ, ne, L = 
          update_fs!(bix, Ψ, idf, αc, σλc, ϵc, llc, dλ, ssλ, Σλ, nλ, ne, L, 
            sns, δt, srδt, snodes!, scond0)
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

  return Ψ, idf, llc, prc, αc, σλc, ϵc, ϵtn, sns
end




"""
     mcmc_gbmce(Ψ       ::Vector{iTgbmct},
                idf     ::Vector{iBffs},
                llc     ::Float64,
                prc     ::Float64,
                αc      ::Float64,
                σλc     ::Float64,
                ϵc      ::Float64,
                ϵtn     ::Float64,
                sns     ::NTuple{3,BitVector},
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
                prints  ::Int64,
                snodes! ::Function,
                scond   ::Function,
                scond0  ::Function)

MCMC chain for `gbmct`.
"""
function mcmc_gbmct(Ψ       ::Vector{iTgbmct},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    ϵc      ::Float64,
                    ϵtn     ::Float64,
                    sns     ::NTuple{3,BitVector},
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
                    prints  ::Int64,
                    snodes! ::Function,
                    scond   ::Function,
                    scond0  ::Function)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # maximum bounds according to uniform priors
  lλxpr = log(λa_prior[2])
  ϵxpr  = ϵ_prior[2]

  L       = treelength(Ψ)            # tree length
  dλ      = deltaλ(Ψ)                # delta change in λ
  ssλ, nλ = sss_gbm(Ψ, αc)           # sum squares in λ
  Σλ      = Σλ_gbm(Ψ)                # sum of λ
  ne      = Float64(ntipsextinct(Ψ)) # number of extinction events
  nin     = lastindex(inodes)        # number of internal nodes
  el      = lastindex(idf)           # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 7)

  # make Ψ vector
  Ψv = iTgbmct[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    ii = 0
    # parameter updates
    for pupi in pup

      if pupi === 1

        llc, prc, αc  = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ψ, αc)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, ϵc, δt, srδt) + scond(Ψ, ϵc, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
        #    return 
        # end

      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, ϵc, δt, srδt) + scond(Ψ, ϵc, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
        #    return 
        # end

      elseif pupi === 3

        llc, ϵc = update_ϵ!(Ψ, llc, ϵc, ϵtn, ne, Σλ, sns, ϵxpr, scond)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, ϵc, δt, srδt) + scond(Ψ, ϵc, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
        #    return 
        # end

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, Σλ = 
          update_gbm!(bix, Ψ, idf, αc, σλc, ϵc, llc, dλ, ssλ, Σλ, δt, 
            srδt, lλxpr)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, ϵc, δt, srδt) + scond(Ψ, ϵc, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ψ
        #    return 
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, Σλ, nλ, ne, L = 
          update_fs!(bix, Ψ, idf, αc, σλc, ϵc, llc, dλ, ssλ, Σλ, nλ, ne, L, 
            sns, δt, srδt, snodes!, scond0)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, ϵc, δt, srδt) + scond(Ψ, ϵc, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-5)
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
        R[lit,7] = ϵc
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
               Ψ      ::Vector{iTgbmct},
               idf    ::Vector{iBffs},
               α      ::Float64,
               σλ     ::Float64,
               ϵ      ::Float64,
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
                    Ψ      ::Vector{iTgbmct},
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
                    sns    ::NTuple{3,BitVector},
                    δt     ::Float64,
                    srδt   ::Float64,
                    snodes!::Function, 
                    scond0 ::Function)

  bi  = idf[bix]
  itb = it(bi) # if is terminal
 
  ψc  = Ψ[bix]
  if !itb
    ψ1  = Ψ[d1(bi)]
    ψ2  = Ψ[d2(bi)]
    λ1 = lλ(ψ1)[end]
    λ2 = lλ(ψ2)[end]
    e1 = e(ψ1)
    e2 = e(ψ2)
  else
    λ1 = 0.0
    λ2 = 0.0
    e1 = 0.0
    e2 = 0.0
  end

  # forward simulate an internal branch
  ψp, ntp, np, λf, lU, acr = 
    fsbi_ct(bi, ψc, λ1, λ2, e1, e2, α, σλ, ϵ, δt, srδt)

  # check for survival or non-exploding simulation
  if np > 0

    ρbi = ρi(bi) # get branch sampling fraction
    nc  = ni(bi) # current ni
    ntc = nt(bi) # current nt

    # if terminal branch
    if itb
      llr   = log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc))
      acr  += llr
      drλ   = 0.0
      ssrλ  = 0.0
      Σrλ   = 0.0
    else
      np  -= 1
      llr  = log((1.0 - ρbi)^(np - nc))
      acr += llr + log(Float64(ntc)/Float64(ntp))
      # change daughters
      if isfinite(acr)

        llrd, acrd, drλ, ssrλ, Σrλ, λ1p, λ2p = 
          _daughters_update!(ψ1, ψ2, λf, α, σλ, ϵ, δt, srδt)

        llr += llrd
        acr += acrd
      else
        return llc, dλ, ssλ, Σλ, nλ, ne, L
      end
    end

    # MH ratio
    if lU < acr

      ll1, dλ1, ssλ1, Σλ1, nλ1 = llik_gbm_ssλ(ψp, α, σλ, ϵ, δt, srδt)
      ll0, dλ0, ssλ0, Σλ0, nλ0 = llik_gbm_ssλ(ψc, α, σλ, ϵ, δt, srδt)

      # if stem or crown conditioned
      scn = (iszero(pa(bi)) && e(bi) > 0.0) || 
            (isone(pa(bi)) && iszero(e(Ψ[1])))
      if scn
        llr += scond0(ψp, ϵ, itb) - scond0(ψc, ϵ, itb)
      end

      # update llr, ssλ, nλ, sns, ne, L,
      llc += llr + ll1  - ll0
      dλ  += dλ1  - dλ0  + drλ
      ssλ += ssλ1 - ssλ0 + ssrλ
      Σλ  += Σλ1  - Σλ0  + Σrλ
      nλ  += nλ1  - nλ0
      ne  += ntipsextinct(ψp) - ntipsextinct(ψc)
      L   += treelength(ψp)   - treelength(ψc)

      Ψ[bix] = ψp          # set new tree
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

  return llc, dλ, ssλ, Σλ, nλ, ne, L
end







"""
    fsbi_ct(bi  ::iBffs,
            λ0  ::Float64,
            λ1  ::Float64,
            λ2  ::Float64,
            e1  ::Float64,
            e2  ::Float64,
            α   ::Float64,
            σλ  ::Float64,
            ϵ   ::Float64,
            δt  ::Float64,
            srδt::Float64)
Forward simulation for branch `bi`
"""
function fsbi_ct(bi  ::iBffs,
                 ψc  ::iTgbmct,
                 λ1  ::Float64,
                 λ2  ::Float64,
                 e1  ::Float64,
                 e2  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 ϵ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

  # times
  tfb = tf(bi)

  # MH uniform
  lU = -randexp()

  # if terminal branch
  if it(bi)

    # forward simulation during branch length
    t0, na, nsp = _sim_gbmct(e(bi), lλ(ψc)[1], α, σλ, ϵ, δt, srδt, 
                    0, 1, 1_000)

    if nsp >= 1_000
      return iTgbmct(), 0, 0, 0.0, Inf, -Inf
    end

    λf = fixrtip!(t0, na, NaN)

    nat = na

    return t0, nat, na, λf, lU, 0.0

  # if internal branch
  else

    # tip rates
    λtsp = Float64[]

    # forward simulation during branch length
    t0, na, nsp = _sim_gbmct(e(bi), lλ(ψc)[1], α, σλ, ϵ, δt, srδt, 
                    0, 1, 1_000, λtsp)

    nat = na

    if nsp >= 1_000
      return iTgbmct(), 0, 0, 0.0, Inf, -Inf
    end

    # get tips -> daughters likelihoods for current
    λtsc = Float64[]
    _λat!(ψc, e(bi), λtsc, 0.0)

    push!(λtsc, λt(bi))

    # current MH `acr`
    acrc = 0.0
    for λi in λtsc
      acrc += exp(λi) * duodnorm(λi, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ)
    end
    acrc = log(acrc)

    # proposal MH `acr`
    wp   = Float64[]
    acrp = 0.0
    for λi in λtsp
      wi    = exp(λi) * duodnorm(λi, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ)
      acrp += wi
      push!(wp, wi)
    end
    acrp = log(acrp)

    # continue simulation only if acr on sum of tip rates is accepted
    acr = acrp - acrc

    if lU < acr

      # sample tip
      wti = sample(wp)

      # fix sampled tip
      lw = lastindex(wp)

      if wti <= div(lw,2)
        fixtip1!(t0, wti, 0)
      else
        fixtip2!(t0, lw - wti + 1, 0)
      end

      # simulated remaining tips until the present
      tx, na = tip_sims!(t0, tfb, α, σλ, ϵ, δt, srδt, na)

      return t0, nat, na, λtsp[wti], lU, acr
    end
  end

  return iTgbmct(), 0, 0, 0.0, Inf, -Inf
end







# """
#     fsbi_ct(bi  ::iBffs,
#             λ0  ::Float64,
#             α   ::Float64,
#             σλ  ::Float64,
#             ϵ   ::Float64,
#             δt  ::Float64,
#             srδt::Float64)

# Forward simulation for branch `bi`
# """
# function fsbi_ct(bi  ::iBffs,
#                  λ0  ::Float64,
#                  α   ::Float64,
#                  σλ  ::Float64,
#                  ϵ   ::Float64,
#                  δt  ::Float64,
#                  srδt::Float64)

#   # times
#   tfb = tf(bi)

#   # forward simulation during branch length
#   t0, na, nsp = _sim_gbmct(e(bi), λ0, α, σλ, ϵ, δt, srδt, 0, 1, 1_000)

#   if nsp >= 1_000
#     return iTgbmct(), 0, 0, 0.0
#   end

#   nat = na

#   if isone(na)
#     f, λf = fixalive!(t0, NaN)

#     return t0, na, nat, λf
#   elseif na > 1
#     # fix random tip
#     λf = fixrtip!(t0, na, NaN)

#     if !it(bi)
#       # add tips until the present
#       tx, na = tip_sims!(t0, tfb, α, σλ, ϵ, δt, srδt, na)
#     end

#     return t0, na, nat, λf
#   end

#   return iTgbmct(), 0, 0, 0.0
# end




"""
    tip_sims!(tree::iTgbmct,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              ϵ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::iTgbmct,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   ϵ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   na  ::Int64)

  if istip(tree) 
    if !isfix(tree) && isalive(tree)

      fdti = fdt(tree)
      lλ0  = lλ(tree)

      # simulate
      stree, na, nsp = 
        _sim_gbmct(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, ϵ, δt, srδt, 
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
    tree.d1, na = tip_sims!(tree.d1, t, α, σλ, ϵ, δt, srδt, na)
    tree.d2, na = tip_sims!(tree.d2, t, α, σλ, ϵ, δt, srδt, na)
  end

  return tree, na
end




"""
    update_gbm!(bix  ::Int64,
                Ψ    ::Vector{iTgbmct},
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
                     Ψ    ::Vector{iTgbmct},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     ϵ    ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     Σλ   ::Float64,
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
  end

  # if crown root
  if iszero(pa(bi)) && iszero(e(bi))
    llc, dλ, ssλ, Σλ = 
      _crown_update!(ψi, ψ1, ψ2, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, lλxpr)
    setλt!(bi, lλ(ψi)[1])
  else
    # if stem branch
    if iszero(pa(bi))
      llc, dλ, ssλ, Σλ = 
        _stem_update!(ψi, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, lλxpr)
    end

    # updates within the stem branch in stem conditioning
    llc, dλ, ssλ, Σλ = 
      _update_gbm!(ψi, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, false)

    # get fixed tip 
    lψi = fixtip(ψi) 

    # make between decoupled trees node update
    llc, dλ, ssλ, Σλ = 
      update_triad!(lλ(lψi), lλ(ψ1), lλ(ψ2), e(lψi), e(ψ1), e(ψ2), 
        fdt(lψi), fdt(ψ1), fdt(ψ2), α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt)

    # set fixed `λ(t)` in branch
    setλt!(bi, lλ(lψi)[end])

    # carry on updates in the daughters
    llc, dλ, ssλ, Σλ = 
      _update_gbm!(ψ1, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, ter1)
    llc, dλ, ssλ, Σλ = 
      _update_gbm!(ψ2, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, ter2)
  end

  return llc, dλ, ssλ, Σλ
end




"""
    update_ϵ!(psi  ::Vector{iTgbmct},
              llc  ::Float64,
              ϵc   ::Float64,
              ϵtn  ::Float64,
              lac  ::Float64,
              ne   ::Float64,
              Σλ   ::Float64,
              sns  ::NTuple{3,BitVector},
              ϵxpr ::Float64,
              scond::Function)

MCMC update for `σ` with acceptance log.
"""
function update_ϵ!(psi  ::Vector{iTgbmct},
                   llc  ::Float64,
                   ϵc   ::Float64,
                   ϵtn  ::Float64,
                   lac  ::Float64,
                   ne   ::Float64,
                   Σλ   ::Float64,
                   sns  ::NTuple{3,BitVector},
                   ϵxpr ::Float64,
                   scond::Function)

  # parameter proposal
  ϵp = mulupt(ϵc, ϵtn)::Float64

  # log likelihood and prior ratio
  ϵr  = log(ϵp/ϵc)
  llr = ne*ϵr + Σλ*(ϵc - ϵp) + scond(psi, ϵp, sns) - scond(psi, ϵc, sns)

  # prior ratio
  prr = ϵp > ϵxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + ϵr)
    ϵc   = ϵp
    llc += llr
    lac += 1.0
  end

  return llc, ϵc, lac
end




"""
    update_ϵ!(psi  ::Vector{iTgbmct},
              llc  ::Float64,
              ϵc   ::Float64,
              ϵtn  ::Float64,
              ne   ::Float64,
              Σλ   ::Float64,
              sns  ::NTuple{3,BitVector},
              ϵxpr ::Float64,
              scond::Function)

MCMC update for `ϵ`.
"""
function update_ϵ!(psi  ::Vector{iTgbmct},
                   llc  ::Float64,
                   ϵc   ::Float64,
                   ϵtn  ::Float64,
                   ne   ::Float64,
                   Σλ   ::Float64,
                   sns  ::NTuple{3,BitVector},
                   ϵxpr ::Float64,
                   scond::Function)

  # parameter proposal
  ϵp = mulupt(ϵc, ϵtn)::Float64

  # log likelihood and prior ratio
  ϵr   = log(ϵp/ϵc)
  llr = ne*ϵr + Σλ*(ϵc - ϵp) + scond(psi, ϵp, sns) - scond(psi, ϵc, sns)

  # prior ratio
  prr = ϵp > ϵxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + ϵr)
    ϵc   = ϵp
    llc += llr
  end

  return llc, ϵc
end

