#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#





"""
    insane_gbmce(tree    ::sTbd, 
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
                 prints  ::Int64             = 5)

Run insane for GBM birth-death.
"""
function insane_gbmce(tree    ::sTbd, 
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
                      prints  ::Int64             = 5)

  # fix tree
  fixtree!(tree)

  # `n` tips, `th` treeheight define δt
  n    = sntn(tree, 0)
  th   = treeheight(tree)
  δt  *= th
  srδt = sqrt(δt)

   # starting parameters (using method of moments)
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # make Ψ current and proposal parameters
  Ψc = iTgbmce(tree, δt, srδt, log(λc), αi, σλi)
  Ψp = deepcopy(Ψc)

  # make fix Ψ directory
  idf = iBffs[]
  bit = BitArray{1}()
  makeiBf!(Ψc, idf, bit)

  # allocate `bb` for each fix branch and their `ts` vectors consisting of 
  # [pe(tree), fdt(tree)]
  bbλc = Array{Float64,1}[]
  tsv  = Array{Float64,1}[]

  makebbv!(Ψc, bbλc, tsv)

  bbλp = deepcopy(bbλc)

  # make trios
  triads, terminus, btotriad = make_triads(idf)

  # make survival conditioning function (stem or crown)
  svf = iszero(e(Ψc)) ? cond_surv_crown : cond_surv_stem

  # parameter updates (1: α, 2: σλ, 3: μ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death gbm with constant μ"

  # burn-in phase
  Ψp, Ψc, llc, prc, αc, σλc, μc, μtn =
    mcmc_burn_gbmce(Ψp, Ψc, bbλp, bbλc, tsv, λa_prior, α_prior, 
      σλ_prior, μ_prior, nburn, tune_int, αi, σλi, μc, μtni, δt, srδt, 
      idf, triads, terminus, btotriad, pup, nlim, prints, scalef, svf)

  # mcmc
  R, Ψv =
    mcmc_gbmce(Ψp, Ψc, llc, prc, αc, σλc, μc, μtn, bbλp, bbλc, tsv,
      λa_prior, α_prior, σλ_prior, μ_prior, niter, nthin, δt, srδt, 
      idf, triads, terminus, btotriad, pup, nlim, prints, svf)

  pardic = Dict(("lambda_root"  => 1,
                 "alpha"        => 2,
                 "sigma_lambda" => 3,
                 "mu"           => 4,
                 "n_extinct"    => 5))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmce(Ψp      ::iTgbmce,
                    Ψc      ::iTgbmce,
                    bbλp    ::Array{Array{Float64,1},1},
                    bbλc    ::Array{Array{Float64,1},1},
                    tsv     ::Array{Array{Float64,1},1},
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
                    idf     ::Array{iBffs,1},
                    triads  ::Array{Array{Int64,1},1},
                    terminus::Array{BitArray{1}},
                    btotriad::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    nlim    ::Int64,
                    prints  ::Int64,
                    scalef  ::Function,
                    svf     ::Function)

MCMC burn-in chain for `gbmce`.
"""
function mcmc_burn_gbmce(Ψp      ::iTgbmce,
                         Ψc      ::iTgbmce,
                         bbλp    ::Array{Array{Float64,1},1},
                         bbλc    ::Array{Array{Float64,1},1},
                         tsv     ::Array{Array{Float64,1},1},
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
                         idf     ::Array{iBffs,1},
                         triads  ::Array{Array{Int64,1},1},
                         terminus::Array{BitArray{1}},
                         btotriad::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         nlim    ::Int64,
                         prints  ::Int64,
                         scalef  ::Function,
                         svf     ::Function)

  ltn = 0
  lup = 0.0
  lac = 0.0

  # crown or stem conditioning
  icr = iszero(e(Ψc))

  llc = llik_gbm(Ψc, αc, σλc, μc, δt, srδt) + svf(Ψc, μc)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])      + 
        logdunif(exp(lλ(Ψc)[1]), λa_prior[1], λa_prior[2]) +
        logdnorm(αc, α_prior[1], α_prior[2]^2)             +
        logdunif(μc, μ_prior[1], μ_prior[2])

  lλmxpr = log(λa_prior[2])
  μmxpr  = log(μ_prior[2])

  # number of branches and of triads
  nbr  = lastindex(idf)
  ntr  = lastindex(triads)

  # make empty triad
  emptytriad = Int64[]
  emptyter   = BitArray([])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      if pupi === 1

        llc, prc, αc  = update_α!(αc, σλc, Ψc, llc, prc, α_prior)

      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, Ψc, llc, prc, σλ_prior)

      elseif pupi === 3

        llc, μc, lac  = update_μ!(μc, Ψc, llc, μtn, lac, μmxpr, svf)

        lup += 1.0

      # gbm update
      elseif pupi === 4

        tix = ceil(Int64,rand()*ntr)

        pr, d1, d2 = triads[tix]

        bi  = idf[pr]
        dri = dr(bi)
        ldr = length(dri)
        ter = terminus[tix]

        wbc = 23
        if iszero(sc(bi))
          wbc = 0
        elseif isone(sc(bi))
          wbc = 1
        end

        llc = lvupdate!(Ψp, Ψc, llc, bbλp, bbλc, tsv, pr, d1, d2,
            αc, σλc, μc, δt, srδt, lλmxpr, icr, wbc, dri, ldr, ter, 0)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*nbr)
        bi  = idf[bix]

        if it(bi)
          triad = emptytriad
          ter   = emptyter
        else
          triad = triads[btotriad[bix]]
          ter   = terminus[btotriad[bix]]
        end

        wbc = 23
        if iszero(sc(bi))
          wbc = 0
        elseif isone(sc(bi))
          wbc = 1
        end

        Ψp, Ψc, llc = fsp(Ψp, Ψc, bi, llc, αc, σλc, μc, tsv, bbλp, bbλc, 
              bix, triad, ter, δt, srδt, nlim, icr, wbc)

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

  return Ψp, Ψc, llc, prc, αc, σλc, μc, μtn
end






"""
     mcmc_gbmce(Ψp      ::iTgbmce,
                Ψc      ::iTgbmce,
                llc     ::Float64,
                prc     ::Float64,
                αc      ::Float64,
                σλc     ::Float64,
                μc      ::Float64,
                μtn     ::Float64,
                bbλp    ::Array{Array{Float64,1},1},
                bbλc    ::Array{Array{Float64,1},1},
                tsv     ::Array{Array{Float64,1},1},
                λa_prior::NTuple{2,Float64},
                α_prior ::NTuple{2,Float64},
                σλ_prior::NTuple{2,Float64},
                μ_prior ::NTuple{2,Float64},
                niter   ::Int64,
                nthin   ::Int64,
                δt      ::Float64,
                srδt    ::Float64,
                idf     ::Array{iBffs,1},
                triads  ::Array{Array{Int64,1},1},
                terminus::Array{BitArray{1}},
                btotriad::Array{Int64,1},
                pup     ::Array{Int64,1},
                nlim    ::Int64,
                prints  ::Int64,
                svf     ::Function)

MCMC chain for `gbmce`.
"""
function mcmc_gbmce(Ψp      ::iTgbmce,
                    Ψc      ::iTgbmce,
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    μc      ::Float64,
                    μtn     ::Float64,
                    bbλp    ::Array{Array{Float64,1},1},
                    bbλc    ::Array{Array{Float64,1},1},
                    tsv     ::Array{Array{Float64,1},1},
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    μ_prior ::NTuple{2,Float64},
                    niter   ::Int64,
                    nthin   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idf     ::Array{iBffs,1},
                    triads  ::Array{Array{Int64,1},1},
                    terminus::Array{BitArray{1}},
                    btotriad::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    nlim    ::Int64,
                    prints  ::Int64,
                    svf     ::Function)

  # crown or stem conditioning
  icr = iszero(e(Ψc))

  lλmxpr = log(λa_prior[2])
  μmxpr  = log(μ_prior[2])

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 8)

  # make Ψ vector
  Ψv = iTgbmce[]

  # number of branches and of triads
  nbr  = lastindex(idf)
  ntr  = lastindex(triads)

  # make empty triad
  emptytriad = Int64[]
  emptyter   = BitArray([])

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update σλ or σμ
      if pupi === 1

        llc, prc, αc  = update_α!(αc, σλc, Ψc, llc, prc, α_prior)

        # ll0 = llik_gbm(Ψc, αc, σλc, μc, δt, srδt) + svf(Ψc, μc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, 1, i, Ψc
        #    return 
        # end

      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, Ψc, llc, prc, σλ_prior)

        # ll0 = llik_gbm(Ψc, αc, σλc, μc, δt, srδt) + svf(Ψc, μc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, 2, i, Ψc
        #    return 
        # end

      elseif pupi === 3

        llc, μc = update_μ!(μc, Ψc, llc, μtn, μmxpr, svf)

        # ll0 = llik_gbm(Ψc, αc, σλc, μc, δt, srδt) + svf(Ψc, μc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, 3, i, Ψc
        #    return 
        # end

      # gbm update
      elseif pupi === 4

        tix = ceil(Int64,rand()*ntr)

        pr, d1, d2 = triads[tix]

        bi  = idf[pr]
        dri = dr(bi)
        ldr = length(dri)
        ter = terminus[tix]

        wbc = 23
        if iszero(sc(bi))
          wbc = 0
        elseif isone(sc(bi))
          wbc = 1
        end

        llc = lvupdate!(Ψp, Ψc, llc, bbλp, bbλc, tsv, pr, d1, d2,
            αc, σλc, μc, δt, srδt, lλmxpr, icr, wbc, dri, ldr, ter, 0)

        # ll0 = llik_gbm(Ψc, αc, σλc, μc, δt, srδt) + svf(Ψc, μc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, 4, i, Ψc
        #    return 
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*nbr)
        bi  = idf[bix]

        if it(bi)
          triad = emptytriad
          ter   = emptyter
        else
          triad = triads[btotriad[bix]]
          ter   = terminus[btotriad[bix]]
        end

        wbc = 23
        if iszero(sc(bi))
          wbc = 0
        elseif isone(sc(bi))
          wbc = 1
        end

        Ψp, Ψc, llc = fsp(Ψp, Ψc, bi, llc, αc, σλc, μc, tsv, bbλp, bbλc, 
              bix, triad, ter, δt, srδt, nlim, icr, wbc)

        # ll0 = llik_gbm(Ψc, αc, σλc, μc, δt, srδt) + svf(Ψc, μc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, 5, i, Ψc
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
        R[lit,4] = exp(lλ(Ψc)[1])
        R[lit,5] = αc
        R[lit,6] = σλc
        R[lit,7] = μc
        R[lit,8] = snen(Ψc, 0)
        push!(Ψv, deepcopy(Ψc))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ψv
end





"""
    fsp(Ψp   ::iTgbmce,
        Ψc   ::iTgbmce,
        bi   ::iBffs,
        llc  ::Float64,
        α    ::Float64, 
        σλ   ::Float64, 
        μ    ::Float64,
        tsv  ::Array{Array{Float64,1},1},
        bbλp ::Array{Array{Float64,1},1}, 
        bbλc ::Array{Array{Float64,1},1}, 
        bix  ::Int64,
        triad::Array{Int64,1},
        ter  ::BitArray{1},
        δt   ::Float64, 
        srδt ::Float64,
        nlim ::Int64,
        icr  ::Bool, 
        wbc  ::Int64)

Forward simulation proposal function for gbm birth-death.
"""
function fsp(Ψp   ::iTgbmce,
             Ψc   ::iTgbmce,
             bi   ::iBffs,
             llc  ::Float64,
             α    ::Float64, 
             σλ   ::Float64, 
             μ    ::Float64,
             tsv  ::Array{Array{Float64,1},1},
             bbλp ::Array{Array{Float64,1},1}, 
             bbλc ::Array{Array{Float64,1},1}, 
             bix  ::Int64,
             triad::Array{Int64,1},
             ter  ::BitArray{1},
             δt   ::Float64, 
             srδt ::Float64,
             nlim ::Int64,
             icr  ::Bool, 
             wbc  ::Int64)

  t0, ret, λf, λf1, dft0 = 
    fsbi_ce(bi, bbλc[bix][1], α, σλ, μ, δt, srδt, nlim)

  # if retain simulation
  if ret

    # get branch information
    dri = dr(bi)
    ldr = lastindex(dri)
    itb = it(bi)

    if !itb
      # make daughter proposal to be concordant with `t0`
      pr, d1, d2 = triad
      llr, acr = ldprop!(Ψp, Ψc, λf, bbλp, bbλc,
        tsv, pr, d1, d2, α, σλ, μ, icr, wbc, δt, srδt, dri, ldr, ter, 0)

      # change last event by speciation for llr
      iλ = 0.5*(λf1 + λf) + log(dft0) + 
           dft0*(exp(0.5*(λf1 + λf)) + μ)

      bbλi = bbλc[pr]
      l    = lastindex(bbλi)

      # acceptance ratio
      acr += 0.5*(λf1 + λf) - 0.5*(bbλi[l-1] + bbλi[l])
    else
      pr  = bix
      iλ  = 0.0
      llr = 0.0
      acr = 0.0
    end


    # mh ratio
    if -randexp() < acr

      llr += llik_gbm( t0, α, σλ, μ, δt, srδt) + iλ - 
             br_ll_gbm(Ψc, α, σλ, μ, δt, srδt, dri, ldr, 0)

      if icr && isone(wbc)
        if dri[1]
          llr += cond_surv_stem_p(t0, μ) - 
                 cond_surv_stem(  Ψc.d1, μ)
        else
          llr += cond_surv_stem_p(t0, μ) -
                 cond_surv_stem(  Ψc.d2, μ)
        end
      elseif iszero(wbc)
        llr += cond_surv_stem_p(t0, μ) -
               cond_surv_stem(  Ψc, μ) 
      end

      llc += llr

      # copy parent to aid vectors
      gbm_copy_f!(t0, bbλc[pr], 0)

      # copy daughters vectors
      if !itb
        pr, d1, d2 = triad

        copyto!(bbλc[d1], bbλp[d1])
        copyto!(bbλc[d2], bbλp[d2])
        gbm_copy_dsf!(Ψc, Ψp, dri, ldr, 0)
      end

      # make combined swap branch
      Ψp, Ψc = swapbranch!(Ψp, Ψc, t0, dri, ldr, itb, 0)
    end
  end

  return Ψp, Ψc, llc
end





"""
    fsbi_ce(bi  ::iBffs, 
            iλ  ::Float64, 
            α   ::Float64, 
            σλ  ::Float64, 
            μ   ::Float64, 
            δt  ::Float64, 
            srδt::Float64,
            nlim::Int64)

Forward gbmce simulation for branch `bi`.
"""
function fsbi_ce(bi  ::iBffs, 
                 iλ  ::Float64, 
                 α   ::Float64, 
                 σλ  ::Float64, 
                 μ   ::Float64, 
                 δt  ::Float64, 
                 srδt::Float64,
                 nlim::Int64)

  # retain the simulation?
  ret = true

  # times
  tfb = tf(bi)

  # simulate tree
  t0, nsp = _sim_gbmce(ti(bi) - tfb, iλ, α, σλ, μ, δt, srδt, 1, nlim)

  na = snan(t0, 0)

  λf, λf1, dft0 = NaN, NaN, NaN

  # if simulation goes extinct or maximum number of species reached
  if iszero(na) || nsp === nlim
    ret = false
  # if one surviving lineage
  elseif isone(na)
    f, λf, λf1, dft0 = fixalive!(t0, NaN, NaN, NaN)
  elseif na > 1
    # if terminal branch
    if it(bi)
      ret = false
    # if continue the simulation
    else
      # fix random tip and return end λ(t)
      λf, λf1, dft0 = fixrtip!(t0, na, NaN, NaN, NaN)

      for j in Base.OneTo(na - 1)
        # get their final λ to continue forward simulation
        ix, λt, fdti = fλ1(t0, NaN, NaN, false)

        for i in Base.OneTo(2)
          st0, nsp = 
            _sim_gbmce(max(δt - fdti, 0.0), tfb, λt, α, σλ, μ, δt, srδt, 1, nlim)
          # if maximum number of species reached.
          if nsp === nlim
            if i === 2
              ret = false
            end
            continue
          end
          # if goes extinct before the present
          if iszero(snan(st0, 0))
            # graft to tip
            addtotip(t0, st0, false)
            break
          end
          if i === 2
            ret = false
          end
        end
        !ret && break
      end
    end
  end

  # speciates at time `tfb`
  if iszero(dft0)
    ret = false
  end

  return t0, ret, λf, λf1, dft0
end




"""
    ldprop!(treep::iTgbmce,
            treec::iTgbmce,
            λf   ::Float64,
            bbλp ::Array{Array{Float64,1},1}, 
            bbλc ::Array{Array{Float64,1},1}, 
            tsv  ::Array{Array{Float64,1},1},
            pr   ::Int64,
            d1   ::Int64,
            d2   ::Int64,
            α    ::Float64,
            σλ   ::Float64, 
            μ    ::Float64, 
            icr  ::Bool, 
            wbc  ::Int64,
            δt   ::Float64, 
            srδt ::Float64, 
            dri  ::BitArray{1},
            ldr  ::Int64,
            ter  ::BitArray{1},
            ix   ::Int64)

Make a gbm update for speciation and extinction for a fixed triad.
"""
function ldprop!(treep::iTgbmce,
                 treec::iTgbmce,
                 λf   ::Float64,
                 bbλp ::Array{Array{Float64,1},1}, 
                 bbλc ::Array{Array{Float64,1},1}, 
                 tsv  ::Array{Array{Float64,1},1},
                 pr   ::Int64,
                 d1   ::Int64,
                 d2   ::Int64,
                 α    ::Float64,
                 σλ   ::Float64, 
                 μ    ::Float64, 
                 icr  ::Bool, 
                 wbc  ::Int64,
                 δt   ::Float64, 
                 srδt ::Float64, 
                 dri  ::BitArray{1},
                 ldr  ::Int64,
                 ter  ::BitArray{1},
                 ix   ::Int64)

  if ix === ldr 

    llr, acr = daughters_lprop!(treep, treec, λf,
      bbλp, bbλc, tsv, pr, d1, d2, ter, α, σλ, μ, icr, wbc, δt, srδt)

  elseif ix < ldr

    ifx1 = isfix(treec.d1)
    if ifx1 && isfix(treec.d2)
      ix += 1
      if dri[ix]
        llr, acr = 
          ldprop!(treep.d1, treec.d1, λf, bbλp, bbλc, 
            tsv, pr, d1, d2, α, σλ, μ, icr, wbc, δt, srδt, dri, ldr, ter, ix)
      else
        llr, acr = 
          ldprop!(treep.d2, treec.d2, λf, bbλp, bbλc, 
            tsv, pr, d1, d2, α, σλ, μ, icr, wbc, δt, srδt, dri, ldr, ter, ix)
      end
    elseif ifx1
      llr, acr = 
        ldprop!(treep.d1, treec.d1, λf, bbλp, bbλc, 
          tsv, pr, d1, d2, α, σλ, μ, icr, wbc, δt, srδt, dri, ldr, ter, ix)
    else
      llr, acr = 
        ldprop!(treep.d2, treec.d2, λf, bbλp, bbλc, 
          tsv, pr, d1, d2, α, σλ, μ, icr, wbc, δt, srδt, dri, ldr, ter, ix)
    end
  end

  return llr, acr
end




"""
    addtotip(tree::T, stree::T, ix::Bool) where {T < iTgbm}

Add `stree` to tip in `tree` given by `it` in `tree.d1` order.
"""
function addtotip(tree::T, stree::T, ix::Bool) where {T <: iTgbm}

  if istip(tree)
    if isalive(tree) && !isfix(tree)

      sete!(tree, e(tree) + e(stree))
      setproperty!(tree, :iμ, isextinct(stree))

      lλ0 = lλ(tree)
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

      ix = true
    end

    return ix
  end

  if !ix
    ix = addtotip(tree.d1, stree, ix)
  end
  if !ix
    ix = addtotip(tree.d2, stree, ix)
  end

  return ix
end




"""
    fixrtip!(tree::T, 
             na  ::Int64, 
             λf  ::Float64, 
             dft0::Float64) where {T <: iTgbm}

Fixes the the path for a random non extinct tip and returns final `λ(t)`.
"""
function fixrtip!(tree::T, 
                  na  ::Int64, 
                  λf  ::Float64, 
                  λf1 ::Float64,
                  dft0::Float64) where {T <: iTgbm}

  fix!(tree)

  if isdefined(tree, :d1)
    if isextinct(tree.d1)
      λf, λf1, dft0 = fixrtip!(tree.d2, na, λf, λf1, dft0)
    elseif isextinct(tree.d2)
      λf, λf1, dft0 = fixrtip!(tree.d1, na, λf, λf1, dft0)
    else
      na1 = snan(tree.d1, 0)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        λf, λf1, dft0 = fixrtip!(tree.d2, na - na1, λf, λf1, dft0)
      else
        λf, λf1, dft0 = fixrtip!(tree.d1, na1,      λf, λf1, dft0)
      end
    end
  else
    λv   = lλ(tree)
    l    = lastindex(λv)
    λf   = λv[l]
    λf1  = λv[l-1]
    dft0 = fdt(tree)
  end

  return λf, λf1, dft0
end




"""
    fixalive!(tree::T, λf::Float64, λf1::Float64, dft0::Float64) where {T <:iTgbm} 

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::T, λf::Float64, λf1::Float64, dft0::Float64) where {T <:iTgbm} 

  if istip(tree) 
    if isalive(tree)
      fix!(tree)
      λv   = lλ(tree)
      l    = lastindex(λv)
      λf   = λv[l]
      λf1  = λv[l-1]
      dft0 = fdt(tree)
      return true, λf, λf1, dft0
    end
  else
    f, λf, λf1, dft0 = fixalive!(tree.d2, λf, λf1, dft0)
    if f 
      fix!(tree)
      return true, λf, λf1, dft0
    end
    f, λf, λf1, dft0 = fixalive!(tree.d1, λf, λf1, dft0)
    if f 
      fix!(tree)
      return true, λf, λf1, dft0
    end
  end

  return false, λf, λf1, dft0
end




"""
    fλ1(tree::T, λt::Float64, fdti::Float64, ix::Bool) where {T <: iTgbm}

Get final `λ(t)` for a `it` tip in `tree` given in `tree.d1` order
not taking into account the fixed tip.
"""
function fλ1(tree::T, λt::Float64, fdti::Float64, ix::Bool) where {T <: iTgbm}

  if istip(tree) 
    if isalive(tree) && !isfix(tree)
      λt   = lλ(tree)[end]
      fdti = fdt(tree)
      ix = true
    end

    return ix, λt, fdti
  end

  if !ix
    ix, λt, fdti = fλ1(tree.d1::T, λt, fdti, ix)
  end
  if !ix
    ix, λt, fdti = fλ1(tree.d2::T, λt, fdti, ix)
  end

  return ix, λt, fdti
end




"""
    lvupdate!(Ψp    ::iTgbmce,
              Ψc    ::iTgbmce,
              llc   ::Float64, 
              bbλp  ::Array{Array{Float64,1},1}, 
              bbλc  ::Array{Array{Float64,1},1}, 
              tsv   ::Array{Array{Float64,1},1},
              pr    ::Int64,
              d1    ::Int64,
              d2    ::Int64,
              α     ::Float64, 
              σλ    ::Float64, 
              μ     ::Float64,
              δt    ::Float64, 
              srδt  ::Float64, 
              lλmxpr::Float64,
              icr   ::Bool,
              wbc   ::Int64,
              dri   ::BitArray{1},
              ldr   ::Int64,
              ter   ::BitArray{1},
              ix    ::Int64)

Make a gbm update for speciation and extinction for a fixed triad.
"""
function lvupdate!(Ψp    ::iTgbmce,
                   Ψc    ::iTgbmce,
                   llc   ::Float64, 
                   bbλp  ::Array{Array{Float64,1},1}, 
                   bbλc  ::Array{Array{Float64,1},1}, 
                   tsv   ::Array{Array{Float64,1},1},
                   pr    ::Int64,
                   d1    ::Int64,
                   d2    ::Int64,
                   α     ::Float64, 
                   σλ    ::Float64, 
                   μ     ::Float64,
                   δt    ::Float64, 
                   srδt  ::Float64, 
                   lλmxpr::Float64,
                   icr   ::Bool,
                   wbc   ::Int64,
                   dri   ::BitArray{1},
                   ldr   ::Int64,
                   ter   ::BitArray{1},
                   ix    ::Int64)

  if ix === ldr 
    # if root
    if ldr === 0
      llc = 
        triad_lupdate_root!(Ψp, Ψc, bbλp, bbλc, 
          tsv, llc, pr, d1, d2, α, σλ, μ, δt, srδt, lλmxpr, icr)
    else
      llc = 
        triad_lvupdate_trio!(Ψp, Ψc, bbλp, bbλc, 
          tsv, llc, pr, d1, d2, α, σλ, μ, δt, srδt, ter, icr, wbc)

    end
  elseif ix < ldr

    ifx1 = isfix(Ψc.d1)
    if ifx1 && isfix(Ψc.d2)
      ix += 1
      if dri[ix]
        llc = 
          lvupdate!(Ψp.d1, Ψc.d1, llc, 
            bbλp, bbλc, tsv, pr, d1, d2, α, σλ, μ, δt, srδt, 
            lλmxpr, icr, wbc, dri, ldr, ter, ix)
      else
        llc = 
          lvupdate!(Ψp.d2, Ψc.d2, llc, 
            bbλp, bbλc, tsv, pr, d1, d2, α, σλ, μ, δt, srδt, 
            lλmxpr, icr, wbc, dri, ldr, ter, ix)
      end
    elseif ifx1
      llc = 
        lvupdate!(Ψp.d1, Ψc.d1, llc, 
          bbλp, bbλc, tsv, pr, d1, d2, α, σλ, μ, δt, srδt, 
          lλmxpr, icr, wbc, dri, ldr, ter, ix)
    else
      llc = 
        lvupdate!(Ψp.d2, Ψc.d2, llc, 
          bbλp, bbλc, tsv, pr, d1, d2, α, σλ, μ, δt, srδt, 
          lλmxpr, icr, wbc, dri, ldr, ter, ix)
    end
  end

  return llc
end






"""
    update_μ!(μc   ::Float64,
              Ψ    ::iTgbmce,
              llc  ::Float64,
              μtn  ::Float64,
              lac  ::Float64,
              μmxpr::Float64,
              svf  ::Function)

MCMC update for `σ` with acceptance log.
"""
function update_μ!(μc   ::Float64,
                   Ψ    ::iTgbmce,
                   llc  ::Float64,
                   μtn  ::Float64,
                   lac  ::Float64,
                   μmxpr::Float64,
                   svf  ::Function)

  # parameter proposal
  μp = mulupt(μc, μtn)::Float64

  # log likelihood and prior ratio
  l, ne = treelength_ne(Ψ, 0.0, 0.0)
  llr = ne*(log(μp) - log(μc)) + l*(μc - μp) + svf(Ψ, μp) - svf(Ψ, μc)

  # prior ratio
  prr = μp > μmxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + log(μp/μc))
    μc   = μp
    llc += llr
    lac += 1.0
  end

  return llc, μc, lac
end




"""
    update_μ!(μc    ::Float64,
              Ψ     ::iTgbmce,
              llc   ::Float64,
              μtn   ::Float64,
              μmxpr::Float64,
              svf  ::Function)

MCMC update for `μ`.
"""
function update_μ!(μc    ::Float64,
                   Ψ     ::iTgbmce,
                   llc   ::Float64,
                   μtn   ::Float64,
                   μmxpr::Float64,
                   svf  ::Function)

  # parameter proposal
  μp = mulupt(μc, μtn)::Float64

  # log likelihood and prior ratio
  l, ne = treelength_ne(Ψ, 0.0, 0.0)
  llr = ne*(log(μp) - log(μc)) + l*(μc - μp) + svf(Ψ, μp) - svf(Ψ, μc)

  # prior ratio
  prr = μp > μmxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + log(μp/μc))
    μc   = μp
    llc += llr
  end

  return llc, μc
end

