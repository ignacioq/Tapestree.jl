#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmce(tree    ::sTbd, 
                 out_file::String;
                 σλ_prior::Float64           = 0.1,
                 σϵ_prior::Float64           = 0.1,
                 λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                 ϵ_prior ::NTuple{2,Float64} = (0.0, 100.0),
                 niter   ::Int64             = 1_000,
                 nthin   ::Int64             = 10,
                 nburn   ::Int64             = 200,
                 tune_int::Int64             = 100,
                 ϵi      ::Float64           = 0.2,
                 λi      ::Float64           = NaN,
                 μi      ::Float64           = NaN,
                 σλi     ::Float64           = 0.01, 
                 σμi     ::Float64           = 0.01,
                 σλtni   ::Float64           = 1.0,
                 σμtni   ::Float64           = 1.0,
                 obj_ar  ::Float64           = 0.4,
                 pupdp   ::NTuple{3,Float64} = (0.3,0.1,0.1),
                 ntry    ::Int64             = 2,
                 nlim    ::Int64             = 500,
                 δt      ::Float64           = 1e-2,
                 prints  ::Int64             = 5)

Run insane for GBM birth-death.
"""
function insane_gbmce(tree    ::sTbd, 
                      out_file::String;
                      λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                      ϵ_prior ::NTuple{2,Float64} = (0.0, 100.0),
                      σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      tune_int::Int64             = 100,
                      ϵi      ::Float64           = 0.2,
                      λi      ::Float64           = NaN,
                      μi      ::Float64           = NaN,
                      σλi     ::Float64           = 0.01,
                      μtni    ::Float64           = 1.0, 
                      obj_ar  ::Float64           = 0.234,
                      pupdp   ::NTuple{3,Float64} = (0.3,0.1,0.1),
                      ntry    ::Int64             = 2,
                      nlim    ::Int64             = 500,
                      δt      ::Float64           = 1e-2,
                      prints  ::Int64             = 5)

  # fix tree
  fixtree!(tree)

  # `n` tips, `th` treeheight define δt
  n    = sntn(tree)
  th   = treeheight(tree)
  δt  *= th
  srδt = sqrt(δt)

   # starting parameters (using method of moments)
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
    ϵc     = ϵi
  else
    λc, μc = λi, μi
  end

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # make Ψ current and proposal parameters
  Ψc = iTgbmct(tree, δt, srδt, log(λc), σλi)
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
  svf = iszero(pe(Ψc)) ? cond_alone_events_crown : cond_alone_events_stem

  # parameter updates (1: σλ & σμ, 2: gbm, 3: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # burn-in phase
  Ψp, Ψc, llc, prc, μc, σλc, μtn =
    mcmc_burn_gbmbd(Ψp, Ψc, bbλp, bbλc, tsv, λa_prior, ϵ_prior, 
      σλ_prior, nburn, tune_int, μc, σλi, μtni, δt, srδt, 
      idf, triads, terminus, btotriad, pup, nlim, prints, scalef, svf)

  # mcmc
  R, Ψv =
    mcmc_gbmbd(Ψp, Ψc, llc, prc, μc, σλc, μtn, bbλp, bbλc, tsv,
      λa_prior, ϵ_prior, σλ_prior, niter, nthin, δt, srδt, 
      idf, triads, terminus, btotriad, pup, nlim, prints, svf)

  pardic = Dict(("lambda_root"  => 1,
                 "mu"           => 2,
                 "sigma_lambda" => 3,
                 "n_extinct"    => 5,
                 "tree_length"  => 6))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmbd(Ψp      ::iTgbmct,
                    Ψc      ::iTgbmct,
                    bbλp    ::Array{Array{Float64,1},1},
                    bbλc    ::Array{Array{Float64,1},1},
                    tsv     ::Array{Array{Float64,1},1},
                    λa_prior::NTuple{2,Float64},
                    ϵ_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    tune_int::Int64,
                    μc      ::Float64,
                    σλc     ::Float64,
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

MCMC burn-in chain for GBM birth-death.
"""
function mcmc_burn_gbmbd(Ψp      ::iTgbmct,
                         Ψc      ::iTgbmct,
                         bbλp    ::Array{Array{Float64,1},1},
                         bbλc    ::Array{Array{Float64,1},1},
                         tsv     ::Array{Array{Float64,1},1},
                         λa_prior::NTuple{2,Float64},
                         ϵ_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         tune_int::Int64,
                         ϵc      ::Float64,
                         σλc     ::Float64,
                         ϵtn     ::Float64,
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
  icr = iszero(pe(Ψc))

  llc = llik_gbm(Ψc, ϵc, σλc, δt, srδt) + svf(Ψc, ϵc)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])      + 
        logdunif(exp(lλ(Ψc)[1]), λa_prior[1], λa_prior[2]) +
        logdunif(ϵc, ϵ_prior[1], ϵ_prior[2])

  lλmxpr = log(λa_prior[2])
  ϵmxpr  = log(ϵ_prior[2])

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

      # update σλ or σμ
      if pupi === 1

        llc, prc, σλc = update_σ!(σλc, Ψc, llc, prc, σλ_prior)

        llc, ϵc, lac  = update_ϵ!(ϵc, Ψc, llc, ϵtn, lac, ϵmxpr, svf)

        lup += 1.0

      # gbm update
      elseif pupi === 2

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
            μc, σλc, δt, srδt, lλmxpr, icr, wbc, dri, ldr, ter, 0)

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

        Ψp, Ψc, llc = fsp(Ψp, Ψc, bi, llc, μc, σλc, tsv, bbλp, bbλc, 
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

  return Ψp, Ψc, llc, prc, μc, σλc, μtn
end






"""
     mcmc_gbmbd(Ψp      ::iTgbmct,
                Ψc      ::iTgbmct,
                llc     ::Float64,
                prc     ::Float64,
                μc      ::Float64,
                σλc     ::Float64,
                μtn     ::Float64,
                bbλp    ::Array{Array{Float64,1},1},
                bbλc    ::Array{Array{Float64,1},1},
                tsv     ::Array{Array{Float64,1},1},
                λa_prior::NTuple{2,Float64},
                ϵ_prior ::NTuple{2,Float64},
                σλ_prior::NTuple{2,Float64},
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
                prints  ::Int64)

MCMC chain for GBM birth-death.
"""
function mcmc_gbmbd(Ψp      ::iTgbmct,
                    Ψc      ::iTgbmct,
                    llc     ::Float64,
                    prc     ::Float64,
                    μc      ::Float64,
                    σλc     ::Float64,
                    μtn     ::Float64,
                    bbλp    ::Array{Array{Float64,1},1},
                    bbλc    ::Array{Array{Float64,1},1},
                    tsv     ::Array{Array{Float64,1},1},
                    λa_prior::NTuple{2,Float64},
                    ϵ_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
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
  icr = iszero(pe(Ψc))

  lλmxpr = log(λa_prior[2])
  μmxpr  = log(ϵ_prior[2])

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 8)

  # make Ψ vector
  Ψv = iTgbmct[]

  # number of branches and of triads
  nbr  = lastindex(idf)
  ntr  = lastindex(triads)

  # make empty triad
  emptytriad = Int64[]
  emptyter   = BitArray([])

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    ii = 0
    # parameter updates
    for pupi in pup

      ii += 1
      # update σλ or σμ
      if pupi === 1

        llc, prc, σλc = 
          update_σ!(σλc, Ψc, llc, prc, σλ_prior)

          μxx = μc
          llxx = llc
        llc, μc = 
          update_μ!(μc, Ψc, llc, μtn, μmxpr, svf)

        # ll0 = llik_gbm(Ψc, μc, σλc, δt, srδt) + svf(Ψc, μc)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, ii, 1, μc
        #    return 
        # end

      # gbm update
      elseif pupi === 2

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
            μc, σλc, δt, srδt, lλmxpr, icr, wbc, dri, ldr, ter, 0)

        # ll0 = llik_gbm(Ψc, μc, σλc, δt, srδt) + svf(Ψc, μc)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, 2
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

        Ψp, Ψc, llc = fsp(Ψp, Ψc, bi, llc, μc, σλc, tsv, bbλp, bbλc, 
              bix, triad, ter, δt, srδt, nlim, icr, wbc)

        # ll0 = llik_gbm(Ψc, μc, σλc, δt, srδt) + svf(Ψc, μc)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, i, 3
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
        R[lit,5] = μc
        R[lit,6] = σλc
        R[lit,7] = snen(Ψc)
        R[lit,8] = treelength(Ψc)
        push!(Ψv, deepcopy(Ψc))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ψv
end





"""
    fsp(Ψp   ::iTgbmct,
        Ψc   ::iTgbmct,
        bi   ::iBffs,
        llc  ::Float64,
        μ    ::Float64,
        σλ   ::Float64, 
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
function fsp(Ψp   ::iTgbmct,
             Ψc   ::iTgbmct,
             bi   ::iBffs,
             llc  ::Float64,
             μ    ::Float64,
             σλ   ::Float64, 
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

  t0, ret, λf = 
    fsbi(bi, bbλc[bix][1], μ, σλ, δt, srδt, nlim)

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
        tsv, pr, d1, d2, μ, σλ, icr, wbc, δt, srδt, dri, ldr, ter, 0)

      # change last event by speciation for llr
      iλ = λf

      # acceptance ratio
      # bbλcpr = bbλc[pr]
      # l = lastindex(bbλcpr)
      # acr += bbλcpr[l] - λf

      # acr += cond_alone_events_stem(Ψc, dri, ldr, 0) -
      #        cond_alone_events_stem_λ(t0)

    else
      pr  = bix
      iλ  = 0.0
      llr = 0.0
      acr = 0.0

      # acr += cond_alone_events_stem(Ψc, dri, ldr, 0) -
      #        cond_alone_events_stem(t0)
    end

    cll = 0.0
    if icr && isone(wbc)
      if dri[1]
        cll += cond_alone_events_stem_λ(t0, μ) - 
               cond_alone_events_stem(Ψc.d1::iTgbmct, μ)
      else
        cll += cond_alone_events_stem_λ(t0, μ) -
               cond_alone_events_stem(Ψc.d2::iTgbmct, μ)
      end
    elseif iszero(wbc)
      cll += cond_alone_events_stem_λ(t0, μ) -
             cond_alone_events_stem(Ψc, μ)
    end


    # if iszero(wbc)
    #   acr += cond_alone_events_stem(Ψc)   - 
    #          cond_alone_events_stem_λ(t0)
    # else
    # acr += cond_alone_events_stem_woλ(Ψc, dri, ldr, 0) -
    #        cond_alone_events_stem(t0)
    # end

    # mh ratio
    if -randexp() < acr #+ cll
      llr += llik_gbm( t0, μ, σλ, δt, srδt) + iλ - 
             br_ll_gbm(Ψc, μ, σλ, δt, srδt, dri, ldr, 0)

      # if icr && isone(wbc)
      #   if dri[1]
      #     llr += cond_alone_events_stem(t0) - 
      #            cond_alone_events_stem(Ψc.d1::iTgbmct)
      #   else
      #     llr += cond_alone_events_stem(t0) -
      #            cond_alone_events_stem(Ψc.d2::iTgbmct)
      #   end
      # elseif iszero(wbc)
      #   llr += cond_alone_events_stem(t0) -
      #          cond_alone_events_stem(Ψc)
      # end

      llc += llr + cll

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
    fsbi(bi  ::iBffs, 
         iλ  ::Float64, 
         fdti::Float64,
         μ   ::Float64, 
         σλ  ::Float64, 
         δt  ::Float64, 
         srδt::Float64,
         nlim::Int64)

Forward gbmce simulation for branch `bi`.
"""
function fsbi(bi  ::iBffs, 
              iλ  ::Float64, 
              μ   ::Float64, 
              σλ  ::Float64, 
              δt  ::Float64, 
              srδt::Float64,
              nlim::Int64)

  # retain the simulation?
  ret = true

  # times
  tfb = tf(bi)

  # simulate tree
  t0, nsp = sim_gbmce(ti(bi) - tfb, iλ, μ, σλ, δt, srδt, 1, nlim)

  na = snan(t0)

  λf, dft0 = NaN, NaN

  # if simulation goes extinct or maximum number of species reached
  if iszero(na) || nsp === nlim
    ret = false
  # if one surviving lineage
  elseif isone(na)
    f, λf, dft0 = fixalive!(t0, NaN, NaN)
  elseif na > 1
    # if terminal branch
    if it(bi)
      ret = false
    # if continue the simulation
    else

      # fix random tip and return end λ(t)
      λf, dft0 = fixrtip!(t0, na, NaN, NaN)

      for j in Base.OneTo(na - 1)

        # get their final λ to continue forward simulation
        ix, λt, fdti = fλ1(t0, NaN, NaN, 1, 0)

        for i in Base.OneTo(2)

          st0, nsp = 
            sim_gbmce(max(δt - fdti, 0.0), tfb, λt, μ, σλ, δt, srδt, 1, nlim)

          # if maximum number of species reached.
          if nsp === nlim
            if i === 2
              ret = false
            end
            continue
          end

          # if goes extinct before the present
          th0 = treeheight(st0)
          if (th0 + 1e-11) < tfb
            # graft to tip
            add1(t0, st0, 1, 0)
            break
          end

          # if not succeeded after 2 tries.
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

  return t0, ret, λf
end




"""
    ldprop!(treep::iTgbmct,
            treec::iTgbmct,
            λf   ::Float64,
            bbλp ::Array{Array{Float64,1},1}, 
            bbλc ::Array{Array{Float64,1},1}, 
            tsv  ::Array{Array{Float64,1},1},
            pr   ::Int64,
            d1   ::Int64,
            d2   ::Int64,
            μ    ::Float64, 
            σλ   ::Float64, 
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
function ldprop!(treep::iTgbmct,
                 treec::iTgbmct,
                 λf   ::Float64,
                 bbλp ::Array{Array{Float64,1},1}, 
                 bbλc ::Array{Array{Float64,1},1}, 
                 tsv  ::Array{Array{Float64,1},1},
                 pr   ::Int64,
                 d1   ::Int64,
                 d2   ::Int64,
                 μ    ::Float64, 
                 σλ   ::Float64, 
                 icr  ::Bool, 
                 wbc  ::Int64,
                 δt   ::Float64, 
                 srδt ::Float64, 
                 dri  ::BitArray{1},
                 ldr  ::Int64,
                 ter  ::BitArray{1},
                 ix   ::Int64)

  if ix === ldr 

    llr, acr = daughters_lprop!(treep::iTgbmct, treec::iTgbmct, λf,
      bbλp, bbλc, tsv, pr, d1, d2, ter, μ, σλ, icr, wbc, δt, srδt)

  elseif ix < ldr

    ifx1 = isfix(treec.d1::iTgbmct)
    if ifx1 && isfix(treec.d2::iTgbmct)
      ix += 1
      if dri[ix]
        llr, acr = 
          ldprop!(treep.d1::iTgbmct, treec.d1::iTgbmct, λf, bbλp, bbλc, 
            tsv, pr, d1, d2, μ, σλ, icr, wbc, δt, srδt, dri, ldr, ter, ix)
      else
        llr, acr = 
          ldprop!(treep.d2::iTgbmct, treec.d2::iTgbmct, λf, bbλp, bbλc, 
            tsv, pr, d1, d2, μ, σλ, icr, wbc, δt, srδt, dri, ldr, ter, ix)
      end
    elseif ifx1
      llr, acr = 
        ldprop!(treep.d1::iTgbmct, treec.d1::iTgbmct, λf, bbλp, bbλc, 
          tsv, pr, d1, d2, μ, σλ, icr, wbc, δt, srδt, dri, ldr, ter, ix)
    else
      llr, acr = 
        ldprop!(treep.d2::iTgbmct, treec.d2::iTgbmct, λf, bbλp, bbλc, 
          tsv, pr, d1, d2, μ, σλ, icr, wbc, δt, srδt, dri, ldr, ter, ix)
    end
  end

  return llr, acr
end




"""
    add1(tree::iTgbmct, stree::iTgbmct, it::Int64, ix::Int64)

Add `stree` to tip in `tree` given by `it` in `tree.d1` order.
"""
function add1(tree::iTgbmct, stree::iTgbmct, it::Int64, ix::Int64) 

  if istip(tree) && !isextinct(tree)
    if !isfix(tree)
      ix += 1
    end

    if ix === it
      pet = pe(tree)
      npe = pet + pe(stree)
      setpe!(tree, npe)

      setproperty!(tree, :iμ, isextinct(stree))

      lλ0 = lλ(tree)

      pop!(lλ0)

      lλs = lλ(stree)

      popfirst!(lλs)

      append!(lλ0, lλs)

      if isone(lastindex(lλs))
        setfdt!(tree, fdt(tree) + fdt(stree))
      else
        setfdt!(tree, fdt(stree))
      end

      tree.d1 = stree.d1
      tree.d2 = stree.d2

      ix += 1
    end

    return ix 
  end

  if ix <= it && !isnothing(tree.d1) 
    ix = add1(tree.d1::iTgbmct, stree, it, ix)
  end
  if ix <= it && !isnothing(tree.d2) 
    ix = add1(tree.d2::iTgbmct, stree, it, ix)
  end

  return ix
end




"""
    fixrtip!(tree::iTgbmct, na::Int64, λf::Float64, dft0::Float64)

Fixes the the path for a random non extinct tip and returns final `λ(t)`.
"""
function fixrtip!(tree::iTgbmct, na::Int64, λf::Float64, dft0::Float64) 

  fix!(tree)

  if !istip(tree)
    if isextinct(tree.d1::iTgbmct)
      λf, dft0 = fixrtip!(tree.d2::iTgbmct, na, λf, dft0)
    elseif isextinct(tree.d2::iTgbmct)
      λf, dft0 = fixrtip!(tree.d1::iTgbmct, na, λf, dft0)
    else
      na1 = snan(tree.d1::iTgbmct)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        λf, dft0 = fixrtip!(tree.d2::iTgbmct, na - na1, λf, dft0)
      else
        λf, dft0 = fixrtip!(tree.d1::iTgbmct, na1,      λf, dft0)
      end
    end
  else
    λf   = lλ(tree)[end]
    dft0 = fdt(tree)
  end

  return λf, dft0
end




"""
    fixalive!(tree::iTgbmct, λf::Float64, dft0::Float64)

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::iTgbmct, λf::Float64, dft0::Float64)

  if istip(tree) && !isextinct(tree)
    fix!(tree)
    λf   = lλ(tree)[end]
    dft0 = fdt(tree)
    return true, λf, dft0
  end

  if !isnothing(tree.d2)
    f, λf, dft0 = fixalive!(tree.d2::iTgbmct, λf, dft0)
    if f 
      fix!(tree)
      return true, λf, dft0
    end
  end

  if !isnothing(tree.d1)
    f, λf, dft0 = fixalive!(tree.d1::iTgbmct, λf, dft0)
    if f 
      fix!(tree)
      return true, λf, dft0
    end
  end

  return false, λf, dft0
end





"""
    fλ1(tree::T, λt::Float64, fdti::Float64, it::Int64, ix::Int64) where {T <: iTgbm}

Get final `λ(t)` for a `it` tip in `tree` given in `tree.d1` order
not taking into account the fixed tip.
"""
function fλ1(tree::T, λt::Float64, fdti::Float64, it::Int64, ix::Int64) where {T <: iTgbm}

  if istip(tree) && !isextinct(tree)
    if !isfix(tree)
      ix += 1
    end

    if ix === it
      λt   = lλ(tree)[end]
      fdti = fdt(tree)
      ix += 1
    end
    return ix, λt, fdti
  end

  if ix <= it && !isnothing(tree.d1)
    ix, λt, fdti = fλ1(tree.d1::T, λt, fdti, it, ix)
  end
  if ix <= it && !isnothing(tree.d2)
    ix, λt, fdti = fλ1(tree.d2::T, λt, fdti, it, ix)
  end

  return ix, λt, fdti
end




"""
    lvupdate!(Ψp    ::iTgbmct,
              Ψc    ::iTgbmct,
              llc   ::Float64, 
              bbλp  ::Array{Array{Float64,1},1}, 
              bbλc  ::Array{Array{Float64,1},1}, 
              tsv   ::Array{Array{Float64,1},1},
              pr    ::Int64,
              d1    ::Int64,
              d2    ::Int64,
              σλ    ::Float64, 
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
function lvupdate!(Ψp    ::iTgbmct,
                   Ψc    ::iTgbmct,
                   llc   ::Float64, 
                   bbλp  ::Array{Array{Float64,1},1}, 
                   bbλc  ::Array{Array{Float64,1},1}, 
                   tsv   ::Array{Array{Float64,1},1},
                   pr    ::Int64,
                   d1    ::Int64,
                   d2    ::Int64,
                   μ     ::Float64,
                   σλ    ::Float64, 
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
        triad_lupdate_root!(Ψp::iTgbmct, Ψc::iTgbmct, bbλp, bbλc, 
          tsv, llc, pr, d1, d2, μ, σλ, δt, srδt, lλmxpr, icr)
    else
      llc = 
        triad_lvupdate_trio!(Ψp::iTgbmct, Ψc::iTgbmct, bbλp, bbλc, 
          tsv, llc, pr, d1, d2, μ, σλ, δt, srδt, ter, icr, wbc)

    end
  elseif ix < ldr

    ifx1 = isfix(Ψc.d1::iTgbmct)
    if ifx1 && isfix(Ψc.d2::iTgbmct)
      ix += 1
      if dri[ix]
        llc = 
          lvupdate!(Ψp.d1::iTgbmct, Ψc.d1::iTgbmct, llc, 
            bbλp, bbλc, tsv, pr, d1, d2, μ, σλ, δt, srδt, 
            lλmxpr, icr, wbc, dri, ldr, ter, ix)
      else
        llc = 
          lvupdate!(Ψp.d2::iTgbmct, Ψc.d2::iTgbmct, llc, 
            bbλp, bbλc, tsv, pr, d1, d2, μ, σλ, δt, srδt, 
            lλmxpr, icr, wbc, dri, ldr, ter, ix)
      end
    elseif ifx1
      llc = 
        lvupdate!(Ψp.d1::iTgbmct, Ψc.d1::iTgbmct, llc, 
          bbλp, bbλc, tsv, pr, d1, d2, μ, σλ, δt, srδt, 
          lλmxpr, icr, wbc, dri, ldr, ter, ix)
    else
      llc = 
        lvupdate!(Ψp.d2::iTgbmct, Ψc.d2::iTgbmct, llc, 
          bbλp, bbλc, tsv, pr, d1, d2, μ, σλ, δt, srδt, 
          lλmxpr, icr, wbc, dri, ldr, ter, ix)
    end
  end

  return llc
end





"""
    update_σ!(σλc     ::Float64,
              Ψ       ::iTgbmct,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64})

Gibbs update for `σλ`.
"""
function update_σ!(σλc     ::Float64,
                   Ψ       ::iTgbmct,
                   llc     ::Float64,
                   prc     ::Float64,
                   σλ_prior::NTuple{2,Float64})

  # standardized sum of squares
  sssλ, n = sss_gbm(Ψ)

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_prior[1] + 0.5 * n, σλ_prior[2] + sssλ)

  # update prior
  prc += llrdinvgamma(σλp2, σλc^2, σλ_prior[1], σλ_prior[2])

  σλp = sqrt(σλp2)

  # update likelihood
  llc += llr_gbm_σp(σλp, σλc, sssλ, n)

  return llc, prc, σλp
end




"""
    update_ϵ!(ϵc    ::Float64,
              Ψ     ::iTgbmct,
              llc   ::Float64,
              ϵtn   ::Float64,
              lac  ::Float64,
              ϵmxpr::Float64,
              svf  ::Function)

MCMC update for `σ` with acceptance log.
"""
function update_ϵ!(ϵc    ::Float64,
                   Ψ     ::iTgbmct,
                   llc   ::Float64,
                   ϵtn   ::Float64,
                   lac  ::Float64,
                   ϵmxpr::Float64,
                   svf  ::Function)

  # parameter proposal
  ϵn = mulupt(ϵc, ϵtn)::Float64

  # log likelihood and prior ratio
  ne   = snenF(Ψ)
  ssλt = sλ_gbm(Ψ)
  llr = ne*(log(ϵn) - log(ϵc)) + 
        ssλt*(ϵc - ϵn) + svf(Ψ, ϵn) - svf(Ψ, ϵc)

  # prior ratio
  prr = ϵn > ϵmxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + log(ϵn/ϵc))
    ϵc   = ϵn
    llc += llr
    lac += 1.0
  end

  return llc, ϵc, lac
end




"""
    update_ϵ!(ϵc    ::Float64,
              Ψ     ::iTgbmct,
              llc   ::Float64,
              ϵtn   ::Float64,
              ϵmxpr::Float64,
              svf  ::Function)

MCMC update for `ϵ`.
"""
function update_ϵ!(ϵc    ::Float64,
                   Ψ     ::iTgbmct,
                   llc   ::Float64,
                   ϵtn   ::Float64,
                   ϵmxpr::Float64,
                   svf  ::Function)

  # parameter proposal
  ϵn = mulupt(ϵc, ϵtn)::Float64

  # log likelihood and prior ratio
  ne   = snenF(Ψ)
  ssλt = sλ_gbm(Ψ)
  llr = ne*(log(ϵn) - log(ϵc)) + 
        ssλt*(ϵc - ϵn) + svf(Ψ, ϵn) - svf(Ψ, ϵc)

  # prior ratio
  prr = ϵn > ϵmxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + log(ϵn/ϵc))
    ϵc   = ϵn
    llc += llr
  end

  return llc, ϵc
end
