#=

Anagenetic `gbmbd` MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmbd(tree    ::sTbd, 
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
                 pupdp   ::NTuple{4,Float64} = (0.1,0.1,0.2,0.2),
                 ntry    ::Int64             = 2,
                 nlim    ::Int64             = 500,
                 δt      ::Float64           = 1e-2,
                 prints  ::Int64             = 5)

Run insane for `gbmbd`.
"""
function insane_gbmbd(tree    ::sTbd, 
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
                      pupdp   ::NTuple{4,Float64} = (0.1,0.1,0.2,0.2),
                      ntry    ::Int64             = 2,
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

  # make Ψ current and proposal parameters
  Ψc = iTgbmbd(tree, δt, srδt, log(λc), log(μc), αi, σλi, σμi)
  Ψp = deepcopy(Ψc)

  # make fix Ψ directory
  idf = iBffs[]
  bit = BitArray{1}()
  makeiBf!(Ψc, idf, bit)

  # allocate `bb` for each fix branch and their `ts` vectors consisting of 
  # [pe(tree), fdt(tree)]
  bbλc = Array{Float64,1}[]
  bbμc = Array{Float64,1}[]
  tsv  = Array{Float64,1}[]

  makebbv!(Ψc, bbλc, bbμc, tsv)

  bbλp = deepcopy(bbλc)
  bbμp = deepcopy(bbμc)

  # make trios
  triads, terminus, btotriad = make_triads(idf)

  # make survival conditioning function (stem or crown)
  svf = iszero(e(Ψc)) ? cond_surv_crown : cond_surv_stem

  # parameter updates (1: α, 2: σλ, 3: σμ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death gbm"

  # burn-in phase
  Ψp, Ψc, llc, prc, αc, σλc, σμc =
    mcmc_burn_gbmbd(Ψp, Ψc, bbλp, bbμp, bbλc, bbμc, tsv, λa_prior, 
     μa_prior, α_prior, σλ_prior, σμ_prior, nburn, αi, σλi, σμi, δt, srδt, 
      idf, triads, terminus, btotriad, pup, nlim, prints, svf)

  # mcmc
  R, Ψv =
    mcmc_gbmbd(Ψp, Ψc, llc, prc, αc, σλc, σμc, bbλp, bbμp, bbλc, bbμc, tsv,
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, niter, nthin, δt, srδt, 
      idf, triads, terminus, btotriad, pup, nlim, prints, svf)

  pardic = Dict(("lambda_root"  => 1,
                 "mu_root"      => 2,
                 "alpha"        => 3,
                 "sigma_lambda" => 4,
                 "sigma_mu"     => 5,
                 "n_extinct"    => 6))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmbd(Ψp      ::iTgbmbd,
                    Ψc      ::iTgbmbd,
                    bbλp    ::Array{Array{Float64,1},1},
                    bbμp    ::Array{Array{Float64,1},1},
                    bbλc    ::Array{Array{Float64,1},1},
                    bbμc    ::Array{Array{Float64,1},1},
                    tsv     ::Array{Array{Float64,1},1},
                    λa_prior::NTuple{2,Float64},
                    μa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc     ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
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

MCMC burn-in chain for `gbmbd`.
"""
function mcmc_burn_gbmbd(Ψp      ::iTgbmbd,
                         Ψc      ::iTgbmbd,
                         bbλp    ::Array{Array{Float64,1},1},
                         bbμp    ::Array{Array{Float64,1},1},
                         bbλc    ::Array{Array{Float64,1},1},
                         bbμc    ::Array{Array{Float64,1},1},
                         tsv     ::Array{Array{Float64,1},1},
                         λa_prior::NTuple{2,Float64},
                         μa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         σμ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc     ::Float64,
                         σλc     ::Float64,
                         σμc     ::Float64,
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

  llc = llik_gbm(Ψc, αc, σλc, σμc, δt, srδt) + svf(Ψc)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])      + 
        logdinvgamma(σμc^2, σμ_prior[1], σμ_prior[2])      + 
        logdnorm(αc, α_prior[1], α_prior[2]^2)             +
        logdunif(exp(lλ(Ψc)[1]), λa_prior[1], λa_prior[2]) +
        logdunif(exp(lμ(Ψc)[1]), μa_prior[1], μa_prior[2])

  lλmxpr = log(λa_prior[2])
  lμmxpr = log(μa_prior[2])

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

        llc, prc, σλc, σμc = 
          update_σ!(σλc, σμc, αc, Ψc, llc, prc, σλ_prior, σμ_prior)

      elseif pupi === 3

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

        llc = lvupdate!(Ψp, Ψc, llc, bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2,
            αc, σλc, σμc, δt, srδt, lλmxpr, lμmxpr, icr, wbc, dri, ldr, ter, 0)

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

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, αc, σλc, σμc, tsv, bbλp, bbμp, bbλc, bbμc, 
              bix, triad, ter, δt, srδt, nlim, icr, wbc)
      end
    end

    next!(pbar)
  end

  return Ψp, Ψc, llc, prc, αc, σλc, σμc
end




"""
     mcmc_gbmbd(Ψp      ::iTgbmbd,
                Ψc      ::iTgbmbd,
                llc     ::Float64,
                prc     ::Float64,
                αc      ::Float64,
                σλc     ::Float64,
                σμc     ::Float64,
                bbλp    ::Array{Array{Float64,1},1},
                bbμp    ::Array{Array{Float64,1},1},
                bbλc    ::Array{Array{Float64,1},1},
                bbμc    ::Array{Array{Float64,1},1},
                tsv     ::Array{Array{Float64,1},1},
                λa_prior::NTuple{2,Float64},
                α_prior ::NTuple{2,Float64},
                μa_prior::NTuple{2,Float64},
                σλ_prior::NTuple{2,Float64},
                σμ_prior::NTuple{2,Float64},
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

MCMC chain for `gbmbd`.
"""
function mcmc_gbmbd(Ψp      ::iTgbmbd,
                    Ψc      ::iTgbmbd,
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    bbλp    ::Array{Array{Float64,1},1},
                    bbμp    ::Array{Array{Float64,1},1},
                    bbλc    ::Array{Array{Float64,1},1},
                    bbμc    ::Array{Array{Float64,1},1},
                    tsv     ::Array{Array{Float64,1},1},
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    μa_prior::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
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
  lμmxpr = log(μa_prior[2])

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 9)

  # make Ψ vector
  Ψv = iTgbmbd[]

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

        llci = llik_gbm(Ψc, αc, σλc, σμc, δt, srδt) + svf(Ψc)
         if !isapprox(llci, llc, atol = 1e-4)
           @show llci, llc, pupi
           return 
        end

      # gbm update
      elseif pupi === 2

        llc, prc, σλc, σμc = 
          update_σ!(σλc, σμc, αc, Ψc, llc, prc, σλ_prior, σμ_prior)

        llci = llik_gbm(Ψc, αc, σλc, σμc, δt, srδt) + svf(Ψc)
         if !isapprox(llci, llc, atol = 1e-4)
           @show llci, llc, pupi
           return 
        end

      # gbm update
      elseif pupi === 3

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

        llc = lvupdate!(Ψp, Ψc, llc, bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2,
            αc, σλc, σμc, δt, srδt, lλmxpr, lμmxpr, icr, wbc, dri, ldr, ter, 0)

        llci = llik_gbm(Ψc, αc, σλc, σμc, δt, srδt) + svf(Ψc)
         if !isapprox(llci, llc, atol = 1e-4)
           @show llci, llc, pupi
           return 
        end

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

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, αc, σλc, σμc, tsv, bbλp, bbμp, bbλc, bbμc, 
              bix, triad, ter, δt, srδt, nlim, icr, wbc)

        llci = llik_gbm(Ψc, αc, σλc, σμc, δt, srδt) + svf(Ψc)
         if !isapprox(llci, llc, atol = 1e-4)
           @show llci, llc, pupi
           return 
        end
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
        R[lit,5] = exp(lμ(Ψc)[1])
        R[lit,6] = αc
        R[lit,7] = σλc
        R[lit,8] = σμc
        R[lit,9] = snen(Ψc, 0)
        push!(Ψv, deepcopy(Ψc))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ψv
end




"""
    fsp(Ψp   ::iTgbmbd,
        Ψc   ::iTgbmbd,
        bi   ::iBffs,
        llc  ::Float64,
        α    ::Float64, 
        σλ   ::Float64, 
        σμ   ::Float64,
        tsv  ::Array{Array{Float64,1},1},
        bbλp ::Array{Array{Float64,1},1}, 
        bbμp ::Array{Array{Float64,1},1}, 
        bbλc ::Array{Array{Float64,1},1}, 
        bbμc ::Array{Array{Float64,1},1}, 
        bix  ::Int64,
        triad::Array{Int64,1},
        ter  ::BitArray{1},
        δt   ::Float64, 
        srδt ::Float64,
        nlim ::Int64,
        icr  ::Bool, 
        wbc  ::Int64)

Forward simulation proposal function for `gbmbd`.
"""
function fsp(Ψp   ::iTgbmbd,
             Ψc   ::iTgbmbd,
             bi   ::iBffs,
             llc  ::Float64,
             α    ::Float64, 
             σλ   ::Float64, 
             σμ   ::Float64,
             tsv  ::Array{Array{Float64,1},1},
             bbλp ::Array{Array{Float64,1},1}, 
             bbμp ::Array{Array{Float64,1},1}, 
             bbλc ::Array{Array{Float64,1},1}, 
             bbμc ::Array{Array{Float64,1},1}, 
             bix  ::Int64,
             triad::Array{Int64,1},
             ter  ::BitArray{1},
             δt   ::Float64, 
             srδt ::Float64,
             nlim ::Int64,
             icr  ::Bool, 
             wbc  ::Int64)

  t0, ret, λf, λf1, μf, μf1, dft0 = 
    fsbi(bi, bbλc[bix][1], bbμc[bix][1], α, σλ, σμ, δt, srδt, nlim)

  # if retain simulation
  if ret

    # get branch information
    dri = dr(bi)
    ldr = lastindex(dri)
    itb = it(bi)

    if !itb
      # make daughter proposal to be concordant with `t0`
      pr, d1, d2 = triad
      llr, acr = ldprop!(Ψp, Ψc, λf, μf, bbλp, bbμp, bbλc, bbμc, 
        tsv, pr, d1, d2, α, σλ, σμ, icr, wbc, δt, srδt, dri, ldr, ter, 0)

      # lambda proposal and current
      bbλi = bbλc[pr]
      # bbμi = bbμc[pr]
      l    = lastindex(bbλi)
      λmp  = 0.5*(λf1 + λf)
      λmc  = 0.5*(bbλi[l-1] + bbλi[l])
      μmp  = 0.5*(μf1 + μf)
      # μmc  = 0.5*(bbμi[l-1] + bbμi[l])
      nep  = -dft0*(exp(λmp) + exp(μmp))
      # nec  = -dft0*(exp(λmc) + exp(μmc))

     # change last event by speciation for llr
      iλ = λmp + log(dft0) - nep

      # acceptance ratio
      acr += λmp - λmc #+ nec - nep

    else
      pr  = bix
      iλ  = 0.0
      llr = 0.0
      acr = 0.0
    end

    # mh ratio
    if -randexp() < acr 

      llr += llik_gbm( t0, α, σλ, σμ, δt, srδt) + iλ - 
             br_ll_gbm(Ψc, α, σλ, σμ, δt, srδt, dri, ldr, 0)

      if icr && isone(wbc)
        if dri[1]
          llr += cond_surv_stem_p(t0) - 
                 cond_surv_stem(Ψc.d1)
        else
          llr += cond_surv_stem_p(t0) -
                 cond_surv_stem(Ψc.d2)
        end
      elseif iszero(wbc)
        llr += cond_surv_stem_p(t0) -
               cond_surv_stem(Ψc)
      end

      llc += llr

      # copy parent to aid vectors
      gbm_copy_f!(t0, bbλc[pr], bbμc[pr], 0)

      # copy daughters vectors
      if !itb
        pr, d1, d2 = triad

        copyto!(bbλc[d1], bbλp[d1])
        copyto!(bbλc[d2], bbλp[d2])
        copyto!(bbμc[d1], bbμp[d1])
        copyto!(bbμc[d2], bbμp[d2])
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
         iμ  ::Float64, 
         α   ::Float64, 
         σλ  ::Float64, 
         σμ  ::Float64, 
         δt  ::Float64, 
         srδt::Float64,
         nlim::Int64)

Forward `gbmbd` simulation for branch `bi`.
"""
function fsbi(bi  ::iBffs, 
              iλ  ::Float64, 
              iμ  ::Float64, 
              α   ::Float64, 
              σλ  ::Float64, 
              σμ  ::Float64, 
              δt  ::Float64, 
              srδt::Float64,
              nlim::Int64)

  # retain the simulation?
  ret = true

  # times
  tfb = tf(bi)

  # simulate tree
  t0, nsp = _sim_gbmbd(ti(bi) - tfb, iλ, iμ, α, σλ, σμ, δt, srδt, 1, nlim)

  na = snan(t0, 0)

  λf, λf1, μf, μf1, dft0 = NaN, NaN, NaN, NaN, NaN

  # if simulation goes extinct or maximum number of species reached
  if iszero(na) || nsp === nlim
    ret = false
  # if one surviving lineage
  elseif isone(na)
    f, λf, λf1, μf, μf1, dft0 = fixalive!(t0, NaN, NaN, NaN, NaN, NaN)
  elseif na > 1
    # if terminal branch
    if it(bi)
      ret = false
    # if continue the simulation
    else
      # fix random tip and return end λ(t) and μ(t) 
      λf, λf1, μf, μf1, dft0 = fixrtip!(t0, na, NaN, NaN, NaN, NaN, NaN)

      for j in Base.OneTo(na - 1)
        # get their final λ and μ to continue forward simulation
        ix, λt, μt, fdti = fλμ1(t0, NaN, NaN, NaN, false)

        for i in Base.OneTo(2)
          st0, nsp = 
            _sim_gbmbd(max(δt - fdti, 0.0), tfb, λt, μt, α, σλ, σμ, 
                      δt, srδt, 1, nlim)
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

  return t0, ret, λf, λf1, μf, μf1, dft0
end




"""
    ldprop!(treep::iTgbmbd,
            treec::iTgbmbd,
            λf   ::Float64,
            μf   ::Float64,
            bbλp ::Array{Array{Float64,1},1}, 
            bbμp ::Array{Array{Float64,1},1}, 
            bbλc ::Array{Array{Float64,1},1}, 
            bbμc ::Array{Array{Float64,1},1}, 
            tsv  ::Array{Array{Float64,1},1},
            pr   ::Int64,
            d1   ::Int64,
            d2   ::Int64,
            α    ::Float64, 
            σλ   ::Float64, 
            σμ   ::Float64, 
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
function ldprop!(treep::iTgbmbd,
                 treec::iTgbmbd,
                 λf   ::Float64,
                 μf   ::Float64,
                 bbλp ::Array{Array{Float64,1},1}, 
                 bbμp ::Array{Array{Float64,1},1}, 
                 bbλc ::Array{Array{Float64,1},1}, 
                 bbμc ::Array{Array{Float64,1},1}, 
                 tsv  ::Array{Array{Float64,1},1},
                 pr   ::Int64,
                 d1   ::Int64,
                 d2   ::Int64,
                 α    ::Float64, 
                 σλ   ::Float64, 
                 σμ   ::Float64, 
                 icr  ::Bool, 
                 wbc  ::Int64,
                 δt   ::Float64, 
                 srδt ::Float64, 
                 dri  ::BitArray{1},
                 ldr  ::Int64,
                 ter  ::BitArray{1},
                 ix   ::Int64)

  if ix === ldr 

    llr, acr = daughters_lprop!(treep, treec, λf, μf, bbλp, bbμp, bbλc, bbμc, 
      tsv, pr, d1, d2, ter, α, σλ, σμ, icr, wbc, δt, srδt)

  elseif ix < ldr

    ifx1 = isfix(treec.d1)
    if ifx1 && isfix(treec.d2)
      ix += 1
      if dri[ix]
        llr, acr = 
          ldprop!(treep.d1, treec.d1, λf, μf, bbλp, bbμp, 
            bbλc, bbμc, tsv, pr, d1, d2, α, σλ, σμ, icr, wbc, δt, srδt, 
            dri, ldr, ter, ix)
      else
        llr, acr = 
          ldprop!(treep.d2, treec.d2, λf, μf, bbλp, bbμp, 
            bbλc, bbμc, tsv, pr, d1, d2, α, σλ, σμ, icr, wbc, δt, srδt, 
            dri, ldr, ter, ix)
      end
    elseif ifx1
      llr, acr = 
        ldprop!(treep.d1, treec.d1, λf, μf, bbλp, bbμp, 
          bbλc, bbμc, tsv, pr, d1, d2, α, σλ, σμ, icr, wbc, δt, srδt, 
          dri, ldr, ter, ix)
    else
      llr, acr = 
        ldprop!(treep.d2, treec.d2, λf, μf, bbλp, bbμp, 
          bbλc, bbμc, tsv, pr, d1, d2, α, σλ, σμ, icr, wbc, δt, srδt, 
          dri, ldr, ter, ix)
    end
  end

  return llr, acr
end




"""
    add1(tree::iTgbmbd, stree::iTgbmbd, it::Int64, ix::Int64)

Add `stree` to tip in `tree` given by `it` in `tree.d1` order.
"""
function addtotip(tree::iTgbmbd, stree::iTgbmbd, ix::Bool) 

  if istip(tree) 
    if isalive(tree) && !isfix(tree)

      sete!(tree, e(tree) + e(stree))
      setproperty!(tree, :iμ, isextinct(stree))

      lλ0 = lλ(tree)
      lμ0 = lμ(tree)
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
    fixrtip!(tree::iTgbmbd, 
             na  ::Int64, 
             λfm1::Float64,
             λf  ::Float64, 
             μf  ::Float64) 

Fixes the the path for a random non extinct tip.
"""
function fixrtip!(tree::iTgbmbd, 
                  na  ::Int64, 
                  λf  ::Float64, 
                  λf1 ::Float64, 
                  μf  ::Float64,
                  μf1 ::Float64,
                  dft0::Float64) 

  fix!(tree)

  if isdefined(tree, :d1)
    if isextinct(tree.d1)
      λf, λf1, μf, μf1, dft0 = 
        fixrtip!(tree.d2, na, λf, λf1, μf, μf1, dft0)
    elseif isextinct(tree.d2)
      λf, λf1, μf, μf1, dft0 = 
        fixrtip!(tree.d1, na, λf, λf1, μf, μf1, dft0)
    else
      na1 = snan(tree.d1, 0)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        λf, λf1, μf, μf1, dft0 = 
          fixrtip!(tree.d2, na - na1, λf, λf1, μf, μf1, dft0)
      else
        λf, λf1, μf, μf1, dft0 = 
          fixrtip!(tree.d1, na1, λf, λf1, μf, μf1, dft0)
      end
    end
  else

    dft0 = fdt(tree)
    λv   = lλ(tree)
    μv   = lμ(tree)
    l    = lastindex(λv)
    λf   = λv[l]
    λf1  = λv[l-1]
    μf   = μv[l]
    μf1  = μv[l-1]
  end

  return λf, λf1, μf, μf1, dft0
end




"""
    ixalive!(tree::iTgbmbd,
             λf  ::Float64,
             λf1 ::Float64,
             μf  ::Float64,
             μf1 ::Float64,
             dft0::Float64)

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::iTgbmbd,
                   λf  ::Float64,
                   λf1 ::Float64,
                   μf  ::Float64,
                   μf1 ::Float64,
                   dft0::Float64)

  if istip(tree) 
    if isalive(tree)
      fix!(tree)
      dft0 = fdt(tree)
      λv   = lλ(tree)
      μv   = lμ(tree)
      l    = lastindex(λv)
      λf   = λv[l]
      λf1  = λv[l-1]
      μf   = μv[l]
      μf1  = μv[l-1]

      return true, λf, λf1, μf, μf1, dft0
    end
  else
    f, λf, λf1, μf, μf1, dft0 = 
      fixalive!(tree.d2, λf, λf1, μf, μf1, dft0)
    if f 
      fix!(tree)
      return true, λf, λf1, μf, μf1, dft0
    end
    f, λf, λf1, μf, μf1, dft0 = 
      fixalive!(tree.d1, λf, λf1, μf, μf1, dft0)
    if f 
      fix!(tree)
      return true, λf, λf1, μf, μf1, dft0
    end
  end

  return false, λf, λf1, μf, μf1, dft0
end





"""
    fλμ1(tree::iTgbmbd, 
         λt  ::Float64, 
         μt  ::Float64, 
         fdti::Float64,
         ix  ::Bool)

Get end `λ` and `μ` for a tip in `tree` given in `tree.d1` order
not taking into account the fixed tip.
"""
function fλμ1(tree::iTgbmbd, 
              λt  ::Float64, 
              μt  ::Float64, 
              fdti::Float64,
              ix  ::Bool)

  if istip(tree) 
    if isalive(tree) &&!isfix(tree)
      @inbounds begin
        lλv  = lλ(tree)
        l    = lastindex(lλv)
        λt   = lλv[l]
        μt   = lμ(tree)[l]
        fdti = fdt(tree)
      end
      ix = true
    end

    return ix, λt, μt, fdti
  end

  if !ix
    ix, λt, μt, fdti = fλμ1(tree.d1, λt, μt, fdti, ix)
  end
  if !ix
    ix, λt, μt, fdti = fλμ1(tree.d2, λt, μt, fdti, ix)
  end

  return ix, λt, μt, fdti
end




"""
    lvupdate!(Ψp    ::iTgbmbd,
              Ψc    ::iTgbmbd,
              llc   ::Float64, 
              bbλp  ::Array{Array{Float64,1},1}, 
              bbμp  ::Array{Array{Float64,1},1}, 
              bbλc  ::Array{Array{Float64,1},1}, 
              bbμc  ::Array{Array{Float64,1},1}, 
              tsv   ::Array{Array{Float64,1},1},
              pr    ::Int64,
              d1    ::Int64,
              d2    ::Int64,
              α     ::Float64, 
              σλ    ::Float64, 
              σμ    ::Float64, 
              δt    ::Float64, 
              srδt  ::Float64, 
              lλmxpr::Float64,
              lμmxpr::Float64,
              icr   ::Bool,
              wbc   ::Int64,
              dri   ::BitArray{1},
              ldr   ::Int64,
              ter   ::BitArray{1},
              ix    ::Int64)

Make a gbm update for speciation and extinction for a fixed triad.
"""
function lvupdate!(Ψp    ::iTgbmbd,
                   Ψc    ::iTgbmbd,
                   llc   ::Float64, 
                   bbλp  ::Array{Array{Float64,1},1}, 
                   bbμp  ::Array{Array{Float64,1},1}, 
                   bbλc  ::Array{Array{Float64,1},1}, 
                   bbμc  ::Array{Array{Float64,1},1}, 
                   tsv   ::Array{Array{Float64,1},1},
                   pr    ::Int64,
                   d1    ::Int64,
                   d2    ::Int64,
                   α     ::Float64, 
                   σλ    ::Float64, 
                   σμ    ::Float64, 
                   δt    ::Float64, 
                   srδt  ::Float64, 
                   lλmxpr::Float64,
                   lμmxpr::Float64,
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
        triad_lupdate_root!(Ψp, Ψc, bbλp, bbμp, bbλc, bbμc, 
          tsv, llc, pr, d1, d2, α, σλ, σμ, δt, srδt, lλmxpr, lμmxpr, icr)
    else
      llc = 
        triad_lvupdate_trio!(Ψp, Ψc, bbλp, bbμp, bbλc, bbμc, 
          tsv, llc, pr, d1, d2, α, σλ, σμ, δt, srδt, ter, icr, wbc)

    end
  elseif ix < ldr

    ifx1 = isfix(Ψc.d1)
    if ifx1 && isfix(Ψc.d2)
      ix += 1
      if dri[ix]
        llc = 
          lvupdate!(Ψp.d1, Ψc.d1, llc, 
            bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2, α, σλ, σμ, δt, srδt, 
            lλmxpr, lμmxpr, icr, wbc, dri, ldr, ter, ix)
      else
        llc = 
          lvupdate!(Ψp.d2, Ψc.d2, llc, 
            bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2, α, σλ, σμ, δt, srδt, 
            lλmxpr, lμmxpr, icr, wbc, dri, ldr, ter, ix)
      end
    elseif ifx1
      llc = 
        lvupdate!(Ψp.d1, Ψc.d1, llc, 
          bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2, α, σλ, σμ, δt, srδt, 
          lλmxpr, lμmxpr, icr, wbc, dri, ldr, ter, ix)
    else
      llc = 
        lvupdate!(Ψp.d2, Ψc.d2, llc, 
          bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2, α, σλ, σμ, δt, srδt, 
          lλmxpr, lμmxpr, icr, wbc, dri, ldr, ter, ix)
    end
  end

  return llc
end





"""
    update_σ!(σλc     ::Float64,
              σμc     ::Float64,
              α       ::Float64,
              Ψ       ::iTgbmbd,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64},
              σμ_prior::NTuple{2,Float64})

Gibbs update for `σλ` and `σμ`.
"""
function update_σ!(σλc     ::Float64,
                   σμc     ::Float64,
                   α       ::Float64,
                   Ψ       ::iTgbmbd,
                   llc     ::Float64,
                   prc     ::Float64,
                   σλ_prior::NTuple{2,Float64},
                   σμ_prior::NTuple{2,Float64})

  # standardized sum of squares
  sssλ, sssμ, n = sss_gbm(Ψ, α)

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_prior[1] + 0.5 * n, σλ_prior[2] + sssλ)
  σμp2 = randinvgamma(σμ_prior[1] + 0.5 * n, σμ_prior[2] + sssμ)

  # update prior
  prc += llrdinvgamma(σλp2, σλc^2, σλ_prior[1], σλ_prior[2]) + 
         llrdinvgamma(σμp2, σμc^2, σμ_prior[1], σμ_prior[2])

  σλp = sqrt(σλp2)
  σμp = sqrt(σμp2)

  # update likelihood
  llc += sssλ*(1.0/σλc^2 - 1.0/σλp^2) - n*(log(σλp/σλc)) + 
         sssμ*(1.0/σμc^2 - 1.0/σμp^2) - n*(log(σμp/σμc))

  return llc, prc, σλp, σμp
end



