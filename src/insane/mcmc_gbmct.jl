#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmct(tree    ::sTbd, 
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
                 prints  ::Int64             = 5)

Run insane for GBM birth-death.
"""
function insane_gbmct(tree    ::sTbd, 
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
                      prints  ::Int64             = 5)

  # fix tree
  fixtree!(tree)

  # `n` tips, `th` treeheight define δt
  n    = sntn(tree, 0)
  th   = treeheight(tree)
  δt  *= th
  srδt = sqrt(δt)

   # starting parameters (using method of moments)
  if isnan(λi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc = λi
  end
  ϵc = ϵi

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # make Ψ current and proposal parameters
  Ψc = iTgbmct(tree, δt, srδt, log(λc), αi, σλi)
  Ψp = deepcopy(Ψc)

  # make fix Ψ directory
  idf = iBffs[]
  bit = BitArray{1}()
  makeiBf!(Ψc, idf, bit)

  # make survival conditioning function (stem or crown)
  svf = iszero(e(Ψc)) ? cond_surv_crown : cond_surv_stem

  # parameter updates (1: α, 2: σλ, 3: ϵ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death gbm with constant ϵ"

  # burn-in phase
  Ψp, Ψc, llc, prc, αc, σλc, ϵc, ϵtn =
    mcmc_burn_gbmct(Ψp, Ψc,λa_prior, α_prior, σλ_prior, ϵ_prior, 
      nburn, tune_int, αi, σλi, ϵc, ϵtni, δt, srδt, 
      idf, pup, nlim, prints, scalef, svf)

  # mcmc
  R, Ψv = 
    mcmc_gbmct(Ψp, Ψc, llc, prc, αc, σλc, ϵc, ϵtn,
      λa_prior, α_prior, σλ_prior, ϵ_prior, niter, nthin, δt, srδt, 
      idf, pup, nlim, prints, svf)

  pardic = Dict(("lambda_root"   => 1,
                  "alpha"        => 2,
                  "sigma_lambda" => 3,
                  "epsilon"      => 4,
                  "n_extinct"    => 5))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmct(Ψp      ::iTgbmct,
                    Ψc      ::iTgbmct,
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
                    δt      ::Float64,
                    srδt    ::Float64,
                    idf     ::Array{iBffs,1},
                    pup     ::Array{Int64,1},
                    nlim    ::Int64,
                    prints  ::Int64,
                    scalef  ::Function,
                    svf     ::Function)

MCMC burn-in chain for `gbmct`.
"""
function mcmc_burn_gbmct(Ψp      ::iTgbmct,
                         Ψc      ::iTgbmct,
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
                         δt      ::Float64,
                         srδt    ::Float64,
                         idf     ::Array{iBffs,1},
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

  llc = llik_gbm(Ψc, αc, σλc, ϵc, δt, srδt) + svf(Ψc, ϵc)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])      + 
        logdunif(exp(lλ(Ψc)[1]), λa_prior[1], λa_prior[2]) +
        logdnorm(αc, α_prior[1], α_prior[2]^2)             +
        logdunif(ϵc, ϵ_prior[1], ϵ_prior[2])

  lλmxpr = log(λa_prior[2])
  ϵmxpr  = ϵ_prior[2]

  # number of branches
  nbr  = lastindex(idf)

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

        llc, ϵc, lac  = update_ϵ!(ϵc, Ψc, llc, ϵtn, lac, ϵmxpr, svf)

        lup += 1.0

      # gbm update
      elseif pupi === 4

        bi = idf[ceil(Int64,rand()*nbr)]

        # if root
        if iszero(sc(bi)) 
          llc = root_update!(Ψp, Ψc, αc, σλc, ϵc, llc, δt, srδt, lλmxpr, icr)
        elseif sc(bi) === 23
          llc = gbm!(Ψp, Ψc, bi, llc, αc, σλc, ϵc, δt, srδt)
        end

      # forward simulation update
      else

        bi  = idf[ceil(Int64,rand()*nbr)]

        if iszero(sc(bi)) 
          llc = root_update!(Ψp, Ψc, αc, σλc, ϵc, llc, δt, srδt, lλmxpr, icr)
        end

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, αc, σλc, ϵc, δt, srδt, nlim, icr)

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

  return Ψp, Ψc, llc, prc, αc, σλc, ϵc, ϵtn
end




"""
     mcmc_gbmct(Ψp      ::iTgbmct,
                Ψc      ::iTgbmct,
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
                idf     ::Array{iBffs,1},
                pup     ::Array{Int64,1},
                nlim    ::Int64,
                prints  ::Int64,
                svf     ::Function)

MCMC chain for `gbmct`.
"""
function mcmc_gbmct(Ψp      ::iTgbmct,
                    Ψc      ::iTgbmct,
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
                    idf     ::Array{iBffs,1},
                    pup     ::Array{Int64,1},
                    nlim    ::Int64,
                    prints  ::Int64,
                    svf     ::Function)

  # crown or stem conditioning
  icr = iszero(e(Ψc))

  lλmxpr = log(λa_prior[2])
  ϵmxpr  = ϵ_prior[2]

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 8)

  # make Ψ vector
  Ψv = iTgbmct[]

  # number of branches
  nbr  = lastindex(idf)

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    ii = 0
    # parameter updates
    for pupi in pup

      if pupi === 1

        llc, prc, αc  = update_α!(αc, σλc, Ψc, llc, prc, α_prior)

        # ll0 = llik_gbm(Ψc, αc, σλc, ϵc, δt, srδt) + svf(Ψc, ϵc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i,
        #    return 
        # end

      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, Ψc, llc, prc, σλ_prior)

        # ll0 = llik_gbm(Ψc, αc, σλc, ϵc, δt, srδt) + svf(Ψc, ϵc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i,
        #    return 
        # end

      elseif pupi === 3

        llc, ϵc  = update_ϵ!(ϵc, Ψc, llc, ϵtn, ϵmxpr, svf)

        # ll0 = llik_gbm(Ψc, αc, σλc, ϵc, δt, srδt) + svf(Ψc, ϵc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i,
        #    return 
        # end

      # gbm update
      elseif pupi === 4

        bi = idf[ceil(Int64,rand()*nbr)]

        # if root
        if iszero(sc(bi)) 
          llc = root_update!(Ψp, Ψc, αc, σλc, ϵc, llc, δt, srδt, lλmxpr, icr)
        elseif sc(bi) === 23
          llc = gbm!(Ψp, Ψc, bi, llc, αc, σλc, ϵc, δt, srδt)
        end

        # ll0 = llik_gbm(Ψc, αc, σλc, ϵc, δt, srδt) + svf(Ψc, ϵc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i,
        #    return 
        # end

      # forward simulation update
      else


        bi  = idf[ceil(Int64,rand()*nbr)]

        if iszero(sc(bi)) 
          llc = root_update!(Ψp, Ψc, αc, σλc, ϵc, llc, δt, srδt, lλmxpr, icr)
        end

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, αc, σλc, ϵc, δt, srδt, nlim, icr)

        # ll0 = llik_gbm(Ψc, αc, σλc, ϵc, δt, srδt) + svf(Ψc, ϵc)
        #  if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i,
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
        R[lit,7] = ϵc
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
    gbm!(Ψp   ::iTgbmct,
         Ψc   ::iTgbmct,
         bi   ::iBffs,
         llc  ::Float64,
         α    ::Float64,
         σλ   ::Float64,
         μ    ::Float64,
         δt   ::Float64,
         srδt ::Float64)

Update the gbm for branch `bi`.
"""
function gbm!(Ψp   ::iTgbmct,
              Ψc   ::iTgbmct,
              bi   ::iBffs,
              llc  ::Float64,
              α    ::Float64,
              σλ   ::Float64,
              ϵ    ::Float64,
              δt   ::Float64,
              srδt ::Float64)

  # get branch information
  dri = dr(bi)
  ldr = lastindex(dri)
  itb = it(bi)

  # go to branch to be updated
  treec, treep = drtree(Ψc, Ψp, dri, ldr, 0)

  llc = gbm_update!(treep, treec, α, σλ, ϵ, llc, δt, srδt)

  return llc
end




"""
    fsp(Ψp   ::iTgbmct,
        Ψc   ::iTgbmct,
        bi   ::iBffs,
        llc  ::Float64,
        α    ::Float64, 
        σλ   ::Float64, 
        ϵ    ::Float64,
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
             α    ::Float64, 
             σλ   ::Float64, 
             ϵ    ::Float64,
             δt   ::Float64, 
             srδt ::Float64,
             nlim ::Int64,
             icr  ::Bool)

  # get branch information
  dri = dr(bi)
  ldr = lastindex(dri)
  itb = it(bi)

  # go to branch to be updated
  treec, treep = drtree(Ψc, Ψp, dri, ldr, 0)

  t0, ret, λf, dft0 = 
    fsbi_ct(bi, lλ(treec)[1], α, σλ, ϵ, δt, srδt, nlim)

  # if retain simulation
  if ret

    if !itb
      # make daughter proposal to be concordant with `t0`
      llr, acr = daughters_lprop!(treep, treec, λf, α, σλ, ϵ, δt, srδt)

      # speciation
      iλ = λf

      # get previous ending λ
      treecd1, treecd2 = fixds(treec)

      # acceptance ratio
      acr += λf - lλ(treecd1)[1]

    else
      iλ  = 0.0
      llr = 0.0
      acr = 0.0
    end

    # mh ratio
    if -randexp() < acr

      llr += llik_gbm(     t0, α, σλ, ϵ, δt, srδt) + iλ - 
             llik_gbm_f(treec, α, σλ, ϵ, δt, srδt)

      if isone(sc(bi)) && icr
        css = itb ? cond_surv_stem : cond_surv_stem_p
        if dri[1]
          llr += css(t0, ϵ) - cond_surv_stem(treec, ϵ)
        else
          llr += css(t0, ϵ) - cond_surv_stem(treec, ϵ)
        end
      elseif iszero(sc(bi))
        llr += cond_surv_stem_p(t0,  ϵ) -
               cond_surv_stem(treec, ϵ) 
      end

      llc += llr

      # copy daughters vectors
      if !itb
        treepd1, treepd2 = fixds(treep)

        copyto!(lλ(treecd1), lλ(treepd1))
        copyto!(lλ(treecd2), lλ(treepd2))
      end

      # make combined swap branch
      Ψp, Ψc = swapbranch!(Ψp, Ψc, t0, dri, ldr, itb, 0)
    end
  end

  return Ψp, Ψc, llc
end





"""
    fsbi_ct(bi  ::iBffs, 
            iλ  ::Float64, 
            α   ::Float64, 
            σλ  ::Float64, 
            ϵ   ::Float64, 
            δt  ::Float64, 
            srδt::Float64,
            nlim::Int64)

Forward gbmct simulation for branch `bi`.
"""
function fsbi_ct(bi  ::iBffs, 
                 iλ  ::Float64, 
                 α   ::Float64, 
                 σλ  ::Float64, 
                 ϵ   ::Float64, 
                 δt  ::Float64, 
                 srδt::Float64,
                 nlim::Int64)

  # retain the simulation?
  ret = true

  # times
  tfb = tf(bi)

  # simulate tree
  t0, nsp = _sim_gbmct(ti(bi) - tfb, iλ, α, σλ, ϵ, δt, srδt, 1, nlim)

  na = snan(t0, 0)

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
        ix, λt, fdti = fλ1(t0, NaN, NaN, false)

        for i in Base.OneTo(2)
          st0, nsp = 
            _sim_gbmct(max(δt-fdti, 0.0), tfb, λt, α, σλ, ϵ, δt, srδt, 1, nlim)
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
          # if not succeeded after 1 tries.
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

  return t0, ret, λf, dft0
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
function update_ϵ!(ϵc   ::Float64,
                   Ψ    ::iTgbmct,
                   llc  ::Float64,
                   ϵtn  ::Float64,
                   lac  ::Float64,
                   ϵmxpr::Float64,
                   svf  ::Function)

  # parameter proposal
  ϵp = mulupt(ϵc, ϵtn)::Float64

  # log likelihood and prior ratio
  ne   = snenF(Ψ, 0.0)
  ssλt = sλ_gbm(Ψ)
  llr = ne*(log(ϵp) - log(ϵc)) + 
        ssλt*(ϵc - ϵp) + svf(Ψ, ϵp) - svf(Ψ, ϵc)

  # prior ratio
  prr = ϵp > ϵmxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + log(ϵp/ϵc))
    ϵc   = ϵp
    llc += llr
    lac += 1.0
  end

  return llc, ϵc, lac
end




"""
    update_ϵ!(ϵc   ::Float64,
              Ψ    ::iTgbmct,
              llc  ::Float64,
              ϵtn  ::Float64,
              ϵmxpr::Float64,
              svf  ::Function)

MCMC update for `ϵ`.
"""
function update_ϵ!(ϵc   ::Float64,
                   Ψ    ::iTgbmct,
                   llc  ::Float64,
                   ϵtn  ::Float64,
                   ϵmxpr::Float64,
                   svf  ::Function)

  # parameter proposal
  ϵp = mulupt(ϵc, ϵtn)::Float64

  # log likelihood and prior ratio
  ne   = snenF(Ψ, 0.0)
  ssλt = sλ_gbm(Ψ)
  llr = ne*(log(ϵp) - log(ϵc)) + 
        ssλt*(ϵc - ϵp) + svf(Ψ, ϵp) - svf(Ψ, ϵc)

  # prior ratio
  prr = ϵp > ϵmxpr ? -Inf : 0.0

  if -randexp() < (llr + prr + log(ϵp/ϵc))
    ϵc   = ϵp
    llc += llr
  end

  return llc, ϵc
end

