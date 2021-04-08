#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmbd(tree    ::sTbd, 
                 out_file::String;
                 σλprior ::Float64           = 0.1,
                 σμprior ::Float64           = 0.1,
                 λa_prior::NTuple{2,Float64} = (0.0,10.0),
                 μa_prior::NTuple{2,Float64} = (0.0,10.0),
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
                 pupdp   ::NTuple{3,Float64} = (0.8,0.2,0.1),
                 ntry    ::Int64             = 2,
                 δt      ::Float64           = 1e-2,
                 prints  ::Int64             = 5)

Run insane for GBM birth-death.
"""
function insane_gbmbd(tree    ::sTbd, 
                      out_file::String;
                      σλprior ::Float64           = 0.1,
                      σμprior ::Float64           = 0.1,
                      λa_prior::NTuple{2,Float64} = (0.0,10.0),
                      μa_prior::NTuple{2,Float64} = (0.0,10.0),
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
                      pupdp   ::NTuple{3,Float64} = (0.8,0.2,0.1),
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
  else
    λc, μc = λi, μi
  end

  # make Ψ current and proposal parameters
  Ψc = iTgbmbd(tree, δt, srδt, log(λc), log(μc), σλi, σμi)
  Ψp = deepcopy(Ψc)

  # make fix Ψ directory
  idf = iBf[]
  bit = BitArray{1}()
  makeiBf!(Ψc, idf, bit)

  # allocate `bb` for each fix branch and their `ts` vectors
  bbλc = Array{Float64,1}[]
  bbμc = Array{Float64,1}[]
  tsv  = Array{Float64,1}[]

  makebbv!(Ψc, bbλc, bbμc, tsv)

  bbλp = deepcopy(bbλc)
  bbμp = deepcopy(bbμc)

  # make trios
  triads, terminus = make_triads(idf)

  # make survival conditioning function (stem or crown)
  # svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
  #                          stem_prob_surv_cbd

  scalef = makescalef(obj_ar)

  # parameter updates (1: σλ & σμ, 2: gbm, 3: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, floor(Int64, 100.0 * pupdp[i]/spup)))
  end

  # burn-in phase
  Ψp, Ψc, llc, prc, σλc, σμc, σλtn, σμtn =
    mcmc_burn_gbmbd(Ψp, Ψc, bbλp, bbμp, bbλc, bbμc, tsv, λa_prior, μa_prior, 
      σλprior, σμprior, nburn, tune_int, σλi, σμi, σλtni, σμtni, 
      δt, srδt, idf, triads, terminus, pup, ntry, nlim, prints, scalef)

  # mcmc
  R, Ψv =
    mcmc_gbmbd(Ψp, Ψc, llc, prc, σλc, σμc, bbλp, bbμp, bbλc, bbμc, tsv,
      λa_prior, μa_prior, σλprior, σμprior, niter, nthin, 
      σλtn, σμtn, δt, srδt, idf, triads, terminus, pup, ntry, nlim, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "mu_root"      => 2,
                 "sigma_lambda" => 3,
                 "sigma_mu"     => 4,
                 "n_tips"       => 5,
                 "n_extinct"    => 6,
                 "tree_length"  => 7))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmbd(Ψp      ::iTgbmbd,
                    Ψc      ::iTgbmbd,
                    λa_prior::Tuple{Float64,Float64},
                    μa_prior::Tuple{Float64,Float64},
                    σλprior ::Float64,
                    σμprior ::Float64,
                    nburn   ::Int64,
                    tune_int::Int64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    σλtni   ::Float64,
                    σμtni   ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idf     ::Array{iBf,1},
                    triads  ::Array{Array{Int64,1},1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    prints  ::Int64,
                    scalef  ::Function

MCMC burn-in chain for GBM birth-death.
"""
function mcmc_burn_gbmbd(Ψp      ::iTgbmbd,
                         Ψc      ::iTgbmbd,
                         bbλp    ::Array{Array{Float64,1},1},
                         bbμp    ::Array{Array{Float64,1},1},
                         bbλc    ::Array{Array{Float64,1},1},
                         bbμc    ::Array{Array{Float64,1},1},
                         tsv     ::Array{Array{Float64,1},1},
                         λa_prior::Tuple{Float64,Float64},
                         μa_prior::Tuple{Float64,Float64},
                         σλprior ::Float64,
                         σμprior ::Float64,
                         nburn   ::Int64,
                         tune_int::Int64,
                         σλc     ::Float64,
                         σμc     ::Float64,
                         σλtni   ::Float64,
                         σμtni   ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         idf     ::Array{iBf,1},
                         triads  ::Array{Array{Int64,1},1},
                         terminus::Array{BitArray{1}},
                         pup     ::Array{Int64,1},
                         ntry    ::Int64,
                         nlim    ::Int64,
                         prints  ::Int64,
                         scalef  ::Function)

  # initialize acceptance log
  ltn = 0
  lup = lλac = lμac = 0.0

  σλtn = σλtni
  σμtn = σμtni

  llc = llik_gbm(Ψc, σλc, σμc, δt, srδt)
  prc = logdexp(σλc, σλprior)                            +
        logdexp(σμc, σμprior)                            +
        logdnorm_tc(lλ(Ψc)[1], λa_prior[1], λa_prior[2]) +
        logdnorm_tc(lμ(Ψc)[1], μa_prior[1], μa_prior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  # number of branches and of triads
  nbr  = lastindex(idf)
  ntr  = lastindex(triads)

  br = false

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      if pupi === 1
        # update σλ or σμ

        ltn += 1
        lup += 1.0

        llc, prc, σλc, lλac = 
          update_σ!(σλc, Ψc, llc, prc, σλtn, lλac, δt, srδt, σλprior, lλ)

        llc, prc, σμc, lμac = 
          update_σ!(σμc, Ψc, llc, prc, σμtn, lμac, δt, srδt, σμprior, lμ)

      elseif pupi === 2
        # gbm update

        tix = ceil(Int64,rand()*ntr)

        pr, d1, d2 = triads[tix]

        bi  = idf[pr]
        dri = dr(bi)
        ldr = length(dri)
        ter = terminus[tix]

        llc, prc = 
          lvupdate!(Ψp, Ψc, llc, prc, bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2,
            σλc, σμc, δt, srδt, λa_prior, μa_prior, dri, ldr, ter, 0)

      else
        # forward simulation update

        bix   = ceil(Int64,rand()*nbr)
        bi    = idf[bix]
        tsi   = tsv[bix]
        bbiλp = bbλp[bix]
        bbiμp = bbμp[bix]
        bbiλc = bbλc[bix]
        bbiμc = bbμc[bix]

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, σλc, σμc, tsi, bbiλp, bbiμp, bbiλc, bbiμc, 
              δt, srδt, ntry, nlim)

        # llik_gbm(Ψc,  σλc, σμc, δt, srδt)

      end

      # tune parameters
      if ltn == tune_int
        σλtn = scalef(σλtn,lλac/lup)
        σμtn = scalef(σμtn,lμac/lup)
        ltn = 0
      end
    end
    if br
      break
    end

    next!(pbar)
  end

  return Ψp, Ψc, llc, prc, σλc, σμc, σλtn, σμtn
end




"""
     mcmc_gbmbd(Ψp      ::iTgbmbd,
                Ψc      ::iTgbmbd,
                llc     ::Float64,
                prc     ::Float64,
                σλc     ::Float64,
                σμc     ::Float64,
                bbλp    ::Array{Array{Float64,1},1},
                bbμp    ::Array{Array{Float64,1},1},
                bbλc    ::Array{Array{Float64,1},1},
                bbμc    ::Array{Array{Float64,1},1},
                tsv     ::Array{Array{Float64,1},1},
                λa_prior::Tuple{Float64,Float64},
                μa_prior::Tuple{Float64,Float64},
                σλprior ::Float64,
                σμprior ::Float64,
                niter   ::Int64,
                nthin   ::Int64,
                σλtn    ::Float64,
                σμtn    ::Float64,
                δt      ::Float64,
                srδt    ::Float64,
                idf     ::Array{iBf,1},
                triads  ::Array{Array{Int64,1},1},
                terminus::Array{BitArray{1}},
                pup     ::Array{Int64,1},
                ntry    ::Int64,
                nlim    ::Int64,
                prints  ::Int64)

MCMC chain for GBM birth-death.
"""
function mcmc_gbmbd(Ψp      ::iTgbmbd,
                    Ψc      ::iTgbmbd,
                    llc     ::Float64,
                    prc     ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    bbλp    ::Array{Array{Float64,1},1},
                    bbμp    ::Array{Array{Float64,1},1},
                    bbλc    ::Array{Array{Float64,1},1},
                    bbμc    ::Array{Array{Float64,1},1},
                    tsv     ::Array{Array{Float64,1},1},
                    λa_prior::Tuple{Float64,Float64},
                    μa_prior::Tuple{Float64,Float64},
                    σλprior ::Float64,
                    σμprior ::Float64,
                    niter   ::Int64,
                    nthin   ::Int64,
                    σλtn    ::Float64,
                    σμtn    ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idf     ::Array{iBf,1},
                    triads  ::Array{Array{Int64,1},1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    ntry    ::Int64,
                    nlim    ::Int64,
                    prints  ::Int64)

  @info "started mcmc"

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 10)

  # make Ψ vector
  Ψv = iTgbmbd[]

  # number of branches and of triads
  nbr  = lastindex(idf)
  ntr  = lastindex(triads)

  pbar = Progress(niter, prints, "running mcmc...", 20)

  br = false
  for it in Base.OneTo(niter)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      if pupi === 1
        # update σλ or σμ

        llc, prc, σλc = 
          update_σ!(σλc, Ψc, llc, prc, σλtn, δt, srδt, σλprior, lλ)

        llc, prc, σμc  = 
          update_σ!(σμc, Ψc, llc, prc, σμtn, δt, srδt, σμprior, lμ)

      elseif pupi === 2
        # gbm update

        tix = ceil(Int64,rand()*ntr)

        pr, d1, d2 = triads[tix]

        bi  = idf[pr]
        dri = dr(bi)
        ldr = length(dri)
        ter = terminus[tix]

        llc, prc = 
          lvupdate!(Ψp, Ψc, llc, prc, bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2,
            σλc, σμc, δt, srδt, λa_prior, μa_prior, dri, ldr, ter, 0)

      else
        # forward simulation update

        bix   = ceil(Int64,rand()*nbr)
        bi    = idf[bix]
        tsi   = tsv[bix]
        bbiλp = bbλp[bix]
        bbiμp = bbμp[bix]
        bbiλc = bbλc[bix]
        bbiμc = bbμc[bix]

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, σλc, σμc, tsi, bbiλp, bbiμp, bbiλc, bbiμc, 
              δt, srδt, ntry, nlim)
      end
    end

    if br
      break
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
        R[lit,6] = σλc
        R[lit,7] = σμc
        R[lit,8] = sntn(Ψc)
        R[lit,9] = snen(Ψc)
        R[lit,10] = treelength(Ψc)
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
        bi   ::iBf,
        llc  ::Float64,
        σλ   ::Float64, 
        σμ   ::Float64,
        tsi  ::Array{Float64,1},
        bbiλp::Array{Float64,1}, 
        bbiμp::Array{Float64,1}, 
        bbiλc::Array{Float64,1}, 
        bbiμc::Array{Float64,1}, 
        δt   ::Float64, 
        srδt ::Float64,
        ntry ::Int64,
        nlim ::Int64)

Forward simulation proposal function for gbm birth-death.
"""
function fsp(Ψp   ::iTgbmbd,
             Ψc   ::iTgbmbd,
             bi   ::iBf,
             llc  ::Float64,
             σλ   ::Float64, 
             σμ   ::Float64,
             tsi  ::Array{Float64,1},
             bbiλp::Array{Float64,1}, 
             bbiμp::Array{Float64,1}, 
             bbiλc::Array{Float64,1}, 
             bbiμc::Array{Float64,1}, 
             δt   ::Float64, 
             srδt ::Float64,
             ntry ::Int64,
             nlim ::Int64)

  # get branch start and end λ & μ
  dri = dr(bi)
  ldr = lastindex(dri)
  λ0, μ0, λ1, μ1 = λμ01(Ψc, dri, ldr, 0, NaN, NaN)

  # make bb given endpoints
  bb!(bbiλp, λ0, λ1, tsi, σλ, srδt)
  bb!(bbiμp, μ0, μ1, tsi, σμ, srδt)

  # forward simulate a branch
  t0, ret = fsbi(bi, bbiλp, bbiμp, tsi, σλ, σμ, δt, srδt, ntry, nlim)

  # if retain simulation
  if ret

    itb = it(bi)

    # if speciation (if branch is internal)
    iλ = itb ? 0.0 : (log(2.0) + λ1)

    # likelihood ratio
    llr = llik_gbm(t0, σλ, σμ, δt, srδt) + iλ - 
          br_ll_gbm(Ψc, σλ, σμ, δt, srδt, dri, ldr, 0)

    if -randexp() <= 0.0
      llc += llr

      # swap branch
      Ψp = swapbranch!(Ψp, deepcopy(t0), dri, ldr, itb, 0)
      Ψc = swapbranch!(Ψc, t0, dri, ldr, itb, 0)

      copyto!(bbiλc, bbiλp)
      copyto!(bbiμc, bbiμp)
    end
  end

  return Ψp, Ψc, llc
end




"""
    fsbi(bi  ::iBf, 
         bbiλ::Array{Float64,1}, 
         bbiμ::Array{Float64,1}, 
         tsi ::Array{Float64,1},
         σλ  ::Float64, 
         σμ  ::Float64, 
         δt  ::Float64, 
         srδt::Float64,
         ntry::Int64)

Forward gbm birth-death simulation for branch `bi`.
"""
function fsbi(bi  ::iBf, 
              bbiλ::Array{Float64,1}, 
              bbiμ::Array{Float64,1}, 
              tsi ::Array{Float64,1},
              σλ  ::Float64, 
              σμ  ::Float64, 
              δt  ::Float64, 
              srδt::Float64,
              ntry::Int64, 
              nlim::Int64)

  # retain the simulation?
  ret = true

  # times
  tfb = tf(bi)

  # gbm length
  tl = lastindex(tsi)

  # simulate tree
  t0, nsp = sim_ov_gbm(ti(bi) - tfb, 1, tl, bbiλ, bbiμ, tsi, 
    σλ, σμ, δt, srδt, 1, nlim)

  # if fix goes extinct
  if nsp === nlim
    return t0, false
  elseif ifxe(t0)
    return t0, false
  else
    nt = sntn(t0)
    ne = snen(t0)

    # remaining time for last non-standard δt for simulation
    nsδt = δt - (tsi[tl] - tsi[tl-1])

    # ntry per unobserved branch to go extinct
    ii = 0
    for i in Base.OneTo(nt - ne - 1)

      ii += 1

      # get their final λ and μ to continue forward simulation
      ix, λt, μt = fλμ1(t0, NaN, NaN, ii, 0)

      for j in Base.OneTo(ntry)
        st0, nsp = sim_gbm(nsδt, tfb, λt, μt, σλ, σμ, δt, srδt, 1, nlim)

        if nsp === nlim
          continue
        end

        th0 = treeheight(st0)

        # if goes extinct before the present
        if (th0 + 1e-10) < tfb
          # graft to tip
          add1(t0, st0, ii, 0)
          ii -= 1
          break
        end
        # if not succeeded after `ntry` tries. 
        if j === ntry
          ret = false
        end
      end
      # if not a successful simulation after `ntry` tries, 
      # stop for further lineages
      if ret === false
        break
      end
    end
    return t0, ret
  end
end




"""
    add1(tree::iTgbmbd, stree::iTgbmbd, it::Int64, ix::Int64)

Add `stree` to tip in `tree` given by `it` in `tree.d1` order.
"""
function add1(tree::iTgbmbd, stree::iTgbmbd, it::Int64, ix::Int64) 

  if istip(tree) && !isextinct(tree)
    if !isfix(tree)
      ix += 1
    end

    if ix === it
      pet = pe(tree)
      npe = pet + pe(stree)
      setpe!(tree, npe)

      setproperty!(tree, :iμ, stree.iμ)

      sts0 = ts(stree)
      ls   = lastindex(sts0)

      @simd for i in Base.OneTo(ls) 
        sts0[i] += pet
      end

      ts0 = ts(tree)
      lλ0 = lλ(tree)
      lμ0 = lμ(tree)

      pop!(ts0)
      pop!(lλ0)
      pop!(lμ0)

      @views append!(ts0, sts0[2:ls])
      @views append!(lλ0, lλ(stree)[2:ls])
      @views append!(lμ0, lμ(stree)[2:ls])

      tree.d1 = stree.d1
      tree.d2 = stree.d2

      ix += 1
    end

    return ix 
  end

  if ix <= it && !isnothing(tree.d1) 
    ix = add1(tree.d1::iTgbmbd, stree, it, ix)
  end
  if ix <= it && !isnothing(tree.d2) 
    ix = add1(tree.d2::iTgbmbd, stree, it, ix)
  end

  return ix
end




"""
    fλμ1(tree::iTgbmbd, λt::Float64, μt::Float64, it::Int64, ix::Int64)

Get end `λ` and `μ` for a `it` tip in `tree` given in `tree.d1` order
not taking into account the fixed tip.
"""
function fλμ1(tree::iTgbmbd, 
              λt   ::Float64, 
              μt   ::Float64, 
              it   ::Int64, 
              ix   ::Int64)

  if istip(tree) && !isextinct(tree)
    if !isfix(tree)
      ix += 1
    end

    if ix === it
      λt = lλ(tree)[end]
      μt = lμ(tree)[end]
      ix += 1
    end
    return ix, λt, μt
  end

  if ix <= it && !isnothing(tree.d1)
    ix, λt, μt = fλμ1(tree.d1::iTgbmbd, λt, μt, it, ix)
  end
  if ix <= it && !isnothing(tree.d2)
    ix, λt, μt = fλμ1(tree.d2::iTgbmbd, λt, μt, it, ix)
  end

  return ix, λt, μt
end





"""
    lvupdate!(Ψp     ::iTgbmbd,
              Ψc     ::iTgbmbd,
              llc    ::Float64, 
              prc    ::Float64,
              bbiλp  ::Array{Float64,1}, 
              bbiμp  ::Array{Float64,1}, 
              bbiλc  ::Array{Float64,1}, 
              bbiμc  ::Array{Float64,1}, 
              σλ     ::Float64, 
              σμ     ::Float64, 
              δt     ::Float64, 
              srδt   ::Float64, 
              λa_prior::Tuple{Float64,Float64},
              μa_prior::Tuple{Float64,Float64},
              dri    ::BitArray{1},
              ldr    ::Int64,
              ter    ::BitArray{1},
              ix     ::Int64)

Make a gbm update for speciation and extinction for a fixed triad.
"""
function lvupdate!(Ψp      ::iTgbmbd,
                   Ψc      ::iTgbmbd,
                   llc     ::Float64, 
                   prc     ::Float64,
                   bbλp    ::Array{Array{Float64,1},1}, 
                   bbμp    ::Array{Array{Float64,1},1}, 
                   bbλc    ::Array{Array{Float64,1},1}, 
                   bbμc    ::Array{Array{Float64,1},1}, 
                   tsv     ::Array{Array{Float64,1},1},
                   pr      ::Int64,
                   d1      ::Int64,
                   d2      ::Int64,
                   σλ      ::Float64, 
                   σμ      ::Float64, 
                   δt      ::Float64, 
                   srδt    ::Float64, 
                   λa_prior::Tuple{Float64,Float64},
                   μa_prior::Tuple{Float64,Float64},
                   dri     ::BitArray{1},
                   ldr     ::Int64,
                   ter     ::BitArray{1},
                   ix      ::Int64)

  if ix == ldr 
    # if root
    if ldr === 0
      llc, prc = 
        triad_lupdate_root!(Ψp::iTgbmbd, Ψc::iTgbmbd, bbλp, bbμp, bbλc, bbμc, 
              tsv, llc, prc, pr, d1, d2, σλ, σμ, δt, srδt,
              λa_prior, μa_prior)
    else
      if ter[1]
        if ter[2]
          # if both are terminal
          llc = 
            triad_lupdate_noded12!(Ψp::iTgbmbd, Ψc::iTgbmbd,
              bbλp, bbμp, bbλc, bbμc, tsv, llc, pr, d1, d2, σλ, σμ, δt, srδt)
        else
          # if d1 is terminal
          llc = 
            triad_lupdate_noded1!(Ψp::iTgbmbd, Ψc::iTgbmbd,
              bbλp, bbμp, bbλc, bbμc, tsv, llc, pr, d1, d2, σλ, σμ, δt, srδt)
        end
      elseif ter[2]
        # if d2 is terminal
        llc = 
          triad_lupdate_noded2!(Ψp::iTgbmbd, Ψc::iTgbmbd,
            bbλp, bbμp, bbλc, bbμc, tsv, llc, pr, d1, d2, σλ, σμ, δt, srδt)
      else
        # if no terminal branches involved
        llc = 
          triad_lupdate_node!(Ψp::iTgbmbd, Ψc::iTgbmbd,
            bbλp, bbμp, bbλc, bbμc, tsv, llc, pr, d1, d2, σλ, σμ, δt, srδt)
      end
    end

  elseif ix < ldr

    ifx1 = isfix(Ψc.d1::iTgbmbd)
    if ifx1 && isfix(Ψc.d2::iTgbmbd)
      ix += 1
      if dri[ix]
        llc, prc = 
          lvupdate!(Ψp.d1::iTgbmbd, Ψc.d1::iTgbmbd, llc, prc, 
            bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2, σλ, σμ, δt, srδt, 
            λa_prior, μa_prior, dri, ldr, ter, ix)
      else
        llc, prc = 
          lvupdate!(Ψp.d2::iTgbmbd, Ψc.d2::iTgbmbd, llc, prc, 
            bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2, σλ, σμ, δt, srδt, 
            λa_prior, μa_prior, dri, ldr, ter, ix)
      end
    elseif ifx1
      llc, prc = 
        lvupdate!(Ψp.d1::iTgbmbd, Ψc.d1::iTgbmbd, llc, prc, 
          bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2, σλ, σμ, δt, srδt, 
          λa_prior, μa_prior, dri, ldr, ter, ix)
    else
      llc, prc = 
        lvupdate!(Ψp.d2::iTgbmbd, Ψc.d2::iTgbmbd, llc, prc, 
          bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2, σλ, σμ, δt, srδt, 
          λa_prior, μa_prior, dri, ldr, ter, ix)
    end
  end

  return llc, prc
end




"""
    triad_lupdate_noded12!(treep::iTgbmbd, 
                           treec::iTgbmbd,
                           bbλp ::Array{Array{Float64,1},1}, 
                           bbμp ::Array{Array{Float64,1},1}, 
                           bbλc ::Array{Array{Float64,1},1}, 
                           bbμc ::Array{Array{Float64,1},1}, 
                           tsv  ::Array{Array{Float64,1},1}, 
                           llc  ::Float64,
                           pr   ::Int64,
                           d1   ::Int64,
                           d2   ::Int64,
                           σλ   ::Float64,
                           σμ   ::Float64,
                           δt   ::Float64, 
                           srδt ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
both daughters are terminal.
"""
function triad_lupdate_noded12!(treep::iTgbmbd, 
                                treec::iTgbmbd,
                                bbλp ::Array{Array{Float64,1},1}, 
                                bbμp ::Array{Array{Float64,1},1}, 
                                bbλc ::Array{Array{Float64,1},1}, 
                                bbμc ::Array{Array{Float64,1},1}, 
                                tsv  ::Array{Array{Float64,1},1}, 
                                llc  ::Float64,
                                pr   ::Int64,
                                d1   ::Int64,
                                d2   ::Int64,
                                σλ   ::Float64,
                                σμ   ::Float64,
                                δt   ::Float64, 
                                srδt ::Float64)

  # speciation vectors
  λprv_p = bbλp[pr]
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λprv_c = bbλc[pr]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  lipr = lastindex(λprv_c)

  # extinction vectors
  μprv_p = bbμp[pr]
  μd1v_p = bbμp[d1]
  μd2v_p = bbμp[d2]
  μprv_c = bbμc[pr]
  μd1v_c = bbμc[d1]
  μd2v_c = bbμc[d2]

  # get fixed daughters
  treecd1 = fixd1(treec)
  treecd2 = fixd2(treec)
  treepd1 = fixd1(treep)
  treepd2 = fixd2(treep)

  # simulate fix tree vector
  bm!(λprv_p, μprv_p, λprv_c[1], μprv_c[1], tsv[pr], σλ, σμ, srδt)
  lλp = λprv_p[lipr]
  lμp = μprv_p[lipr]
  bm!(λd1v_p, μd1v_p, lλp, lμp, tsv[d1], σλ, σμ, srδt)
  bm!(λd2v_p, μd2v_p, lλp, lμp, tsv[d2], σλ, σμ, srδt)

  # fill fix and simulate unfix tree
  bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
  bm!(treepd1, λd1v_p, μd1v_p, 1, lastindex(λd1v_p), σλ, σμ, srδt)
  bm!(treepd2, λd2v_p, μd2v_p, 1, lastindex(λd2v_p), σλ, σμ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       σλ, σμ, δt, srδt)

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
    copyto!(μprv_c, μprv_p)
    copyto!(μd1v_c, μd1v_p)
    copyto!(μd2v_c, μd2v_p)
    gbm_copy_f!(treec,   treep)
    gbm_copy_f!(treecd1, treepd1)
    gbm_copy_f!(treecd2, treepd2)
  end

  return llc
end




"""
    triad_lupdate_noded1!(treep  ::iTgbmbd, 
                          treec  ::iTgbmbd,
                          bbλp::Array{Array{Float64,1},1}, 
                          bbμp::Array{Array{Float64,1},1}, 
                          bbλc::Array{Array{Float64,1},1}, 
                          bbμc::Array{Array{Float64,1},1}, 
                          tsv ::Array{Array{Float64,1},1}, 
                          llc ::Float64,
                          pr  ::Int64,
                          d1  ::Int64,
                          d2  ::Int64,
                          σλ  ::Float64,
                          σμ  ::Float64,
                          δt  ::Float64, 
                          srδt::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 1 is terminal.
"""
function triad_lupdate_noded1!(treep  ::iTgbmbd, 
                               treec  ::iTgbmbd,
                               bbλp::Array{Array{Float64,1},1}, 
                               bbμp::Array{Array{Float64,1},1}, 
                               bbλc::Array{Array{Float64,1},1}, 
                               bbμc::Array{Array{Float64,1},1}, 
                               tsv ::Array{Array{Float64,1},1}, 
                               llc ::Float64,
                               pr  ::Int64,
                               d1  ::Int64,
                               d2  ::Int64,
                               σλ  ::Float64,
                               σμ  ::Float64,
                               δt  ::Float64, 
                               srδt::Float64)

  ## get reference vectors
  # time vectors
  tprv = tsv[pr]
  td1v = tsv[d1]
  td2v = tsv[d2]
  lipr = lastindex(tprv)
  lid1 = lastindex(td1v)
  lid2 = lastindex(td2v)

  # speciation vectors
  λprv_p = bbλp[pr]
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λprv_c = bbλc[pr]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  λpr    = λprv_c[1]
  λd2    = λd2v_c[lid2]

  # extinction vectors
  μprv_p = bbμp[pr]
  μd1v_p = bbμp[d1]
  μd2v_p = bbμp[d2]
  μprv_c = bbμc[pr]
  μd1v_c = bbμc[d1]
  μd2v_c = bbμc[d2]
  μpr    = μprv_c[1]
  μd2    = μd2v_c[lid2]

  # get fixed daughters
  treecd1 = fixd1(treec)
  treecd2 = fixd2(treec)
  treepd1 = fixd1(treep)
  treepd2 = fixd2(treep)

  # pendant edges
  pepr = tprv[lipr]
  ped1 = td1v[lid1]
  ped2 = td2v[lid2]

  # node proposal
  lλp = duoprop(λpr, λd2, pepr, ped2, σλ)
  lμp = duoprop(μpr, μd2, pepr, ped2, σμ)

  # simulate fix tree vector
  bb!(λprv_p, λpr, lλp, μprv_p, μpr, lμp, tprv, σλ, σμ, srδt)
  bm!(λd1v_p, μd1v_p, lλp, lμp, td1v, σλ, σμ, srδt)
  bb!(λd2v_p, lλp, λd2, μd2v_p, lμp, μd2, td2v, σλ, σμ, srδt)

  # fill fix and simulate unfix tree
  bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
  bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
  bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       σλ, σμ, δt, srδt)

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
    copyto!(μprv_c, μprv_p)
    copyto!(μd1v_c, μd1v_p)
    copyto!(μd2v_c, μd2v_p)
    gbm_copy_f!(treec, treep)
    gbm_copy_f!(treecd1, treepd1)
    gbm_copy_f!(treecd2, treepd2)
  end

  return llc
end




"""
    triad_lupdate_noded2!(treep  ::iTgbmbd, 
                          treec  ::iTgbmbd,
                          bbλp::Array{Array{Float64,1},1}, 
                          bbμp::Array{Array{Float64,1},1}, 
                          bbλc::Array{Array{Float64,1},1}, 
                          bbμc::Array{Array{Float64,1},1}, 
                          tsv ::Array{Array{Float64,1},1}, 
                          llc ::Float64,
                          pr  ::Int64,
                          d1  ::Int64,
                          d2  ::Int64,
                          σλ  ::Float64,
                          σμ  ::Float64,
                          δt  ::Float64, 
                          srδt::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 2 is terminal.
"""
function triad_lupdate_noded2!(treep  ::iTgbmbd, 
                               treec  ::iTgbmbd,
                               bbλp::Array{Array{Float64,1},1}, 
                               bbμp::Array{Array{Float64,1},1}, 
                               bbλc::Array{Array{Float64,1},1}, 
                               bbμc::Array{Array{Float64,1},1}, 
                               tsv ::Array{Array{Float64,1},1}, 
                               llc ::Float64,
                               pr  ::Int64,
                               d1  ::Int64,
                               d2  ::Int64,
                               σλ  ::Float64,
                               σμ  ::Float64,
                               δt  ::Float64, 
                               srδt::Float64)

  ## get reference vectors
  # time vectors
  tprv = tsv[pr]
  td1v = tsv[d1]
  td2v = tsv[d2]
  lipr = lastindex(tprv)
  lid1 = lastindex(td1v)
  lid2 = lastindex(td2v)

  # speciation vectors
  λprv_p = bbλp[pr]
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λprv_c = bbλc[pr]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  λpr    = λprv_c[1]
  λd1    = λd1v_c[lid1]

  # extinction vectors
  μprv_p = bbμp[pr]
  μd1v_p = bbμp[d1]
  μd2v_p = bbμp[d2]
  μprv_c = bbμc[pr]
  μd1v_c = bbμc[d1]
  μd2v_c = bbμc[d2]
  μpr    = μprv_c[1]
  μd1    = μd1v_c[lid1]

  # get fixed daughters
  treecd1 = fixd1(treec)
  treecd2 = fixd2(treec)
  treepd1 = fixd1(treep)
  treepd2 = fixd2(treep)

  # pendant edges
  pepr = tprv[lipr]
  ped1 = td1v[lid1]
  ped2 = td2v[lid2]

  # node proposal
  lλp = duoprop(λpr, λd1, pepr, ped1, σλ)
  lμp = duoprop(μpr, μd1, pepr, ped1, σμ)

  # simulate fix tree vector
  bb!(λprv_p, λpr, lλp, μprv_p, μpr, lμp, tprv, σλ, σμ, srδt)
  bb!(λd1v_p, lλp, λd1, μd1v_p, lμp, μd1, td1v, σλ, σμ, srδt)
  bm!(λd2v_p, μd2v_p, lλp, lμp, td2v, σλ, σμ, srδt)

  # fill fix and simulate unfix tree
  bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
  bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
  bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       σλ, σμ, δt, srδt)

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
    copyto!(μprv_c, μprv_p)
    copyto!(μd1v_c, μd1v_p)
    copyto!(μd2v_c, μd2v_p)
    gbm_copy_f!(treec, treep)
    gbm_copy_f!(treecd1, treepd1)
    gbm_copy_f!(treecd2, treepd2)
  end

  return llc
end




"""
    triad_lupdate_node!(treep  ::iTgbmbd, 
                        treec  ::iTgbmbd,
                        bbλp::Array{Array{Float64,1},1}, 
                        bbμp::Array{Array{Float64,1},1}, 
                        bbλc::Array{Array{Float64,1},1}, 
                        bbμc::Array{Array{Float64,1},1}, 
                        tsv ::Array{Array{Float64,1},1}, 
                        llc ::Float64,
                        pr  ::Int64,
                        d1  ::Int64,
                        d2  ::Int64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        δt  ::Float64, 
                        srδt::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
no daughters are terminal.
"""
function triad_lupdate_node!(treep  ::iTgbmbd, 
                             treec  ::iTgbmbd,
                             bbλp::Array{Array{Float64,1},1}, 
                             bbμp::Array{Array{Float64,1},1}, 
                             bbλc::Array{Array{Float64,1},1}, 
                             bbμc::Array{Array{Float64,1},1}, 
                             tsv ::Array{Array{Float64,1},1}, 
                             llc ::Float64,
                             pr  ::Int64,
                             d1  ::Int64,
                             d2  ::Int64,
                             σλ  ::Float64,
                             σμ  ::Float64,
                             δt  ::Float64, 
                             srδt::Float64)

  ## get reference vectors
  # time vectors
  tprv = tsv[pr]
  td1v = tsv[d1]
  td2v = tsv[d2]
  lipr = lastindex(tprv)
  lid1 = lastindex(td1v)
  lid2 = lastindex(td2v)

  # speciation vectors
  λprv_p = bbλp[pr]
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λprv_c = bbλc[pr]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  λpr    = λprv_c[1]
  λd1    = λd1v_c[lid1]
  λd2    = λd2v_c[lid2]

  # extinction vectors
  μprv_p = bbμp[pr]
  μd1v_p = bbμp[d1]
  μd2v_p = bbμp[d2]
  μprv_c = bbμc[pr]
  μd1v_c = bbμc[d1]
  μd2v_c = bbμc[d2]
  μpr    = μprv_c[1]
  μd1    = μd1v_c[lid1]
  μd2    = μd2v_c[lid2]

  # get fixed daughters
  treecd1 = fixd1(treec)
  treecd2 = fixd2(treec)
  treepd1 = fixd1(treep)
  treepd2 = fixd2(treep)

  # pendant edges
  pepr = tprv[lipr]
  ped1 = td1v[lid1]
  ped2 = td2v[lid2]

  # node proposal
  lλp  = trioprop(λpr, λd1, λd2, pepr, ped1, ped2, σλ)
  lμp  = trioprop(μpr, μd1, μd2, pepr, ped1, ped2, σμ)

  # simulate fix tree vector
  bb!(λprv_p, λpr, lλp, μprv_p, μpr, lμp, tprv, σλ, σμ, srδt)
  bb!(λd1v_p, lλp, λd1, μd1v_p, lμp, μd1, td1v, σλ, σμ, srδt)
  bb!(λd2v_p, lλp, λd2, μd2v_p, lμp, μd2, td2v, σλ, σμ, srδt)

  # fill fix and simulate unfix tree
  bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
  bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
  bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       σλ, σμ, δt, srδt)

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
    copyto!(μprv_c, μprv_p)
    copyto!(μd1v_c, μd1v_p)
    copyto!(μd2v_c, μd2v_p)
    gbm_copy_f!(treec, treep)
    gbm_copy_f!(treecd1, treepd1)
    gbm_copy_f!(treecd2, treepd2)
  end

  return llc
end




"""
    triad_lupdate_root!(treep  ::iTgbmbd, 
                        treec  ::iTgbmbd,
                        bbλp::Array{Array{Float64,1},1}, 
                        bbμp::Array{Array{Float64,1},1}, 
                        bbλc::Array{Array{Float64,1},1}, 
                        bbμc::Array{Array{Float64,1},1}, 
                        tsv ::Array{Array{Float64,1},1}, 
                        llc ::Float64,
                        prc ::Float64,
                        pr  ::Int64,
                        d1  ::Int64,
                        d2  ::Int64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        δt  ::Float64, 
                        srδt::Float64,
                        λa_prior::Tuple{Float64, Float64},
                        μa_prior::Tuple{Float64, Float64})

Make a trio of Brownian motion MCMC updates when the root is involved.
"""
function triad_lupdate_root!(treep  ::iTgbmbd, 
                             treec  ::iTgbmbd,
                             bbλp::Array{Array{Float64,1},1}, 
                             bbμp::Array{Array{Float64,1},1}, 
                             bbλc::Array{Array{Float64,1},1}, 
                             bbμc::Array{Array{Float64,1},1}, 
                             tsv ::Array{Array{Float64,1},1}, 
                             llc ::Float64,
                             prc ::Float64,
                             pr  ::Int64,
                             d1  ::Int64,
                             d2  ::Int64,
                             σλ  ::Float64,
                             σμ  ::Float64,
                             δt  ::Float64, 
                             srδt::Float64,
                             λa_prior::Tuple{Float64, Float64},
                             μa_prior::Tuple{Float64, Float64})

  ## get reference vectors
  # time vectors
  tprv = tsv[pr]
  td1v = tsv[d1]
  td2v = tsv[d2]
  lipr = lastindex(tprv)
  lid1 = lastindex(td1v)
  lid2 = lastindex(td2v)

  # speciation vectors
  λprv_p = bbλp[pr]
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λprv_c = bbλc[pr]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  λpr    = λprv_c[1]
  λd1    = λd1v_c[lid1]
  λd2    = λd2v_c[lid2]

  # extinction vectors
  μprv_p = bbμp[pr]
  μd1v_p = bbμp[d1]
  μd2v_p = bbμp[d2]
  μprv_c = bbμc[pr]
  μd1v_c = bbμc[d1]
  μd2v_c = bbμc[d2]
  μpr    = μprv_c[1]
  μd1    = μd1v_c[lid1]
  μd2    = μd2v_c[lid2]

  # get fixed daughters
  treecd1 = fixd1(treec)
  treecd2 = fixd2(treec)
  treepd1 = fixd1(treep)
  treepd2 = fixd2(treep)

  # pendant edges
  pepr = tprv[lipr]
  ped1 = td1v[lid1]
  ped2 = td2v[lid2]

  # proposal given daughters
  lλp = duoprop(λd1, λd2, ped1, ped2, σλ)
  lμp = duoprop(μd1, μd2, ped1, ped2, σμ)

  # propose for root
  srpepr = sqrt(pepr)
  lλrp = rnorm(lλp, srpepr*σλ)
  lμrp = rnorm(lμp, srpepr*σμ)

  # simulate fix tree vector
  bb!(λprv_p, lλrp, lλp, μprv_p, lμrp, lμp, tprv, σλ, σμ, srδt)
  bb!(λd1v_p, lλp,  λd1, μd1v_p,  lμp, μd1, td1v, σλ, σμ, srδt)
  bb!(λd2v_p, lλp,  λd2, μd2v_p,  lμp, μd2, td2v, σλ, σμ, srδt)

  # fill fix and simulate unfix tree
  bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
  bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
  bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       σλ, σμ, δt, srδt)

  # prior ratio
  prr = llrdnorm_x(lλrp, λpr, λa_prior[1], λa_prior[2]) +
        llrdnorm_x(lμrp, μpr, μa_prior[1], μa_prior[2])

  # acceptance ratio
  acr += prr

  if -randexp() < acr
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
    copyto!(μprv_c, μprv_p)
    copyto!(μd1v_c, μd1v_p)
    copyto!(μd2v_c, μd2v_p)
    gbm_copy_f!(treec, treep)
    gbm_copy_f!(treecd1, treepd1)
    gbm_copy_f!(treecd2, treepd2)
  end

  return llc, prc
end





"""
    llr_propr(treep  ::iTgbmbd,
              treepd1::iTgbmbd,
              treepd2::iTgbmbd,
              treec  ::iTgbmbd,
              treecd1::iTgbmbd,
              treecd2::iTgbmbd,
              σλ     ::Float64,
              σμ     ::Float64,
              δt     ::Float64,
              srδt   ::Float64)

Return the likelihood and proposal ratio for birth-death gbm.
"""
function llr_propr(treep  ::iTgbmbd,
                   treepd1::iTgbmbd,
                   treepd2::iTgbmbd,
                   treec  ::iTgbmbd,
                   treecd1::iTgbmbd,
                   treecd2::iTgbmbd,
                   σλ     ::Float64,
                   σμ     ::Float64,
                   δt     ::Float64,
                   srδt   ::Float64)

  llrbm_pr, llrpb_pr = llr_gbm_sep_f(treep, treec, σλ, σμ, δt, srδt)
  llrbm_d1, llrpb_d1 = llr_gbm_sep_f(treepd1, treecd1, σλ, σμ, δt, srδt)
  llrbm_d2, llrpb_d2 = llr_gbm_sep_f(treepd2, treecd2, σλ, σμ, δt, srδt)

  llr = llrbm_pr + llrpb_pr +
        llrbm_d1 + llrpb_d1 +
        llrbm_d2 + llrpb_d2

  acr = llr - llrbm_pr - llrbm_d1 - llrbm_d2 

  return llr, acr
end




"""
    update_σ!(σc    ::Float64,
              Ψ     ::iTgbmbd,
              llc   ::Float64,
              prc   ::Float64,
              σtn   ::Float64,
              ltn   ::Int64,
              lup   ::Float64,
              lac   ::Float64,
              δt    ::Float64,
              srδt  ::Float64,
              σprior::Float64)

MCMC update for `σ` with acceptance log.
"""
function update_σ!(σc    ::Float64,
                   Ψ     ::iTgbmbd,
                   llc   ::Float64,
                   prc   ::Float64,
                   σtn   ::Float64,
                   lac   ::Float64,
                   δt    ::Float64,
                   srδt  ::Float64,
                   σprior::Float64,
                   lf    ::Function)

  # parameter proposals
  σp = mulupt(σc, σtn)::Float64

  # log likelihood and prior ratio
  llr = llr_gbm_bm(Ψ, σp, σc, srδt, lf)
  prr = llrdexp_x(σp, σc, σprior)

  if -randexp() < (llr + prr + log(σp/σc))
    σc   = σp
    llc += llr
    prc += prr
    lac += 1.0
  end

  return llc, prc, σc, lac
end





"""
    update_σ!(σc    ::Float64,
              Ψ     ::iTgbmbd,
              llc   ::Float64,
              prc   ::Float64,
              σtn   ::Float64,
              δt    ::Float64,
              srδt  ::Float64,
              σprior::Float64,
              lf    ::Function)


MCMC update for `σ`.
"""
function update_σ!(σc    ::Float64, 
                   Ψ     ::iTgbmbd,
                   llc   ::Float64,
                   prc   ::Float64,
                   σtn   ::Float64,
                   δt    ::Float64,
                   srδt  ::Float64,
                   σprior::Float64,
                   lf    ::Function)

  # parameter proposals
  σp = mulupt(σc, σtn)::Float64

  # log likelihood and prior ratio
  llr = llr_gbm_bm(Ψ, σp, σc, srδt, lf)
  prr = llrdexp_x(σp, σc, σprior)

  if -randexp() < (llr + prr + log(σp/σc))
    σc   = σp
    llc += llr
    prc += prr
  end

  return llc, prc, σc
end





