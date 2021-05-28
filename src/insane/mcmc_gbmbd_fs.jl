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
                 pupdp   ::NTuple{3,Float64} = (0.1,0.4,0.4),
                 ntry    ::Int64             = 2,
                 nlim    ::Int64             = 500,
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
                      pupdp   ::NTuple{3,Float64} = (0.1,0.4,0.4),
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
  triads, terminus, btotriad = make_triads(idf)

  # make survival conditioning function (stem or crown)
  # svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
  #                          stem_prob_surv_cbd

  scalef = makescalef(obj_ar)

  # parameter updates (1: σλ & σμ, 2: gbm, 3: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, ceil(Int64, Float64(n) * pupdp[i]/spup)))
  end

  # burn-in phase
  Ψp, Ψc, llc, prc, σλc, σμc, σλtn, σμtn =
    mcmc_burn_gbmbd(Ψp, Ψc, bbλp, bbμp, bbλc, bbμc, tsv, λa_prior, μa_prior, 
      σλprior, σμprior, nburn, tune_int, σλi, σμi, σλtni, σμtni, 
      δt, srδt, idf, triads, terminus, btotriad, pup, nlim, prints, scalef)

  # mcmc
  R, Ψv =
    mcmc_gbmbd(Ψp, Ψc, llc, prc, σλc, σμc, bbλp, bbμp, bbλc, bbμc, tsv,
      λa_prior, μa_prior, σλprior, σμprior, niter, nthin, 
      σλtn, σμtn, δt, srδt, idf, triads, terminus, btotriad, pup, nlim, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "mu_root"      => 2,
                 "sigma_lambda" => 3,
                 "sigma_mu"     => 4,
                 "n_extinct"    => 5,
                 "tree_length"  => 6))

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
                         btotriad::Array{Int64,1},
                         pup     ::Array{Int64,1},
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

        ltn += 1
        lup += 1.0

        llc, prc, σλc, lλac = 
          update_σ!(σλc, Ψc, llc, prc, σλtn, lλac, δt, srδt, σλprior, lλ)

        llc, prc, σμc, lμac = 
          update_σ!(σμc, Ψc, llc, prc, σμtn, lμac, δt, srδt, σμprior, lμ)

      # gbm update
      elseif pupi === 2

        tix = ceil(Int64,rand()*ntr)

        pr, d1, d2 = triads[tix]

        bi  = idf[pr]
        dri = dr(bi)
        ldr = length(dri)
        ter = terminus[tix]

        llc, prc = 
          lvupdate!(Ψp, Ψc, llc, prc, bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2,
            σλc, σμc, δt, srδt, λa_prior, μa_prior, dri, ldr, ter, 0)

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

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, σλc, σμc, tsv, bbλp, bbμp, bbλc, bbμc, 
              bix, triad, ter, δt, srδt, nlim)

      end

      # tune parameters
      if ltn === tune_int
        σλtn = scalef(σλtn,lλac/lup)
        σμtn = scalef(σμtn,lμac/lup)
        ltn = 0
      end
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
                    btotriad::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    nlim    ::Int64,
                    prints  ::Int64)

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

        llc, prc, σλc = 
          update_σ!(σλc, Ψc, llc, prc, σλtn, δt, srδt, σλprior, lλ)
        llc, prc, σμc  = 
          update_σ!(σμc, Ψc, llc, prc, σμtn, δt, srδt, σμprior, lμ)

      # gbm update
      elseif pupi === 2

        tix = ceil(Int64,rand()*ntr)

        pr, d1, d2 = triads[tix]

        bi  = idf[pr]
        dri = dr(bi)
        ldr = length(dri)
        ter = terminus[tix]

        llc, prc = 
          lvupdate!(Ψp, Ψc, llc, prc, bbλp, bbμp, bbλc, bbμc, tsv, pr, d1, d2,
            σλc, σμc, δt, srδt, λa_prior, μa_prior, dri, ldr, ter, 0)

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

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, σλc, σμc, tsv, bbλp, bbμp, bbλc, bbμc, 
              bix, triad, ter, δt, srδt, nlim)

          llik_gbm(Ψc, σλc, σμc, δt, srδt)

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
        R[lit,6] = σλc
        R[lit,7] = σμc
        R[lit,8] = snen(Ψc)
        R[lit,9] = treelength(Ψc)
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
        tsv  ::Array{Array{Float64,1},1},
        bbλp ::Array{Array{Float64,1},1}, 
        bbμp ::Array{Array{Float64,1},1}, 
        bbλc ::Array{Array{Float64,1},1}, 
        bbμc ::Array{Array{Float64,1},1}, 
        bix  ::Int64,
        triad::Array{Int64,1},
        δt   ::Float64, 
        srδt ::Float64,
        nlim ::Int64)


Forward simulation proposal function for gbm birth-death.
"""
function fsp(Ψp   ::iTgbmbd,
             Ψc   ::iTgbmbd,
             bi   ::iBf,
             llc  ::Float64,
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
             nlim ::Int64)

  tsi  = tsv[bix]
  tl   = lastindex(tsi)
  nsδt = δt - (tsi[tl] - tsi[tl-1])

  t0, ret, λf, μf = fsbi(bi, bbλc[bix][1], bbμc[bix][1], 
    nsδt, σλ, σμ, δt, srδt, nlim)

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
        tsv, d1, d2, σλ, σμ, δt, srδt, dri, ldr, ter, 0)

      # estimate for `llr` and `acr`
      iλ   = 0.6931471805599453 + λf
      acr += λf - bbλc[pr][tl]

      if ter[1]
        if ter[2]
          # if both are terminal
        else
          # if d1 is terminal
          acr += duoldnorm(λf, bbλc[pr][1], bbλc[d2][end], tsv[pr][end],tsv[d2][end], σλ) - 
                 duoldnorm(bbλc[pr][tl], bbλc[pr][1], bbλc[d2][end], tsv[pr][end],tsv[d2][end], σλ) + 
                 duoldnorm(μf, bbμc[pr][end], bbμc[d2][end], tsv[pr][end], tsv[d2][end], σμ) - 
                 duoldnorm(bbμc[pr][tl], bbμc[pr][end], bbμc[d2][end], tsv[pr][end], tsv[d2][end], σμ)
        end
      elseif ter[2]
        # if d2 is terminal
        acr += duoldnorm(λf, bbλc[pr][1], bbλc[d1][end], tsv[pr][end],tsv[d1][end], σλ) - 
               duoldnorm(bbλc[pr][tl], bbλc[pr][1], bbλc[d1][end], tsv[pr][end],tsv[d1][end], σλ) + 
               duoldnorm(μf, bbμc[pr][end], bbμc[d1][end], tsv[pr][end], tsv[d1][end], σμ) - 
               duoldnorm(bbμc[pr][tl], bbμc[pr][end], bbμc[d1][end], tsv[pr][end], tsv[d1][end], σμ)

      else
        # if no terminal branches involved
        acr += trioldnorm(λf, bbλc[pr][1], bbλc[d1][end], bbλc[d2][end], tsv[pr][end], tsv[d1][end], tsv[d2][end], σλ) - 
               trioldnorm(bbλc[pr][tl], bbλc[pr][1], bbλc[d1][end], bbλc[d2][end], tsv[pr][end], tsv[d1][end], tsv[d2][end], σλ) + 
               trioldnorm(μf, bbμc[pr][end], bbμc[d1][end], bbμc[d2][end], tsv[pr][end], tsv[d1][end], tsv[d2][end], σμ) - 
               trioldnorm(bbμc[pr][tl], bbμc[pr][end], bbμc[d1][end], bbμc[d2][end], tsv[pr][end], tsv[d1][end], tsv[d2][end], σμ)
      end

    else
      pr  = bix
      iλ  = 0.0
      llr = 0.0
      acr = 0.0
    end

    # mh ratio
    if -randexp() < acr
      llc += llr + 
             llik_gbm( t0, σλ, σμ, δt, srδt) + iλ - 
             br_ll_gbm(Ψc, σλ, σμ, δt, srδt, dri, ldr, 0)

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
    ldprop!(treep   ::iTgbmbd,
            treec   ::iTgbmbd,
            bbλp    ::Array{Array{Float64,1},1}, 
            bbμp    ::Array{Array{Float64,1},1}, 
            bbλc    ::Array{Array{Float64,1},1}, 
            bbμc    ::Array{Array{Float64,1},1}, 
            tsv     ::Array{Array{Float64,1},1},
            d1      ::Int64,
            d2      ::Int64,
            σλ      ::Float64, 
            σμ      ::Float64, 
            δt      ::Float64, 
            srδt    ::Float64, 
            dri     ::BitArray{1},
            ldr     ::Int64,
            ix      ::Int64)

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
                 d1   ::Int64,
                 d2   ::Int64,
                 σλ   ::Float64, 
                 σμ   ::Float64, 
                 δt   ::Float64, 
                 srδt ::Float64, 
                 dri  ::BitArray{1},
                 ldr  ::Int64,
                 ter  ::BitArray{1},
                 ix   ::Int64)

  if ix === ldr 

    llr, acr = daughters_lprop!(treep::iTgbmbd, treec::iTgbmbd,
        λf, μf, bbλp, bbμp, bbλc, bbμc, tsv, d1, d2, ter, σλ, σμ, δt, srδt)

  elseif ix < ldr

    ifx1 = isfix(treec.d1::iTgbmbd)
    if ifx1 && isfix(treec.d2::iTgbmbd)
      ix += 1
      if dri[ix]
        llr, acr = 
          ldprop!(treep.d1::iTgbmbd, treec.d1::iTgbmbd, λf, μf, bbλp, bbμp, 
            bbλc, bbμc, tsv, d1, d2, σλ, σμ, δt, srδt, dri, ldr, ter, ix)
      else
        llr, acr = 
          ldprop!(treep.d2::iTgbmbd, treec.d2::iTgbmbd, λf, μf, bbλp, bbμp, 
            bbλc, bbμc, tsv, d1, d2, σλ, σμ, δt, srδt, dri, ldr, ter, ix)
      end
    elseif ifx1
      llr, acr = 
        ldprop!(treep.d1::iTgbmbd, treec.d1::iTgbmbd, λf, μf, bbλp, bbμp, 
          bbλc, bbμc, tsv, d1, d2, σλ, σμ, δt, srδt, dri, ldr, ter, ix)
    else
      llr, acr = 
        ldprop!(treep.d2::iTgbmbd, treec.d2::iTgbmbd, λf, μf, bbλp, bbμp, 
          bbλc, bbμc, tsv, d1, d2, σλ, σμ, δt, srδt, dri, ldr, ter, ix)
    end
  end

  return llr, acr
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
              iλ  ::Float64, 
              iμ  ::Float64, 
              nsδt::Float64,
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
  t0, nsp = sim_gbm(ti(bi) - tfb, iλ, iμ, σλ, σμ, δt, srδt, 1, nlim)
  na = snan(t0)

  λf, μf = NaN, NaN

  # if simulation goes extinct
  if iszero(na)
    ret = false
  # if simulation reached the maximum limit of species
  elseif nsp === nlim
    ret = false
  # if one surviving lineage
  elseif isone(na)
    f, λf, μf = fixalive!(t0, NaN, NaN)
  elseif na > 1
    # if terminal branch
    if it(bi)
      ret = false
    # if continue the simulation
    else
      # fix random tip and return end λ(t) and μ(t) 
      λf, μf = fixrtip!(t0, na, NaN, NaN)

      for j in Base.OneTo(na - 1)

        # get their final λ and μ to continue forward simulation
        ix, λt, μt = fλμ1(t0, NaN, NaN, 1, 0)

        for i in Base.OneTo(2)

          st0, nsp = sim_gbm(nsδt, tfb, λt, μt, σλ, σμ, δt, srδt, 1, nlim)

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

  return t0, ret, λf, μf
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
    fixrtip!(tree::iTgbmbd, 
             na  ::Int64, 
             λf  ::Float64, 
             μf  ::Float64) 

Fixes the the path for a random non extinct tip.
"""
function fixrtip!(tree::iTgbmbd, 
                  na  ::Int64, 
                  λf  ::Float64, 
                  μf  ::Float64) 

  fix!(tree)

  if !istip(tree)
    if isextinct(tree.d1::iTgbmbd)
      λf, μf = fixrtip!(tree.d2::iTgbmbd, na, λf, μf)
    elseif isextinct(tree.d2::iTgbmbd)
      λf, μf = fixrtip!(tree.d1::iTgbmbd, na, λf, μf)
    else
      na1 = snan(tree.d1::iTgbmbd)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        λf, μf = fixrtip!(tree.d2::iTgbmbd, na - na1, λf, μf)
      else
        λf, μf = fixrtip!(tree.d1::iTgbmbd, na1, λf, μf)
      end
    end
  else
    λf = lλ(tree)[end]
    μf = lμ(tree)[end]
  end

  return λf, μf
end




"""
    fixalive!(tree::iTgbmbd, λf::Float64, μf::Float64)

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::iTgbmbd, λf::Float64, μf::Float64)

  if istip(tree) && !isextinct(tree)
    fix!(tree)
    λf = lλ(tree)[end]
    μf = lμ(tree)[end]
    return true, λf, μf
  end

  if !isnothing(tree.d2)
    f, λf, μf = fixalive!(tree.d2::iTgbmbd, λf, μf)
    if f 
      fix!(tree)
      return true, λf, μf
    end
  end

  if !isnothing(tree.d1)
    f, λf, μf = fixalive!(tree.d1::iTgbmbd, λf, μf)
    if f 
      fix!(tree)
      return true, λf, μf
    end
  end

  return false, λf, μf
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

  if ix === ldr 
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
  σp = mulupt(σc, rand() < 0.3 ? σtn*4.0 : σtn)::Float64

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





