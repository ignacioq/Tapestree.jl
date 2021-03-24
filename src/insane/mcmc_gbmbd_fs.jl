#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#

function insane_gbmbd(tree    ::sTbd, 
                      out_file::String;
                      λprior  ::Float64           = 0.1,
                      μprior  ::Float64           = 0.1,
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      tune_int::Int64             = 100,
                      ϵi      ::Float64           = 0.2,
                      λi      ::Float64           = NaN,
                      μi      ::Float64           = NaN,
                      λtni    ::Float64           = 1.0,
                      μtni    ::Float64           = 1.0,
                      obj_ar  ::Float64           = 0.4,
                      pupdp   ::NTuple{3,Float64} = (0.8,0.2,0.1),
                      prints  ::Int64              = 5)

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
  tsv = Array{Float64,1}[]

  makebbv!(Ψc, bbλc, bbμc, tsv)

  bbλp = deepcopy(bbλc)
  bbμp = deepcopy(bbμc)

  # make trios
  triads, terminus = make_triads(idf)


  # make survival conditioning function (stem or crown)
  # svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
  #                          stem_prob_surv_cbd

  scalef = makescalef(obj_ar)

  # parameter updates (1: σλ & σμ, 2: gbm, 3: forward simulation,)
  pup = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, Int64(100.0 * pupdp[i])))
  end

  # initialize acceptance log
  lλtn = 0
  lλup = lλac = 0.0
  σλtn = σλtni
  lμtn = 0
  lμup = lμac = 0.0
  σμtn = σμtni

  # starting parameters
  σλc = σλi
  σμc = σμi

  llc = llik_gbm(Ψc, σλc, σμc, δt, srδt)
  prc = logdexp(σλc, σλprior)                            +
        logdexp(σμc, σμprior)                            +
        logdnorm_tc(lλ(Ψc)[1], λa_prior[1], λa_prior[2]) +
        logdnorm_tc(lμ(Ψc)[1], μa_prior[1], μa_prior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  # number of branches
  nbr  = lastindex(idf)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup


      ## parameter updates
      # update σλ or σμ
      if pupi === 1

        llc, prc, σλc, lλtn, lλup, lλac = 
          update_σ!(σλc, Ψc, llc, prc, 
            σλtn, lλtn, lλup, lλac, δt, srδt, σλprior, lλ)

        llc, prc, σμc, lμtn, lμup, lμac = 
          update_σ!(σμc, Ψc, llc, prc, 
            σμtn, lμtn, lμup, lμac, δt, srδt, σμprior, lμ)

      # gbm update
      elseif pupi === 2



        """
        here! check if only GBM or 3 branch forward simulation update
        or both!
        """




        # ter  = terminus[pupi]
        # inod = inodes[pupi]

        # dri = dr(idv[inod])
        # ldr = lastindex(dri)

        # llc, prc = 
        #   lvupdate!(Ψp, Ψc, llc, prc, σλc, δt, srδt, 
        #     λa_prior, dri, ldr, ter, 0)


      # forward simulation update
      else
        bix   = ceil(Int64,rand()*nbr)
        bi    = idf[bix]
        tsi   = tsv[bix]
        bbiλp = bbλp[bix]
        bbiμp = bbμp[bix]
        bbiλc = bbλc[bix]
        bbiμc = bbμc[bix]

        Ψc, llc = 
          fsp(Ψc, bi, llc, σλc, σμc, tsi, bbiλp, bbiμp, bbiλc, bbiμc, 
              δt, srδt, ntry)
      end


      if ltn == tune_int
        σλtn = scalef(σλtn,lac/lup)
        ltn = 0
      end

    end

    next!(pbar)
  end

  return llc
end











"""
    fsp(tree::sTbd,
        bi  ::iBf,
        llc ::Float64,
        λc  ::Float64, 
        μc  ::Float64,
        ntry::Int64)

Forward simulation proposal function for gbm birth-death.
"""
function fsp(Ψc   ::iTgbmbd,
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
             ntry ::Int64)

  # get branch start and end λ & μ
  dri = dr(bi)
  ldr = lastindex(dri)
  λ0, μ0, λ1, μ1 = λμ01(Ψc, dri, ldr, 0, NaN, NaN)

  # make bb given endpoints
  bb!(bbiλp, λ0, λ1, tsi, σλ, srδt)
  bb!(bbiμp, μ0, μ1, tsi, σμ, srδt)

  # forward simulate a branch
  t0, ret = fsbi(bi, bbiλp, bbiμp, tsi, σλ, σμ, δt, srδt, ntry)

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
      Ψc = swapbranch!(Ψc, t0, dri, ldr, itb, 0)

      copy!(bbiλc, bbiλp)
      copy!(bbiμc, bbiμp)
    end
  end

  return Ψc, llc
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
              ntry::Int64)

  # retain the simulation?
  ret = true

  # times
  tfb = tf(bi)

  # gbm length
  tl = lastindex(tsi)

  # simulate tree
  t0 = sim_ov_gbm(ti(bi) - tfb, 1, tl, bbiλ, bbiμ, tsi, σλ, σμ, δt, srδt)

  # if fix goes extinct
  if ifxe(t0)
    ret = false
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
        st0 = sim_gbm(nsδt, tfb, λt, μt, σλ, σμ, δt, srδt)
        th0 = treeheight(st0)

        # if goes extinct before the present
        if (th0 + 1e-10) < tfb
          # graft to tip
          add1(t0, st0, ii, 0)
          ii -= 1
          break
        end
        if j === ntry
          ret = false
        end
      end
    end
  end

  return t0, ret
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








pr, d1, d2 = triads[3]

bi = idf[pr]
dri = dr(bi)
ldr = length(dri)



"""
    lvupdate!(Ψp      ::iTgbmbd,
              Ψc      ::iTgbmbd,
              llc     ::Float64, 
              prc     ::Float64,
              σ      ::Float64, 
              δt      ::Float64, 
              srδt    ::Float64, 
              a_prior::Tuple{Float64,Float64},
              dri     ::BitArray{1},
              ldr     ::Int64,
              ter     ::BitArray{1},
              ix      ::Int64)

Make a gbm update for a triad.
"""
function lvupdate!(Ψp     ::iTgbmbd,
                   Ψc     ::iTgbmbd,
                   llc    ::Float64, 
                   prc    ::Float64,
                   bbiλp::Array{Float64,1}, 
                   bbiμp::Array{Float64,1}, 
                   bbiλc::Array{Float64,1}, 
                   bbiμc::Array{Float64,1}, 


                   σ      ::Float64, 
                   δt     ::Float64, 
                   srδt   ::Float64, 
                   a_prior::Tuple{Float64,Float64},
                   dri    ::BitArray{1},
                   ldr    ::Int64,
                   ter    ::BitArray{1},
                   ix     ::Int64,
                   lf     ::Function)




  if ix == ldr 
    # if root
    if ldr == 0
      llc, prc = triad_lupdate_root!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                   llc, prc, σ, δt, srδt, a_prior)
    else
      if ter[1]
        if ter[2]
          # if both are terminal




          llc = triad_lupdate_noded12!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                        llc, σ, δt, srδt)
        else
          # if d1 is terminal
          llc = triad_lupdate_noded1!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                       llc, σ, δt, srδt)
        end
      elseif ter[2]
        # if d2 is terminal
        llc = triad_lupdate_noded2!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                     llc, σ, δt, srδt)
      else

        """
        *********
        Make all the rest as respective "copies" for the rest from this one
        *********
        """

        # if no terminal branches involved
        llc = triad_lλupdate_node!(Ψp::iTgbmbd, Ψc::iTgbmbd,
                bbλp, bbμp, bbλc, bbμc, tsv, llc, pr, d1, d2, σλ, δt, srδt)
      end
    end

  elseif ix < ldr
    ix += 1
    if dri[ix]
      llc, prc = 
        lvupdate!(Ψp.d1::iTgbmbd, Ψc.d1::iTgbmbd, 
          llc, prc, σ, δt, srδt, a_prior, dri, ldr, ter, ix)
    else
      llc, prc = 
        lvupdate!(Ψp.d2::iTgbmbd, Ψc.d2::iTgbmbd, 
          llc, prc, σ, δt, srδt, a_prior, dri, ldr, ter, ix)
    end
  end

  return llc, prc
end






"""
    triad_lλupdate_noded12!(Ψc  ::iTgbmpb, 
                            Ψp  ::iTgbmpb,
                            llc ::Float64,
                            σλ  ::Float64, 
                            δt  ::Float64)
                            srδt::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
both daughters are terminal.
"""
function triad_lλupdate_noded12!(Ψp  ::iTgbmpb, 
                                 Ψc  ::iTgbmpb,
                                 llc ::Float64,
                   bbλp::Array{Array{Float64,1},1}, 
                   bbμp::Array{Array{Float64,1},1}, 
                   bbλc::Array{Array{Float64,1},1}, 
                   bbμc::Array{Array{Float64,1},1}, 
                   tsi  ::Array{Array{Float64,1},1}, 




                                 σλ  ::Float64, 
                                 δt  ::Float64,
                                 srδt::Float64,
                                 lf  ::Function)

  # speciation vectors
  λprv_p = bbλp[pr]
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λprv_c = bbλc[pr]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  λpr    = λprv_c[1]

  # extinction vectors
  μprv_p = bbμp[pr]
  μd1v_p = bbμp[d1]
  μd2v_p = bbμp[d2]
  μprv_c = bbμc[pr]
  μd1v_c = bbμc[d1]
  μd2v_c = bbμc[d2]
  μpr    = μprv_c[1]

  # time vectors
  tprv = tsv[pr]
  td1v = tsv[d1]
  td2v = tsv[d2]


  # fill with Brownian motion
  bm!(λprv_p, λpr, tprv, σλ, srδt)
  lλp = λprv_p[end]
  bm!(λd1v_p, lλp, td1v, σλ, srδt)
  bm!(λd2v_p, lλp, td2v, σλ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_noded1!(Ψc::iTgbmpb, 
                           Ψp::iTgbmpb,
                           llc  ::Float64,
                           σλ   ::Float64, 
                           δt   ::Float64)
                           srδt ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 1 is terminal.
"""
function triad_lλupdate_noded1!(Ψp  ::iTgbmpb, 
                                Ψc  ::iTgbmpb,
                                llc ::Float64,
                                σλ  ::Float64, 
                                δt  ::Float64,
                                srδt::Float64)

  # get reference vectors
  λprv_p = lλ(Ψp)
  λd1v_p = lλ(Ψp.d1)
  λd2v_p = lλ(Ψp.d2)
  λprv_c = lλ(Ψc)
  λd1v_c = lλ(Ψc.d1)
  λd2v_c = lλ(Ψc.d2)
  λpr = λprv_c[1]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(Ψc)
  td1v = ts(Ψc.d1)
  td2v = ts(Ψc.d2)

  # node proposal
  lλp = duoprop(λpr, λd2, pe(Ψc), pe(Ψc.d2), σλ)

  bb!(λprv_p, λpr, lλp, tprv, σλ, srδt)
  bm!(λd1v_p, lλp, td1v, σλ, srδt)
  bb!(λd2v_p, lλp, λd2, td2v, σλ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_noded2!(Ψc::iTgbmpb, 
                           Ψp::iTgbmpb,
                           llc  ::Float64,
                           σλ   ::Float64, 
                           δt   ::Float64,
                           srδt ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 2 is terminal.
"""
function triad_lλupdate_noded2!(Ψp  ::iTgbmpb, 
                                Ψc  ::iTgbmpb,
                                llc ::Float64,
                                σλ  ::Float64, 
                                δt  ::Float64,
                                srδt::Float64)

  # get reference vectors
  λprv_p = lλ(Ψp)
  λd1v_p = lλ(Ψp.d1)
  λd2v_p = lλ(Ψp.d2)
  λprv_c = lλ(Ψc)
  λd1v_c = lλ(Ψc.d1)
  λd2v_c = lλ(Ψc.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]

  # time vectors
  tprv = ts(Ψc)
  td1v = ts(Ψc.d1)
  td2v = ts(Ψc.d2)

  # node proposal
  lλp = duoprop(λpr, λd1, pe(Ψc), pe(Ψc.d1), σλ)

  bb!(λprv_p, λpr, lλp, tprv, σλ, srδt)
  bb!(λd1v_p, lλp, λd1, td1v, σλ, srδt)
  bm!(λd2v_p, lλp, td2v, σλ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  if -randexp() < acr 
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end






"""
    triad_lλupdate_node!(Ψp  ::iTgbmpb, 
                         Ψc  ::iTgbmpb,
                         llc ::Float64,
                         σλ  ::Float64,
                         δt  ::Float64, 
                         srδt::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
no daughters are terminal.
"""
function triad_lλupdate_node!(treep  ::iTgbmbd, 
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

  # node proposal
  pepr = pe(treec)
  ped1 = pe(treecd1) 
  ped2 = pe(treecd1)
  lλp = trioprop(λpr, λd1, λd2, pepr, ped1, ped2, σλ)
  lμp = trioprop(μpr, μd1, μd2, pepr, ped1, ped2, σμ)

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
    triad_lλupdate_root!(Ψc      ::iTgbmpb, 
                         Ψp      ::iTgbmpb,
                         llc     ::Float64,
                         prc     ::Float64,
                         σλ      ::Float64, 
                         δt      ::Float64,
                         σλ      ::Float64, 
                         δt      ::Float64,
                         srδt    ::Float64,
                         λa_prior::Tuple{Float64, Float64})

Make a trio of Brownian motion MCMC updates when the root is involved.
"""
function triad_lλupdate_root!(Ψp      ::iTgbmpb, 
                              Ψc      ::iTgbmpb,
                              llc     ::Float64,
                              prc     ::Float64,
                              σλ      ::Float64, 
                              δt      ::Float64,
                              srδt    ::Float64,
                              λa_prior::Tuple{Float64, Float64})

  # get reference vectors
  λprv_p = lλ(Ψp)
  λd1v_p = lλ(Ψp.d1)
  λd2v_p = lλ(Ψp.d2)
  λprv_c = lλ(Ψc)
  λd1v_c = lλ(Ψc.d1)
  λd2v_c = lλ(Ψc.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(Ψc)
  td1v = ts(Ψc.d1)
  td2v = ts(Ψc.d2)

  # proposal given daughters
  lλp = duoprop(λd1, λd2, pe(Ψc.d1), pe(Ψc.d2), σλ)

  # propose for root
  lλrp = rnorm(lλp, sqrt(pe(Ψc))*σλ)

  # make Brownian bridge proposals
  bb!(λprv_p, lλrp, lλp, tprv, σλ, srδt)
  bb!(λd1v_p, lλp,  λd1, td1v, σλ, srδt)
  bb!(λd2v_p, lλp,  λd2, td2v, σλ, srδt)

  ## make acceptance ratio 
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  # prior ratio
  prr = llrdnorm_x(lλrp, λpr, λa_prior[1], λa_prior[2])

  # acceptance ratio
  acr += prr

  if -randexp() < acr 
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
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
  #     lλp - lλc

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
                   ltn   ::Int64,
                   lup   ::Float64,
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

  ltn += 1
  lup += 1.0

  if -randexp() < (llr + prr + log(σp/σc))
    σc   = σp
    llc += llr
    prc += prr
    lac += 1.0
  end

  return llc, prc, σc, ltn, lup, lac
end







