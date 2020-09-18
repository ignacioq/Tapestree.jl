#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#





"""
    insane_gbmpb(tree    ::iTpb, 
                 out_file::String;
                 λa_prior::Tuple{Float64,Float64} = (0.0,10.0),
                 σλprior ::Float64  = 0.1,
                 δt      ::Float64  = 1e-2,
                 niter   ::Int64    = 1_000,
                 nthin   ::Int64    = 10,
                 nburn   ::Int64    = 200,
                 tune_int::Int64    = 100,
                 σλtni   ::Float64  = 1.0,
                 obj_ar  ::Float64  = 0.234,
                 prints  ::Int64    = 5)

Run insane for constant pure-birth.
"""
function insane_gbmpb(tree    ::iTpb, 
                      out_file::String;
                      λa_prior::Tuple{Float64,Float64} = (0.0,10.0),
                      σλprior ::Float64  = 0.1,
                      δt      ::Float64  = 1e-2,
                      niter   ::Int64    = 1_000,
                      nthin   ::Int64    = 10,
                      nburn   ::Int64    = 200,
                      tune_int::Int64    = 100,
                      σλtni   ::Float64  = 1.0,
                      obj_ar  ::Float64  = 0.234,
                      prints  ::Int64    = 5)

  δt  *= treeheight(tree)
  srδt = sqrt(δt)

  # lλ root node
  lλa = log(λmle_cpb(tree))

  # make tree current and proposal parameters
  treec = iTgbmpb(tree, δt, srδt, lλa, 0.1)
  treep = deepcopy(treec)

  # make fix tree directory
  idv = iDir[]
  bit = BitArray{1}()
  makeiDir!(treec, idv, bit)

  # make parent node directory to `iDir`
  inodes, terminus = make_inodes(idv)

  # parameter update vector
  pup  = Random.SamplerRangeFast(0:lastindex(inodes))
  pups = rand(pup, lastindex(inodes))

  # make scaling function
  scalef = makescalef(obj_ar)

  # burn-in phase
  llc, prc, σλc, σλtn =
    mcmc_burn_gbmpb(treec, treep, λa_prior, σλprior, nburn, tune_int, σλtni, 
      δt, srδt, idv, inodes, terminus, pups, pup, prints, scalef)

  # mcmc
  R, treev = mcmc_gbmpb(treec, treep, llc, prc, σλc, λa_prior, σλprior, 
        niter, nthin, σλtn, δt, srδt, idv, inodes, terminus, pups, pup, prints)

  pardic = Dict(("lambda_root" => 1,
                 "sigma_lambda" => 2))

  write_ssr(R, pardic, out_file)

  return R, treev
end




"""
    mcmc_burn_gbmpb(treec   ::iTgbmpb,
                    treep   ::iTgbmpb,
                    λa_prior::Tuple{Float64,Float64},
                    σλprior::Float64,
                    nburn   ::Int64,
                    tune_int::Int64,
                    σλtni  ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    pups    ::Array{Int64,1},
                    pup     ::Random.SamplerRangeFast{UInt64,Int64},
                    scalef  ::Function)

MCMC chain for constant pure-birth.
"""
function mcmc_burn_gbmpb(treec   ::iTgbmpb,
                         treep   ::iTgbmpb,
                         λa_prior::Tuple{Float64,Float64},
                         σλprior::Float64,
                         nburn   ::Int64,
                         tune_int::Int64,
                         σλtni  ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         idv     ::Array{iDir,1},
                         inodes  ::Array{Int64,1},
                         terminus::Array{BitArray{1}},
                         pups    ::Array{Int64,1},
                         pup     ::Random.SamplerRangeFast{UInt64,Int64},
                         prints  ::Int64,
                         scalef  ::Function)

  # initialize acceptance log
  ltn = 0
  lup = 0.0
  lac = 0.0
  σλtn = σλtni

  # starting parameters
  σλc = 1.0
  llc  = llik_gbm(treec, σλc, δt, srδt)
  prc  = logdexp(σλc, σλprior) + 
         logdnorm_tc(lλ(treec)[1], λa_prior[1], λa_prior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    rand!(pups, pup)

    for pupi in pups
      ## parameter updates
      if iszero(pupi)
        # `λ` diffusion rate updates
        llc, prc, σλc, ltn, lup, lac = 
          update_σλ!(σλc, treec, llc, prc, 
            σλtn, ltn, lup, lac, δt, srδt, σλprior)
      else
        # gbm updates
        ter  = terminus[pupi]
        inod = inodes[pupi]

        dri = dr(idv[inod])
        ldr = lastindex(dri)

        llc, prc = 
          lλupdate!(treec, treep, llc, prc, σλc, δt, srδt, 
            λa_prior, dri, ldr, ter, 0)
      end

      if ltn == tune_int
        σλtn = scalef(σλtn,lac/lup)
        ltn = 0
      end

    end

    next!(pbar)
  end

  return llc, prc, σλc, σλtn
end




"""
    mcmc_gbmpb(treec   ::iTgbmpb,
               treep   ::iTgbmpb,
               llc     ::Float64,
               prc     ::Float64,
               σλc    ::Float64,
               λa_prior::Tuple{Float64,Float64},
               σλprior ::Float64,
               niter   ::Int64,
               nthin   ::Int64,
               σλtn    ::Float64,
               δt      ::Float64,
               srδt    ::Float64,
               idv     ::Array{iDir,1},
               inodes  ::Array{Int64,1},
               terminus::Array{BitArray{1}},
               pups    ::Array{Int64,1},
               pup     ::Random.SamplerRangeFast{UInt64,Int64},
               prints  ::Int64)

MCMC chain for constant pure-birth.
"""
function mcmc_gbmpb(treec   ::iTgbmpb,
                    treep   ::iTgbmpb,
                    llc     ::Float64,
                    prc     ::Float64,
                    σλc    ::Float64,
                    λa_prior::Tuple{Float64,Float64},
                    σλprior ::Float64,
                    niter   ::Int64,
                    nthin   ::Int64,
                    σλtn    ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idv     ::Array{iDir,1},
                    inodes  ::Array{Int64,1},
                    terminus::Array{BitArray{1}},
                    pups    ::Array{Int64,1},
                    pup     ::Random.SamplerRangeFast{UInt64,Int64},
                    prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  R = Array{Float64,2}(undef, nlogs, 5)

  # make tree vector
  treev = iTgbmpb[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    rand!(pups, pup)

    for pupi in pups
      ## parameter updates
      if iszero(pupi)
        # `λ` diffusion rate updates
        llc, prc, σλc = 
          update_σλ!(σλc, treec, llc, prc, σλtn, δt, srδt, σλprior)
      else
        # gbm updates
        ter  = terminus[pupi]
        inod = inodes[pupi]

        dri = dr(idv[inod])
        ldr = lastindex(dri)

        llc, prc = 
          lλupdate!(treec, treep, llc, prc, σλc, δt, srδt, 
            λa_prior, dri, ldr, ter, 0)
      end

    end

    # log parameters
    lthin += 1
    if lthin == nthin
      lit += 1
      @inbounds begin
        R[lit,1] = Float64(lit)
        R[lit,2] = llc
        R[lit,3] = prc
        R[lit,4] = exp(lλ(treec)[1])
        R[lit,5] = σλc
        push!(treev, deepcopy(treec))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev
end




"""
    lλupdate!(treec   ::iTgbmpb,
              treep   ::iTgbmpb,
              llc     ::Float64, 
              prc     ::Float64,
              σλ      ::Float64, 
              δt      ::Float64, 
              srδt    ::Float64, 
              λa_prior::Tuple{Float64,Float64},
              dri     ::BitArray{1},
              ldr     ::Int64,
              ter     ::BitArray{1},
              ix      ::Int64)

Make a λgbm update for a triad.
"""
function lλupdate!(treec   ::iTgbmpb,
                   treep   ::iTgbmpb,
                   llc     ::Float64, 
                   prc     ::Float64,
                   σλ      ::Float64, 
                   δt      ::Float64, 
                   srδt    ::Float64, 
                   λa_prior::Tuple{Float64,Float64},
                   dri     ::BitArray{1},
                   ldr     ::Int64,
                   ter     ::BitArray{1},
                   ix      ::Int64)

  if ix == ldr 
    # if root
    if ldr == 0
      llc, prc = triad_lλupdate_root!(treec::iTgbmpb, treep::iTgbmpb, 
                   llc, prc, σλ, δt, srδt, λa_prior)
    else
      if ter[1]
        if ter[2]
          # if both are terminal
          llc = triad_lλupdate_noded12!(treec::iTgbmpb, treep::iTgbmpb, 
                        llc, σλ, δt, srδt)
        else
          # if d1 is terminal
          llc = triad_lλupdate_noded1!(treec::iTgbmpb, treep::iTgbmpb, 
                       llc, σλ, δt, srδt)
        end
      elseif ter[2]
        # if d2 is terminal
        llc = triad_lλupdate_noded2!(treec::iTgbmpb, treep::iTgbmpb, 
                     llc, σλ, δt, srδt)
      else
        # if no terminal branches involved
        llc = triad_lλupdate_node!(treec::iTgbmpb, treep::iTgbmpb, 
                     llc, σλ, δt, srδt)
      end
    end

  elseif ix < ldr
    ix += 1
    if dri[ix]
      llc, prc = 
        lλupdate!(treec.d1::iTgbmpb, treep.d1::iTgbmpb, 
          llc, prc, σλ, δt, srδt, λa_prior, dri, ldr, ter, ix)
    else
      llc, prc = 
        lλupdate!(treec.d2::iTgbmpb, treep.d2::iTgbmpb, 
          llc, prc, σλ, δt, srδt, λa_prior, dri, ldr, ter, ix)
    end
  end

  return llc, prc
end




"""
    triad_lλupdate_noded12!(treec::iTgbmpb, 
                            treep::iTgbmpb,
                            llc  ::Float64,
                            σλ   ::Float64, 
                            δt   ::Float64)
                            srδt ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
both daughters are terminal.
"""
function triad_lλupdate_noded12!(treec::iTgbmpb, 
                                 treep::iTgbmpb,
                                 llc  ::Float64,
                                 σλ   ::Float64, 
                                 δt   ::Float64,
                                 srδt ::Float64)

  # get reference vectors
  λprv_p = lλ(treep)
  λd1v_p = lλ(treep.d1)
  λd2v_p = lλ(treep.d2)
  λprv_c = lλ(treec)
  λd1v_c = lλ(treec.d1)
  λd2v_c = lλ(treec.d2)
  λpr = λprv_c[1]

  # time vectors
  tprv = ts(treec)
  td1v = ts(treec.d1)
  td2v = ts(treec.d2)

  # make sigma proposal
  σλϕ = randexp()*10.0

  # fill with Brownian motion
  bm!(λprv_p, λpr, tprv, σλϕ, srδt)
  lλp = λprv_p[end]
  bm!(λd1v_p, lλp, td1v, σλϕ, srδt)
  bm!(λd2v_p, lλp, td2v, σλϕ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr = ll_gbm_b(tprv, λprv_p, σλ, δt, srδt) + 
        ll_gbm_b(td1v, λd1v_p, σλ, δt, srδt) + 
        ll_gbm_b(td2v, λd2v_p, σλ, δt, srδt) -
        ll_gbm_b(tprv, λprv_c, σλ, δt, srδt) -
        ll_gbm_b(td1v, λd1v_c, σλ, δt, srδt) -
        ll_gbm_b(td2v, λd2v_c, σλ, δt, srδt) +
        2.0*lλp - 2.0*λd1v_c[1]

  # proposal ratio
  propr = ll_bm(λprv_c, tprv, σλϕ, srδt) +
          ll_bm(λd1v_c, td1v, σλϕ, srδt) +
          ll_bm(λd2v_c, td2v, σλϕ, srδt) +
          ll_bm(λprv_p, tprv, σλϕ, srδt) -
          ll_bm(λd1v_p, td1v, σλϕ, srδt) -
          ll_bm(λd2v_p, td2v, σλϕ, srδt)

  # acceptance ratio
  acr = llr + propr

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_noded1!(treec::iTgbmpb, 
                           treep::iTgbmpb,
                           llc  ::Float64,
                           σλ   ::Float64, 
                           δt   ::Float64)
                           srδt ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 1 is terminal.
"""
function triad_lλupdate_noded1!(treec::iTgbmpb, 
                                treep::iTgbmpb,
                                llc  ::Float64,
                                σλ   ::Float64, 
                                δt   ::Float64,
                                srδt ::Float64)

  # get reference vectors
  λprv_p = lλ(treep)
  λd1v_p = lλ(treep.d1)
  λd2v_p = lλ(treep.d2)
  λprv_c = lλ(treec)
  λd1v_c = lλ(treec.d1)
  λd2v_c = lλ(treec.d2)
  λpr = λprv_c[1]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(treec)
  td1v = ts(treec.d1)
  td2v = ts(treec.d2)

  # make sigma proposal
  σλϕ = randexp()*10.0

  # node proposal
  lλp = duoprop(λpr, λd2, pe(treec), pe(treec.d2), σλϕ)

  bb!(λprv_p, λpr, lλp, tprv, σλϕ, srδt)
  bm!(λd1v_p, lλp, td1v, σλϕ, srδt)
  bb!(λd2v_p, lλp, λd2, td2v, σλϕ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr = ll_gbm_b(tprv, λprv_p, σλ, δt, srδt) + 
        ll_gbm_b(td1v, λd1v_p, σλ, δt, srδt) + 
        ll_gbm_b(td2v, λd2v_p, σλ, δt, srδt) -
        ll_gbm_b(tprv, λprv_c, σλ, δt, srδt) -
        ll_gbm_b(td1v, λd1v_c, σλ, δt, srδt) -
        ll_gbm_b(td2v, λd2v_c, σλ, δt, srδt) +
        2.0*lλp - 2.0*λd1v_c[1]

  # proposal ratio
  propr = ll_bm(λprv_c, tprv, σλϕ, srδt) +
          ll_bm(λd1v_c, td1v, σλϕ, srδt) +
          ll_bm(λd2v_c, td2v, σλϕ, srδt) +
          ll_bm(λprv_p, tprv, σλϕ, srδt) -
          ll_bm(λd1v_p, td1v, σλϕ, srδt) -
          ll_bm(λd2v_p, td2v, σλϕ, srδt)

  # acceptance ratio
  acr = llr + propr

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_noded2!(treec::iTgbmpb, 
                           treep::iTgbmpb,
                           llc  ::Float64,
                           σλ   ::Float64, 
                           δt   ::Float64,
                           srδt ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 2 is terminal.
"""
function triad_lλupdate_noded2!(treec::iTgbmpb, 
                                treep::iTgbmpb,
                                llc  ::Float64,
                                σλ   ::Float64, 
                                δt   ::Float64,
                                srδt ::Float64)

  # get reference vectors
  λprv_p = lλ(treep)
  λd1v_p = lλ(treep.d1)
  λd2v_p = lλ(treep.d2)
  λprv_c = lλ(treec)
  λd1v_c = lλ(treec.d1)
  λd2v_c = lλ(treec.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]

  # time vectors
  tprv = ts(treec)
  td1v = ts(treec.d1)
  td2v = ts(treec.d2)

  # make sigma proposal
  σλϕ = randexp()*10.0

  # node proposal
  lλp = duoprop(λpr, λd1, pe(treec), pe(treec.d1), σλϕ)

  bb!(λprv_p, λpr, lλp, tprv, σλϕ, srδt)
  bb!(λd1v_p, lλp, λd1, td1v, σλϕ, srδt)
  bm!(λd2v_p, lλp, td2v, σλϕ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr = ll_gbm_b(tprv, λprv_p, σλ, δt, srδt) + 
        ll_gbm_b(td1v, λd1v_p, σλ, δt, srδt) + 
        ll_gbm_b(td2v, λd2v_p, σλ, δt, srδt) -
        ll_gbm_b(tprv, λprv_c, σλ, δt, srδt) -
        ll_gbm_b(td1v, λd1v_c, σλ, δt, srδt) -
        ll_gbm_b(td2v, λd2v_c, σλ, δt, srδt) +
        2.0*lλp - 2.0*λd1v_c[1]

  # proposal ratio
  propr = ll_bm(λprv_c, tprv, σλϕ, srδt) +
          ll_bm(λd1v_c, td1v, σλϕ, srδt) +
          ll_bm(λd2v_c, td2v, σλϕ, srδt) +
          ll_bm(λprv_p, tprv, σλϕ, srδt) -
          ll_bm(λd1v_p, td1v, σλϕ, srδt) -
          ll_bm(λd2v_p, td2v, σλϕ, srδt)

  # acceptance ratio
  acr = llr + propr

  if -randexp() < acr 
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_node!(treec     ::iTgbmpb, 
                         treep    ::iTgbmpb,
                         llc     ::Float64,
                         σλ      ::Float64,
                         δt      ::Float64, 
                         srδt    ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
no daughters are terminal.
"""
function triad_lλupdate_node!(treec   ::iTgbmpb, 
                              treep   ::iTgbmpb,
                              llc     ::Float64,
                              σλ      ::Float64,
                              δt      ::Float64, 
                              srδt    ::Float64)

  # get reference vectors
  λprv_p = lλ(treep)
  λd1v_p = lλ(treep.d1)
  λd2v_p = lλ(treep.d2)
  λprv_c = lλ(treec)
  λd1v_c = lλ(treec.d1)
  λd2v_c = lλ(treec.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(treec)
  td1v = ts(treec.d1)
  td2v = ts(treec.d2)

  # make sigma proposal
  σλϕ = randexp()*10.0

  # node proposal
  lλp = trioprop(λpr, λd1, λd2, pe(treep), pe(treep.d1), pe(treep.d1), σλϕ)

  bb!(λprv_p, λpr, lλp, tprv, σλϕ, srδt)
  bb!(λd1v_p, lλp, λd1, td1v, σλϕ, srδt)
  bb!(λd2v_p, lλp, λd2, td2v, σλϕ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr = ll_gbm_b(tprv, λprv_p, σλ, δt, srδt) + 
        ll_gbm_b(td1v, λd1v_p, σλ, δt, srδt) + 
        ll_gbm_b(td2v, λd2v_p, σλ, δt, srδt) -
        ll_gbm_b(tprv, λprv_c, σλ, δt, srδt) -
        ll_gbm_b(td1v, λd1v_c, σλ, δt, srδt) -
        ll_gbm_b(td2v, λd2v_c, σλ, δt, srδt) +
        2.0*lλp - 2.0*λd1v_c[1]

  # proposal ratio
  propr = ll_bm(λprv_c, tprv, σλϕ, srδt) +
          ll_bm(λd1v_c, td1v, σλϕ, srδt) +
          ll_bm(λd2v_c, td2v, σλϕ, srδt) +
          ll_bm(λprv_p, tprv, σλϕ, srδt) -
          ll_bm(λd1v_p, td1v, σλϕ, srδt) -
          ll_bm(λd2v_p, td2v, σλϕ, srδt)

  # acceptance ratio
  acr = llr + propr

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_root!(treec   ::iTgbmpb, 
                         treep   ::iTgbmpb,
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
function triad_lλupdate_root!(treec   ::iTgbmpb, 
                              treep   ::iTgbmpb,
                              llc     ::Float64,
                              prc     ::Float64,
                              σλ      ::Float64, 
                              δt      ::Float64,
                              srδt    ::Float64,
                              λa_prior::Tuple{Float64, Float64})

  # get reference vectors
  λprv_p = lλ(treep)
  λd1v_p = lλ(treep.d1)
  λd2v_p = lλ(treep.d2)
  λprv_c = lλ(treec)
  λd1v_c = lλ(treec.d1)
  λd2v_c = lλ(treec.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(treec)
  td1v = ts(treec.d1)
  td2v = ts(treec.d2)

  # make sigma proposal
  σλϕ = randexp()*10.0

  # proposal given daughters
  lλp = duoprop(λd1, λd2, pe(treec.d1), pe(treec.d2), σλϕ)

  # propose for root
  lλrp = rnorm(lλp, sqrt(pe(treec))*σλϕ)

  # make Brownian bridge proposals
  bb!(λprv_p, lλrp, lλp, tprv, σλϕ, srδt)
  bb!(λd1v_p, lλp,  λd1, td1v, σλϕ, srδt)
  bb!(λd2v_p, lλp,  λd2, td2v, σλϕ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr = ll_gbm_b(tprv, λprv_p, σλ, δt, srδt) + 
        ll_gbm_b(td1v, λd1v_p, σλ, δt, srδt) + 
        ll_gbm_b(td2v, λd2v_p, σλ, δt, srδt) -
        ll_gbm_b(tprv, λprv_c, σλ, δt, srδt) -
        ll_gbm_b(td1v, λd1v_c, σλ, δt, srδt) -
        ll_gbm_b(td2v, λd2v_c, σλ, δt, srδt) +
        2.0*lλp - 2.0*λd1v_c[1]

  # proposal ratio
  propr = ll_bm(λprv_c, tprv, σλϕ, srδt) +
          ll_bm(λd1v_c, td1v, σλϕ, srδt) +
          ll_bm(λd2v_c, td2v, σλϕ, srδt) +
          ll_bm(λprv_p, tprv, σλϕ, srδt) -
          ll_bm(λd1v_p, td1v, σλϕ, srδt) -
          ll_bm(λd2v_p, td2v, σλϕ, srδt)

  # prior ratio
  prr = llrdnorm_x(lλrp, λpr, λa_prior[1], λa_prior[2])

  # acceptance ratio
  acr = llr + propr + prr

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
    update_σλ!(σλc ::Float64,
               llc  ::Float64,
               prc  ::Float64,
               σλtn::Float64,
               δt   ::Float64,
               σλprior::Float64)

MCMC update for σλ.
"""
function update_σλ!(σλc    ::Float64,
                    tree   ::iTgbmpb,
                    llc    ::Float64,
                    prc    ::Float64,
                    σλtn   ::Float64,
                    δt     ::Float64,
                    srδt   ::Float64,
                    σλprior::Float64)

  # parameter proposals
  σλp = mulupt(σλc, σλtn)::Float64

  # one could make a ratio likelihood function
  llp = llik_gbm(tree, σλp, δt, srδt)
  prr = llrdexp_x(σλp, σλc, σλprior)

  if -randexp() < (llp - llc + prr + log(σλp/σλc))
    llc  = llp::Float64
    prc += prr::Float64
    σλc  = σλp::Float64
  end

  return llc, prc, σλc
end




"""
    update_σλ!(σλc    ::Float64,
               tree   ::iTgbmpb,
               llc    ::Float64,
               prc    ::Float64,
               σλtn   ::Float64,
               ltn    ::Int64,
               lup    ::Float64,
               lac    ::Float64,
               δt     ::Float64,
               srδt   ::Float64,
               σλprior::Float64)

MCMC update for σλ with acceptance log.
"""
function update_σλ!(σλc    ::Float64,
                    tree   ::iTgbmpb,
                    llc    ::Float64,
                    prc    ::Float64,
                    σλtn   ::Float64,
                    ltn    ::Int64,
                    lup    ::Float64,
                    lac    ::Float64,
                    δt     ::Float64,
                    srδt   ::Float64,
                    σλprior::Float64)

  # parameter proposals
  σλp = mulupt(σλc, σλtn)::Float64

  # one could make a ratio likelihood function
  llp = llik_gbm(tree, σλp, δt, srδt)
  prr = llrdexp_x(σλp, σλc, σλprior)

  ltn += 1
  lup += 1.0

  if -randexp() < (llp - llc + prr + log(σλp/σλc))
    llc  = llp::Float64
    prc += prr::Float64
    σλc  = σλp::Float64
    lac += 1.0
  end

  return llc, prc, σλc, ltn, lup, lac
end

