#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#





"""
    daughters_lprop!(treep::iTgbmct, 
                     treec::iTgbmct,
                     λf   ::Float64,
                     bbλp ::Array{Array{Float64,1},1}, 
                     bbλc ::Array{Array{Float64,1},1}, 
                     tsv  ::Array{Array{Float64,1},1}, 
                     pr   ::Int64,
                     d1   ::Int64,
                     d2   ::Int64,
                     ter  ::BitArray{1},
                     ϵ    ::Float64,
                     σλ   ::Float64,
                     icr  ::Bool, 
                     wbc  ::Int64,
                     δt   ::Float64, 
                     srδt ::Float64)

Make a `gbm-bd` proposal for daughters when node is internal and both
daughters are terminal.
"""
function daughters_lprop!(treep::iTgbmct, 
                          treec::iTgbmct,
                          λf   ::Float64,
                          bbλp ::Array{Array{Float64,1},1}, 
                          bbλc ::Array{Array{Float64,1},1}, 
                          tsv  ::Array{Array{Float64,1},1}, 
                          pr   ::Int64,
                          d1   ::Int64,
                          d2   ::Int64,
                          ter  ::BitArray{1},
                          ϵ    ::Float64,
                          σλ   ::Float64,
                          icr  ::Bool, 
                          wbc  ::Int64,
                          δt   ::Float64, 
                          srδt ::Float64)

  ## get reference vectors
  # time vectors
  td1v = tsv[d1]
  td2v = tsv[d2]
  # speciation vectors
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  lid1   = lastindex(λd1v_c)
  lid2   = lastindex(λd2v_c)
  λd1    = λd1v_c[lid1]
  λd2    = λd2v_c[lid2]

  # get fixed daughters
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

  # pendant edges
  ped1 = td1v[1]
  ped2 = td2v[2]

  bb!(λd1v_p, λf, λd1, td1v[2], σλ, δt, srδt)
  bb!(λd2v_p, λf, λd2, td2v[2], σλ, δt, srδt)

  normprop = duoldnorm(λf,        λd1, λd2, ped1, ped2, σλ) -
             duoldnorm(λd1v_c[1], λd1, λd2, ped1, ped2, σλ)

  llrcond = 0.0
  # only if crown conditioning
  if icr && iszero(wbc)
    llrcond += cond_alone_events_crown(treep, ϵ) -
               cond_alone_events_crown(treec, ϵ)
  end

  # fill fix and simulate unfix tree
  bm!(treepd1, λd1v_p, 1, lid1, σλ, srδt)
  bm!(treepd2, λd2v_p, 1, lid2, σλ, srδt)

  llrbm_d1, llrbd_d1 = llr_gbm_sep_f(treepd1, treecd1, ϵ, σλ, δt, srδt)
  llrbm_d2, llrbd_d2 = llr_gbm_sep_f(treepd2, treecd2, ϵ, σλ, δt, srδt)

  acr  = llrbd_d1 + llrbd_d2 + llrcond
  llr  = llrbm_d1 + llrbm_d2 + acr
  acr += normprop

  return llr, acr
end




"""
    triad_lvupdate_trio!(treep::iTgbmct, 
                         treec::iTgbmct,
                         bbλp ::Array{Array{Float64,1},1}, 
                         bbλc ::Array{Array{Float64,1},1}, 
                         tsv  ::Array{Array{Float64,1},1}, 
                         llc  ::Float64,
                         pr   ::Int64,
                         d1   ::Int64,
                         d2   ::Int64,
                         σλ   ::Float64,
                         δt   ::Float64, 
                         srδt ::Float64,
                         ter  ::BitArray{1},
                         icr  ::Bool,
                         wbc  ::Int64)

Make a trio of Brownian motion MCMC updates when node is internal and 
no daughters are terminal.
"""
function triad_lvupdate_trio!(treep::iTgbmct, 
                              treec::iTgbmct,
                              bbλp ::Array{Array{Float64,1},1}, 
                              bbλc ::Array{Array{Float64,1},1}, 
                              tsv  ::Array{Array{Float64,1},1}, 
                              llc  ::Float64,
                              pr   ::Int64,
                              d1   ::Int64,
                              d2   ::Int64,
                              ϵ    ::Float64,
                              σλ   ::Float64,
                              δt   ::Float64, 
                              srδt ::Float64,
                              ter  ::BitArray{1},
                              icr  ::Bool,
                              wbc  ::Int64)

  ## get reference vectors
  # time vectors
  tprv = tsv[pr]
  td1v = tsv[d1]
  td2v = tsv[d2]

  # speciation vectors
  λprv_p = bbλp[pr]
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λprv_c = bbλc[pr]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  lipr = lastindex(λprv_c)
  lid1 = lastindex(λd1v_c)
  lid2 = lastindex(λd2v_c)
  λpr    = λprv_c[1]
  λd1    = λd1v_c[lid1]
  λd2    = λd2v_c[lid2]

  # get fixed daughters
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

  # pendant edges
  pepr = tprv[1]
  ped1 = td1v[1]
  ped2 = td2v[1]

  if ter[1]
    if ter[2]
      # if both are terminal
      bm!(λprv_p, λpr, tprv[2], σλ, srδt)
      lλp = λprv_p[lipr]
      bm!(λd1v_p, lλp, td1v[2], σλ, srδt)
      bm!(λd2v_p, lλp, td2v[2], σλ, srδt)

      # fill fix and simulate unfix tree
      bm!(treep,   λprv_p, 1, lipr, σλ, srδt)
      bm!(treepd1, λd1v_p, 1, lid1, σλ, srδt)
      bm!(treepd2, λd2v_p, 1, lid2, σλ, srδt)

    else
      # if d1 is terminal
      lλp = duoprop(λpr, λd2, pepr, ped2, σλ)

      # simulate fix tree vector
      bb!(λprv_p, λpr, lλp, tprv[2], σλ, δt, srδt)
      bm!(λd1v_p, lλp, td1v[2], σλ, srδt)
      bb!(λd2v_p, lλp, λd2, td2v[2], σλ, δt, srδt)

      # fill fix and simulate unfix tree
      bm!(treep,   λprv_p, 1, lipr, σλ, srδt)
      bm!(treepd1, λd1v_p, 1, lid1, σλ, srδt)
      bm!(treepd2, λd2v_p, 1, lid2, σλ, srδt)

    end
  elseif ter[2]
    # if d2 is terminal
    # node proposal
    lλp = duoprop(λpr, λd1, pepr, ped1, σλ)

    # simulate fix tree vector
    bb!(λprv_p, λpr, lλp, tprv[2], σλ, δt, srδt)
    bb!(λd1v_p, lλp, λd1, td1v[2], σλ, δt, srδt)
    bm!(λd2v_p, lλp, td2v[2], σλ, srδt)

    # fill fix and simulate unfix tree
    bm!(treep,   λprv_p, 1, lipr, σλ, srδt)
    bm!(treepd1, λd1v_p, 1, lid1, σλ, srδt)
    bm!(treepd2, λd2v_p, 1, lid2, σλ, srδt)

  else
    # if no terminal branches involved
    # node proposal
    lλp  = trioprop(λpr, λd1, λd2, pepr, ped1, ped2, σλ)

    # simulate fix tree vector
    bb!(λprv_p, λpr, lλp, tprv[2], σλ, δt, srδt)
    bb!(λd1v_p, lλp, λd1, td1v[2], σλ, δt, srδt)
    bb!(λd2v_p, lλp, λd2, td2v[2], σλ, δt, srδt)

    # fill fix and simulate unfix tree
    bm!(treep,   λprv_p, 1, lipr, σλ, srδt)
    bm!(treepd1, λd1v_p, 1, lid1, σλ, srδt)
    bm!(treepd2, λd2v_p, 1, lid2, σλ, srδt)
  end

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       ϵ, σλ, δt, srδt, icr, wbc)

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
    gbm_copy_f!(treec, treep)
    gbm_copy_f!(treecd1, treepd1)
    gbm_copy_f!(treecd2, treepd2)
  end

  return llc
end




"""
    triad_lupdate_root!(treep ::iTgbmct, 
                        treec ::iTgbmct,
                        bbλp  ::Array{Array{Float64,1},1}, 
                        bbλc  ::Array{Array{Float64,1},1}, 
                        tsv   ::Array{Array{Float64,1},1}, 
                        llc   ::Float64,
                        pr    ::Int64,
                        d1    ::Int64,
                        d2    ::Int64,
                        ϵ     ::Float64,
                        σλ    ::Float64,
                        δt    ::Float64, 
                        srδt  ::Float64,
                        lλmxpr::Float64,
                        icr   ::Bool)

Make a trio of Brownian motion MCMC updates when the root is involved.
"""
function triad_lupdate_root!(treep ::iTgbmct, 
                             treec ::iTgbmct,
                             bbλp  ::Array{Array{Float64,1},1}, 
                             bbλc  ::Array{Array{Float64,1},1}, 
                             tsv   ::Array{Array{Float64,1},1}, 
                             llc   ::Float64,
                             pr    ::Int64,
                             d1    ::Int64,
                             d2    ::Int64,
                             ϵ     ::Float64,
                             σλ    ::Float64,
                             δt    ::Float64, 
                             srδt  ::Float64,
                             lλmxpr::Float64,
                             icr   ::Bool)

  ## get reference vectors
  # time vectors
  tprv = tsv[pr]
  td1v = tsv[d1]
  td2v = tsv[d2]

  # speciation vectors
  λprv_p = bbλp[pr]
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λprv_c = bbλc[pr]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  lipr = lastindex(λprv_c)
  lid1 = lastindex(λd1v_c)
  lid2 = lastindex(λd2v_c)
  λpr    = λprv_c[1]
  λd1    = λd1v_c[lid1]
  λd2    = λd2v_c[lid2]

  # get fixed daughters
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

  # pendant edges
  pepr = tprv[1]
  ped1 = td1v[1]
  ped2 = td2v[1]

  # proposal given daughters
  lλp = duoprop(λd1, λd2, ped1, ped2, σλ)

  # propose for root
  srpepr = sqrt(pepr)
  lλrp = rnorm(lλp, srpepr*σλ)

  # simulate fix tree vector
  bb!(λprv_p, lλrp, lλp, tprv[2], σλ, δt, srδt)
  bb!(λd1v_p, lλp,  λd1, td1v[2], σλ, δt, srδt)
  bb!(λd2v_p, lλp,  λd2, td2v[2], σλ, δt, srδt)

  # fill fix and simulate unfix tree
  bm!(treep,   λprv_p, 1, lipr, σλ, srδt)
  bm!(treepd1, λd1v_p, 1, lid1, σλ, srδt)
  bm!(treepd2, λd2v_p, 1, lid2, σλ, srδt)

  ## make acceptance ratio 

  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       ϵ, σλ, δt, srδt, icr, 0)

  # prior ratio
  if lλrp > lλmxpr
    acr += -Inf
  end

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
    gbm_copy_f!(treec, treep)
    gbm_copy_f!(treecd1, treepd1)
    gbm_copy_f!(treecd2, treepd2)
  end

  return llc
end




"""
    llr_propr(treep  ::iTgbmct,
              treepd1::iTgbmct,
              treepd2::iTgbmct,
              treec  ::iTgbmct,
              treecd1::iTgbmct,
              treecd2::iTgbmct,
              ϵ      ::Float64,
              σλ     ::Float64,
              δt     ::Float64,
              srδt   ::Float64,
              icr    ::Bool,
              wbc    ::Int64)

Return the likelihood and proposal ratio for birth-death gbm.
"""
function llr_propr(treep  ::iTgbmct,
                   treepd1::iTgbmct,
                   treepd2::iTgbmct,
                   treec  ::iTgbmct,
                   treecd1::iTgbmct,
                   treecd2::iTgbmct,
                   ϵ      ::Float64,
                   σλ     ::Float64,
                   δt     ::Float64,
                   srδt   ::Float64,
                   icr    ::Bool,
                   wbc    ::Int64)

  llrbm_pr, llrbd_pr = llr_gbm_sep_f(treep,   treec,   ϵ, σλ, δt, srδt)
  llrbm_d1, llrbd_d1 = llr_gbm_sep_f(treepd1, treecd1, ϵ, σλ, δt, srδt)
  llrbm_d2, llrbd_d2 = llr_gbm_sep_f(treepd2, treecd2, ϵ, σλ, δt, srδt)

  llrcond = 0.0

  if iszero(wbc)
    if icr 
      llrcond += cond_alone_events_crown(treep, ϵ) -
                 cond_alone_events_crown(treec, ϵ)
    else
      llrcond += cond_alone_events_stem(treep, ϵ) -
                 cond_alone_events_stem(treec, ϵ)
    end
  elseif isone(wbc) && icr
    llrcond += cond_alone_events_stem(treep, ϵ) -
               cond_alone_events_stem(treec, ϵ)
  end

  acr = llrbd_pr + llrbd_d1 + llrbd_d2 + llrcond
  llr = acr + llrbm_pr + llrbm_d1 + llrbm_d2

  return llr, acr
end

