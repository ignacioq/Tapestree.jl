#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#





"""
    daughters_lprop12!(treep::iTgbmbd, 
                       treec::iTgbmbd,
                       bbλp::Array{Array{Float64,1},1}, 
                       bbμp::Array{Array{Float64,1},1}, 
                       bbλc::Array{Array{Float64,1},1}, 
                       bbμc::Array{Array{Float64,1},1}, 
                       tsv ::Array{Array{Float64,1},1}, 
                       d1  ::Int64,
                       d2  ::Int64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       δt  ::Float64, 
                       srδt::Float64)

Make a `gbm-bd` proposal for daughters when node is internal and both
daughters are terminal.
"""
function daughters_lprop!(treep::iTgbmbd, 
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
                          ter  ::BitArray{1},
                          σλ   ::Float64,
                          σμ   ::Float64,
                          δt   ::Float64, 
                          srδt ::Float64)

  ## get reference vectors
  # time vectors
  td1v = tsv[d1]
  td2v = tsv[d2]
  lid1 = lastindex(td1v)
  lid2 = lastindex(td2v)

  # speciation vectors
  λd1v_p = bbλp[d1]
  λd2v_p = bbλp[d2]
  λd1v_c = bbλc[d1]
  λd2v_c = bbλc[d2]
  λd1    = λd1v_c[lid1]
  λd2    = λd2v_c[lid2]

  # extinction vectors
  μd1v_p = bbμp[d1]
  μd2v_p = bbμp[d2]
  μd1v_c = bbμc[d1]
  μd2v_c = bbμc[d2]
  μd1    = μd1v_c[lid1]
  μd2    = μd2v_c[lid2]

  # get fixed daughters
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

  # pendant edges
  ped1 = td1v[lid1]
  ped2 = td2v[lid2]

  # simulate fix tree vector
  if ter[1]
    if ter[2]
      # if both are terminal
      bm!(λd1v_p, μd1v_p, λf, μf, td1v, σλ, σμ, srδt)
      bm!(λd2v_p, μd2v_p, λf, μf, td2v, σλ, σμ, srδt)
    else
      # if d1 is terminal
      bm!(λd1v_p, μd1v_p, λf, μf, td1v, σλ, σμ, srδt)
      bb!(λd2v_p, λf, λd2, μd2v_p, μf, μd2, td2v, σλ, σμ, srδt)
    end
  elseif ter[2]
    # if d2 is terminal
    bb!(λd1v_p, λf, λd1, μd1v_p, μf, μd1, td1v, σλ, σμ, srδt)
    bm!(λd2v_p, μd2v_p, λf, μf, td2v, σλ, σμ, srδt)
  else
    # if no terminal branches involved
    bb!(λd1v_p, λf, λd1, μd1v_p, μf, μd1, td1v, σλ, σμ, srδt)
    bb!(λd2v_p, λf, λd2, μd2v_p, μf, μd2, td2v, σλ, σμ, srδt)
  end

  # fill fix and simulate unfix tree
  bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
  bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

  llrbm_d1, llrbd_d1 = llr_gbm_sep_f(treepd1, treecd1, σλ, σμ, δt, srδt)
  llrbm_d2, llrbd_d2 = llr_gbm_sep_f(treepd2, treecd2, σλ, σμ, δt, srδt)

  acr = llrbd_d1 + llrbd_d2 
  llr = llrbm_d1 + llrbm_d2 + acr

  return llr, acr
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
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

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
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

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
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

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
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

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
  treecd1, treecd2 = fixds(treec)
  treepd1, treepd2 = fixds(treep)

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

  llrbm_pr, llrbd_pr = llr_gbm_sep_f(treep,   treec,   σλ, σμ, δt, srδt)
  llrbm_d1, llrbd_d1 = llr_gbm_sep_f(treepd1, treecd1, σλ, σμ, δt, srδt)
  llrbm_d2, llrbd_d2 = llr_gbm_sep_f(treepd2, treecd2, σλ, σμ, δt, srδt)

  acr = llrbd_pr + llrbd_d1 + llrbd_d2 
  llr = llrbm_pr + llrbm_d1 + llrbm_d2 + acr

  return llr, acr
end

