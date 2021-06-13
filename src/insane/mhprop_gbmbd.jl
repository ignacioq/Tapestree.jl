#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#





"""
    daughters_lprop!(treep::iTgbmbd, 
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
                     ter  ::BitArray{1},
                     σλ   ::Float64,
                     σμ   ::Float64,
                     icr  ::Bool, 
                     wbc  ::Int64,
                     δt   ::Float64, 
                     srδt ::Float64)

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
                          pr   ::Int64,
                          d1   ::Int64,
                          d2   ::Int64,
                          ter  ::BitArray{1},
                          σλ   ::Float64,
                          σμ   ::Float64,
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
  lid1 = lastindex(λd1v_c)
  lid2 = lastindex(λd2v_c)
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
  ped1 = td1v[1]
  ped2 = td2v[2]

  # acceptance rate
  normprop = 0.0

  # simulate fix tree vector
  if ter[1]
    if ter[2]
      # if both are terminal
      bm!(λd1v_p, μd1v_p, λf, μf, td1v[2], σλ, σμ, srδt)
      bm!(λd2v_p, μd2v_p, λf, μf, td2v[2], σλ, σμ, srδt)
    else
      # if d1 is terminal
      bm!(λd1v_p, μd1v_p, λf, μf, td1v[2], σλ, σμ, srδt)
      bb!(λd2v_p, λf, λd2, μd2v_p, μf, μd2, td2v[2], σλ, σμ, δt, srδt)

      # acceptance ratio
      λpr1_c = bbλc[pr][1]
      μpr1_c = bbμc[pr][1]
      pepr   = tsv[pr][1]

      # normprop += 
      #   trioldnorm(λf, λpr1_c, λd1, λd2, pepr, ped1, ped2, σλ)        -
      #   trioldnorm(λd1v_c[1], λpr1_c, λd1, λd2, pepr, ped1, ped2, σλ) +
      #   trioldnorm(μf, μpr1_c, μd1, μd2, pepr, ped1, ped2, σμ)        -
      #   trioldnorm(μd1v_c[1], μpr1_c, μd1, μd2, pepr, ped1, ped2, σμ)

      # normprop += 
      #       duoldnorm(λf, λpr1_c, λd2, pepr, ped2, σλ)        -
      #       duoldnorm(λd2v_c[1], λpr1_c, λd2, pepr, ped2, σλ) +
      #       duoldnorm(μf, μpr1_c, μd2, pepr, ped2, σμ)        -
      #       duoldnorm(μd2v_c[1], μpr1_c, μd2, pepr, ped2, σμ)

      # normprop += 
      #        ldnorm_bm(λf, λd2, sqrt(ped2)*σλ)        - 
      #        ldnorm_bm(λd2v_c[1], λd2, sqrt(ped2)*σλ) + 
      #        ldnorm_bm(μf, μd2, sqrt(ped2)*σμ)        - 
      #        ldnorm_bm(μd2v_c[1], μd2, sqrt(ped2)*σμ)
  
      normprop += 
            duoldnorm(λf, λpr1_c, λd2, pepr, ped2, σλ) + 
            duoldnorm(μf, μpr1_c, μd2, pepr, ped2, σμ)

    end
  elseif ter[2]
    # if d2 is terminal
    bb!(λd1v_p, λf, λd1, μd1v_p, μf, μd1, td1v[2], σλ, σμ, δt, srδt)
    bm!(λd2v_p, μd2v_p, λf, μf, td2v[2], σλ, σμ, srδt)

    # acceptance ratio
    λpr1_c = bbλc[pr][1]
    μpr1_c = bbμc[pr][1]
    pepr   = tsv[pr][1]

    # normprop += 
    #   trioldnorm(λf, λpr1_c, λd1, λd2, pepr, ped1, ped2, σλ)        -
    #   trioldnorm(λd1v_c[1], λpr1_c, λd1, λd2, pepr, ped1, ped2, σλ) +
    #   trioldnorm(μf, μpr1_c, μd1, μd2, pepr, ped1, ped2, σμ)        -
    #   trioldnorm(μd1v_c[1], μpr1_c, μd1, μd2, pepr, ped1, ped2, σμ)

    # normprop += 
    #       duoldnorm(λf, λpr1_c, λd1, pepr, ped1, σλ)        -
    #       duoldnorm(λd1v_c[1], λpr1_c, λd1, pepr, ped1, σλ) +
    #       duoldnorm(μf, μpr1_c, μd1, pepr, ped1, σμ)        -
    #       duoldnorm(μd1v_c[1], μpr1_c, μd1, pepr, ped1, σμ)

    # normprop += 
    #    ldnorm_bm(λf, λd1, sqrt(ped1)*σλ)        - 
    #    ldnorm_bm(λd1v_c[1], λd1, sqrt(ped1)*σλ) + 
    #    ldnorm_bm(μf, μd1, sqrt(ped1)*σμ)        - 
    #    ldnorm_bm(μd1v_c[1], μd1, sqrt(ped1)*σμ)

    normprop += 
          duoldnorm(λf, λpr1_c, λd1, pepr, ped1, σλ) +
          duoldnorm(μf, μpr1_c, μd1, pepr, ped1, σμ)

  else
    # if no terminal branches involved
    bb!(λd1v_p, λf, λd1, μd1v_p, μf, μd1, td1v[2], σλ, σμ, δt, srδt)
    bb!(λd2v_p, λf, λd2, μd2v_p, μf, μd2, td2v[2], σλ, σμ, δt, srδt)

    # acceptance ratio
    λpr1_c = bbλc[pr][1]
    μpr1_c = bbμc[pr][1]
    pepr   = tsv[pr][1]

    # normprop += 
    #       trioldnorm(λf, λpr1_c, λd1, λd2, pepr, ped1, ped2, σλ)        -
    #       trioldnorm(λd1v_c[1], λpr1_c, λd1, λd2, pepr, ped1, ped2, σλ) +
    #       trioldnorm(μf, μpr1_c, μd1, μd2, pepr, ped1, ped2, σμ)        -
    #       trioldnorm(μd1v_c[1], μpr1_c, μd1, μd2, pepr, ped1, ped2, σμ)

    # normprop += 
    #        duoldnorm(λf, λd1, λd2, ped1, ped2, σλ)        -
    #        duoldnorm(λd1v_c[1], λd1, λd2, ped1, ped2, σλ) +
    #        duoldnorm(μf, μd1, μd2, ped1, ped2, σμ)        -
    #        duoldnorm(μd1v_c[1], μd1, μd2, ped1, ped2, σμ)

    normprop += 
          trioldnorm(λf, λpr1_c, λd1, λd2, pepr, ped1, ped2, σλ) +
          trioldnorm(μf, μpr1_c, μd1, μd2, pepr, ped1, ped2, σμ)

  end

  llrcond = 0.0
  # only if crown conditioning
  if icr && iszero(wbc)
    llrcond += cond_alone_events_crown(treep) -
               cond_alone_events_crown(treec)
  end

  # fill fix and simulate unfix tree
  bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
  bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

  llrbm_d1, llrbd_d1 = llr_gbm_sep_f(treepd1, treecd1, σλ, σμ, δt, srδt)
  llrbm_d2, llrbd_d2 = llr_gbm_sep_f(treepd2, treecd2, σλ, σμ, δt, srδt)

  acr  = llrbd_d1 + llrbd_d2 + llrcond
  llr  = llrbm_d1 + llrbm_d2 + acr
  acr += normprop

  return llr, acr
end







"""
    triad_lvupdate_trio!(treep::iTgbmbd, 
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
                         srδt ::Float64,
                         ter  ::BitArray{1},
                         wbc  ::Bool,
                         svf  ::Function)

Make a trio of Brownian motion MCMC updates when node is internal and 
no daughters are terminal.
"""
function triad_lvupdate_trio!(treep::iTgbmbd, 
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
  pepr = tprv[1]
  ped1 = td1v[1]
  ped2 = td2v[1]

  if ter[1]
    if ter[2]
      # if both are terminal
      bm!(λprv_p, μprv_p, λpr, μpr, tprv[2], σλ, σμ, srδt)
      lλp = λprv_p[lipr]
      lμp = μprv_p[lipr]
      bm!(λd1v_p, μd1v_p, lλp, lμp, td1v[2], σλ, σμ, srδt)
      bm!(λd2v_p, μd2v_p, lλp, lμp, td2v[2], σλ, σμ, srδt)

      # fill fix and simulate unfix tree
      bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
      bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
      bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

    else
      # if d1 is terminal
      lλp = duoprop(λpr, λd2, pepr, ped2, σλ)
      lμp = duoprop(μpr, μd2, pepr, ped2, σμ)

      # simulate fix tree vector
      bb!(λprv_p, λpr, lλp, μprv_p, μpr, lμp, tprv[2], σλ, σμ, δt, srδt)
      bm!(λd1v_p, μd1v_p, lλp, lμp, td1v[2], σλ, σμ, srδt)
      bb!(λd2v_p, lλp, λd2, μd2v_p, lμp, μd2, td2v[2], σλ, σμ, δt, srδt)

      # fill fix and simulate unfix tree
      bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
      bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
      bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

    end
  elseif ter[2]
    # if d2 is terminal
    # node proposal
    lλp = duoprop(λpr, λd1, pepr, ped1, σλ)
    lμp = duoprop(μpr, μd1, pepr, ped1, σμ)

    # simulate fix tree vector
    bb!(λprv_p, λpr, lλp, μprv_p, μpr, lμp, tprv[2], σλ, σμ, δt, srδt)
    bb!(λd1v_p, lλp, λd1, μd1v_p, lμp, μd1, td1v[2], σλ, σμ, δt, srδt)
    bm!(λd2v_p, μd2v_p, lλp, lμp, td2v[2], σλ, σμ, srδt)

    # fill fix and simulate unfix tree
    bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
    bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
    bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

  else
    # if no terminal branches involved
    # node proposal
    lλp  = trioprop(λpr, λd1, λd2, pepr, ped1, ped2, σλ)
    lμp  = trioprop(μpr, μd1, μd2, pepr, ped1, ped2, σμ)

    # simulate fix tree vector
    bb!(λprv_p, λpr, lλp, μprv_p, μpr, lμp, tprv[2], σλ, σμ, δt, srδt)
    bb!(λd1v_p, lλp, λd1, μd1v_p, lμp, μd1, td1v[2], σλ, σμ, δt, srδt)
    bb!(λd2v_p, lλp, λd2, μd2v_p, lμp, μd2, td2v[2], σλ, σμ, δt, srδt)

    # fill fix and simulate unfix tree
    bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
    bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
    bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)
  end

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       σλ, σμ, δt, srδt, icr, wbc)

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
function triad_lupdate_root!(treep ::iTgbmbd, 
                             treec ::iTgbmbd,
                             bbλp  ::Array{Array{Float64,1},1}, 
                             bbμp  ::Array{Array{Float64,1},1}, 
                             bbλc  ::Array{Array{Float64,1},1}, 
                             bbμc  ::Array{Array{Float64,1},1}, 
                             tsv   ::Array{Array{Float64,1},1}, 
                             llc   ::Float64,
                             pr    ::Int64,
                             d1    ::Int64,
                             d2    ::Int64,
                             σλ    ::Float64,
                             σμ    ::Float64,
                             δt    ::Float64, 
                             srδt  ::Float64,
                             lλmxpr::Float64,
                             lμmxpr::Float64,
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
  pepr = tprv[1]
  ped1 = td1v[1]
  ped2 = td2v[1]

  # proposal given daughters
  lλp = duoprop(λd1, λd2, ped1, ped2, σλ)
  lμp = duoprop(μd1, μd2, ped1, ped2, σμ)

  # propose for root
  srpepr = sqrt(pepr)
  lλrp = rnorm(lλp, srpepr*σλ)
  lμrp = rnorm(lμp, srpepr*σμ)

  # simulate fix tree vector
  bb!(λprv_p, lλrp, lλp, μprv_p, lμrp, lμp, tprv[2], σλ, σμ, δt, srδt)
  bb!(λd1v_p, lλp,  λd1, μd1v_p, lμp,  μd1, td1v[2], σλ, σμ, δt, srδt)
  bb!(λd2v_p, lλp,  λd2, μd2v_p, lμp,  μd2, td2v[2], σλ, σμ, δt, srδt)

  # fill fix and simulate unfix tree
  bm!(treep,   λprv_p, μprv_p, 1, lipr, σλ, σμ, srδt)
  bm!(treepd1, λd1v_p, μd1v_p, 1, lid1, σλ, σμ, srδt)
  bm!(treepd2, λd2v_p, μd2v_p, 1, lid2, σλ, σμ, srδt)

  ## make acceptance ratio 

  # estimate likelihoods
  llr, acr = llr_propr(treep, treepd1, treepd2, 
                       treec, treecd1, treecd2, 
                       σλ, σμ, δt, srδt, icr, 0)

  # prior ratio
  if lλrp > lλmxpr
    acr += -Inf
  end
  if lμrp > lμmxpr
    acr += -Inf
  end

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
                   srδt   ::Float64,
                   icr    ::Bool,
                   wbc    ::Int64)

  llrbm_pr, llrbd_pr = llr_gbm_sep_f(treep,   treec,   σλ, σμ, δt, srδt)
  llrbm_d1, llrbd_d1 = llr_gbm_sep_f(treepd1, treecd1, σλ, σμ, δt, srδt)
  llrbm_d2, llrbd_d2 = llr_gbm_sep_f(treepd2, treecd2, σλ, σμ, δt, srδt)

  llrcond = 0.0

  if iszero(wbc)
    if icr 
      llrcond += cond_alone_events_crown(treep) -
                 cond_alone_events_crown(treec)
    else
      llrcond += cond_alone_events_stem(treep) -
                 cond_alone_events_stem(treec)
    end
  elseif isone(wbc) && icr
    llrcond += cond_alone_events_stem(treep) -
               cond_alone_events_stem(treec)
  end

  acr = llrbd_pr + llrbd_d1 + llrbd_d2 + llrcond
  llr = llrbm_pr + llrbm_d1 + llrbm_d2 + acr

  return llr, acr
end

