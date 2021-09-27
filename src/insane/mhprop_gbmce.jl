#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#



# """
# do this until both fixed
# """



# triad_lvupdate!(treep, treec, αc, σλc, δt, srδt)

# function xx()


#   if isdefined(treec, :d1)
#     triad_lvupdate!(treep, treec, αc, σλc, δt, srδt)

#   else
#     return nothing
#   end

#     """
#     here
#     """


#     ifx1 = 
#     if isfix(treec.d1) && isfix(treec.d2)
#       triad_lvupdate!(treep, treec, αc, σλc, δt, srδt)

#     elseif ifx1
#       triad_lvupdate!(treep.d1, treec.d2, αc, σλc, δt, srδt)
#     else
#       triad_lvupdate!(treep, treec, αc, σλc, δt, srδt)
#     end
#   end
# end





# """
#     triad_lvupdate!(treep::iTgbmce,
#                     treec::iTgbmce,
#                     α    ::Float64,
#                     σλ   ::Float64,
#                     δt   ::Float64,
#                     srδt ::Float64)

# Make a `gbmce` trio proposal.
# """
# function triad_lvupdate!(treep::iTgbmce,
#                          treec::iTgbmce,
#                          α    ::Float64,
#                          σλ   ::Float64,
#                          δt   ::Float64,
#                          srδt ::Float64)

#   @inbounds begin

#     it2 = istip(treec.d2)

#     if istip(treec.d1)
#       # if both daughters are terminal
#       if it2

#         lλpr = lλ(treep)

#         bm!(lλpr, lλ(treec)[1], α, σλ, δt, fdt(treec), srδt)
#         lλp = lλpr[end]
#         bm!(lλ(treep.d1), lλp, α, σλ, δt, fdt(treec.d1), srδt)
#         bm!(lλ(treep.d2), lλp, α, σλ, δt, fdt(treec.d2), srδt)

#       # if d1 is terminal
#       else

#         λpr  = lλ(treec)[1]
#         λd2  = lλ(treec.d2)[end]
#         epr  = e(treec)
#         ed2  = e(treec.d2)

#         # node proposal
#         lλp = duoprop(λpr + α*epr, λd2 - α*ed2, epr, ed2, σλ)

#         # simulate fix tree vector
#         bb!(lλ(treep),    λpr, lλp, σλ, δt, fdt(treec),    srδt)
#         bm!(lλ(treep.d1), lλp,   α, σλ, δt, fdt(treec.d1), srδt)
#         bb!(lλ(treep.d2), lλp, λd2, σλ, δt, fdt(treec.d2), srδt)

#       end
    
#     # if d2 is terminal
#     elseif it2

#       λpr  = lλ(treec)[1]
#       λd1  = lλ(treec.d1)[end]
#       epr  = e(treec)
#       ed1  = e(treec.d1)

#       # node proposal
#       lλp = duoprop(λpr + α*epr, λd1 - α*ed1, epr, ed1, σλ)

#       # simulate fix tree vector
#       bb!(lλ(treep),    λpr, lλp, σλ, δt, fdt(treec),    srδt)
#       bb!(lλ(treep.d1), lλp, λd1, σλ, δt, fdt(treec.d1), srδt)
#       bm!(lλ(treep.d2), lλp,   α, σλ, δt, fdt(treec.d2), srδt)

#     # if no terminal branches involved
#     else

#       λpr  = lλ(treec)[1]
#       λd1  = lλ(treec.d1)[end]
#       λd2  = lλ(treec.d2)[end]
#       epr  = e(treec)
#       ed1  = e(treec.d1)
#       ed2  = e(treec.d2)

#       # node proposal
#       lλp  = trioprop(λpr + α*epr, λd1 - α*ed1, λd2 - α*ed2, 
#                epr, ed1, ed2, σλ)

#       # simulate fix tree vector
#       bb!(lλ(treep),    λpr, lλp, σλ, δt, fdt(treec),    srδt)
#       bb!(lλ(treep.d1), lλp, λd1, σλ, δt, fdt(treec.d1), srδt)
#       bb!(lλ(treep.d2), lλp, λd2, σλ, δt, fdt(treec.d2), srδt)
#     end

#   end

#   """
#   here: return likelihood
#   """



#   return nothing
# end





"""
    daughters_lprop!(treep::iTgbmce, 
                     treec::iTgbmce,
                     λf   ::Float64,
                     α    ::Float64,
                     σλ   ::Float64,
                     μ    ::Float64,
                     δt   ::Float64, 
                     srδt ::Float64)

Make a `gbmce` proposal for daughters of forwards simulated branch.
"""
function daughters_lprop!(treep::iTgbmce, 
                          treec::iTgbmce,
                          λf   ::Float64,
                          α    ::Float64,
                          σλ   ::Float64,
                          μ    ::Float64,
                          δt   ::Float64, 
                          srδt ::Float64)
  @inbounds begin

    # get fixed daughters
    treecd1, treecd2 = fixds(treec)
    treepd1, treepd2 = fixds(treep)

    # edges
    ed1 = e(treecd1)
    ed2 = e(treecd2)

    lλ1 = lλ(treecd1)
    λn  = lλ1[1]
    λd1 = lλ1[end]
    λd2 = lλ(treecd2)[end]

    bb!(lλ(treepd1), λf, λd1, σλ, δt, fdt(treecd1), srδt)
    bb!(lλ(treepd2), λf, λd2, σλ, δt, fdt(treecd2), srδt)

    # acceptance rate
    normprop = 
      duoldnorm(λf, λd1 - α*ed1, λd2 - α*ed2, ed1, ed2, σλ) -
      duoldnorm(λn, λd1 - α*ed1, λd2 - α*ed2, ed1, ed2, σλ)

    llrbm1, llrbd1 = 
      llr_gbm_b_sep(lλ(treepd1), lλ(treecd1), α, σλ, δt, 
        fdt(treecd1), srδt, !istip(treecd1))
    llrbm2, llrbd2 = 
      llr_gbm_b_sep(lλ(treepd2), lλ(treecd2), α, σλ, δt, 
        fdt(treecd2), srδt, !istip(treecd2))

    acr  = llrbd1 + llrbd2
    llr  = llrbm1 + llrbm2 + acr
    acr += normprop
  end

  return llr, acr
end





"""
    root_update!(treep ::iTgbmce, 
                 treec ::iTgbmce,
                 α     ::Float64,
                 σλ    ::Float64,
                 μ     ::Float64,
                 δt    ::Float64, 
                 srδt  ::Float64,
                 lλmxpr::Float64,
                 icr   ::Bool)

Make a trio of Brownian motion MCMC updates when the root is involved.
"""
function root_update!(treep ::iTgbmce, 
                      treec ::iTgbmce,
                      α     ::Float64,
                      σλ    ::Float64,
                      μ     ::Float64,
                      llc   ::Float64,
                      δt    ::Float64, 
                      srδt  ::Float64,
                      lλmxpr::Float64,
                      icr   ::Bool)

  @inbounds begin

    λpp = lλ(treep)
    λ1p = lλ(treep.d1)
    λ2p = lλ(treep.d2)
    λpc = lλ(treec)
    λ1c = lλ(treec.d1)
    λ2c = lλ(treec.d2)

    epr  = e(treec)
    ed1  = e(treec.d1)
    ed2  = e(treec.d2)
    fdtp = fdt(treec)
    fdt1 = fdt(treec.d1)
    fdt2 = fdt(treec.d2)

    λd1 = λ1c[end]
    λd2 = λ2c[end]

    # proposal given daughters
    lλp = duoprop(λd1 - α*ed1, λd2 - α*ed2, ed1, ed2, σλ)

    # propose for root
    lλrp = rnorm(lλp - α*epr, sqrt(epr)*σλ)

    # simulate fix tree vector
    bb!(λpp, lλrp, lλp, σλ, δt, fdtp, srδt)
    bb!(λ1p, lλp,  λd1, σλ, δt, fdt1, srδt)
    bb!(λ2p, lλp,  λd2, σλ, δt, fdt2, srδt)

    llrbm, llrbd = 
      llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtp, srδt, true)
    llrbm1, llrbd1 = 
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, !istip(treec.d1))
    llrbm2, llrbd2 = 
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, !istip(treec.d2))

    if icr 
      llrcond = cond_surv_crown(treep, μ) -
                cond_surv_crown(treec, μ)
    else
      llrcond = cond_surv_stem(treep, μ) -
                cond_surv_stem(treec, μ)
    end

    acr = llrbd + llrbd1 + llrbd2 + llrcond
    llr = acr + llrbm + llrbm1 + llrbm2

    # prior ratio
    if lλrp > lλmxpr
      acr += -Inf
    end

    if -randexp() < acr
      llc += llr
      copyto!(λpc, λpp)
      copyto!(λ1c, λ1p)
      copyto!(λ2c, λ2p)
    end 
  end
  
  return llc
end









































"""
    daughters_lprop!(treep::iTgbmce, 
                     treec::iTgbmce,
                     λf   ::Float64,
                     bbλp ::Array{Array{Float64,1},1}, 
                     bbλc ::Array{Array{Float64,1},1}, 
                     tsv  ::Array{Array{Float64,1},1}, 
                     pr   ::Int64,
                     d1   ::Int64,
                     d2   ::Int64,
                     ter  ::BitArray{1},
                     α    ::Float64,
                     σλ   ::Float64,
                     μ    ::Float64,
                     icr  ::Bool, 
                     wbc  ::Int64,
                     δt   ::Float64, 
                     srδt ::Float64)

Make a `gbmce` proposal for daughters of forwards simulated branch.
"""
function daughters_lprop!(treep::iTgbmce, 
                          treec::iTgbmce,
                          λf   ::Float64,
                          bbλp ::Array{Array{Float64,1},1}, 
                          bbλc ::Array{Array{Float64,1},1}, 
                          tsv  ::Array{Array{Float64,1},1}, 
                          pr   ::Int64,
                          d1   ::Int64,
                          d2   ::Int64,
                          ter  ::BitArray{1},
                          α    ::Float64,
                          σλ   ::Float64,
                          μ    ::Float64,
                          icr  ::Bool, 
                          wbc  ::Int64,
                          δt   ::Float64, 
                          srδt ::Float64)
  @inbounds begin

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

    # acceptance rate
    normprop = 
      duoldnorm(λf,        λd1 - α*ped1, λd2 - α*ped2, ped1, ped2, σλ) -
      duoldnorm(λd1v_c[1], λd1 - α*ped1, λd2 - α*ped2, ped1, ped2, σλ)

    # fill fix and simulate unfix tree
    bm!(treepd1, λd1v_p, 1, lid1, α, σλ, δt, srδt)
    bm!(treepd2, λd2v_p, 1, lid2, α, σλ, δt, srδt)

    llrbm_d1, llrbd_d1 = llr_gbm_sep_f(treepd1, treecd1, α, σλ, δt, srδt)
    llrbm_d2, llrbd_d2 = llr_gbm_sep_f(treepd2, treecd2, α, σλ, δt, srδt)

    acr  = llrbd_d1 + llrbd_d2
    llr  = llrbm_d1 + llrbm_d2 + acr
    acr += normprop
  end

  return llr, acr
end




"""
    triad_lvupdate_trio!(treep::iTgbmce, 
                         treec::iTgbmce,
                         bbλp ::Array{Array{Float64,1},1}, 
                         bbλc ::Array{Array{Float64,1},1}, 
                         tsv  ::Array{Array{Float64,1},1}, 
                         llc  ::Float64,
                         pr   ::Int64,
                         d1   ::Int64,
                         d2   ::Int64,
                         α    ::Float64,
                         σλ   ::Float64,
                         μ    ::Float64,
                         δt   ::Float64, 
                         srδt ::Float64,
                         ter  ::BitArray{1},
                         icr  ::Bool,
                         wbc  ::Int64)

Make a `gbmce` trio proposal.
"""
function triad_lvupdate_trio!(treep::iTgbmce, 
                              treec::iTgbmce,
                              bbλp ::Array{Array{Float64,1},1}, 
                              bbλc ::Array{Array{Float64,1},1}, 
                              tsv  ::Array{Array{Float64,1},1}, 
                              llc  ::Float64,
                              pr   ::Int64,
                              d1   ::Int64,
                              d2   ::Int64,
                              α    ::Float64,
                              σλ   ::Float64,
                              μ    ::Float64,
                              δt   ::Float64, 
                              srδt ::Float64,
                              ter  ::BitArray{1},
                              icr  ::Bool,
                              wbc  ::Int64)
  @inbounds begin

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
        bm!(λprv_p, λpr, α, σλ, δt, tprv[2], srδt)
        lλp = λprv_p[lipr]
        bm!(λd1v_p, lλp, α, σλ, δt, td1v[2], srδt)
        bm!(λd2v_p, lλp, α, σλ, δt, td2v[2], srδt)

      else
        # if d1 is terminal
        lλp = duoprop(λpr + α*pepr, λd2 - α*ped2, pepr, ped2, σλ)

        # simulate fix tree vector
        bb!(λprv_p, λpr, lλp, tprv[2], σλ, δt, srδt)
        bm!(λd1v_p, lλp, α, σλ, δt, td1v[2], srδt)
        bb!(λd2v_p, lλp, λd2, td2v[2], σλ, δt, srδt)

      end
    elseif ter[2]
      # if d2 is terminal
      # node proposal
      lλp = duoprop(λpr + α*pepr, λd1 - α*ped1, pepr, ped1, σλ)

      # simulate fix tree vector
      bb!(λprv_p, λpr, lλp, tprv[2], σλ, δt, srδt)
      bb!(λd1v_p, lλp, λd1, td1v[2], σλ, δt, srδt)
      bm!(λd2v_p, lλp, α, σλ, δt, td2v[2], srδt)

    else
      # if no terminal branches involved
      # node proposal
      lλp  = trioprop(λpr + α*pepr, λd1 - α*ped1, λd2 - α*ped2, 
               pepr, ped1, ped2, σλ)

      # simulate fix tree vector
      bb!(λprv_p, λpr, lλp, tprv[2], σλ, δt, srδt)
      bb!(λd1v_p, lλp, λd1, td1v[2], σλ, δt, srδt)
      bb!(λd2v_p, lλp, λd2, td2v[2], σλ, δt, srδt)
    end

    # fill fix and simulate unfix tree
    bm!(treep,   λprv_p, 1, lipr, α, σλ, δt, srδt)
    bm!(treepd1, λd1v_p, 1, lid1, α, σλ, δt, srδt)
    bm!(treepd2, λd2v_p, 1, lid2, α, σλ, δt, srδt)

    ## make acceptance ratio 
    # estimate likelihoods
    llr, acr = llr_propr(treep, treepd1, treepd2, 
                         treec, treecd1, treecd2, 
                         α, σλ, μ, δt, srδt, icr, wbc)

    if -randexp() < acr
      llc += llr
      unsafe_copyto!(λprv_c, 1, λprv_p, 1, lipr)
      unsafe_copyto!(λd1v_c, 1, λd1v_p, 1, lid1)
      unsafe_copyto!(λd2v_c, 1, λd2v_p, 1, lid2)
      gbm_copy_f!(treec, treep)
      gbm_copy_f!(treecd1, treepd1)
      gbm_copy_f!(treecd2, treepd2)
    end
  end

  return llc
end




"""
    llr_propr(treep  ::iTgbmce,
              treepd1::iTgbmce,
              treepd2::iTgbmce,
              treec  ::iTgbmce,
              treecd1::iTgbmce,
              treecd2::iTgbmce,
              α      ::Float64,
              σλ     ::Float64,
              μ      ::Float64,
              δt     ::Float64,
              srδt   ::Float64,
              icr    ::Bool,
              wbc    ::Int64)

Return the likelihood and proposal ratio for birth-death gbm.
"""
function llr_propr(treep  ::iTgbmce,
                   treepd1::iTgbmce,
                   treepd2::iTgbmce,
                   treec  ::iTgbmce,
                   treecd1::iTgbmce,
                   treecd2::iTgbmce,
                   α      ::Float64,
                   σλ     ::Float64,
                   μ      ::Float64,
                   δt     ::Float64,
                   srδt   ::Float64,
                   icr    ::Bool,
                   wbc    ::Int64)

  llrbm_pr, llrbd_pr = llr_gbm_sep_f(treep,   treec,   α, σλ, δt, srδt)
  llrbm_d1, llrbd_d1 = llr_gbm_sep_f(treepd1, treecd1, α, σλ, δt, srδt)
  llrbm_d2, llrbd_d2 = llr_gbm_sep_f(treepd2, treecd2, α, σλ, δt, srδt)

  llrcond = 0.0

  if iszero(wbc)
    if icr 
      llrcond += cond_surv_crown(treep, μ) -
                 cond_surv_crown(treec, μ)
    else
      llrcond += cond_surv_stem(treep, μ) -
                 cond_surv_stem(treec, μ)
    end
  elseif isone(wbc) && icr
    llrcond += cond_surv_stem(treep, μ) -
               cond_surv_stem(treec, μ)
  end

  acr = llrbd_pr + llrbd_d1 + llrbd_d2 + llrcond
  llr = acr + llrbm_pr + llrbm_d1 + llrbm_d2

  return llr, acr
end



