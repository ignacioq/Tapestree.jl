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
                     α    ::Float64,
                     σλ   ::Float64,
                     σμ   ::Float64,
                     δt   ::Float64, 
                     srδt ::Float64)

Make a `gbmbd` proposal for daughters of forwards simulated branch.
"""
function daughters_lprop!(treep::iTgbmbd, 
                          treec::iTgbmbd,
                          λf   ::Float64,
                          μf   ::Float64,
                          α    ::Float64,
                          σλ   ::Float64,
                          σμ   ::Float64,
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
    lμ1 = lμ(treecd1)
    λn  = lλ1[1]
    μn  = lμ1[1]
    λd1 = lλ1[end]
    λd2 = lλ(treecd2)[end]
    μd1 = lμ1[end]
    μd2 = lμ(treecd2)[end]

    bb!(lλ(treepd1), λf, λd1, lμ(treepd1), μf, μd1, 
        σλ, σμ, δt, fdt(treecd1), srδt)
    bb!(lλ(treepd2), λf, λd2, lμ(treepd2), μf, μd2, 
        σλ, σμ, δt, fdt(treecd2), srδt)

    # acceptance rate
    normprop = 
      duoldnorm(λf, λd1 - α*ed1, λd2 - α*ed2, ed1, ed2, σλ) -
      duoldnorm(λn, λd1 - α*ed1, λd2 - α*ed2, ed1, ed2, σλ) +
      duoldnorm(μf, μd1, μd2, ed1, ed2, σμ)                 -
      duoldnorm(μn, μd1, μd2, ed1, ed2, σμ)

    llrbm1, llrbd1 = 
      llr_gbm_b_sep(lλ(treepd1), lμ(treepd1), lλ1, lμ1, 
        α, σλ, σμ, δt, fdt(treecd1), srδt, !istip(treec.d1), isextinct(treecd1))
    llrbm2, llrbd2 = 
      llr_gbm_b_sep(lλ(treepd2), lμ(treepd2), lλ(treecd2), lμ(treecd2), 
        α, σλ, σμ, δt, fdt(treecd2), srδt, !istip(treec.d2), isextinct(treecd2))

    acr  = llrbd1 + llrbd2
    llr  = llrbm1 + llrbm2 + acr
    acr += normprop
  end

  return llr, acr
end




"""
    gbm_update!(treep::iTgbmbd,
                treec::iTgbmbd,
                α    ::Float64,
                σλ   ::Float64,
                σμ   ::Float64,
                llc  ::Float64,
                δt   ::Float64,
                srδt ::Float64)

Do gbm updates for all internal nodes across a fixed branch.
"""
function gbm_update!(treep::iTgbmbd,
                     treec::iTgbmbd,
                     α    ::Float64,
                     σλ   ::Float64,
                     σμ   ::Float64,
                     llc  ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)

  if isdefined(treec, :d1)
    llc = triad_lvupdate!(treep, treec, α, σλ, σμ, llc, δt, srδt)

    if isfix(treec.d1) && isfix(treec.d2)
      return llc
    else
      llc = gbm_update!(treep.d1, treec.d1, α, σλ, σμ, llc, δt, srδt)
      llc = gbm_update!(treep.d2, treec.d2, α, σλ, σμ, llc, δt, srδt)
    end
  else
    llc = fs_lvupdate!(treep, treec, α, σλ, σμ, llc, δt, srδt)
  end

  return llc
end




"""
    fs_lvupdate!(treep::iTgbmbd,
                 treec::iTgbmbd,
                 α    ::Float64,
                 σλ   ::Float64,
                 σμ   ::Float64,
                 δt   ::Float64,
                 srδt ::Float64)

Make a `gbmbd` trio proposal.
"""
function fs_lvupdate!(treep::iTgbmbd,
                      treec::iTgbmbd,
                      α    ::Float64,
                      σλ   ::Float64,
                      σμ   ::Float64,
                      llc  ::Float64,
                      δt   ::Float64,
                      srδt ::Float64)

  @inbounds begin

    λpp  = lλ(treep)
    μpp  = lμ(treep)
    λpc  = lλ(treec)
    μpc  = lμ(treec)
    fdtp = fdt(treec)

    bm!(λpp, μpp, λpc[1], μpc[1], α, σλ, σμ, δt, fdtp, srδt)

    llrbm, llrbd = llr_gbm_b_sep(λpp, μpp, λpc, μpc, α, σλ, σμ, δt, fdtp, srδt, 
        false, isextinct(treec))

    acr = llrbd
    llr = acr + llrbm

    if -randexp() < acr
      llc += llr
      copyto!(λpc, λpp)
      copyto!(μpc, μpp)
    end 
  end

  return llc
end




"""
    triad_lvupdate!(treep::iTgbmbd,
                    treec::iTgbmbd,
                    α    ::Float64,
                    σλ   ::Float64,
                    σμ   ::Float64,
                    llc  ::Float64,
                    δt   ::Float64,
                    srδt ::Float64)

Make a `gbmbd` trio proposal.
"""
function triad_lvupdate!(treep::iTgbmbd,
                         treec::iTgbmbd,
                         α    ::Float64,
                         σλ   ::Float64,
                         σμ   ::Float64,
                         llc  ::Float64,
                         δt   ::Float64,
                         srδt ::Float64)

  @inbounds begin

    λpp = lλ(treep)
    λ1p = lλ(treep.d1)
    λ2p = lλ(treep.d2)
    λpc = lλ(treec)
    λ1c = lλ(treec.d1)
    λ2c = lλ(treec.d2)
    μpp = lμ(treep)
    μ1p = lμ(treep.d1)
    μ2p = lμ(treep.d2)
    μpc = lμ(treec)
    μ1c = lμ(treec.d1)
    μ2c = lμ(treec.d2)

    epr  = e(treec)
    ed1  = e(treec.d1)
    ed2  = e(treec.d2)
    fdtp = fdt(treec)
    fdt1 = fdt(treec.d1)
    fdt2 = fdt(treec.d2)

    λpr  = λpc[1]
    λd1  = λ1c[end]
    λd2  = λ2c[end]
    μpr  = μpc[1]
    μd1  = μ1c[end]
    μd2  = μ2c[end]

    # node proposal
    lλp = trioprop(λpr + α*epr, λd1 - α*ed1, λd2 - α*ed2, epr, ed1, ed2, σλ)
    lμp = trioprop(μpr, μd1, μd2, epr, ed1, ed2, σμ)

    # simulate fix tree vector
    bb!(λpp, λpr, lλp, μpp, μpr, lμp, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, lλp, λd1, μ1p, lμp, μd1, σλ, σμ, δt, fdt1, srδt)
    bb!(λ2p, lλp, λd2, μ2p, lμp, μd2, σλ, σμ, δt, fdt2, srδt)

    llrbm, llrbd = 
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, 
        α, σλ, σμ, δt, fdtp, srδt, true, false)
    llrbm1, llrbd1 = 
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, 
        α, σλ, σμ, δt, fdt1, srδt, !istip(treec.d1), isextinct(treec.d1))
    llrbm2, llrbd2 = 
      llr_gbm_b_sep(λ2p, μ2p, λ2c, μ2c, 
        α, σλ, σμ, δt, fdt2, srδt, !istip(treec.d2), isextinct(treec.d2))

    acr = llrbd + llrbd1 + llrbd2
    llr = acr + llrbm + llrbm1 + llrbm2

    if -randexp() < acr
      llc += llr
      copyto!(λpc, λpp)
      copyto!(λ1c, λ1p)
      copyto!(λ2c, λ2p)
      copyto!(μpc, μpp)
      copyto!(μ1c, μ1p)
      copyto!(μ2c, μ2p)
    end
  end

  return llc
end




"""
    root_update!(treep ::iTgbmbd, 
                 treec ::iTgbmbd,
                 α     ::Float64,
                 σλ    ::Float64,
                 σμ    ::Float64,
                 δt    ::Float64, 
                 srδt  ::Float64,
                 lλmxpr::Float64,
                 icr   ::Bool)

Make a trio of Brownian motion MCMC updates when the root is involved.
"""
function root_update!(treep ::iTgbmbd, 
                      treec ::iTgbmbd,
                      α     ::Float64,
                      σλ    ::Float64,
                      σμ    ::Float64,
                      llc   ::Float64,
                      δt    ::Float64, 
                      srδt  ::Float64,
                      lλmxpr::Float64,
                      lμmxpr::Float64,
                      icr   ::Bool)

  @inbounds begin


    λpp = lλ(treep)
    λ1p = lλ(treep.d1)
    λ2p = lλ(treep.d2)
    λpc = lλ(treec)
    λ1c = lλ(treec.d1)
    λ2c = lλ(treec.d2)
    μpp = lμ(treep)
    μ1p = lμ(treep.d1)
    μ2p = lμ(treep.d2)
    μpc = lμ(treec)
    μ1c = lμ(treec.d1)
    μ2c = lμ(treec.d2)

    epr  = e(treec)
    ed1  = e(treec.d1)
    ed2  = e(treec.d2)
    fdtp = fdt(treec)
    fdt1 = fdt(treec.d1)
    fdt2 = fdt(treec.d2)

    λpr  = λpc[1]
    λd1  = λ1c[end]
    λd2  = λ2c[end]
    μpr  = μpc[1]
    μd1  = μ1c[end]
    μd2  = μ2c[end]

    if icr
      # node proposal
      lλp  = duoprop(λd1 - α*ed1, λd2 - α*ed2, ed1, ed2, σλ)
      lμp  = duoprop(μd1, μd2, ed1, ed2, σμ)
      lλrp = lλp
      lμrp = lμp
    else
      # node proposal
      lλp = trioprop(λpr + α*epr, λd1 - α*ed1, λd2 - α*ed2, 
             epr, ed1, ed2, σλ)
      lμp = trioprop(μpr, μd1, μd2, epr, ed1, ed2, σμ)

      # propose for root
      lλrp = rnorm(lλp - α*epr, sqrt(epr)*σλ)
      lμrp = rnorm(lμp, sqrt(epr)*σμ)
    end

    # simulate fix tree vector
    bb!(λpp, lλrp, lλp, μpp, lμrp, lμp, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, lλp, λd1, μ1p, lμp, μd1, σλ, σμ, δt, fdt1, srδt)
    bb!(λ2p, lλp, λd2, μ2p, lμp, μd2, σλ, σμ, δt, fdt2, srδt)

    llrbm, llrbd = 
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, 
        α, σλ, σμ, δt, fdtp, srδt, true, false)
    llrbm1, llrbd1 = 
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, 
        α, σλ, σμ, δt, fdt1, srδt, !istip(treec.d1), isextinct(treec.d1))
    llrbm2, llrbd2 = 
      llr_gbm_b_sep(λ2p, μ2p, λ2c, μ2c, 
        α, σλ, σμ, δt, fdt2, srδt, !istip(treec.d2), isextinct(treec.d2))

    if icr
      llrcond = cond_surv_crown(treep) -
                cond_surv_crown(treec)
    else
      llrcond = cond_surv_stem(treep) -
                cond_surv_stem(treec)
    end

    acr = llrbd + llrbd1 + llrbd2 + llrcond
    llr = acr + llrbm + llrbm1 + llrbm2

    # prior ratio
    if lλrp > lλmxpr
      acr += -Inf
    end
    if lμrp > lμmxpr
      acr += -Inf
    end

    if -randexp() < acr
      llc += llr
      copyto!(λpc, λpp)
      copyto!(λ1c, λ1p)
      copyto!(λ2c, λ2p)
      copyto!(μpc, μpp)
      copyto!(μ1c, μ1p)
      copyto!(μ2c, μ2p)
    end 
  end
  
  return llc
end



