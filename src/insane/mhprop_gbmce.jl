#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#




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
    gbm_update!(treep::iTgbmce,
                treec::iTgbmce,
                α    ::Float64,
                σλ   ::Float64,
                μ    ::Float64,
                llc  ::Float64,
                δt   ::Float64,
                srδt ::Float64,
                scond::Bool)

Do gbm updates for all internal nodes across a fixed branch.
"""
function gbm_update!(treep::iTgbmce,
                     treec::iTgbmce,
                     α    ::Float64,
                     σλ   ::Float64,
                     μ    ::Float64,
                     llc  ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)

  if isdefined(treec, :d1)
    llc = triad_lvupdate!(treep, treec, α, σλ, μ, llc, δt, srδt)

    if isfix(treec.d1) && isfix(treec.d2)
      return llc
    else
      llc = gbm_update!(treep.d1, treec.d1, α, σλ, μ, llc, δt, srδt)
      llc = gbm_update!(treep.d2, treec.d2, α, σλ, μ, llc, δt, srδt)
    end
  else
    llc = fs_lvupdate!(treep, treec, α, σλ, μ, llc, δt, srδt)
  end

  return llc
end




"""
    fs_lvupdate!(treep::iTgbmce,
                 treec::iTgbmce,
                 α    ::Float64,
                 σλ   ::Float64,
                 μ    ::Float64,
                 δt   ::Float64,
                 srδt ::Float64,
                 scond::Bool)

Make a `gbmce` trio proposal.
"""
function fs_lvupdate!(treep::iTgbmce,
                      treec::iTgbmce,
                      α    ::Float64,
                      σλ   ::Float64,
                      μ    ::Float64,
                      llc  ::Float64,
                      δt   ::Float64,
                      srδt ::Float64)

  @inbounds begin

    λpp  = lλ(treep)
    λpc  = lλ(treec)
    fdtp = fdt(treec)

    bm!(λpp, λpc[1], α, σλ, δt, fdtp, srδt)

    llrbm, llrbd = llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtp, srδt, false)

    acr = llrbd
    llr = acr + llrbm

    if -randexp() < acr
      llc += llr
      copyto!(λpc, λpp)
    end 
  end

  return llc
end




"""
    triad_lvupdate!(treep::iTgbmce,
                    treec::iTgbmce,
                    α    ::Float64,
                    σλ   ::Float64,
                    μ    ::Float64,
                    δt   ::Float64,
                    srδt ::Float64,
                    scond::Bool)

Make a `gbmce` trio proposal.
"""
function triad_lvupdate!(treep::iTgbmce,
                         treec::iTgbmce,
                         α    ::Float64,
                         σλ   ::Float64,
                         μ    ::Float64,
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

    epr  = e(treec)
    ed1  = e(treec.d1)
    ed2  = e(treec.d2)
    fdtp = fdt(treec)
    fdt1 = fdt(treec.d1)
    fdt2 = fdt(treec.d2)

    λpr  = λpc[1]
    λd1  = λ1c[end]
    λd2  = λ2c[end]

    # node proposal
    lλp = trioprop(λpr + α*epr, λd1 - α*ed1, λd2 - α*ed2, 
             epr, ed1, ed2, σλ)

    # simulate fix tree vector
    bb!(λpp, λpr, lλp, σλ, δt, fdtp, srδt)
    bb!(λ1p, lλp, λd1, σλ, δt, fdt1, srδt)
    bb!(λ2p, lλp, λd2, σλ, δt, fdt2, srδt)

    llrbm, llrbd = 
      llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtp, srδt, true)
    llrbm1, llrbd1 = 
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, !istip(treec.d1))
    llrbm2, llrbd2 = 
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, !istip(treec.d2))

    acr = llrbd + llrbd1 + llrbd2
    llr = acr + llrbm + llrbm1 + llrbm2

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

    λpr = λpc[1]
    λd1 = λ1c[end]
    λd2 = λ2c[end]

    if icr
      # node proposal
      lλp  = duoprop(λd1 - α*ed1, λd2 - α*ed2, ed1, ed2, σλ)
      lλrp = lλp
    else
      # node proposal
      lλp = trioprop(λpr + α*epr, λd1 - α*ed1, λd2 - α*ed2, 
             epr, ed1, ed2, σλ)
      # propose for root
      lλrp = rnorm(lλp - α*epr, sqrt(epr)*σλ)
    end

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





