#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#




"""
    _daughters_update!(ξ1  ::iTct,
                       ξ2  ::iTct,
                       λf  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       ϵ   ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

Make a `gbm-ct` proposal for daughters from forwards simulated branch.
"""
function _daughters_update!(ξ1  ::iTct,
                            ξ2  ::iTct,
                            λf  ::Float64,
                            α   ::Float64,
                            σλ  ::Float64,
                            ϵ   ::Float64,
                            δt  ::Float64,
                            srδt::Float64)
  @inbounds begin

    λ1c  = lλ(ξ1)
    λ2c  = lλ(ξ2)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    λi   = λ1c[1]
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]
    e1   = e(ξ1)
    e2   = e(ξ2)
    fdt1 = fdt(ξ1)
    fdt2 = fdt(ξ2)

    bb!(λ1p, λf, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λf, λ2, σλ, δt, fdt2, srδt)

    # acceptance rate
    gp = duoldnorm(λf, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ) -
         duoldnorm(λi, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ)

    # log likelihood ratios
    llrbm1, llrct1, ssrλ1, Σrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, ϵ, δt, fdt1, srδt, false, false)
    llrbm2, llrct2, ssrλ2, Σrλ2 =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, ϵ, δt, fdt2, srδt, false, false)

    acr  = llrct1 + llrct2 + λf - λi 
    llr  = llrbm1 + llrbm2 + acr
    acr += gp
    drλ  = 2.0*(λi - λf)
    ssrλ = ssrλ1 + ssrλ2
    Σrλ  = Σrλ1 + Σrλ2
  end

  return llr, acr, drλ, ssrλ, Σrλ, λ1p, λ2p
end




"""
    _stem_update!(ξi   ::iTct,
                  α    ::Float64,
                  σλ   ::Float64,
                  ϵ    ::Float64,
                  llc  ::Float64,
                  dλ   ::Float64,
                  ssλ  ::Float64,
                  Σλ   ::Float64,
                  mc   ::Float64,
                  th   ::Float64,
                  δt   ::Float64,
                  srδt ::Float64,
                  lλxpr::Float64)

Do gbm update for stem root.
"""
function _stem_update!(ξi   ::iTct,
                       α    ::Float64,
                       σλ   ::Float64,
                       ϵ    ::Float64,
                       llc  ::Float64,
                       dλ   ::Float64,
                       ssλ  ::Float64,
                       Σλ   ::Float64,
                       mc   ::Float64,
                       th   ::Float64,
                       δt   ::Float64,
                       srδt ::Float64,
                       lλxpr::Float64)

  @inbounds begin
    λc   = lλ(ξi)
    l    = lastindex(λc)
    λp   = Vector{Float64}(undef,l)
    λn   = λc[l]
    el   = e(ξi)
    fdtp = fdt(ξi)

    # node proposal
    λr = rnorm(λn - α*el, σλ*sqrt(el))

    # prior ratio
    if λr > lλxpr
      return llc, dλ, ssλ, Σλ, mc
    end

    # simulate fix tree vector
    bb!(λp, λr, λn, σλ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ, Σrλ = llr_gbm_b_sep(λp, λc, α, σλ, ϵ, δt, fdtp, srδt,
      false, false)

    # survival
    mp  = m_surv_gbmct(th, λr, α, σλ, ϵ, δt, srδt, 1_000, true)
    llr = log(mp/mc)

    acr = llrbd + llr

    if -randexp() < acr
      llc += acr + llrbm
      dλ  += λc[1] - λr
      ssλ += ssrλ
      Σλ  += Σrλ
      mc   = mp
      unsafe_copyto!(λc, 1, λp, 1, l)
    end
  end

  return llc, dλ, ssλ, Σλ, mc
end




"""
    _crown_update!(ξi   ::iTct,
                   ξ1   ::iTct,
                   ξ2   ::iTct,
                   α    ::Float64,
                   σλ   ::Float64,
                   ϵ    ::Float64,
                   llc  ::Float64,
                   dλ   ::Float64,
                   ssλ  ::Float64,
                   Σλ   ::Float64,
                   mc   ::Float64,
                   th   ::Float64,
                   δt   ::Float64,
                   srδt ::Float64,
                   lλxpr::Float64)

Do gbm update for crown root.
"""
function _crown_update!(ξi   ::iTct,
                        ξ1   ::iTct,
                        ξ2   ::iTct,
                        α    ::Float64,
                        σλ   ::Float64,
                        ϵ    ::Float64,
                        llc  ::Float64,
                        dλ   ::Float64,
                        ssλ  ::Float64,
                        Σλ   ::Float64,
                        mc   ::Float64,
                        th   ::Float64,
                        δt   ::Float64,
                        srδt ::Float64,
                        lλxpr::Float64)

  @inbounds begin
    λpc  = lλ(ξi)
    λi   = λpc[1]
    λ1c  = lλ(ξ1)
    λ2c  = lλ(ξ2)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]
    e1   = e(ξ1)
    e2   = e(ξ2)
    fdt1 = fdt(ξ1)
    fdt2 = fdt(ξ2)

    σS = exp(randn())

    # node proposal
    λr = duoprop(λ1, λ2, e1, e2, σS)

    # prior ratio
    if λr > lλxpr
      return llc, dλ, ssλ, Σλ, mc
    end

    # simulate fix tree vector
    bb!(λ1p, λr, λ1, σS, δt, fdt1, srδt)
    bb!(λ2p, λr, λ2, σS, δt, fdt2, srδt)

    llr1, prr1, ssrλ1, Σrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, ϵ, σS, δt, fdt1, srδt, false, false)
    llr2, prr2, ssrλ2, Σrλ2 =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, ϵ, σS, δt, fdt2, srδt, false, false)

    # survival
    # mp  = m_surv_gbmct(th, λr, α, σλ, ϵ, δt, srδt, 1_000, false)
    mp = 1.0
    llr = log(mp/mc)

    llr += llr1 + llr2 + λr - λi

    if -randexp() < llr + prr1 + prr2
      llc += llr
      dλ  += 2.0*(λi - λr)
      ssλ += ssrλ1 + ssrλ2
      Σλ  += Σrλ1 + Σrλ2
      mc   = mp
      fill!(λpc, λr)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, dλ, ssλ, Σλ, mc
end




"""
    _update_gbm!(tree::iTct,
                 α   ::Float64,
                 σλ  ::Float64,
                 ϵ   ::Float64,
                 llc ::Float64,
                 dλ  ::Float64,
                 ssλ ::Float64,
                 Σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 ter ::Bool)

Do gbm updates on a decoupled tree recursively.
"""
function _update_gbm!(tree::iTct,
                      α   ::Float64,
                      σλ  ::Float64,
                      ϵ   ::Float64,
                      llc ::Float64,
                      dλ  ::Float64,
                      ssλ ::Float64,
                      Σλ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64)

  if def1(tree)
    llc, dλ, ssλ, Σλ = update_triad!(tree, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt)

    llc, dλ, ssλ, Σλ =
      _update_gbm!(tree.d1, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt)
    llc, dλ, ssλ, Σλ =
      _update_gbm!(tree.d2, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt)
  end

  return llc, dλ, ssλ, Σλ
end




# """
#     update_tip!(tree::iTct,
#                 α   ::Float64,
#                 σλ  ::Float64,
#                 ϵ   ::Float64,
#                 llc ::Float64,
#                 dλ  ::Float64,
#                 ssλ ::Float64,
#                 Σλ  ::Float64,
#                 δt  ::Float64,
#                 srδt::Float64)

# Make a `gbm` tip proposal.
# """
# function update_tip!(tree::iTct,
#                      α   ::Float64,
#                      σλ  ::Float64,
#                      ϵ   ::Float64,
#                      llc ::Float64,
#                      dλ  ::Float64,
#                      ssλ ::Float64,
#                      Σλ  ::Float64,
#                      δt  ::Float64,
#                      srδt::Float64)

#   @inbounds begin

#     λc   = lλ(tree)
#     l    = lastindex(λc)
#     fdtp = fdt(tree)
#     λp   = Vector{Float64}(undef, l)

#     bm!(λp, λc[1], α, σλ, δt, fdtp, srδt)

#     llrbm, llrbd, ssrλ, Σrλ = llr_gbm_b_sep(λp, λc, α, σλ, ϵ, δt, fdtp, srδt,
#       false, isextinct(tree))

#     acr = llrbd

#     if -randexp() < acr
#       llc += llrbm + acr
#       dλ  += λp[l] - λc[l]
#       ssλ += ssrλ
#       Σλ  += Σrλ
#       unsafe_copyto!(λc, 1, λp, 1, l)
#     end
#   end

#   return llc, dλ, ssλ, Σλ
# end




"""
    update_triad!(λpc ::Vector{Float64},
                  λ1c ::Vector{Float64},
                  λ2c ::Vector{Float64},
                  ep  ::Float64,
                  e1  ::Float64,
                  e2  ::Float64,
                  fdtp::Float64,
                  fdt1::Float64,
                  fdt2::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  llc ::Float64,
                  dλ  ::Float64,
                  ssλ ::Float64,
                  Σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` trio proposal for observed speciation event.
"""
function update_triad!(λpc ::Vector{Float64},
                       λ1c ::Vector{Float64},
                       λ2c ::Vector{Float64},
                       ep  ::Float64,
                       e1  ::Float64,
                       e2  ::Float64,
                       fdtp::Float64,
                       fdt1::Float64,
                       fdt2::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       ϵ   ::Float64,
                       llc ::Float64,
                       dλ  ::Float64,
                       ssλ ::Float64,
                       Σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  @inbounds begin

    lp  = lastindex(λpc)
    l1  = lastindex(λ1c)
    l2  = lastindex(λ2c)
    λpp = Vector{Float64}(undef,lp)
    λ1p = Vector{Float64}(undef,l1)
    λ2p = Vector{Float64}(undef,l2)
    λp  = λpc[1]
    λ1  = λ1c[l1]
    λ2  = λ2c[l2]

    σS = exp(randn())

    λn = trioprop(λp, λ1, λ2, ep, e1, e2, σS)

    # simulate fix tree vector
    bb!(λpp, λp, λn, σS, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σS, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, σS, δt, fdt2, srδt)

    llr, prr, ssrλ, Σrλ = llr_propr(λpp, λ1p, λ2p, λpc, λ1c, λ2c,
      α, σλ, ϵ, σS, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < llr + prr
      llc += llr
      dλ  += (λ1c[1] - λn)
      ssλ += ssrλ
      Σλ  += Σrλ
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, dλ, ssλ, Σλ
end




"""
    update_triad!(tree::iTct,
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  llc ::Float64,
                  dλ  ::Float64,
                  ssλ ::Float64,
                  Σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` trio proposal.
"""
function update_triad!(tree::iTct,
                       α   ::Float64,
                       σλ  ::Float64,
                       ϵ   ::Float64,
                       llc ::Float64,
                       dλ  ::Float64,
                       ssλ ::Float64,
                       Σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  @inbounds begin

    λpc  = lλ(tree)
    λ1c  = lλ(tree.d1)
    λ2c  = lλ(tree.d2)
    lp   = lastindex(λpc)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λpp  = Vector{Float64}(undef,lp)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    λp   = λpc[1]
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]
    ep   = e(tree)
    e1   = e(tree.d1)
    e2   = e(tree.d2)
    fdtp = fdt(tree)
    fdt1 = fdt(tree.d1)
    fdt2 = fdt(tree.d2)

    σS = exp(randn())

    λn = trioprop(λp, λ1, λ2, ep, e1, e2, σS)

    # simulate fix tree vector
    bb!(λpp, λp, λn, σS, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σS, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, σS, δt, fdt2, srδt)

    llr, prr, ssrλ, Σrλ = llr_propr(λpp, λ1p, λ2p, λpc, λ1c, λ2c,
      α, σλ, ϵ, σS, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < llr + prr
      llc += llr
      dλ  += (λ1c[1] - λn)
      ssλ += ssrλ
      Σλ  += Σrλ
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, dλ, ssλ, Σλ
end




"""
    llr_propr(λpp  ::Array{Float64,1},
              λ1p  ::Array{Float64,1},
              λ2p  ::Array{Float64,1},
              λpc  ::Array{Float64,1},
              λ1c  ::Array{Float64,1},
              λ2c  ::Array{Float64,1},
              α    ::Float64,
              σλ   ::Float64,
              ϵ    ::Float64,
              δt   ::Float64,
              fdtp::Float64,
              fdt1::Float64,
              fdt2::Float64,
              srδt ::Float64)

Return the likelihood and proposal ratio for gbm-ct.
"""
function llr_propr(λpp  ::Array{Float64,1},
                   λ1p  ::Array{Float64,1},
                   λ2p  ::Array{Float64,1},
                   λpc  ::Array{Float64,1},
                   λ1c  ::Array{Float64,1},
                   λ2c  ::Array{Float64,1},
                   α    ::Float64,
                   σλ   ::Float64,
                   ϵ    ::Float64,
                   σS   ::Float64,
                   δt   ::Float64,
                   fdtp::Float64,
                   fdt1::Float64,
                   fdt2::Float64,
                   srδt ::Float64)

  # log likelihood ratios
  llrp, prrp, ssrλp, Σrλp =
    llr_gbm_b_sep(λpp, λpc, α, σλ, ϵ, σS, δt, fdtp, srδt, true, false)
  llr1, prr1, ssrλ1, Σrλ1 =
    llr_gbm_b_sep(λ1p, λ1c, α, σλ, ϵ, σS, δt, fdt1, srδt, false, false)
  llr2, prr2, ssrλ2, Σrλ2 =
    llr_gbm_b_sep(λ2p, λ2c, α, σλ, ϵ, σS, δt, fdt2, srδt, false, false)

  llr  = llrp + llr1 + llr2
  prr  = prrp + prr1 + prr2
  ssrλ = ssrλp  + ssrλ1  + ssrλ2
  Σrλ  = Σrλp   + Σrλ1   + Σrλ2

  return llr, prr, ssrλ, Σrλ
end
