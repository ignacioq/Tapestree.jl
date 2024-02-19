#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#




"""
    _daughter_update!(ξ1  ::iTct,
                      λf  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      ϵ   ::Float64,
                      δt  ::Float64,
                      srδt::Float64)

Make a `gbm-ct` proposal for daughters from forwards simulated branch.
"""
function _daughter_update!(ξ1  ::iTct,
                           λf  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           ϵ   ::Float64,
                           δt  ::Float64,
                           srδt::Float64)
  @inbounds begin

    λ1c  = lλ(ξ1)
    l1   = lastindex(λ1c)
    λ1p  = Vector{Float64}(undef,l1)
    λi   = λ1c[1]
    λ1   = λ1c[l1]
    e1   = e(ξ1)
    fdt1 = fdt(ξ1)

    bb!(λ1p, λf, λ1, σλ, δt, fdt1, srδt)

    # log likelihood ratios
    llrbm1, llrct1, ssrλ1, Σrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, ϵ, δt, fdt1, srδt, false, false)

    acr  = llrct1
    llr  = llrbm1 + acr
    acr += lrdnorm_bm_x(λf, λi, λ1 - α*e1, σλ * sqrt(e1))
    drλ  = λi - λf
  end

  return llr, acr, drλ, ssrλ1, Σrλ1, λ1p
end




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

    # log likelihood ratios
    llrbm1, llrct1, ssrλ1, Σrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, ϵ, δt, fdt1, srδt, false, false)
    llrbm2, llrct2, ssrλ2, Σrλ2 =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, ϵ, δt, fdt2, srδt, false, false)

    acr  = llrct1 + llrct2 + λf - λi
    llr  = llrbm1 + llrbm2 + acr
    acr += duoldnorm(λf, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ) -
           duoldnorm(λi, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ)
    drλ  = 2.0*(λi - λf)
    ssrλ = ssrλ1 + ssrλ2
    Σrλ  = Σrλ1 + Σrλ2
  end

  return llr, acr, drλ, ssrλ, Σrλ, λ1p, λ2p
end




"""
    _stem_update!(ξi      ::iTct,
                  α       ::Float64,
                  σλ      ::Float64,
                  ϵ       ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  dλ      ::Float64,
                  ssλ     ::Float64,
                  Σλ      ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  surv    ::Int64,
                  δt      ::Float64,
                  srδt    ::Float64,
                  λa_prior::NTuple{2,Float64})

Do gbm update for stem root.
"""
function _stem_update!(ξi      ::iTct,
                       α       ::Float64,
                       σλ      ::Float64,
                       ϵ       ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       dλ      ::Float64,
                       ssλ     ::Float64,
                       Σλ      ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       surv    ::Int64,
                       δt      ::Float64,
                       srδt    ::Float64,
                       λa_prior::NTuple{2,Float64})

  @inbounds begin
    lλc   = lλ(ξi)
    l    = lastindex(lλc)
    lλp   = Vector{Float64}(undef,l)
    lλn   = lλc[l]
    el   = e(ξi)
    fdtp = fdt(ξi)

    # node proposal
    lλr = rnorm(lλn - α*el, σλ*sqrt(el))

    # simulate fix tree vector
    bb!(lλp, lλr, lλn, σλ, δt, fdtp, srδt)

    llrbm, llrct, ssrλ, Σrλ = 
      llr_gbm_b_sep(lλp, lλc, α, σλ, ϵ, δt, fdtp, srδt, false, false)

    # log probability
    lU = -randexp()

    llr = llrct
    prr = llrdgamma(exp(lλp[1]), exp(lλc[1]), λa_prior[1], λa_prior[2])

    if lU < llr + prr + log(5_000.0/mc)

      # survival
      mp   = m_surv_gbmct(th, lλr, α, σλ, ϵ, δt, srδt, 5_000, surv)
      llr += log(mp/mc)

      if lU < llr + prr
        llc += llrbm + llr
        prc += prr
        dλ  += lλc[1] - lλr
        ssλ += ssrλ
        Σλ  += Σrλ
        mc   = mp
        unsafe_copyto!(lλc, 1, lλp, 1, l)
      end
    end
  end

  return llc, prc, dλ, ssλ, Σλ, mc
end




"""
    _crown_update!(ξi      ::iTct,
                   ξ1      ::iTct,
                   ξ2      ::iTct,
                   α       ::Float64,
                   σλ      ::Float64,
                   ϵ       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   dλ      ::Float64,
                   ssλ     ::Float64,
                   Σλ      ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   surv    ::Int64,
                   δt      ::Float64,
                   srδt    ::Float64,
                   λa_prior::NTuple{2,Float64})

Do gbm update for crown root.
"""
function _crown_update!(ξi      ::iTct,
                        ξ1      ::iTct,
                        ξ2      ::iTct,
                        α       ::Float64,
                        σλ      ::Float64,
                        ϵ       ::Float64,
                        llc     ::Float64,
                        prc     ::Float64,
                        dλ      ::Float64,
                        ssλ     ::Float64,
                        Σλ      ::Float64,
                        mc      ::Float64,
                        th      ::Float64,
                        surv    ::Int64,
                        δt      ::Float64,
                        srδt    ::Float64,
                        λa_prior::NTuple{2,Float64})

  @inbounds begin
    lλpc  = lλ(ξi)
    lλi   = lλpc[1]
    lλ1c  = lλ(ξ1)
    lλ2c  = lλ(ξ2)
    l1   = lastindex(lλ1c)
    l2   = lastindex(lλ2c)
    lλ1p  = Vector{Float64}(undef,l1)
    lλ2p  = Vector{Float64}(undef,l2)
    lλ1   = lλ1c[l1]
    lλ2   = lλ2c[l2]
    e1   = e(ξ1)
    e2   = e(ξ2)
    fdt1 = fdt(ξ1)
    fdt2 = fdt(ξ2)

    # node proposal
    lλr = duoprop(lλ1 - α*e1, lλ2 - α*e2, e1, e2, σλ)

    # simulate fix tree vector
    bb!(lλ1p, lλr, lλ1, σλ, δt, fdt1, srδt)
    bb!(lλ2p, lλr, lλ2, σλ, δt, fdt2, srδt)

    # log likelihood ratios
    llrbm1, llrct1, ssrλ1, Σrλ1 =
      llr_gbm_b_sep(lλ1p, lλ1c, α, σλ, ϵ, δt, fdt1, srδt, false, false)
    llrbm2, llrct2, ssrλ2, Σrλ2 =
      llr_gbm_b_sep(lλ2p, lλ2c, α, σλ, ϵ, δt, fdt2, srδt, false, false)

    # log probability
    lU = -randexp()

    llr = llrct1 + llrct2
    prr = llrdgamma(exp(lλr), exp(lλi), λa_prior[1], λa_prior[2])

    if lU < llr + prr + log(5_000.0/mc)

      # survival
      mp   = m_surv_gbmct(th, lλr, α, σλ, ϵ, δt, srδt, 5_000, surv)
      llr += log(mp/mc)

      if lU < llr + prr
        llc += llrbm1 + llrbm2 + llr
        prc += prr
        dλ  += 2.0*(lλi - lλr)
        ssλ += ssrλ1 + ssrλ2
        Σλ  += Σrλ1 + Σrλ2
        mc   = mp
        fill!(lλpc, lλr)
        unsafe_copyto!(lλ1c, 1, lλ1p, 1, l1)
        unsafe_copyto!(lλ2c, 1, lλ2p, 1, l2)
      end
    end
  end

  return llc, prc, dλ, ssλ, Σλ, mc
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
                      srδt::Float64,
                      ter  ::Bool)


  if def1(tree)
    llc, dλ, ssλ, Σλ =
      update_triad!(tree, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt)

    llc, dλ, ssλ, Σλ =
      _update_gbm!(tree.d1, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, ter)
    llc, dλ, ssλ, Σλ =
      _update_gbm!(tree.d2, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt, ter)
  elseif !isfix(tree) || ter

    llc, dλ, ssλ, Σλ =
      update_tip!(tree, α, σλ, ϵ, llc, dλ, ssλ, Σλ, δt, srδt)
  end

  return llc, dλ, ssλ, Σλ
end




"""
    update_tip!(tree::iTct,
                α   ::Float64,
                σλ  ::Float64,
                ϵ   ::Float64,
                llc ::Float64,
                dλ  ::Float64,
                ssλ ::Float64,
                Σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a `gbm` tip proposal.
"""
function update_tip!(tree::iTct,
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

    λc   = lλ(tree)
    l    = lastindex(λc)
    fdtp = fdt(tree)
    λp   = Vector{Float64}(undef, l)

    bm!(λp, λc[1], α, σλ, δt, fdtp, srδt)

    llrbm, llrct, ssrλ, Σrλ = llr_gbm_b_sep(λp, λc, α, σλ, ϵ, δt, fdtp, srδt,
      false, isextinct(tree))

    if -randexp() < llrct
      llc += llrbm + llrct
      dλ  += λp[l] - λc[l]
      ssλ += ssrλ
      Σλ  += Σrλ
      unsafe_copyto!(λc, 1, λp, 1, l)
    end
  end

  return llc, dλ, ssλ, Σλ
end




"""
    update_duo_ϵ!(λpc ::Vector{Float64},
                  λ1c ::Vector{Float64},
                  ep  ::Float64,
                  e1  ::Float64,
                  fdtp::Float64,
                  fdt1::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  llc ::Float64,
                  ssλ ::Float64,
                  Σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` duo proposal.
"""
function update_duo_ϵ!(λpc ::Vector{Float64},
                       λ1c ::Vector{Float64},
                       ep  ::Float64,
                       e1  ::Float64,
                       fdtp::Float64,
                       fdt1::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       ϵ   ::Float64,
                       llc ::Float64,
                       ssλ ::Float64,
                       Σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  @inbounds begin

    lp  = lastindex(λpc)
    l1  = lastindex(λ1c)
    λpp = Vector{Float64}(undef,lp)
    λ1p = Vector{Float64}(undef,l1)
    λp  = λpc[1]
    λ1  = λ1c[l1]

    # node proposal
    λn = duoprop(λp + α*ep, λ1 - α*e1, ep, e1, σλ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, σλ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σλ, δt, fdt1, srδt)

    llrbmp, llrctp, ssrλp, Σrλp =
      llr_gbm_b_sep(λpp, λpc, α, σλ, ϵ, δt, fdtp, srδt,
        false, false)
    llrbm1, llrct1, ssrλ1, Σrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, ϵ, δt, fdt1, srδt,
        false, false)

    acr = llrctp + llrct1

    if -randexp() < acr
      llc += llrbmp + llrbm1 + acr
      ssλ += ssrλp + ssrλ1
      Σλ  += Σrλp  + Σrλ1
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
    end
  end

  return llc, ssλ, Σλ
end




"""
    update_triad_ϵ!(λpc ::Vector{Float64},
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
                    srδt::Float64,
                    mσλ ::Float64)

Make a `gbm` trio proposal for observed speciation event.
"""
function update_triad_ϵ!(λpc ::Vector{Float64},
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

    # node proposal
    λn = trioprop(λp + α*ep, λ1 - α*e1, λ2 - α*e2, ep, e1, e2, σλ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, σλ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, σλ, δt, fdt2, srδt)

    llr, acr, ssrλ, Σrλ = llr_propr(λpp, λ1p, λ2p, λpc, λ1c, λ2c,
      α, σλ, ϵ, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < acr
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

    # node proposal
    λn = trioprop(λp + α*ep, λ1 - α*e1, λ2 - α*e2, ep, e1, e2, σλ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, σλ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, σλ, δt, fdt2, srδt)

    llrbmp, llrctp, ssrλp, Σrλp =
      llr_gbm_b_sep(λpp, λpc, α, σλ, ϵ, δt, fdtp, srδt,
        true, false)
    llrbm1, llrct1, ssrλ1, Σrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, ϵ, δt, fdt1, srδt,
        false, isextinct(tree.d1))
    llrbm2, llrct2, ssrλ2, Σrλ2 =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, ϵ, δt, fdt2, srδt,
        false, isextinct(tree.d2))

    acr = llrctp + llrct1 + llrct2

    if -randexp() < acr
      llc += llrbmp + llrbm1 + llrbm2 + acr
      dλ  += (λ1c[1] - λn)
      ssλ += ssrλp  + ssrλ1  + ssrλ2
      Σλ  += Σrλp   + Σrλ1   + Σrλ2
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
                   δt   ::Float64,
                   fdtp::Float64,
                   fdt1::Float64,
                   fdt2::Float64,
                   srδt ::Float64)

  # log likelihood ratios
  llrbmp, llrctp, ssrλp, Σrλp =
    llr_gbm_b_sep(λpp, λpc, α, σλ, ϵ, δt, fdtp, srδt, true, false)
  llrbm1, llrct1, ssrλ1, Σrλ1 =
    llr_gbm_b_sep(λ1p, λ1c, α, σλ, ϵ, δt, fdt1, srδt, false, false)
  llrbm2, llrct2, ssrλ2, Σrλ2 =
    llr_gbm_b_sep(λ2p, λ2c, α, σλ, ϵ, δt, fdt2, srδt, false, false)

  acr  = llrctp + llrct1 + llrct2
  llr  = llrbmp + llrbm1 + llrbm2 + acr
  ssrλ = ssrλp  + ssrλ1  + ssrλ2
  Σrλ  = Σrλp   + Σrλ1   + Σrλ2

  return llr, acr, ssrλ, Σrλ
end

