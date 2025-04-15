#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#





"""
    _daughter_update!(ξ1  ::T,
                      λf  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64) where {T <: iTbdU}

Make a `fbdd` proposal for daughter for forward simulated fossil branch.
"""
function _daughter_update!(ξ1  ::T,
                           λf  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           δt  ::Float64,
                           srδt::Float64) where {T <: iTbdU}
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
    llrbm, llrbd, ssrλ =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)

    acr  = llrbd
    llr  = llrbm + acr
    srt  = sqrt(e1)
    acr += lrdnorm_bm_x(λf, λi, λ1 - α*e1, σλ * srt)
    drλ  = λi - λf
  end

  return llr, acr, drλ, ssrλ, λ1p
end




"""
    _daughters_update!(ξ1  ::T,
                       ξ2  ::T,
                       λf  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64) where {T <: iTbdU}

Make a `gbm-bd` proposal for daughters from forward simulated branch.
"""
function _daughters_update!(ξ1  ::T,
                            ξ2  ::T,
                            λf  ::Float64,
                            α   ::Float64,
                            σλ  ::Float64,
                            δt  ::Float64,
                            srδt::Float64) where {T <: iTbdU}
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
    llrbm1, llrbd1, ssrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
    llrbm2, llrbd2, ssrλ2 =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

    acr  = llrbd1 + llrbd2 + λf - λi
    llr  = llrbm1 + llrbm2 + acr
    acr += duoldnorm(λf, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ) -
           duoldnorm(λi, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ)
    drλ  = 2.0*(λi - λf)
    ssrλ = ssrλ1 + ssrλ2
  end

  return llr, acr, drλ, ssrλ, λ1p, λ2p
end




"""
    _stem_update!(ξi   ::T,
                  α    ::Float64,
                  σλ   ::Float64,
                  σμ   ::Float64,
                  llc  ::Float64,
                  dλ   ::Float64,
                  ssλ  ::Float64,
                  ssμ  ::Float64,
                  mc   ::Float64,
                  th   ::Float64,
                  δt   ::Float64,
                  srδt ::Float64,
                  lλxpr::Float64,
                  lμxpr::Float64)  where {T <: iTbdU}

Do gbm update for stem root.
"""
function _stem_update!(ξi      ::iTbd,
                       α       ::Float64,
                       σλ      ::Float64,
                       σμ      ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       dλ      ::Float64,
                       ssλ     ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       δt      ::Float64,
                       srδt    ::Float64,
                       λ0_prior::NTuple{2,Float64},
                       surv    ::Int64)

  @inbounds begin
    λc   = lλ(ξi)
    l    = lastindex(λc)
    λp   = Vector{Float64}(undef,l)
    λn   = λc[l]
    μr   = lμ(ξi)[1]
    el   = e(ξi)
    fdtp = fdt(ξi)

    # node proposal
    λr = duoprop(λn - α*el, λ0_prior[1], σλ^2*el, λ0_prior[2])

    # simulate fix tree vector
    bb!(λp, λr, λn, σλ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ =
      llr_gbm_b_sep(λp, λc, α, σλ, δt, fdtp, srδt, false)

    # log probability
    lU = -randexp()

    llr = llrbd

    if lU < llr + log(1000.0/mc)

      #survival
      mp   = m_surv_gbmbd(th, λr, μr, α, σλ, σμ, δt, srδt, 1_000, surv)
      llr += log(mp/mc)

      if lU < llr
        llc += llrbm + llr
        prc += llrdnorm_x(λr, λc[1], λ0_prior[1], λ0_prior[2])
        dλ  += λc[1] - λr
        ssλ += ssrλ
        mc   = mp
        unsafe_copyto!(λc, 1, λp, 1, l)
      end
    end
  end

  return llc, prc, dλ, ssλ, mc
end




"""
    _crown_update!(ξi      ::iTbd,
                   ξ1      ::iTbd,
                   ξ2      ::iTbd,
                   α       ::Float64,
                   σλ      ::Float64,
                   σμ      ::Float64,
                   llc     ::Float64,
                   dλ      ::Float64,
                   ssλ     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   δt      ::Float64,
                   srδt    ::Float64,
                   λ0_prior::NTuple{2,Float64},
                   surv    ::Int64)

Do `gbm-bd` update for crown root.
"""
function _crown_update!(ξi      ::iTbd,
                        ξ1      ::iTbd,
                        ξ2      ::iTbd,
                        α       ::Float64,
                        σλ      ::Float64,
                        σμ      ::Float64,
                        llc     ::Float64,
                        prc     ::Float64,
                        dλ      ::Float64,
                        ssλ     ::Float64,
                        mc      ::Float64,
                        th      ::Float64,
                        δt      ::Float64,
                        srδt    ::Float64,
                        λ0_prior::NTuple{2,Float64},
                        surv    ::Int64)

  @inbounds begin
    λpc  = lλ(ξi)
    λi   = λpc[1]
    μr   = lμ(ξi)[1]
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

    # node proposal
    λr = trioprop(λ1 - α*e1, λ2 - α*e2, λ0_prior[1], 
                  e1*σλ^2,     e2*σλ^2, λ0_prior[2])

    # simulate fix tree vector
    bb!(λ1p, λr, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λr, λ2, σλ, δt, fdt2, srδt)

    # log likelihood ratios
    llrbm1, llrbd1, ssrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
    llrbm2, llrbd2, ssrλ2  =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

    # log probability
    lU = -randexp()

    llr = llrbd1 + llrbd2

    if lU < llr + log(1000.0/mc)

      #survival
      mp   = m_surv_gbmbd(th, λr, μr, α, σλ, σμ, δt, srδt, 1_000, surv)
      llr += log(mp/mc)

      if lU < llr
        llc += llrbm1 + llrbm2 + llr
        prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
        dλ  += 2.0*(λi - λr)
        ssλ += ssrλ1 + ssrλ2
        mc   = mp
        fill!(λpc, λr)
        unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
        unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
      end
    end
  end

  return llc, prc, dλ, ssλ, mc
end




"""
    _update_gbm!(tree::T,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 llc ::Float64,
                 dλ  ::Float64,
                 ssλ ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 mσλ ::Float64,
                 mσμ ::Float64) where {T <: iTbdU}

Do `gbm-bd` updates on a decoupled tree recursively.
"""
function _update_gbm!(tree::T,
                      α   ::Float64,
                      σλ  ::Float64,
                      llc ::Float64,
                      dλ  ::Float64,
                      ssλ ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ter ::Bool) where {T <: iTbdU}

  if def1(tree)
    llc, dλ, ssλ =
      update_triad!(tree, α, σλ, llc, dλ, ssλ, δt, srδt)

    llc, dλ, ssλ =
      _update_gbm!(tree.d1, α, σλ, llc, dλ, ssλ, δt, srδt, ter)
    llc, dλ, ssλ =
      _update_gbm!(tree.d2, α, σλ, llc, dλ, ssλ, δt, srδt, ter)
  elseif !isfix(tree) || ter
    llc, dλ, ssλ = 
      update_tip!(tree, α, σλ, llc, dλ, ssλ, δt, srδt)
  end

  return llc, dλ, ssλ
end




"""
    update_tip!(tree::T,
                α   ::Float64,
                σλ  ::Float64,
                llc ::Float64,
                dλ  ::Float64,
                ssλ ::Float64,
                δt  ::Float64,
                srδt::Float64) where {T <: iTbdU}

Make a `gbm` tip proposal.
"""
function update_tip!(tree::T,
                     α   ::Float64,
                     σλ  ::Float64,
                     llc ::Float64,
                     dλ  ::Float64,
                     ssλ ::Float64,
                     δt  ::Float64,
                     srδt::Float64) where {T <: iTbdU}

  @inbounds begin

    λc   = lλ(tree)
    l    = lastindex(λc)
    fdtp = fdt(tree)
    λp   = Vector{Float64}(undef, l)

    bm!(λp, λc[1], α, σλ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ =
      llr_gbm_b_sep(λp, λc, α, σλ, δt, fdtp, srδt, false)

    if -randexp() < llrbd
      llc += llrbm + llrbd
      dλ  += λp[l] - λc[l]
      ssλ += ssrλ
      unsafe_copyto!(λc, 1, λp, 1, l)
    end
  end

  return llc, dλ, ssλ
end




"""
    update_duo!(λpc ::Vector{Float64},
                λ1c ::Vector{Float64},
                ep  ::Float64,
                e1  ::Float64,
                fdtp::Float64,
                fdt1::Float64,
                α   ::Float64,
                σλ  ::Float64,
                llc ::Float64,
                ssλ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a `gbm` duo proposal.
"""
function update_duo!(λpc ::Vector{Float64},
                     λ1c ::Vector{Float64},
                     ep  ::Float64,
                     e1  ::Float64,
                     fdtp::Float64,
                     fdt1::Float64,
                     α   ::Float64,
                     σλ  ::Float64,
                     llc ::Float64,
                     ssλ ::Float64,
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

    bb!(λpp, λp, λn, σλ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σλ, δt, fdt1, srδt)

    # log likelihood ratios
    llrbmp, llrbdp, ssrλp =
      llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtp, srδt, false)
    llrbm1, llrbd1, ssrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)

    acr = llrbdp + llrbd1

    if -randexp() < acr
      llc += llrbmp + llrbm1 + acr
      ssλ += ssrλp + ssrλ1
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
    end
  end

  return llc, ssλ
end




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
                  llc ::Float64,
                  dλ  ::Float64,
                  ssλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` trio proposal.
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
                       llc ::Float64,
                       dλ  ::Float64,
                       ssλ ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  @inbounds begin

    lp   = lastindex(λpc)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λpp  = Vector{Float64}(undef,lp)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    λp   = λpc[1]
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]

    # node proposal
    λn = trioprop(λp + α*ep, λ1 - α*e1, λ2 - α*e2, ep, e1, e2, σλ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, σλ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, σλ, δt, fdt2, srδt)

    llrbmp, llrbdp, ssrλp, =
      llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtp, srδt, true)
    llrbm1, llrbd1, ssrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
    llrbm2, llrbd2, ssrλ2 = 
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

    acr = llrbdp + llrbd1 + llrbd2

    if -randexp() < acr
      llc += llrbmp + llrbm1 + llrbm2 + acr
      dλ  += (λ1c[1] - λn)
      ssλ += ssrλp + ssrλ1 + ssrλ2
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, dλ, ssλ
end




"""
    update_triad!(tree::T,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  llc ::Float64,
                  dλ  ::Float64,
                  ssλ ::Float64,
                  ssμ ::Float64,
                  δt  ::Float64,
                  srδt::Float64) where {T <: iTbdU}

Make a `gbm` trio proposal.
"""
function update_triad!(tree::T,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       dλ  ::Float64,
                       ssλ ::Float64,
                       δt  ::Float64,
                       srδt::Float64) where {T <: iTbdU}

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

    llrbmp, llrbdp, ssrλp, =
      llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtp, srδt, true)
    llrbm1, llrbd1, ssrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
    llrbm2, llrbd2, ssrλ2 = 
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

    acr = llrbdp + llrbd1 + llrbd2

    if -randexp() < acr
      llc += llrbmp + llrbm1 + llrbm2 + acr
      dλ  += (λ1c[1] - λn)
      ssλ += ssrλp + ssrλ1 + ssrλ2
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, dλ, ssλ
end




