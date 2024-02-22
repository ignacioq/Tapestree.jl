#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#



"""
    _daughter_update!(ξ1  ::T,
                      λf  ::Float64,
                      μf  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64) where {T <: iTbdU}

Make a gbm proposal for one daughter from forward simulated branch.
"""
function _daughter_update!(ξ1  ::T,
                           λf  ::Float64,
                           μf  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           σμ  ::Float64,
                           δt  ::Float64,
                           srδt::Float64) where {T <: iTbdU}
  @inbounds begin

    λ1c  = lλ(ξ1)
    μ1c  = lμ(ξ1)
    l1   = lastindex(λ1c)
    λ1p  = Vector{Float64}(undef,l1)
    μ1p  = Vector{Float64}(undef,l1)
    λi   = λ1c[1]
    λ1   = λ1c[l1]
    μi   = μ1c[1]
    μ1   = μ1c[l1]
    e1   = e(ξ1)
    fdt1 = fdt(ξ1)

    bb!(λ1p, λf, λ1, μ1p, μf, μ1, σλ, σμ, δt, fdt1, srδt)

    # log likelihood ratios
    llrbm1, llrbd1, ssrλ1, ssrμ1, irrλ1, irrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt, false, false)

    acr = llrbd1
    llr = llrbm1 + acr
    srt = sqrt(e1)
    acr += lrdnorm_bm_x(λf, λi, λ1 - α*e1, σλ * srt) + 
           lrdnorm_bm_x(μf, μi, μ1,        σμ * srt)
    drλ  = λi - λf
  end

  return llr, acr, drλ, ssrλ1, ssrμ1, irrλ, irrμ, λ1p, μ1p
end




"""
    _daughters_update!(ξ1  ::T,
                       ξ2  ::T,
                       λf  ::Float64,
                       μf  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64) where {T <: iTbdU}

Make a `gbm` proposal for daughters from forward simulated branch.
"""
function _daughters_update!(ξ1  ::T,
                            ξ2  ::T,
                            λf  ::Float64,
                            μf  ::Float64,
                            α   ::Float64,
                            σλ  ::Float64,
                            σμ  ::Float64,
                            δt  ::Float64,
                            srδt::Float64) where {T <: iTbdU}
  @inbounds begin

    λ1c  = lλ(ξ1)
    λ2c  = lλ(ξ2)
    μ1c  = lμ(ξ1)
    μ2c  = lμ(ξ2)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    μ1p  = Vector{Float64}(undef,l1)
    μ2p  = Vector{Float64}(undef,l2)
    λi   = λ1c[1]
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]
    μi   = μ1c[1]
    μ1   = μ1c[l1]
    μ2   = μ2c[l2]
    e1   = e(ξ1)
    e2   = e(ξ2)
    fdt1 = fdt(ξ1)
    fdt2 = fdt(ξ2)

    bb!(λ1p, λf, λ1, μ1p, μf, μ1, σλ, σμ, δt, fdt1, srδt)
    bb!(λ2p, λf, λ2, μ2p, μf, μ2, σλ, σμ, δt, fdt2, srδt)

    # log likelihood ratios
    llrbm1, llrbd1, ssrλ1, ssrμ1, irrλ1, irrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt, false, false)
    llrbm2, llrbd2, ssrλ2, ssrμ2, irrλ2, irrμ2 =
      llr_gbm_b_sep(λ2p, μ2p, λ2c, μ2c, α, σλ, σμ, δt, fdt2, srδt, false, false)

    acr  = llrbd1 + llrbd2 + λf - λi
    llr  = llrbm1 + llrbm2 + acr
    acr += duoldnorm(λf, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ) -
           duoldnorm(λi, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ) +
           duoldnorm(μf, μ1, μ2, e1, e2, σμ)               -
           duoldnorm(μi, μ1, μ2, e1, e2, σμ)
    drλ  = 2.0*(λi - λf)
    ssrλ = ssrλ1 + ssrλ2
    ssrμ = ssrμ1 + ssrμ2
    irrλ = irrλ1 + irrλ2
    irrμ = irrμ1 + irrμ2
  end

  return llr, acr, drλ, ssrλ, ssrμ, irrλ, irrμ, λ1p, λ2p, μ1p, μ2p
end




"""
    _stem_update!(ξi   ::T,
                  α    ::Float64,
                  σλ   ::Float64,
                  σμ   ::Float64,
                  llc  ::Float64,
                  ddλ  ::Float64,
                  ssλ  ::Float64,
                  ssμ  ::Float64,
                  irλ  ::Float64,
                  irμ  ::Float64,
                  mc   ::Float64,
                  th   ::Float64,
                  δt   ::Float64,
                  srδt ::Float64,
                  lλxpr::Float64,
                  lμxpr::Float64,
                  surv ::Int64) where {T <: iTbdU}

Do gbm update for stem root.
"""
function _stem_update!(ξi   ::T,
                       α    ::Float64,
                       σλ   ::Float64,
                       σμ   ::Float64,
                       llc  ::Float64,
                       ddλ  ::Float64,
                       ssλ  ::Float64,
                       ssμ  ::Float64,
                       irλ  ::Float64,
                       irμ  ::Float64,
                       mc   ::Float64,
                       th   ::Float64,
                       δt   ::Float64,
                       srδt ::Float64,
                       lλxpr::Float64,
                       lμxpr::Float64,
                       surv ::Int64) where {T <: iTbdU}

  @inbounds begin
    λc   = lλ(ξi)
    μc   = lμ(ξi)
    l    = lastindex(λc)
    λp   = Vector{Float64}(undef,l)
    μp   = Vector{Float64}(undef,l)
    λn   = λc[l]
    μn   = μc[l]
    el   = e(ξi)
    sqre = sqrt(el)
    fdtp = fdt(ξi)

    # node proposal
    λr = rnorm(λn - α*el, σλ*sqre)
    μr = rnorm(μn, σμ*sqre)

    # prior ratio
    if λr > lλxpr || μr > lμxpr
      return llc, ddλ, ssλ, ssμ, mc
    end

    # simulate fix tree vector
    bb!(λp, λr, λn, μp, μr, μn, σλ, σμ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ, ssrμ, irrλ, irrμ =
      llr_gbm_b_sep(λp, μp, λc, μc, α, σλ, σμ, δt, fdtp, srδt, false, false)

    # log probability
    lU = -randexp()

    llr = llrbd

    if lU < llr + log(1000.0/mc)

      #survival
      mp   = m_surv_gbmbd(th, λr, μr, α, σλ, σμ, δt, srδt, 1_000, surv)
      llr += log(mp/mc)

      if lU < llr
        llc += llrbm + llr
        ddλ += λc[1] - λr
        ssλ += ssrλ
        ssμ += ssrμ
        irλ += irrλ
        irμ += irrμ
        mc   = mp
        unsafe_copyto!(λc, 1, λp, 1, l)
        unsafe_copyto!(μc, 1, μp, 1, l)
      end
    end
  end

  return llc, ddλ, ssλ, ssμ, irλ, irμ, mc
end




"""
   _crown_update!(ξi   ::T,
                  ξ1   ::T,
                  ξ2   ::T,
                  α    ::Float64,
                  σλ   ::Float64,
                  σμ   ::Float64,
                  llc  ::Float64,
                  ddλ  ::Float64,
                  ssλ  ::Float64,
                  ssμ  ::Float64,
                  irλ  ::Float64, 
                  irμ  ::Float64,
                  mc   ::Float64,
                  th   ::Float64,
                  δt   ::Float64,
                  srδt ::Float64,
                  lλxpr::Float64,
                  lμxpr::Float64,
                  surv ::Int64) where {T <: iTbdU}

Do gbm update for crown root.
"""
function _crown_update!(ξi   ::T,
                        ξ1   ::T,
                        ξ2   ::T,
                        α    ::Float64,
                        σλ   ::Float64,
                        σμ   ::Float64,
                        llc  ::Float64,
                        ddλ  ::Float64,
                        ssλ  ::Float64,
                        ssμ  ::Float64,
                        irλ  ::Float64, 
                        irμ  ::Float64,
                        mc   ::Float64,
                        th   ::Float64,
                        δt   ::Float64,
                        srδt ::Float64,
                        lλxpr::Float64,
                        lμxpr::Float64,
                        surv ::Int64) where {T <: iTbdU}

  @inbounds begin
    λpc  = lλ(ξi)
    μpc  = lμ(ξi)
    λi   = λpc[1]
    μi   = μpc[1]
    λ1c  = lλ(ξ1)
    λ2c  = lλ(ξ2)
    μ1c  = lμ(ξ1)
    μ2c  = lμ(ξ2)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    μ1p  = Vector{Float64}(undef,l1)
    μ2p  = Vector{Float64}(undef,l2)
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]
    μ1   = μ1c[l1]
    μ2   = μ2c[l2]
    e1   = e(ξ1)
    e2   = e(ξ2)
    fdt1 = fdt(ξ1)
    fdt2 = fdt(ξ2)

    # node proposal
    λr = duoprop(λ1 - α*e1, λ2 - α*e2, e1, e2, σλ)
    μr = duoprop(μ1, μ2, e1, e2, σμ)

    # prior ratio
    if λr > lλxpr || μr > lμxpr 
      return llc, ddλ, ssλ, ssμ, irλ, irμ, mc
    end

    # simulate fix tree vector
    bb!(λ1p, λr, λ1, μ1p, μr, μ1, σλ, σμ, δt, fdt1, srδt)
    bb!(λ2p, λr, λ2, μ2p, μr, μ2, σλ, σμ, δt, fdt2, srδt)

    # log likelihood ratios
    llrbm1, llrbd1, ssrλ1, ssrμ1, irrλ1, irrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt, false, false)
    llrbm2, llrbd2, ssrλ2, ssrμ2, irrλ2, irrμ2 =
      llr_gbm_b_sep(λ2p, μ2p, λ2c, μ2c, α, σλ, σμ, δt, fdt2, srδt, false, false)

    # log probability
    lU = -randexp()

    llr = llrbd1 + llrbd2

    if lU < llr + log(1000.0/mc)

      #survival
      mp   = m_surv_gbmbd(th, λr, μr, α, σλ, σμ, δt, srδt, 1_000, surv)
      llr += log(mp/mc)

      if lU < llr
        llc += llrbm1 + llrbm2 + llr
        ddλ += 2.0*(λi - λr)
        ssλ += ssrλ1 + ssrλ2
        ssμ += ssrμ1 + ssrμ2
        irλ += irrλ1 + irrλ2
        irμ += irrμ1 + irrμ2
        mc   = mp
        fill!(λpc, λr)
        fill!(μpc, μr)
        unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
        unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
        unsafe_copyto!(μ1c, 1, μ1p, 1, l1)
        unsafe_copyto!(μ2c, 1, μ2p, 1, l2)
      end
    end
  end

  return llc, ddλ, ssλ, ssμ, irλ, irμ, mc
end




"""
    _update_gbm!(tree::iTbd,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 llc ::Float64,
                 ddλ ::Float64,
                 ssλ ::Float64,
                 ssμ ::Float64,
                 irλ ::Float64, 
                 irμ ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 ter ::Bool)

Do gbm updates on a decoupled tree recursively.
"""
function _update_gbm!(tree::iTbd,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      llc ::Float64,
                      ddλ ::Float64,
                      ssλ ::Float64,
                      ssμ ::Float64,
                      irλ ::Float64, 
                      irμ ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ter ::Bool)

  if def1(tree)
    llc, ddλ, ssλ, ssμ, irλ, irμ =
      update_triad!(tree, α, σλ, σμ, llc, ddλ, ssλ, ssμ, irλ, irμ, δt, srδt)

    llc, ddλ, ssλ, ssμ, irλ, irμ =
      _update_gbm!(tree.d1, α, σλ, σμ, llc, 
        ddλ, ssλ, ssμ, irλ, irμ, δt, srδt, ter)
    llc, ddλ, ssλ, ssμ, irλ, irμ =
      _update_gbm!(tree.d2, α, σλ, σμ, llc, 
        ddλ, ssλ, ssμ, irλ, irμ, δt, srδt, ter)
  elseif !isfix(tree) || ter
    llc, ddλ, ssλ, ssμ, irλ, irμ = 
      update_tip!(tree, α, σλ, σμ, llc, ddλ, ssλ, ssμ, irλ, irμ, δt, srδt)
  end

  return llc, ddλ, ssλ, ssμ, irλ, irμ
end




"""
    update_tip!(tree::T,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                llc ::Float64,
                ddλ ::Float64,
                ssλ ::Float64,
                ssμ ::Float64,
                irλ ::Float64, 
                irμ ::Float64,
                δt  ::Float64,
                srδt::Float64) where {T <: iTbdU}

Make a gbm tip proposal.
"""
function update_tip!(tree::T,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     llc ::Float64,
                     ddλ ::Float64,
                     ssλ ::Float64,
                     ssμ ::Float64,
                     irλ ::Float64, 
                     irμ ::Float64,
                     δt  ::Float64,
                     srδt::Float64) where {T <: iTbdU}

  @inbounds begin

    λc   = lλ(tree)
    μc   = lμ(tree)
    l    = lastindex(λc)
    fdtp = fdt(tree)
    λp   = Vector{Float64}(undef, l)
    μp   = Vector{Float64}(undef, l)

    bm!(λp, μp, λc[1], μc[1], α, σλ, σμ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ, ssrμ, irrλ, irrμ =
      llr_gbm_b_sep(λp, μp, λc, μc, α, σλ, σμ, δt, fdtp, srδt,
        false, isextinct(tree))

    if -randexp() < llrbd
      llc += llrbm + llrbd
      ddλ  += λp[l] - λc[l]
      ssλ += ssrλ
      ssμ += ssrμ
      irλ += irrλ 
      irμ += irrμ
      unsafe_copyto!(λc, 1, λp, 1, l)
      unsafe_copyto!(μc, 1, μp, 1, l)
    end
  end

  return llc, ddλ, ssλ, ssμ, irλ, irμ
end




"""
    update_duo!(λpc ::Vector{Float64},
                λ1c ::Vector{Float64},
                μpc ::Vector{Float64},
                μ1c ::Vector{Float64},
                ep  ::Float64,
                e1  ::Float64,
                fdtp::Float64,
                fdt1::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                llc ::Float64,
                ssλ ::Float64,
                ssμ ::Float64,
                irλ ::Float64, 
                irμ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a gbm duo proposal.
"""
function update_duo!(λpc ::Vector{Float64},
                     λ1c ::Vector{Float64},
                     μpc ::Vector{Float64},
                     μ1c ::Vector{Float64},
                     ep  ::Float64,
                     e1  ::Float64,
                     fdtp::Float64,
                     fdt1::Float64,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     llc ::Float64,
                     ssλ ::Float64,
                     ssμ ::Float64,
                     irλ ::Float64, 
                     irμ ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  @inbounds begin

    lp  = lastindex(λpc)
    l1  = lastindex(λ1c)
    λpp = Vector{Float64}(undef,lp)
    λ1p = Vector{Float64}(undef,l1)
    μpp = Vector{Float64}(undef,lp)
    μ1p = Vector{Float64}(undef,l1)
    λp  = λpc[1]
    λ1  = λ1c[l1]
    μp  = μpc[1]
    μ1  = μ1c[l1]

    # node proposal
    λn = duoprop(λp + α*ep, λ1 - α*e1, ep, e1, σλ)
    μn = duoprop(μp, μ1, ep, e1, σμ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, μpp, μp, μn, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, μ1p, μn, μ1, σλ, σμ, δt, fdt1, srδt)

    # log likelihood ratios
    llrbmp, llrbdp, ssrλp, ssrμp, irrλp, irrμp =
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, α, σλ, σμ, δt, fdtp, srδt,
        false, false)
    llrbm1, llrbd1, ssrλ1, ssrμ1, irrλ1, irrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt,
        false, false)

    acr = llrbdp + llrbd1

    if -randexp() < acr
      llc += llrbmp + llrbm1 + acr
      ssλ += ssrλp + ssrλ1
      ssμ += ssrμp + ssrμ1
      irλ += irrλp + irrλ1
      irμ += irrμp + irrμ1 
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(μpc, 1, μpp, 1, lp)
      unsafe_copyto!(μ1c, 1, μ1p, 1, l1)
    end
  end

  return llc, ssλ, ssμ, irλ, irμ
end




"""
    update_triad!(λpc ::Vector{Float64},
                  λ1c ::Vector{Float64},
                  λ2c ::Vector{Float64},
                  μpc ::Vector{Float64},
                  μ1c ::Vector{Float64},
                  μ2c ::Vector{Float64},
                  ep  ::Float64,
                  e1  ::Float64,
                  e2  ::Float64,
                  fdtp::Float64,
                  fdt1::Float64,
                  fdt2::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64,
                  ssμ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a gbm trio proposal.
"""
function update_triad!(λpc ::Vector{Float64},
                       λ1c ::Vector{Float64},
                       λ2c ::Vector{Float64},
                       μpc ::Vector{Float64},
                       μ1c ::Vector{Float64},
                       μ2c ::Vector{Float64},
                       ep  ::Float64,
                       e1  ::Float64,
                       e2  ::Float64,
                       fdtp::Float64,
                       fdt1::Float64,
                       fdt2::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       ssμ ::Float64,
                       irλ ::Float64, 
                       irμ ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  @inbounds begin

    lp  = lastindex(λpc)
    l1  = lastindex(λ1c)
    l2  = lastindex(λ2c)
    λpp = Vector{Float64}(undef,lp)
    λ1p = Vector{Float64}(undef,l1)
    λ2p = Vector{Float64}(undef,l2)
    μpp = Vector{Float64}(undef,lp)
    μ1p = Vector{Float64}(undef,l1)
    μ2p = Vector{Float64}(undef,l2)
    λp  = λpc[1]
    λi  = λ1c[1]
    λ1  = λ1c[l1]
    λ2  = λ2c[l2]
    μp  = μpc[1]
    μi  = μ1c[1]
    μ1  = μ1c[l1]
    μ2  = μ2c[l2]

    # node proposal
    λn = trioprop(λp + α*ep, λ1 - α*e1, λ2 - α*e2, ep, e1, e2, σλ)
    μn = trioprop(μp, μ1, μ2, ep, e1, e2, σμ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, μpp, μp, μn, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, μ1p, μn, μ1, σλ, σμ, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, μ2p, μn, μ2, σλ, σμ, δt, fdt2, srδt)

    # log likelihood ratios
    llrbmp, llrbdp, ssrλp, ssrμp, irrλp, irrμp =
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, α, σλ, σμ, δt, fdtp, srδt,
        true, false)
    llrbm1, llrbd1, ssrλ1, ssrμ1, irrλ1, irrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt,
        false, false)
    llrbm2, llrbd2, ssrλ2, ssrμ2, irrλ2, irrμ2 =
      llr_gbm_b_sep(λ2p, μ2p, λ2c, μ2c, α, σλ, σμ, δt, fdt2, srδt,
        false, false)

    acr = llrbdp + llrbd1 + llrbd2

    if -randexp() < acr
      llc += llrbmp + llrbm1 + llrbm2 + acr
      ddλ  += λi - λn
      ssλ += ssrλp + ssrλ1 + ssrλ2
      ssμ += ssrμp + ssrμ1 + ssrμ2
      irλ += irrλp + irrλ1 + irrλ2
      irμ += irrμp + irrμ1 + irrμ2
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
      unsafe_copyto!(μpc, 1, μpp, 1, lp)
      unsafe_copyto!(μ1c, 1, μ1p, 1, l1)
      unsafe_copyto!(μ2c, 1, μ2p, 1, l2)
      λi = λn
    end
  end

  return llc, ddλ, ssλ, ssμ, irλ, irμ, λi
end




"""
    update_triad!(tree::T,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  llc ::Float64,
                  ddλ  ::Float64,
                  ssλ ::Float64,
                  ssμ ::Float64,
                  irλ ::Float64,
                  irμ ::Float64,
                  δt  ::Float64,
                  srδt::Float64) where {T <: iTbdU}

Make a gbm trio proposal.
"""
function update_triad!(tree::T,
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       llc ::Float64,
                       ddλ  ::Float64,
                       ssλ ::Float64,
                       ssμ ::Float64,
                       irλ ::Float64,
                       irμ ::Float64,
                       δt  ::Float64,
                       srδt::Float64) where {T <: iTbdU}

  @inbounds begin

    λpc  = lλ(tree)
    λ1c  = lλ(tree.d1)
    λ2c  = lλ(tree.d2)
    μpc  = lμ(tree)
    μ1c  = lμ(tree.d1)
    μ2c  = lμ(tree.d2)
    lp   = lastindex(λpc)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λpp  = Vector{Float64}(undef,lp)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    μpp  = Vector{Float64}(undef,lp)
    μ1p  = Vector{Float64}(undef,l1)
    μ2p  = Vector{Float64}(undef,l2)
    λp   = λpc[1]
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]
    μp   = μpc[1]
    μ1   = μ1c[l1]
    μ2   = μ2c[l2]
    ep   = e(tree)
    e1   = e(tree.d1)
    e2   = e(tree.d2)
    fdtp = fdt(tree)
    fdt1 = fdt(tree.d1)
    fdt2 = fdt(tree.d2)

    # node proposal
    λn = trioprop(λp + α*ep, λ1 - α*e1, λ2 - α*e2, ep, e1, e2, σλ)
    μn = trioprop(μp, μ1, μ2, ep, e1, e2, σμ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, μpp, μp, μn, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, μ1p, μn, μ1, σλ, σμ, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, μ2p, μn, μ2, σλ, σμ, δt, fdt2, srδt)

    llrbmp, llrbdp, ssrλp, ssrμp, irrλp, irrμp =
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, α, σλ, σμ, δt, fdtp, srδt,
        true, false)
    llrbm1, llrbd1, ssrλ1, ssrμ1, irrλ1, irrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt,
        false, isextinct(tree.d1))
    llrbm2, llrbd2, ssrλ2, ssrμ2, irrλ2, irrμ2 =
      llr_gbm_b_sep(λ2p, μ2p, λ2c, μ2c, α, σλ, σμ, δt, fdt2, srδt,
        false, isextinct(tree.d2))

    acr = llrbdp + llrbd1 + llrbd2

    if -randexp() < acr
      llc += llrbmp + llrbm1 + llrbm2 + acr
      ddλ  += (λ1c[1] - λn)
      ssλ += ssrλp + ssrλ1 + ssrλ2
      ssμ += ssrμp + ssrμ1 + ssrμ2
      irλ += irrλp + irrλ1 + irrλ2
      irμ += irrμp + irrμ1 + irrμ2
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
      unsafe_copyto!(μpc, 1, μpp, 1, lp)
      unsafe_copyto!(μ1c, 1, μ1p, 1, l1)
      unsafe_copyto!(μ2c, 1, μ2p, 1, l2)
    end
  end

  return llc, ddλ, ssλ, ssμ, irλ, irμ
end



