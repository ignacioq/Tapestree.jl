#=

Anagenetic GBM pure-birth MCMC MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 14 11 2021
=#




"""
    _daughters_update!(ξ1  ::iTxb,
                       ξ2  ::iTxb,
                       λf  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

Make a `gbmb` proposal for daughters from forwards simulated branch.
"""
function _daughters_update!(ξ1  ::iTxb,
                            ξ2  ::iTxb,
                            λf  ::Float64,
                            α   ::Float64,
                            σλ  ::Float64,
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
    llrbm1, llrb1, ssrλ1, irrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
    llrbm2, llrb2, ssrλ2, irrλ2 =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

    acr  = llrb1 + llrb2 + λf - λi
    llr  = llrbm1 + llrbm2 + acr
    acr += gp
    drλ  = 2.0*(λi - λf)
    ssrλ = ssrλ1 + ssrλ2
    irrλ = irrλ1 + irrλ2
  end

  return llr, acr, drλ, ssrλ, irrλ, λ1p, λ2p
end




"""
    _stem_update!(ξi      ::iTxb,
                  ασc     ::Float64, 
                  σσc     ::Float64, 
                  αλc     ::Float64, 
                  βλc     ::Float64, 
                  σλc     ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  dxs     ::Float64,
                  dxl     ::Float64,
                  ddx     ::Float64,
                  ddσ     ::Float64,
                  ssσ     ::Float64,
                  ddλ     ::Float64,
                  ssλ     ::Float64,
                  irλ     ::Float64,
                  δt      ::Float64,
                  srδt    ::Float64,
                  λ0_prior::NTuple{2,Float64})

Do diffusions' stem update.
"""
function _stem_update!(ξi      ::iTxb,
                       ασ      ::Float64, 
                       σσ      ::Float64, 
                       αλ      ::Float64, 
                       βλ      ::Float64, 
                       σλ      ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       dxs     ::Float64,
                       dxl     ::Float64,
                       ddx     ::Float64,
                       ddσ     ::Float64,
                       ssσ     ::Float64,
                       ddλ     ::Float64,
                       ssλ     ::Float64,
                       irλ     ::Float64,
                       δt      ::Float64,
                       srδt    ::Float64,
                       λ0_prior::NTuple{2,Float64})
  @inbounds begin
    xc   = xv(ξi)
    lσ2c = lσ2(ξi)
    lλc  = lλ(ξi)
    l    = lastindex(lλc)
    xp   = Vector{Float64}(undef,l)
    lσ2p = Vector{Float64}(undef,l)
    lλp  = Vector{Float64}(undef,l)
    x1   = xc[1]
    xn   = xc[l]
    lσ2n = lσ2c[l]
    lλn  = lλc[l]
    el   = e(ξi)
    fdtp = fdt(ξi)

    # rate path sample
    lσ2r = rnorm(lσ2n - ασ*el, γ*sqrt(el))
    bb!(lσ2p, lσ2r, lσ2n, σσ, δt, fdtp, srδt)

    llr, ssσr = llr_xb_σ(xc, ασ, lσ2p, lσ2c, δt, fdtp)

    if -randexp() < llr
      llc += llr
      ddσ += lσ2c[1] - lσ2r 
      ssσ += ssσr
      unsafe_copyto!(lσ2c, 1, lσ2p, 1, l)
    end

    # trait and speciation rate path sample
    xr  = rnorm(xn, sqrt(intσ2(lσ2c, δt, fdtp)))
    lλr = duoprop(lλn - βλ*(xn - xr), λ0_prior[1], σλ^2*el, λ0_prior[2])
    cbb!(xp, xr, xn, lσ2c, lλp, lλr, lλn, βλ, σλ, δt, fdt, srδt)

    llbmr, llbr, dxsr, dxlr, ssλr, irλr = 
      llr_xb_b_sep(vxp, vxc, vlσ2, lλp, lλc, 
        ασ, σσ, αλ, βλ, σλ, δt, fdtp, srδt, false)

    if -randexp() < llbr
      llc += llbmr + llbr
      prc += llrdnorm_x(λr, lλc[1], λ0_prior[1], λ0_prior[2])
      dxs += dxsr
      dxl += dxlr
      ddx += xc[1]  - xr
      ddλ += lλc[1] - λr
      ssλ += ssλr
      irλ += irλr
      unsafe_copyto!(xc,  1, xp,  1, l)
      unsafe_copyto!(lλc, 1, lλp, 1, l)
    end
  end

  return llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ
end




"""
    _crown_update!(ξi      ::iTxb,
                   ξ1      ::iTxb,
                   ξ2      ::iTxb,
                   ασ      ::Float64, 
                   σσ      ::Float64, 
                   αλ      ::Float64, 
                   βλ      ::Float64, 
                   σλ      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   dxs     ::Float64,
                   dxl     ::Float64,
                   ddx     ::Float64,
                   ddσ     ::Float64,
                   ssσ     ::Float64,
                   ddλ     ::Float64,
                   ssλ     ::Float64,
                   irλ     ::Float64,
                   δt      ::Float64,
                   srδt    ::Float64,
                   λ0_prior::NTuple{2,Float64})

Do diffusions' crown update.
"""
function _crown_update!(ξi      ::iTxb,
                        ξ1      ::iTxb,
                        ξ2      ::iTxb,
                        ασ      ::Float64, 
                        σσ      ::Float64, 
                        αλ      ::Float64, 
                        βλ      ::Float64, 
                        σλ      ::Float64,
                        llc     ::Float64,
                        prc     ::Float64,
                        dxs     ::Float64,
                        dxl     ::Float64,
                        ddx     ::Float64,
                        ddσ     ::Float64,
                        ssσ     ::Float64,
                        ddλ     ::Float64,
                        ssλ     ::Float64,
                        irλ     ::Float64,
                        δt      ::Float64,
                        srδt    ::Float64,
                        λ0_prior::NTuple{2,Float64})

  @inbounds begin
    xac   = xv(ξi)
    x1c   = xv(ξ1)
    x2c   = xv(ξ2)
    lσ2ac = lσ2(ξi)
    lσ21c = lσ2(ξi)
    lσ22c = lσ2(ξi)
    lλac  = lλ(ξi)
    lλ1c  = lλ(ξ1)
    lλ2c  = lλ(ξ2)
    l1    = lastindex(lλ1c)
    l2    = lastindex(lλ2c)
    x1p   = Vector{Float64}(undef,l1)
    x2p   = Vector{Float64}(undef,l2)
    lσ21p = Vector{Float64}(undef,l1)
    lσ22p = Vector{Float64}(undef,l2)
    lλ1p  = Vector{Float64}(undef,l1)
    lλ2p  = Vector{Float64}(undef,l2)
    x1f   = x1c[l1]
    x2f   = x2c[l2]
    lσ21f = lσ21c[l1]
    lσ22f = lσ22c[l2]
    lλ1f  = lλ1c[l1]
    lλ2f  = lλ2c[l2]
    e1    = e(ξ1)
    e2    = e(ξ2)
    fdt1  = fdt(ξ1)
    fdt2  = fdt(ξ2)

    # rate path sample
    lσ2n = duoprop(lσ21f - ασ*e1, lσ22f - ασ*e2, σσ^2*e1, σσ^2*e2)
    bb!(lσ21p, lσ2n, lσ21f, σσ, δt, fdt1, srδt)
    bb!(lσ22p, lσ2n, lσ22f, σσ, δt, fdt2, srδt)

    ll1r, ssσ1r = llr_xb_σ(x1c, ασ, lσ21p, lσ21c, δt, fdt1)
    ll2r, ssσ2r = llr_xb_σ(x2c, ασ, lσ22p, lσ22c, δt, fdt2)

    llr = ll1r + ll2r

    if -randexp() < llr
      llc += llr
      ssσ += ssσ1r + ssσ2r
      ddσ += 2.0*(lσ2ac[1] - lσ2n)
      unsafe_copyto!(lσ21c, 1, lσ21p, 1, l1)
      unsafe_copyto!(lσ22c, 1, lσ22p, 1, l2)
      fill!(lσ2ac, lσ2n)
    end




    """
    here
    """






    # node proposal
    λr = trioprop(λ1 - α*e1, λ2 - α*e2, λ0_prior[1], 
                  e1*σλ^2,     e2*σλ^2, λ0_prior[2])

    # simulate fix tree vector
    bb!(λ1p, λr, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λr, λ2, σλ, δt, fdt2, srδt)

    # log likelihood ratios
    llrbm1, llrb1, ssrλ1, irrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
    llrbm2, llrb2, ssrλ2, irrλ2 =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

    acr  = llrb1 + llrb2

    if -randexp() < acr
      llc += llrbm1 + llrbm2 + acr
      prc += llrdnorm_x(λr, λ1c[1], λ0_prior[1], λ0_prior[2])
      ddλ += 2.0*(λ1c[1] - λr)
      ssλ += ssrλ1 + ssrλ2
      irλ += irrλ1 + irrλ2
      fill!(λpc, λr)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, prc, ddλ, ssλ, irλ
end




"""
    _update_gbm!(tree::iTxb,
                 α   ::Float64,
                 σλ  ::Float64,
                 llc ::Float64,
                 ddλ ::Float64,
                 ssλ ::Float64,
                 irλ ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 ter ::Bool)

Do gbm updates on a decoupled tree recursively.
"""
function _update_gbm!(tree::iTxb,
                      α   ::Float64,
                      σλ  ::Float64,
                      llc ::Float64,
                      ddλ ::Float64,
                      ssλ ::Float64,
                      irλ ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ter ::Bool)

  if def1(tree)

    llc, ddλ, ssλ, irλ = 
      update_triad!(tree, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt)

    llc, ddλ, ssλ, irλ =
      _update_gbm!(tree.d1, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt, ter)
    llc, ddλ, ssλ, irλ =
      _update_gbm!(tree.d2, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt, ter)
  elseif !isfix(tree) || ter
    llc, ddλ, ssλ, irλ = 
      update_tip!(tree, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt)
  end

  return llc, ddλ, ssλ, irλ
end




"""
    update_tip!(tree::iTxb,
                α   ::Float64,
                σλ  ::Float64,
                llc ::Float64,
                ddλ ::Float64,
                ssλ ::Float64,
                irλ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a `gbm` tip proposal.
"""
function update_tip!(tree::iTxb,
                     α   ::Float64,
                     σλ  ::Float64,
                     llc ::Float64,
                     ddλ ::Float64,
                     ssλ ::Float64,
                     irλ ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  @inbounds begin

    λc   = lλ(tree)
    l    = lastindex(λc)
    fdtp = fdt(tree)
    λp   = Vector{Float64}(undef, l)

    bm!(λp, λc[1], α, σλ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ, irrλ = 
      llr_gbm_b_sep(λp, λc, α, σλ, δt, fdtp, srδt, false)

    acr = llrbd

    if -randexp() < acr
      llc += llrbm + acr
      ddλ += λp[l] - λc[l]
      ssλ += ssrλ
      irλ += irrλ
      unsafe_copyto!(λc, 1, λp, 1, l)
    end
  end

  return llc, ddλ, ssλ, irλ
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
                  ddλ ::Float64,
                  ssλ ::Float64,
                  irλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` trio proposal.
"""
function update_triad_b!(λpc ::Vector{Float64},
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
                       ddλ ::Float64,
                       ssλ ::Float64,
                       irλ ::Float64,
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

    llr, acr, ssrλ, irrλ = llr_propr_b( λpp, λ1p, λ2p, λpc, λ1c, λ2c,
      α, σλ, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < acr
      llc += llr
      ddλ += (λ1c[1] - λn)
      ssλ += ssrλ
      irλ += irrλ
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, ddλ, ssλ, irλ
end




"""
    update_triad!(tree::iTxb,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64,
                  irλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` trio proposal.
"""
function update_triad!(tree::iTxb,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       irλ ::Float64,
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

    llr, acr, ssrλ, irrλ = llr_propr_b( λpp, λ1p, λ2p, λpc, λ1c, λ2c,
      α, σλ, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < acr
      llc += llr
      ddλ += (λ1c[1] - λn)
      ssλ += ssrλ
      irλ += irrλ
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, ddλ, ssλ, irλ
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
              δt   ::Float64,
              fdtpr::Float64,
              fdtd1::Float64,
              fdtd2::Float64,
              srδt ::Float64)

Return the likelihood and proposal ratio for pure-birth gbm.
"""
function llr_propr_b( λpp  ::Array{Float64,1},
                   λ1p  ::Array{Float64,1},
                   λ2p  ::Array{Float64,1},
                   λpc  ::Array{Float64,1},
                   λ1c  ::Array{Float64,1},
                   λ2c  ::Array{Float64,1},
                   α    ::Float64,
                   σλ   ::Float64,
                   δt   ::Float64,
                   fdtp::Float64,
                   fdt1::Float64,
                   fdt2::Float64,
                   srδt ::Float64)

  # log likelihood ratios
  llrbmp, llrbp, ssrλp, irrλp = 
    llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtp, srδt, true)
  llrbm1, llrb1, ssrλ1, irrλ1 = 
    llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
  llrbm2, llrb2, ssrλ2, irrλ2 = 
    llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

  acr  = llrbp + llrb1 + llrb2
  llr  = llrbmp + llrbm1 + llrbm2 + acr
  ssrλ = ssrλp + ssrλ1 + ssrλ2
  irrλ = irrλp + irrλ1 + irrλ2

  return llr, acr, ssrλ, irrλ
end