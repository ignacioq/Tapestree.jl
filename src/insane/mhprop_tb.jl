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

Make a `xb` proposal for daughters from forwards simulated branch.
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
    _update_stem!(ξi      ::iTxb,
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

Do diffusions' stem update.
"""
function _update_stem!(ξi      ::iTxb,
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
    xf   = xc[l]
    lσ2n = lσ2c[l]
    lλf  = lλc[l]
    el   = e(ξi)
    fdtp = fdt(ξi)

    # rate path sample
    lσ2r = rnorm(lσ2n - ασ*el, σσ*sqrt(el))
    bb!(lσ2p, lσ2r, lσ2n, σσ, δt, fdtp, srδt)

    llr, ssσr = llr_xb_σ(xc, ασ, lσ2p, lσ2c, δt, fdtp)

    if -randexp() < llr
      llc += llr
      ddσ += lσ2c[1] - lσ2r 
      ssσ += ssσr
      unsafe_copyto!(lσ2c, 1, lσ2p, 1, l)
    end

    # trait and speciation rate path sample
    xn  = rnorm(xf, sqrt(intσ2(lσ2c, δt, fdtp)))
    lλn = duoprop(lλf - αλ*el - βλ*(xf - xn), λ0_prior[1], σλ^2*el, λ0_prior[2])
    cbb!(xp, xn, xf, lσ2c, lλp, lλn, lλf, βλ, σλ, δt, fdtp, srδt)

    llbmr, llbr, dxsr, dxlr, ssλr, irλr = 
      llr_xb_b_sep(xp, xc, lσ2c, lλp, lλc, 
        ασ, σσ, αλ, βλ, σλ, δt, fdtp, srδt, false)

    if -randexp() < llbr
      llc += llbmr + llbr
      prc += llrdnorm_x(lλn, lλc[1], λ0_prior[1], λ0_prior[2])
      dxs += dxsr
      dxl += dxlr
      ddx += xc[1]  - xn
      ddλ += lλc[1] - lλn
      ssλ += ssλr
      irλ += irλr
      unsafe_copyto!(xc,  1, xp,  1, l)
      unsafe_copyto!(lλc, 1, lλp, 1, l)
    end
  end

  return llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ
end




"""
    _update_crown!(ξi      ::iTxb,
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
function _update_crown!(ξi      ::iTxb,
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
    xac,     x1c,   x2c =  xv(ξi),  xv(ξ1),  xv(ξ2)
    lσ2ac, lσ21c, lσ22c = lσ2(ξi), lσ2(ξ1), lσ2(ξ2)
    lλac,   lλ1c,  lλ2c =  lλ(ξi),  lλ(ξ1),  lλ(ξ2)
    l1, l2 = lastindex(lλ1c), lastindex(lλ2c)
    x1p,     x2p = Vector{Float64}(undef,l1), Vector{Float64}(undef,l2)
    lσ21p, lσ22p = Vector{Float64}(undef,l1), Vector{Float64}(undef,l2)
    lλ1p,   lλ2p = Vector{Float64}(undef,l1), Vector{Float64}(undef,l2)
    lσ21f, lσ22f = lσ21c[l1], lσ22c[l2]
    x1f,     x2f =  x1c[l1],  x2c[l2]
    lλ1f,   lλ2f = lλ1c[l1], lλ2c[l2]
    e1, e2, fdt1, fdt2  = e(ξ1), e(ξ2), fdt(ξ1), fdt(ξ2)

    # trait rate path sample
    lσ2n = duoprop(lσ21f - ασ*e1, lσ22f - ασ*e2, σσ^2*e1, σσ^2*e2)
    bb!(lσ21p, lσ2n, lσ21f, σσ, δt, fdt1, srδt)
    bb!(lσ22p, lσ2n, lσ22f, σσ, δt, fdt2, srδt)

    llσx1r, llσσ1r, ssσ1r = llr_xb_σ(x1c, ασ, σσ, lσ21p, lσ21c, δt, fdt1)
    llσx2r, llσσ2r, ssσ2r = llr_xb_σ(x2c, ασ, σσ, lσ22p, lσ22c, δt, fdt2)

    llr = llσx1r + llσx2r

    if -randexp() < llr
      llc += llr + llσσ1r + llσσ2r
      ssσ += ssσ1r + ssσ2r
      ddσ += 2.0*(lσ2ac[1] - lσ2n)
      unsafe_copyto!(lσ21c, 1, lσ21p, 1, l1)
      unsafe_copyto!(lσ22c, 1, lσ22p, 1, l2)
      fill!(lσ2ac, lσ2n)
    end

    # trait and speciation path samples
    xn  = duoprop(x1f, x2f, intσ2(lσ21c, δt, fdt1), intσ2(lσ22c, δt, fdt2))
    lλn = trioprop(lλ1f - αλ*e1 - βλ*(x1f - xn), 
                   lλ2f - αλ*e2 - βλ*(x2f - xn), 
                   λ0_prior[1],
                   σλ^2*e1, σλ^2*e2, λ0_prior[2])

    cbb!(x1p, xn, x1f, lσ21c, lλ1p, lλn, lλ1f, βλ, σλ, δt, fdt1, srδt)
    cbb!(x2p, xn, x2f, lσ22c, lλ2p, lλn, lλ2f, βλ, σλ, δt, fdt2, srδt)

    # likelihood ratio
    llbm1r, llb1r, dxs1r, dxl1r, ssλ1r, irλ1r = 
      llr_xb_b_sep(x1p, x1c, lσ21c, lλ1p, lλ1c, 
        ασ, σσ, αλ, βλ, σλ, δt, fdt1, srδt, false)
    llbm2r, llb2r, dxs2r, dxl2r, ssλ2r, irλ2r = 
      llr_xb_b_sep(x2p, x2c, lσ22c, lλ2p, lλ2c, 
        ασ, σσ, αλ, βλ, σλ, δt, fdt2, srδt, false)

    llr = llb1r + llb2r

    if -randexp() < llr
      llc += llbm1r + llbm2r + llr
      prc += llrdnorm_x(lλn, lλ1c[1], λ0_prior[1], λ0_prior[2])
      dxs += dxs1r + dxs2r
      dxl += dxl1r + dxl2r
      ddx += 2.0*(x1c[1]  - xn)
      ddλ += 2.0*(lλ1c[1] - lλn)
      ssλ += ssλ1r + ssλ2r
      irλ += irλ1r + irλ2r
      fill!(xac, xn)
      unsafe_copyto!(x1c,  1, x1p,  1, l1)
      unsafe_copyto!(x2c,  1, x2p,  1, l2)
      fill!(lλac, lλn)
      unsafe_copyto!(lλ1c, 1, lλ1p, 1, l1)
      unsafe_copyto!(lλ2c, 1, lλ2p, 1, l2)
    end
  end

  return llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ
end




"""
    _update_node!(tree::iTxb,
                  ασ  ::Float64, 
                  σσ  ::Float64, 
                  αλ  ::Float64, 
                  βλ  ::Float64, 
                  σλ  ::Float64,
                  llc ::Float64,
                  dxs ::Float64,
                  dxl ::Float64,
                  ddx ::Float64,
                  ddσ ::Float64,
                  ssσ ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64,
                  irλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  ter ::Bool)

Perform xb node updates recursively.
"""
function _update_node!(tree::iTxb,
                       xavg::Float64,
                       xstd::Float64,
                       ασ  ::Float64, 
                       σσ  ::Float64, 
                       αλ  ::Float64, 
                       βλ  ::Float64, 
                       σλ  ::Float64,
                       llc ::Float64,
                       dxs ::Float64,
                       dxl ::Float64,
                       ddx ::Float64,
                       ddσ ::Float64,
                       ssσ ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       irλ ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       ter ::Bool)

  if def1(tree)
    llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ = 
      update_triad!(tree, ασ, σσ, αλ, βλ, σλ, 
        llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt)

    llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ =
      _update_node!(tree.d1, xavg, xstd, ασ, σσ, αλ, βλ, σλ,
        llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, ter)
    llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ =
      _update_node!(tree.d2, xavg, xstd, ασ, σσ, αλ, βλ, σλ,
        llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, ter)
  else
    if !isfix(tree)
      llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ = 
        update_tip!(tree, NaN, NaN, ασ, σσ, αλ, βλ, σλ,
          llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt)
    else
      if ter
        llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ = 
          update_tip!(tree, xavg, xstd, ασ, σσ, αλ, βλ, σλ,
            llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt)
      end
    end
  end

  return llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ
end




"""
    update_tip!(tree::iTxb,
                xavg::Float64, 
                xstd::Float64,
                ασ  ::Float64, 
                σσ  ::Float64, 
                αλ  ::Float64, 
                βλ  ::Float64, 
                σλ  ::Float64,
                llc ::Float64,
                dxs ::Float64,
                dxl ::Float64,
                ddx ::Float64,
                ddσ ::Float64,
                ssσ ::Float64,
                ddλ ::Float64,
                ssλ ::Float64,
                irλ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Perform xb tip updates.
"""
function update_tip!(tree::iTxb,
                     xavg::Float64, 
                     xstd::Float64,
                     ασ  ::Float64, 
                     σσ  ::Float64, 
                     αλ  ::Float64, 
                     βλ  ::Float64, 
                     σλ  ::Float64,
                     llc ::Float64,
                     dxs ::Float64,
                     dxl ::Float64,
                     ddx ::Float64,
                     ddσ ::Float64,
                     ssσ ::Float64,
                     ddλ ::Float64,
                     ssλ ::Float64,
                     irλ ::Float64,
                     δt  ::Float64,
                     srδt::Float64)
  @inbounds begin

    xc   = xv(tree)
    lσ2c = lσ2(tree)
    lλc  = lλ(tree)
    l    = lastindex(lλc)
    xic,   xfc = xc[1],   xc[l]
    lλic, lλfc = lλc[1], lλc[l]
    xp   = Vector{Float64}(undef, l)
    lλp  = Vector{Float64}(undef, l)
    lσ2p = Vector{Float64}(undef, l)
    ei   = e(tree)
    fdti = fdt(tree)

    # trait rate path sample
    bm!(lσ2p, lσ2c[1], ασ, σσ, δt, fdti, srδt)

    llσxr, llσσr, ssσr = llr_xb_σ(xc, ασ, σσ, lσ2p, lσ2c, δt, fdti)

    if -randexp() < llσxr
      llc += llσxr + llσσr
      ssσ += ssσr
      ddσ += lσ2p[l] - lσ2c[l] 
      unsafe_copyto!(lσ2c, 1, lσ2p, 1, l)
    end

    # trait and speciation path samples
    xfp = xfc
    if isnan(xavg)
      xfp  = rnorm(xic, intσ2(lσ2c, δt, fdti))
    elseif xstd > 0.0
      xfp = duoprop(xavg, xic, xstd^2, intσ2(lσ2c, δt, fdti))
    end

    lλfp = rnorm(lλic + αλ*ei + βλ*(xic - xfp), sqrt(ei)*σλ)

    cbb!(xp, xic, xfp, lσ2c, lλp, lλic, lλfp, βλ, σλ, δt, fdti, srδt)

    llbmr, llbr, dxsr, dxlr, ssλr, irλr = 
      llr_xb_b_sep(xp, xc, lσ2c, lλp, lλc, 
        ασ, σσ, αλ, βλ, σλ, δt, fdti, srδt, false)

    if -randexp() < llbr
      llc += llbmr + llbr
      dxs += dxsr
      dxl += dxlr
      ddx += xfp  - xfc
      ddλ += lλfp - lλfc
      ssλ += ssλr
      irλ += irλr
      unsafe_copyto!(xc,  1, xp,  1, l)
      unsafe_copyto!(lλc, 1, lλp, 1, l)
    end
  end

  return llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ
end




"""
    update_triad!(tree::iTxb,
                       ασ  ::Float64, 
                       σσ  ::Float64, 
                       αλ  ::Float64, 
                       βλ  ::Float64, 
                       σλ  ::Float64,
                       llc ::Float64,
                       dxs ::Float64,
                       dxl ::Float64,
                       ddx ::Float64,
                       ddσ ::Float64,
                       ssσ ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       irλ ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

Perform xb trio updates.
"""
function update_triad!(tree::iTxb,
                       ασ  ::Float64, 
                       σσ  ::Float64, 
                       αλ  ::Float64, 
                       βλ  ::Float64, 
                       σλ  ::Float64,
                       llc ::Float64,
                       dxs ::Float64,
                       dxl ::Float64,
                       ddx ::Float64,
                       ddσ ::Float64,
                       ssσ ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       irλ ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ =
    update_triad_xb!(tree, tree.d1, tree.d2, ασ, σσ, αλ, βλ, σλ, 
      llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt)

  return llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ
end




"""
    update_triad!(ξa  ::iTxb,
                  ξ1  ::iTxb,
                  ξ2  ::iTxb,
                  ασ  ::Float64, 
                  σσ  ::Float64, 
                  αλ  ::Float64, 
                  βλ  ::Float64, 
                  σλ  ::Float64,
                  llc ::Float64,
                  dxs ::Float64,
                  dxl ::Float64,
                  ddx ::Float64,
                  ddσ ::Float64,
                  ssσ ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64,
                  irλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Perform xb trio updates.
"""
function update_triad!(ξa  ::iTxb,
                       ξ1  ::iTxb,
                       ξ2  ::iTxb,
                       ασ  ::Float64, 
                       σσ  ::Float64, 
                       αλ  ::Float64, 
                       βλ  ::Float64, 
                       σλ  ::Float64,
                       llc ::Float64,
                       dxs ::Float64,
                       dxl ::Float64,
                       ddx ::Float64,
                       ddσ ::Float64,
                       ssσ ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       irλ ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  @inbounds begin
    xac,     x1c,   x2c =  xv(ξa),  xv(ξ1),  xv(ξ2)
    lσ2ac, lσ21c, lσ22c = lσ2(ξa), lσ2(ξ1), lσ2(ξ2)
    lλac,   lλ1c,  lλ2c =  lλ(ξa),  lλ(ξ1),  lλ(ξ2)
    la, l1, l2 = lastindex(lλac), lastindex(lλ1c), lastindex(lλ2c)
    xap,     x1p,   x2p = Vector{Float64}(undef,la), Vector{Float64}(undef,l1), Vector{Float64}(undef,l2)
    lσ2ap, lσ21p, lσ22p = Vector{Float64}(undef,la), Vector{Float64}(undef,l1), Vector{Float64}(undef,l2)
    lλap,   lλ1p,  lλ2p = Vector{Float64}(undef,la), Vector{Float64}(undef,l1), Vector{Float64}(undef,l2)
    ea, e1, e2, fdta, fdt1, fdt2 = e(ξa), e(ξ1), e(ξ2), fdt(ξa), fdt(ξ1), fdt(ξ2)

    lσ2ai, lσ21f, lσ22f = lσ2ac[1], lσ21c[l1], lσ22c[l2]
    xai,   x1f,   x2f   = xac[1],   x1c[l1],   x2c[l2]
    lλai,  lλ1f,  lλ2f  = lλac[1],  lλ1c[l1],  lλ2c[l2]

    # rate path sample
    lσ2n = trioprop(lσ2ai + ασ*ea, lσ21f - ασ*e1, lσ22f - ασ*e2, ea, e1, e2, σσ)

    bb!(lσ2ap, lσ2ai, lσ2n, σσ, δt, fdta, srδt)
    bb!(lσ21p, lσ2n, lσ21f, σσ, δt, fdt1, srδt)
    bb!(lσ22p, lσ2n, lσ22f, σσ, δt, fdt2, srδt)

    llσxar, llσσar, ssσar = llr_xb_σ(xac, ασ, σσ, lσ2ap, lσ2ac, δt, fdta)
    llσx1r, llσσ1r, ssσ1r = llr_xb_σ(x1c, ασ, σσ, lσ21p, lσ21c, δt, fdt1)
    llσx2r, llσσ2r, ssσ2r = llr_xb_σ(x2c, ασ, σσ, lσ22p, lσ22c, δt, fdt2)

    llr = llσxar + llσx1r + llσx2r

    if -randexp() < llr
      llc += llr + llσσar + llσσ1r + llσσ2r
      ssσ += ssσar + ssσ1r + ssσ2r
      ddσ += (lσ21c[1] - lσ2n)
      unsafe_copyto!(lσ2ac, 1, lσ2ap, 1, la)
      unsafe_copyto!(lσ21c, 1, lσ21p, 1, l1)
      unsafe_copyto!(lσ22c, 1, lσ22p, 1, l2)
    end

    # trait and speciation path samples
    xn  = trioprop(xai, x1f, x2f, intσ2(lσ2ac, δt, fdta), 
                                  intσ2(lσ21c, δt, fdt1), 
                                  intσ2(lσ22c, δt, fdt2))
    lλn = trioprop(lλai + αλ*ea + βλ*(xn - xai), 
                   lλ1f - αλ*e1 - βλ*(x1f - xn), 
                   lλ2f - αλ*e2 - βλ*(x2f - xn), 
                   σλ^2*ea, σλ^2*e1, σλ^2*e2)

    cbb!(xap, xai, xn, lσ2ac, lλap, lλai, lλn, βλ, σλ, δt, fdta, srδt)
    cbb!(x1p, xn, x1f, lσ21c, lλ1p, lλn, lλ1f, βλ, σλ, δt, fdt1, srδt)
    cbb!(x2p, xn, x2f, lσ22c, lλ2p, lλn, lλ2f, βλ, σλ, δt, fdt2, srδt)

    # likelihood ratio
    llbmar, llbar, dxsar, dxlar, ssλar, irλar = 
      llr_xb_b_sep(xap, xac, lσ2ac, lλap, lλac, 
        ασ, σσ, αλ, βλ, σλ, δt, fdta, srδt, true)
    llbm1r, llb1r, dxs1r, dxl1r, ssλ1r, irλ1r = 
      llr_xb_b_sep(x1p, x1c, lσ21c, lλ1p, lλ1c, 
        ασ, σσ, αλ, βλ, σλ, δt, fdt1, srδt, false)
    llbm2r, llb2r, dxs2r, dxl2r, ssλ2r, irλ2r = 
      llr_xb_b_sep(x2p, x2c, lσ22c, lλ2p, lλ2c, 
        ασ, σσ, αλ, βλ, σλ, δt, fdt2, srδt, false)

    llr = llbar + llb1r + llb2r

    if -randexp() < llr
      llc += llbmar + llbm1r + llbm2r + llr
      dxs += dxsar + dxs1r + dxs2r
      dxl += dxlar + dxl1r + dxl2r
      ddx += x1c[1]  -  xn
      ddλ += lλ1c[1] - lλn
      ssλ += ssλar + ssλ1r + ssλ2r
      irλ += irλar + irλ1r + irλ2r
      unsafe_copyto!(xac,  1, xap,  1, la)
      unsafe_copyto!(x1c,  1, x1p,  1, l1)
      unsafe_copyto!(x2c,  1, x2p,  1, l2)
      unsafe_copyto!(lλac, 1, lλap, 1, la)
      unsafe_copyto!(lλ1c, 1, lλ1p, 1, l1)
      unsafe_copyto!(lλ2c, 1, lλ2p, 1, l2)
    end
  end

  return llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ
end



