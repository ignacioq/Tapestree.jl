#=

Diffused Brownian motion proposals

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#





"""
    _fstem_update!(ξi       ::sTxs,
                   ξ1       ::sTxs,
                   αx       ::Float64,
                   ασ       ::Float64,
                   γ        ::Float64,
                   δt       ::Float64,
                   srδt     ::Float64)

Do dbm update for fossil stem root.
"""
function _fstem_update!(ξi       ::sTxs,
                        ξ1       ::sTxs,
                        αx       ::Float64,
                        ασ       ::Float64,
                        γ        ::Float64,
                        δt       ::Float64,
                        srδt     ::Float64)

  @inbounds begin
    lσ21 = lσ2(ξ1)
    x1   = xv(ξ1)
    l    = lastindex(lσ21)
    lσ2p = Vector{Float64}(undef,l)
    lσ2n = lσ21[l]
    xn   = x1[l]
    el   = e(ξ1)
    fdt1 = fdt(ξ1)

    # rate path sample
    lσ2r = rnorm(lσ2n - ασ*el, γ*sqrt(el))
    bb!(lσ2p, lσ2r, lσ2n, γ, δt, fdt1, srδt)

    llr = llr_dbm_σ(x1, αx, lσ2p, lσ21, δt, fdt1)

    if -randexp() < llr
      unsafe_copyto!(lσ21, 1, lσ2p, 1, l)
      fill!(lσ2(ξi), lσ2r)
    end

    # trait path sample
    xr = rnorm(xn - αx*el, sqrt(intσ2(lσ21, δt, fdt1)))
    dbb!(x1, xr, xn, lσ21, δt, fdt1, srδt)
    fill!(xv(ξi), xr)

    # likelihood
    ll, Ls, Xs, dd, ss = ll_dbm_dd_ss_b(x1, αx, lσ21, ασ, γ, δt, fdt1)
  end

  return ll, Ls, Xs, dd, ss
end




"""
    _stem_update!(ξi       ::sTxs,
                  αx       ::Float64,
                  ασ       ::Float64,
                  γ        ::Float64,
                  δt       ::Float64,
                  srδt     ::Float64)

Do dbm update for stem root.
"""
function _stem_update!(ξi       ::sTxs,
                       αx       ::Float64,
                       ασ       ::Float64,
                       γ        ::Float64,
                       δt       ::Float64,
                       srδt     ::Float64)

  @inbounds begin
    lσ2c = lσ2(ξi)
    xc   = xv(ξi)
    l    = lastindex(lσ2c)
    lσ2p = Vector{Float64}(undef,l)
    lσ2n = lσ2c[l]
    xn   = xc[l]
    el   = e(ξi)
    fdtp = fdt(ξi)

    # rate path sample
    lσ2r = rnorm(lσ2n - ασ*el, γ*sqrt(el))
    bb!(lσ2p, lσ2r, lσ2n, γ, δt, fdtp, srδt)

    llr = llr_dbm_σ(xc, αx, lσ2p, lσ2c, δt, fdtp)

    if -randexp() < llr
      unsafe_copyto!(lσ2c, 1, lσ2p, 1, l)
    end

    # trait path sample
    xr = rnorm(xn - αx*el, sqrt(intσ2(lσ2c, δt, fdtp)))
    dbb!(xc, xr, xn, lσ2c, δt, fdtp, srδt)

    # likelihood
    ll, Ls, Xs, dd, ss = ll_dbm_dd_ss_b(xc, αx, lσ2c, ασ, γ, δt, fdtp)
  end

  return ll, Ls, Xs, dd, ss
end




"""
    _crown_update!(ξi       ::sTxs,
                   ξ1       ::sTxs,
                   ξ2       ::sTxs,
                   αx       ::Float64,
                   ασ       ::Float64,
                   γ        ::Float64,
                   δt       ::Float64,
                   srδt     ::Float64)

Do dbm update for crown root.
"""
function _crown_update!(ξi       ::sTxs,
                        ξ1       ::sTxs,
                        ξ2       ::sTxs,
                        αx       ::Float64,
                        ασ       ::Float64,
                        γ        ::Float64,
                        δt       ::Float64,
                        srδt     ::Float64)
  @inbounds begin
    lσ2c  = lσ2(ξi)
    lσ21  = lσ2(ξ1)
    lσ22  = lσ2(ξ2)
    xc    = xv(ξi)
    x1    = xv(ξ1)
    x2    = xv(ξ2)
    l1    = lastindex(x1)
    l2    = lastindex(x2)
    lσ21p = Vector{Float64}(undef,l1)
    lσ22p = Vector{Float64}(undef,l2)
    lσ21f = lσ21[l1]
    lσ22f = lσ22[l2]
    x1f   = x1[l1]
    x2f   = x2[l2]
    e1    = e(ξ1)
    e2    = e(ξ2)
    fdt1  = fdt(ξ1)
    fdt2  = fdt(ξ2)

    # rate path sample
    lσ2n = duoprop(lσ21f - ασ*e1, lσ22f - ασ*e2, γ^2*e1, γ^2*e2)
    bb!(lσ21p, lσ2n, lσ21f, γ, δt, fdt1, srδt)
    bb!(lσ22p, lσ2n, lσ22f, γ, δt, fdt2, srδt)

    llr = llr_dbm_σ(x1, αx, lσ21p, lσ21, δt, fdt1) +
          llr_dbm_σ(x2, αx, lσ22p, lσ22, δt, fdt2)

    if -randexp() < llr
      unsafe_copyto!(lσ21, 1, lσ21p, 1, l1)
      unsafe_copyto!(lσ22, 1, lσ22p, 1, l2)
      fill!(lσ2c, lσ2n)
    end

    # trait path sample
    xn = duoprop(x1f - αx*e1, x2f - αx*e2, 
           intσ2(lσ21, δt, fdt1), intσ2(lσ22, δt, fdt2))
    dbb!(x1, xn, x1f, lσ21, δt, fdt1, srδt)
    dbb!(x2, xn, x2f, lσ22, δt, fdt2, srδt)

    # fill root
    fill!(xc, xn)

    # log likelihood ratios
    ll1, Ls1, Xs1, dd1, ss1 = ll_dbm_dd_ss_b(x1, αx, lσ21, ασ, γ, δt, fdt1)
    ll2, Ls2, Xs2, dd2, ss2 = ll_dbm_dd_ss_b(x2, αx, lσ22, ασ, γ, δt, fdt2)
  end

  return ll1, ll2, Ls1, Ls2, Xs1, Xs2, dd1, dd2, ss1, ss2
end




"""
    _update_leaf_x!(ξi  ::sTxs,
                    xavg::Float64,
                    xstd::Float64,
                    αx  ::Float64,
                    ασ  ::Float64,
                    γ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

Make a `dbm` **fixed** tip proposal.
"""
function _update_leaf_x!(ξi  ::sTxs,
                         xavg::Float64,
                         xstd::Float64,
                         αx  ::Float64,
                         ασ  ::Float64,
                         γ   ::Float64,
                         δt  ::Float64,
                         srδt::Float64)
  lσ2c = lσ2(ξi)
  xc   = xv(ξi)
  l    = lastindex(lσ2c)
  lσ2p = Vector{Float64}(undef,l)
  xi   = xc[1]
  xn   = xc[l]
  el   = e(ξi)
  fdtp = fdt(ξi)
  
  # rate path sample
  bm!(lσ2p, lσ2c[1], ασ, γ, δt, fdtp, srδt)

  llr = llr_dbm_σ(xc, αx, lσ2p, lσ2c, δt, fdtp)

  if -randexp() < llr
    unsafe_copyto!(lσ2c, 1, lσ2p, 1, l)
  end

  # trait path sample
  if !iszero(xstd)
    xn = duoprop(xavg, xi + αx*el, xstd^2, intσ2(lσ2c, δt, fdtp))
  end
  dbb!(xc, xi, xn, lσ2c, δt, fdtp, srδt)

  # likelihood
  ll, Ls, Xs, dd, ss = ll_dbm_dd_ss_b(xc, αx, lσ2c, ασ, γ, δt, fdtp)

  return ll, Ls, Xs, dd, ss
end




"""
    _update_leaf_x!(ξi  ::sTxs,
                    αx  ::Float64,
                    ασ  ::Float64,
                    γ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

Make a `dbm` **unfixed** tip proposal.
"""
function _update_leaf_x!(ξi  ::sTxs,
                         αx  ::Float64,
                         ασ  ::Float64,
                         γ   ::Float64,
                         δt  ::Float64,
                         srδt::Float64)

  lσ2c = lσ2(ξi)
  xc   = xv(ξi)
  fdtp = fdt(ξi)

  # trait and rate path sample
  dbm!(xc, xc[1], αx, lσ2c, lσ2c[1], ασ, γ, δt, fdtp, srδt)

  # likelihood
  ll, Ls, Xs, dd, ss = ll_dbm_dd_ss_b(xc, αx, lσ2c, ασ, γ, δt, fdtp)

  return ll, Ls, Xs, dd, ss
end




"""
    _update_duo_x!(ξi  ::sTxs,
                   ξ1  ::sTxs,
                   xavg::Float64,
                   xstd::Float64,
                   αx  ::Float64,
                   ασ  ::Float64,
                   γ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

Do duo `dbm` update for **fixed** node.
"""
function _update_duo_x!(ξi  ::sTxs,
                        ξ1  ::sTxs,
                        xavg::Float64,
                        xstd::Float64,
                        αx  ::Float64,
                        ασ  ::Float64,
                        γ   ::Float64,
                        δt  ::Float64,
                        srδt::Float64)

  @inbounds begin
    lσ2a  = lσ2(ξi)
    lσ21  = lσ2(ξ1)
    xa    = xv(ξi)
    x1    = xv(ξ1)
    la    = lastindex(xa)
    l1    = lastindex(x1)
    lσ2ap = Vector{Float64}(undef,la)
    lσ21p = Vector{Float64}(undef,l1)
    lσ2ai = lσ2a[1]
    lσ21f = lσ21[l1]
    xai   = xa[1]
    xn    = xa[la]
    x1f   = x1[l1]
    ei    = e(ξi)
    e1    = e(ξ1)
    fdta  = fdt(ξi)
    fdt1  = fdt(ξ1)

    # rate path sample
    lσ2n = duoprop(lσ2ai + ασ*ei, lσ21f - ασ*e1, ei, e1, γ)
    bb!(lσ2ap, lσ2ai, lσ2n, γ, δt, fdta, srδt)
    bb!(lσ21p, lσ2n, lσ21f, γ, δt, fdt1, srδt)

    llr = llr_dbm_σ(xa, αx, lσ2ap, lσ2a, δt, fdta) + 
          llr_dbm_σ(x1, αx, lσ21p, lσ21, δt, fdt1)

    if -randexp() < llr
      unsafe_copyto!(lσ2a, 1, lσ2ap, 1, la)
      unsafe_copyto!(lσ21, 1, lσ21p, 1, l1)
    end

    # trait path sample
    if !iszero(xstd)
      xn = trioprop(xavg, xai + αx*ei, x1f - αx*e1, 
             xstd^2, intσ2(lσ2a, δt, fdta), intσ2(lσ21, δt, fdt1))
    end
    dbb!(xa, xai, xn, lσ2a, δt, fdta, srδt)
    dbb!(x1, xn, x1f, lσ21, δt, fdt1, srδt)

    # log likelihood ratios
    lla, Lsa, Xsa, dda, ssa = ll_dbm_dd_ss_b(xa, αx, lσ2a, ασ, γ, δt, fdta)
    ll1, Ls1, Xs1, dd1, ss1 = ll_dbm_dd_ss_b(x1, αx, lσ21, ασ, γ, δt, fdt1)
  end

  return lla, ll1, Lsa, Ls1, Xsa, Xs1, dda, dd1, ssa, ss1
end




"""
    _update_duo_x!(ξi  ::sTxs,
                   ξ1  ::sTxs,
                   αx  ::Float64,
                   ασ  ::Float64,
                   γ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

Do duo `dbm` update for **unfixed** node.
"""
function _update_duo_x!(ξi  ::sTxs,
                        ξ1  ::sTxs,
                        αx  ::Float64,
                        ασ  ::Float64,
                        γ   ::Float64,
                        δt  ::Float64,
                        srδt::Float64)

  @inbounds begin
    lσ2a  = lσ2(ξi)
    lσ21  = lσ2(ξ1)
    xa    = xv(ξi)
    x1    = xv(ξ1)
    la    = lastindex(xa)
    l1    = lastindex(x1)
    lσ2ap = Vector{Float64}(undef,la)
    lσ21p = Vector{Float64}(undef,l1)
    lσ2ai = lσ2a[1]
    lσ21f = lσ21[l1]
    xai   = xa[1]
    x1f   = x1[l1]
    ei    = e(ξi)
    e1    = e(ξ1)
    fdta  = fdt(ξi)
    fdt1  = fdt(ξ1)

    # rate path sample
    lσ2n = duoprop(lσ2ai + ασ*ei, lσ21f - ασ*e1, ei, e1, γ)
    bb!(lσ2ap, lσ2ai, lσ2n, γ, δt, fdta, srδt)
    bb!(lσ21p, lσ2n, lσ21f, γ, δt, fdt1, srδt)

    llr = llr_dbm_σ(xa, αx, lσ2ap, lσ2a, δt, fdta) + 
          llr_dbm_σ(x1, αx, lσ21p, lσ21, δt, fdt1)

    if -randexp() < llr
      unsafe_copyto!(lσ2a, 1, lσ2ap, 1, la)
      unsafe_copyto!(lσ21, 1, lσ21p, 1, l1)
    end

    # trait path sample
    xn = duoprop(xai + αx*ei, x1f - αx*e1, 
          intσ2(lσ2a, δt, fdta), intσ2(lσ21, δt, fdt1))
    dbb!(xa, xai, xn, lσ2a, δt, fdta, srδt)
    dbb!(x1, xn, x1f, lσ21, δt, fdt1, srδt)

    # log likelihood ratios
    lla, Lsa, Xsa, dda, ssa = ll_dbm_dd_ss_b(xa, αx, lσ2a, ασ, γ, δt, fdta)
    ll1, Ls1, Xs1, dd1, ss1 = ll_dbm_dd_ss_b(x1, αx, lσ21, ασ, γ, δt, fdt1)
  end

  return lla, ll1, Lsa, Ls1, Xsa, Xs1, dda, dd1, ssa, ss1
end




"""
    _update_triad_x!(ξi   ::sTxs,
                     ξ1   ::sTxs,
                     ξ2   ::sTxs,
                     αx  ::Float64,
                     ασ  ::Float64,
                     γ    ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)

Make a `gbm` trio proposal.
"""
function _update_triad_x!(ξi   ::sTxs,
                          ξ1   ::sTxs,
                          ξ2   ::sTxs,
                          αx  ::Float64,
                          ασ  ::Float64,
                          γ    ::Float64,
                          δt   ::Float64,
                          srδt ::Float64)

  @inbounds begin
    lσ2a  = lσ2(ξi)
    lσ21  = lσ2(ξ1)
    lσ22  = lσ2(ξ2)
    xa    = xv(ξi)
    x1    = xv(ξ1)
    x2    = xv(ξ2)
    la    = lastindex(xa)
    l1    = lastindex(x1)
    l2    = lastindex(x2)
    lσ2ap = Vector{Float64}(undef,la)
    lσ21p = Vector{Float64}(undef,l1)
    lσ22p = Vector{Float64}(undef,l2)
    lσ2ai = lσ2a[1]
    lσ21f = lσ21[l1]
    lσ22f = lσ22[l2]
    xai   = xa[1]
    x1f   = x1[l1]
    x2f   = x2[l2]
    ei    = e(ξi) 
    e1    = e(ξ1) 
    e2    = e(ξ2)
    fdta  = fdt(ξi)
    fdt1  = fdt(ξ1)
    fdt2  = fdt(ξ2)

    # rate path sample
    lσ2n = trioprop(lσ2ai + ασ*ei, lσ21f - ασ*e1, lσ22f - ασ*e2, ei, e1, e2, γ)
    bb!(lσ2ap, lσ2ai, lσ2n, γ, δt, fdta, srδt)
    bb!(lσ21p, lσ2n, lσ21f, γ, δt, fdt1, srδt)
    bb!(lσ22p, lσ2n, lσ22f, γ, δt, fdt2, srδt)

    llr = llr_dbm_σ(xa, αx, lσ2ap, lσ2a, δt, fdta) +
          llr_dbm_σ(x1, αx, lσ21p, lσ21, δt, fdt1) +
          llr_dbm_σ(x2, αx, lσ22p, lσ22, δt, fdt2)

    if -randexp() < llr
      unsafe_copyto!(lσ2a, 1, lσ2ap, 1, la)
      unsafe_copyto!(lσ21, 1, lσ21p, 1, l1)
      unsafe_copyto!(lσ22, 1, lσ22p, 1, l2)
    end

    # trait path sample
    xn = trioprop(xai + αx*ei, x1f - αx*e1, x2f - αx*e2, 
           intσ2(lσ2a, δt, fdta), intσ2(lσ21, δt, fdt1), intσ2(lσ22, δt, fdt2))
    dbb!(xa, xai, xn, lσ2a, δt, fdta, srδt)
    dbb!(x1, xn, x1f, lσ21, δt, fdt1, srδt)
    dbb!(x2, xn, x2f, lσ22, δt, fdt2, srδt)

    # log likelihood ratios
    lla, Lsa, Xsa, dda, ssa = ll_dbm_dd_ss_b(xa, αx, lσ2a, ασ, γ, δt, fdta)
    ll1, Ls1, Xs1, dd1, ss1 = ll_dbm_dd_ss_b(x1, αx, lσ21, ασ, γ, δt, fdt1)
    ll2, Ls2, Xs2, dd2, ss2 = ll_dbm_dd_ss_b(x2, αx, lσ22, ασ, γ, δt, fdt2)
  end

  return lla, ll1, ll2, Lsa, Ls1, Ls2, Xsa, Xs1, Xs2, 
         dda, dd1, dd2, ssa, ss1, ss2
end




