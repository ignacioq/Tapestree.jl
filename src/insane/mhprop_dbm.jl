#=

Diffused Brownian motion proposals

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#




"""
    _stem_update!(ξi   ::sTxs,
                  γ    ::Float64,
                  δt   ::Float64,
                  srδt ::Float64)

Do dbm update for stem root.
"""
function _stem_update!(ξi   ::sTxs,
                       γ    ::Float64,
                       δt   ::Float64,
                       srδt ::Float64)
  @inbounds begin
    σc   = lσ(ξi)
    xc   = xv(ξi)
    l    = lastindex(σc)
    σp   = Vector{Float64}(undef,l)
    σn   = σc[l]
    xn   = xc[l]
    fdtp = fdt(ξi)

    # rate path sample
    σr = rnorm(σn, γ*sqrt(e(ξi)))
    bb!(σp, σr, σn, γ, δt, fdtp, srδt)

    llr = llr_dbm(xc, σp, σc, δt, fdtp)

    if -randexp() < llr
      unsafe_copyto!(σc, 1, σp, 1, l)
    end

    # trait path sample
    xr = rnorm(xn, intσ(σc, δt, fdtp))
    dbb!(xc, xr, xn, σc, δt, fdtp, srδt)

    # likelihood
    ll, ss = ll_dbm_ss_b(xc, σc, γ, δt, fdtp)
  end

  return ll, ss
end




"""
    _crown_update!(ξi   ::sTxs,
                   ξ1   ::sTxs,
                   ξ2   ::sTxs,
                   γ    ::Float64,
                   δt   ::Float64,
                   srδt ::Float64)

Do dbm update for crown root.
"""
function _crown_update!(ξi   ::sTxs,
                        ξ1   ::sTxs,
                        ξ2   ::sTxs,
                        γ    ::Float64,
                        δt   ::Float64,
                        srδt ::Float64)

  @inbounds begin
    σc   = lσ(ξi)
    σ1   = lσ(ξ1)
    σ2   = lσ(ξ2)
    xc   = xv(ξi)
    x1   = xv(ξ1)
    x2   = xv(ξ2)
    l1   = lastindex(x1)
    l2   = lastindex(x2)
    σ1p  = Vector{Float64}(undef,l1)
    σ2p  = Vector{Float64}(undef,l2)
    σ1f  = σ1[l1]
    σ2f  = σ2[l2]
    x1f  = x1[l1]
    x2f  = x2[l2]
    fdt1 = fdt(ξ1)
    fdt2 = fdt(ξ2)

    # rate path sample
    σn = duoprop(σ1f, σ2f,  e(ξ1), e(ξ2), γ)
    bb!(σ1p, σn, σ1f, γ, δt, fdt1, srδt)
    bb!(σ2p, σn, σ2f, γ, δt, fdt2, srδt)

    llr = llr_dbm(x1, σ1p, σ1, δt, fdt1) +
          llr_dbm(x2, σ2p, σ2, δt, fdt2)

    if -randexp() < llr
      unsafe_copyto!(σ1, 1, σ1p, 1, l1)
      unsafe_copyto!(σ2, 1, σ2p, 1, l2)
      fill!(σc, σn)
    end

    # trait path sample
    xn = duoprop(x1f, x2f, intσ(σ1, δt, fdt1), intσ(σ2, δt, fdt2))
    dbb!(x1, xn, x1f, σ1, δt, fdt1, srδt)
    dbb!(x2, xn, x2f, σ2, δt, fdt2, srδt)

    # fill root
    fill!(xc, xn)

    # log likelihood ratios
    ll1, ss1 = ll_dbm_ss_b(x1, σ1, γ, δt, fdt1)
    ll2, ss2 = ll_dbm_ss_b(x2, σ2, γ, δt, fdt2)
  end

  return ll1, ll2, ss1, ss2
end




"""
    _update_leaf_x!(ξi  ::sTxs,
                    xavg::Float64,
                    xstd::Float64,
                    γ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

Make a `dbm` **fixed** tip proposal.
"""
function _update_leaf_x!(ξi  ::sTxs,
                         xavg::Float64,
                         xstd::Float64,
                         γ   ::Float64,
                         δt  ::Float64,
                         srδt::Float64)
  σc   = lσ(ξi)
  xc   = xv(ξi)
  l    = lastindex(σc)
  σp   = Vector{Float64}(undef,l)
  xn   = xc[l]
  fdtp = fdt(ξi)

  # rate path sample
  bm!(σp, σc[1], γ, δt, fdtp, srδt)

  llr = llr_dbm(xc, σp, σc, δt, fdtp)

  if -randexp() < llr
    unsafe_copyto!(σc, 1, σp, 1, l)
  end

  # trait path sample
  xi = xc[1]
  if !iszero(xstd)
    xp = rnorm(xi, intσ(σc, δt, fdtp))
    if -randexp() < llrdnorm_x(xp, xn, xavg, xstd^2)
      xn = xp
    end
  end

  dbb!(xc, xi, xn, σc, δt, fdtp, srδt)

  # likelihood
  ll, ss = ll_dbm_ss_b(xc, σc, γ, δt, fdtp)

  return ll, ss
end




"""
    _update_leaf_x!(ξi  ::sTxs,
                    γ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

Make a `dbm` **unfixed** tip proposal.
"""
function _update_leaf_x!(ξi  ::sTxs,
                         γ   ::Float64,
                         δt  ::Float64,
                         srδt::Float64)

  σc   = lσ(ξi)
  xc   = xv(ξi)
  fdtp = fdt(ξi)

  # trait and rate path sample
  dbm!(xc, xc[1], σc, σc[1], γ, fdtp, srδt)

  # likelihood
  ll, ss = ll_dbm_ss_b(xc, σc, γ, δt, fdtp)

  return ll, ss
end




"""
    _update_duo_x!(ξi  ::sTxs,
                   ξ1  ::sTxs,
                   xavg::Float64,
                   xstd::Float64,
                   γ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

Do duo `dbm` update for **fixed** node.
"""
function _update_duo_x!(ξi  ::sTxs,
                        ξ1  ::sTxs,
                        xavg::Float64,
                        xstd::Float64,
                        γ   ::Float64,
                        δt  ::Float64,
                        srδt::Float64)

  @inbounds begin
    σa   = lσ(ξi)
    σ1   = lσ(ξ1)
    xa   = xv(ξi)
    x1   = xv(ξ1)
    la   = lastindex(xa)
    l1   = lastindex(x1)
    σap  = Vector{Float64}(undef,la)
    σ1p  = Vector{Float64}(undef,l1)
    σai  = σa[1]
    σ1f  = σ1[l1]
    xai  = xa[1]
    xn   = xa[la]
    x1f  = x1[l1]
    fdta = fdt(ξi)
    fdt1 = fdt(ξ1)

    # rate path sample
    σn = duoprop(σai, σ1f, e(ξi), e(ξ1), γ)
    bb!(σap, σai, σn, γ, δt, fdta, srδt)
    bb!(σ1p, σn, σ1f, γ, δt, fdt1, srδt)

    llr = llr_dbm(xa, σap, σa, δt, fdta) + 
          llr_dbm(x1, σ1p, σ1, δt, fdt1)

    if -randexp() < llr
      unsafe_copyto!(σa, 1, σap, 1, la)
      unsafe_copyto!(σ1, 1, σ1p, 1, l1)
    end

    # trait path sample
    if !iszero(xstd)
      xp = duoprop(xai, x1f, intσ(σa, δt, fdta), intσ(σ1, δt, fdt1))
      if -randexp() < llrdnorm_x(xp, xn, xavg, xstd^2)
        xn = xp
      end
    end
    dbb!(xa, xai, xn, σa, δt, fdta, srδt)
    dbb!(x1, xn, x1f, σ1, δt, fdt1, srδt)

    # log likelihood ratios
    lla, ssa = ll_dbm_ss_b(xa, σa, γ, δt, fdta)
    ll1, ss1 = ll_dbm_ss_b(x1, σ1, γ, δt, fdt1)
  end

  return lla, ll1, ssa, ss1
end





"""
    _update_duo_x!(ξi   ::sTxs,
                   ξ1   ::sTxs,
                   γ    ::Float64,
                   δt   ::Float64,
                   srδt ::Float64)

Do duo `dbm` update for **unfixed** node.
"""
function _update_duo_x!(ξi   ::sTxs,
                        ξ1   ::sTxs,
                        γ    ::Float64,
                        δt   ::Float64,
                        srδt ::Float64)

  @inbounds begin
    σa   = lσ(ξi)
    σ1   = lσ(ξ1)
    xa   = xv(ξi)
    x1   = xv(ξ1)
    la   = lastindex(xa)
    l1   = lastindex(x1)
    σap  = Vector{Float64}(undef,la)
    σ1p  = Vector{Float64}(undef,l1)
    σai  = σa[1]
    σ1f  = σ1[l1]
    xai  = xa[1]
    x1f  = x1[l1]
    fdta = fdt(ξi)
    fdt1 = fdt(ξ1)

    # rate path sample
    σn = duoprop(σai, σ1f, e(ξi), e(ξ1), γ)
    bb!(σap, σai, σn, γ, δt, fdta, srδt)
    bb!(σ1p, σn, σ1f, γ, δt, fdt1, srδt)

    llr = llr_dbm(xa, σap, σa, δt, fdta) + 
          llr_dbm(x1, σ1p, σ1, δt, fdt1)

    if -randexp() < llr
      unsafe_copyto!(σa, 1, σap, 1, la)
      unsafe_copyto!(σ1, 1, σ1p, 1, l1)
    end

    # trait path sample
    xn = duoprop(xai, x1f, intσ(σa, δt, fdta), intσ(σ1, δt, fdt1))
    dbb!(xa, xai, xn, σa, δt, fdta, srδt)
    dbb!(x1, xn, x1f, σ1, δt, fdt1, srδt)

    # log likelihood ratios
    lla, ssa = ll_dbm_ss_b(xa, σa, γ, δt, fdta)
    ll1, ss1 = ll_dbm_ss_b(x1, σ1, γ, δt, fdt1)
  end

  return lla, ll1, ssa, ss1
end




"""
    _update_triad_x!(ξi   ::sTxs,
                     ξ1   ::sTxs,
                     ξ2   ::sTxs,
                     γ    ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)

Make a `gbm` trio proposal.
"""
function _update_triad_x!(ξi   ::sTxs,
                          ξ1   ::sTxs,
                          ξ2   ::sTxs,
                          γ    ::Float64,
                          δt   ::Float64,
                          srδt ::Float64)

  @inbounds begin
    σa   = lσ(ξi)
    σ1   = lσ(ξ1)
    σ2   = lσ(ξ2)
    xa   = xv(ξi)
    x1   = xv(ξ1)
    x2   = xv(ξ2)
    la   = lastindex(xa)
    l1   = lastindex(x1)
    l2   = lastindex(x2)
    σap  = Vector{Float64}(undef,la)
    σ1p  = Vector{Float64}(undef,l1)
    σ2p  = Vector{Float64}(undef,l2)
    σai  = σa[1]
    σ1f  = σ1[l1]
    σ2f  = σ2[l2]
    xai  = xa[1]
    x1f  = x1[l1]
    x2f  = x2[l2]
    fdta = fdt(ξi)
    fdt1 = fdt(ξ1)
    fdt2 = fdt(ξ2)

    # rate path sample
    σn = trioprop(σai, σ1f, σ2f, e(ξi), e(ξ1), e(ξ2), γ)
    bb!(σap, σai, σn, γ, δt, fdta, srδt)
    bb!(σ1p, σn, σ1f, γ, δt, fdt1, srδt)
    bb!(σ2p, σn, σ2f, γ, δt, fdt2, srδt)

    llr = llr_dbm(xa, σap, σa, δt, fdta) +
          llr_dbm(x1, σ1p, σ1, δt, fdt1) +
          llr_dbm(x2, σ2p, σ2, δt, fdt2)

    if -randexp() < llr
      unsafe_copyto!(σa, 1, σap, 1, la)
      unsafe_copyto!(σ1, 1, σ1p, 1, l1)
      unsafe_copyto!(σ2, 1, σ2p, 1, l2)
    end

    # trait path sample
    xn = trioprop(xai, x1f, x2f, 
           intσ(σa, δt, fdta), intσ(σ1, δt, fdt1), intσ(σ2, δt, fdt2))
    dbb!(xa, xai, xn, σa, δt, fdta, srδt)
    dbb!(x1, xn, x1f, σ1, δt, fdt1, srδt)
    dbb!(x2, xn, x2f, σ2, δt, fdt2, srδt)

    # log likelihood ratios
    lla, ssc = ll_dbm_ss_b(xa, σa, γ, δt, fdta)
    ll1, ss1 = ll_dbm_ss_b(x1, σ1, γ, δt, fdt1)
    ll2, ss2 = ll_dbm_ss_b(x2, σ2, γ, δt, fdt2)
  end

  return lla, ll1, ll2, ssc, ss1, ss2
end
