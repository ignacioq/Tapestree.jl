#=

Anagenetic GBM pure-birth MCMC MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 14 11 2021
=#




"""
    _daughters_update!(ξ1  ::iTbX,
                       ξ2  ::iTbX,
                       λf  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       xp  ::Float64,
                       σx  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

Make a `gbmbX` proposal for daughters from forwards simulated branch.
"""
function _daughters_update!(ξ1  ::iTbX,
                            ξ2  ::iTbX,
                            λf  ::Float64,
                            α   ::Float64,
                            σλ  ::Float64,
                            xp  ::Float64,
                            βλ  ::Float64,
                            σx  ::Float64,
                            δt  ::Float64,
                            srδt::Float64)
  @inbounds begin

    λ1c  = lλ(ξ1)
    λ2c  = lλ(ξ2)
    x1c  = xv(ξ1)
    x2c  = xv(ξ2)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    x1p  = Vector{Float64}(undef,l1)
    x2p  = Vector{Float64}(undef,l2)
    λi   = λ1c[1]
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]
    e1   = e(ξ1)
    e2   = e(ξ2)
    fdt1 = fdt(ξ1)
    fdt2 = fdt(ξ2)

    bb!(λ1p, λf, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λf, λ2, σλ, δt, fdt2, srδt)
    bb!(x1p, xp, x1c[l1], σx, δt, fdt1, srδt)
    bb!(x2p, xp, x2c[l2], σx, δt, fdt2, srδt)

    # log likelihood ratios (llrbm contains both λ and x llr)
    llrbm1, llrb1, ssrλ1, ssrx1 =
      llr_gbm_b_sep(λ1p, λ1c, x1p, x1c, α, σλ, βλ, σx, δt, fdt1, srδt, false)
    llrbm2, llrb2, ssrλ2, ssrx2 =
      llr_gbm_b_sep(λ2p, λ2c, x2p, x2c, α, σλ, βλ, σx, δt, fdt2, srδt, false)

    acr  = llrb1 + llrb2
    llr  = llrbm1 + llrbm2 + acr  + λf - λi
    drλ  = 2.0*(λi - λf)
    ssrλ = ssrλ1 + ssrλ2
    ssrx = ssrx1 + ssrx2
  end

  return llr, acr, drλ, ssrλ, ssrx, λ1p, λ2p, x1p, x2p
end




"""
    _stem_update!(ξi   ::iTbX,
                  α    ::Float64,
                  σλ   ::Float64,
                  llc  ::Float64,
                  dλ   ::Float64,
                  ssλ  ::Float64,
                  δt   ::Float64,
                  srδt ::Float64)

Do gbm update for crown root.
"""
function _stem_update!(ξi   ::iTbX,
                       α    ::Float64,
                       σλ   ::Float64,
                       llc  ::Float64,
                       dλ   ::Float64,
                       ssλ  ::Float64,
                       δt   ::Float64,
                       srδt ::Float64)

  @inbounds begin
    λc   = lλ(ξi)
    l    = lastindex(λc)
    λp   = Vector{Float64}(undef,l)
    λn   = λc[l]
    el   = e(ξi)
    fdtp = fdt(ξi)

    # node proposal
    λr = rnorm(λn - α*el, σλ*sqrt(el))

    # simulate fix tree vector
    bb!(λp, λr, λn, σλ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ = llr_gbm_b_sep(λp, λc, α, σλ, δt, fdtp, srδt, false)

    acr = llrbd

    if -randexp() < acr
      llc += acr + llrbm
      dλ  += λc[1] - λr
      ssλ += ssrλ
      unsafe_copyto!(λc, 1, λp, 1, l)
    end
  end

  return llc, dλ, ssλ
end




"""
    _crown_update!(ξi   ::iTbX,
                   ξ1   ::iTbX,
                   ξ2   ::iTbX,
                   α    ::Float64,
                   σλ   ::Float64,
                   llc  ::Float64,
                   dλ   ::Float64,
                   ssλ  ::Float64,
                   δt   ::Float64,
                   srδt ::Float64)

Do gbm update for crown root.
"""
function _crown_update!(ξi   ::iTbX,
                        ξ1   ::iTbX,
                        ξ2   ::iTbX,
                        α    ::Float64,
                        σλ   ::Float64,
                        llc  ::Float64,
                        dλ   ::Float64,
                        ssλ  ::Float64,
                        δt   ::Float64,
                        srδt ::Float64)

  @inbounds begin
    λpc  = lλ(ξi)
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
    λr = duoprop(λ1 - α*e1, λ2 - α*e2, e1, e2, σλ)

    # simulate fix tree vector
    bb!(λ1p, λr, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λr, λ2, σλ, δt, fdt2, srδt)

    # log likelihood ratios
    llrbm1, llrb1, ssrλ1 =
      llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
    llrbm2, llrb2, ssrλ2 =
      llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

    acr  = llrb1 + llrb2

    if -randexp() < acr
      llc += llrbm1 + llrbm2 + acr
      dλ  += 2.0*(λ1c[1] - λr)
      ssλ += ssrλ1 + ssrλ2
      fill!(λpc, λr)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, dλ, ssλ
end




"""
    _update_gbm!(tree::iTbX,
                 α   ::Float64,
                 σλ  ::Float64,
                 llc ::Float64,
                 dλ  ::Float64,
                 ssλ ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 ter ::Bool)

Do gbm updates on a decoupled tree recursively.
"""
function _update_gbm!(tree::iTbX,
                      α   ::Float64,
                      σλ  ::Float64,
                      llc ::Float64,
                      dλ  ::Float64,
                      ssλ ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ter ::Bool)

  if def1(tree)
    llc, dλ, ssλ = update_triad!(tree, α, σλ, llc, dλ, ssλ, δt, srδt)

    llc, dλ, ssλ =
      _update_gbm!(tree.d1, α, σλ, llc, dλ, ssλ, δt, srδt, ter)
    llc, dλ, ssλ =
      _update_gbm!(tree.d2, α, σλ, llc, dλ, ssλ, δt, srδt, ter)
  else
    # if !isfix(tree) || ter
    #   llc, dλ, ssλ = update_tip!(tree, α, σλ, llc, dλ, ssλ, δt, srδt)
    # end
  end

  return llc, dλ, ssλ
end




"""
    update_tip!(tree::iTbX,
                α   ::Float64,
                σλ  ::Float64,
                llc ::Float64,
                dλ  ::Float64,
                ssλ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a `gbm` tip proposal.
"""
function update_tip!(tree::iTbX,
                     α   ::Float64,
                     σλ  ::Float64,
                     llc ::Float64,
                     dλ  ::Float64,
                     ssλ ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  @inbounds begin

    λc   = lλ(tree)
    l    = lastindex(λc)
    fdtp = fdt(tree)
    λp   = Vector{Float64}(undef, l)

    bm!(λp, λc[1], α, σλ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ = llr_gbm_b_sep(λp, λc, α, σλ, δt, fdtp, srδt, false)

    acr = llrbd

    if -randexp() < acr
      llc += llrbm + acr
      dλ  += λp[l] - λc[l]
      ssλ += ssrλ
      unsafe_copyto!(λc, 1, λp, 1, l)
    end
  end

  return llc, dλ, ssλ
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

    llr, acr, ssrλ = llr_propr(λpp, λ1p, λ2p, λpc, λ1c, λ2c,
      α, σλ, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < acr
      llc += llr
      dλ  += (λ1c[1] - λn)
      ssλ += ssrλ
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, dλ, ssλ
end




"""
    update_triad!(tree::iTbX,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ssλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` trio proposal.
"""
function update_triad!(tree::iTbX,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       dλ  ::Float64,
                       ssλ ::Float64,
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

    llr, acr, ssrλ = llr_propr(λpp, λ1p, λ2p, λpc, λ1c, λ2c,
      α, σλ, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < acr
      llc += llr
      dλ  += (λ1c[1] - λn)
      ssλ += ssrλ
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(λ2c, 1, λ2p, 1, l2)
    end
  end

  return llc, dλ, ssλ
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
function llr_propr(λpp  ::Array{Float64,1},
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
  llrbmp, llrbp, ssrλp = llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtp, srδt, true)
  llrbm1, llrb1, ssrλ1 = llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdt1, srδt, false)
  llrbm2, llrb2, ssrλ2 = llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdt2, srδt, false)

  acr  = llrbp + llrb1 + llrb2
  llr  = llrbmp + llrbm1 + llrbm2 + acr
  ssrλ = ssrλp + ssrλ1 + ssrλ2

  return llr, acr, ssrλ
end


