#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#




"""
    _daughters_update!(ξ1  ::iTfbd,
                       ξ2  ::iTfbd,
                       λf  ::Float64,
                       μf  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

Make a `gbm-bd` proposal for daughters from forwards simulated branch.
"""
function _daughters_update!(ξ1  ::iTfbd,
                            ξ2  ::iTfbd,
                            λf  ::Float64,
                            μf  ::Float64,
                            α   ::Float64,
                            σλ  ::Float64,
                            σμ  ::Float64,
                            δt  ::Float64,
                            srδt::Float64)
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

    # acceptance rate
    normprop =
      duoldnorm(λf, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ) -
      duoldnorm(λi, λ1 - α*e1, λ2 - α*e2, e1, e2, σλ) +
      duoldnorm(μf, μ1, μ2, e1, e2, σμ)               -
      duoldnorm(μi, μ1, μ2, e1, e2, σμ)

    # log likelihood ratios
    llrbm1, llrbd1, ssrλ1, ssrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt, false, false)
    llrbm2, llrbd2, ssrλ2, ssrμ2 =
      llr_gbm_b_sep(λ2p, μ2p, λ2c, μ2c, α, σλ, σμ, δt, fdt2, srδt, false, false)

    acr  = llrbd1 + llrbd2 + λf - λi
    llr  = llrbm1 + llrbm2 + acr
    acr += normprop
    drλ  = 2.0*(λi - λf)
    ssrλ = ssrλ1 + ssrλ2
    ssrμ = ssrμ1 + ssrμ2
  end

  return llr, acr, drλ, ssrλ, ssrμ, λ1p, λ2p, μ1p, μ2p
end




"""
    _update_gbm!(tree::iTfbd,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 llc ::Float64,
                 dλ  ::Float64,
                 ssλ ::Float64,
                 ssμ ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 ter ::Bool)

Do `gbm-bd` updates on a decoupled tree recursively.
"""
function _update_gbm!(tree::iTfbd,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      llc ::Float64,
                      dλ  ::Float64,
                      ssλ ::Float64,
                      ssμ ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ter ::Bool)

  if def1(tree)
    if def2(tree)
      llc, dλ, ssλ, ssμ =
        update_triad!(tree, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt)
      llc, dλ, ssλ, ssμ =
        _update_gbm!(tree.d1, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter)
      llc, dλ, ssλ, ssμ =
        _update_gbm!(tree.d2, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter)
    else
      llc, dλ, ssλ, ssμ =
        update_duo!(tree, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt)
      llc, dλ, ssλ, ssμ =
        _update_gbm!(tree.d1, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter)
    end
  else
    if !isfix(tree) || ter
      llc, dλ, ssλ, ssμ =
        update_tip!(tree, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt)
    end
  end

  return llc, dλ, ssλ, ssμ
end




"""
    update_tip!(tree::iTfbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                llc ::Float64,
                dλ  ::Float64,
                ssλ ::Float64,
                ssμ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a `gbm` tip proposal.
"""
function update_tip!(tree::iTfbd,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     llc ::Float64,
                     dλ  ::Float64,
                     ssλ ::Float64,
                     ssμ ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  @inbounds begin

    λc   = lλ(tree)
    μc   = lμ(tree)
    l    = lastindex(λc)
    fdtp = fdt(tree)
    λp   = Vector{Float64}(undef, l)
    μp   = Vector{Float64}(undef, l)

    bm!(λp, μp, λc[1], μc[1], α, σλ, σμ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ, ssrμ =
      llr_gbm_b_sep(λp, μp, λc, μc, α, σλ, σμ, δt, fdtp, srδt,
        false, isextinct(tree))

    acr = llrbd

    if -randexp() < acr
      llc += llrbm + acr
      dλ  += λp[l] - λc[l]
      ssλ += ssrλ
      ssμ += ssrμ
      unsafe_copyto!(λc, 1, λp, 1, l)
      unsafe_copyto!(μc, 1, μp, 1, l)
    end
  end

  return llc, dλ, ssλ, ssμ
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
                dλ  ::Float64,
                ssλ ::Float64,
                ssμ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a `gbm` duo proposal.
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
                     dλ  ::Float64,
                     ssλ ::Float64,
                     ssμ ::Float64,
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
    λi  = λ1c[1]
    λ1  = λ1c[l1]
    μp  = μpc[1]
    μi  = μ1c[1]
    μ1  = μ1c[l1]

    # node proposal
    λn = duoprop(λp + α*ep, λ1 - α*e1, ep, e1, σλ)
    μn = duoprop(μp, μ1, ep, e1, σμ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, μpp, μp, μn, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, μ1p, μn, μ1, σλ, σμ, δt, fdt1, srδt)

    # log likelihood ratios
    llrbmp, llrbdp, ssrλp, ssrμp =
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, α, σλ, σμ, δt, fdtp, srδt,
        true, false)
    llrbm1, llrbd1, ssrλ1, ssrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt,
        false, false)

    acr = llrbdp + llrbd1

    if -randexp() < acr
      llc += llrbmp + llrbm1 + acr
      dλ  += λi - λn
      ssλ += ssrλp + ssrλ1
      ssμ += ssrμp + ssrμ1
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(μpc, 1, μpp, 1, lp)
      unsafe_copyto!(μ1c, 1, μ1p, 1, l1)
      λi = λn
    end
  end

  return llc, dλ, ssλ, ssμ, λi
end




"""
    update_duo!(tree::iTfbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                llc ::Float64,
                dλ  ::Float64,
                ssλ ::Float64,
                ssμ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a `gbm` trio proposal.
"""
function update_duo!(tree::iTfbd,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     llc ::Float64,
                     dλ  ::Float64,
                     ssλ ::Float64,
                     ssμ ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  @inbounds begin

    λpc  = lλ(tree)
    λ1c  = lλ(tree.d1)
    μpc  = lμ(tree)
    μ1c  = lμ(tree.d1)
    lp   = lastindex(λpc)
    l1   = lastindex(λ1c)
    λpp  = Vector{Float64}(undef,lp)
    λ1p  = Vector{Float64}(undef,l1)
    μpp  = Vector{Float64}(undef,lp)
    μ1p  = Vector{Float64}(undef,l1)
    λp   = λpc[1]
    λ1   = λ1c[l1]
    μp   = μpc[1]
    μ1   = μ1c[l1]
    ep   = e(tree)
    e1   = e(tree.d1)
    fdtp = fdt(tree)
    fdt1 = fdt(tree.d1)

    # node proposal
    λn = trioprop(λp + α*ep, λ1 - α*e1, ep, e1, σλ)
    μn = trioprop(μp, μ1, ep, e1, σμ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, μpp, μp, μn, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, μ1p, μn, μ1, σλ, σμ, δt, fdt1, srδt)

    llrbmp, llrbdp, ssrλp, ssrμp =
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, α, σλ, σμ, δt, fdtp, srδt,
        true, false)
    llrbm1, llrbd1, ssrλ1, ssrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt,
        false, isextinct(tree.d1))

    acr = llrbdp + llrbd1

    if -randexp() < acr
      llc += llrbmp + llrbm1 + acr
      dλ  += (λ1c[1] - λn)
      ssλ += ssrλp + ssrλ1
      ssμ += ssrμp + ssrμ1
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(μpc, 1, μpp, 1, lp)
      unsafe_copyto!(μ1c, 1, μ1p, 1, l1)
    end
  end

  return llc, dλ, ssλ, ssμ
end



