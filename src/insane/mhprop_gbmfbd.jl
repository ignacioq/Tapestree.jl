#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#



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
      llc, ssλ, ssμ =
        update_duo!(tree, α, σλ, σμ, llc, ssλ, ssμ, δt, srδt)
      llc, dλ, ssλ, ssμ =
        _update_gbm!(tree.d1, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter)
    end
  elseif !isfix(tree) || ter

    llc, dλ, ssλ, ssμ = 
      update_tip!(tree, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt)
  end

  return llc, dλ, ssλ, ssμ
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

Make a `gbm` duo proposal.
"""
function update_duo!(tree::iTfbd,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     llc ::Float64,
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
    λn = duoprop(λp + α*ep, λ1 - α*e1, ep, e1, σλ)
    μn = duoprop(μp, μ1, ep, e1, σμ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, μpp, μp, μn, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, μ1p, μn, μ1, σλ, σμ, δt, fdt1, srδt)

    llrbmp, llrbdp, ssrλp, ssrμp =
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, α, σλ, σμ, δt, fdtp, srδt,
        false, false)
    llrbm1, llrbd1, ssrλ1, ssrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt,
        false, false)

    acr = llrbdp + llrbd1

    if -randexp() < acr
      llc += llrbmp + llrbm1 + acr
      ssλ += ssrλp + ssrλ1
      ssμ += ssrμp + ssrμ1
      unsafe_copyto!(λpc, 1, λpp, 1, lp)
      unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
      unsafe_copyto!(μpc, 1, μpp, 1, lp)
      unsafe_copyto!(μ1c, 1, μ1p, 1, l1)
    end
  end

  return llc, ssλ, ssμ
end

