#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#




"""
    _fstem_update!(ξi   ::iTfbd,
                   ξ1   ::iTfbd,
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
                   surv ::Int64)

Do `gbm-bd` update for fossil stem root.
"""
function _fstem_update!(ξi   ::iTfbd,
                        ξ1   ::iTfbd,
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
                        surv ::Int64)

  @inbounds begin
    λpc  = lλ(ξi)
    μpc  = lμ(ξi)
    λi   = λpc[1]
    μi   = μpc[1]
    λ1c  = lλ(ξ1)
    μ1c  = lμ(ξ1)
    l1   = lastindex(λ1c)
    λ1p  = Vector{Float64}(undef,l1)
    μ1p  = Vector{Float64}(undef,l1)
    λ1   = λ1c[l1]
    μ1   = μ1c[l1]
    el   = e(ξ1)
    sqre = sqrt(el)
    fdt1 = fdt(ξ1)

    # node proposal
    λr = rnorm(λ1 - α*el, σλ*sqre)
    μr = rnorm(μ1, σμ*sqre)

    # simulate fix tree vector
    bb!(λ1p, λr, λ1, μ1p, μr, μ1, σλ, σμ, δt, fdt1, srδt)

    # log likelihood ratios
    llrbm, llrbd, ssrλ, ssrμ, irrλ, irrμ =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt, false, false)

    # log probability
    lU = -randexp()

    llr = llrbd

    if lU < llr + log(1000.0/mc)

      #survival
      mp   = m_surv_gbmbd(th, λr, μr, α, σλ, σμ, δt, srδt, 1_000, surv)
      llr += log(mp/mc)

      if lU < llr
        llc += llrbm + llr
        ddλ  += λ1c[1] - λr
        ssλ += ssrλ
        ssμ += ssrμ
        irλ += irrλ
        irμ += irrμ
        mc   = mp
        fill!(λpc, λr)
        fill!(μpc, μr)
        unsafe_copyto!(λ1c, 1, λ1p, 1, l1)
        unsafe_copyto!(μ1c, 1, μ1p, 1, l1)
      end
    end
  end

  return llc, ddλ, ssλ, ssμ, irλ, irμ, mc
end




"""
    _update_gbm!(tree::iTfbd,
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

Do `gbm-bd` updates on a decoupled tree recursively.
"""
function _update_gbm!(tree::iTfbd,
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
    if def2(tree)
      llc, ddλ, ssλ, ssμ, irλ, irμ =
        update_triad!(tree, α, σλ, σμ, llc, ddλ, ssλ, ssμ, irλ, irμ, δt, srδt)

      llc, ddλ, ssλ, ssμ, irλ, irμ =
        _update_gbm!(tree.d1, α, σλ, σμ, llc, ddλ, ssλ, ssμ, irλ, irμ, 
          δt, srδt, ter)
      llc, ddλ, ssλ, ssμ, irλ, irμ =
        _update_gbm!(tree.d2, α, σλ, σμ, llc, ddλ, ssλ, ssμ, irλ, irμ, 
          δt, srδt, ter)
    else

      if e(tree.d1) < (δt + accerr) && isextinct(tree.d1)
        llc, ddλ, ssλ, ssμ, irλ, irμ =
          update_fduo!(tree, α, σλ, σμ, llc, ddλ, ssλ, ssμ, irλ, irμ, δt, srδt)
      else
        llc, ssλ, ssμ, irλ, irμ =
          update_duo!(tree, α, σλ, σμ, llc, ssλ, ssμ, irλ, irμ, δt, srδt)
      end

      llc, ddλ, ssλ, ssμ, irλ, irμ =
        _update_gbm!(tree.d1, α, σλ, σμ, llc, ddλ, ssλ, ssμ, irλ, irμ, 
          δt, srδt, ter)
    end
  elseif !isfix(tree) || ter
    llc, ddλ, ssλ, ssμ, irλ, irμ = 
      update_tip!(tree, α, σλ, σμ, llc, ddλ, ssλ, ssμ, irλ, irμ, δt, srδt)
  end

  return llc, ddλ, ssλ, ssμ, irλ, irμ
end




"""
    update_duo!(tree::iTfbd,
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

Make a `gbm` duo proposal.
"""
function update_duo!(tree::iTfbd,
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
    update_fduo!(tree::iTfbd,
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

Make a `gbm` duo proposal.
"""
function update_fduo!(tree::iTfbd,
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
    ec   = ep + e1
    fdtp = fdt(tree)
    fdt1 = fdt(tree.d1)
    sqre = sqrt(ec)

    # node proposal
    λi = rnorm(λp + α*ec, σλ*sqre)
    μi = rnorm(μp, σμ*sqre)
    λn = duoprop(λp + α*ep, λi - α*e1, ep, e1, σλ)
    μn = duoprop(μp, μi, ep, e1, σμ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, μpp, μp, μn, σλ, σμ, δt, fdtp, srδt)
    bb!(λ1p, λn, λi, μ1p, μn, μi, σλ, σμ, δt, fdt1, srδt)

    llrbmp, llrbdp, ssrλp, ssrμp, irrλp, irrμp =
      llr_gbm_b_sep(λpp, μpp, λpc, μpc, α, σλ, σμ, δt, fdtp, srδt,
        false, false)
    llrbm1, llrbd1, ssrλ1, ssrμ1, irrλ1, irrμ1 =
      llr_gbm_b_sep(λ1p, μ1p, λ1c, μ1c, α, σλ, σμ, δt, fdt1, srδt,
        false, true)

    acr = llrbdp + llrbd1

    if -randexp() < acr
      llc += llrbmp + llrbm1 + acr
      ddλ  += λi - λ1
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

  return llc, ddλ, ssλ, ssμ, irλ, irμ
end


