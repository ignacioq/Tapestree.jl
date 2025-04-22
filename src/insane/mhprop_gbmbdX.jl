#=

Anagenetic GBM birth-death MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 27 05 2020
=#




"""
    _daughters_update!(ξ1  ::T,
                       ξ2  ::T,
                       λf  ::Float64,
                       μf  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       xp  ::Float64,
                       βλ  ::Float64,
                       σx  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

Make a `gbm-bd` proposal for daughters from forwards simulated branch.
"""
function _daughters_update!(ξ1  ::T,
                            ξ2  ::T,
                            λf  ::Float64,
                            μf  ::Float64,
                            α   ::Float64,
                            σλ  ::Float64,
                            σμ  ::Float64,
                            xp  ::Float64,
                            βλ  ::Float64,
                            σx  ::Float64,
                            δt  ::Float64,
                            srδt::Float64) where {T <: Tx}
  @inbounds begin

    λ1c  = lλ(ξ1)
    λ2c  = lλ(ξ2)
    μ1c  = lμ(ξ1)
    μ2c  = lμ(ξ2)
    x1c  = xv(ξ1)
    x2c  = xv(ξ2)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    μ1p  = Vector{Float64}(undef,l1)
    μ2p  = Vector{Float64}(undef,l2)
    x1p  = Vector{Float64}(undef,l1)
    x2p  = Vector{Float64}(undef,l2)
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
    bb!(x1p, xp, x1c[l1], σx, δt, fdt1, srδt)
    bb!(x2p, xp, x2c[l2], σx, δt, fdt2, srδt)

    # log likelihood ratios
    llrbm1, llrbd1, ssrλ1, ssrμ1, ssrx1 =
      llr_gbm_b_sep(λ1p, μ1p, x1p, λ1c, μ1c, x1c, α, σλ, σμ, βλ, σx, 
        δt, fdt1, srδt, false, false)
    llrbm2, llrbd2, ssrλ2, ssrμ2, ssrx2 =
      llr_gbm_b_sep(λ2p, μ2p, x2p, λ2c, μ2c, x2c, α, σλ, σμ, βλ, σx, 
        δt, fdt2, srδt, false, false)

    acr  = llrbd1 + llrbd2
    llr  = llrbm1 + llrbm2 + λf - λi + acr
    drλ  = 2.0*(λi - λf)
    ssrλ = ssrλ1 + ssrλ2
    ssrμ = ssrμ1 + ssrμ2
    ssrx = ssrx1 + ssrx2
  end

  return llr, acr, drλ, ssrλ, ssrμ, ssrx, λ1p, λ2p, μ1p, μ2p, x1p, x2p
end






