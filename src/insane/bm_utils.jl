#=

Brownian motion utilities

Ignacio Quintero Mächler

t(-_-t)

Created 10 09 2020
=#




"""
    bm!(tree::T,
        λt  ::Float64,
        α   ::Float64,
        σλ  ::Float64,
        δt  ::Float64,
        srδt::Float64) where {T <: iT}

Simulate birth-death geometric Brownian motion in place.
"""
function bm!(tree::T,
             λt  ::Float64,
             α   ::Float64,
             σλ  ::Float64,
             δt  ::Float64,
             srδt::Float64) where {T <: iT}

  λv  = lλ(tree)

  bm!(λv, λt, α, σλ, δt, fdt(tree), srδt)

  l = lastindex(λv)

  if def1(tree)
    bm!(tree.d1::T, λv[l], α, σλ, δt, srδt)
    bm!(tree.d2::T, λv[l], α, σλ, δt, srδt)
  end

  return nothing
end




"""
    bm!(tree::iTbd,
        λt  ::Float64,
        μt  ::Float64,
        α   ::Float64,
        σλ  ::Float64,
        σμ  ::Float64,
        δt  ::Float64,
        srδt::Float64)

Simulate birth-death geometric Brownian motion in place.
"""
@inline function bm!(tree::iTbd,
                     λt  ::Float64,
                     μt  ::Float64,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  @inbounds begin

    λv  = lλ(tree)
    μv  = lμ(tree)

    bm!(λv, μv, λt, μt, α, σλ, σμ, δt, fdt(tree), srδt)

    l = lastindex(λv)

    if def1(tree)
      bm!(tree.d1, λv[l], μv[l], α, σλ, σμ, δt, srδt)
      bm!(tree.d2, λv[l], μv[l], α, σλ, σμ, δt, srδt)
    end
  end

  return nothing
end





"""
    ll_bm(x  ::Array{Float64,1},
          α  ::Float64,
          σ  ::Float64,
          δt ::Float64,
          fdt::Float64,
          srδt::Float64)

Returns the log-likelihood for Brownian motion with drift.
"""
function ll_bm(x  ::Array{Float64,1},
               α  ::Float64,
               σ  ::Float64,
               δt ::Float64,
               fdt::Float64,
               srδt::Float64)

  # estimate standard `δt` likelihood
  nI = lastindex(x)-2

  ll = 0.0
  @avx for i in Base.OneTo(nI)
    ll += (x[i+1] - x[i] - α*δt)^2
  end

  # add to global likelihood
  ll *= (-0.5/((σ*srδt)^2))
  ll -= Float64(nI)*(log(σ*srδt) + 0.5*log(2.0π))

  # add final non-standard `δt`
  ll += ldnorm_bm(x[nI+2], x[nI+1] + α*fdt, sqrt(fdt)*σ)

  return ll
end




"""
    llr_bm(xp  ::Array{Float64,1},
           xc  ::Array{Float64,1},
           α   ::Float64,
           σ   ::Float64,
           δt  ::Float64,
           fdt ::Float64,
           srδt::Float64)

Returns the log-likelihood ratio for Brownian motion.
"""
function llr_bm(xp  ::Array{Float64,1},
                xc  ::Array{Float64,1},
                α   ::Float64,
                σ   ::Float64,
                δt  ::Float64,
                fdt ::Float64,
                srδt::Float64)

  # estimate standard `δt` likelihood
  nI = lastindex(xp) - 2

  llr = 0.0
  @avx for i in Base.OneTo(nI)
    llr += (xp[i+1] - xp[i] - α*δt)^2 -
           (xc[i+1] - xc[i] - α*δt)^2
  end

  # add to global likelihood
  llr *= -0.5/((σ*srδt)^2)

  # add final non-standard `δt`
  llr += lrdnorm_bm_x(xp[nI+2], xp[nI+1] + α*fdt,
                      xc[nI+2], xc[nI+1] + α*fdt, sqrt(fdt)*σ)

  return llr
end




"""
    bm!(x0  ::Array{Float64,1},
        x1  ::Array{Float64,1},
        x0i ::Float64,
        x1i ::Float64,
        α   ::Float64,
        σ0  ::Float64,
        σ1  ::Float64,
        δt  ::Float64,
        fdt ::Float64,
        srδt::Float64)

Brownian motion simulation function for updating a branch for two
vectors that share times and x0 follows drift α.
"""
@inline function bm!(x0  ::Array{Float64,1},
                     x1  ::Array{Float64,1},
                     x0i ::Float64,
                     x1i ::Float64,
                     α   ::Float64,
                     σ0  ::Float64,
                     σ1  ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin
    l = lastindex(x0)

    randn!(x0)
    randn!(x1)

    # for standard δt
    x0[1] = x0i
    x1[1] = x1i
    @simd for i = Base.OneTo(l-2)
      x0[i+1] *= srδt*σ0
      x0[i+1] += α*δt
      x1[i+1] *= srδt*σ1
    end
    srfdt  = sqrt(fdt)
    x0[l] *= srfdt*σ0
    x0[l] += α*fdt
    x1[l] *= srfdt*σ1

    cumsum!(x0, x0)
    cumsum!(x1, x1)
  end

  return nothing
end




"""
    bm!(x   ::Array{Float64,1},
        xi  ::Float64,
        α   ::Float64,
        σ   ::Float64,
        δt  ::Float64,
        fdt ::Float64,
        srδt::Float64)

Brownian motion simulation function for updating a branch in place.
"""
@inline function bm!(x   ::Array{Float64,1},
                     xi  ::Float64,
                     α   ::Float64,
                     σ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(x)

    # for standard δt
    x[1] = xi
    @simd for i = Base.OneTo(l-2)
      x[i+1] *= srδt*σ
      x[i+1] += α*δt
    end
    x[l] *= sqrt(fdt)*σ
    x[l] += α*fdt
    cumsum!(x, x)
  end

  return nothing
end




"""
    bb!(x   ::Array{Float64,1},
        xi  ::Float64,
        xf  ::Float64,
        σ   ::Float64,
        δt  ::Float64,
        fdt::Float64,
        srδt::Float64)

Brownian bridge simulation function for updating a branch in place.
"""
@inline function bb!(x   ::Array{Float64,1},
                     xi  ::Float64,
                     xf  ::Float64,
                     σ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin

    l = lastindex(x)

    randn!(x)

    # for standard δt
    x[1] = xi
    @simd for i = Base.OneTo(l-2)
      x[i+1] *= srδt*σ
    end
    x[l] *= sqrt(fdt)*σ
    cumsum!(x, x)

    # make bridge
    if l > 2
      ite = 1.0/(Float64(l-2) * δt + fdt)
      xdf = (x[l] - xf)

      @avx for i = Base.OneTo(l-1)
        x[i] -= (Float64(i-1) * δt * ite * xdf)
      end
    end
    # for last non-standard δt
    x[l] = xf
  end

  return nothing
end




"""
    bb!(x0  ::Array{Float64,1},
        x0i ::Float64,
        x0f ::Float64,
        x1  ::Array{Float64,1},
        x1i ::Float64,
        x1f ::Float64,
        σ0  ::Float64,
        σ1  ::Float64,
        δt  ::Float64,
        fdt ::Float64,
        srδt::Float64)

Brownian bridge simulation function for updating two vectors
(`0` & `1`) with shared times in place.
"""
@inline function bb!(x0  ::Array{Float64,1},
                     x0i ::Float64,
                     x0f ::Float64,
                     x1  ::Array{Float64,1},
                     x1i ::Float64,
                     x1f ::Float64,
                     σ0  ::Float64,
                     σ1  ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin
    l = lastindex(x0)

    randn!(x0)
    randn!(x1)

    # for standard δt
    x0[1] = x0i
    x1[1] = x1i
    @simd for i = Base.OneTo(l-2)
      x0[i+1] *= srδt*σ0
      x1[i+1] *= srδt*σ1
    end
    srlt  = sqrt(fdt)
    x0[l] *= srlt*σ0
    x1[l] *= srlt*σ1
    cumsum!(x0, x0)
    cumsum!(x1, x1)

    # make bridge
    if l > 2
      ite = 1.0/(Float64(l-2) * δt + fdt)
      x0df = (x0[l] - x0f)
      x1df = (x1[l] - x1f)

      @avx for i = Base.OneTo(l-1)
        iti    = Float64(i-1) * δt * ite
        x0[i] -= (iti * x0df)
        x1[i] -= (iti * x1df)
      end
    end

    # for last non-standard δt
    x0[l] = x0f
    x1[l] = x1f
  end

  return nothing
end




"""
    sim_bm(xa  ::Float64,
           α   ::Float64,
           σ   ::Float64,
           δt  ::Float64,
           fdt ::Float64,
           srδt::Float64,
           nt  ::Int64)

Returns a Brownian motion vector starting in `xa`, with diffusion rate
`σ` and times `t`.
"""
@inline function sim_bm(xa  ::Float64,
                        α   ::Float64,
                        σ   ::Float64,
                        δt  ::Float64,
                        fdt ::Float64,
                        srδt::Float64,
                        nt  ::Int64)
  @inbounds begin
    l = nt + 2
    x = randn(l)
    # for standard δt
    x[1] = xa
    @simd for i in Base.OneTo(nt)
      x[i+1] *= srδt*σ
      x[i+1] += α*δt
    end
    x[l] *= sqrt(fdt)*σ
    x[l] += α*fdt
    cumsum!(x, x)
  end

  return x
end




"""
    duoprop(xd1::Float64,
            xd2::Float64,
            td1::Float64,
            td2::Float64,
            σ  ::Float64)

Proposal for a duo of Gaussians.
"""
function duoprop(xd1::Float64,
                 xd2::Float64,
                 td1::Float64,
                 td2::Float64,
                 σ  ::Float64)

  invt = 1.0/(td1 + td2)
  return rnorm((td2 * invt * xd1 + td1 * invt * xd2),
               sqrt(td1 * td2 * invt)*σ)
end




"""
    trioprop(xpr::Float64,
             xd1::Float64,
             xd2::Float64,
             tpr::Float64,
             td1::Float64,
             td2::Float64,
             σ  ::Float64)

Proposal for a trio of Gaussians.
"""
function trioprop(xpr::Float64,
                  xd1::Float64,
                  xd2::Float64,
                  tpr::Float64,
                  td1::Float64,
                  td2::Float64,
                  σ  ::Float64)

    t = 1.0/(1.0/tpr + 1.0/td1 + 1.0/td2)
    return rnorm((xpr/tpr + xd1/td1 + xd2/td2)*t,
                 sqrt(t)*σ)
end




"""
    duodnorm(x  ::Float64,
             xd1::Float64,
             xd2::Float64,
             td1::Float64,
             td2::Float64,
             σ  ::Float64)
Likelihood for a duo of Gaussians.
"""
function duodnorm(x  ::Float64,
                  xd1::Float64,
                  xd2::Float64,
                  td1::Float64,
                  td2::Float64,
                  σ  ::Float64)

  invt = 1.0/(td1 + td2)
  return dnorm_bm(x, td2 * invt * xd1 + td1 * invt * xd2,
    sqrt(td1 * td2 * invt)*σ)
end




"""
    duoldnorm(x  ::Float64,
              xd1::Float64,
              xd2::Float64,
              td1::Float64,
              td2::Float64,
              σ  ::Float64)

Likelihood for a duo of Gaussians.
"""
function duoldnorm(x  ::Float64,
                   xd1::Float64,
                   xd2::Float64,
                   td1::Float64,
                   td2::Float64,
                   σ  ::Float64)

  invt = 1.0/(td1 + td2)
  return ldnorm_bm(x, invt * (td2 * xd1 + td1 * xd2), sqrt(td1 * td2 * invt)*σ)
end




"""
    trioldnorm(x  ::Float64,
               xpr::Float64,
               xd1::Float64,
               xd2::Float64,
               tpr::Float64,
               td1::Float64,
               td2::Float64,
               σ  ::Float64)

Likelihood for a trio of Gaussians.
"""
function trioldnorm(x  ::Float64,
                    xpr::Float64,
                    xd1::Float64,
                    xd2::Float64,
                    tpr::Float64,
                    td1::Float64,
                    td2::Float64,
                    σ  ::Float64)

    t = 1.0/(1.0/tpr + 1.0/td1 + 1.0/td2)
    return ldnorm_bm(x, (xpr/tpr + xd1/td1 + xd2/td2)*t, sqrt(t)*σ)
end




"""
    ldnorm_bm(x::Float64, μ::Float64, σsrt::Float64)

Compute the **Normal** density in logarithmic scale with
mean `μ` and standard density `σsrt` for `x`.
"""
ldnorm_bm(x::Float64, μ::Float64, σsrt::Float64) =
  -0.5*log(2.0π) - log(σsrt) - 0.5*((x - μ)/σsrt)^2




"""
    ldnorm_bm(x::Float64, μ::Float64, σsrt::Float64)

Compute the **Normal** density in with
mean `μ` and standard density `σsrt` for `x`.
"""
dnorm_bm(x::Float64, μ::Float64, σsrt::Float64) =
  1.0/(σsrt*sqrt(2.0π)) * exp(-0.5*((x - μ)/σsrt)^2)




"""
    lrdnorm_bm_x(xp::Float64,
                 μp::Float64,
                 xc::Float64,
                 μc::Float64,
                 σsrt::Float64)

Compute the **Normal** density ratio in logarithmic scale with
standard density `σ,` proposal mean `μp` and current mean `μc` for `xp`
and `xc`, respectively.
"""
lrdnorm_bm_x(xp::Float64, μp::Float64, xc::Float64, μc::Float64, σsrt::Float64) =
  -0.5*((xp - μp)^2 - (xc - μc)^2)/σsrt^2




"""
    lrdnorm_bm_x(xp::Float64,
                 xc::Float64,
                 μ ::Float64,
                 σsrt::Float64)

Compute the **Normal** density ratio in logarithmic scale with
standard density `σ` and mean `μ` for `xp`
and `xc`, respectively.
"""
lrdnorm_bm_x(xp::Float64, xc::Float64, μ::Float64, σsrt::Float64) =
  -0.5*((xp - μ)^2 - (xc - μ)^2)/σsrt^2




"""
    pnorm(x::Float64, y::Float64, μ::Float64, σ::Float64)

Cumulative probability between `x` and `y`, with `x < y` for a
**Normal** Distribution with mean `μ` and standard deviation `σ`.
"""
function pnorm(x::Float64, y::Float64, μ::Float64, σ::Float64)
  iσ = 1.0/(σ*sqrt(2.0))
  0.5*(erf((x - μ)*iσ, (y - μ)*iσ))
end




"""
    ncrep!(xp::Array{Float64,1},
           xc::Array{Float64,1},
           t ::Array{Float64,1},
           σ ::Float64)

Non-centered reparametization of data augmentation for `σ`, based on
Roberts and Stramer (2001).
"""
function ncrep!(xp::Array{Float64,1},
                xc::Array{Float64,1},
                t ::Array{Float64,1},
                σ ::Float64)

  @inbounds begin
    l   = lastindex(xc)
    xi  = xc[1]
    xf  = xc[l]
    itf = 1.0/t[l]
    iσ  = 1.0/σ

    for i in 2:(l-1)
      xp[i] = (xc[i] - xi + t[i]*(xi - xf)*itf)*iσ
    end

    xp[1] = xp[l] = 0.0
  end

  return nothing
end


