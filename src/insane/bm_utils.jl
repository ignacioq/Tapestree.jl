#=

Brownian motion utilities

Ignacio Quintero Mächler

t(-_-t)

Created 10 09 2020
=#




"""
    bm!(tree::T,
        bbiλ::Array{Float64,1},
        ii  ::Int64,
        tl  ::Int64,
        σλ  ::Float64,
        srδt::Float64) where {T <: iTgbm}

Fill fix with previously simulated geometric Brownian motion and simulate
geometric Brownian motion in place for the unfixed trees.
"""
@inline function bm!(tree::T,
                     bbiλ::Array{Float64,1},
                     ii  ::Int64,
                     tl  ::Int64,
                     σλ  ::Float64,
                     srδt::Float64) where {T <: iTgbm}

  @inbounds begin

    λv  = lλ(tree)
    l   = lastindex(λv)
    fi  = ii + l - 1
    cix = ii:fi

    @simd for i in Base.OneTo(l)
      λv[i] = bbiλ[cix[i]]
    end

    if fi < tl
      if isfix(tree.d1)
        bm!(tree.d1::T, bbiλ, fi, tl, σλ, srδt)
        bm!(tree.d2::T, λv[l], σλ, srδt)
      elseif isfix(tree.d2)
        bm!(tree.d1::T, λv[l], σλ, srδt)
        bm!(tree.d2::T, bbiλ, fi, tl, σλ, srδt)
      end
    end

  end

  return nothing
end




"""
  bm!(tree::T,
      λt  ::Float64,
      σλ  ::Float64,
      srδt::Float64) where {T <: iTgbm}

Simulate birth-death geometric Brownian motion in place.
"""
function bm!(tree::T,
             λt  ::Float64,
             σλ  ::Float64,
             srδt::Float64) where {T <: iTgbm}

  λv  = lλ(tree)

  bm!(λv, λt, fdt(tree), σλ, srδt)

  l = lastindex(λv)

  if !isnothing(tree.d1)
    bm!(tree.d1::T, λv[l], σλ, srδt)
  end

  if !isnothing(tree.d2)
    bm!(tree.d2::T, λv[l], σλ, srδt)
  end

  return nothing
end




"""
    bm!(tree::iTgbmbd,
        bbiλ::Array{Float64,1},
        bbiμ::Array{Float64,1},
        ii  ::Int64,
        tl  ::Int64,
        σλ  ::Float64,
        σμ  ::Float64,
        srδt::Float64)

Fill fix with previously simulated geometric Brownian motion and simulate
geometric Brownian motion in place for the unfixed trees.
"""
@inline function bm!(tree::iTgbmbd,
                     bbiλ::Array{Float64,1},
                     bbiμ::Array{Float64,1},
                     ii  ::Int64,
                     tl  ::Int64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     srδt::Float64)

  @inbounds begin

    λv  = lλ(tree)
    μv  = lμ(tree)
    l   = lastindex(λv)
    fi  = ii + l - 1
    cix = ii:fi

    @simd for i in Base.OneTo(l)
      λv[i] = bbiλ[cix[i]]
      μv[i] = bbiμ[cix[i]]
    end

    if fi < tl
      if isfix(tree.d1)
        bm!(tree.d1::iTgbmbd, bbiλ, bbiμ, fi, tl, σλ, σμ, srδt)
        bm!(tree.d2::iTgbmbd, λv[l], μv[l], σλ, σμ, srδt)
      elseif isfix(tree.d2)
        bm!(tree.d1::iTgbmbd, λv[l], μv[l], σλ, σμ, srδt)
        bm!(tree.d2::iTgbmbd, bbiλ, bbiμ, fi, tl, σλ, σμ, srδt)
      end
    end

  end

  return nothing
end




"""
  bm!(tree::iTgbmbd,
      λt  ::Float64,
      μt  ::Float64,
      σλ  ::Float64,
      σμ  ::Float64,
      srδt::Float64)

Simulate birth-death geometric Brownian motion in place.
"""
function bm!(tree::iTgbmbd,
             λt  ::Float64,
             μt  ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             srδt::Float64)

  λv  = lλ(tree)
  μv  = lμ(tree)

  bm!(λv, μv, λt, μt, fdt(tree), σλ, σμ, srδt)

  l = lastindex(λv)

  if !isnothing(tree.d1)
    bm!(tree.d1::iTgbmbd, λv[l], μv[l], σλ, σμ, srδt)
  end

  if !isnothing(tree.d2)
    bm!(tree.d2::iTgbmbd, λv[l], μv[l], σλ, σμ, srδt)
  end

  return nothing
end




"""
    ll_bm(x   ::Array{Float64,1},
          fdt ::Float64,
          σ   ::Float64, 
          srδt::Float64)

Returns the log-likelihood for Brownian motion.
"""
@inline function ll_bm(x   ::Array{Float64,1},
                       fdt ::Float64,
                       σ   ::Float64, 
                       srδt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(x)-2

    ll = 0.0
    xi = x[1]
    @simd for i in Base.OneTo(nI)
      xi1 = x[i+1]
      ll += (xi1 - xi)^2
      xi  = xi1
    end

    # add to global likelihood
    ll *= (-0.5/((σ*srδt)^2))
    ll -= Float64(nI)*(log(σ*srδt) + 0.5*log(2.0π))

    # add final non-standard `δt`
    ll += ldnorm_bm(x[nI+2], x[nI+1], sqrt(fdt)*σ)
  end

  return ll
end





"""
    llr_bm(xp  ::Array{Float64,1},
           xc  ::Array{Float64,1},
           fdt::Float64,
           σ   ::Float64, 
           srδt::Float64)

Returns the log-likelihood ratio for Brownian motion.
"""
@inline function llr_bm(xp  ::Array{Float64,1},
                        xc  ::Array{Float64,1},
                        fdt::Float64,
                        σ   ::Float64, 
                        srδt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(xp)-2

    llr = 0.0
    xpi = xp[1]
    xci = xc[1]
    @simd for i in Base.OneTo(nI)
      xpi1 = xp[i+1]
      xci1 = xc[i+1]
      llr += (xpi1 - xpi)^2 - (xci1 - xci)^2
      xpi  = xpi1
      xci  = xci1
    end

    # add to global likelihood
    llr *= (-0.5/((σ*srδt)^2))

    # add final non-standard `δt`
    llr += lrdnorm_bm_x(xp[nI+2], xp[nI+1], xc[nI+2], xc[nI+1], sqrt(fdt)*σ)
  end

  return llr
end




"""
    bm!(x0   ::Array{Float64,1},
        x1   ::Array{Float64,1},
        x0i  ::Float64,
        x1i  ::Float64,
        fdt ::Float64,
        σ0   ::Float64,
        σ1   ::Float64,
        srδt::Float64)

Brownian motion simulation function for updating a branch for two 
vectors that share times in place.
"""
@inline function bm!(x0   ::Array{Float64,1},
                     x1   ::Array{Float64,1},
                     x0i  ::Float64,
                     x1i  ::Float64,
                     fdt ::Float64,
                     σ0   ::Float64,
                     σ1   ::Float64,
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
    srfdt  = sqrt(fdt)
    x0[l] *= srfdt*σ0
    x1[l] *= srfdt*σ1

    cumsum!(x0, x0)
    cumsum!(x1, x1)
  end

  return nothing
end




"""
    bm!(x   ::Array{Float64,1},
        xi  ::Float64,
        fdt::Float64,
        σ   ::Float64,
        srδt::Float64)

Brownian motion simulation function for updating a branch in place.
"""
@inline function bm!(x   ::Array{Float64,1},
                     xi  ::Float64,
                     fdt::Float64,
                     σ   ::Float64,
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
  end

  return nothing
end




"""
    bb!(x   ::Array{Float64,1},
        xi  ::Float64,
        xf  ::Float64,
        fdt::Float64,
        σ   ::Float64,
        δt  ::Float64,
        srδt::Float64)

Brownian bridge simulation function for updating a branch in place.
"""
@inline function bb!(x   ::Array{Float64,1},
                     xi  ::Float64,
                     xf  ::Float64,
                     fdt ::Float64,
                     σ   ::Float64,
                     δt  ::Float64,
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
    ite = 1.0/(Float64(l-2) * δt + fdt)
    xdf = (x[l] - xf)

    @simd for i = Base.OneTo(l-1)
      x[i] -= (Float64(i-1) * δt * ite * xdf)
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
        fdt::Float64,
        σ0  ::Float64,
        σ1  ::Float64,
        δt  ::Float64,
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
                     fdt::Float64,
                     σ0  ::Float64,
                     σ1  ::Float64,
                     δt  ::Float64,
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
    ite = 1.0/(Float64(l-2) * δt + fdt)
    x0df = (x0[l] - x0f)
    x1df = (x1[l] - x1f)

    @simd for i = Base.OneTo(l-1)
      iti    = Float64(i-1) * δt * ite
      x0[i] -= (iti * x0df)
      x1[i] -= (iti * x1df)
    end

    # for last non-standard δt
    x0[l] = x0f
    x1[l] = x1f
  end

  return nothing
end





"""
   sim_bm(xa  ::Float64, 
          σ   ::Float64, 
          srδt::Float64, 
          nt  ::Int64, 
          fdt::Float64)

Returns a Brownian motion vector starting in `xa`, with diffusion rate
`σ` and times `t`. 
"""
@inline function sim_bm(xa  ::Float64, 
                        σ   ::Float64, 
                        srδt::Float64, 
                        nt  ::Int64, 
                        fdt::Float64)
  @inbounds begin
    l = nt + 2
    x = randn(l)
    # for standard δt
    x[1] = xa
    @simd for i in Base.OneTo(l-2)
      x[i+1] *= srδt*σ
    end
    x[l] *= sqrt(fdt)*σ
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
  return ldnorm_bm(x, td2 * invt * xd1 + td1 * invt * xd2, 
    sqrt(td1 * td2 * invt)*σ)
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
    lrdnorm_bm_σ(x::Float64, μ::Float64, σpsrt::Float64, σcsrt::Float64)

Compute the **Normal** density ratio in logarithmic scale with mean `μ` 
and proposal standard density `σpsrt` and current `σcsrt` for `x`.
"""
lrdnorm_bm_σ(x::Float64, μ::Float64, σpsrt::Float64, σcsrt::Float64) =
  - log(σpsrt/σcsrt) - 0.5*(x - μ)^2*(1.0/σpsrt^2 - 1.0/σcsrt^2)




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


