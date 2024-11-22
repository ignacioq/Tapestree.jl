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
  if nI > 0
    @turbo for i in Base.OneTo(nI)
      ll += (x[i+1] - x[i] - α*δt)^2
    end

    # add to global likelihood
    ll *= (-0.5/((σ*srδt)^2))
    ll -= Float64(nI)*(log(σ*srδt) + 0.5*log(2.0π))
  end

  # add final non-standard `δt`
  if fdt > 0.0
    ll += ldnorm_bm(x[nI+2], x[nI+1] + α*fdt, sqrt(fdt)*σ)
  end

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
  if nI > 0
    @turbo for i in Base.OneTo(nI)
      llr += (xp[i+1] - xp[i] - α*δt)^2 -
             (xc[i+1] - xc[i] - α*δt)^2
    end

    # add to global likelihood
    llr *= -0.5/((σ*srδt)^2)
  end

  # add final non-standard `δt`
  if fdt > 0.0
    llr += lrdnorm_bm_x(xp[nI+2], xp[nI+1] + α*fdt,
                        xc[nI+2], xc[nI+1] + α*fdt, sqrt(fdt)*σ)
  end

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
    bm!(x0  ::Array{Float64,1},
        x1  ::Array{Float64,1},
        x0i ::Float64,
        x1i ::Float64,
        α0  ::Float64,
        α1  ::Float64,
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
                     α0  ::Float64,
                     α1  ::Float64,
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
      x0[i+1] += α0*δt
      x1[i+1] *= srδt*σ1
      x1[i+1] += α1*δt
    end
    srfdt  = sqrt(fdt)
    x0[l] *= srfdt*σ0
    x0[l] += α0*fdt
    x1[l] *= srfdt*σ1
    x1[l] += α1*fdt

    cumsum!(x0, x0)
    cumsum!(x1, x1)
  end

  return nothing
end




"""
    bm!(x   ::Array{Float64,1},
        xi  ::Float64,
        σ   ::Float64,
        δt  ::Float64,
        fdt ::Float64,
        srδt::Float64)

Brownian motion without drift simulation function for updating a branch 
in place.
"""
@inline function bm!(x   ::Array{Float64,1},
                     xi  ::Float64,
                     σ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(x)
    # for standard δt
    x[1] = xi
    if l > 2
      @turbo for i = Base.OneTo(l-2)
        x[i+1] *= srδt*σ
      end
    end
    x[l] *= sqrt(fdt)*σ
    cumsum!(x, x)
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
    if l > 2
      @turbo for i = Base.OneTo(l-2)
        x[i+1] *= srδt*σ
        x[i+1] += α*δt
      end
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
    x[1] = xi
    if l > 2
      for i = Base.OneTo(l-2)
        x[i+1] *= srδt*σ
        x[i+1] += x[i]
      end
      x[l] *= sqrt(fdt)*σ
      x[l] += x[l-1]

      # make bridge
      ite = 1.0/(Float64(l-2) * δt + fdt)
      xdf = (x[l] - xf)

      @turbo for i = Base.OneTo(l-1)
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
    if l > 2
      for i = Base.OneTo(l-2)
        x0[i+1] *= srδt*σ0
        x0[i+1] += x0[i]
        x1[i+1] *= srδt*σ1
        x1[i+1] += x1[i]
      end
      srlt  = sqrt(fdt)
      x0[l] *= srlt*σ0
      x0[l] += x0[l-1]
      x1[l] *= srlt*σ1
      x1[l] += x1[l-1]

      # make bridge
      ite = 1.0/(Float64(l-2) * δt + fdt)
      x0df = (x0[l] - x0f)
      x1df = (x1[l] - x1f)

      for i = Base.OneTo(l-1)
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
    bm(xa  ::Float64,
       α   ::Float64,
       σ   ::Float64,
       δt  ::Float64,
       fdt ::Float64,
       srδt::Float64,
       nt  ::Int64)

Returns a Brownian motion vector starting in `xa`, with diffusion rate
`σ` and times `t`.
"""
@inline function bm(xa  ::Float64,
                    α   ::Float64,
                    σ   ::Float64,
                    δt  ::Float64,
                    fdt ::Float64,
                    srδt::Float64,
                    n   ::Int64)

  @inbounds begin
    l = n + 2
    x = randn(l)
    # for standard δt
    x[1] = xa
    if n > 0
      @turbo for i in Base.OneTo(n)
        x[i+1] *= srδt*σ
        x[i+1] += α*δt
      end
    end
    x[l] *= sqrt(fdt)*σ
    x[l] += α*fdt
    cumsum!(x, x)
  end

  return x
end




"""
    bb(xi  ::Float64,
       xf  ::Float64,
       σ   ::Float64,
       δt  ::Float64,
       fdt ::Float64,
       srδt::Float64,
       n   ::Int64)

Brownian bridge simulation.
"""
@inline function bb(xi  ::Float64,
                    xf  ::Float64,
                    σ   ::Float64,
                    δt  ::Float64,
                    fdt ::Float64,
                    srδt::Float64,
                    n   ::Int64)

  @inbounds begin

    l = n + 2
    x = randn(l)

    # for standard δt
    x[1] = xi
    if l > 2
      for i = Base.OneTo(l-2)
        x[i+1] *= srδt*σ
        x[i+1] += x[i]
      end
      x[l] *= sqrt(fdt)*σ
      x[l] += x[l-1]

      # make bridge
      ite = 1.0/(Float64(l-2) * δt + fdt)
      xdf = (x[l] - xf)

      @turbo for i = Base.OneTo(l-1)
        x[i] -= (Float64(i-1) * δt * ite * xdf)
      end
    end
    # for last non-standard δt
    x[l] = xf
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
function duoprop(x1::Float64,
                 x2::Float64,
                 t1::Float64,
                 t2::Float64,
                 σ  ::Float64)

  it = 1.0/(t1 + t2)
  return rnorm((t2*x1 + t1*x2) * it, sqrt(t1*t2*it)*σ)
end




"""
    trioprop(xp::Float64,
             x1::Float64,
             x2::Float64,
             tp::Float64,
             t1::Float64,
             t2::Float64,
             σ ::Float64)

Proposal for a trio of Gaussians.
"""
function trioprop(xp::Float64,
                  x1::Float64,
                  x2::Float64,
                  tp::Float64,
                  t1::Float64,
                  t2::Float64,
                  σ ::Float64)

    t = 1.0/(1.0/tp + 1.0/t1 + 1.0/t2)
    return rnorm((xp/tp + x1/t1 + x2/t2)*t, sqrt(t)*σ)
end




"""
    duodnorm(x ::Float64,
             x1::Float64,
             x2::Float64,
             t1::Float64,
             t2::Float64,
             σ  ::Float64)

Likelihood for a duo of Gaussians.
"""
function duodnorm(x  ::Float64,
                  x1::Float64,
                  x2::Float64,
                  t1::Float64,
                  t2::Float64,
                  σ  ::Float64)

  it = 1.0/(t1 + t2)
  return dnorm_bm(x, (t2*x1 + t1*x2)*it, sqrt(t1*t2*it)*σ)
end




"""
    duoldnorm(x  ::Float64,
              x1::Float64,
              x2::Float64,
              t1::Float64,
              t2::Float64,
              σ  ::Float64)

Likelihood for a duo of Gaussians.
"""
function duoldnorm(x  ::Float64,
                   x1::Float64,
                   x2::Float64,
                   t1::Float64,
                   t2::Float64,
                   σ  ::Float64)

  it = 1.0/(t1 + t2)
  return ldnorm_bm(x, (t2*x1 + t1*x2)*it, sqrt(t1*t2*it)*σ)
end




"""
    trioldnorm(x  ::Float64,
               xp::Float64,
               x1::Float64,
               x2::Float64,
               tp::Float64,
               t1::Float64,
               t2::Float64,
               σ  ::Float64)

Likelihood for a trio of Gaussians.
"""
function trioldnorm(x  ::Float64,
                    xp::Float64,
                    x1::Float64,
                    x2::Float64,
                    tp::Float64,
                    t1::Float64,
                    t2::Float64,
                    σ  ::Float64)

    t = 1.0/(1.0/tp + 1.0/t1 + 1.0/t2)
    return ldnorm_bm(x, (xp/tp + x1/t1 + x2/t2)*t, sqrt(t)*σ)
end




"""
    ldnorm_bm(x::Float64, μ::Float64, σsrt::Float64)

Compute the **Normal** density in logarithmic scale with
mean `μ` and standard density `σsrt` for `x`.
"""
ldnorm_bm(x::Float64, μ::Float64, σsrt::Float64) =
  -0.5*log(2.0π) - log(σsrt) - 0.5*((x - μ)/σsrt)^2




"""
    dnorm_bm(x::Float64, μ::Float64, σsrt::Float64)

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
    _sss(tree::T,
         α   ::Float64,
         f   ::Function,
         ss  ::Float64,
         n   ::Float64) where {T <: iTree}

Returns the standardized sum of squares of a diffusion with drift `α`.
"""
function _sss(tree::T,
              α   ::Float64,
              f   ::Function,
              ss  ::Float64,
              n   ::Float64) where {T <: iTree}

  ss0, n0 = _sss_b(f(tree), α, dt(tree), fdt(tree))

  ss += ss0
  n  += n0

  if def1(tree)
    ss, n = _sss(tree.d1, α, f, ss, n)
    if def2(tree)
      ss, n = _sss(tree.d2, α, f, ss, n)
    end
  end

  return ss, n
end




"""
    _sss_b(v::Array{Float64,1},
           α  ::Float64,
           δt ::Float64,
           fdt::Float64)

Returns the standardized sum of squares of a diffusion with drift `α`.
"""
function _sss_b(v::Array{Float64,1},
                α  ::Float64,
                δt ::Float64,
                fdt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    n = lastindex(v)-2

    ss = 0.0
    if n > 0
      @turbo for i in Base.OneTo(n)
        ss  += (v[i+1] - v[i] - α*δt)^2
      end
        # standardize
      ss *= 1.0/(2.0*δt)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      ss += (v[n+2] - v[n+1] - α*fdt)^2/(2.0*fdt)
      n = Float64(n + 1)
    else
      n = Float64(n)
    end
  end

  return ss, n
end




"""
    _sss(tree::T,
         f   ::Function,
         ss  ::Float64,
         n   ::Float64) where {T <: iTree}

Returns the standardized sum of squares of a diffusion without drift.
"""
function _sss(tree::T,
              f   ::Function,
              ss  ::Float64,
              n   ::Float64) where {T <: iTree}

  ss0, n0 = _sss_b(f(tree), dt(tree), fdt(tree))

  ss += ss0
  n  += n0

  if def1(tree)
    ss, n = _sss(tree.d1, f, ss, n)
    if def2(tree)
      ss, n = _sss(tree.d2, f, ss, n)
    end
  end

  return ss, n
end




"""
    _sss_b(v  ::Array{Float64,1},
           δt ::Float64,
           fdt::Float64)

Returns the standardized sum of squares of a diffusion without drift.
"""
function _sss_b(v  ::Array{Float64,1},
                δt ::Float64,
                fdt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    n = lastindex(v)-2

    ss = 0.0
    if n > 0
      @turbo for i in Base.OneTo(n)
        ss += (v[i+1] - v[i])^2
      end
      # standardize
      ss *= 1.0/(2.0*δt)
    end

    nF = Float64(n)
    # add final non-standard `δt`
    if fdt > 0.0
      ss += (v[n+2] - v[n+1])^2/(2.0*fdt)
      nF  += 1.0
    end
  end

  return ss, nF
end




"""
    dα(Ξ::Vector{T}) where {T <: iT}

Returns the overall drift of a diffusion for `α` proposal.
"""
function dα(Ξ::Vector{T}, f::Function) where {T <: iTree}

  da = 0.0
  for ξi in Ξ
    da += _dα(ξi, f)
  end

  return da
end




"""
    _dα(tree::T) where {T <: iT}

Returns the log-likelihood ratio for a `iTpb` according
to GBM birth-death for a `α` proposal.
"""
function _dα(tree::T, f::Function) where {T <: iT}

  v = f(tree)

  if def1(tree)
    if def2(tree)
      fv[end] - fv[1] + _dα(tree.d1, f) + _dα(tree.d2, f)
    else
      fv[end] - fv[1] + _dα(tree.d1, f)
    end
  else
    fv[end] - fv[1]
  end
end



