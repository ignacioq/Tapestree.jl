#=

Brownian motion utilities

Ignacio Quintero Mächler

t(-_-t)

Created 10 09 2020
=#




"""
    ll_bm(x ::Array{Float64,1},
          t ::Array{Float64,1},
          σ ::Float64, 
          srδt::Float64)

Returns the log-likelihood for a brownian motion.
"""
function ll_bm(x ::Array{Float64,1},
               t ::Array{Float64,1},
               σ ::Float64, 
               srδt::Float64)

 @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

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
    ll += ldnorm_bm(x[nI+2], x[nI+1], sqrt(t[nI+2] - t[nI+1])*σ)
  end

  return ll
end




"""
    bm!(x   ::Array{Float64,1},
        xi  ::Float64,
        t   ::Array{Float64,1},
        σ   ::Float64,
        srδt::Float64)

Brownian motion simulation function for updating a branch in place.
"""
function bm!(x   ::Array{Float64,1},
             xi  ::Float64,
             t   ::Array{Float64,1},
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

    cumsum!(x, x)

    # for last non-standard δt
    x[l] = rnorm(x[l-1], sqrt(t[l] - t[l-1])*σ)
  end

  return nothing
end





"""
    bb!(x   ::Array{Float64,1},
        xi  ::Float64,
        xf  ::Float64,
        t   ::Array{Float64,1},
        σ   ::Float64,
        srδt::Float64)

Brownian bridge simulation function for updating a branch in place.
"""
function bb!(x   ::Array{Float64,1},
             xi  ::Float64,
             xf  ::Float64,
             t   ::Array{Float64,1},
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

    cumsum!(x, x)

    # for last non-standard δt
    x[l] = rnorm(x[l-1], sqrt(t[l] - t[l-1])*σ)

    # make bridge
    ite = 1.0/t[l]
    xdf = (x[l] - xf)

    @simd for i = Base.OneTo(l)
      x[i] -= (t[i] * ite * xdf)
    end
  end

  return nothing
end




"""
    sim_bm(xa::Float64, σ::Float64, srδt ::Float64, t::Array{Float64,1})

Returns a Brownian motion vector starting in `xa`, with diffusion rate
`σ` and times `t`. 
"""
function sim_bm(xa::Float64, σ::Float64, srδt::Float64, t::Array{Float64,1})

  @inbounds begin

    l = lastindex(t)
    x = randn(l)
    # for standard δt
    x[1] = xa
    @simd for i in Base.OneTo(l-2)
      x[i+1] *= srδt*σ
    end
    cumsum!(x, x)

    # for last non-standard δt
    x[l] = rnorm(x[l-1], sqrt(t[l] - t[l-1])*σ)
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
    ldnorm_bm(x::Float64, μ::Float64, σsrt::Float64)
  
Compute the logarithmic transformation of the 
**Normal** density with mean `μ` and standard density `σ` for `x`.
"""
ldnorm_bm(x::Float64, μ::Float64, σsrt::Float64) =
  -0.5*log(2.0π) - log(σsrt) - 0.5*((x - μ)/σsrt)^2




