#=

Brownian motion utilities

Ignacio Quintero Mächler

t(-_-t)

Created 10 09 2020
=#




"""
    ll_bm(t ::Array{Float64,1},
          x ::Array{Float64,1},
          σ ::Float64, 
          srδt::Float64)

Returns the log-likelihood for Brownian motion.
"""
function ll_bm(t ::Array{Float64,1},
               x ::Array{Float64,1},
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
    llr_bm(t   ::Array{Float64,1},
           xp  ::Array{Float64,1},
           xc  ::Array{Float64,1},
           σ   ::Float64, 
           srδt::Float64)

Returns the log-likelihood ratio for Brownian motion.
"""
function llr_bm(t   ::Array{Float64,1},
                xp  ::Array{Float64,1},
                xc  ::Array{Float64,1},
                σ   ::Float64, 
                srδt::Float64)

 @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

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
    llr += lrdnorm_bm_x(xp[nI+2], xp[nI+1], xc[nI+2], xc[nI+1], 
              sqrt(t[nI+2] - t[nI+1])*σ)
  end

  return llr
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
    sim_bm(xa::Float64, σ::Float64, srδt::Float64, t::Array{Float64,1})

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

Compute the **Normal** density in logarithmic scale with 
mean `μ` and standard density `σ` for `x`.
"""
ldnorm_bm(x::Float64, μ::Float64, σsrt::Float64) =
  -0.5*log(2.0π) - log(σsrt) - 0.5*((x - μ)/σsrt)^2




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


