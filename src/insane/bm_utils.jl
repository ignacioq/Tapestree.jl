#=

Brownian motion utilities

Ignacio Quintero Mächler

t(-_-t)

Created 10 09 2020
=#




"""
    ll_bm(t  ::Array{Float64,1},
             lλv::Array{Float64,1},
             σ²λ::Float64, 
             δt  ::Float64)

Returns the log-likelihood for a brownian motion.
"""
function ll_bm(t ::Array{Float64,1},
               x ::Array{Float64,1},
               σ²::Float64, 
               δt::Float64)
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
    ll *= (-1.0/(2.0*σ²*δt))
    ll -= 0.5*Float64(nI)*log(σ²*δt)

    # add final non-standard `δt`
    ll += logdnorm_tc(x[nI+2], x[nI+1], (t[nI+2] - t[nI+1])*σ²)
  end

  return ll
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

    # for standard δt
    x[1] = xi
    for i = Base.OneTo(l-2)
      x[i+1] = rnorm(x[i], srδt*σ)
    end

    # for last non-standard δt
    x[l] = rnorm(x[l-1], sqrt(t[l] - t[l-1])*σ)

    # make bridge
    ite = 1.0/t[l]
    xdf = (x[l] - xf)

    for i = Base.OneTo(l)
      x[i] -= (t[i] * ite * xdf)
    end
  end

  return nothing
end




"""
    sim_bm(xa::Float64, σ::Float64, tsv::Array{Float64,1})

Returns a Brownian motion vector starting in `xa`, with diffusion rate
`σ` and times `tsv`. 
"""
function sim_bm(xa::Float64, σ::Float64, tsv::Array{Float64,1})

  @inbounds begin
    x    = similar(tsv)
    x[1] = xa
    for i in Base.OneTo(lastindex(x)-1)
      x[i+1] = rnorm(x[i], sqrt(tsv[i+1] - tsv[i])*σ)
    end
  end

  return x
end



"""
    duoprop(xd1::Float64,
           xd2::Float64,
           td1::Float64, 
           td2::Float64,
           σ² ::Float64)

Proposal for a duo of Gaussians.
"""
function duoprop(xd1::Float64,
                xd2::Float64,
                td1::Float64, 
                td2::Float64,
                σ² ::Float64)
  invt = 1.0/(td1 + td2)
  return rnorm((td2 * invt * xd1 + td1 * invt * xd2),
               sqrt(td1 * td2 * invt * σ²))
end




"""
    trioprop(xpr::Float64,
            xd1::Float64,
            xd2::Float64,
            tpr::Float64, 
            td1::Float64, 
            td2::Float64,
            σ² ::Float64)

Proposal for a trio of Gaussians.
"""
function trioprop(xpr::Float64,
                 xd1::Float64,
                 xd2::Float64,
                 tpr::Float64, 
                 td1::Float64, 
                 td2::Float64,
                 σ² ::Float64)

    t = 1.0/(1.0/tpr + 1.0/td1 + 1.0/td2)
    return rnorm((xpr/tpr + xd1/td1 + xd2/td2)*t,
                 sqrt(t*σ²))
end

