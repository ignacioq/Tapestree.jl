#=

MCMC utility functions for Tapestree

Ignacio Quintero Mächler

February 6 2017

t(-_-t)

=#



"""
    uniupt(p::Float64, tn::Float64)

Uniform parameter window move.
"""
uniupt(p::Float64, tn::Float64) = abs(p + (rand()-0.5) * tn)




"""
    addupt(p::Float64, tn::Float64)

Gaussian parameter window move.
"""
addupt(p::Float64, tn::Float64) = p + randn() * tn




"""
    addupt_lims(p::Float64, tn::Float64, xmin::Float64, xmax::Float64)

Gaussian parameter window move within a region of interest using rejection.
"""
function addupt_lims(p::Float64, tn::Float64, xmin::Float64, xmax::Float64)

  s = p + randn() * tn
  if !(xmin < s < xmax)
    s = p
  end

  return s
end




"""
    addupt(p::Float64, tn::Float64)

Gaussian parameter window move to vector.
"""
function addupt!(p::Vector{Float64}, tn::Vector{Float64}, j::Int64, i::Int64)

  @inbounds p[j] += randn() * tn[i]

  return nothing
end




"""
    duoupd(xd1::Float64,
           xd2::Float64,
           td1::Float64, 
           td2::Float64,
           σ²ϕ::Float64)

Duo of Gaussians parameter update.
"""
function duoupd(xd1::Float64,
                xd2::Float64,
                td1::Float64, 
                td2::Float64,
                σ²ϕ::Float64)
  invt = 1.0/(td1 + td2)
  return randn()*sqrt(td1 * td2 * invt * σ²ϕ) + 
    (td2 * invt * xd1 + td1 * invt * xd2)
end






"""
    trioupd(xpr::Float64,
            xd1::Float64,
            xd2::Float64,
            tpr::Float64, 
            td1::Float64, 
            td2::Float64,
            σ²ϕ::Float64)

Trio of Gaussians parameter update.
"""
function trioupd(xpr::Float64,
                 xd1::Float64,
                 xd2::Float64,
                 tpr::Float64, 
                 td1::Float64, 
                 td2::Float64,
                 σ²ϕ::Float64)

    t = 1.0/(1.0/tpr + 1.0/td1 + 1.0/td2)
    return randn()*sqrt(t*σ²ϕ) + (xpr/tpr + xd1/td1 + xd2/td2)*t
end






"""
    absaddupt(p::Float64, tn::Float64)

Non-negative Gaussian parameter window move.
"""
absaddupt(p::Float64, tn::Float64) = abs(p + randn() * tn)




"""
    mulupt(p::Float64, tn::Float64)

Multiplicative parameter window move.
"""
mulupt(p::Float64, tn::Float64) = p * exp((rand() - 0.5) * tn)





"""
    makescalef(obj_ar::Float64)

Make scaling function given the objective acceptance rates.
"""
function makescalef(obj_ar::Float64)
  nar::Float64 = 1.0 - obj_ar

  function f(window::Float64, rate::Float64)
    if rate > obj_ar
      window *= (1.0 + (rate - obj_ar) / nar)::Float64
    else
      window /= (2.0 - rate / obj_ar)::Float64
    end
    
    return window::Float64
  end

  return f
end




"""
    globalscalef(λ::Float64, grate::Float64, stepsize::Float64, obj_ar::Float64)

Estimate global scaling factor.
"""
function globalscalef(λ       ::Float64, 
                      grate   ::Float64, 
                      stepsize::Float64, 
                      obj_ar  ::Float64)
  return (exp(log(λ) + stepsize * (grate - obj_ar)))::Float64
end





"""
    adaptiveupd!(Σ::Array{Float64,2}, psam::Array{Float64,1}, pmean::Array{Float64,1}, stepsize::Float64)

Adaptive update for parameter mean and Σ in place.
"""
function adaptiveupd!(Σ       ::Array{Float64,2},
                      psam    ::Array{Float64,1},
                      pmean   ::Array{Float64,1},
                      stepsize::Float64)

  @inbounds begin
    for i in Base.OneTo(length(pmean))
      psam[i]  -= pmean[i]::Float64
      pmean[i] += (stepsize * psam[i])::Float64
    end

    BLAS.axpy!(stepsize,
               BLAS.gemm!('N', 'T', 1.0, psam, psam, -1.0, copy(Σ)),
               Σ)
  end
end





"""
    makestepsize(C::Float64, η::Float64)

Make function for the stepsize for the adaptive update.
"""
function makestepsize(η::Float64)
  
  β::Float64 = rand(range((1.0/(1.0 + η)),1))

  function f(t::Float64, C::Float64)
    return (C/(t^β))::Float64
  end

  return f
end




"""
    makemvnproposal(Σ::Array{Float64,2})

Make the multivariate update given the covariance matrix.
"""
function makemvnproposal(Σ::Array{Float64,2})

  spde = *(eigvecs(Σ),sqrt.(diagm(eigvals(Σ))))
  ln   = size(spde,1)

  function f(pvec::Array{Float64,1})
    (pvec .+ 
     BLAS.gemv('N', spde, randn(ln))::Array{Float64,1}
     )::Array{Float64,1}
  end

  return f
end

