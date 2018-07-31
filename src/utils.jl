#=

Utility functions for Tapestree

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
    rowind(x::Int64, nrow::Int64)

Get row indexing from matrix indexing.
"""
rowind(x::Int64, nrow::Int64) = mod1(x,nrow)




"""
    colind(x::Int64, nrow::Int64)

Get column indexing from matrix indexing
"""
colind(x::Int64, nrow::Int64) = cld(x, nrow)





"""
    vecind(row::Int64, col::Int64, nrow::Int64)

Get vector indexing from column and row.
"""
vecind(row::Int64, col::Int64, nrow::Int64) = row + nrow*(col - 1)





"""
    Pc(λi::Float64, λj::Float64, δt::Float64)

Estimate probability of collision.
"""
function Pc(λ1::Float64, λ0::Float64, δt::Float64)
  @fastmath begin
    λt = (λ1 + λ0)*δt
    er = exp(-λt)
    return 1.0 - er - λt*er
  end
end




"""
    makescalef(obj_ar::Float64)

Make scaling function given the objective acceptance rates.
"""
function makescalef(obj_ar::Float64)
  const nar::Float64 = 1.0 - obj_ar

  function f(window::Float64, rate::Float64)
    if rate > obj_ar
      window *= (1. + (rate - obj_ar) / nar)::Float64
    else
      window /= (2. - rate / obj_ar)::Float64
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
  return @fastmath (exp(log(λ) + stepsize * (grate - obj_ar)))::Float64
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
  
  β::Float64 = rand(linspace((1./(1. + η)),1))

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

  const spde = *(eigvecs(Σ),sqrt.(diagm(eigvals(Σ))))
  const ln   = size(spde,1)

  function f(pvec::Array{Float64,1})
    @fastmath (pvec .+ 
               BLAS.gemv('N', spde, randn(ln))::Array{Float64,1}
               )::Array{Float64,1}
  end

  return f
end





"""
    indmindif(x::Array{Float64,1}, val::Float64)

Return index for closest value 
in non-sorted arrays.
"""
function indmindif(x::Array{Float64,1}, val::Float64) 
  ibest  = start(eachindex(x)) 
  dxbest = abs(x[ibest]-val) 

  for j in eachindex(x) 
    dx = abs(x[j]-val) 
    if dx < dxbest 
        dxbest = dx 
        ibest  = j
    end 
  end 
  return ibest::Int64
end 





"""
    uniroot(f, approx = 1e-8, a = 0.0, b = 0.1)

Find the root of function between `0.0` and `b`.
"""
function uniroot(f; approx = 1e-8, a = 0.0, b = 0.1) 
  # choose b
  while sign(f(a)::Float64)::Float64 == sign(f(b)::Float64)::Float64
    b += 0.1
  end
  m::Float64 = (a + b)/2.0::Float64

  while abs(f(m)::Float64)::Float64 > approx
    if sign(f(a)::Float64)::Float64 == sign(f(m)::Float64)::Float64
      a = m::Float64
    else 
      b = m::Float64
    end
    m = (a + b)/2.0::Float64
  end
  return m::Float64
end 





"""
    int_λt(t::Float64, x::Array{Float64,1}, y::Array{Float64,1})

Cumulative pdf of λ(t) from `0` to `t`.
"""
function int_λt(t     ::Float64, 
                cumδts::Array{Float64,1}, 
                Δx    ::Array{Float64,1},
                λ     ::Float64,
                ω     ::Float64)

  d::Int64 = idxlessthan(cumδts, t)

  # riemman sums
  s::Float64 = 0.0
  for i in Base.OneTo(d-1)
    s += λ*exp(ω*Δx[i])*exp(-λ*exp(ω*Δx[i])*(cumδts[i]))*(cumδts[i+1] - cumδts[i])
  end

  # last piece to sum
  s += λ*exp(ω*Δx[d])*exp(-λ*exp(ω*Δx[d])*(cumδts[d]))*(t - cumδts[d])

  return s
end







