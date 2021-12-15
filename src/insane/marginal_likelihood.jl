#=

Marginal Likelihood estimation given power posteriors

Ignacio Quintero Mächler

t(-_-t)

Created 01 12 2021
=#




""" 
    gss(pp::Vector{Vector{Float64}}, βs::Vector{Float64})

Return marginal likelihood estimate using generalized stepping stone.
"""
function gss(pp::Vector{Vector{Float64}}, βs::Vector{Float64})

  K = lastindex(pp)

  # exponentiate
  mxs = Vector{Float64}(undef, K)
  for k in Base.OneTo(K)
    mx = -Inf
    ppk = pp[k]
    for i in Base.OneTo(lastindex(ppk))
      mx = max(ppk[i], mx)
    end
    mxs[k] = mx
  end

  ml = 0.0
  for k in Base.OneTo(K-1)
    ml  += (βs[k+1] - βs[k])*mxs[k]
    ppk  = pp[k]
    mxk  = mxs[k]
    n    = lastindex(ppk)
    ssk  = 0.0
    dβk = (βs[k+1] - βs[k])
    for i in Base.OneTo(n)
      ssk += exp(dβk * (ppk[i] - mxk))
    end
    ml += log(ssk/Float64(n))
  end

  return ml
end




""" 
    path_sampling(pp::Vector{Vector{Float64}}, βs::Vector{Float64})

Return marginal likelihood estimate using path sampling based on the 
power posterior.
"""
function path_sampling(pp::Vector{Vector{Float64}}, βs::Vector{Float64})

  K = lastindex(pp)

  avgs = Vector{Float64}(undef, K)
  for k in Base.OneTo(K)
    avgs[k] = mean(pp[k])
  end

  ml = 0.0
  for k in Base.OneTo(K-1)
    ml += (avgs[k] + avgs[k+1]) * (βs[k+1] - βs[k]) * 0.5
  end

  return ml
end




""" 
    stepping_stone(pp::Vector{Vector{Float64}}, βs::Vector{Float64})

Return marginal likelihood estimate using path sampling based on the 
power posterior.
"""
function stepping_stone(pp::Vector{Vector{Float64}}, βs::Vector{Float64})

  K = lastindex(pp)

  ml = 0.0
  for k in Base.OneTo(K-1)
    ppk = pp[k]
    mxs = maximum(ppk)
    dβ  = βs[k+1] - βs[k]
    n   = lastindex(ppk)

    ss = 0.0
    for i in Base.OneTo(lastindex(ppk))
      ss += exp((ppk[i] - mxs) * dβ)
    end
    ss /= n
    
    ml += log(ss) + dβ * mxs
  end

  return ml
end