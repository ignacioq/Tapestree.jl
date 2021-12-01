#=

Marginal Likelihood estimation given power posteriors

Ignacio Quintero Mächler

t(-_-t)

Created 01 12 2021
=#




""" 
    path_sampling(PP::Vector{Vector{Float64}}, βs::Vector{Float64})

Return marginal likelihood estimate using path sampling based on the 
power posterior.
"""
function path_sampling(PP::Vector{Vector{Float64}}, βs::Vector{Float64})

  K = lastindex(PP)

  avgs = Vector{Float64}(undef, K)
  for k in Base.OneTo(K)
    avgs[k] = mean(PP[k])
  end

  ml = 0.0
  for k in Base.OneTo(K-1)
    ml += (avgs[k] + avgs[k+1]) * (βs[k+1] - βs[k]) * 0.5
  end

  return ml
end




""" 
    stepping_stone(PP::Vector{Vector{Float64}}, βs::Vector{Float64})

Return marginal likelihood estimate using path sampling based on the 
power posterior.
"""
function stepping_stone(PP::Vector{Vector{Float64}}, βs::Vector{Float64})

  K = lastindex(PP)

  ml = 0.0
  for k in Base.OneTo(K-1)
    ppk = PP[k]
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