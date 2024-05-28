#=
Random variable generators

Ignacio Quintero Mächler

t(-_-t)

Created 08 06 2021
=#




"""
    sample(weights::Vector{Float64})

Sample one from weights.
"""
function sample(weights::Vector{Float64})
  U  = rand()
  s  = sum(weights)
  ss = 0.0
  @inbounds begin
    for i in Base.OneTo(lastindex(weights))
      ss += weights[i]/s
      if ss > U
        return i
      end
    end
  end
  return 0
end




"""
    sample(items::Vector{Int64}, weights::Vector{Float64})

Sample from `items` using weights.
"""
function sample(items::Vector{Int64}, weights::Vector{Float64})
  U  = rand()
  s  = sum(weights)
  ss = 0.0
  @inbounds begin
    for i in Base.OneTo(lastindex(weights))
      ss += weights[i]/s
      if ss > U
        return items[i]
      end
    end
  end
end




"""
  randinvgamma(α::Float64, β::Float64)

Generate a random sample from an **Inverse Gamma** distribution with
shape `α` and rate `β` based on Tanizaki, H. (2008) "A Simple 
Gamma Random Number Generator for Arbitrary Shape Parameters".
"""
randinvgamma(α::Float64, β::Float64) = β / randgamma(α)




"""
  randgamma(α::Float64, β::Float64)

Generate a random sample from a ** Gamma** distribution with
shape `α` and rate `β` based on Tanizaki, H. (2008) "A Simple 
Gamma Random Number Generator for Arbitrary Shape Parameters".
"""
randgamma(α::Float64, β::Float64) = randgamma(α) * 1.0/β




"""
  randgamma(α::Float64)

Generate a random sample from a ** Gamma** distribution with
shape `α` and rate `1` based on Tanizaki, H. (2008) "A Simple 
Gamma Random Number Generator for Arbitrary Shape Parameters".
"""
function randgamma(α::Float64)
  
  c1 = 0.0
  if α <= 0.4
    n = 1.0/α
  elseif α <= 4.0
    n = (1.0/α)*(1.0 + (α - 0.4)/3.6)
  else
    n = 1.0/sqrt(α)
  end
  invn = 1.0/n
  b1 = α - invn
  b2 = α + invn
  if 0.4 < α  
    c1 = b1*(log(b1) - 1.0)*0.5
  end
  c2 = b2*(log(b2) - 1.0)*0.5

  return _gamma_rej(n, b1, b2, c1, c2)
end




"""
  _gamma_rej(n::Float64, b1::Float64, b2::Float64, c1::Float64, c2::Float64)

Internal function for rejection sampling for gamma generation.
"""
function _gamma_rej(n::Float64, b1::Float64, b2::Float64, c1::Float64, c2::Float64)

  w1, w2, y = _gamma_try(n, b1, b2, c1, c2)
  while y < 0.0
    w1, w2, y = _gamma_try(n, b1, b2, c1, c2)
  end
  x = n*(w2 - w1)

  while log(y) < x
    w1, w2, y = _gamma_try(n, b1, b2, c1, c2)
    while y < 0.0
      w1, w2, y = _gamma_try(n, b1, b2, c1, c2)
    end
    x = n*(w2 - w1)
  end

  return exp(x)
end




"""
  _gamma_rej(n::Float64, b1::Float64, b2::Float64, c1::Float64, c2::Float64)

Internal function for rejection sampling for gamma generation.
"""
function _gamma_try(n::Float64, b1::Float64, b2::Float64, c1::Float64, c2::Float64)
  w1 = c1 - randexp()
  w2 = c2 - randexp()
  y  = n*(b1*w2 - b2*w1)
  return w1, w2, y
end

