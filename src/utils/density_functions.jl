#=

Denisity functions for Tapestree

Ignacio Quintero Mächler

t(-_-t)

September 23 2017

=#





"""
    logdexp(x::Float64, λ::Float64)

Compute the logarithmic transformation of the 
**Exponential** density with mean `λ` for `x`.
"""
logdexp(x::Float64, λ::Float64) = (log(λ) - λ * x)::Float64




"""
    logdunifU(x::Float64)

Standard uniform distribution (`[0.0,1.0]`).
"""
logdunifU(x::Float64) = 0.0





"""
    logdunif(x::Float64, a::Float64, b::Float64)

Make function to compute the logarithmic transformation of the 
**Uniform** density with lower bound `a` and upper bound `b` for `x`.
"""
function logdunif(x::Float64, a::Float64, b::Float64)
  if x < a 
    return -Inf
  elseif x <= b 
    return -log(b-a)::Float64
  else 
    return -Inf
  end
end





"""
    llrdexp_x(xp::Float64, xc::Float64, λ::Float64)

Compute the loglik ratio of the 
**Exponential** density for proposal 
`xp` given current `xc` both with mean `λ`.
"""
llrdexp_x(xp::Float64, xc::Float64, λ::Float64) = 
  λ * (xc - xp)





"""
    logdbeta(x::Float64, α::Float64, β::Float64)

Compute the logarithmic transformation of the 
**Beta** density with shape `α` and shape `β` for `x`.
"""
logdbeta(x::Float64, α::Float64, β::Float64) = 
  ((α-1.0) * log(x)                 +
  (β-1.0) * log(1.0 - x)           +
  log(gamma(α + β)/(gamma(α)*gamma(β))))





"""
    llrdbeta_x(xp::Float64, xc::Float64, α::Float64, β::Float64)

Compute the logarithmic ratio for the **Beta** density 
with shape `α` and shape `β` between `xp` and `xc`.
"""
function llrdbeta_x(xp::Float64, xc::Float64, α::Float64, β::Float64) 
  if !(0.0 < xp < 1.0)
    return -Inf
  else
    return ((α-1.0) * log(xp/xc) +
            (β-1.0) * log((1.0 - xp)/(1.0 - xc)))
  end
end




"""
    logdnorm(x::Float64, μ::Float64, σ²::Float64)
  
Compute the logarithmic transformation of the 
**Normal** density with mean `μ` and variance `σ²` for `x`.
"""
logdnorm(x::Float64, μ::Float64, σ²::Float64) = 
  -(0.5*log(2.0π) + 0.5*log(σ²) + (x - μ)^2/(2.0σ²))





"""
    logdnorm_tc(x::Float64, μ::Float64, σ²::Float64)

Compute the logarithmic transformation of the 
**Normal** density with mean `μ` and variance `σ²` for `x`, up to a constant
"""
logdnorm_tc(x::Float64, μ::Float64, σ²::Float64) =
  -0.5*log(σ²) - (x - μ)^2/(2.0σ²)::Float64





"""
    llrdnorm_ωx(x::Float64, xi::Float64, μp::Float64, μc::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `ωx` updates
"""
llrdnorm_ωx(x::Float64, xi::Float64, μp::Float64, μc::Float64, σ²::Float64) =
  (-(x - xi - μp)^2 + (x - xi - μc)^2)/(2.0σ²)





"""
    llrdnorm_σ²(x::Float64, μ::Float64, σ²p::Float64, σ²c::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `σ²` updates
"""
llrdnorm_σ²(x::Float64, μ::Float64, σ²p::Float64, σ²c::Float64) = 
  -0.5*(log(σ²p/σ²c) + (x - μ)^2*(1.0/σ²p - 1.0/σ²c))





"""
    llrdnorm_μ(x::Float64, μp::Float64, μc::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `μ` updates
"""
llrdnorm_μ(x::Float64, μp::Float64, μc::Float64, σ²::Float64) =
  ((x - μc)^2 - (x - μp)^2)/(2.0σ²)





"""
    llrdnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `x` updates
"""
llrdnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64) =
  ((xc - μ)^2 - (xp - μ)^2)/(2.0σ²)




"""
    llrdnorm_xμ(xp::Float64, xc::Float64, μp::Float64, μc::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `x` and `μ` updates
"""
llrdnorm_xμ(xp::Float64, xc::Float64, μp::Float64, μc::Float64, σ²::Float64) =
  ((xc - μc)^2 - (xp - μp)^2)/(2.0σ²)





"""
    logdtnorm(x::Float64, σ²::Float64)

Compute the log-likelihood density for the **Truncated Normal** density
with `a = -1` and `b = Inf`
"""
function logdtnorm(x::Float64, σ²::Float64)
  if x < -1.0 
    return -Inf
  else
    return (-x^2/(2.0*σ²) - 0.5*log(2.0π) - log(0.5*sqrt(σ²)) -
            log(1.0 - erf_custom(-1.0/(sqrt(2.0σ²)))))
  end
end





"""
    llrtnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Truncated Normal** density
for `x` with `a = -1`, `b = Inf`
"""
function llrdtnorm_x(xp::Float64, xc::Float64, σ²::Float64)
  if xp < -1.0 
    return -Inf
  else
    return (xc^2 - xp^2)/(2.0σ²)
  end
end






"""
    erf(x::Float64)

Compute the error function
"""
function erf_custom(x::Float64)
  z = abs(x)
  t = 1.0/(1.0+0.5*z)
  r = t*exp(-z*z-1.26551223 +
        t*(1.00002368 +
            t*(0.37409196 + 
              t*(0.09678418 +
                t*(-0.18628806 +
                  t*(0.27886807 +
                    t*(-1.13520398 +
                      t*(1.48851587 +
                        t*(-0.82215223 + 
                          t*0.17087277)))))))))

  return x >= 0.0 ? (1.0-r) : (r -1.0)
end





"""
    logdhcau(x::Float64, scl::Float64)

Compute the logarithmic transformation of the 
**Half-Cauchy** density with scale `scl` for `x`.
"""
logdhcau(x::Float64, scl::Float64) = 
  log(2.0 * scl/(π *(x * x + scl * scl)))





"""
    logdhcau1(x::Float64)
  
Compute the logarithmic transformation of the 
**Half-Cauchy** density with scale of 1 for `x`.
"""
logdhcau1(x::Float64) = 
  log(2.0/(π * (x * x + 1.)))

