#=

birth-death with traits likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cbd(tree::sTbd, λ::Float64, μ::Float64)

Log-likelihood up to a constant for constant birth-death and trait evolution
given a complete `iTree` recursively.
"""
function llik_cbd(tree::sTbdX, λ::Float64, μ::Float64, σx::Float64)

  el  = e(tree)
  bml = iszero(el) ? 0.0 : ldnorm_bm(xf(tree), xi(tree), sqrt(el)*σx)

  if istip(tree)
    - el*(λ + μ) + (isextinct(tree) ? log(μ) : 0.0) + bml
  else
    log(λ) - el*(λ + μ) + bml         +
    llik_cbd(tree.d1::sTbdX, λ, μ, σx) +
    llik_cbd(tree.d2::sTbdX, λ, μ, σx)
  end
end




"""
    llik_cbd(Ξ::Vector{sTbdX}, λ::Float64, μ::Float64, σx::Float64)

Log-likelihood up to a constant for constant birth-death and trait evolution
given a complete `iTree` for decoupled trees.
"""
function llik_cbd(Ξ::Vector{sTbdX}, λ::Float64, μ::Float64, σx::Float64)

  ll = 0.0
  for ξ in Ξ
    ll += llik_cbd(ξ, λ, μ, σx)
  end

  ll += Float64(lastindex(Ξ) - 1) * 0.5 * log(λ)

  return ll
end



