#=

birth-death likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cbd(tree::iTree, λ::Float64, μ::Float64)

Log-likelihood up to a constant for constant birth-death 
given a complete `iTree`.
"""
function llik_cbd(tree::iTree, λ::Float64, μ::Float64)

  if istip(tree) 
    - pe(tree)*(λ + μ) + (isextinct(tree) ? log(μ) : 0.0)
  else
    llik_cbd(tree.d1, λ, μ) +
    llik_cbd(tree.d2, λ, μ) + 
    log(λ) - pe(tree)*(λ + μ)
  end
end





