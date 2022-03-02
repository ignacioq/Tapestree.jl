#=

fossilized birth-death likelihoods

Jérémy Andréoletti
Adapted from birth-death likelihoods by Ignacio Quintero Mächler

v(^-^v)

Created 16 12 2021
=#




"""
    llik_cfbd(tree::sTfbd, λ::Float64, μ::Float64, ψ::Float64)

Log-likelihood up to a constant for constant fossilized birth-death 
given a complete `iTree` recursively.
"""
function llik_cfbd(tree::sTfbd, λ::Float64, μ::Float64, ψ::Float64)
  if istip(tree)
    # fossil tips are labelled extinct despite no actual extinction event
    if isfossil(tree)
      return - e(tree)*(λ + μ + ψ) + log(ψ)
    elseif isextinct(tree)
      return - e(tree)*(λ + μ + ψ) + log(μ)
    else
      return - e(tree)*(λ + μ + ψ) 
    end
  elseif issampledancestor(tree)
    return - e(tree)*(λ + μ + ψ) + log(ψ) +
             llik_cfbd(tree.d1::sTfbd, λ, μ, ψ)
  else
    return - e(tree)*(λ + μ + ψ) + log(λ) +
             llik_cfbd(tree.d1::sTfbd, λ, μ, ψ) + 
             llik_cfbd(tree.d2::sTfbd, λ, μ, ψ)
  end
end




"""
    llik_cfbd(Ξ::Vector{sTfbd}, 
              λ::Float64, 
              μ::Float64,
              ψ::Float64)

Log-likelihood up to a constant for constant fossilized birth-death 
given a complete `iTree` for decoupled trees.
"""
function llik_cfbd(Ξ::Vector{sTfbd}, 
                   λ::Float64, 
                   μ::Float64,
                   ψ::Float64)

  ll  = 0.0
  nsa = 0.0 # number of sampled ancestors
  for ξ in Ξ
    nsa += issampledancestorfix(ξ)
    ll  += llik_cfbd(ξ, λ, μ, ψ)
  end

  ll += Float64(lastindex(Ξ) - nsa - 1) * 0.5 * log(λ)

  return ll
end



