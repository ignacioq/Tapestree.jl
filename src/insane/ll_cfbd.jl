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
  if def1(tree)
    if def2(tree)
      - e(tree)*(λ + μ + ψ) + log(λ)       +
        llik_cfbd(tree.d1::sTfbd, λ, μ, ψ) +
        llik_cfbd(tree.d2::sTfbd, λ, μ, ψ)
    else
      - e(tree)*(λ + μ + ψ) + log(ψ)       +
        llik_cfbd(tree.d1::sTfbd, λ, μ, ψ)
    end
  else
    - e(tree)*(λ + μ + ψ)              + 
      (isextinct(tree) ? log(μ) : 0.0) +
      (isfossil(tree)  ? log(ψ) : 0.0)
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


  ll = 0.0
  nf = 0 # number of sampled ancestors
  for ξ in Ξ
    nf += isinternalfossil(ξ)
    ll += llik_cfbd(ξ, λ, μ, ψ)
  end

  ll += Float64(lastindex(Ξ) - nf - 1) * 0.5 * log(λ)

  return ll
end



