#=

pure-birth likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cpb(tree::sTpb, λ::Float64)

Log-likelihood up to a constant for constant pure-birth 
given a complete `iTree`.
"""
function llik_cpb(tree::sTpb, λ::Float64)

  if istip(tree)
    - e(tree) * λ
  else
    log(λ) - e(tree) * λ       + 
    llik_cpb(tree.d1::sTpb, λ) + 
    llik_cpb(tree.d2::sTpb, λ)
  end
end




"""
    llik_cpb(Ξ::Vector{sTpb}, λ::Float64)

Log-likelihood up to a constant for constant pure-birth 
given a complete `iTree` for decoupled trees.
"""
function llik_cpb(Ξ::Vector{sTpb}, λ::Float64)

  ll = 0.0
  for ξ in Ξ
    ll += llik_cpb(ξ, λ)
  end

  ll += ((Float64(lastindex(Ξ) - 1) * 0.5) - 1.0) * log(λ)

  return ll
end




"""
    λmle_cpb(tree::T)

Returns the maximum likelihood estimate for `λ` according
to a constant pure-birth process.
"""
λmle_cpb(tree::T) where {T <: iTree} = 
  Float64(ntips(tree)-2)/treelength(tree)
