#=

pure-birth likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cb(Ξ::Vector{sTb}, λ::Float64)

Log-likelihood up to a constant for constant pure-birth
given a complete `iTree` for decoupled trees.
"""
function llik_cb(Ξ::Vector{sTb}, λ::Float64)

  ll = 0.0
  for ξ in Ξ
    ll += llik_cb(ξ, λ)
  end

  ll += (Float64(lastindex(Ξ) - 1) * 0.5) * log(λ)

  return ll
end




"""
    llik_cb(tree::sTb, λ::Float64)

Log-likelihood up to a constant for constant pure-birth
given a complete `iTree`.
"""
function llik_cb(tree::sTb, λ::Float64)

  if istip(tree)
    - e(tree) * λ
  else
    log(λ) - e(tree) * λ       +
    llik_cb(tree.d1::sTb, λ) +
    llik_cb(tree.d2::sTb, λ)
  end
end




"""
    λmle_cb(tree::T)

Returns the maximum likelihood estimate for `λ` according
to a constant pure-birth process.
"""
λmle_cb(tree::T) where {T <: iTree} =
  Float64(ntips(tree) - 2)/treelength(tree)
