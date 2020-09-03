#=

pure-birth likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cpb(tree::iTree, λ::Float64)

Log-likelihood up to a constant for constant pure-birth 
given a complete `iTree`.
"""
function llik_cpb(tree::T, λ::Float64) where {T <: iTree}

  if istip(tree) 
    - pe(tree)*λ
  else
    llik_cpb(tree.d1::T, λ) + llik_cpb(tree.d2::T, λ) + log(λ) - pe(tree)*λ
  end
end
