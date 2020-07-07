#=

pure birth likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cpb(tree::iTree, λ::Float64)

Log-likelihood up to a constant for constant pure-birth 
given a complete `iTree`.
"""
function llik_cpb(tree::iTree, λ::Float64)

  if istip(tree) 
    - pe(tree)*λ
  else
    llik_cpb(tree.d1, λ) +
    llik_cpb(tree.d2, λ) + 
    log(λ) - pe(tree)*λ
  end
end

