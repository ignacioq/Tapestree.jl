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
    - pe(tree)*λ
  else
    llik_cpb(tree.d1::sTpb, λ) + llik_cpb(tree.d2::sTpb, λ) + log(λ) - pe(tree)*λ
  end
end



"""
    λmle_cpb(tree::T)

Returns the maximum likelihood estimate for `λ` according
to a constant pure-birth process.
"""
λmle_cpb(tree::T) where {T <: iTree} = 
  Float64(sntn(tree)-1)/treelength(tree)
