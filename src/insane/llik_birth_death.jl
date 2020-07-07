#=

birth-death likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#


"""
    llik_cbd(nλ        ::Int64, 
             nμ        ::Int64, 
             treeheight::Float64, 
             λ         ::Float64, 
             μ         ::Float64)

Log-likelihood for constant birth death given a complete tree up to a constant
"""
function llik_cbd(nλ        ::Int64, 
                  nμ        ::Int64, 
                  treeheight::Float64, 
                  λ         ::Float64, 
                  μ         ::Float64)

  return nλ*log(λ) + nμ*log(μ) - (λ+μ)*treeheight
end

