#=

fossilized birth-death with traits likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cfbd(tree::sTfbdX,
              λ   ::Float64,
              μ   ::Float64,
              ψ   ::Float64,
              σx  ::Float64)

Log-likelihood up to a constant for constant fossilized birth-death and 
trait evolution given a complete `iTree` recursively.
"""
function llik_cfbd(tree::sTfbdX,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   σx  ::Float64)

  el  = e(tree)
  bml = iszero(el) ? 0.0 : ldnorm_bm(xf(tree), xi(tree), sqrt(el)*σx)

  if istip(tree)
    if isfossil(tree)
      return - e(tree)*(λ + μ + ψ) + log(ψ) + bml
    elseif isextinct(tree)
      return - e(tree)*(λ + μ + ψ) + log(μ) + bml
    else
      return - e(tree)*(λ + μ + ψ) + bml
    end
  elseif isfossil(tree)
    return - e(tree)*(λ + μ + ψ) + log(ψ) + bml     +
             llik_cfbd(tree.d1::sTfbdX, λ, μ, ψ, σx)
  else
    return - e(tree)*(λ + μ + ψ) + log(λ) + bml      +
             llik_cfbd(tree.d1::sTfbdX, λ, μ, ψ, σx) +
             llik_cfbd(tree.d2::sTfbdX, λ, μ, ψ, σx)
  end
end




"""
    llik_cfbd(Ξ ::Vector{sTfbdX},
              λ ::Float64,
              μ ::Float64,
              ψ ::Float64,
              σx::Float64)

Log-likelihood up to a constant for constant fossilized birth-death and 
trait evolution given a complete `iTree` for decoupled trees.
"""
function llik_cfbd(Ξ ::Vector{sTfbdX},
                   λ ::Float64,
                   μ ::Float64,
                   ψ ::Float64,
                   σx::Float64)

  ll = 0.0
  for ξ in Ξ
    ll += llik_cfbd(ξ, λ, μ, ψ, σx)
  end

  ll += Float64(lastindex(Ξ) - 1) * 0.5 * log(λ)

  return ll
end



