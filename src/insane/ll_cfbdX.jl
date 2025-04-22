#=

fossilized birth-death with traits likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cfbd(tree::sTfbdx,
              λ   ::Float64,
              μ   ::Float64,
              ψ   ::Float64,
              σx  ::Float64)

Log-likelihood up to a constant for constant fossilized birth-death and 
trait evolution given a complete `iTree` recursively.
"""
function llik_cfbd(tree::sTfbdx,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   σx  ::Float64)

  el  = e(tree)
  bml = iszero(el) ? 0.0 : ldnorm_bm(xf(tree), xi(tree), sqrt(el)*σx)

  if istip(tree)
    if isfossil(tree)
      return - el*(λ + μ + ψ) + log(ψ) + bml
    elseif isextinct(tree)
      return - el*(λ + μ + ψ) + log(μ) + bml
    else
      return - el*(λ + μ + ψ) + bml
    end
  elseif isfossil(tree)
    return - el*(λ + μ + ψ) + log(ψ) + bml     +
             llik_cfbd(tree.d1::sTfbdx, λ, μ, ψ, σx)
  else
    return - el*(λ + μ + ψ) + log(λ) + bml      +
             llik_cfbd(tree.d1::sTfbdx, λ, μ, ψ, σx) +
             llik_cfbd(tree.d2::sTfbdx, λ, μ, ψ, σx)
  end
end




"""
    llik_cfbd(Ξ ::Vector{sTfbdx},
              λ ::Float64,
              μ ::Float64,
              ψ ::Float64,
              σx::Float64)

Log-likelihood up to a constant for constant fossilized birth-death and 
trait evolution given a complete `iTree` for decoupled trees.
"""
function llik_cfbd(Ξ ::Vector{sTfbdx},
                   λ ::Float64,
                   μ ::Float64,
                   ψ ::Float64,
                   σx::Float64)

  ll = 0.0
  nf = 0
  for ξ in Ξ
    nf += isinternalfossil(ξ)
    ll += llik_cfbd(ξ, λ, μ, ψ, σx)
  end
  ll += Float64(lastindex(Ξ) - nf - 1) * 0.5 * log(λ)

  return ll
end




"""
    _sdeltaX(tree::T, sdX::Float64, nX::Int64) where {T <: Tx}

Returns time standardized trait differences.
"""
function _sdeltaX(tree::T, sdX::Float64, nX::Float64) where {T <: sTfbdx}

  el = e(tree)
  if !iszero(el)
    sdX += (xf(tree) - xi(tree))^2 / (2.0*el)
    nX  += 1.0
  end

  if def1(tree)
    sdX, nX = _sdeltaX(tree.d1, sdX, nX)
    if def2(tree)
      sdX, nX = _sdeltaX(tree.d2, sdX, nX)
    end
  end

  return sdX, nX
end


