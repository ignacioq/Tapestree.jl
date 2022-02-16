#=

pure-birth likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cpb(tree::sTpbX, λ::Float64, σx::Float64)

Log-likelihood for constant pure-birth and trait evolution
given a complete `iTree`.
"""
function llik_cpb(tree::sTpbX, λ::Float64, σx::Float64)

  el  = e(tree)
  bml = iszero(el) ? 0.0 : ldnorm_bm(xf(tree), xi(tree), sqrt(el)*σx)

  if istip(tree)
    - el * λ + bml
  else
    log(λ) - el * λ + bml           +
    llik_cpb(tree.d1::sTpbX, λ, σx) +
    llik_cpb(tree.d2::sTpbX, λ, σx)
  end
end




"""
    llik_cpb(Ξ::Vector{sTpbX}, λ::Float64, σx::Float64)

Log-likelihood up to a constant for constant pure-birth and trait evolution
given a complete `iTree` for decoupled trees.
"""
function llik_cpb(Ξ::Vector{sTpbX}, λ::Float64, σx::Float64)

  ll = 0.0
  for ξ in Ξ
    ll += llik_cpb(ξ, λ, σx)
  end

  ll += (Float64(lastindex(Ξ) - 1) * 0.5) * log(λ)

  return ll
end




"""
    λmle_cpb(tree::T)

Returns the maximum likelihood estimate for `λ` according
to a constant pure-birth process.
"""
λmle_cpb(tree::T) where {T <: iTree} = 
  Float64(ntips(tree)-2)/treelength(tree)




"""
    sdX(Ξ::Vector{T}) where {T <: iTgbm}

Returns the time standardized trait difference.
"""
function sdeltaX(Ξ::Vector{T}) where {T <: iTreeX}

  sdX = 0.0
  nX  = 0.0
  for ξ in Ξ
    sdX, nX = _sdeltaX(ξ, sdX, nX)
  end

  return sdX, nX
end




"""
    _sdX(tree::T, sdX::Float64, nX::Int64) where {T <: iTreeX}

Returns time standardized trait differences.
"""
function _sdeltaX(tree::T, sdX::Float64, nX::Float64) where {T <: iTreeX}

  el = e(tree)
  if !iszero(el)
    sdX += (xf(tree) - xi(tree))^2 / (2.0*el)
    nX  += 1.0
  end

  if isdefined(tree, :d1) 
    sdX, nX = _sdeltaX(tree.d1, sdX, nX)
    sdX, nX = _sdeltaX(tree.d2, sdX, nX)
  end

  return sdX, nX
end

