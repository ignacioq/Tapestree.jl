#=

constant punkeek likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 22 11 2024
=#




"""
    llik_pe(Ξ ::Vector{peT}, 
            λ ::Float64, 
            μ ::Float64, 
            σa::Float64,
            σk::Float64,
            nλ::Float64)

Log-likelihood up to a constant for constant birth-death punctuated equilibrium
given a complete `iTree` for decoupled trees.
"""
function llik_pe(Ξ ::Vector{peT}, 
                 λ ::Float64, 
                 μ ::Float64, 
                 σa::Float64,
                 σk::Float64,
                 nλ::Float64)

  ll = 0.0
  for ξ in Ξ
    ll += llik_pe(ξ, λ, μ, σa, σk)
  end

  ll += nλ * log(λ)

  return ll
end




"""
    llik_pe(tree::peT, λ::Float64, μ::Float64, σa::Float64, σk::Float64)

Log-likelihood up to a constant for constant birth-death punctuated equilibrium
given a complete `iTree` recursively.
"""
function llik_pe(tree::peT, λ::Float64, μ::Float64, σa::Float64, σk::Float64)

  ei = e(tree)
  if istip(tree)
    - ei*(λ + μ) + (isextinct(tree) ? log(μ) : 0.0) + 
    ldnorm_bm(xf(tree), xi(tree), σa*sqrt(ei))
  else
    xfi = xf(tree)
    log(λ) - ei*(λ + μ)                                      +
    ldnorm_bm(xfi, xi(tree), σa*sqrt(ei))                    +
    ldnorm_bm(sh(tree) ? xi(tree.d1) : xi(tree.d2), xfi, σk) +
    llik_pe(tree.d1, λ, μ, σa, σk)                           +
    llik_pe(tree.d2, λ, μ, σa, σk)
  end
end

