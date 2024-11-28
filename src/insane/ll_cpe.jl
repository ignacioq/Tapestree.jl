#=

constant punkeek likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 22 11 2024
=#




"""
    llik_cpe(Ξ  ::Vector{sTpe}, 
            idf::Vector{iBffs},
            λ  ::Float64, 
            μ  ::Float64, 
            σa ::Float64,
            σk ::Float64,
            nλ ::Float64)

Log-likelihood up to a constant for constant birth-death punctuated equilibrium
given a complete `iTree` for decoupled trees.
"""
function llik_cpe(Ξ  ::Vector{sTpe}, 
                 idf::Vector{iBffs},
                 λ  ::Float64, 
                 μ  ::Float64, 
                 σa ::Float64,
                 σk ::Float64,
                 nλ ::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi = idf[i]
      ξi = Ξ[i]

      if d2(bi) > 0
        if sh(ξi)
          ξd = Ξ[d1(bi)]
        else
          ξd = Ξ[d2(bi)]
        end
        ll += ldnorm_bm(xi(ξd), xf(tree), σk)
      end
      iszero(ei) && continue
      ll += llik_cpe(ξi, λ, μ, σa, σk)
    end
    ll += nλ * log(λ)
  end

  return ll
end




"""
    llik_cpe(tree::sTpe, λ::Float64, μ::Float64, σa::Float64, σk::Float64)

Log-likelihood up to a constant for constant birth-death punctuated equilibrium
given a complete `iTree` recursively.
"""
function llik_cpe(tree::sTpe, λ::Float64, μ::Float64, σa::Float64, σk::Float64)

  ei = e(tree)
  if istip(tree)
    - ei*(λ + μ) + (isextinct(tree) ? log(μ) : 0.0) + 
    ldnorm_bm(xf(tree), xi(tree), σa*sqrt(ei))
  else
    xfi = xf(tree)
    log(λ) - ei*(λ + μ)                                      +
    ldnorm_bm(xfi, xi(tree), σa*sqrt(ei))                    +
    ldnorm_bm(sh(tree) ? xi(tree.d1) : xi(tree.d2), xfi, σk) +
    llik_cpe(tree.d1, λ, μ, σa, σk)                           +
    llik_cpe(tree.d2, λ, μ, σa, σk)
  end
end




"""
    llik_quartet(xa ::Float64,
                 xi ::Float64,
                 xk ::Float64,
                 xad::Float64,
                 xkd::Float64,
                 ei ::Float64,
                 ea ::Float64,
                 ek ::Float64,
                 σa2::Float64,
                 σk2::Float64)

Likelihood for a `quartet` under constant punctuated equilibrium.
"""
function llik_quartet(xa ::Float64,
                      xi ::Float64,
                      xk ::Float64,
                      xad::Float64,
                      xkd::Float64,
                      ei ::Float64,
                      ea ::Float64,
                      ek ::Float64,
                      σa2::Float64,
                      σk2::Float64)

  return logdnorm(xi,  xa, ei*σa2) + # anagenetic ancestor
         logdnorm(xad, xi, ea*σa2) + # anagenetic daughter
         logdnorm(xk,  xi,    σk2) + # cladogenetic shift
         logdnorm(xkd, xk, ek*σa2)   # cladogenetic daughter
end




"""
    llik_trio(xi ::Float64,
              xk ::Float64,
              xad::Float64,
              xkd::Float64,
              ea ::Float64,
              ek ::Float64,
              σa2::Float64,
              σk2::Float64)

Likelihood for a `trio` (crown) under constant punctuated equilibrium.
"""
function llik_trio(xi ::Float64,
                   xk ::Float64,
                   xad::Float64,
                   xkd::Float64,
                   ea ::Float64,
                   ek ::Float64,
                   σa2::Float64,
                   σk2::Float64)

  return logdnorm(xad, xi, ea*σa2) + # anagenetic daughter
         logdnorm(xk,  xi,    σk2) + # cladogenetic shift
         logdnorm(xkd, xk, ek*σa2)   # cladogenetic daughter
end



