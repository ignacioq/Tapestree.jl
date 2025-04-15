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
        lξi = fixtip(ξi)
        ξd  = if sh(lξi) Ξ[d1(bi)] else Ξ[d2(bi)] end
        ll += ldnorm_bm(xi(ξd), xf(lξi), σk)
      end

      iszero(e(bi)) && continue
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
    llik_cpe(tree.d1, λ, μ, σa, σk)                          +
    llik_cpe(tree.d2, λ, μ, σa, σk)
  end
end




"""
    llik_cpe_track(tree::sTpe,
                   λ   ::Float64, 
                   μ   ::Float64, 
                   σa  ::Float64, 
                   σk  ::Float64,
                   ns  ::Float64,
                   ne  ::Float64,
                   L   ::Float64,
                   sσa ::Float64, 
                   sσk ::Float64)

Log-likelihood up to a constant for constant birth-death punctuated equilibrium
given a complete `iTree` recursively.
"""
function llik_cpe_track(tree::sTpe,
                        λ   ::Float64, 
                        μ   ::Float64, 
                        σa2 ::Float64, 
                        σk2 ::Float64,
                        ll  ::Float64,
                        ns  ::Float64,
                        ne  ::Float64,
                        L   ::Float64,
                        sσa ::Float64, 
                        sσk ::Float64)

  ei = e(tree)
  L += ei
  if istip(tree)
    sqi  = (xf(tree) - xi(tree))^2
    ll  -= ei*(λ + μ)                                                        +
           0.5*log(6.28318530717958623199592693708837032318115234375*σa2*ei) + 
           sqi/(2.0*σa2*ei)
    sσa += sqi/ei
    if isextinct(tree)
      ll += log(μ)
      ne += 1.0
    end
  else
    xfi  = xf(tree)
    sqa  = (xfi - xi(tree))^2
    sqk  = ((sh(tree) ? xi(tree.d1) : xi(tree.d2)) -  xfi)^2
    ll  += log(λ) - ei*(λ + μ)                                               -
           0.5*log(6.28318530717958623199592693708837032318115234375*σa2*ei) -
           sqa/(2.0*σa2*ei)                                                  -
           0.5*log(6.28318530717958623199592693708837032318115234375*σk2)    -
           sqk/(2.0*σk2)
    sσa += sqa/ei
    sσk += sqk
    ns  += 1.0
    ll, ns, ne, L, sσa, sσk = 
      llik_cpe_track(tree.d1, λ, μ, σa2, σk2, ll, ns, ne, L, sσa, sσk)
    ll, ns, ne, L, sσa, sσk = 
      llik_cpe_track(tree.d2, λ, μ, σa2, σk2, ll, ns, ne, L, sσa, sσk)
  end

  return ll, ns, ne, L, sσa, sσk
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

Likelihood for a `trio` under constant punctuated equilibrium.
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



"""
    ssσak(Ξ::Vector{T}, idf::Vector{iBffs}) where T <: sT

Estimate the anagenetic and cladogenetic sum of squared differences, 
`sσa` and `sσk`.
"""
function ssσak(Ξ::Vector{T}, idf::Vector{iBffs}) where T <: sT

  @inbounds begin
    sσa = sσk = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi = idf[i]
      ξi = Ξ[i]

      if d2(bi) > 0
        lξi = fixtip(ξi)
        ξd  = if sh(lξi) Ξ[d1(bi)] else Ξ[d2(bi)] end
        sσk += (xi(ξd) - xf(lξi))^2
      end

      iszero(e(bi)) && continue
      sσa, sσk = ssσak(ξi, sσa, sσk)
    end
  end

  return sσa, sσk
end




"""
    ssσak(tree::sTpe, sσa::Float64, sσk::Float64)

Estimate the anagenetic and cladogenetic sum of squared differences, 
`sσa` and `sσk`.
"""
function ssσak(tree::sTpe, sσa::Float64, sσk::Float64)

  ei   = e(tree)
  sσa += (xf(tree) - xi(tree))^2/ei

  if def1(tree)
    xk   = sh(tree) ? xi(tree.d1) : xi(tree.d2)
    sσk += (xf(tree) - xk)^2
    sσa, sσk = ssσak(tree.d1, sσa, sσk)
    sσa, sσk = ssσak(tree.d2, sσa, sσk)
  end

  return sσa, sσk
end





