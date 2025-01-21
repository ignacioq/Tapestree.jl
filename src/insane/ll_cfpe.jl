#=

constant punkeek likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 22 11 2024
=#




"""
    llik_cfpe(Ξ  ::Vector{sTfpe}, 
              idf::Vector{iBffs},
              λ  ::Float64, 
              μ  ::Float64, 
              ψ  ::Vector{Float64},
              σa ::Float64,
              σk ::Float64,
              nλ ::Float64,
              ψts::Vector{Float64},
              bst::Vector{Float64},
              eix::Vector{Int64})

Log-likelihood up to a constant for constant birth-death punctuated equilibrium
given a complete `iTree` for decoupled trees.
"""
function llik_cfpe(Ξ  ::Vector{sTfpe}, 
                   idf::Vector{iBffs},
                   λ  ::Float64, 
                   μ  ::Float64, 
                   ψ  ::Vector{Float64},
                   σa ::Float64,
                   σk ::Float64,
                   nλ ::Float64,
                   ψts::Vector{Float64},
                   fex::Vector{Int64},
                   bst::Vector{Float64},
                   eix::Vector{Int64})

  @inbounds begin

    nep = lastindex(ψts) + 1
    ll  = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi = idf[i]
      ξi = Ξ[i]

      if d2(bi) > 0
        lξi = fixtip(ξi)
        ξd  = if sh(lξi) Ξ[d1(bi)] else Ξ[d2(bi)] end
        ll += ldnorm_bm(xi(ξd), xf(lξi), σk)
      end

      iszero(e(bi)) && continue
      ll += llik_cfpe(ξi, λ, μ, ψ, σa, σk, bst[i], ψts, eix[i], nep)
    end


    for i in Base.OneTo(nep)
      ll += Float64(fex[i]) * log(ψ[i])
    end

    ll += nλ * log(λ)
  end

  return ll
end




"""
    llik_cfpe(tree::sTfpe, 
              λ   ::Float64, 
              μ   ::Float64, 
              ψ   ::Vector{Float64},
              σa  ::Float64, 
              σk  ::Float64,
              t   ::Float64,
              ψts ::Vector{Float64},
              ix  ::Int64,
              nep ::Int64)

Log-likelihood up to a constant for constant fossil birth-death 
punctuated equilibrium given a complete `iTree` recursively.
"""
function llik_cfpe(tree::sTfpe, 
                   λ   ::Float64, 
                   μ   ::Float64, 
                   ψ   ::Vector{Float64},
                   σa  ::Float64, 
                   σk  ::Float64,
                   t   ::Float64,
                   ψts ::Vector{Float64},
                   ix  ::Int64,
                   nep ::Int64)

  @inbounds begin

    ei = e(tree)
    ll = ldnorm_bm(xf(tree), xi(tree), σa*sqrt(ei))
    # if epoch change
    while ix < nep && t - ei < ψts[ix]
      li  = t - ψts[ix]
      ll -= li*(λ + μ + ψ[ix])
      ei -= li
      t   = ψts[ix]
      ix += 1
    end

    ll -= ei*(λ + μ + ψ[ix])
    t  -= ei

    if def1(tree)
      if def2(tree)
        ll += log(λ)                                                        +
              ldnorm_bm(sh(tree) ? xi(tree.d1) : xi(tree.d2), xf(tree), σk) +
              llik_cfpe(tree.d1::sTfpe, λ, μ, ψ, σa, σk, t, ψts, ix, nep)   +
              llik_cfpe(tree.d2::sTfpe, λ, μ, ψ, σa, σk, t, ψts, ix, nep)
      else
        ll += log(ψ[ix])                                                  +
              llik_cfpe(tree.d1::sTfpe, λ, μ, ψ, σa, σk, t, ψts, ix, nep)
      end
    else
      ll += (isextinct(tree) ? log(μ)     : 0.0) +
            (isfossil(tree)  ? log(ψ[ix]) : 0.0)
    end
  end

  return ll
end




"""
    llik_cfpe_track(tree::sTfpe,
                    λ   ::Float64, 
                    μ   ::Float64, 
                    ψ   ::Vector{Float64},
                    σa2 ::Float64, 
                    σk2 ::Float64,
                    ll  ::Float64,
                    ns  ::Float64,
                    ne  ::Float64,
                    L   ::Vector{Float64},
                    sσa ::Float64, 
                    sσk ::Float64,
                    t   ::Float64,
                    ψts ::Vector{Float64},
                    ix  ::Int64,
                    nep ::Int64, 
                    sos ::Function)

Log-likelihood up to a constant for constant fossil birth-death 
punctuated equilibrium given a complete `iTree` recursively.
"""
function llik_cfpe_track!(tree::sTfpe,
                          λ   ::Float64, 
                          μ   ::Float64, 
                          ψ   ::Vector{Float64},
                          σa2 ::Float64, 
                          σk2 ::Float64,
                          ll  ::Float64,
                          ns  ::Float64,
                          ne  ::Float64,
                          L   ::Vector{Float64},
                          sσa ::Float64, 
                          sσk ::Float64,
                          t   ::Float64,
                          ψts ::Vector{Float64},
                          ix  ::Int64,
                          nep ::Int64, 
                          sos ::Function)
  @inbounds begin
    ei = e(tree)

    # anagenetic squares
    sqa = (xf(tree) - xi(tree))^2
    sσa = sos(sσa, sqa/ei)
    ll  = sos(ll, 
            -(0.5*log(6.28318530717958623199592693708837032318115234375*σa2*ei) +
            sqa/(2.0*σa2*ei)))

    # if epoch change
    while ix < nep && t - ei < ψts[ix]
      li    = t - ψts[ix]
      L[ix] = sos(L[ix], li)
      ll    = sos(ll, -(li*(λ + μ + ψ[ix])))
      ei   -= li
      t     = ψts[ix]
      ix   += 1
    end
    ll    = sos(ll, -(ei*(λ + μ + ψ[ix])))
    t    -= ei
    L[ix] = sos(L[ix], ei)

    if def1(tree)
      if def2(tree)
        xfi = xf(tree)
        sqk = ((sh(tree) ? xi(tree.d1) : xi(tree.d2)) -  xfi)^2
        ll  = sos(ll, 
               (log(λ) -
               0.5*log(6.28318530717958623199592693708837032318115234375*σk2) -
               sqk/(2.0*σk2)))
        sσk = sos(sσk, sqk)
        ns  = sos(ns, 1.0)

        ll, ns, ne, sσa, sσk = 
          llik_cfpe_track!(tree.d1, λ, μ, ψ, σa2, σk2, ll, ns, ne, 
            L, sσa, sσk, t, ψts, ix, nep, sos)
        ll, ns, ne, sσa, sσk = 
          llik_cfpe_track!(tree.d2, λ, μ, ψ, σa2, σk2, ll, ns, ne, 
            L, sσa, sσk, t, ψts, ix, nep, sos)
      else
        ll = sos(ll, log(ψ[ix]))
        ll, ns, ne, sσa, sσk = 
          llik_cfpe_track!(tree.d1, λ, μ, ψ, σa2, σk2, ll, ns, ne, 
            L, sσa, sσk, t, ψts, ix, nep, sos)
      end
    else
      if isextinct(tree)
        ll = sos(ll, log(μ))
        ne = sos(ne, 1.0)
      end
    end
  end

  return ll, ns, ne, sσa, sσk
end




"""
    ssσak(tree::sTfpe, sσa::Float64, sσk::Float64)

Estimate the anagenetic and cladogenetic sum of squared differences, 
`sσa` and `sσk`.
"""
function ssσak(tree::sTfpe, sσa::Float64, sσk::Float64)

  ei   = e(tree)
  sσa += (xf(tree) - xi(tree))^2/ei

  if def1(tree)
    sσa, sσk = ssσak(tree.d1, sσa, sσk)
    if def2(tree)
      xk   = sh(tree) ? xi(tree.d1) : xi(tree.d2)
      sσk += (xf(tree) - xk)^2
      sσa, sσk = ssσak(tree.d2, sσa, sσk)
    end
  end

  return sσa, sσk
end





