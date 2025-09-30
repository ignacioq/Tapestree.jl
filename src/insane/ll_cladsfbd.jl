#=

Clads fossilized birth-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 28 07 2025
=#





"""
    llik_clads(Ξ  ::Vector{cTfbd},
               idf::Vector{iBffs},
               αλ  ::Float64,
               αμ  ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               ψ   ::Vector{Float64},
               ψts ::Vector{Float64},
               bst ::Vector{Float64},
               eix ::Vector{Int64})

Returns the log-likelihood for a `cTfbd` according to fclads.
"""
function llik_clads(Ξ  ::Vector{cTfbd},
                    idf::Vector{iBffs},
                    αλ  ::Float64,
                    αμ  ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Vector{Float64},
                    ψts ::Vector{Float64},
                    bst ::Vector{Float64},
                    eix ::Vector{Int64})


  @inbounds begin
    ll = 0.0
    nep = lastindex(ψts) + 1
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      ll += llik_clads(Ξ[i], αλ, αμ, σλ, σμ, ψ, bst[i], ψts, eix[i], nep)

      bi2 = d2(bi)
      if bi2 > 0
        lλi = λt(bi)
        lμi = μt(bi)
        ξ1  = Ξ[d1(bi)]
        ξ2  = Ξ[bi2]

        ll += lλi + logdnorm2(lλ(ξ1), lλ(ξ2), lλi + αλ, σλ) +
                    logdnorm2(lμ(ξ1), lμ(ξ2), lμi + αμ, σμ) 
      end
    end
  end

  return ll
end




"""
    llik_clads(tree::cTfbd,
               αλ  ::Float64,
               αμ  ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               ψ   ::Vector{Float64},
               t   ::Float64,
               ψts ::Vector{Float64},
               ix  ::Int64,
               nep ::Int64)
 
Returns the log-likelihood for a `cTfbd` according to clads.
"""
function llik_clads(tree::cTfbd,
                    αλ  ::Float64,
                    αμ  ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Vector{Float64},
                    t   ::Float64,
                    ψts ::Vector{Float64},
                    ix  ::Int64,
                    nep ::Int64)
  @inbounds begin

    ei = e(tree)
    ll = 0.0

    lλi, lμi = lλ(tree), lμ(tree)
    λi, μi   = exp(lλi), exp(lμi)

    # if epoch change
    while ix < nep && t - ei < ψts[ix]
      li  = t - ψts[ix]
      ll -= li*(λi + μi + ψ[ix])
      ei -= li
      t   = ψts[ix]
      ix += 1
    end

    ll -= ei*(λi + μi + ψ[ix])
    t  -= ei

    if def1(tree)
      if def2(tree)
        td1 = tree.d1
        td2 = tree.d2

        ll += lλi                                                 +
              logdnorm2(lλ(td1), lλ(td2), lλi + αλ, σλ)           +
              logdnorm2(lμ(td1), lμ(td2), lμi + αμ, σμ)           +
              llik_clads(td1, αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep) +
              llik_clads(td2, αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep)
      else
        ll += log(ψ[ix])                                          +
              llik_clads(tree.d1, αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep)
      end
    else
      ll += (isextinct(tree) ? lμi        : 0.0) +
            (isfossil(tree)  ? log(ψ[ix]) : 0.0)
    end
  end

  return ll
end




"""
    llik_cladsfbd_track!(tree::cTfbd,
                         αλ  ::Float64,
                         αμ  ::Float64,
                         σλ  ::Float64,
                         σμ  ::Float64,
                         ψ   ::Vector{Float64},
                         t   ::Float64,
                         ψts ::Vector{Float64},
                         ix  ::Int64,
                         nep ::Int64,
                         ll  ::Float64,
                         L   ::Vector{Float64},
                         ddλ ::Float64,
                         ddμ ::Float64,
                         ssλ ::Float64,
                         ssμ ::Float64,
                         ns  ::Float64,
                         ne  ::Float64,
                         sos ::Function)

Returns the log-likelihood for a `cTfbd` according to clads.
"""
function llik_cladsfbd_track!(tree::cTfbd,
                              αλ  ::Float64,
                              αμ  ::Float64,
                              σλ  ::Float64,
                              σμ  ::Float64,
                              ψ   ::Vector{Float64},
                              t   ::Float64,
                              ψts ::Vector{Float64},
                              ix  ::Int64,
                              nep ::Int64,
                              ll  ::Float64,
                              L   ::Vector{Float64},
                              ddλ ::Float64,
                              ddμ ::Float64,
                              ssλ ::Float64,
                              ssμ ::Float64,
                              ns  ::Float64,
                              ne  ::Float64,
                              sos ::Function)

  @inbounds begin

    ei = e(tree)
    ll = 0.0

    lλi, lμi = lλ(tree), lμ(tree)
    λi, μi   = exp(lλi), exp(lμi)

    # if epoch change
    while ix < nep && t - ei < ψts[ix]
      li    = t - ψts[ix]
      L[ix] = sos(L[ix], li)
      ll    = sos(ll, - li*(λi + μi + ψ[ix]))
      ei   -= li
      t     = ψts[ix]
      ix   += 1
    end

    ll    = sos(ll, - ei*(λi + μi + ψ[ix]))
    t    -= ei
    L[ix] = sos(L[ix], ei)

    if def1(tree)
      if def2(tree)
        ns  = sos(ns, 1.0)
        td1 = tree.d1
        td2 = tree.d2
        lλ1, lλ2  = lλ(td1), lλ(td2)
        lμ1, lμ2  = lμ(td1), lμ(td2)
        sqλ = 0.5*((lλ1 - lλi - αλ)^2 + (lλ2 - lλi - αλ)^2)
        sqμ = 0.5*((lμ1 - lμi - αμ)^2 + (lμ2 - lμi - αμ)^2)
        ll  = sos(ll, 
                  lλi - 2.0 * log(σλ) - 
                  1.83787706640934533908193770912475883960723876953125 - 
                  sqλ/σλ^2 - 
                  2.0 * log(σμ) - 
                  1.83787706640934533908193770912475883960723876953125 - 
                  sqμ/σμ^2)
        ssλ = sos(ssλ, sqλ)
        ssμ = sos(ssμ, sqμ)
        ddλ = sos(ddλ, lλ1 + lλ2 - 2.0*lλi)
        ddμ = sos(ddμ, lμ1 + lμ2 - 2.0*lμi)

        ll, ddλ, ddμ, ssλ, ssμ, ns, ne = 
          llik_cladsfbd_track!(td1, αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep, 
            ll, L, ddλ, ddμ, ssλ, ssμ, ns, ne, sos)
        ll, ddλ, ddμ, ssλ, ssμ, ns, ne = 
          llik_cladsfbd_track!(td2, αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep, 
            ll, L, ddλ, ddμ, ssλ, ssμ, ns, ne, sos)
      else
        ll = sos(ll, log(ψ[ix]))

        ll, ddλ, ddμ, ssλ, ssμ, ns, ne = 
          llik_cladsfbd_track!(tree.d1, αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep, 
            ll, L, ddλ, ddμ, ssλ, ssμ, ns, ne, sos)
      end
    else
      if isextinct(tree)
        ne = sos(ne, 1.0)
        ll = sos(ll, lμi)
      end
      if isfossil(tree)
        ll = sos(ll, log(ψ[ix]))
      end
    end
  end

  return ll, ddλ, ddμ, ssλ, ssμ, ns, ne
end




"""
    _dd_ss(tree::cTfbd,
           αλ  ::Float64,
           αμ  ::Float64,
           ddλ ::Float64,
           ddμ ::Float64,
           ssλ ::Float64,
           ssμ ::Float64)

Returns the standardized sum of squares & the delta drifts.
"""
function _dd_ss(tree::cTfbd,
                αλ  ::Float64,
                αμ  ::Float64,
                ddλ ::Float64,
                ddμ ::Float64,
                ssλ ::Float64,
                ssμ ::Float64)

  if def1(tree)
    td1 = tree.d1
    ddλ, ddμ, ssλ, ssμ = _dd_ss(td1, αλ, αμ, ddλ, ddμ, ssλ, ssμ)
    if def2(tree)
      td2 = tree.d2
      ddλ, ddμ, ssλ, ssμ = _dd_ss(td2, αλ, αμ, ddλ, ddμ, ssλ, ssμ)

      lλi, lλ1, lλ2 = lλ(tree), lλ(td1), lλ(td2) 
      lμi, lμ1, lμ2 = lμ(tree), lμ(td1), lμ(td2) 
      ddλ += lλ1 + lλ2 - 2.0*lλi
      ddμ += lμ1 + lμ2 - 2.0*lμi
      ssλ += 0.5*((lλ1 - lλi - αλ)^2 + (lλ2 - lλi - αλ)^2)
      ssμ += 0.5*((lμ1 - lμi - αμ)^2 + (lμ2 - lμi - αμ)^2)
    end
  end

  return ddλ, ddμ, ssλ, ssμ
end




"""
    _ss(tree::cTfbd, α::Float64, ss::Float64, f::Function)

Returns the standardized sum of squares for `f`. 
"""
function _ss(tree::cTfbd, α::Float64, ss::Float64, f::Function)

  if def1(tree)
    td1 = tree.d1
    ss = _ss(td1, α, ss, f)
    if def2(tree)
      td2 = tree.d2
      ss = _ss(td2, α, ss, f)
      fi = f(tree)
      f1 = f(td1)
      f2 = f(td2)
      ss += 0.5*((f1 - fi - α)^2 + (f2 - fi - α)^2)
    end
  end

  return ss
end




