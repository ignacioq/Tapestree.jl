#=

Clads birth-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 28 07 2025
=#




"""
    llik_clads(Ξ  ::Vector{cTbd},
               idf::Vector{iBffs},
               α  ::Float64,
               σλ ::Float64,
               σμ ::Float64)

Returns the log-likelihood for a `cTbd` according to clads.
"""
function llik_clads(Ξ  ::Vector{cTbd},
                    idf::Vector{iBffs},
                    α  ::Float64,
                    σλ ::Float64,
                    σμ ::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      ll += llik_clads(Ξ[i], α, σλ, σμ)

      bi2 = d2(bi)
      if bi2 > 0
        lλi = λt(bi)
        ξ1  = Ξ[d1(bi)]
        ξ2  = Ξ[bi2]

        ll += lλi + logdnorm2(lλ(ξ1), lλ(ξ2), lλi + α, σλ) +
                    logdnorm2(lμ(ξ1), lμ(ξ2), μt(bi),  σμ) 
      end
    end
  end

  return ll
end




"""
    llik_clads(tree::cTbd,
               α   ::Float64,
               σλ  ::Float64,
               σμ  ::Float64)

Returns the log-likelihood for a `cTbd` according to clads.
"""
function llik_clads(tree::cTbd,
                    α   ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64)

  if istip(tree)
    lμi = lμ(tree)
    - e(tree) * (exp(lλ(tree)) + exp(lμi)) + 
      (isextinct(tree) ? lμi : 0.0)
  else
    td1 = tree.d1
    td2 = tree.d2
    lλi = lλ(tree)
    lμi = lμ(tree)

    lλi - e(tree) * (exp(lλi) + exp(lμi))    +
    logdnorm2(lλ(td1), lλ(td2), lλi + α, σλ) +
    logdnorm2(lμ(td1), lμ(td2), lμi, σμ)     +
    llik_clads(td1, α, σλ, σμ)               +
    llik_clads(td2, α, σλ, σμ)
  end
end





"""
    llik_cladsbd_track!(tree::cTbd,
                        α   ::Float64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        ll  ::Float64,
                        dd  ::Float64,
                        ssλ ::Float64,
                        ssμ ::Float64,
                        ns  ::Float64,
                        ne  ::Float64,
                        sos ::Function)

Returns the log-likelihood for a `cTbd` according to clads.
"""
function llik_cladsbd_track!(tree::cTbd,
                             α   ::Float64,
                             σλ  ::Float64,
                             σμ  ::Float64,
                             ll  ::Float64,
                             dd  ::Float64,
                             ssλ ::Float64,
                             ssμ ::Float64,
                             ns  ::Float64,
                             ne  ::Float64,
                             sos ::Function)

  lλi = lλ(tree)
  μi = lμ(tree)
  ei = e(tree)
  ll = sos(ll, - ei * (exp(lλi) + exp(μi)))
 
  if def1(tree)
    ns  = sos(ns, 1.0)
    td1, td2 = tree.d1, tree.d2
    lλ1, lλ2 = lλ(td1), lλ(td2)
    sqλ = 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
    sqμ = 0.5*((lμ(td1) - μi)^2 + (lμ(td2) - μi)^2)
    ll  = sos(ll, 
              lλi - 2.0 * log(σλ) - 
              1.83787706640934533908193770912475883960723876953125 - sqλ/σλ^2 - 
              2.0 * log(σμ) - 
              1.83787706640934533908193770912475883960723876953125 - sqμ/σμ^2)
    ssλ = sos(ssλ, sqλ)
    ssμ = sos(ssμ, sqμ)
    dd  = sos(dd, lλ1 + lλ2 - 2.0*lλi)

    ll, dd, ssλ, ssμ, ns, ne = 
      llik_cladsbd_track!(td1, α, σλ, σμ, ll, dd, ssλ, ssμ, ns, ne, sos)
    ll, dd, ssλ, ssμ, ns, ne = 
      llik_cladsbd_track!(td2, α, σλ, σμ, ll, dd, ssλ, ssμ, ns, ne, sos)
  elseif isextinct(tree)
    ne = sos(ne, 1.0)
    ll = sos(ll, μi)
  end

  return ll, dd, ssλ, ssμ, ns, ne
end




"""
    _dd_ss(tree::cTbd,
           α   ::Float64,
           dd  ::Float64,
           ssλ ::Float64,
           ssμ ::Float64)

Returns the standardized sum of squares and the delta drift.
"""
function _dd_ss(tree::cTbd,
                α   ::Float64,
                dd  ::Float64,
                ssλ ::Float64,
                ssμ ::Float64)

  if def1(tree)
    td1 = tree.d1
    dd, ssλ, ssμ = _dd_ss(td1, α, dd, ssλ, ssμ)
    if def2(tree)
      td2 = tree.d2
      dd, ssλ, ssμ = _dd_ss(td2, α, dd, ssλ, ssμ)

      lλi  = lλ(tree)
      lμi  = lμ(tree)
      lλ1, lλ2  = lλ(td1), lλ(td2)
      dd  += lλ1 + lλ2 - 2.0*lλi
      ssλ += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
      ssμ += 0.5*((lμ(td1) - lμi)^2 + (lμ(td2) - lμi)^2)
    end
  end

  return dd, ssλ, ssμ
end




"""
    _ir(tree::cTbd, irλ::Float64, irμ::Float64)

Returns the the integrated speciation and extinction rate `irλ` and `irμ`.
"""
function _ir(tree::cT, irλ::Float64, irμ::Float64)

  ei   = e(tree)
  irλ += ei * exp(lλ(tree))
  irμ += ei * exp(lμ(tree))

  if def1(tree)
    irλ, irμ = _ir(tree.d1, irλ, irμ)
    if def2(tree)
      irλ, irμ = _ir(tree.d2, irλ, irμ)
    end
  end

  return irλ, irμ
end



