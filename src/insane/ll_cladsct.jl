#=

Clads constant-extinction likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    llik_clads(Ξ  ::Vector{cTct},
               idf::Vector{iBffs},
               α  ::Float64,
               σλ ::Float64,
               ϵ  ::Float64)

Returns the log-likelihood for a `cTct` according to clads.
"""
function llik_clads(Ξ  ::Vector{cTct},
                    idf::Vector{iBffs},
                    α  ::Float64,
                    σλ ::Float64,
                    ϵ  ::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      ll += llik_clads(Ξ[i], α, σλ, ϵ)

      bi2 = d2(bi)
      if bi2 > 0
        lλi = λt(bi)
        lλ1 = lλ(Ξ[d1(bi)])
        lλ2 = lλ(Ξ[bi2])

        ll += lλi + logdnorm2(lλ1, lλ2, lλi + α, σλ)
      end
    end
  end

  return ll
end




"""
    llik_clads(tree::cTct,
               α   ::Float64,
               σλ  ::Float64,
               ϵ   ::Float64)

Returns the log-likelihood for a `cTct` according to clads.
"""
function llik_clads(tree::cTct,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64)

  lλi = lλ(tree)

  if istip(tree)
    - e(tree) * exp(lλi)*(1.0 + ϵ) + (isextinct(tree) ? log(ϵ) + lλi : 0.0)
  else
    td1 = tree.d1
    td2 = tree.d2
    lλi = lλ(tree)

    lλi - e(tree) * exp(lλi)*(1.0 + ϵ)       +
    logdnorm2(lλ(td1), lλ(td2), lλi + α, σλ) +
    llik_clads(td1, α, σλ, ϵ)                +
    llik_clads(td2, α, σλ, ϵ)
  end
end




"""
    llik_cladsct_track!(tree::cTct,
                        α   ::Float64,
                        σλ  ::Float64,
                        ϵ   ::Float64,
                        ll  ::Float64,
                        dd  ::Float64,
                        ss  ::Float64,
                        seλ ::Float64,
                        ns  ::Float64,
                        ne  ::Float64,
                        sos ::Function)

Returns the log-likelihood for a `cTct` according to clads.
"""
function llik_cladsct_track!(tree::cTct,
                             α   ::Float64,
                             σλ  ::Float64,
                             ϵ   ::Float64,
                             ll  ::Float64,
                             dd  ::Float64,
                             ss  ::Float64,
                             seλ ::Float64,
                             ns  ::Float64,
                             ne  ::Float64,
                             sos ::Function)

  λi  = lλ(tree)
  ei  = e(tree)
  eλ  = ei * exp(λi)
  ll  = sos(ll, - eλ * (1.0 + ϵ))
  seλ = sos(seλ, eλ)

  if def1(tree)
    td1 = tree.d1
    td2 = tree.d2
    λ1  = lλ(td1)
    λ2  = lλ(td2)

    ns = sos(ns, 1.0)
    sq = 0.5*((λ1 - λi - α)^2 + (λ2 - λi - α)^2)
    ll = sos(ll, λi - 2.0 * log(σλ) - 
                 1.83787706640934533908193770912475883960723876953125 - sq/σλ^2)
    dd = sos(dd, λ1 + λ2 - 2.0*λi)
    ss = sos(ss, sq)

    ll, dd, ss, seλ, ns, ne = 
      llik_cladsct_track!(td1, α, σλ, ϵ, ll, dd, ss, seλ, ns, ne, sos)
    ll, dd, ss, seλ, ns, ne = 
      llik_cladsct_track!(td2, α, σλ, ϵ, ll, dd, ss, seλ, ns, ne, sos)
  elseif isextinct(tree)
    ne = sos(ne, 1.0)
    ll = sos(ll, log(ϵ) + λi)
  end

  return ll, dd, ss, seλ, ns, ne
end




"""
    _dd_ss_seλ(tree::cTct,
              α   ::Float64,
              dd  ::Float64,
              ss  ::Float64,
              seλ  ::Float64)

Returns the standardized sum of squares `ss`, the delta drift `dd` and the rate
sum `seλ`.
"""
function _dd_ss_seλ(tree::cTct,
                    α    ::Float64,
                    dd   ::Float64,
                    ss   ::Float64,
                    seλ  ::Float64)

  lλi = lλ(tree)
  seλ += e(tree) * exp(lλi)

  if def1(tree)
    td1 = tree.d1
    dd, ss, seλ = _dd_ss_seλ(td1, α, dd, ss, seλ)
    if def2(tree)
      td2 = tree.d2
      dd, ss, seλ = _dd_ss_seλ(td2, α, dd, ss, seλ)

      lλ1 = lλ(td1)
      lλ2 = lλ(td2)
      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
      dd += lλ1 + lλ2 - 2.0*lλi
    end
  end

  return dd, ss, seλ
end



