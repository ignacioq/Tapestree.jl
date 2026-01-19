#=

Clads constant-extinction likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    llik_clads(Ξ  ::Vector{cTce},
               idf::Vector{iBffs},
               α  ::Float64,
               σλ ::Float64,
               μ  ::Float64)

Returns the log-likelihood for a `cTce` according to clads.
"""
function llik_clads(Ξ  ::Vector{cTce},
                    idf::Vector{iBffs},
                    α  ::Float64,
                    σλ ::Float64,
                    μ  ::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      ll += llik_clads(Ξ[i], α, σλ, μ)

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
    llik_clads(tree::cTce,
               α   ::Float64,
               σλ  ::Float64,
               μ   ::Float64)

Returns the log-likelihood for a `cTce` according to clads.
"""
function llik_clads(tree::cTce,
                    α   ::Float64,
                    σλ  ::Float64,
                    μ   ::Float64)

  if istip(tree)
    - e(tree) * (exp(lλ(tree)) + μ) + (isextinct(tree) ? log(μ) : 0.0)
  else
    td1 = tree.d1
    td2 = tree.d2
    lλi = lλ(tree)

    lλi - e(tree) * (exp(lλi) + μ)           +
    logdnorm2(lλ(td1), lλ(td2), lλi + α, σλ) +
    llik_clads(td1, α, σλ, μ)                +
    llik_clads(td2, α, σλ, μ)
  end
end





"""
    llik_cladsce_track!(tree::cTce,
                        α   ::Float64,
                        σλ  ::Float64,
                        μ   ::Float64,
                        ll  ::Float64,
                        dd  ::Float64,
                        ss  ::Float64,
                        ns  ::Float64,
                        ne  ::Float64,
                        L   ::Float64,
                        sos ::Function)

Returns the log-likelihood for a `cTce` according to clads.
"""
function llik_cladsce_track!(tree::cTce,
                             α   ::Float64,
                             σλ  ::Float64,
                             μ   ::Float64,
                             ll  ::Float64,
                             dd  ::Float64,
                             ss  ::Float64,
                             ns  ::Float64,
                             ne  ::Float64,
                             L   ::Float64,
                             sos ::Function)

  λi = lλ(tree)
  ei = e(tree)
  ll = sos(ll, - ei * (exp(λi) + μ))
  L  = sos(L, ei)

  if def1(tree)
    td1 = tree.d1
    td2 = tree.d2
    λ1  = lλ(td1)
    λ2  = lλ(td2)

    ns = sos(ns, 1.0)
    sq = 0.5*((λ1 - λi - α)^2 + (λ2 - λi - α)^2)
    ll = sos(ll, λi - 2.0 * log(σλ) - 
                 1.83787706640934533908193770912475883960723876953125 - sq/σλ^2)
    ss = sos(ss, sq)
    dd = sos(dd, λ1 + λ2 - 2.0*λi)

    ll, dd, ss, ns, ne, L = 
      llik_cladsce_track!(td1, α, σλ, μ, ll, dd, ss, ns, ne, L, sos)
    ll, dd, ss, ns, ne, L = 
      llik_cladsce_track!(td2, α, σλ, μ, ll, dd, ss, ns, ne, L, sos)
  elseif isextinct(tree)
    ne = sos(ne, 1.0)
    ll = sos(ll, log(μ))
  end

  return ll, dd, ss, ns, ne, L
end



