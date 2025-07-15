#=

GBM pure-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_clads(Ξ  ::Vector{cTpb},
               idf::Vector{iBffs},
               α  ::Float64,
               σλ ::Float64)

Returns the log-likelihood for a `cTpb` according to clads
"""
function llik_clads(Ξ  ::Vector{cTpb},
                    idf::Vector{iBffs},
                    α  ::Float64,
                    σλ ::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      ll += llik_clads(Ξ[i], α, σλ)

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
    llik_clads(tree::cTpb,
             α   ::Float64,
             σλ  ::Float64)

Returns the log-likelihood for a `cTpb` according to clads.
"""
function llik_clads(tree::cTpb,
                    α   ::Float64,
                    σλ  ::Float64)

  if istip(tree)
    - e(tree) * exp(lλ(tree))
  else
    td1 = tree.d1
    td2 = tree.d2
    lλi = lλ(tree)

    lλi - e(tree) * exp(lλi)                 +
    logdnorm2(lλ(td1), lλ(td2), lλi + α, σλ) +
    llik_clads(td1, α, σλ)                   +
    llik_clads(td2, α, σλ)
  end
end





"""
    llik_clads_track!(tree::cTpb,
                      α   ::Float64,
                      σλ  ::Float64,
                      ll  ::Float64,
                      dd  ::Float64,
                      ss  ::Float64,
                      ns  ::Float64,
                      sos ::Function)

Returns the log-likelihood for a `cTpb` according to clads.
"""
function llik_cladspb_track!(tree::cTpb,
                             α   ::Float64,
                             σλ  ::Float64,
                             ll  ::Float64,
                             dd  ::Float64,
                             ss  ::Float64,
                             ns  ::Float64,
                             sos ::Function)

  λi  = lλ(tree)
  ll  = sos(ll, - e(tree) * exp(λi))

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

    ll, dd, ss, ns = 
      llik_cladspb_track!(td1, α, σλ, ll, dd, ss, ns, sos)
    ll, dd, ss, ns = 
      llik_cladspb_track!(td2, α, σλ, ll, dd, ss, ns, sos)
  end

  return ll, dd, ss, ns
end




"""
    _ss_dd(tree::cTpb,
           f   ::Function,
           α   ::Float64,
           dd  ::Float64,
           ss  ::Float64)

Returns the standardized sum of squares for rate `v`, the path number `n`,
and the delta drift `dd`.
"""
function _ss_dd(tree::cTpb,
                f   ::Function,
                α   ::Float64,
                dd  ::Float64,
                ss  ::Float64)

  lλi = lλ(tree)

  if def1(tree)
    td1 = tree.d1
    dd, ss = _ss_dd(td1, f, α, dd, ss)
    if def2(tree)
      td2 = tree.d2
      dd, ss = _ss_dd(td2, f, α, dd, ss)

      lλ1 = lλ(td1)
      lλ2 = lλ(td2)
      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
      dd += lλ1 + lλ2 - 2.0*lλi
    end
  end

  return dd, ss
end




"""
    _ss(tree::cTpb, f::Function, α::Float64, ss::Float64)

Returns the standardized sum of squares for rate `v`.
"""
function _ss(tree::cTpb, f::Function, α::Float64, ss::Float64)

  if def1(tree)
    td1 = tree.d1
    ss = _ss(td1, f, α, ss)
    if def2(tree)
      td2 = tree.d2
      ss = _ss(td2, f, α, ss)
      lλi = lλ(tree)
      lλ1 = lλ(td1)
      lλ2 = lλ(td2)
      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
    end
  end

  return ss
end




"""
    _ir(tree::cTpb, ir::Float64)

Returns the the integrated rate `ir`.
"""
function _ir(tree::cTpb, ir::Float64)

  ir += exp(lλ(tree)) * e(tree)

  if def1(tree)
    ir = _ir(tree.d1, ir)
    if def2(tree)
      ir = _ir(tree.d2, ir)
    end
  end

  return ir
end



