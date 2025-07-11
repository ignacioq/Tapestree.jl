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
                      ir  ::Float64,
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
                             ir  ::Float64,
                             ns  ::Float64,
                             sos ::Function)

  λi  = lλ(tree)
  irr = e(tree) * exp(λi) 
  ll  = sos(ll, - irr)
  ir  = sos(ir, irr)

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

    ll, dd, ss, ir, ns = 
      llik_cladspb_track!(td1, α, σλ, ll, dd, ss, ir, ns, sos)
    ll, dd, ss, ir, ns = 
      llik_cladspb_track!(td2, α, σλ, ll, dd, ss, ir, ns, sos)
  end

  return ll, dd, ss, ir, ns
end




"""
    _ss_ir_dd(tree::cTpb,
              f   ::Function,
              α   ::Float64,
              dd  ::Float64,
              ss  ::Float64,
              ir  ::Float64)

Returns the standardized sum of squares for rate `v`, the path number `n`,
the integrated rate `ir` and the delta drift `dd`.
"""
function _ss_ir_dd(tree::cTpb,
                   f   ::Function,
                   α   ::Float64,
                   dd  ::Float64,
                   ss  ::Float64,
                   ir  ::Float64)

  lλi = lλ(tree)
  ir += exp(lλi) * e(tree)

  if def1(tree)
    td1 = tree.d1
    dd, ss, ir = _ss_ir_dd(td1, f, α, dd, ss, ir)
    if def2(tree)
      td2 = tree.d2

      dd, ss, ir = _ss_ir_dd(td2, f, α, dd, ss, ir)

      lλ1 = lλ(td1)
      lλ2 = lλ(td2)

      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
      dd += lλ1 + lλ2 - 2.0*lλi
    end
  end

  return dd, ss, ir
end



"""
    _ss(tree::cTpb, f::Function, α::Float64, ss::Float64)

Returns the standardized sum of squares for rate `v`.
"""
function _ss(tree::cTpb, f::Function, α::Float64, ss::Float64)

  if def1(tree)
    td1 = tree.d1
    ss += _ss(td1, f, α, ss)
    if def2(tree)
      td2 = tree.d2
      ss += _ss(td2, f, α, ss)
      lλi = lλ(tree)
      lλ1 = lλ(td1)
      lλ2 = lλ(td2)
      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
    end
  end

  return ss
end



