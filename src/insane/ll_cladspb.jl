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
             σλ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

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
    _ss_ir_dd(tree::cTpb,
              f   ::Function,
              α   ::Float64,
              dd  ::Float64,
              ss  ::Float64,
              n   ::Float64)

Returns the standardized sum of squares for rate `v`, the path number `n`,
the integrated rate `ir` and the delta drift `dd`.
"""
function _ss_ir_dd(tree::cTpb,
                   f   ::Function,
                   α   ::Float64,
                   dd  ::Float64,
                   ss  ::Float64,
                   n   ::Float64,
                   ir  ::Float64)
  lλi = lλ(tree)
  ir += exp(lλi) * e(tree)

  if def1(tree)
    td1 = tree.d1
    dd, ss, n, ir = _ss_ir_dd(td1, f, α, dd, ss, n, ir)
    if def2(tree)
      td2 = tree.d2

      dd, ss, n, ir = _ss_ir_dd(td2, f, α, dd, ss, n, ir)

      lλ1 = lλ(td1)
      lλ2 = lλ(td2)

      n  += 2.0
      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
      dd += lλ1 + lλ2 - 2.0*lλi
    end
  end

  return dd, ss, n, ir
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






















"""
    llik_clads_ssλ(tree::cTpb,
                   α   ::Float64,
                   σλ  ::Float64,
                   ns  ::Float64)

Returns the log-likelihood for a `cTpb` according to clads.
"""
function llik_clads_ssλ(tree::cTpb,
                        α   ::Float64,
                        σλ  ::Float64,
                        ns  ::Float64)


  if def1(tree)
    ll1, dλ1, ssλ1, nλ1, irλ1, ns = 
      llik_clads_ssλ(tree.d1, α, σλ, δt, srδt, ns)

    if def2(tree)



          ll1, dλ1, ssλ1, nλ1, irλ1, ns = 
      
      llik_clads_ssλ(tree.d1, α, σλ, δt, srδt, ns)


    end
  end

    ll, dλ, ssλ, nλ, irλ = 
      ll_clads_b_ssλ(lλ(tree), α, σλ, δt, fdt(tree), srδt, false)
  else



    td1 = tree.d1
    td2 = tree.d2
    lλi = lλ(tree)

    lλi - e(tree) * exp(lλi)                 +
    logdnorm2(lλ(td1), lλ(td2), lλi + α, σλ) +


    ns += 1.0

    ll, dλ, ssλ, nλ, irλ = 
      ll_clads_b_ssλ(lλ(tree), α, σλ, δt, fdt(tree), srδt, true)

    ll1, dλ1, ssλ1, nλ1, irλ1, ns = 
      llik_clads_ssλ(tree.d1, α, σλ, δt, srδt, ns)
    ll2, dλ2, ssλ2, nλ2, irλ2, ns = 
      llik_clads_ssλ(tree.d2, α, σλ, δt, srδt, ns)

    ll  += ll1  + ll2
    dλ  += dλ1  + dλ2
    ssλ += ssλ1 + ssλ2
    nλ  += nλ1  + nλ2
    irλ += irλ1 + irλ2
  end

  return ll, dλ, ssλ, nλ, irλ, ns
end




"""
    ll_clads_b_ssλ(lλv ::Array{Float64,1},
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64,
                 λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_clads_b_ssλ(lλv ::Array{Float64,1},
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64,
                      λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2
  nλ = Float64(nI)

  ll = llbm = llpb = ssλ = irλ = 0.0
  if nI > 0
    @turbo for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lλvi1 = lλv[i+1]
      llbm += (lλvi1 - lλvi - α*δt)^2
      llpb += exp(0.5*(lλvi + lλvi1))
    end

    # standardized sum of squares
    ssλ += llbm/(2.0*δt)
    irλ += llpb*δt

    # add to global likelihood
    ll += llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) - 
          llpb*δt
  end

  lλvi1 = lλv[nI+2]
  dλ    = lλvi1 - lλv[1]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi = lλv[nI+1]
    iri  = fdt*exp(0.5*(lλvi + lλvi1))
    ll  += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
           iri
    ssλ += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
    nλ  += 1.0
    irλ += iri
  end
  if λev
    ll += lλvi1
  end

  return ll, dλ, ssλ, nλ, irλ
end



"""
    llr_clads_b_sep(lλp ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth
separately for the Brownian motion and the pure-birth
"""
function llr_clads_b_sep(lλp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       srδt::Float64,
                       λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλp)-2

  llrbm = llrpb = ssrλ = 0.0
  if nI > 0
    @turbo for i in Base.OneTo(nI)
      lλpi   = lλp[i]
      lλci   = lλc[i]
      lλpi1  = lλp[i+1]
      lλci1  = lλc[i+1]
      llrbm += (lλpi1 - lλpi - α*δt)^2 - (lλci1 - lλci - α*δt)^2
      llrpb += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
    end

    # standardized sum of squares
    ssrλ  += llrbm/(2.0*δt)
    # add to global likelihood
    llrbm *= (-0.5/((σλ*srδt)^2))
    llrpb *= (-δt)
  end

  lλpi1 = lλp[nI+2]
  lλci1 = lλc[nI+2]

 # add final non-standard `δt`
  if fdt > 0.0
    lλpi   = lλp[nI+1]
    lλci   = lλc[nI+1]
    ssrλ  += ((lλpi1 - lλpi - α*fdt)^2 - (lλci1 - lλci - α*fdt)^2)/(2.0*fdt)
    llrbm += lrdnorm_bm_x(lλpi1, lλpi + α*fdt,
                          lλci1, lλci + α*fdt, sqrt(fdt)*σλ)
    llrpb -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
  end

  irrλ = -llrpb
  #if speciation
  if λev
    llrpb += lλpi1 - lλci1
  end

  return llrbm, llrpb, ssrλ, irrλ
end






"""
    int_rate(tree::iTree, f::Function)

Integrate rate given by `f`.
"""
function int_rate(tree::iTree, f::Function)

  if def1(tree)
    if def2(tree)
      int_rate(f(tree), dt(tree), fdt(tree)) +
      int_rate(tree.d1, f)         +
      int_rate(tree.d2, f)
    else
      int_rate(f(tree), dt(tree), fdt(tree)) +
      int_rate(tree.d1, f)
    end
  else
    int_rate(f(tree), dt(tree), fdt(tree))
  end
end




"""
    int_rate(v::Vector{Float64}, δt::Float64, fdt::Float64)

Integrate `v` rate.
"""
function int_rate(v::Vector{Float64}, δt::Float64, fdt::Float64)

  # estimate standard `δt` likelihood
  nI = lastindex(v)-2

  ir = 0.0
  if nI > 0
    @turbo for i in Base.OneTo(nI)
      ir += exp(0.5*(v[i] + v[i+1]))
    end
    ir *= δt
  end

  # add final non-standard `δt`
  if fdt > 0.0
    ir += fdt*exp(0.5*(v[nI+1] + v[nI+2]))
  end

  return ir
end



