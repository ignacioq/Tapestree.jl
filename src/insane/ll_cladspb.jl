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
      if d2(bi) > 0
        ll += λt(bi)
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
here
"""








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
    _ss_ir_dd(tree::T,
              f   ::Function,
              α   ::Float64,
              dd  ::Float64,
              ss  ::Float64,
              n   ::Float64,
              ir  ::Float64) where {T <: iTree}

Returns the standardized sum of squares for rate `v`, the path number `n`,
the integrated rate `ir` and the delta drift `dd`.
"""
function _ss_ir_dd(tree::T,
                   f   ::Function,
                   α   ::Float64,
                   dd  ::Float64,
                   ss  ::Float64,
                   n   ::Float64,
                   ir  ::Float64) where {T <: iTree}

  dd0, ss0, n0, ir0 = _ss_ir_dd_b(f(tree), α, dt(tree), fdt(tree))

  dd += dd0
  ss += ss0
  n  += n0
  ir += ir0

  if def1(tree)
    dd, ss, n, ir = _ss_ir_dd(tree.d1, f, α, dd, ss, n, ir)
    if def2(tree)
      dd, ss, n, ir = _ss_ir_dd(tree.d2, f, α, dd, ss, n, ir)
    end
  end

  return dd, ss, n, ir
end



"""
    _ss_ir_dd_b(v  ::Array{Float64,1},
                α  ::Float64,
                δt ::Float64,
                fdt::Float64)

Returns the standardized sum of squares for rate `v`, the path number `n`,
the integrated rate `ir` and the delta drift `dd`.
"""
function _ss_ir_dd_b(v  ::Array{Float64,1},
                     α  ::Float64,
                     δt ::Float64,
                     fdt::Float64)


    # estimate standard `δt` likelihood
    nI = lastindex(v)-2

    ss = ir = n = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        vi  = v[i]
        vi1 = v[i+1]
        ss += (vi1 - vi - α*δt)^2
        ir += exp(0.5*(vi + vi1))
      end
    
      # standardize
      ss *= 1.0/(2.0*δt)
      ir *= δt
      n  += Float64(nI)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      vi  = v[nI+1]
      vi1 = v[nI+2]
      ss += (vi1 - vi - α*fdt)^2/(2.0*fdt)
      n  += 1.0
      ir += fdt*exp(0.5*(vi + vi1))
    end

  return (v[nI+2] - v[1]), ss, n, ir
end



"""
    _ss(tree::T, f::Function, α::Float64) where {T <: iTree}

Returns the standardized sum of squares for rate `v`.
"""
function _ss(tree::T, f::Function, α::Float64) where {T <: iTree}

  ss = _ss_b(f(tree), α, dt(tree), fdt(tree))

  if def1(tree)
    ss += _ss(tree.d1, f, α)
    if def2(tree)
      ss += _ss(tree.d2, f, α)
    end
  end

  return ss
end




"""
    _ss_b(v::Array{Float64,1},
          α  ::Float64,
          δt ::Float64,
          fdt::Float64)

Returns the standardized sum of squares for rate `v`.
"""
function _ss_b(v::Array{Float64,1},
               α  ::Float64,
               δt ::Float64,
               fdt::Float64)

    # estimate standard `δt` likelihood
    nI = lastindex(v)-2

    ss = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        ss += (v[i+1] - v[i] - α*δt)^2
      end

      # standardize
      ss *= 1.0/(2.0*δt)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      ss += (v[nI+2] - v[nI+1] - α*fdt)^2/(2.0*fdt)
    end

  return ss
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



