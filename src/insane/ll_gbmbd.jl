#=

`gbmbd` likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(Ξ   ::Vector{iTbd},
             idf ::Vector{iBffs},
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTbd` according to `gbm-bd`.
"""
function llik_gbm(Ξ   ::Vector{iTbd},
                  idf ::Vector{iBffs},
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      ll += llik_gbm(Ξ[i], α, σλ, σμ, δt, srδt)
      if d2(idf[i]) > 0
        ll += λt(Ξ[i])
      end
    end
  end

  return ll
end




"""
    llik_gbm(tree::iTbd,
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTbd` according to `gbmbd`.
"""
function llik_gbm(tree::iTbd,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  if istip(tree)
    ll_gbm_b(lλ(tree), lμ(tree), α, σλ, σμ, δt, fdt(tree), srδt,
      false, isextinct(tree))
  else
    ll_gbm_b(lλ(tree), lμ(tree), α, σλ, σμ, δt, fdt(tree), srδt, true, false) +
    llik_gbm(tree.d1, α, σλ, σμ, δt, srδt)                                    +
    llik_gbm(tree.d2, α, σλ, σμ, δt, srδt)
  end
end




"""
    ll_gbm_b(lλv ::Array{Float64,1},
             lμv ::Array{Float64,1},
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64,
             fdt ::Float64,
             srδt::Float64,
             λev ::Bool,
             μev ::Bool)

Returns the log-likelihood for a branch according to `gbmbd`.
"""
function ll_gbm_b(lλv ::Array{Float64,1},
                  lμv ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool,
                  μev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ  = 0.0
    llμ  = 0.0
    llbd = 0.0
    @turbo for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi - α*δt)^2
      llμ  += (lμvi1 - lμvi)^2
      llbd += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1))
    end

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) +
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    # add to global likelihood
    ll -= llbd*δt

    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi = lλv[nI+1]
      lμvi = lμv[nI+1]
      srfdt = sqrt(fdt)
      ll += ldnorm_bm(lλvi1, lλvi + α*fdt, srfdt*σλ)
      ll += ldnorm_bm(lμvi1, lμvi, srfdt*σμ)
      ll -= fdt*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)))
    end

    #if speciation
    if λev
      ll += lλvi1
    #if extinction
    elseif μev
      ll += lμvi1
    end

  end

  return ll
end




"""
    llik_gbm_ss(tree::iTbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ns  ::Float64,
                ne  ::Float64)

Returns the log-likelihood for a `iTbd` according to `gbmbd`.
"""
function llik_gbm_ss(tree::iTbd,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     δt  ::Float64,
                     srδt::Float64,
                     ns  ::Float64,
                     ne  ::Float64)

  if istip(tree)
    ie  = isextinct(tree)
    ne += Float64(ie)

    ll, dλ, ssλ, ssμ, nλ, irλ, irμ =
      ll_gbm_b_ss(lλ(tree), lμ(tree), α, σλ, σμ, δt, fdt(tree), srδt,
        false, ie)
  else
    ns += 1.0

    ll, dλ, ssλ, ssμ, nλ, irλ, irμ =
      ll_gbm_b_ss(lλ(tree), lμ(tree), α, σλ, σμ, δt, fdt(tree), srδt,
        true, false)

    ll1, dλ1, ssλ1, ssμ1, nλ1, irλ1, irμ1, ns, ne =
      llik_gbm_ss(tree.d1, α, σλ, σμ, δt, srδt, ns, ne)
    ll2, dλ2, ssλ2, ssμ2, nλ2, irλ2, irμ2, ns, ne =
      llik_gbm_ss(tree.d2, α, σλ, σμ, δt, srδt, ns, ne)

    ll  += ll1  + ll2
    dλ  += dλ1  + dλ2
    ssλ += ssλ1 + ssλ2
    ssμ += ssμ1 + ssμ2
    nλ  += nλ1  + nλ2
    irλ += irλ1 + irλ2
    irμ += irμ1 + irμ2
  end

  return ll, dλ, ssλ, ssμ, nλ, irλ, irμ, ns, ne
end




"""
    ll_gbm_b_ss(lλv ::Array{Float64,1},
                lμv ::Array{Float64,1},
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64,
                fdt ::Float64,
                srδt::Float64,
                λev ::Bool,
                μev ::Bool)

Returns the log-likelihood for a branch according to `gbmbd`.
"""
function ll_gbm_b_ss(lλv ::Array{Float64,1},
                     lμv ::Array{Float64,1},
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64,
                     λev ::Bool,
                     μev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ = llμ = llbdλ = llbdμ = 0.0
    @turbo for i in Base.OneTo(nI)
      lλvi   = lλv[i]
      lμvi   = lμv[i]
      lλvi1  = lλv[i+1]
      lμvi1  = lμv[i+1]
      llλ   += (lλvi1 - lλvi - α*δt)^2
      llμ   += (lμvi1 - lμvi)^2
      llbdλ += exp(0.5*(lλvi + lλvi1))
      llbdμ += exp(0.5*(lμvi + lμvi1))
    end

    # standardized sum of squares
    ssλ  = llλ/(2.0*δt)
    ssμ  = llμ/(2.0*δt)
    nλ   = Float64(nI)

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) +
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    # add to global likelihood
    irλ = llbdλ*δt
    irμ = llbdμ*δt
    ll -= (llbdλ + llbdμ)*δt

    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]
    ddλ = lλvi1 - lλv[1]

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi  = lλv[nI+1]
      lμvi  = lμv[nI+1]
      srfdt = sqrt(fdt)
      ll   += ldnorm_bm(lλvi1, lλvi + α*fdt, srfdt*σλ) +
              ldnorm_bm(lμvi1, lμvi, srfdt*σμ)
      ssλ  += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
      ssμ  += (lμvi1 - lμvi)^2/(2.0*fdt)
      nλ   += 1.0
      iriλ  = fdt*(exp(0.5*(lλvi + lλvi1)))
      iriμ  = fdt*(exp(0.5*(lμvi + lμvi1)))
      ll   -= (iriλ + iriμ)
      irλ  += iriλ 
      irμ  += iriμ
    end

    #if speciation
    if λev
      ll += lλvi1
    #if extinction
    elseif μev
      ll += lμvi1
    end
  end
  return ll, ddλ, ssλ, ssμ, nλ, irλ, irμ
end




"""
    _ss_ir_dd(tree::T,
              α   ::Float64,
              dd  ::Float64,
              ssλ ::Float64,
              ssμ ::Float64,
              n   ::Float64,
              irλ ::Float64,
              irμ ::Float64) where {T <: iTbdU}

Returns the standardized sum of squares for rate `v`, the path number `n`,
the integrated rate `ir` and the delta drift `dd`.
"""
function _ss_ir_dd(tree::T,
                   α   ::Float64,
                   dd  ::Float64,
                   ssλ ::Float64,
                   ssμ ::Float64,
                   n   ::Float64,
                   irλ ::Float64,
                   irμ ::Float64) where {T <: iTbdU}

  dd0, ssλ0, ssμ0, n0, irλ0, irμ0 = 
    _ss_ir_dd_b(lλ(tree), lμ(tree), α, dt(tree), fdt(tree))

  dd  += dd0
  ssλ += ssλ0
  ssμ += ssμ0
  n   += n0
  irλ += irλ0
  irμ += irμ0

  if def1(tree)
    dd, ssλ, ssμ, n, irλ, irμ = 
      _ss_ir_dd(tree.d1, α, dd, ssλ, ssμ, n, irλ, irμ)
    if def2(tree)
      dd, ssλ, ssμ, n, irλ, irμ = 
        _ss_ir_dd(tree.d2, α, dd, ssλ, ssμ, n, irλ, irμ)
    end
  end

  return dd, ssλ, ssμ, n, irλ, irμ
end




"""
    _ss_ir_dd_b(lλv::Array{Float64,1},
                lμv::Array{Float64,1},
                α  ::Float64,
                δt ::Float64,
                fdt::Float64)

Returns the standardized sum of squares for rate `v`, the path number `n`,
the integrated rate `ir` and the delta drift `dd`.
"""
function _ss_ir_dd_b(lλv::Array{Float64,1},
                     lμv::Array{Float64,1},
                     α  ::Float64,
                     δt ::Float64,
                     fdt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ = ssμ = irλ = irμ = 0.0
    @turbo for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      ssλ  += (lλvi1 - lλvi - α*δt)^2
      ssμ  += (lμvi1 - lμvi)^2
      irλ  += exp(0.5*(lλvi + lλvi1))
      irμ  += exp(0.5*(lμvi + lμvi1))
    end

    # standardize
    invt = 1.0/(2.0*δt)
    ssλ *= invt
    ssμ *= invt
    irλ *= δt
    irμ *= δt

    n = Float64(nI)
    # add final non-standard `δt`
    if fdt > 0.0
      invt = 1.0/(2.0*fdt)
      lλvi  = lλv[nI+1]
      lμvi  = lμv[nI+1]
      lλvi1 = lλv[nI+2]
      lμvi1 = lμv[nI+2]
      ssλ += invt * (lλvi1 - lλvi - α*fdt)^2
      ssμ += invt * (lμvi1 - lμvi)^2
      irλ += fdt*exp(0.5*(lλvi + lλvi1))
      irμ += fdt*exp(0.5*(lμvi + lμvi1))
      n += 1.0
    end
  end

  return (lλv[nI+2] - lλv[1]), ssλ, ssμ, n, irλ, irμ
end




"""
    _ss(tree::T, α::Float64, ssλ::Float64, ssμ::Float64) where {T <: iTbdU}

Returns the standardized sum of squares for the gbm part of a branch
for `bdd`.
"""
function _ss(tree::T, α::Float64, ssλ::Float64, ssμ::Float64) where {T <: iTbdU}

  ssλ0, ssμ0 = _ss_b(lλ(tree), lμ(tree), α, dt(tree), fdt(tree))
  ssλ += ssλ0
  ssμ += ssμ0

  if def1(tree)
    ssλ, ssμ = _ss(tree.d1, α, ssλ, ssμ)
    if def2(tree)
      ssλ, ssμ = _ss(tree.d2, α, ssλ, ssμ)
    end
  end

  return ssλ, ssμ
end




"""
    _ss_b(lλv::Array{Float64,1},
          lμv::Array{Float64,1},
          α  ::Float64,
          δt ::Float64,
          fdt::Float64)

Returns the standardized sum of squares for the gbm part of a branch
for `bdd`.
"""
function _ss_b(lλv::Array{Float64,1},
               lμv::Array{Float64,1},
               α  ::Float64,
               δt ::Float64,
               fdt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ = ssμ = 0.0
    @turbo for i in Base.OneTo(nI)
      ssλ  += (lλv[i+1] - lλv[i] - α*δt)^2
      ssμ  += (lμv[i+1] - lμv[i])^2
    end

    # add to global likelihood
    invt = 1.0/(2.0*δt)
    ssλ *= invt
    ssμ *= invt

    # add final non-standard `δt`
    if fdt > 0.0
      invt = 1.0/(2.0*fdt)
      ssλ += invt * (lλv[nI+2] - lλv[nI+1] - α*fdt)^2
      ssμ += invt * (lμv[nI+2] - lμv[nI+1])^2
    end
  end

  return ssλ, ssμ
end




"""
    llr_gbm_b_sep(lλp ::Array{Float64,1},
                  lμp ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  lμc ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool,
                  μev ::Bool)

Returns the log-likelihood for a branch according to `gbmbd`
separately (for gbm and bd).
"""
function llr_gbm_b_sep(lλp ::Array{Float64,1},
                       lμp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       lμc ::Array{Float64,1},
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       srδt::Float64,
                       λev ::Bool,
                       μev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλc)-2

    llrbmλ = llrbmμ = llrbdλ = llrbdμ = 0.0
    @turbo for i in Base.OneTo(nI)
      lλpi    = lλp[i]
      lλci    = lλc[i]
      lμpi    = lμp[i]
      lμci    = lμc[i]
      lλpi1   = lλp[i+1]
      lλci1   = lλc[i+1]
      lμpi1   = lμp[i+1]
      lμci1   = lμc[i+1]
      llrbmλ += (lλpi1 - lλpi - α*δt)^2 - (lλci1 - lλci - α*δt)^2
      llrbmμ += (lμpi1 - lμpi)^2 - (lμci1 - lμci)^2
      llrbdλ += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
      llrbdμ += exp(0.5*(lμpi + lμpi1)) - exp(0.5*(lμci + lμci1))
    end

    # standardized sum of squares
    ssrλ = llrbmλ/(2.0*δt)
    ssrμ = llrbmμ/(2.0*δt)

    # overall
    llrbmλ *= (-0.5/((σλ*srδt)^2))
    llrbmμ *= (-0.5/((σμ*srδt)^2))
    llrbm   = llrbmλ + llrbmμ
    llrbdλ *= (-δt)
    llrbdμ *= (-δt)

    lλpi1 = lλp[nI+2]
    lμpi1 = lμp[nI+2]
    lλci1 = lλc[nI+2]
    lμci1 = lμc[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lλpi    = lλp[nI+1]
      lλci    = lλc[nI+1]
      lμpi    = lμp[nI+1]
      lμci    = lμc[nI+1]
      ssrλ   += ((lλpi1 - lλpi - α*fdt)^2 - (lλci1 - lλci - α*fdt)^2)/(2.0*fdt)
      ssrμ   += ((lμpi1 - lμpi)^2 - (lμci1 - lμci)^2)/(2.0*fdt)
      srfdt   = sqrt(fdt)
      llrbm  += lrdnorm_bm_x(lλpi1, lλpi + α*fdt,
                             lλci1, lλci + α*fdt, srfdt*σλ) +
                lrdnorm_bm_x(lμpi1, lμpi, lμci1, lμci, srfdt*σμ)
      llrbdλ -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
      llrbdμ -= fdt*(exp(0.5*(lμpi + lμpi1)) - exp(0.5*(lμci + lμci1)))
    end
    irrλ  = -llrbdλ
    irrμ  = -llrbdμ
    llrbd = llrbdλ + llrbdμ
    if λev
      llrbd += lλpi1 - lλci1
    elseif μev
      llrbd += lμpi1 - lμci1
    end
  end

  return llrbm, llrbd, ssrλ, ssrμ, irrλ, irrμ
end



