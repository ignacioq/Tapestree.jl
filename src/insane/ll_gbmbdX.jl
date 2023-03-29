#=

`gbmbd` likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(Ξ   ::Vector{iTbdX},
             idf ::Vector{iBffs},
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             βλ  ::Float64,
             σx  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTbdX` according to `gbm-bd`.
"""
function llik_gbm(Ξ   ::Vector{iTbdX},
                  idf ::Vector{iBffs},
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      ll += llik_gbm(Ξ[i], α, σλ, σμ, βλ, σx, δt, srδt)
      if !iszero(d1(idf[i]))
        ll += λt(Ξ[i])
      end
    end
  end

  return ll
end




"""
    llik_gbm(tree::iTbdX,
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             βλ  ::Float64,
             σx  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTbdX` according to `gbmbd`.
"""
function llik_gbm(tree::iTbdX,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  if istip(tree)
    ll_gbm_b(lλ(tree), lμ(tree), xv(tree), α, σλ, σμ, βλ, σx, 
      δt, fdt(tree), srδt, false, isextinct(tree))
  else
    ll_gbm_b(lλ(tree), lμ(tree), xv(tree), α, σλ, σμ, βλ, σx,
      δt, fdt(tree), srδt, true, false)             +
    llik_gbm(tree.d1, α, σλ, σμ, βλ, σx, δt, srδt)  +
    llik_gbm(tree.d2, α, σλ, σμ, βλ, σx, δt, srδt)
  end
end




"""
    ll_gbm_b(lλv ::Array{Float64,1},
             lμv ::Array{Float64,1},
             xv  ::Array{Float64,1},
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             βλ  ::Float64,
             σx  ::Float64,
             δt  ::Float64,
             fdt ::Float64,
             srδt::Float64,
             λev ::Bool,
             μev ::Bool)

Returns the log-likelihood for a branch according to `gbmbd`.
"""
function ll_gbm_b(lλv ::Array{Float64,1},
                  lμv ::Array{Float64,1},
                  xv  ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
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
    llx  = 0.0
    llbd = 0.0
    @avx for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      xvi   = xv[i]
      llx  += (xv[i+1] - xvi)^2
      llλ  += (lλvi1 - lλvi - (α + βλ*xvi)*δt)^2
      llμ  += (lμvi1 - lμvi)^2
      llbd += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1))
    end

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) +
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    # add to global likelihood
    ll -= llbd*δt
    ll += llx*(-0.5/((σx*srδt)^2)) - Float64(nI)*(log(σx*srδt) + 0.5*log(2.0π))

    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi = lλv[nI+1]
      lμvi = lμv[nI+1]
      xvi  = xv[nI+1]
      srfdt = sqrt(fdt)
      ll += ldnorm_bm(lλvi1, lλvi + (α + βλ*xvi)*fdt, srfdt*σλ)     +
            ldnorm_bm(lμvi1, lμvi, srfdt*σμ)                        -
            fdt*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1))) +
            ldnorm_bm(xv[nI+2], xvi,  srfdt*σx)
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
    llik_gbm_ss(tree::iTbdX,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64)

Returns the log-likelihood for a `iTbdX` according to `gbmbd`.
"""
function llik_gbm_ss(tree::iTbdX,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     βλ  ::Float64,
                     σx  ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  if istip(tree)
    ll, dλ, ssλ, ssμ, ssx, nx =
      ll_gbm_b_ss(lλ(tree), lμ(tree), xv(tree), α, σλ, σμ, βλ, σx, 
        δt, fdt(tree), srδt, false, isextinct(tree))
  else
    ll, dλ, ssλ, ssμ, ssx, nx =
      ll_gbm_b_ss(lλ(tree), lμ(tree), xv(tree), α, σλ, σμ, βλ, σx, 
        δt, fdt(tree), srδt, true, false)

    ll1, dλ1, ssλ1, ssμ1, ssx1, nx1 =
      llik_gbm_ss(tree.d1, α, σλ, σμ, βλ, σx, δt, srδt)
    ll2, dλ2, ssλ2, ssμ2, ssx2, nx2 =
      llik_gbm_ss(tree.d2, α, σλ, σμ, βλ, σx, δt, srδt)

    ll  += ll1  + ll2
    dλ  += dλ1  + dλ2
    ssλ += ssλ1 + ssλ2
    ssμ += ssμ1 + ssμ2
    ssx += ssx1 + ssx2
    nx  += nx1  + nx2
  end

  return ll, dλ, ssλ, ssμ, ssx, nx
end




"""
    ll_gbm_b_ss(lλv ::Array{Float64,1},
                lμv ::Array{Float64,1},
                x   ::Array{Float64,1},
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                fdt ::Float64,
                srδt::Float64,
                λev ::Bool,
                μev ::Bool)

Returns the log-likelihood for a branch according to `gbmbd`.
"""
function ll_gbm_b_ss(lλv ::Array{Float64,1},
                     lμv ::Array{Float64,1},
                     x   ::Array{Float64,1},
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     βλ  ::Float64,
                     σx  ::Float64,
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
    llx  = 0.0
    llbd = 0.0
    @avx for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      xvi   = x[i]
      llx  += (x[i+1] - xvi)^2
      llλ  += (lλvi1 - lλvi - (α + βλ*xvi)*δt)^2
      llμ  += (lμvi1 - lμvi)^2
      llbd += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1))
    end

    # standardized sum of squares
    ssλ = llλ/(2.0*δt)
    ssμ = llμ/(2.0*δt)
    ssx = llx/(2.0*δt)

    nx  = Float64(nI)

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) +
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    # add to global likelihood
    ll -= llbd*δt
    ll += llx*(-0.5/((σx*srδt)^2)) - Float64(nI)*(log(σx*srδt) + 0.5*log(2.0π))

    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]

    dλ = lλvi1 - lλv[1]

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi = lλv[nI+1]
      lμvi = lμv[nI+1]
      xvi  = x[nI+1]
      xvi1 = x[nI+2]
      srfdt = sqrt(fdt)
      ll  += ldnorm_bm(lλvi1, lλvi + (α + βλ*xvi)*fdt, srfdt*σλ)     +
             ldnorm_bm(lμvi1, lμvi, srfdt*σμ)                        -
             fdt*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1))) +
             ldnorm_bm(xvi1, xvi, srfdt*σx)
      ssλ += (lλvi1 - lλvi - (α + βλ*xvi)*fdt)^2/(2.0*fdt)
      ssμ += (lμvi1 - lμvi)^2/(2.0*fdt)
      ssx += (xvi1 - xvi)^2/(2.0*fdt)
      nx  += 1.0
    end

    #if speciation
    if λev
      ll += lλvi1
    #if extinction
    elseif μev
      ll += lμvi1
    end
  end
  return ll, dλ, ssλ, ssμ, ssx, nx
end




"""
    llr_gbm_b_sep(lλp ::Array{Float64,1},
                  lμp ::Array{Float64,1},
                  xp  ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  lμc ::Array{Float64,1},
                  xc  ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
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
                       xp  ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       lμc ::Array{Float64,1},
                       xc  ::Array{Float64,1},
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       βλ  ::Float64,
                       σx  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       srδt::Float64,
                       λev ::Bool,
                       μev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλc)-2

    llrλ  = 0.0
    llrμ  = 0.0
    llrx  = 0.0
    llrbd = 0.0
    @avx for i in Base.OneTo(nI)
      lλpi   = lλp[i]
      lλci   = lλc[i]
      lμpi   = lμp[i]
      lμci   = lμc[i]
      lλpi1  = lλp[i+1]
      lλci1  = lλc[i+1]
      lμpi1  = lμp[i+1]
      lμci1  = lμc[i+1]
      xpi    = xp[i]
      xci    = xc[i]
      llrx  += (xp[i+1] - xpi)^2 - (xc[i+1] - xci)^2
      llrλ  += (lλpi1 - lλpi - (α + βλ*xpi)*δt)^2 - 
               (lλci1 - lλci - (α + βλ*xci)*δt)^2
      llrμ  += (lμpi1 - lμpi)^2 - (lμci1 - lμci)^2
      llrbd += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)) +
               exp(0.5*(lμpi + lμpi1)) - exp(0.5*(lμci + lμci1))
    end

    # standardized sum of squares
    ssrλ = llrλ/(2.0*δt)
    ssrμ = llrμ/(2.0*δt)
    ssrx = llrx/(2.0*δt)

    # overall
    llrλ  *= (-0.5/((σλ*srδt)^2))
    llrμ  *= (-0.5/((σμ*srδt)^2))
    llrbm  = llrλ + llrμ + llrx * (-0.5/((σx*srδt)^2))
    llrbd *= (-δt)

    lλpi1 = lλp[nI+2]
    lμpi1 = lμp[nI+2]
    lλci1 = lλc[nI+2]
    lμci1 = lμc[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lλpi  = lλp[nI+1]
      lλci  = lλc[nI+1]
      lμpi  = lμp[nI+1]
      lμci  = lμc[nI+1]
      xpi    = xp[nI+1]
      xpi1   = xp[nI+2]
      xci    = xc[nI+1]
      xci1   = xc[nI+2]
      srfdt  = sqrt(fdt)
      λd     = (lλpi1 - lλpi - (α + βλ*xpi)*fdt)^2 - 
               (lλci1 - lλci - (α + βλ*xci)*fdt)^2
      μd     = (lμpi1 - lμpi)^2 - (lμci1 - lμci)^2
      xd     = ((xpi1 - xpi)^2 - (xci1 - xci)^2)
      ssrλ  += λd/(2.0*fdt)
      ssrμ  += μd/(2.0*fdt)
      ssrx  += xd/(2.0*fdt)
      llrbm += λd*(-0.5/((σλ*srfdt)^2)) + 
               μd*(-0.5/((σμ*srfdt)^2)) + 
               xd*(-0.5/((σx*srfdt)^2))
      llrbd -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)) +
                    exp(0.5*(lμpi + lμpi1)) - exp(0.5*(lμci + lμci1)))
    end

    if λev
      llrbd += lλpi1 - lλci1
    elseif μev
      llrbd += lμpi1 - lμci1
    end
  end

  return llrbm, llrbd, ssrλ, ssrμ, ssrx
end




"""
    _sss_gbm(tree::T,
             α   ::Float64,
             βλ  ::Float64,
             ssλ ::Float64,
             ssμ ::Float64,
             ssx ::Float64,
             n   ::Float64) where {T <: iTX}

Returns the standardized sum of squares a `iT` according
to `gbm-bd` for a `σ` proposal.
"""
function _sss_gbm(tree::T,
                  α   ::Float64,
                  βλ  ::Float64,
                  ssλ ::Float64,
                  ssμ ::Float64,
                  ssx ::Float64,
                  n   ::Float64) where {T <: iTX}

  ssλ0, ssμ0, ssx0, n0 = _sss_gbm_b(lλ(tree), lμ(tree), xv(tree), α, βλ, 
    dt(tree), fdt(tree))

  ssλ += ssλ0
  ssμ += ssμ0
  ssx += ssx0
  n   += n0

  if def1(tree)
    ssλ, ssμ, ssx, n = _sss_gbm(tree.d1, α, βλ, ssλ, ssμ, ssx, n)
    if def2(tree)
      ssλ, ssμ, ssx, n = _sss_gbm(tree.d2, α, βλ, ssλ, ssμ, ssx, n)
    end
  end

  return ssλ, ssμ, ssx, n
end




"""
    _sss_gbm_b(lλv::Array{Float64,1},
               lμv::Array{Float64,1},
               xv ::Array{Float64,1},
               α  ::Float64,
               βλ ::Float64,
               δt ::Float64,
               fdt::Float64)

Returns the standardized sum of squares for the GBM part of a branch
for `gbmbd`.
"""
function _sss_gbm_b(lλv::Array{Float64,1},
                    lμv::Array{Float64,1},
                    xv ::Array{Float64,1},
                    α  ::Float64,
                    βλ ::Float64,
                    δt ::Float64,
                    fdt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ  = 0.0
    ssμ  = 0.0
    ssx = 0.0
    @avx for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      xvi   = xv[i]
      ssλ  += (lλvi1 - lλvi - (α + βλ*xvi)*δt)^2
      ssμ  += (lμvi1 - lμvi)^2
      ssx  += (xv[i+1] - xvi)^2
    end

    # add to global likelihood
    invt = 1.0/(2.0*δt)
    ssλ *= invt
    ssμ *= invt
    ssx *= invt

    # add final non-standard `δt`
    if fdt > 0.0
      invt = 1.0/(2.0*fdt)
      xvi  = xv[nI+1]
      ssλ += invt * (lλv[nI+2] - lλv[nI+1] - (α + βλ*xvi)*fdt)^2
      ssμ += invt * (lμv[nI+2] - lμv[nI+1])^2
      ssx += (xv[nI+2] - xvi)^2/(2.0*fdt)
      n = Float64(nI + 1)
    else
      n = Float64(nI)
    end
  end
  return ssλ, ssμ, ssx, n
end



