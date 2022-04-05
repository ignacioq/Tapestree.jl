#=

GBM pure-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    llik_gbm(tree::iTpbX,
             α   ::Float64,
             σλ  ::Float64,
             βλ  ::Float64,
             σx  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTpbX` according to GBM birth-death.
"""
function llik_gbm(tree::iTpbX,
                  α   ::Float64,
                  σλ  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  if istip(tree)
    ll_gbm_b(lλ(tree), xv(tree), α, σλ, βλ, σx, δt, fdt(tree), srδt, false)
  else
    ll_gbm_b(lλ(tree), xv(tree), α, σλ, βλ, σx, δt, fdt(tree), srδt, true) +
    llik_gbm(tree.d1::iTpbX, α, σλ, βλ, σx, δt, srδt)                      +
    llik_gbm(tree.d2::iTpbX, α, σλ, βλ, σx, δt, srδt)
  end
end




"""
    llik_gbm(Ξ   ::Vector{iTpbX},
             idf ::Vector{iBffs},
             α   ::Float64,
             σλ  ::Float64,
             βλ  ::Float64,
             σx  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTpbX` according to GBM birth-death.
"""
function llik_gbm(Ξ   ::Vector{iTpbX},
                  idf ::Vector{iBffs},
                  α   ::Float64,
                  σλ  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      ll += llik_gbm(Ξ[i], α, σλ, βλ, σx, δt, srδt)

      if !it(bi)
        ll += λt(bi)
      end
    end
  end

  return ll
end




"""
    ll_gbm_b(lλv ::Array{Float64,1},
             xv  ::Array{Float64,1},
             α   ::Float64,
             σλ  ::Float64,
             σx  ::Float64,
             δt  ::Float64,
             fdt ::Float64,
             srδt::Float64,
             λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_gbm_b(lλv ::Array{Float64,1},
                  xv  ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2

  llx  = 0.0
  llbm = 0.0
  llpb = 0.0
  @avx for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    xvi   = xv[i]
    llx  += (xv[i+1] - xvi)^2
    llbm += (lλvi1 - lλvi - (α + βλ*xvi)*δt)^2
    llpb += exp(0.5*(lλvi + lλvi1))
  end

  # add to global likelihood
  ll  = llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
  ll -= llpb*δt
  ll += llx *(-0.5/((σx*srδt)^2)) - Float64(nI)*(log(σx*srδt) + 0.5*log(2.0π))

  lλvi1 = lλv[nI+2]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi = lλv[nI+1]
    xvi  = xv[nI+1]
    ll  += ldnorm_bm(lλvi1, lλvi + (α + βλ*xvi)*fdt, sqrt(fdt)*σλ) -
           fdt*exp(0.5*(lλvi + lλvi1))                             +
           ldnorm_bm(xv[nI+2], xvi,  sqrt(fdt)*σx)
  end
  if λev
    ll += lλvi1
  end

  return ll
end




"""
    llik_gbm_ss(tree::iTpbX,
                α   ::Float64,
                σλ  ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64)

Returns the log-likelihood for a `iTpbX` according to GBM birth-death.
"""
function llik_gbm_ss(tree::iTpbX,
                     α   ::Float64,
                     σλ  ::Float64,
                     βλ  ::Float64,
                     σx  ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  if istip(tree)
    ll, dλ, ssλ, ssx, nd =
      ll_gbm_b_ss(lλ(tree), α, σλ, βλ, σx, δt, fdt(tree), srδt, false)
  else
    ll, dλ, ssλ, ssx, nd =
      ll_gbm_b_ss(lλ(tree), α, σλ, βλ, σx, δt, fdt(tree), srδt, true)

    ll1, dλ1, ssλ1, ssx1, nd1 =
      llik_gbm_ss(tree.d1::iTpbX, α, σλ, βλ, σx, δt, srδt)
    ll2, dλ2, ssλ2, ssx2, nd2 =
      llik_gbm_ss(tree.d2::iTpbX, α, σλ, βλ, σx, δt, srδt)

    ll  += ll1  + ll2
    dλ  += dλ1  + dλ2
    ssλ += ssλ1 + ssλ2
    ssx += ssx1 + ssx2
    nd  += nd1  + nd2
  end

  return ll, dλ, ssλ, ssx, nd
end




"""
    ll_gbm_b_ss(lλv ::Array{Float64,1},
                α   ::Float64,
                σλ  ::Float64,
                βλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                fdt ::Float64,
                srδt::Float64,
                λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_gbm_b_ss(lλv ::Array{Float64,1},
                     α   ::Float64,
                     σλ  ::Float64,
                     βλ  ::Float64,
                     σx  ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64,
                     λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2
  nd = Float64(nI)

  llbm = 0.0
  llpb = 0.0
  llx  = 0.0
  @avx for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    xvi   = xv[i]
    llx  += (xv[i+1] - xvi)^2
    llbm += (lλvi1 - lλvi - (α + βλ*xvi)*δt)^2
    llpb += exp(0.5*(lλvi + lλvi1))
  end

  # standardized sum of squares
  ssλ = llbm/(2.0*δt)
  ssx = llx /(2.0*δt)

  # add to global likelihood
  ll  = llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
  ll -= llpb*δt
  ll += llx *(-0.5/((σx*srδt)^2)) - Float64(nI)*(log(σx*srδt) + 0.5*log(2.0π))

  lλvi1 = lλv[nI+2]

  dλ = lλvi1 - lλv[1]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi = lλv[nI+1]
    xvi  = xv[nI+1]
    xvi1 = xv[nI+2]
    ll  += ldnorm_bm(lλvi1, lλvi + (α + βλ*xvi)*fdt, sqrt(fdt)*σλ) -
           fdt*exp(0.5*(lλvi + lλvi1))                             +
           ldnorm_bm(xvi1, xvi, sqrt(fdt)*σx)
    ssλ += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
    ssx += (xvi1 - xvi)^2/(2.0*fdt)
    nd  += 1.0
  end
  if λev
    ll += lλvi1
  end

  return ll, dλ, ssλ, ssx, nd
end





"""
    llr_xv(lxp ::Array{Float64,1},
           lxc ::Array{Float64,1},
           σx  ::Float64,
           δt  ::Float64,
           fdt ::Float64,
           srδt::Float64)

Returns the log-likelihood ratio for a branch according to GBM pure-birth
separately for the Brownian motion and the pure-birth
"""
function llr_ssx(xp  ::Array{Float64,1},
                 xc  ::Array{Float64,1},
                 σx  ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64)

  # estimate standard `δt` likelihood
  n = lastindex(xc)

  ssrx = 0.0
  @avx for i in Base.OneTo(n-2)
    ssrx += (xp[i+1] - xp[i])^2 - (xc[i+1] - xc[i])^2
  end

  # likelihood ratio
  llr  = ssrx * (-0.5/((σx*srδt)^2))

  # standardized sum of squares
  ssrx *= 1.0/(2.0*δt)

 # add final non-standard `δt`
  if fdt > 0.0
    xd    = (xp[n] - xp[n-1])^2 - (xc[n] - xc[n-1])^2
    ssrx += xd/(2.0*fdt)
    llr  += xd*(-0.5/((σx*sqrt(fdt))^2))
  end

  return llr, ssrx
end




"""
    llr_gbm_b_sep(lλp ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  xp  ::Array{Float64,1},
                  xc  ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool)

Returns the log-likelihood ratio for a branch according to GBM pure-birth
separately for the Brownian motion and the pure-birth
"""
function llr_gbm_b_sep(lλp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       xp  ::Array{Float64,1},
                       xc  ::Array{Float64,1},
                       α   ::Float64,
                       σλ  ::Float64,
                       βλ  ::Float64,
                       σx  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       srδt::Float64,
                       λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλp)-2

  llrbm = 0.0
  llrpb = 0.0
  llrx  = 0.0
  @avx for i in Base.OneTo(nI)
    lλpi   = lλp[i]
    lλci   = lλc[i]
    lλpi1  = lλp[i+1]
    lλci1  = lλc[i+1]
    xpi    = xp[i]
    xci    = xc[i]
    llrx  += (xp[i+1] - xpi)^2 - (xc[i+1] - xci)^2
    llrbm += (lλpi1 - lλpi - (α + βλ*xpi)*δt)^2 - 
             (lλci1 - lλci - (α + βλ*xci)*δt)^2
    llrpb += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
  end

  # standardized sum of squares
  ssrλ  = llrbm/(2.0*δt)
  ssrx  = llrx/(2.0*δt)

  # likelihood ratio
  # add to global likelihood
  llrbm *= (-0.5/((σλ*srδt)^2))
  llrbm += llrx * (-0.5/((σx*srδt)^2))
  llrpb *= (-δt)

  lλpi1 = lλp[nI+2]
  lλci1 = lλc[nI+2]

 # add final non-standard `δt`
  if fdt > 0.0
    lλpi   = lλp[nI+1]
    lλci   = lλc[nI+1]
    xci    = xv[nI+1]
    xpi    = xv[nI+1]
    xpi1   = xp[nI+2]
    xci1   = xc[nI+2]
    srfdt  = sqrt(fdt)
    λd    += (lλpi1 - lλpi - (α + βλ*xpi)*fdt)^2 - 
             (lλci1 - lλci - (α + βλ*xci)*fdt)^2
    xd     = ((xpi1 - xpi)^2 - (xci1 - xci)^2)
    ssrλ  += λd/(2.0*fdt)
    ssrx  += xd/(2.0*fdt)
    llrbm += λd*(-0.5/((σλ*srfdt)^2)) + xd*(-0.5/((σx*srfdt)^2))
    llrpb -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
  end
  #if speciation
  if λev
    llrpb += lλpi1 - lλci1
  end

  return llrbm, llrpb, ssrλ, ssrx
end




"""
    _sss_gbm(tree::T,
             α   ::Float64,
             βλ  ::Float64,
             ssλ ::Float64,
             ssx ::Float64,
             n   ::Float64)

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function _sss_gbm(tree::T,
                  α   ::Float64,
                  βλ  ::Float64,
                  ssλ ::Float64,
                  ssx ::Float64,
                  n   ::Float64) where {T <: iT}

  ssλ0, ssx0, n0 = _sss_gbm_b(lλ(tree), xv(tree), α, βλ, dt(tree), fdt(tree))

  ssλ += ssλ0
  ssx += ssx0
  n   += n0

  if def1(tree)
    ssλ, ssx, n = _sss_gbm(tree.d1, α, βλ, ssλ, ssx, n)
    ssλ, ssx, n = _sss_gbm(tree.d2, α, βλ, ssλ, ssx, n)
  end

  return ssλ, ssx, n
end




"""
    _sss_gbm_b(lλv::Array{Float64,1},
               xv ::Array{Float64,1},
               α  ::Float64,
               βλ ::Float64,
               δt ::Float64,
               fdt::Float64)

Returns the standardized sum of squares for the stochastic GBM part of a branch
for GBM birth-death.
"""
function _sss_gbm_b(lλv::Array{Float64,1},
                    xv ::Array{Float64,1},
                    α  ::Float64,
                    βλ ::Float64,
                    δt ::Float64,
                    fdt::Float64)

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ = 0.0
    ssx = 0.0
    @avx for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lλvi1 = lλv[i+1]
      xvi   = xv[i]
      ssλ  += (lλvi1 - lλvi - (α + βλ*xvi)*δt)^2
      ssx  += (xv[i+1] - xvi)^2
    end

    # standardize
    ssλ *= 1.0/(2.0*δt)
    ssx *= 1.0/(2.0*δt)

    # add final non-standard `δt`
    if fdt > 0.0
      xvi  = xv[nI+1]
      ssλ += (lλv[nI+2] - lλv[nI+1] - (α + βλ*xvi)*fdt)^2/(2.0*fdt)
      ssx += (xv[nI+2] - xvi)^2/(2.0*fdt)
      n = Float64(nI + 1)
    else
      n = Float64(nI)
    end

  return ssλ, ssx, n
end



