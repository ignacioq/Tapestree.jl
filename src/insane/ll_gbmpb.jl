#=

GBM pure-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(Ξ   ::Vector{iTpb},
             idf ::Vector{iBffs},
             α   ::Float64,
             σλ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTpb` according to GBM birth-death.
"""
function llik_gbm(Ξ   ::Vector{iTpb},
                  idf ::Vector{iBffs},
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      ll += llik_gbm(Ξ[i], α, σλ, δt, srδt)
      if d2(idf[i]) > 0
        ll += λt(Ξ[i])
      end
    end
  end

  return ll
end




"""
    llik_gbm(tree::iTpb,
             α   ::Float64,
             σλ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTpb` according to GBM birth-death.
"""
function llik_gbm(tree::iTpb,
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  if istip(tree)
    ll_gbm_b(lλ(tree), α, σλ, δt, fdt(tree), srδt, false)
  else
    ll_gbm_b(lλ(tree), α, σλ, δt, fdt(tree), srδt, true) +
    llik_gbm(tree.d1::iTpb, α, σλ, δt, srδt)          +
    llik_gbm(tree.d2::iTpb, α, σλ, δt, srδt)
  end
end




"""
    ll_gbm_b(lλv ::Array{Float64,1},
             α   ::Float64,
             σλ  ::Float64,
             δt  ::Float64,
             fdt ::Float64,
             srδt::Float64,
             λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_gbm_b(lλv ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2

  llbm = 0.0
  llpb = 0.0
  @turbo for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    llbm += (lλvi1 - lλvi - α*δt)^2
    llpb += exp(0.5*(lλvi + lλvi1))
  end

  # add to global likelihood
  ll  = llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
  ll -= llpb*δt

  lλvi1 = lλv[nI+2]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi = lλv[nI+1]
    ll   += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
            fdt*exp(0.5*(lλvi + lλvi1))
  end
  if λev
    ll += lλvi1
  end

  return ll
end




"""
    llik_gbm_ssλ(tree::iTpb,
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

Returns the log-likelihood for a `iTpb` according to GBM birth-death.
"""
function llik_gbm_ssλ(tree::iTpb,
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64)

  if istip(tree)
    ll, dλ, ssλ, nλ = ll_gbm_b_ssλ(lλ(tree), α, σλ, δt, fdt(tree), srδt, false)
  else
    ll, dλ, ssλ, nλ = ll_gbm_b_ssλ(lλ(tree), α, σλ, δt, fdt(tree), srδt, true)

    ll1, dλ1, ssλ1, nλ1 = llik_gbm_ssλ(tree.d1::iTpb, α, σλ, δt, srδt)
    ll2, dλ2, ssλ2, nλ2 = llik_gbm_ssλ(tree.d2::iTpb, α, σλ, δt, srδt)

    ll  += ll1  + ll2
    dλ  += dλ1  + dλ2
    ssλ += ssλ1 + ssλ2
    nλ  += nλ1  + nλ2
  end

  return ll, dλ, ssλ, nλ
end




"""
    ll_gbm_b_ssλ(lλv ::Array{Float64,1},
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64,
                 λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_gbm_b_ssλ(lλv ::Array{Float64,1},
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64,
                      λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2
  nλ = Float64(nI)

  llbm = 0.0
  llpb = 0.0
  @turbo for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    llbm += (lλvi1 - lλvi - α*δt)^2
    llpb += exp(0.5*(lλvi + lλvi1))
  end

  # standardized sum of squares
  ssλ  = llbm/(2.0*δt)

  # add to global likelihood
  ll  = llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
  ll -= llpb*δt

  lλvi1 = lλv[nI+2]

  dλ = lλvi1 - lλv[1]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi = lλv[nI+1]
    ll  += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
           fdt*exp(0.5*(lλvi + lλvi1))
    ssλ += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
    nλ  += 1.0
  end
  if λev
    ll += lλvi1
  end

  return ll, dλ, ssλ, nλ
end





"""
    llr_gbm_b_sep(lλp ::Array{Float64,1},
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
function llr_gbm_b_sep(lλp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       srδt::Float64,
                       λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλp)-2

  llrbm = 0.0
  llrpb = 0.0
  @turbo for i in Base.OneTo(nI)
    lλpi   = lλp[i]
    lλci   = lλc[i]
    lλpi1  = lλp[i+1]
    lλci1  = lλc[i+1]
    llrbm += (lλpi1 - lλpi - α*δt)^2 - (lλci1 - lλci - α*δt)^2
    llrpb += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
  end

  # standardized sum of squares
  ssrλ  = llrbm/(2.0*δt)
  # add to global likelihood
  llrbm *= (-0.5/((σλ*srδt)^2))
  llrpb *= (-δt)

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
  #if speciation
  if λev
    llrpb += lλpi1 - lλci1
  end

  return llrbm, llrpb, ssrλ
end




"""
    _sss_gbm(tree::T,
             α   ::Float64,
             ssλ ::Float64,
             n   ::Float64) where {T <: iT}

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function _sss_gbm(tree::T,
                  α   ::Float64,
                  ssλ ::Float64,
                  n   ::Float64) where {T <: iT}

  ssλ0, n0 = _sss_gbm_b(lλ(tree), α, dt(tree), fdt(tree))

  ssλ += ssλ0
  n   += n0

  if def1(tree)
    ssλ, n = _sss_gbm(tree.d1, α, ssλ, n)
    ssλ, n = _sss_gbm(tree.d2, α, ssλ, n)
  end

  return ssλ, n
end




"""
    sss_gbm_b(lλv::Array{Float64,1},
              α  ::Float64,
              δt ::Float64,
              fdt::Float64)

Returns the standardized sum of squares for the stochastic GBM part of a branch
for GBM birth-death.
"""
function _sss_gbm_b(lλv::Array{Float64,1},
                    α  ::Float64,
                    δt ::Float64,
                    fdt::Float64)


    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ  = 0.0
    @turbo for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lλvi1 = lλv[i+1]
      ssλ  += (lλvi1 - lλvi - α*δt)^2
    end

    # standardize
    ssλ *= 1.0/(2.0*δt)

    # add final non-standard `δt`
    if fdt > 0.0
      ssλ += (lλv[nI+2] - lλv[nI+1] - α*fdt)^2/(2.0*fdt)
      n = Float64(nI + 1)
    else
      n = Float64(nI)
    end

  return ssλ, n
end




"""
    deltaλ(Ξ::Vector{T}) where {T <: iT}

Returns the log-likelihood ratio for according to GBM
for a drift `α` proposal.
"""
function deltaλ(Ξ::Vector{T}) where {T <: iT}

  dλ = 0.0

  for ξi in Ξ
    dλ += _deltaλ(ξi)
  end

  return dλ
end




"""
    _deltaλ(tree::T) where {T <: iT}

Returns the log-likelihood ratio for a `iTpb` according
to GBM birth-death for a `α` proposal.
"""
function _deltaλ(tree::T) where {T <: iT}

  lλv = lλ(tree)

  if def1(tree)
    lλv[end] - lλv[1] + _deltaλ(tree.d1) + _deltaλ(tree.d2)
  else
    lλv[end] - lλv[1]
  end
end




"""
    llik_gbm_lλshift(Ξ      ::Vector{T},
                     δt     ::Float64,
                     lλshift::Float64) where {T <: iT}

Returns the exponential term of the birth-death log-likelihood ratio 
for a lλshift on an `iT`.
"""
function llr_gbm_lλshift(Ξ      ::Vector{T},
                         δt     ::Float64,
                         lλshift::Float64) where {T <: iT}
  @inbounds begin
    explλshiftm1 = exp(lλshift)-1
    llr = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      llr += _llr_gbm_lλshift(Ξ[i], δt, lλshift, explλshiftm1)
    end
  end

  return llr
end




"""
    _llr_gbm_lλshift(tree         ::T,
                     δt           ::Float64,
                     lλshift      ::Float64,
                     explλshiftm1 ::Float64) where {T <: iT}

Returns the exponential term of the birth-death log-likelihood ratio 
for a lλshift on an `iT`.
"""
function _llr_gbm_lλshift(tree         ::T,
                          δt           ::Float64,
                          lλshift      ::Float64,
                          explλshiftm1 ::Float64) where {T <: iT}
  @inbounds begin
    llr = 0.0
    lλtree = lλ(tree)
    nI = lastindex(lλtree)-2
    fdti = fdt(tree)

    for i in Base.OneTo(nI)
      llr -= exp(0.5*(lλtree[i] + lλtree[i+1]))
    end
    llr *= explλshiftm1*δt

    # add final non-standard `δt`
    if fdti > 0.0
      llr -= (exp(0.5*(lλtree[nI+1] + lλtree[nI+2])))*explλshiftm1*fdti
    end

    if !istip(tree)
      llr += _llr_gbm_lλshift(tree.d1, δt, lλshift, explλshiftm1)
      llr += _llr_gbm_lλshift(tree.d2, δt, lλshift, explλshiftm1)
    end
  end

  return llr
end




