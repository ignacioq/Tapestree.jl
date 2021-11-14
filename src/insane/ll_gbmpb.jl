#=

GBM pure-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    llik_gbm(tree::iTgbmpb, 
             α   ::Float64,
             σλ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmpb` according to GBM birth-death.
"""
function llik_gbm(tree::iTgbmpb, 
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  if istip(tree) 
    ll_gbm_b(lλ(tree), α, σλ, δt, fdt(tree), srδt, false)
  else
    ll_gbm_b(lλ(tree), α, σλ, δt, fdt(tree), srδt, true) +
    llik_gbm(tree.d1::iTgbmpb, α, σλ, δt, srδt)          +
    llik_gbm(tree.d2::iTgbmpb, α, σλ, δt, srδt)
  end
end




"""
    llik_gbm(psi ::Vector{iTgbmpb},
             idf ::Vector{iBffs}, 
             α   ::Float64,
             σλ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmpb` according to GBM birth-death.
"""
function llik_gbm(psi ::Vector{iTgbmpb},
                  idf ::Vector{iBffs}, 
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(psi))
      bi  = idf[i]
      ll += llik_gbm(psi[i], α, σλ, δt, srδt)

      if !it(bi)
        ll += λt(bi)
      end
    end
  end

  return ll
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

  @inbounds @fastmath begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llbm = 0.0
    llpb = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      llbm += (lλvi1 - lλvi - α*δt)^2
      llpb += exp(0.5*(lλvi + lλvi1))
      lλvi  = lλvi1
    end

    # add to global likelihood
    ll  = llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
    ll -= llpb*δt

    lλvi1 = lλv[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      ll += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
            fdt*exp(0.5*(lλvi + lλvi1))
    end
    if λev
      ll += lλvi1
    end
  end

  return ll
end




"""
    llik_gbm(tree::iTgbmpb, 
             α   ::Float64,
             σλ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmpb` according to GBM birth-death.
"""
function llik_gbm_ssλ(tree::iTgbmpb, 
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64)

  if istip(tree) 
    ll, dλ, ssλ, nλ = ll_gbm_b_ssλ(lλ(tree), α, σλ, δt, fdt(tree), srδt, false)
  else
    ll, dλ, ssλ, nλ = ll_gbm_b_ssλ(lλ(tree), α, σλ, δt, fdt(tree), srδt, true)

    ll1, dλ1, ssλ1, nλ1 = llik_gbm_ssλ(tree.d1::iTgbmpb, α, σλ, δt, srδt)
    ll2, dλ2, ssλ2, nλ2 = llik_gbm_ssλ(tree.d2::iTgbmpb, α, σλ, δt, srδt)
  
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

  @inbounds @fastmath begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llbm = 0.0
    llpb = 0.0
    ssλ  = 0.0
    nλ   = Float64(nI)
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      llbm += (lλvi1 - lλvi - α*δt)^2
      llpb += exp(0.5*(lλvi + lλvi1))
      lλvi  = lλvi1
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
      ll  += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
             fdt*exp(0.5*(lλvi + lλvi1))
      ssλ += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
      nλ  += 1.0
    end
    if λev
      ll += lλvi1
    end
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
@inline function llr_gbm_b_sep(lλp ::Array{Float64,1},
                               lλc ::Array{Float64,1},
                               α   ::Float64,
                               σλ  ::Float64, 
                               δt  ::Float64,
                               fdt ::Float64,
                               srδt::Float64,
                               λev ::Bool)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλp)-2

    llrbm = 0.0
    llrpb = 0.0
    lλpi = lλp[1]
    lλci = lλc[1]
    @simd for i in Base.OneTo(nI)
      lλpi1  = lλp[i+1]
      lλci1  = lλc[i+1]
      llrbm += (lλpi1 - lλpi - α*δt)^2 - (lλci1 - lλci - α*δt)^2
      llrpb += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
      lλpi   = lλpi1
      lλci   = lλci1
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
      ssrλ  += ((lλpi1 - lλpi - α*fdt)^2 - (lλci1 - lλci - α*fdt)^2)/(2.0*fdt) 
      llrbm += lrdnorm_bm_x(lλpi1, lλpi + α*fdt, 
                            lλci1, lλci + α*fdt, sqrt(fdt)*σλ)
      llrpb -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
    end
    #if speciation
    if λev
      llrpb += lλpi1 - lλci1
    end
  end

  return llrbm, llrpb, ssrλ
end




"""
    sss_gbm(psi::Vector{T}, α::Float64) where {T <: iTgbm}

Returns the standardized sum of squares a `iTgbm` according 
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(psi::Vector{T}, α::Float64) where {T <: iTgbm}

  n   = 0.0
  ssλ = 0.0
  for ψi in psi
    ssλ, n = _sss_gbm(ψi, α, ssλ, n)
  end

  return ssλ, n
end




"""
    sss_gbm(tree::T, α::Float64, ssλ::Float64, n::Float64)

Returns the standardized sum of squares a `iTgbm` according 
to GBM birth-death for a `σ` proposal.
"""
function _sss_gbm(tree::T, α::Float64, ssλ::Float64, n::Float64) where {T <: iTgbm}

  ssλ0, n0 = _sss_gbm_b(lλ(tree), α, dt(tree), fdt(tree))

  ssλ += ssλ0
  n   += n0

  if isdefined(tree, :d1) 
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
@inline function _sss_gbm_b(lλv::Array{Float64,1},
                            α  ::Float64,
                            δt ::Float64, 
                            fdt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ  = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      ssλ  += (lλvi1 - lλvi - α*δt)^2
      lλvi  = lλvi1
    end

    # add to global likelihood
    ssλ *= 1.0/(2.0*δt)

    # add final non-standard `δt`
    if fdt > 0.0
      ssλ += (lλv[nI+2] - lλvi - α*fdt)^2/(2.0*fdt)
      n = Float64(nI + 1)
    else
      n = Float64(nI)
    end
  end

  return ssλ, n
end




"""
    deltaλ(tree::T, dλ::Float64, l::Float64) where {T <: iTgbm}

Returns the log-likelihood ratio for according to GBM 
for a drift `α` proposal.
"""
function deltaλ(psi::Vector{T}) where {T <: iTgbm}

  dλ = 0.0

  for ψi in psi
    dλ += _deltaλ(ψi)
  end

  return dλ
end




"""
    _deltaλ(tree::T) where {T <: iTgbm}

Returns the log-likelihood ratio for a `iTgbmpb` according 
to GBM birth-death for a `α` proposal.
"""
function _deltaλ(tree::T) where {T <: iTgbm}

  lλv = lλ(tree)

  if isdefined(tree, :d1) 
    lλv[end] - lλv[1] + _deltaλ(tree.d1) + _deltaλ(tree.d2)
  else
    lλv[end] - lλv[1]
  end
end


