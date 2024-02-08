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
      bi  = idf[i]
      ll += llik_gbm(Ξ[i], α, σλ, δt, srδt)
      if d2(bi) > 0
        ll += λt(bi)
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





