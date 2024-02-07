#=

`gbmct` likelihood for constant turnover

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    llik_gbm(Ξ   ::Vector{iTct},
             idf ::Vector{iBffs},
             α   ::Float64,
             σλ  ::Float64,
             ϵ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTct` according to `gbmct`.
"""
function llik_gbm(Ξ   ::Vector{iTct},
                  idf ::Vector{iBffs},
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      ll += llik_gbm(Ξ[i], α, σλ, ϵ, δt, srδt)
      if d2(idf[i]) > 0
        ll += λt(Ξ[i])
      end
    end
  end

  return ll
end




"""
    llik_gbm(tree::iTct,
             α   ::Float64,
             σλ  ::Float64,
             ϵ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTct` according to `gbmct`.
"""
function llik_gbm(tree::iTct,
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  if istip(tree)
    ll_gbm_b_ϵ(lλ(tree), α, σλ, ϵ, δt, fdt(tree), srδt, false, isextinct(tree))
  else
    ll_gbm_b_ϵ(lλ(tree), α, σλ, ϵ, δt, fdt(tree), srδt, true, false) +
    llik_gbm(tree.d1, α, σλ, ϵ, δt, srδt)                            +
    llik_gbm(tree.d2, α, σλ, ϵ, δt, srδt)
  end
end




"""
    ll_gbm_b_ϵ(lλv ::Array{Float64,1},
               α   ::Float64,
               σλ  ::Float64,
               ϵ   ::Float64,
               δt  ::Float64,
               fdt ::Float64,
               srδt::Float64,
               λev ::Bool,
               μev ::Bool)

Returns the log-likelihood for a branch according to `gbmct`.
"""
function ll_gbm_b_ϵ(lλv ::Array{Float64,1},
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    δt  ::Float64,
                    fdt ::Float64,
                    srδt::Float64,
                    λev ::Bool,
                    μev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2

  llbm = 0.0
  llct = 0.0
  @turbo for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    llbm += (lλvi1 - lλvi - α*δt)^2
    llct += exp(0.5*(lλvi + lλvi1))
  end

  # add to global likelihood
  ll  = llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
 
  # birth-death likelihood
  ll -= llct*δt*(1.0 + ϵ)

  lλvi1 = lλv[nI+2]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi = lλv[nI+1]
    ll  += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
           fdt*exp(0.5*(lλvi + lλvi1))*(1.0 + ϵ)
  end

  # if speciation
  if λev
    ll += lλvi1
  # if extinction
  elseif μev
    ll += lλvi1 + log(ϵ)
  end

  return ll
end




"""
    llik_gbm_ssλ(tree::iTct,
                 α   ::Float64,
                 σλ  ::Float64,
                 ϵ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

Returns the log-likelihood for a `iTct` according to `gbmct`.
"""
function llik_gbm_ssλ(tree::iTct,
                      α   ::Float64,
                      σλ  ::Float64,
                      ϵ   ::Float64,
                      δt  ::Float64,
                      srδt::Float64)

  if istip(tree)
    ll, dλ, ssλ, Σλ, nλ =
      ll_gbm_b_ϵ_ssλ(lλ(tree), α, σλ, ϵ, δt, fdt(tree), srδt,
        false, isextinct(tree))
  else
    ll, dλ, ssλ, Σλ, nλ =
      ll_gbm_b_ϵ_ssλ(lλ(tree), α, σλ, ϵ, δt, fdt(tree), srδt,
        true, false)

    ll1, dλ1, ssλ1, Σλ1, nλ1 =
      llik_gbm_ssλ(tree.d1, α, σλ, ϵ, δt, srδt)
    ll2, dλ2, ssλ2, Σλ2, nλ2 =
      llik_gbm_ssλ(tree.d2, α, σλ, ϵ, δt, srδt)

    ll  += ll1  + ll2
    dλ  += dλ1  + dλ2
    ssλ += ssλ1 + ssλ2
    Σλ  += Σλ1  + Σλ2
    nλ  += nλ1  + nλ2
  end

  return ll, dλ, ssλ, Σλ, nλ
end




"""
    ll_gbm_b_ϵ_ssλ(lλv ::Array{Float64,1},
                   α   ::Float64,
                   σλ  ::Float64,
                   ϵ   ::Float64,
                   δt  ::Float64,
                   fdt ::Float64,
                   srδt::Float64,
                   λev ::Bool,
                   μev ::Bool)

Returns the log-likelihood for a branch according to `gbmct`.
"""
function ll_gbm_b_ϵ_ssλ(lλv ::Array{Float64,1},
                        α   ::Float64,
                        σλ  ::Float64,
                        ϵ   ::Float64,
                        δt  ::Float64,
                        fdt ::Float64,
                        srδt::Float64,
                        λev ::Bool,
                        μev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2

  llbm = 0.0
  llct = 0.0
  @turbo for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    llbm += (lλvi1 - lλvi - α*δt)^2
    llct += exp(0.5*(lλvi + lλvi1))
  end

  # standardized sum of squares
  ssλ = llbm/(2.0*δt)
  nλ  = Float64(nI)

  # add to global likelihood
  ll    = llbm *
          (-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
  llct *= δt
  Σλ    = llct
  ll   -= llct*(1.0 + ϵ)

  lλvi1 = lλv[nI+2]

  dλ = lλvi1 - lλv[1]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi = lλv[nI+1]
    lli  = fdt * exp(0.5*(lλvi + lλvi1))
    ll  += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) - lli*(1.0 + ϵ)
    ssλ += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
    nλ  += 1.0
    Σλ  += lli
  end
  # if speciation
  if λev
    ll += lλvi1
  end
  # if extinction
  if μev
    ll += lλvi1 + log(ϵ)
  end

  return ll, dλ, ssλ, Σλ, nλ
end




"""
    Σλ_gbm(tree::iTct)

Returns the sum of `λ` rates for a `iTct` according
to `gbmct` for a `ϵ` proposal.
"""
function Σλ_gbm(tree::iTct)

  if istip(tree)
    Σλ_gbm_b(lλ(tree), dt(tree), fdt(tree))
  else
    Σλ_gbm_b(lλ(tree), dt(tree), fdt(tree)) +
    Σλ_gbm(tree.d1::iTct) +
    Σλ_gbm(tree.d2::iTct)
  end
end




"""
    Σλ_gbm_b(lλv::Array{Float64,1},
             δt ::Float64,
             fdt::Float64)

Returns the sum of `λ` rates for a `iTct` branch according
to `gbmct` for a `ϵ` proposal.
"""
function Σλ_gbm_b(lλv::Array{Float64,1},
                  δt ::Float64,
                  fdt::Float64)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2

  Σλ   = 0.0
  @turbo for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    Σλ   += exp(0.5*(lλvi + lλvi1))
  end

  # add to global sum
  Σλ *= δt

  # add final non-standard `δt`
  if fdt > 0.0
    Σλ += fdt * exp(0.5*(lλv[nI+2] + lλv[nI+1]))
  end

  return Σλ
end



"""
    llr_gbm_b_sep(lλp ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  δt  ::Float64,
                  fdt::Float64,
                  srδt::Float64,
                  λev ::Bool,
                  μev ::Bool)
Returns the log-likelihood for a branch according to `gbmct`
separately (for gbm and bd).
"""
function llr_gbm_b_sep(lλp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       α   ::Float64,
                       σλ  ::Float64,
                       ϵ   ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       srδt::Float64,
                       λev ::Bool,
                       μev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλc)-2

  llrbm = 0.0
  llrct = 0.0
  @turbo for i in Base.OneTo(nI)
    lλpi   = lλp[i]
    lλci   = lλc[i]
    lλpi1  = lλp[i+1]
    lλci1  = lλc[i+1]
    llrbm += (lλpi1 - lλpi - α*δt)^2 - (lλci1 - lλci - α*δt)^2
    llrct += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
  end

  # standardized sum of squares
  ssrλ = llrbm/(2.0*δt)
  Σrλ  = llrct * δt

  # overall
  llrbm *= (-0.5/((σλ*srδt)^2))
  llrct *= -δt*(1.0 + ϵ)

  lλpi1 = lλp[nI+2]
  lλci1 = lλc[nI+2]

  # add final non-standard `δt`
  if fdt > 0.0
    lλpi  = lλp[nI+1]
    lλci  = lλc[nI+1]
    ssrλ  += ((lλpi1 - lλpi - α*fdt)^2 - (lλci1 - lλci - α*fdt)^2)/(2.0*fdt)
    llri   = fdt * (exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
    Σrλ   += llri
    llrbm += lrdnorm_bm_x(lλpi1, lλpi + α*fdt,
                          lλci1, lλci + α*fdt, sqrt(fdt)*σλ)
    llrct -= (1.0 + ϵ) * llri
  end
  # if speciation or extinction
  if λev || μev
    llrct += lλpi1 - lλci1
  end

  return llrbm, llrct, ssrλ, Σrλ
end




"""
    llik_gbm_lλshift(Ξ      ::Vector{iTct},
                     δt     ::Float64,
                     lλshift::Float64,
                     ϵ      ::Float64)

Returns the exponential term of the birth-death log-likelihood ratio 
for a lλshift on a `iTct`.
"""
function llr_gbm_lλshift(Ξ      ::Vector{iTct},
                         δt     ::Float64,
                         lλshift::Float64,
                         ϵ      ::Float64)
  @inbounds begin
    explλshiftm1 = exp(lλshift)-1
    llr = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      llr += _llr_gbm_lλshift(Ξ[i], δt, lλshift, explλshiftm1, ϵ)
    end
  end

  return llr
end




"""
    _llr_gbm_lλshift(tree        ::iTct,
                     δt          ::Float64,
                     lλshift     ::Float64,
                     explλshiftm1::Float64,
                     ϵ           ::Float64)

Returns the exponential term of the birth-death log-likelihood ratio 
for a lλshift on a `iTct`.
"""
function _llr_gbm_lλshift(tree        ::iTct,
                          δt          ::Float64,
                          lλshift     ::Float64,
                          explλshiftm1::Float64,
                          ϵ           ::Float64)
  @inbounds begin
    llr = 0.0
    lλtree = lλ(tree)
    nI = lastindex(lλtree)-2
    fdti = fdt(tree)

    for i in Base.OneTo(nI)
      llr -= exp(0.5*(lλtree[i] + lλtree[i+1]))
    end
    llr *= explλshiftm1*δt*(1.0 + ϵ)

    # add final non-standard `δt`
    if fdti > 0.0
      llr -= (exp(0.5*(lλtree[nI+1] + lλtree[nI+2])))*explλshiftm1*fdti*(1.0 + ϵ)
    end

    if !istip(tree)
      llr += _llr_gbm_lλshift(tree.d1, δt, lλshift, explλshiftm1, ϵ)
      llr += _llr_gbm_lλshift(tree.d2, δt, lλshift, explλshiftm1, ϵ)
    end
  end

  return llr
end



