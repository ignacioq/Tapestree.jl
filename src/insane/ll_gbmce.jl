#=

`gbmce` likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(Ξ   ::Vector{iTce},
             idf ::Vector{iBffs},
             α   ::Float64,
             σλ  ::Float64,
             μ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTce` according to `gbm-ce`.
"""
function llik_gbm(Ξ   ::Vector{iTce},
                  idf ::Vector{iBffs},
                  α   ::Float64,
                  σλ  ::Float64,
                  μ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      ll += llik_gbm(Ξ[i], α, σλ, μ, δt, srδt)
      if !it(bi)
        ll += λt(bi)
      end
    end
  end

  return ll
end




"""
    llik_gbm(tree::iTce,
             α   ::Float64,
             σλ  ::Float64,
             μ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTce` according to `gbmce`.
"""
function llik_gbm(tree::iTce,
                  α   ::Float64,
                  σλ  ::Float64,
                  μ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  if istip(tree)
    ll_gbm_b(lλ(tree), α, σλ, μ, δt, fdt(tree), srδt, false, isextinct(tree))
  else
    ll_gbm_b(lλ(tree), α, σλ, μ, δt, fdt(tree), srδt, true, false) +
    llik_gbm(tree.d1, α, σλ, μ, δt, srδt)                          +
    llik_gbm(tree.d2, α, σλ, μ, δt, srδt)
  end
end




"""
    ll_gbm_b(lλv ::Array{Float64,1},
             α   ::Float64,
             σλ  ::Float64,
             μ   ::Float64,
             δt  ::Float64,
             fdt ::Float64,
             srδt::Float64,
             λev ::Bool,
             μev ::Bool)

Returns the log-likelihood for a branch according to `gbmce`.
"""
@inline function ll_gbm_b(lλv ::Array{Float64,1},
                          α   ::Float64,
                          σλ  ::Float64,
                          μ   ::Float64,
                          δt  ::Float64,
                          fdt ::Float64,
                          srδt::Float64,
                          λev ::Bool,
                          μev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2

  llλ  = 0.0
  llbd = 0.0
  @avx for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    llλ  += (lλvi1 - lλvi - α*δt)^2
    llbd += exp(0.5*(lλvi + lλvi1))
  end

  # add to global likelihood
  ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))

  # add to global likelihood
  llbd += Float64(nI)*μ
  ll   -= llbd*δt

  lλvi1 = lλv[nI+2]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi  = lλv[nI+1]
    ll += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
          fdt*(exp(0.5*(lλvi + lλvi1)) + μ)
  end
  # if speciation
  if λev
    ll += lλvi1
  end
  # if extinction
  if μev
    ll += log(μ)
  end

  return ll
end




"""
    llik_gbm_ssλ(tree::iTce,
                 α   ::Float64,
                 σλ  ::Float64,
                 μ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

Returns the log-likelihood for a `iTce` according to GBM birth-death.
"""
function llik_gbm_ssλ(tree::iTce,
                      α   ::Float64,
                      σλ  ::Float64,
                      μ   ::Float64,
                      δt  ::Float64,
                      srδt::Float64)

  if istip(tree)
    ll, dλ, ssλ, nλ =
      ll_gbm_b_ssλ(lλ(tree), α, σλ, μ, δt, fdt(tree), srδt,
        false, isextinct(tree))
  else
    ll, dλ, ssλ, nλ =
      ll_gbm_b_ssλ(lλ(tree), α, σλ, μ, δt, fdt(tree), srδt,
        true, false)

    ll1, dλ1, ssλ1, nλ1 =
      llik_gbm_ssλ(tree.d1::iTce, α, σλ, μ, δt, srδt)
    ll2, dλ2, ssλ2, nλ2 =
      llik_gbm_ssλ(tree.d2::iTce, α, σλ, μ, δt, srδt)

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
                 μ   ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64,
                 λev ::Bool,
                 μev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_gbm_b_ssλ(lλv ::Array{Float64,1},
                      α   ::Float64,
                      σλ  ::Float64,
                      μ   ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64,
                      λev ::Bool,
                      μev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2

  llbm = 0.0
  llbd = 0.0
  @avx for i in Base.OneTo(nI)
    lλvi  = lλv[i]
    lλvi1 = lλv[i+1]
    llbm += (lλvi1 - lλvi - α*δt)^2
    llbd += exp(0.5*(lλvi + lλvi1))
  end

  # standardized sum of squares
  ssλ  = llbm/(2.0*δt)
  nλ   = Float64(nI)

  # add to global likelihood
  ll    = llbm *
          (-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
  llbd += Float64(nI)*μ
  ll   -= llbd*δt

  lλvi1 = lλv[nI+2]

  dλ = lλvi1 - lλv[1]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi  = lλv[nI+1]
    ll  += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
           fdt*(exp(0.5*(lλvi + lλvi1)) + μ)
    ssλ += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
    nλ  += 1.0
  end
  if λev
    ll += lλvi1
  end
  if μev
    ll += log(μ)
  end

  return ll, dλ, ssλ, nλ
end
