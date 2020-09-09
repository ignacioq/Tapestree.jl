#=

GBM birth-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(tree::iTgbmbd, 
             σ_λ ::Float64, 
             σ_μ ::Float64)

Returns the log-likelihood for a `iTgbmbd` according to GBM birth-death.
"""
function llik_gbm(tree::iTgbmbd, 
                  σ_λ ::Float64, 
                  σ_μ ::Float64)

  tsb = ts(tree)
  lλb = lλ(tree)
  lμb = lμ(tree)

  if istip(tree) 
    ll_gbm_b(tsb, lλb, lμb, σ_λ, σ_μ) + (isextinct(tree) ? lμb[end] : 0.0)
  else
    ll_gbm_b(tsb, lλb, lμb, σ_λ, σ_μ)    + (2.0*lλb[end]) + 
    llik_gbm(tree.d1::iTgbmbd, σ_λ, σ_μ) + 
    llik_gbm(tree.d2::iTgbmbd, σ_λ, σ_μ)
  end
end




"""
    ll_gbm_b(t  ::Array{Float64,1},
             lλv::Array{Float64,1},
             lμv::Array{Float64,1},
             σ_λ::Float64,
             σ_μ::Float64)

Returns the log-likelihood for a branch according to GBM birth-death.
"""
function ll_gbm_b(t  ::Array{Float64,1},
                  lλv::Array{Float64,1},
                  lμv::Array{Float64,1},
                  σ_λ::Float64,
                  σ_μ::Float64)

  ll   = 0.0
  ti   = t[1]
  lλvi = lλv[1]
  lμvi = lμv[1]
  for i in 2:lastindex(t)
    ti1   = t[i]
    δt    = ti1-ti
    lλvi1 = lλv[i]
    lμvi1 = lμv[i]
    ll   += logdnorm_tc(lλvi1, lλvi, δt*σ_λ) + 
            logdnorm_tc(lμvi1, lμvi, δt*σ_μ) - 
            δt*(exp(geoavg2(lλvi, lλvi1)) + exp(geoavg2(lμvi, lμvi1)))
    lλvi = lλvi1
    lμvi = lμvi1
    ti   = ti1 
  end

  return ll
end

