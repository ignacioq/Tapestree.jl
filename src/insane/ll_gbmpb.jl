#=

GBM pure-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(tree::iTgbmpb, 
             σ²λ ::Float64)

Returns the log-likelihood for a `iTgbmpb` according to GBM birth-death.
"""
function llik_gbm(tree::iTgbmpb, 
                  σ²λ ::Float64,
                  δt  ::Float64)

  tsb = ts(tree)
  lλb = lλ(tree)

  if istip(tree) 
    ll_gbm_b(tsb, lλb, σ²λ, δt)
  else
    ll_gbm_b(tsb, lλb, σ²λ, δt)     +
    2.0*lλb[end]                    +
    llik_gbm(tree.d1::iTgbmpb, σ²λ, δt) +
    llik_gbm(tree.d2::iTgbmpb, σ²λ, δt)
  end
end




"""
    ll_gbm_b(t  ::Array{Float64,1},
             lλv::Array{Float64,1},
             σ²λ::Float64, 
             δt  ::Float64)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_gbm_b(t  ::Array{Float64,1},
                  lλv::Array{Float64,1},
                  σ²λ::Float64, 
                  δt ::Float64)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    llbm = 0.0
    llpb = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      llbm += (lλvi1 - lλvi)^2
      llpb += exp(0.5*(lλvi + lλvi1))
      lλvi  = lλvi1
    end

    # add to global likelihood
    ll  = llbm*(-1.0/(2.0*σ²λ*δt)) - 0.5*Float64(nI)*log(σ²λ*δt)
    ll -= llpb*δt

    # add final non-standard `δt`
    δtf   = t[nI+2] - t[nI+1]
    lλvi1 = lλv[nI+2]
    ll += logdnorm_tc(lλvi1, lλvi, δtf*σ²λ)
    ll -= δtf*exp(0.5*(lλvi + lλvi1))
  end

  return ll
end
