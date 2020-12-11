#=

GBM birth-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(tree::iTgbmbd, 
             σλ  ::Float64, 
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmbd` according to GBM birth-death.
"""
function llik_gbm(tree::iTgbmbd, 
                  σλ  ::Float64, 
                  σμ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  tsb = ts(tree)
  lλb = lλ(tree)
  lμb = lμ(tree)

  if istip(tree) 
    ll_gbm_b(tsb, lλb, lμb, σλ, σμ, δt, srδt) + 
    (isextinct(tree) ? lμb[end] : 0.0)
  else
    ll_gbm_b(tsb, lλb, lμb, σλ, σμ, δt, srδt)    + 
    log(2.0) + lλb[end]                          + 
    llik_gbm(tree.d1::iTgbmbd, σλ, σμ, δt, srδt) + 
    llik_gbm(tree.d2::iTgbmbd, σλ, σμ, δt, srδt)
  end
end




"""
    ll_gbm_b(t   ::Array{Float64,1},
             lλv ::Array{Float64,1},
             lμv ::Array{Float64,1},
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64, 
             srδt::Float64)

Returns the log-likelihood for a branch according to GBM birth-death.
"""
function ll_gbm_b(t   ::Array{Float64,1},
                  lλv ::Array{Float64,1},
                  lμv ::Array{Float64,1},
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64, 
                  srδt::Float64)

  @inbounds @fastmath begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    llλ  = 0.0
    llμ  = 0.0
    llbd = 0.0
    lλvi = lλv[1]
    lμvi = lμv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi)^2
      llμ  += (lμvi1 - lμvi)^2
      llbd += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)) 
      lλvi  = lλvi1
      lμvi  = lμvi1
    end

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) + 
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    # add to global likelihood
    ll -= llbd*δt

    # add final non-standard `δt`
    δtf   = t[nI+2] - t[nI+1]
    srδtf = sqrt(δtf)
    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]
    ll += ldnorm_bm(lλvi1, lλvi, srδtf*σλ)
    ll += ldnorm_bm(lμvi1, lμvi, srδtf*σμ)
    ll -= δtf*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)))
  end

  return ll
end






"""
estimate survival probability
"""





