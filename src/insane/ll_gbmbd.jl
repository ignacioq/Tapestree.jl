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
    llr_gbm_bm(tree::iTgbmbd, 
               σp  ::Float64,
               σc  ::Float64,
               srδt::Float64,
               lf  ::Function)

Returns the log-likelihood ration for a `iTgbmpb` according 
to GBM birth-death for a `σ` proposal.
"""
function llr_gbm_bm(tree::iTgbmbd, 
                    σp  ::Float64,
                    σc  ::Float64,
                    srδt::Float64,
                    lf  ::Function)

  if istip(tree) 
    llr_gbm_bm(ts(tree), lf(tree), σp, σc, srδt)
  else
    llr_gbm_bm(ts(tree), lf(tree), σp, σc, srδt)   +
    llr_gbm_bm(tree.d1::iTgbmbd, σp, σc, srδt, lf) +
    llr_gbm_bm(tree.d2::iTgbmbd, σp, σc, srδt, lf)
  end
end




"""
    llik_gbm_f(tree::iTgbmbd,
               σλ  ::Float64,
               σμ  ::Float64,
               δt  ::Float64,
               srδt::Float64)

Estimate gbm birth-death likelihood for the tree in a branch.
"""
function llik_gbm_f(tree::iTgbmbd,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  tsb = ts(tree)
  lλb = lλ(tree)
  lμb = lμ(tree)

  if istip(tree)
    ll = ll_gbm_b(tsb, lλb, lμb, σλ, σμ, δt, srδt) + 
         (isextinct(tree) ? lμb[end] : 0.0)
  else
    ll = ll_gbm_b(tsb, lλb, lμb, σλ, σμ, δt, srδt) + 
         log(2.0) + lλb[end]

    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      return ll
    elseif ifx1
      ll += llik_gbm_f(tree.d1::iTgbmbd, σλ, σμ, δt, srδt) +
            llik_gbm(  tree.d2::iTgbmbd, σλ, σμ, δt, srδt)
    else
      ll += llik_gbm(  tree.d1::iTgbmbd, σλ, σμ, δt, srδt) + 
            llik_gbm_f(tree.d2::iTgbmbd, σλ, σμ, δt, srδt)
    end
  end

  return ll
end





"""
    br_ll_cbd(tree::sTbd,
              λc  ::Float64, 
              μc  ::Float64,
              dri ::BitArray{1}, 
              ldr ::Int64,
              ix  ::Int64)

Returns gbm birth-death likelihood for whole branch `br`.
"""
function br_ll_gbm(tree::iTgbmbd,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   dri ::BitArray{1},
                   ldr ::Int64,
                   ix  ::Int64)

  if ix === ldr
    return llik_gbm_f(tree::iTgbmbd, σλ, σμ, δt, srδt)
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTgbmbd)
    if ifx1 && isfix(tree.d2::iTgbmbd)
      ix += 1
      if dri[ix]
        ll = br_ll_gbm(tree.d1::iTgbmbd, σλ, σμ, δt, srδt, dri, ldr, ix)
      else
        ll = br_ll_gbm(tree.d2::iTgbmbd, σλ, σμ, δt, srδt, dri, ldr, ix)
      end
    elseif ifx1
      ll = br_ll_gbm(tree.d1::iTgbmbd, σλ, σμ, δt, srδt, dri, ldr, ix)
    else
      ll = br_ll_gbm(tree.d2::iTgbmbd, σλ, σμ, δt, srδt, dri, ldr, ix)
    end
  end

  return ll
end




