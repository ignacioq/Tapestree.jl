#=

GBM pure-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(tree::iTgbmpb, 
             σλ  ::Float64,
             δt  ::Float64
             srδt::Float64)

Returns the log-likelihood for a `iTgbmpb` according to GBM birth-death.
"""
function llik_gbm(tree::iTgbmpb, 
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  tsb = ts(tree)
  lλb = lλ(tree)

  if istip(tree) 
    ll_gbm_b(tsb, lλb, σλ, δt, srδt)
  else
    ll_gbm_b(tsb, lλb, σλ, δt, srδt)         +
    lλb[end]                             +
    llik_gbm(tree.d1::iTgbmpb, σλ, δt, srδt) +
    llik_gbm(tree.d2::iTgbmpb, σλ, δt, srδt)
  end
end




"""
    ll_gbm_b(t   ::Array{Float64,1},
             lλv ::Array{Float64,1},
             σλ  ::Float64, 
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_gbm_b(t   ::Array{Float64,1},
                  lλv ::Array{Float64,1},
                  σλ  ::Float64, 
                  δt  ::Float64,
                  srδt::Float64)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    lλvi = lλv[1]
    llbm = 0.0
    llpb = 0.5*exp(lλvi)
    @simd for i in Base.OneTo(nI-1)
      lλvi1 = lλv[i+1]
      llbm += (lλvi1 - lλvi)^2
      llpb += exp(lλvi1)
      lλvi  = lλvi1
    end

    # last interval for the end 0.5 in pure-birth
    lλvi1 = lλv[nI+1]
    llbm += (lλvi1 - lλvi)^2
    elλv1 = exp(lλvi1)
    llpb += 0.5*elλv1
    lλvi  = lλvi1

    # add to global likelihood
    ll  = llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
    ll -= llpb*δt

    # add final non-standard `δt`
    δtf   = t[nI+2] - t[nI+1]
    lλvi1 = lλv[nI+2]
    ll += ldnorm_bm(lλvi1, lλvi, sqrt(δtf)*σλ)
    ll -= δtf*0.5*(elλv1 + exp(lλvi1))
  end

  return ll
end





"""
    llr_gbm_b(t   ::Array{Float64,1},
              lλp ::Array{Float64,1},
              lλc ::Array{Float64,1},
              σλ  ::Float64, 
              δt  ::Float64,
              srδt::Float64)

Returns the log-likelihood ratio for a branch according to GBM pure-birth.
"""
function llr_gbm_b(t   ::Array{Float64,1},
                   lλp ::Array{Float64,1},
                   lλc ::Array{Float64,1},
                   σλ  ::Float64, 
                   δt  ::Float64,
                   srδt::Float64)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    lλpi = lλp[1]
    lλci = lλc[1]
    llbm = 0.0
    llpb = 0.5*exp(lλpi) - 0.5*exp(lλci)
    @simd for i in Base.OneTo(nI-1)
      lλpi1 = lλp[i+1]
      lλci1 = lλc[i+1]
      llbm += (lλpi1 - lλpi)^2 - (lλci1 - lλci)^2
      llpb += exp(lλpi1) - exp(lλci1)
      lλpi  = lλpi1
      lλci  = lλci1
    end

    # last interval for the end 0.5 in pure-birth
    lλpi1 = lλp[nI+1]
    lλci1 = lλc[nI+1]
    llbm += (lλpi1 - lλpi)^2 - (lλci1 - lλci)^2
    elλp1 = exp(lλpi1)
    elλc1 = exp(lλci1)
    llpb += 0.5*elλp1 - 0.5*elλc1
    lλpi  = lλpi1
    lλci  = lλci1

    # add to global likelihood
    llr  = llbm*(-0.5/((σλ*srδt)^2))
    llr -= δt*llpb

    # add final non-standard `δt`
    δtf   = t[nI+2] - t[nI+1]
    lλpi1 = lλp[nI+2]
    lλci1 = lλc[nI+2]
    llr += lrdnorm_bm_x(lλpi1, lλpi, lλci1, lλci, sqrt(δtf)*σλ)
    llr -= δtf*0.5*((elλp1 + exp(lλpi1)) - (elλc1 + exp(lλci1)))
  end

  return llr
end




"""
    ll_gbm_b_sep(t   ::Array{Float64,1},
                 lλv ::Array{Float64,1},
                 σλ  ::Float64, 
                 δt  ::Float64,
                 srδt::Float64)

Returns the log-likelihood for a branch according to GBM pure-birth 
separately for the Brownian motion and the pure-birth
"""
function ll_gbm_b_sep(t   ::Array{Float64,1},
                      lλv ::Array{Float64,1},
                      σλ  ::Float64, 
                      δt  ::Float64,
                      srδt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    lλvi = lλv[1]
    llbm = 0.0
    llpb = 0.5*exp(lλvi)
    @simd for i in Base.OneTo(nI-1)
      lλvi1 = lλv[i+1]
      llbm += (lλvi1 - lλvi)^2
      llpb += exp(lλvi1)
      lλvi  = lλvi1
    end

    # last interval for the end 0.5 in pure-birth
    lλvi1 = lλv[nI+1]
    llbm += (lλvi1 - lλvi)^2
    elλv1 = exp(lλvi1)
    llpb += 0.5*elλv1
    lλvi  = lλvi1

    # add to global likelihood
    llbm *= (-0.5/((σλ*srδt)^2))
    llbm -= Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
    llpb *= (-δt)

    # add final non-standard `δt`
    δtf   = t[nI+2] - t[nI+1]
    lλvi1 = lλv[nI+2]
    llbm += ldnorm_bm(lλvi1, lλvi, sqrt(δtf)*σλ)
    llpb -= δtf*0.5*(elλv1 + exp(lλvi1))
  end

  return llbm, llpb
end




"""
    llr_gbm_b_sep(t   ::Array{Float64,1},
                 lλv ::Array{Float64,1},
                 σλ  ::Float64, 
                 δt  ::Float64,
                 srδt::Float64)

Returns the log-likelihood for a branch according to GBM pure-birth 
separately for the Brownian motion and the pure-birth
"""
function llr_gbm_b_sep(t   ::Array{Float64,1},
                       lλp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       σλ  ::Float64, 
                       δt  ::Float64,
                       srδt::Float64)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    lλpi = lλp[1]
    lλci = lλc[1]
    llrbm = 0.0
    llrpb = 0.5*exp(lλpi) - 0.5*exp(lλci)
    @simd for i in Base.OneTo(nI-1)
      lλpi1 = lλp[i+1]
      lλci1 = lλc[i+1]
      llrbm += (lλpi1 - lλpi)^2 - (lλci1 - lλci)^2
      llrpb += exp(lλpi1) - exp(lλci1)
      lλpi  = lλpi1
      lλci  = lλci1
    end

    # last interval for the end 0.5 in pure-birth
    lλpi1 = lλp[nI+1]
    lλci1 = lλc[nI+1]
    llrbm += (lλpi1 - lλpi)^2 - (lλci1 - lλci)^2
    elλp1 = exp(lλpi1)
    elλc1 = exp(lλci1)
    llrpb += 0.5*elλp1 - 0.5*elλc1
    lλpi  = lλpi1
    lλci  = lλci1

    # add to global likelihood
    llrbm *= (-0.5/((σλ*srδt)^2))
    llrpb *= (-δt)

    # add final non-standard `δt`
    δtf   = t[nI+2] - t[nI+1]
    lλpi1 = lλp[nI+2]
    lλci1 = lλc[nI+2]
    llrbm += lrdnorm_bm_x(lλpi1, lλpi, lλci1, lλci, sqrt(δtf)*σλ)
    llrpb -= δtf*0.5*((elλp1 + exp(lλpi1)) - (elλc1 + exp(lλci1)))
  end

  return llrbm, llrpb
end




"""
    ll_gbm_b_pb(t  ::Array{Float64,1},
                lλv::Array{Float64,1},
                δt ::Float64)

Returns the log-likelihood for a branch according to GBM pure-birth 
separately for the Brownian motion and the pure-birth
"""
function ll_gbm_b_pb(t  ::Array{Float64,1},
                     lλv::Array{Float64,1},
                     δt ::Float64)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    lλvi = lλv[1]
    ll = 0.5*exp(lλvi)
    @simd for i in Base.OneTo(nI-1)
      lλvi1 = lλv[i+1]
      ll   += exp(lλvi1)
      lλvi  = lλvi1
    end

    # last interval for the end 0.5 in pure-birth
    lλvi1 = lλv[nI+1]
    elλv1 = exp(lλvi1)
    ll   += 0.5*elλv1
    lλvi  = lλvi1

    # add to global likelihood
    ll *= (-δt)

    # add final non-standard `δt`
    ll -= (t[nI+2] - t[nI+1])*0.5*(elλv1 + exp(lλv[nI+2]))
  end

  return ll
end




"""
    llr_gbm_bm(tree::iTgbmpb, 
               σλp  ::Float64,
               σλc  ::Float64,
               srδt::Float64)

Returns the log-likelihood for a `iTgbmpb` according to GBM birth-death.
"""
function llr_gbm_bm(tree::iTgbmpb, 
                    σλp  ::Float64,
                    σλc  ::Float64,
                    srδt::Float64)

  if istip(tree) 
    llr_gbm_bm(ts(tree), lλ(tree), σλp, σλc, srδt)
  else
    llr_gbm_bm(ts(tree), lλ(tree), σλp, σλc, srδt) +
    llr_gbm_bm(tree.d1::iTgbmpb, σλp, σλc, srδt)   +
    llr_gbm_bm(tree.d2::iTgbmpb, σλp, σλc, srδt)
  end
end




"""
    llr_gbm_bm(t   ::Array{Float64,1},
               lλv ::Array{Float64,1},
               σλp  ::Float64, 
               σλc  ::Float64, 
               srδt::Float64)

Returns the log-likelihood ratio for a branch according to BM for `σ` proposal.
"""
function llr_gbm_bm(t   ::Array{Float64,1},
                    lλv ::Array{Float64,1},
                    σλp  ::Float64, 
                    σλc  ::Float64, 
                    srδt::Float64)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    ss   = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      ss   += (lλvi1 - lλvi)^2
      lλvi  = lλvi1
    end

    # likelihood ratio
    llr = ss*(0.5/((σλc*srδt)^2) - 0.5/((σλp*srδt)^2)) - 
          Float64(nI)*(log(σλp/σλc))

    # add final non-standard `δt`
    srδtf = sqrt(t[nI+2] - t[nI+1])
    lλvi1 = lλv[nI+2]
    llr  += lrdnorm_bm_σ(lλvi1, lλvi, srδtf*σλp, srδtf*σλc)
  end

  return llr
end





