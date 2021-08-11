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

    # add final non-standard `δt`
    if !iszero(fdt)
      lλvi1 = lλv[nI+2]
      ll += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ)
      ll -= fdt*exp(0.5*(lλvi + lλvi1))
      if λev
        ll += lλvi1
      end
    end
  end

  return ll
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

    # add to global likelihood
    llrbm *= (-0.5/((σλ*srδt)^2))
    llrpb *= (-δt)

    # add final non-standard `δt`
    if !iszero(fdt)
      lλpi1 = lλp[nI+2]
      lλci1 = lλc[nI+2]
      llrbm += lrdnorm_bm_x(lλpi1, lλpi + α*δt, 
                            lλci1, lλci + α*δt, sqrt(fdt)*σλ)
      llrpb -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
      if λev
        llrpb += lλpi1 - lλci1
      end
    end
  end

  return llrbm, llrpb
end




"""
    sss_gbm(tree::iTgbmpb, α::Float64)

Returns the log-likelihood ratio for a `iTgbmpb` according 
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(tree::iTgbmpb, α::Float64)

  ssλ, n = sss_gbm_b(lλ(tree), α, dt(tree), fdt(tree))

  if isdefined(tree, :d1) 
    ssλ1, n1 = 
      sss_gbm(tree.d1, α)
    ssλ2, n2 = 
      sss_gbm(tree.d2, α)

    ssλ += ssλ1 + ssλ2
    n   += n1 + n2
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
@inline function sss_gbm_b(lλv::Array{Float64,1},
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
    if !iszero(fdt)
      ssλ += 1.0/(2.0*fdt) * (lλv[nI+2] - lλvi - α*fdt)^2
      n = Float64(nI + 1)
    else
      n = Float64(nI)
    end
  end

  return ssλ, n
end



"""
    llr_gbm_σp(σλp::Float64,
               σλc::Float64,
               ssλ::Float64,
               n  ::Float64)

Returns the log-likelihood ratio according to GBM for `σ` proposal
after computing standard sum of squares `sss`.
"""
function llr_gbm_σp(σλp::Float64,
                    σλc::Float64,
                    ssλ::Float64,
                    n  ::Float64)

  llr = ssλ*(1.0/σλc^2 - 1.0/σλp^2) - n*(log(σλp/σλc))

  return llr
end





"""
    treelength(tree::T) where {T <: iTree}

Return the branch length sum of `tree`.
"""
function treelength(tree::T, l::Float64) where {T <: iTree}
  l += e(tree)
  if isdefined(tree, :d1)
    l = treelength(tree.d1, l)::Float64
    l = treelength(tree.d2, l)::Float64
  end
  return l
end




"""
    treelength_dλ(tree::iTgbmpb, dλ::Float64, l::Float64)

Returns the log-likelihood ratio for a `iTgbmpb` according 
to GBM birth-death for a `σ` proposal.
"""
function treelength_dλ(tree::iTgbmpb, dλ::Float64, l::Float64)

  lλv = lλ(tree)
  dλ += lλv[end] - lλv[1]
  l  += e(tree)

  if isdefined(tree, :d1) 
    dλ, l = treelength_dλ(tree.d1, dλ, l)
    dλ, l = treelength_dλ(tree.d2, dλ, l)
  end

  return dλ, l
end



