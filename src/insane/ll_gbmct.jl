#=

`gbmct` likelihood for constant turnover

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    llik_cbd(xi::Vector{iTgbmct}, 
             idf::Vector{iBffs},
             α   ::Float64,
             σλ  ::Float64,
             ϵ   ::Float64, 
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmct` according to `gbmce`.
"""
function llik_gbm(xi::Vector{iTgbmct}, 
                  idf::Vector{iBffs},
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64, 
                  δt  ::Float64,
                  srδt::Float64)
  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(xi))
      bi  = idf[i]
      ll += llik_gbm(xi[i], α, σλ, ϵ, δt, srδt)
      if !it(bi)
        ll += λt(bi)
      end
    end
  end

  return ll
end




"""
    llik_gbm(tree::iTgbmct,
             α   ::Float64, 
             σλ  ::Float64, 
             ϵ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmct` according to `gbmct`.
"""
function llik_gbm(tree::iTgbmct,
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
@inline function ll_gbm_b_ϵ(lλv ::Array{Float64,1},
                            α   ::Float64,
                            σλ  ::Float64,
                            ϵ   ::Float64,
                            δt  ::Float64, 
                            fdt ::Float64,
                            srδt::Float64,
                            λev ::Bool,
                            μev ::Bool)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ  = 0.0
    llbd = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      llλ  += (lλvi1 - lλvi - α*δt)^2
      llbd += exp(0.5*(lλvi + lλvi1))
      lλvi  = lλvi1
    end

    # add to global likelihood
    ll  = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
    ll -= llbd*δt*(1.0 + ϵ)

    lλvi1 = lλv[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      ll += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
            fdt*exp(0.5*(lλvi + lλvi1))*(1.0 + ϵ)
    end
    # if speciation
    if λev
      ll += lλvi1
    end
    # if extinction
    if μev
      ll += lλvi1 + log(ϵ)
    end
  end

  return ll
end




"""
    llik_gbm_ssλ(tree::iTgbmct,
                 α   ::Float64, 
                 σλ  ::Float64, 
                 ϵ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

Returns the log-likelihood for a `iTgbmct` according to `gbmct`.
"""
function llik_gbm_ssλ(tree::iTgbmct,
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
@inline function ll_gbm_b_ϵ_ssλ(lλv ::Array{Float64,1},
                                α   ::Float64,
                                σλ  ::Float64,
                                ϵ   ::Float64,
                                δt  ::Float64, 
                                fdt ::Float64,
                                srδt::Float64,
                                λev ::Bool,
                                μev ::Bool)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llbm = 0.0
    llbd = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      llbm += (lλvi1 - lλvi - α*δt)^2
      llbd += exp(0.5*(lλvi + lλvi1))
      lλvi  = lλvi1
    end

    # standardized sum of squares
    ssλ = llbm/(2.0*δt)
    nλ  = Float64(nI)

    # add to global likelihood
    ll    = llbm * 
            (-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))
    llbd *= δt
    Σλ    = llbd
    ll   -= llbd*(1.0 + ϵ)

    lλvi1 = lλv[nI+2]

    dλ = lλvi1 - lλv[1]

    # add final non-standard `δt`
    if fdt > 0.0
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
  end

  return ll, dλ, ssλ, Σλ, nλ 
end




"""
    Σλ_gbm(tree::iTgbmct)
Returns the sum of `λ` rates for a `iTgbmct` according 
to `gbmct` for a `ϵ` proposal.
"""
function Σλ_gbm(tree::iTgbmct)

  if istip(tree)
    Σλ_gbm_b(lλ(tree), dt(tree), fdt(tree))
  else
    Σλ_gbm_b(lλ(tree), dt(tree), fdt(tree)) +
    Σλ_gbm(tree.d1::iTgbmct) + 
    Σλ_gbm(tree.d2::iTgbmct)
  end
end




"""
    Σλ_gbm_b(lλv::Array{Float64,1},
             δt ::Float64, 
             fdt::Float64)

Returns the sum of `λ` rates for a `iTgbmct` branch according 
to `gbmct` for a `ϵ` proposal.
"""
@inline function Σλ_gbm_b(lλv::Array{Float64,1},
                          δt ::Float64, 
                          fdt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    Σλ   = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      Σλ   += exp(0.5*(lλvi + lλvi1))
      lλvi  = lλvi1
    end

    # add to global sum
    Σλ *= δt

    # add final non-standard `δt`
    if fdt > 0.0
      Σλ += fdt * exp(0.5*(lλv[nI+2] + lλvi))
    end
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

  @inbounds @fastmath begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλc)-2

    llrbm = 0.0
    llrbd = 0.0
    lλpi = lλp[1]
    lλci = lλc[1]
    @simd for i in Base.OneTo(nI)
      lλpi1  = lλp[i+1]
      lλci1  = lλc[i+1]
      llrbm += (lλpi1 - lλpi - α*δt)^2 - (lλci1 - lλci - α*δt)^2
      llrbd += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
      lλpi   = lλpi1
      lλci   = lλci1
    end

    # standardized sum of squares
    ssrλ = llrbm/(2.0*δt)
    Σrλ  = llrbd * δt

    # overall
    llrbm *= (-0.5/((σλ*srδt)^2))
    llrbd *= -δt*(1.0 + ϵ)

    lλpi1 = lλp[nI+2]
    lλci1 = lλc[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      ssrλ  += ((lλpi1 - lλpi - α*fdt)^2 - (lλci1 - lλci - α*fdt)^2)/(2.0*fdt) 
      llri   = fdt * (exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
      Σrλ   += llri
      llrbm += lrdnorm_bm_x(lλpi1, lλpi + α*fdt, 
                            lλci1, lλci + α*fdt, sqrt(fdt)*σλ)
      llrbd -= (1.0 + ϵ) * llri
    end
    # if speciation or extinction
    if λev || μev
      llrbd += lλpi1 - lλci1
    end
  end

  return llrbm, llrbd, ssrλ, Σrλ
end
 



"""
    make_scond(idf::Vector{iBffs}, stem::Bool, ::Type{iTgbmct})

Return closure for log-likelihood for conditioning
"""
function make_scond(idf::Vector{iBffs}, stem::Bool, ::Type{iTgbmct})

  b1  = idf[1]
  d1i = d1(b1)
  d2i = d2(b1)

  if stem
    # for whole likelihood
    f = let d1i = d1i, d2i = d2i
      function (xi::Vector{iTgbmct}, ϵ::Float64, sns::NTuple{3,BitVector})
        sn1 = sns[1]
        cond_ll(xi[1], 0.0, ϵ, sn1, lastindex(sn1), 1)
      end
    end
    # for new proposal
    f0 = (xi::iTgbmct, ϵ::Float64, ter::Bool) -> 
            sum_alone_stem_p(xi, 0.0, 0.0, ϵ)
  else
    # for whole likelihood
    f = let d1i = d1i, d2i = d2i
      function (xi::Vector{iTgbmct}, ϵ::Float64, sns::NTuple{3,BitVector})

        sn2 = sns[2]
        sn3 = sns[3]

        cond_ll(xi[d1i], 0.0, ϵ, sn2, lastindex(sn2), 1) +
        cond_ll(xi[d2i], 0.0, ϵ, sn3, lastindex(sn3), 1) + 
        (1.0 + ϵ)
      end
    end
    # for new proposal
    f0 = function (xi::iTgbmct, ϵ::Float64, ter::Bool)
      if ter
        sum_alone_stem(  xi, 0.0, 0.0, ϵ)
      else
        sum_alone_stem_p(xi, 0.0, 0.0, ϵ)
      end
    end
  end

  return f, f0
end




"""
    cond_ll(tree::iTgbmct,
            ll  ::Float64, 
            ϵ   ::Float64,
            sn  ::BitVector,
            lsn ::Int64,
            ix  ::Int64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function cond_ll(tree::iTgbmct,
                 ll  ::Float64, 
                 ϵ   ::Float64,
                 sn  ::BitVector,
                 lsn ::Int64,
                 ix  ::Int64)

  if lsn >= ix
    if sn[ix]
      ll += log(1.0 + ϵ)
    end

    if lsn === ix
      return ll
    else
      ix += 1
      if isfix(tree.d1::iTgbmct)
        ll = cond_ll(tree.d1::iTgbmct, ll, ϵ, sn, lsn, ix)
      else
        ll = cond_ll(tree.d2::iTgbmct, ll, ϵ, sn, lsn, ix)
      end
    end
  end

  return ll
end



