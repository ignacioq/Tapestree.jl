#=

`gbmct` likelihood for constant turnover

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    cond_surv_crown(tree::iTgbmct, ϵ::Float64)

Condition events when there is only one alive lineage in the stem branch
to only be speciation events.
"""
cond_surv_crown(tree::iTgbmct, ϵ::Float64) = 
  sum_alone_stem(tree.d1, 0.0, 0.0, ϵ) +
  sum_alone_stem(tree.d2, 0.0, 0.0, ϵ) - 
  lλ(tree)[1]




"""
    cond_surv_stem(tree::iTgbmct, ϵ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_surv_stem(tree::iTgbmct, ϵ::Float64) = 
  sum_alone_stem(tree, 0.0, 0.0, ϵ)




"""
    sum_alone_stem(tree::iTgbmct, 
                   tna ::Float64, 
                   ll  ::Float64, 
                   ϵ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem(tree::iTgbmct, 
                        tna ::Float64, 
                        ll  ::Float64, 
                        ϵ   ::Float64)

  if istip(tree)
    return ll
  end

  if tna < e(tree)
    ll += log(1.0 + ϵ)
  end
  tna -= e(tree)

  if isfix(tree.d1::iTgbmct)
    if isfix(tree.d2::iTgbmct)
      return ll
    else
      tnx = treeheight(tree.d2::iTgbmct)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem(tree.d1::iTgbmct, tna, ll, ϵ)
    end
  else
    tnx = treeheight(tree.d1::iTgbmct)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem(tree.d2::iTgbmct, tna, ll, ϵ)
  end
end





"""
    cond_surv_stem_p(tree::iTgbmct, ϵ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_surv_stem_p(tree::iTgbmct, ϵ::Float64) = 
  sum_alone_stem_p(tree, 0.0, 0.0, ϵ)




"""
    sum_alone_stem_p(tree::iTgbmct, 
                     tna ::Float64, 
                     ll  ::Float64, 
                     ϵ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem_p(tree::iTgbmct, 
                          tna ::Float64, 
                          ll  ::Float64, 
                          ϵ   ::Float64)

  if tna < e(tree)
    ll += log(1.0 + ϵ)
  end
  tna -= e(tree)

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmct)
    if isfix(tree.d2::iTgbmct)
      return ll
    else
      tnx = treeheight(tree.d2::iTgbmct)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem_p(tree.d1::iTgbmct, tna, ll, ϵ)
    end
  else
    tnx = treeheight(tree.d1::iTgbmct)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem_p(tree.d2::iTgbmct, tna, ll, ϵ)
  end
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
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))

    # add to global likelihood
    ll   -= llbd*δt*(1.0 + ϵ)

    lλvi1 = lλv[nI+2]

    if fdt > 0.0
      ll += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ)

      if λev
        ll += log(fdt) + 0.5*(lλvi + lλvi1)
      elseif μev
        ll += 0.5*(lλvi + lλvi1) + log(ϵ * fdt)
      else
        ll -= fdt*exp(0.5*(lλvi + lλvi1))*(1.0 + ϵ)
      end
    elseif λev
      ll += lλvi1
    end
  end

  return ll
end




"""
    sλ_gbm(tree::iTgbmct)
Returns the sum of `λ` rates for a `iTgbmct` according 
to `gbmct` for a `ϵ` proposal.
"""
function sλ_gbm(tree::iTgbmct)

  if istip(tree)
    sλt = sλ_gbm_b(lλ(tree), dt(tree), fdt(tree), isextinct(tree))
  else
    sλt  = sλ_gbm_b(lλ(tree), dt(tree), fdt(tree), true)
    sλt += sλ_gbm(tree.d1::iTgbmct) + 
           sλ_gbm(tree.d2::iTgbmct)
  end

  return sλt
end




"""
    sλ_gbm_b(lλv::Array{Float64,1},
             δt ::Float64, 
             fdt::Float64)

Returns the sum of `λ` rates for a `iTgbmct` branch according 
to `gbmct` for a `ϵ` proposal.
"""
@inline function sλ_gbm_b(lλv::Array{Float64,1},
                          δt ::Float64, 
                          fdt::Float64,
                          ev ::Bool)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    sλt  = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      sλt  += exp(0.5*(lλvi + lλvi1))
      lλvi  = lλvi1
    end

    # add to global sum
    sλt *= δt

    # add final non-standard `δt`
    if fdt > 0.0 && !ev
      sλt += fdt * exp(0.5*(lλv[nI+2] + lλvi))
    end
  end

  return sλt
end




"""
    llik_gbm_f(tree::iTgbmct,
               α    ::Float64,
               σλ  ::Float64,
               ϵ   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Estimate `gbmct` likelihood for the tree in a fix branch.
"""
function llik_gbm_f(tree::iTgbmct,
                    α    ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  if istip(tree)
    ll = ll_gbm_b_ϵ(lλ(tree), α, σλ, ϵ, δt, fdt(tree), srδt, false, false)
  else
    ll = ll_gbm_b_ϵ(lλ(tree), α, σλ, ϵ, δt, fdt(tree), srδt, true, false)

    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      return ll
    elseif ifx1
      ll += llik_gbm_f(tree.d1, α, σλ, ϵ, δt, srδt) +
            llik_gbm(  tree.d2, α, σλ, ϵ, δt, srδt)
    else
      ll += llik_gbm(  tree.d1, α, σλ, ϵ, δt, srδt) + 
            llik_gbm_f(tree.d2, α, σλ, ϵ, δt, srδt)
    end
  end

  return ll
end




"""
    br_ll_gbm(tree::iTgbmct,
              α    ::Float64,
              σλ  ::Float64,
              ϵ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              dri ::BitArray{1},
              ldr ::Int64,
              ix  ::Int64)

Returns `gbmct` likelihood for whole branch `br`.
"""
function br_ll_gbm(tree::iTgbmct,
                   α    ::Float64,
                   σλ  ::Float64,
                   ϵ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   dri ::BitArray{1},
                   ldr ::Int64,
                   ix  ::Int64)

  if ix === ldr
    return llik_gbm_f(tree, α, σλ, ϵ, δt, srδt)
  elseif ix < ldr
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      ix += 1
      if dri[ix]
        ll = br_ll_gbm(tree.d1, α, σλ, ϵ, δt, srδt, dri, ldr, ix)
      else
        ll = br_ll_gbm(tree.d2, α, σλ, ϵ, δt, srδt, dri, ldr, ix)
      end
    elseif ifx1
      ll = br_ll_gbm(tree.d1, α, σλ, ϵ, δt, srδt, dri, ldr, ix)
    else
      ll = br_ll_gbm(tree.d2, α, σλ, ϵ, δt, srδt, dri, ldr, ix)
    end
  end

  return ll
end





"""
    llr_gbm_sep_f(treep::iTgbmct,
                  treec::iTgbmct,
                  α    ::Float64,
                  σλ   ::Float64,
                  ϵ    ::Float64,
                  δt   ::Float64,
                  srδt ::Float64)

Returns the log-likelihood for a branch according to `gbmct` 
separately (for gbm and bd).
"""
function llr_gbm_sep_f(treep::iTgbmct,
                       treec::iTgbmct,
                       α    ::Float64,
                       σλ   ::Float64,
                       ϵ    ::Float64,
                       δt   ::Float64,
                       srδt ::Float64)

  if istip(treec)
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), α, σλ, ϵ, δt, fdt(treec), srδt, 
        false, isextinct(treec))
  else
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), α, σλ, ϵ, δt, fdt(treec), srδt, 
        true, false) 

    ifx1 = isfix(treec.d1)
    if ifx1 && isfix(treec.d2)
      return llrbm, llrbd
    elseif ifx1
      llrbm0, llrbd0 = llr_gbm_sep_f(treep.d1, treec.d1, α, σλ, ϵ, δt, srδt)
      llrbm1, llrbd1 = llr_gbm_sep(  treep.d2, treec.d2, α, σλ, ϵ, δt, srδt)
      llrbm += llrbm0 + llrbm1
      llrbd += llrbd0 + llrbd1
    else
      llrbm0, llrbd0 = llr_gbm_sep(  treep.d1, treec.d1, α, σλ, ϵ, δt, srδt)
      llrbm1, llrbd1 = llr_gbm_sep_f(treep.d2, treec.d2, α, σλ, ϵ, δt, srδt)
      llrbm += llrbm0 + llrbm1
      llrbd += llrbd0 + llrbd1
    end
  end

  return llrbm, llrbd
end




"""
    llr_gbm_sep(treep::iTgbmct, 
                treec::iTgbmct,
                α    ::Float64,
                σλ   ::Float64,
                ϵ    ::Float64,
                δt   ::Float64,
                srδt::Float64)

Returns the log-likelihood ratio for a tree according to `gbmct` 
separately (for gbm and bd).
"""
function llr_gbm_sep(treep::iTgbmct, 
                     treec::iTgbmct,
                     α    ::Float64,
                     σλ   ::Float64,
                     ϵ    ::Float64,
                     δt   ::Float64,
                     srδt::Float64)

  if istip(treec) 
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), α, σλ, ϵ, δt, fdt(treec), srδt, 
        false, isextinct(treec))
  else
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), α, σλ, ϵ, δt, fdt(treec), srδt, 
        true, false) 

    llrbm0, llrbd0 = llr_gbm_sep(treep.d1, treec.d1, α, σλ, ϵ, δt, srδt) 
    llrbm1, llrbd1 = llr_gbm_sep(treep.d2, treec.d2, α, σλ, ϵ, δt, srδt)

    llrbm += llrbm0 + llrbm1
    llrbd += llrbd0 + llrbd1
  end
  
  return llrbm, llrbd
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
    # overall
    llrbm *= (-0.5/((σλ*srδt)^2))
    llrbd *= -δt*(1.0 + ϵ)

    lλpi1 = lλp[nI+2]
    lλci1 = lλc[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      srfdt = sqrt(fdt)
      llrbm += lrdnorm_bm_x(lλpi1, lλpi + α*fdt, 
                            lλci1, lλci + α*fdt, sqrt(fdt)*σλ)

      if λev
        llrbd += 0.5*(lλpi + lλpi1) - 0.5*(lλci + lλci1)
      elseif μev
        llrbd += 0.5*(lλpi + lλpi1) - 0.5*(lλci + lλci1)
      else
        llrbd -= fdt*(1.0 + ϵ)*
                 (exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
      end
    elseif λev
      llrbd += lλpi1 - lλci1
    end
  end

  return llrbm, llrbd
end





