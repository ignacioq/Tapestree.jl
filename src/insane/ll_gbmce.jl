#=

`gbmce` likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    cond_surv_crown(tree::iTgbmce, μ::Float64)

Condition events when there is only one alive lineage in the stem branch
to only be speciation events.
"""
cond_surv_crown(tree::iTgbmce, μ::Float64) = 
  cond_surv_stem(tree.d1, 0.0, 0.0, μ) +
  cond_surv_stem(tree.d2, 0.0, 0.0, μ) - 
  lλ(tree)[1]




"""
    cond_surv_stem(tree::iTgbmce, μ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_surv_stem(tree::iTgbmce, μ::Float64) = 
  sum_alone_stem(tree, 0.0, 0.0, μ)




"""
    sum_alone_stem(tree::iTgbmce, 
                   tna ::Float64, 
                   ll  ::Float64,
                   μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem(tree::iTgbmce, 
                        tna ::Float64, 
                        ll  ::Float64,
                        μ   ::Float64)

  if istip(tree)
    return ll
  end

  if tna < e(tree)
    λi  = lλ(tree)[end]
    ll += log((exp(λi) + μ)) - λi
  end
  tna -= e(tree)

  if isfix(tree.d1::iTgbmce)
    if isfix(tree.d2::iTgbmce)
      return ll
    else
      tnx = treeheight(tree.d2::iTgbmce)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem(tree.d1::iTgbmce, tna, ll, μ)
    end
  else
    tnx = treeheight(tree.d1::iTgbmce)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem(tree.d2::iTgbmce, tna, ll, μ)
  end
end







"""
    cond_surv_stem_p(tree::iTgbmce, μ::Float64) = 

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_surv_stem_p(tree::iTgbmce, μ::Float64) = 
  sum_alone_stem_p(tree, 0.0, 0.0, μ)




"""
    sum_alone_stem_p(tree::iTgbmce, 
                     tna ::Float64, 
                     ll  ::Float64, 
                     μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem_p(tree::iTgbmce, 
                          tna ::Float64, 
                          ll  ::Float64, 
                          μ   ::Float64)

  if tna < e(tree)
    λi  = lλ(tree)[end]
    ll += log((exp(λi) + μ)) - λi
  end
  tna -= e(tree)

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmce)
    if isfix(tree.d2::iTgbmce)
      return ll
    else
      tnx = treeheight(tree.d2::iTgbmce)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem_p(tree.d1::iTgbmce, tna, ll, μ)
    end
  else
    tnx = treeheight(tree.d1::iTgbmce)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem_p(tree.d2::iTgbmce, tna, ll, μ)
  end
end





"""
    llik_gbm(tree::iTgbmce, 
             α   ::Float64,
             σλ  ::Float64, 
             μ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmce` according to `gbmce`.
"""
function llik_gbm(tree::iTgbmce, 
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
    llbd += Float64(nI)*μ
    ll   -= llbd*δt

    lλvi1 = lλv[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      ll += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ)
      ll -= fdt*(exp(0.5*(lλvi + lλvi1)) + μ)
    end

    # if speciation
    if λev
      ll += lλvi1
    end
    # if extinction
    if μev
      ll += log(μ)
    end
  end

  return ll
end




"""
    llik_gbm_f(tree::iTgbmce,
               α   ::Float64,
               σλ  ::Float64,
               μ   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Estimate `gbmce` likelihood for the tree in a fix branch.
"""
function llik_gbm_f(tree::iTgbmce,
                    α   ::Float64,
                    σλ  ::Float64,
                    μ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  if istip(tree)
    ll = ll_gbm_b(lλ(tree), α, σλ, μ, δt, fdt(tree), srδt, false, false)
  else
    ll = ll_gbm_b(lλ(tree), α, σλ, μ, δt, fdt(tree), srδt, true, false)

    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      return ll
    elseif ifx1
      ll += llik_gbm_f(tree.d1, α, σλ, μ, δt, srδt) +
            llik_gbm(  tree.d2, α, σλ, μ, δt, srδt)
    else
      ll += llik_gbm(  tree.d1, α, σλ, μ, δt, srδt) + 
            llik_gbm_f(tree.d2, α, σλ, μ, δt, srδt)
    end
  end

  return ll
end




"""
    br_ll_gbm(tree::iTgbmce,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              dri ::BitArray{1},
              ldr ::Int64,
              ix  ::Int64)

Returns `gbmce` likelihood for whole branch `br`.
"""
function br_ll_gbm(tree::iTgbmce,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   dri ::BitArray{1},
                   ldr ::Int64,
                   ix  ::Int64)

  if ix === ldr
    return llik_gbm_f(tree, α, σλ, μ, δt, srδt)
  elseif ix < ldr
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      ix += 1
      if dri[ix]
        ll = br_ll_gbm(tree.d1, α, σλ, μ, δt, srδt, dri, ldr, ix)
      else
        ll = br_ll_gbm(tree.d2, α, σλ, μ, δt, srδt, dri, ldr, ix)
      end
    elseif ifx1
      ll = br_ll_gbm(tree.d1, α, σλ, μ, δt, srδt, dri, ldr, ix)
    else
      ll = br_ll_gbm(tree.d2, α, σλ, μ, δt, srδt, dri, ldr, ix)
    end
  end

  return ll
end





"""
    llr_gbm_sep_f(treep::iTgbmce,
                  treec::iTgbmce,
                  α    ::Float64,
                  σλ   ::Float64,
                  δt   ::Float64,
                  srδt ::Float64)

Returns the log-likelihood for a branch according to `gbmce` 
separately (for gbm and bd).
"""
function llr_gbm_sep_f(treep::iTgbmce,
                       treec::iTgbmce,
                       α    ::Float64,
                       σλ   ::Float64,
                       δt   ::Float64,
                       srδt ::Float64)

  if istip(treec)
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), α, σλ, δt, fdt(treec), srδt, false)
  else
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), α, σλ, δt, fdt(treec), srδt, true) 

    ifx1 = isfix(treec.d1)
    if ifx1 && isfix(treec.d2)
      return llrbm, llrbd
    elseif ifx1
      llrbm0, llrbd0 = 
        llr_gbm_sep_f(treep.d1, treec.d1, α, σλ, δt, srδt)
      llrbm1, llrbd1 = 
        llr_gbm_sep(treep.d2, treec.d2, α, σλ, δt, srδt)
      llrbm += llrbm0 + llrbm1
      llrbd += llrbd0 + llrbd1
    else
      llrbm0, llrbd0 = 
        llr_gbm_sep(treep.d1, treec.d1, α, σλ, δt, srδt)
      llrbm1, llrbd1 = 
        llr_gbm_sep_f(treep.d2, treec.d2, α, σλ, δt, srδt)
      llrbm += llrbm0 + llrbm1
      llrbd += llrbd0 + llrbd1
    end
  end

  return llrbm, llrbd
end




"""
    llr_gbm_sep(treep::iTgbmce, 
                treec::iTgbmce,
                α    ::Float64,
                σλ   ::Float64,
                δt   ::Float64,
                srδt::Float64)

Returns the log-likelihood ratio for a tree according to `gbmce` 
separately (for gbm and bd).
"""
function llr_gbm_sep(treep::iTgbmce, 
                     treec::iTgbmce,
                     α    ::Float64,
                     σλ   ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)

  if istip(treec) 
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), α, σλ, δt, fdt(treec), srδt, false)
  else
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), α, σλ, δt, fdt(treec), srδt, true) 

    llrbm0, llrbd0 = 
      llr_gbm_sep(treep.d1, treec.d1, α, σλ, δt, srδt) 
    llrbm1, llrbd1 = 
      llr_gbm_sep(treep.d2, treec.d2, α, σλ, δt, srδt)

    llrbm += llrbm0 + llrbm1
    llrbd += llrbd0 + llrbd1
  end
  
  return llrbm, llrbd
end




