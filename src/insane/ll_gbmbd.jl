#=

GBM birth-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    cond_alone_events_crown(tree::iTgbmbd, tna::Float64, ll::Int64)

Condition events when there is only one alive lineage in the stem branch
to only be speciation events.
"""
function cond_alone_events_crown(tree::iTgbmbd)

  ll0 = cond_alone_events_ll(tree.d1::iTgbmbd, 0.0, 0.0)
  ll1 = cond_alone_events_ll(tree.d2::iTgbmbd, 0.0, 0.0)

  return (ll0 + ll1 - lλ(tree.d1::iTgbmbd)[1])
end






"""
    cond_alone_events_stem(tree::iTgbmbd,
                           dri ::BitArray{1},
                           ldr ::Int64,
                           ix  ::Int64)

Returns gbm birth-death likelihood for whole branch `br`.
"""
function cond_alone_events_stem(tree::iTgbmbd,
                                dri ::BitArray{1},
                                ldr ::Int64,
                                ix  ::Int64)

  if ix === ldr
    return cond_alone_events_ll(tree, 0.0, 0.0)
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTgbmbd)
    if ifx1 && isfix(tree.d2::iTgbmbd)
      ix += 1
      if dri[ix]
        cond_alone_events_stem(tree.d1::iTgbmbd, dri, ldr, ix)
      else
        cond_alone_events_stem(tree.d2::iTgbmbd, dri, ldr, ix)
      end
    elseif ifx1
      cond_alone_events_stem(tree.d1::iTgbmbd, dri, ldr, ix)
    else
      cond_alone_events_stem(tree.d2::iTgbmbd, dri, ldr, ix)
    end
  end

end





"""
    cond_alone_events_stem(tree::iTgbmbd, tna::Float64, ll::Int64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_alone_events_stem(tree::iTgbmbd) = 
  cond_alone_events_ll(tree, 0.0, 0.0)





#"""
#    cond_alone_events_tip_ll(tree::iTgbmbd, tna::Float64, ll::Int64)

#Condition events when there is only one alive lineage in the crown subtrees 
#to only be speciation events.
#"""
#function cond_alone_events_tip_ll(tree::iTgbmbd, tna::Float64, ll::Float64)

#  if tna < pe(tree)
#    @inbounds begin
#      lλv = lλ(tree)
#      lv  = lastindex(lλv)
#      λi  = lλv[lv]
#      μi  = lμ(tree)[lv]
#    end
#    ll += log(exp(λi) + exp(μi)) - λi
#  end
#  tna -= pe(tree)

#  if istip(tree)
#    return ll
#  end

#  if isfix(tree.d1::iTgbmbd)
#    if isfix(tree.d2::iTgbmbd)
#      return ll
#    else
#      tnx = treeheight(tree.d2::iTgbmbd)
#      tna = tnx > tna ? tnx : tna
#      cond_alone_events_ll(tree.d1::iTgbmbd, tna, ll)
#    end
#  else
#    tnx = treeheight(tree.d1::iTgbmbd)
#    tna = tnx > tna ? tnx : tna
#    cond_alone_events_ll(tree.d2::iTgbmbd, tna, ll)
#  end
#end




"""
    cond_alone_events_ll(tree::iTgbmbd, tna::Float64, ll::Int64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function cond_alone_events_ll(tree::iTgbmbd, tna::Float64, ll::Float64)

  if tna < pe(tree)
    @inbounds begin
      lλv = lλ(tree)
      lv  = lastindex(lλv)
      λi  = lλv[lv]
      μi  = lμ(tree)[lv]
    end
    ll += log(exp(λi) + exp(μi)) - λi
  end
  tna -= pe(tree)

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmbd)
    if isfix(tree.d2::iTgbmbd)
      return ll
    else
      tnx = treeheight(tree.d2::iTgbmbd)
      tna = tnx > tna ? tnx : tna
      cond_alone_events_ll(tree.d1::iTgbmbd, tna, ll)
    end
  else
    tnx = treeheight(tree.d1::iTgbmbd)
    tna = tnx > tna ? tnx : tna
    cond_alone_events_ll(tree.d2::iTgbmbd, tna, ll)
  end
end




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

  lλb  = lλ(tree)
  lμb  = lμ(tree)

  if istip(tree) 
    ll_gbm_b(lλb, lμb, σλ, σμ, δt, fdt(tree), srδt) + 
    (isextinct(tree) ? lμb[end] : 0.0)
  else
    ll_gbm_b(lλb, lμb, σλ, σμ, δt, fdt(tree), srδt) + lλb[end] + 
    llik_gbm(tree.d1::iTgbmbd, σλ, σμ, δt, srδt)    + 
    llik_gbm(tree.d2::iTgbmbd, σλ, σμ, δt, srδt)
  end
end




"""
    ll_gbm_b(lλv ::Array{Float64,1},
             lμv ::Array{Float64,1},
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64, 
             fdt::Float64,
             srδt::Float64)

Returns the log-likelihood for a branch according to GBM birth-death.
"""
@inline function ll_gbm_b(lλv ::Array{Float64,1},
                          lμv ::Array{Float64,1},
                          σλ  ::Float64,
                          σμ  ::Float64,
                          δt  ::Float64, 
                          fdt ::Float64,
                          srδt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

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
    srfdt = sqrt(fdt)
    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]
    ll += ldnorm_bm(lλvi1, lλvi, srfdt*σλ)
    ll += ldnorm_bm(lμvi1, lμvi, srfdt*σμ)
    ll -= fdt*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)))
  end

  return ll
end




"""
    ll_gbm_b_bm(lλv ::Array{Float64,1},
                lμv ::Array{Float64,1},
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64, 
                fdt ::Float64,
                srδt::Float64)

Returns the log-likelihood for the GBM part of a branch for GBM birth-death.
"""
@inline function ll_gbm_b_bm(lλv ::Array{Float64,1},
                             lμv ::Array{Float64,1},
                             σλ  ::Float64,
                             σμ  ::Float64,
                             δt  ::Float64, 
                             fdt ::Float64,
                             srδt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ  = 0.0
    llμ  = 0.0
    lλvi = lλv[1]
    lμvi = lμv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi)^2
      llμ  += (lμvi1 - lμvi)^2
      lλvi  = lλvi1
      lμvi  = lμvi1
    end

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) + 
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))

    # add final non-standard `δt`
    srfdt = sqrt(fdt)
    ll   += ldnorm_bm(lλv[nI+2], lλvi, srfdt*σλ) + 
            ldnorm_bm(lμv[nI+2], lμvi, srfdt*σμ)
  end

  return ll
end






"""
    ll_gbm_b_bd(lλv ::Array{Float64,1},
                lμv ::Array{Float64,1},
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64, 
                fdt ::Float64,
                srδt::Float64)


Returns the log-likelihood for the birth-death part of a branch for GBM 
birth-death.
"""
function ll_gbm_b_bd(lλv ::Array{Float64,1},
                     lμv ::Array{Float64,1},
                     σλ  ::Float64,
                     σμ  ::Float64,
                     δt  ::Float64, 
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ll   = 0.0
    lλvi = lλv[1]
    lμvi = lμv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      ll   += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)) 
      lλvi  = lλvi1
      lμvi  = lμvi1
    end

    # global likelihood
    ll *= (-δt)

    # add final non-standard `δt`
    ll -= fdt * (exp(0.5*(lλvi + lλv[nI+2])) + exp(0.5*(lμvi + lμv[nI+2])))
  end

  return ll
end





"""
    ll_gbm_b_bm(lλv ::Array{Float64,1},
                lμv ::Array{Float64,1},
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64, 
                fdt ::Float64,
                srδt::Float64)

Returns the log-likelihood for the GBM part of a branch for GBM birth-death.
"""
@inline function ll_gbm_b_bm(lλv ::Array{Float64,1},
                             lμv ::Array{Float64,1},
                             σλ  ::Float64,
                             σμ  ::Float64,
                             δt  ::Float64, 
                             fdt ::Float64,
                             srδt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ  = 0.0
    llμ  = 0.0
    lλvi = lλv[1]
    lμvi = lμv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi)^2
      llμ  += (lμvi1 - lμvi)^2
      lλvi  = lλvi1
      lμvi  = lμvi1
    end

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) + 
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))

    # add final non-standard `δt`
    srfdt = sqrt(fdt)
    ll   += ldnorm_bm(lλv[nI+2], lλvi, srfdt*σλ) + 
            ldnorm_bm(lμv[nI+2], lμvi, srfdt*σμ)
  end

  return ll
end



"""
    llr_gbm_bm(tree::iTgbmbd, 
               σp  ::Float64,
               σc  ::Float64,
               srδt::Float64,
               lf  ::Function)

Returns the log-likelihood ratio for a `iTgbmbd` according 
to GBM birth-death for a `σ` proposal.
"""
function llr_gbm_bm(tree::iTgbmbd, 
                    σp  ::Float64,
                    σc  ::Float64,
                    srδt::Float64,
                    lf  ::Function)

  if istip(tree) 
    llr_gbm_bm(lf(tree), σp, σc, fdt(tree), srδt)
  else
    llr_gbm_bm(lf(tree), σp, σc, fdt(tree), srδt)   +
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

Estimate gbm birth-death likelihood for the tree in a fix branch.
"""
function llik_gbm_f(tree::iTgbmbd,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  lλb = lλ(tree)
  lμb = lμ(tree)

  if istip(tree)
    ll = ll_gbm_b(lλb, lμb, σλ, σμ, δt, fdt(tree), srδt)
  else
    ll = ll_gbm_b(lλb, lμb, σλ, σμ, δt, fdt(tree), srδt) + lλb[end]

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
    br_ll_gbm(tree::iTgbmbd,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
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




"""
    br_llr_gbm(treep::iTgbmbd,
               treec::iTgbmbd,
               σλ   ::Float64,
               σμ   ::Float64,
               δt   ::Float64,
               srδt ::Float64,
               dri  ::BitArray{1},
               ldr  ::Int64,
               ix   ::Int64)

Returns gbm birth-death likelihood ratio for whole branch `br`.
"""
function br_llr_gbm(treep::iTgbmbd,
                    treec::iTgbmbd,
                    σλ   ::Float64,
                    σμ   ::Float64,
                    δt   ::Float64,
                    srδt ::Float64,
                    dri  ::BitArray{1},
                    ldr  ::Int64,
                    ix   ::Int64)

  if ix === ldr
    llrbm, llrbd = 
      llr_gbm_sep_f(treep::iTgbmbd, treec::iTgbmbd, σλ, σμ, δt, srδt)
    return llrbm, llrbd
  elseif ix < ldr
    ifx1 = isfix(treec.d1::iTgbmbd)
    if ifx1 && isfix(treec.d2::iTgbmbd)
      ix += 1
      if dri[ix]
        llrbm, llrbd = 
          br_llr_gbm(treep.d1::iTgbmbd, treec.d1::iTgbmbd, 
            σλ, σμ, δt, srδt, dri, ldr, ix)
      else
        llrbm, llrbd = 
          br_llr_gbm(treep.d2::iTgbmbd, treec.d2::iTgbmbd, 
            σλ, σμ, δt, srδt, dri, ldr, ix)
      end
    elseif ifx1
      llrbm, llrbd = 
        br_llr_gbm(treep.d1::iTgbmbd, treec.d1::iTgbmbd, 
          σλ, σμ, δt, srδt, dri, ldr, ix)
    else
      llrbm, llrbd = 
        br_llr_gbm(treep.d2::iTgbmbd, treec.d2::iTgbmbd, 
          σλ, σμ, δt, srδt, dri, ldr, ix)
    end
  end

  return llrbm, llrbd
end




"""
    llr_gbm_sep_f(treep::iTgbmbd,
                  treec::iTgbmbd,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Returns the log-likelihood for a branch according to GBM birth-death 
separately (for gbm and bd).
"""
function llr_gbm_sep_f(treep::iTgbmbd,
                       treec::iTgbmbd,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  lλvp = lλ(treep)
  lμvp = lμ(treep)
  lλvc = lλ(treec)
  lμvc = lμ(treec)
  lv   = lastindex(lλvc)

  if istip(treec)
    llrbm, llrbd = 
      llr_gbm_b_sep(lλvp, lμvp, lλvc, lμvc, σλ, σμ, δt, fdt(treec), srδt)
    llrbd += (isextinct(treec) ? (lμvp[lv] - lμvc[lv]) : 0.0)
  else
    llrbm, llrbd = 
      llr_gbm_b_sep(lλvp, lμvp, lλvc, lμvc, σλ, σμ, δt, fdt(treec), srδt) 
    llrbd += lλvp[lv] - lλvc[lv]

    ifx1 = isfix(treec.d1)
    if ifx1 && isfix(treec.d2)
      return llrbm, llrbd
    elseif ifx1
      llrbm0, llrbd0 = 
        llr_gbm_sep_f(treep.d1::iTgbmbd, treec.d1::iTgbmbd, σλ, σμ, δt, srδt)
      llrbm1, llrbd1 = 
        llr_gbm_sep(treep.d2::iTgbmbd, treec.d2::iTgbmbd, σλ, σμ, δt, srδt)
      llrbm += llrbm0 + llrbm1
      llrbd += llrbd0 + llrbd1
    else
      llrbm0, llrbd0 = 
        llr_gbm_sep(treep.d1::iTgbmbd, treec.d1::iTgbmbd, σλ, σμ, δt, srδt)
      llrbm1, llrbd1 = 
        llr_gbm_sep_f(treep.d2::iTgbmbd, treec.d2::iTgbmbd, σλ, σμ, δt, srδt)
      llrbm += llrbm0 + llrbm1
      llrbd += llrbd0 + llrbd1
    end
  end

  return llrbm, llrbd
end





"""
    llr_gbm_sep(treep::iTgbmbd, 
                treec::iTgbmbd,
                σλ   ::Float64,
                σμ   ::Float64,
                δt   ::Float64,
                srδt::Float64)

Returns the log-likelihood ratio for a tree according to GBM birth-death 
separately (for gbm and bd).
"""
function llr_gbm_sep(treep::iTgbmbd, 
                     treec::iTgbmbd,
                     σλ   ::Float64,
                     σμ   ::Float64,
                     δt   ::Float64,
                     srδt::Float64)

  lλvp = lλ(treep)
  lμvp = lμ(treep)
  lλvc = lλ(treec)
  lμvc = lμ(treec)
  lv   = lastindex(lλvc)

  if istip(treec) 
    llrbm, llrbd = 
      llr_gbm_b_sep(lλvp, lμvp, lλvc, lμvc, σλ, σμ, δt, fdt(treec), srδt)
    llrbd += (isextinct(treec) ? (lμvp[lv] - lμvc[lv]) : 0.0)
  else
    llrbm, llrbd = 
      llr_gbm_b_sep(lλvp, lμvp, lλvc, lμvc, σλ, σμ, δt, fdt(treec), srδt) 
    llrbd += lλvp[lv] - lλvc[lv]

    llrbm0, llrbd0 = 
      llr_gbm_sep(treep.d1::iTgbmbd, treec.d1::iTgbmbd, σλ, σμ, δt, srδt) 
    llrbm1, llrbd1 = 
      llr_gbm_sep(treep.d2::iTgbmbd, treec.d2::iTgbmbd, σλ, σμ, δt, srδt)

    llrbm += llrbm0 + llrbm1
    llrbd += llrbd0 + llrbd1
  end
  
  return llrbm, llrbd
end




"""
    llr_gbm_b_sep(lλp ::Array{Float64,1},
                  lμp ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  lμc ::Array{Float64,1},
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64, 
                  fdt::Float64,
                  srδt::Float64)

Returns the log-likelihood for a branch according to GBM birth-death 
separately (for gbm and bd).
"""
function llr_gbm_b_sep(lλp ::Array{Float64,1},
                       lμp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       lμc ::Array{Float64,1},
                       σλ  ::Float64,
                       σμ  ::Float64,
                       δt  ::Float64, 
                       fdt::Float64,
                       srδt::Float64)

  @inbounds @fastmath begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλc)-2

    llrbmλ = 0.0
    llrbmμ = 0.0
    llrbd  = 0.0
 
    lλpi = lλp[1]
    lλci = lλc[1]
    lμpi = lμp[1]
    lμci = lμc[1]

    @simd for i in Base.OneTo(nI)
      lλpi1   = lλp[i+1]
      lλci1   = lλc[i+1]
      lμpi1   = lμp[i+1]
      lμci1   = lμc[i+1]
      llrbmλ += (lλpi1 - lλpi)^2 - (lλci1 - lλci)^2
      llrbmμ += (lμpi1 - lμpi)^2 - (lμci1 - lμci)^2
      llrbd  += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)) +
                exp(0.5*(lμpi + lμpi1)) - exp(0.5*(lμci + lμci1))
      lλpi    = lλpi1
      lλci    = lλci1
      lμpi    = lμpi1
      lμci    = lμci1
    end
    # overall
    llrbmλ *= (-0.5/((σλ*srδt)^2))
    llrbmμ *= (-0.5/((σμ*srδt)^2))
    llrbd  *= (-δt)

    # add final non-standard `δt`
    srfdt = sqrt(fdt)
    lλpi1 = lλp[nI+2]
    lμpi1 = lμp[nI+2]
    lλci1 = lλc[nI+2]
    lμci1 = lμc[nI+2]

    llrbm  = llrbmλ + lrdnorm_bm_x(lλpi1, lλpi, lλci1, lλci, srfdt*σλ) +
             llrbmμ + lrdnorm_bm_x(lμpi1, lμpi, lμci1, lμci, srfdt*σμ)
    llrbd -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)) +
                  exp(0.5*(lμpi + lμpi1)) - exp(0.5*(lμci + lμci1)))
  end

  return llrbm, llrbd
end





