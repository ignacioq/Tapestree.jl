#=

GBM birth-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    cond_alone_events_crown(tree::iTgbmce, tna::Float64, ll::Int64)

Condition events when there is only one alive lineage in the stem branch
to only be speciation events.
"""
function cond_alone_events_crown(tree::iTgbmce, μ::Float64)

  ll0 = cond_alone_events_ll(tree.d1::iTgbmce, 0.0, 0.0, μ)
  ll1 = cond_alone_events_ll(tree.d2::iTgbmce, 0.0, 0.0, μ)

  return (ll0 + ll1 - lλ(tree.d1::iTgbmce)[1])
end





"""
    cond_alone_events_stem(tree::iTgbmce,
                           dri ::BitArray{1},
                           ldr ::Int64,
                           ix  ::Int64)

Returns gbm birth-death likelihood for whole branch `br`.
"""
function cond_alone_events_stem(tree::iTgbmce,
                                μ   ::Float64,
                                dri ::BitArray{1},
                                ldr ::Int64,
                                ix  ::Int64)

  if ix === ldr
    return cond_alone_events_ll(tree, 0.0, 0.0, μ)
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTgbmce)
    if ifx1 && isfix(tree.d2::iTgbmce)
      ix += 1
      if dri[ix]
        cond_alone_events_stem(tree.d1::iTgbmce, μ, dri, ldr, ix)
      else
        cond_alone_events_stem(tree.d2::iTgbmce, μ, dri, ldr, ix)
      end
    elseif ifx1
      cond_alone_events_stem(tree.d1::iTgbmce, μ, dri, ldr, ix)
    else
      cond_alone_events_stem(tree.d2::iTgbmce, μ, dri, ldr, ix)
    end
  end

end



"""
    cond_alone_events_stem_woλ(tree::iTgbmce,
                               μ   ::Float64,
                               dri ::BitArray{1},
                               ldr ::Int64,
                               ix  ::Int64)

Returns gbm birth-death likelihood for whole branch `br`.
"""
function cond_alone_events_stem_woλ(tree::iTgbmce,
                                    μ   ::Float64,
                                    dri ::BitArray{1},
                                    ldr ::Int64,
                                    ix  ::Int64)

  if ix === ldr
    return cond_alone_events_ll_woλ(tree, 0.0, 0.0, μ)
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTgbmce)
    if ifx1 && isfix(tree.d2::iTgbmce)
      ix += 1
      if dri[ix]
        cond_alone_events_stem_woλ(tree.d1::iTgbmce, μ, dri, ldr, ix)
      else
        cond_alone_events_stem_woλ(tree.d2::iTgbmce, μ, dri, ldr, ix)
      end
    elseif ifx1
      cond_alone_events_stem_woλ(tree.d1::iTgbmce, μ, dri, ldr, ix)
    else
      cond_alone_events_stem_woλ(tree.d2::iTgbmce, μ, dri, ldr, ix)
    end
  end

end





"""
    cond_alone_events_stem_woλ(tree::iTgbmce, μ::Float64) = 

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_alone_events_stem_woλ(tree::iTgbmce, μ::Float64) = 
  cond_alone_events_ll_woλ(tree, 0.0, 0.0, μ)


"""
    cond_alone_events_ll_woλ(tree::iTgbmce, 
                             tna ::Float64, 
                             ll  ::Float64,
                             μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function cond_alone_events_ll_woλ(tree::iTgbmce, 
                                  tna ::Float64, 
                                  ll  ::Float64,
                                  μ   ::Float64)

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmce) && isfix(tree.d2::iTgbmce)
    return ll
  end

  if tna < pe(tree)
    λi  = lλ(tree)[end]
    ll += log((exp(λi) + μ)) - λi
  end
  tna -= pe(tree)

  if isfix(tree.d1::iTgbmce)
    if isfix(tree.d2::iTgbmce)
      return ll
    else
      tnx = treeheight(tree.d2::iTgbmce)
      tna = tnx > tna ? tnx : tna
      cond_alone_events_ll_woλ(tree.d1::iTgbmce, tna, ll, μ)
    end
  else
    tnx = treeheight(tree.d1::iTgbmce)
    tna = tnx > tna ? tnx : tna
    cond_alone_events_ll_woλ(tree.d2::iTgbmce, tna, ll, μ)
  end
end




"""
    cond_alone_events_stem_λ(tree::iTgbmce, μ::Float64) = 

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_alone_events_stem_λ(tree::iTgbmce, μ::Float64) = 
  cond_alone_events_ll_λ(tree, 0.0, 0.0, μ)




"""
    cond_alone_events_ll_λ(tree::iTgbmce, 
                           tna ::Float64, 
                           ll  ::Float64, 
                           μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function cond_alone_events_ll_λ(tree::iTgbmce, 
                                tna ::Float64, 
                                ll  ::Float64, 
                                μ   ::Float64)


  if tna < pe(tree)
    λi  = lλ(tree)[end]
    ll += log((exp(λi) + μ)) - λi
  end
  tna -= pe(tree)

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmce)
    if isfix(tree.d2::iTgbmce)
      return ll
    else
      tnx = treeheight(tree.d2::iTgbmce)
      tna = tnx > tna ? tnx : tna
      cond_alone_events_ll_λ(tree.d1::iTgbmce, tna, ll, μ)
    end
  else
    tnx = treeheight(tree.d1::iTgbmce)
    tna = tnx > tna ? tnx : tna
    cond_alone_events_ll_λ(tree.d2::iTgbmce, tna, ll, μ)
  end
end





"""
    cond_alone_events_stem(tree::iTgbmce, μ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_alone_events_stem(tree::iTgbmce, μ::Float64) = 
  cond_alone_events_ll(tree, 0.0, 0.0, μ)




"""
    cond_alone_events_ll(tree::iTgbmce, 
                         tna ::Float64, 
                         ll  ::Float64, 
                         μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function cond_alone_events_ll(tree::iTgbmce, 
                              tna ::Float64, 
                              ll  ::Float64, 
                              μ   ::Float64)

  if istip(tree)
    return ll
  end

  if tna < pe(tree)
    λi  = lλ(tree)[end]
    ll += log((exp(λi) + μ)) - λi
  end
  tna -= pe(tree)

  if isfix(tree.d1::iTgbmce)
    if isfix(tree.d2::iTgbmce)
      return ll
    else
      tnx = treeheight(tree.d2::iTgbmce)
      tna = tnx > tna ? tnx : tna
      cond_alone_events_ll(tree.d1::iTgbmce, tna, ll, μ)
    end
  else
    tnx = treeheight(tree.d1::iTgbmce)
    tna = tnx > tna ? tnx : tna
    cond_alone_events_ll(tree.d2::iTgbmce, tna, ll, μ)
  end
end







"""
    llik_gbm(tree::iTgbmce, 
             μ   ::Float64
             σλ  ::Float64, 
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmce` according to GBM birth-death.
"""
function llik_gbm(tree::iTgbmce, 
                  μ   ::Float64,
                  σλ  ::Float64, 
                  δt  ::Float64,
                  srδt::Float64)
  if istip(tree) 
    ll_gbm_b(lλ(tree), μ, σλ, δt, fdt(tree), srδt, false, isextinct(tree))
  else
    ll_gbm_b(lλ(tree), μ, σλ, δt, fdt(tree), srδt, true, false) +
    llik_gbm(tree.d1::iTgbmce, μ, σλ, δt, srδt)                 +
    llik_gbm(tree.d2::iTgbmce, μ, σλ, δt, srδt)
  end
end





"""
    ll_gbm_b(lλv ::Array{Float64,1},
             μ   ::Float64,
             σλ  ::Float64,
             δt  ::Float64, 
             fdt ::Float64,
             srδt::Float64,
             λev ::Bool,
             μev ::Bool)

Returns the log-likelihood for a branch according to GBM birth-death.
"""
@inline function ll_gbm_b(lλv ::Array{Float64,1},
                          μ   ::Float64,
                          σλ  ::Float64,
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
      llλ  += (lλvi1 - lλvi)^2
      llbd += exp(0.5*(lλvi + lλvi1))
      lλvi  = lλvi1
    end

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π))

    # add to global likelihood
    llbd += Float64(nI)*μ
    ll   -= llbd*δt

    # add final non-standard `δt`
    if !iszero(fdt)
      lλvi1 = lλv[nI+2]
      ll += ldnorm_bm(lλvi1, lλvi, sqrt(fdt)*σλ)
      ll -= fdt*(exp(0.5*(lλvi + lλvi1)) + μ)

      if λev
        ll += lλvi1
      elseif μev
        ll += log(μ)
      end
    end
  end

  return ll
end





"""
    sss_gbm(tree::iTgbmce)

Returns the log-likelihood ratio for a `iTgbmce` according 
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(tree::iTgbmce)

  if istip(tree) 
    ssλ, n = 
      sss_gbm_b(lλ(tree), dt(tree), fdt(tree))
  else
    ssλ, n = 
      sss_gbm_b(lλ(tree), dt(tree), fdt(tree))
    ssλ1, n1 = 
      sss_gbm(tree.d1::iTgbmce)
    ssλ2, n2 = 
      sss_gbm(tree.d2::iTgbmce)

    ssλ += ssλ1 + ssλ2
    n   += n1 + n2
  end

  return ssλ, n
end




"""
    sss_gbm_b(lλv::Array{Float64,1},
              lμv::Array{Float64,1},
              δt ::Float64, 
              fdt::Float64)

Returns the standardized sum of squares for the GBM part of a branch 
for GBM birth-death.
"""
@inline function sss_gbm_b(lλv::Array{Float64,1},
                           δt ::Float64, 
                           fdt::Float64)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ  = 0.0
    lλvi = lλv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      ssλ  += (lλvi1 - lλvi)^2
      lλvi  = lλvi1
    end

    # add to global sum
    invt = 1.0/(2.0*δt)
    ssλ *= invt

    # add final non-standard `δt`
    if !iszero(fdt)
      invt = 1.0/(2.0*fdt)
      ssλ += invt * (lλv[nI+2] - lλvi)^2
      n = Float64(nI + 1)
    else
      n = Float64(nI)
    end
  end

  return ssλ, n
end




"""
    llr_gbm_bm(tree::iTgbmce, 
               σp  ::Float64,
               σc  ::Float64,
               srδt::Float64,
               lf  ::Function)

Returns the log-likelihood ratio for a `iTgbmce` according 
to GBM birth-death for a `σ` proposal.
"""
function llr_gbm_bm(tree::iTgbmce, 
                    σp  ::Float64,
                    σc  ::Float64,
                    srδt::Float64,
                    lf  ::Function)

  if istip(tree) 
    llr_gbm_bm(lf(tree), σp, σc, fdt(tree), srδt)
  else
    llr_gbm_bm(lf(tree), σp, σc, fdt(tree), srδt)   +
    llr_gbm_bm(tree.d1::iTgbmce, σp, σc, srδt, lf) +
    llr_gbm_bm(tree.d2::iTgbmce, σp, σc, srδt, lf)
  end
end




"""
    llik_gbm_f(tree::iTgbmce,
               μ   ::Float64,
               σλ  ::Float64,
               δt  ::Float64,
               srδt::Float64)

Estimate gbm birth-death likelihood for the tree in a fix branch.
"""
function llik_gbm_f(tree::iTgbmce,
                    μ   ::Float64,
                    σλ  ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  if istip(tree)
    ll = ll_gbm_b(lλ(tree), μ, σλ, δt, fdt(tree), srδt, false, false)
  else
    ll = ll_gbm_b(lλ(tree), μ, σλ, δt, fdt(tree), srδt, true, false)

    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      return ll
    elseif ifx1
      ll += llik_gbm_f(tree.d1::iTgbmce, μ, σλ, δt, srδt) +
            llik_gbm(  tree.d2::iTgbmce, μ, σλ, δt, srδt)
    else
      ll += llik_gbm(  tree.d1::iTgbmce, μ, σλ, δt, srδt) + 
            llik_gbm_f(tree.d2::iTgbmce, μ, σλ, δt, srδt)
    end
  end

  return ll
end




"""
    br_ll_gbm(tree::iTgbmce,
              μ   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              dri ::BitArray{1},
              ldr ::Int64,
              ix  ::Int64)

Returns gbm birth-death likelihood for whole branch `br`.
"""
function br_ll_gbm(tree::iTgbmce,
                   μ   ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   dri ::BitArray{1},
                   ldr ::Int64,
                   ix  ::Int64)

  if ix === ldr
    return llik_gbm_f(tree::iTgbmce, μ, σλ, δt, srδt)
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTgbmce)
    if ifx1 && isfix(tree.d2::iTgbmce)
      ix += 1
      if dri[ix]
        ll = br_ll_gbm(tree.d1::iTgbmce, μ, σλ, δt, srδt, dri, ldr, ix)
      else
        ll = br_ll_gbm(tree.d2::iTgbmce, μ, σλ, δt, srδt, dri, ldr, ix)
      end
    elseif ifx1
      ll = br_ll_gbm(tree.d1::iTgbmce, μ, σλ, δt, srδt, dri, ldr, ix)
    else
      ll = br_ll_gbm(tree.d2::iTgbmce, μ, σλ, δt, srδt, dri, ldr, ix)
    end
  end

  return ll
end




"""
    br_llr_gbm(treep::iTgbmce,
               treec::iTgbmce,
               σλ   ::Float64,
               σμ   ::Float64,
               δt   ::Float64,
               srδt ::Float64,
               dri  ::BitArray{1},
               ldr  ::Int64,
               ix   ::Int64)

Returns gbm birth-death likelihood ratio for whole branch `br`.
"""
function br_llr_gbm(treep::iTgbmce,
                    treec::iTgbmce,
                    σλ   ::Float64,
                    σμ   ::Float64,
                    δt   ::Float64,
                    srδt ::Float64,
                    dri  ::BitArray{1},
                    ldr  ::Int64,
                    ix   ::Int64)

  if ix === ldr
    llrbm, llrbd = 
      llr_gbm_sep_f(treep::iTgbmce, treec::iTgbmce, σλ, σμ, δt, srδt)
    return llrbm, llrbd
  elseif ix < ldr
    ifx1 = isfix(treec.d1::iTgbmce)
    if ifx1 && isfix(treec.d2::iTgbmce)
      ix += 1
      if dri[ix]
        llrbm, llrbd = 
          br_llr_gbm(treep.d1::iTgbmce, treec.d1::iTgbmce, 
            σλ, σμ, δt, srδt, dri, ldr, ix)
      else
        llrbm, llrbd = 
          br_llr_gbm(treep.d2::iTgbmce, treec.d2::iTgbmce, 
            σλ, σμ, δt, srδt, dri, ldr, ix)
      end
    elseif ifx1
      llrbm, llrbd = 
        br_llr_gbm(treep.d1::iTgbmce, treec.d1::iTgbmce, 
          σλ, σμ, δt, srδt, dri, ldr, ix)
    else
      llrbm, llrbd = 
        br_llr_gbm(treep.d2::iTgbmce, treec.d2::iTgbmce, 
          σλ, σμ, δt, srδt, dri, ldr, ix)
    end
  end

  return llrbm, llrbd
end




"""
    llr_gbm_sep_f(treep::iTgbmce,
                  treec::iTgbmce,
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Returns the log-likelihood for a branch according to GBM birth-death 
separately (for gbm and bd).
"""
function llr_gbm_sep_f(treep::iTgbmce,
                       treec::iTgbmce,
                       σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  if istip(treec)
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), σλ, δt, fdt(treec), srδt, false)
  else
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), σλ, δt, fdt(treec), srδt, true) 

    ifx1 = isfix(treec.d1)
    if ifx1 && isfix(treec.d2)
      return llrbm, llrbd
    elseif ifx1
      llrbm0, llrbd0 = 
        llr_gbm_sep_f(treep.d1::iTgbmce, treec.d1::iTgbmce, σλ, δt, srδt)
      llrbm1, llrbd1 = 
        llr_gbm_sep(treep.d2::iTgbmce, treec.d2::iTgbmce, σλ, δt, srδt)
      llrbm += llrbm0 + llrbm1
      llrbd += llrbd0 + llrbd1
    else
      llrbm0, llrbd0 = 
        llr_gbm_sep(treep.d1::iTgbmce, treec.d1::iTgbmce, σλ, δt, srδt)
      llrbm1, llrbd1 = 
        llr_gbm_sep_f(treep.d2::iTgbmce, treec.d2::iTgbmce, σλ, δt, srδt)
      llrbm += llrbm0 + llrbm1
      llrbd += llrbd0 + llrbd1
    end
  end

  return llrbm, llrbd
end




"""
    llr_gbm_sep(treep::iTgbmce, 
                treec::iTgbmce,
                σλ   ::Float64,
                δt   ::Float64,
                srδt::Float64)

Returns the log-likelihood ratio for a tree according to GBM birth-death 
separately (for gbm and bd).
"""
function llr_gbm_sep(treep::iTgbmce, 
                     treec::iTgbmce,
                     σλ   ::Float64,
                     δt   ::Float64,
                     srδt::Float64)

  if istip(treec) 
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), σλ, δt, fdt(treec), srδt, false)
  else
    llrbm, llrbd = 
      llr_gbm_b_sep(lλ(treep), lλ(treec), σλ, δt, fdt(treec), srδt, true) 

    llrbm0, llrbd0 = 
      llr_gbm_sep(treep.d1::iTgbmce, treec.d1::iTgbmce, σλ, δt, srδt) 
    llrbm1, llrbd1 = 
      llr_gbm_sep(treep.d2::iTgbmce, treec.d2::iTgbmce, σλ, δt, srδt)

    llrbm += llrbm0 + llrbm1
    llrbd += llrbd0 + llrbd1
  end
  
  return llrbm, llrbd
end





"""
    llr_gbm_b_sep(lλp ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  σλ  ::Float64,
                  δt  ::Float64, 
                  fdt::Float64,
                  srδt::Float64,
                  λev ::Bool)

Returns the log-likelihood for a branch according to GBM birth-death 
separately (for gbm and bd).
"""
function llr_gbm_b_sep(lλp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       σλ  ::Float64,
                       δt  ::Float64, 
                       fdt::Float64,
                       srδt::Float64,
                       λev ::Bool)

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
      llrbm += (lλpi1 - lλpi)^2 - (lλci1 - lλci)^2
      llrbd += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
      lλpi   = lλpi1
      lλci   = lλci1
    end
    # overall
    llrbm *= (-0.5/((σλ*srδt)^2))
    llrbd *= (-δt)

    # add final non-standard `δt`
    if !iszero(fdt)
      lλpi1 = lλp[nI+2]
      lλci1 = lλc[nI+2]
      llrbm += lrdnorm_bm_x(lλpi1, lλpi, lλci1, lλci, sqrt(fdt)*σλ)
      llrbd -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))

      if λev
        llrbd += lλpi1 - lλci1 
      end
    end
  end

  return llrbm, llrbd
end





