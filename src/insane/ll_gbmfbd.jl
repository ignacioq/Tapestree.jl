#=

`fbdd` likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 10 03 2022
=#




"""
    llik_gbm(Ξ   ::Vector{iTfbd},
             idf ::Vector{iBffs},
             αλ  ::Float64,
             αμ  ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             ψ   ::Vector{Float64},
             ψts ::Vector{Float64},
             bst ::Vector{Float64},
             eix ::Vector{Int64},
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTfbd` according to `fbdd`.
"""
function llik_gbm(Ξ   ::Vector{iTfbd},
                  idf ::Vector{iBffs},
                  αλ  ::Float64,
                  αμ  ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  ψ   ::Vector{Float64},
                  ψts ::Vector{Float64},
                  bst ::Vector{Float64},
                  eix ::Vector{Int64},
                  δt  ::Float64,
                  srδt::Float64)
  @inbounds begin

    nep = lastindex(ψts) + 1
    ll  = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      lli, ix = 
        llik_gbm(Ξ[i], αλ, αμ, σλ, σμ, ψ, bst[i], ψts, eix[i], δt, srδt, nep)
      ll += lli
      if d2(bi) > 0
        ll += λt(bi)
      end
    end
  end

  return ll
end




"""
    llik_gbm(tree::iTfbd,
             αλ  ::Float64,
             αμ  ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             ψ   ::Vector{Float64},
             t   ::Float64,
             ψts ::Vector{Float64},
             ix  ::Int64,
             δt  ::Float64,
             srδt::Float64,
             nep ::Int64)

Returns the log-likelihood for a `iTfbd` according to `fbdd`.
"""
function llik_gbm(tree::iTfbd,
                  αλ  ::Float64,
                  αμ  ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  ψ   ::Vector{Float64},
                  t   ::Float64,
                  ψts ::Vector{Float64},
                  ix  ::Int64,
                  δt  ::Float64,
                  srδt::Float64,
                  nep ::Int64)

  if def1(tree)
    if def2(tree)
      ei = e(tree)
      ll, ix = 
        ll_gbm_b(lλ(tree), lμ(tree), αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep,
          δt, fdt(tree), srδt, true, false, false)

      ll1, ix1 = 
        llik_gbm(tree.d1, αλ, αμ, σλ, σμ, ψ, t - ei, ψts, ix, δt, srδt, nep)
      ll2, ix1 = 
        llik_gbm(tree.d2, αλ, αμ, σλ, σμ, ψ, t - ei, ψts, ix, δt, srδt, nep)

      ll += ll1 + ll2
    else
      ll, ix = ll_gbm_b(lλ(tree), lμ(tree), αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep,
        δt, fdt(tree), srδt, false, false, true)
      ll1, ix = llik_gbm(tree.d1, αλ, αμ, σλ, σμ, ψ, t - e(tree), ψts, ix, 
                  δt, srδt, nep)

      ll += ll1
    end
  else
    ll, ix = ll_gbm_b(lλ(tree), lμ(tree), αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep,
               δt, fdt(tree), srδt, false, isextinct(tree), isfossil(tree))
  end

  return ll, ix
end





"""
    ll_gbm_b(lλv ::Array{Float64,1},
             lμv ::Array{Float64,1},
             αλ  ::Float64,
             αμ  ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             ψ   ::Vector{Float64},
             t   ::Float64,
             ψts ::Vector{Float64},
             ix  ::Int64,
             δt  ::Float64,
             fdt ::Float64,
             srδt::Float64,
             λev ::Bool,
             μev ::Bool,
             ψev ::Bool)

Returns the log-likelihood for a branch according to `fbdd`.
"""
function ll_gbm_b(lλv ::Array{Float64,1},
                  lμv ::Array{Float64,1},
                  αλ  ::Float64,
                  αμ  ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  ψ   ::Vector{Float64},
                  t   ::Float64,
                  ψts ::Vector{Float64},
                  ix  ::Int64,
                  nep ::Int64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool,
                  μev ::Bool,
                  ψev ::Bool)

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ = llμ = llbd = 0.0
    @turbo for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi - αλ*δt)^2
      llμ  += (lμvi1 - lμvi - αμ*δt)^2
      llbd += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1))
    end

    # global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) +
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π)) -
         δt*llbd

    # ψ likelihood
    ψi  = ψ[ix]
    et  = ix < nep ? ψts[ix] : -Inf
    te  = max(t - nI*δt - fdt, 0.0)
    while t >= et > te
      ll -= ψi*(t - et)
      t   = et
      ix += 1
      ψi  = ψ[ix]
      et  = ix < nep ? ψts[ix] : -Inf
    end
    ll -= ψi*(t - te)

    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi  = lλv[nI+1]
      lμvi  = lμv[nI+1]
      srfdt = sqrt(fdt)
      ll   += ldnorm_bm(lλvi1, lλvi + αλ*fdt, srfdt*σλ)               +
              ldnorm_bm(lμvi1, lμvi + αμ*fdt, srfdt*σμ)               -
              fdt*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)))
    end

    #if speciation
    if λev
      ll += lλvi1
    # if extinction
    elseif μev
      ll += lμvi1
    # if fossilization
    elseif ψev
      ll += log(ψi)
    end
  end

  return ll, ix
end




"""
    llik_gbm_ss(tree::iTfbd,
                αλ  ::Float64,
                αμ  ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                t   ::Float64,
                ψts ::Vector{Float64},
                ix  ::Int64,
                δt  ::Float64,
                srδt::Float64,
                nep ::Int64,
                ns  ::Float64,
                ne  ::Float64)

Returns the log-likelihood for a `iTfbd` according to `fbdd`.
"""
function llik_gbm_ss(tree::iTfbd,
                     αλ  ::Float64,
                     αμ  ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     ψ   ::Vector{Float64},
                     t   ::Float64,
                     ψts ::Vector{Float64},
                     ix  ::Int64,
                     δt  ::Float64,
                     srδt::Float64,
                     nep ::Int64,
                     ns  ::Float64,
                     ne  ::Float64)

  if def1(tree)
    if def2(tree)
      ns += 1.0
      ei  = e(tree)
      ll, ix, ddλ, ddμ, ssλ, ssμ, nλ =
        ll_gbm_ss_b(lλ(tree), lμ(tree), αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep,
          δt, fdt(tree), srδt, true, false, false)

      ll1, ix1, ddλ1, ddμ1, ssλ1, ssμ1, nλ1, ns, ne =
        llik_gbm_ss(tree.d1, αλ, αμ, σλ, σμ, ψ, t - ei, ψts, ix, 
          δt, srδt, nep, ns, ne)
      ll2, ix1, ddλ2, ddμ2, ssλ2, ssμ2, nλ2, ns, ne =
        llik_gbm_ss(tree.d2, αλ, αμ, σλ, σμ, ψ, t - ei, ψts, ix, 
          δt, srδt, nep, ns, ne)

      ll  += ll1  + ll2
      ddλ += ddλ1 + ddλ2
      ddμ += ddμ1 + ddμ2
      ssλ += ssλ1 + ssλ2
      ssμ += ssμ1 + ssμ2
      nλ  += nλ1  + nλ2
    else
      ll, ix, ddλ, ddμ, ssλ, ssμ, nλ =
        ll_gbm_ss_b(lλ(tree), lμ(tree), αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep,
          δt, fdt(tree), srδt, false, false, true)

      ll1, ix, ddλ1, ddμ1, ssλ1, ssμ1, nλ1, ns, ne =
        llik_gbm_ss(tree.d1, αλ, αμ, σλ, σμ, ψ, t - e(tree), ψts, ix, 
          δt, srδt, nep, ns, ne)

      ll  += ll1
      ddλ += ddλ1
      ddμ += ddμ1
      ssλ += ssλ1
      ssμ += ssμ1
      nλ  += nλ1
    end
  else
    ie  = isextinct(tree)
    ne += Float64(ie)

    ll, ix, ddλ, ddμ, ssλ, ssμ, nλ =
      ll_gbm_ss_b(lλ(tree), lμ(tree), αλ, αμ, σλ, σμ, ψ, t, ψts, ix, nep,
        δt, fdt(tree), srδt, false, ie, isfossil(tree))
  end

  return ll, ix, ddλ, ddμ, ssλ, ssμ, nλ, ns, ne
end




"""
    ll_gbm_ss_b(lλv ::Array{Float64,1},
                lμv ::Array{Float64,1},
                αλ  ::Float64,
                αμ  ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                t   ::Float64,
                ψts ::Vector{Float64},
                ix  ::Int64,
                nep ::Int64,
                δt  ::Float64,
                fdt ::Float64,
                srδt::Float64,
                λev ::Bool,
                μev ::Bool,
                ψev ::Bool)

Returns the log-likelihood for a branch according to `fbdd`.
"""
function ll_gbm_ss_b(lλv ::Array{Float64,1},
                     lμv ::Array{Float64,1},
                     αλ  ::Float64,
                     αμ  ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     ψ   ::Vector{Float64},
                     t   ::Float64,
                     ψts ::Vector{Float64},
                     ix  ::Int64,
                     nep ::Int64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64,
                     λev ::Bool,
                     μev ::Bool,
                     ψev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ = llμ = llbdλ = llbdμ = 0.0
    @turbo for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi - αλ*δt)^2
      llμ  += (lμvi1 - lμvi - αμ*δt)^2
      llbdλ += exp(0.5*(lλvi + lλvi1))
      llbdμ += exp(0.5*(lμvi + lμvi1))
    end

    # standardized sum of squares
    ssλ = llλ/(2.0*δt)
    ssμ = llμ/(2.0*δt)
    nλ  = Float64(nI)

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) +
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    ll -= (llbdλ + llbdμ)*δt

    # ψ likelihood
    ψi  = ψ[ix]
    et  = ix < nep ? ψts[ix] : -Inf
    te  = max(t - nI*δt - fdt, 0.0)
    while t >= et > te
      ll -= ψi*(t - et)
      t   = et
      ix += 1
      ψi  = ψ[ix]
      et  = ix < nep ? ψts[ix] : -Inf
    end
    ll -= ψi*(t - te)

    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]
    ddλ = lλvi1 - lλv[1]
    ddμ = lμvi1 - lμv[1]

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi  = lλv[nI+1]
      lμvi  = lμv[nI+1]
      srfdt = sqrt(fdt)
      ll   += ldnorm_bm(lλvi1, lλvi + αλ*fdt, srfdt*σλ) +
              ldnorm_bm(lμvi1, lμvi + αμ*fdt, srfdt*σμ)
      ssλ  += (lλvi1 - lλvi - αλ*fdt)^2/(2.0*fdt)
      ssμ  += (lμvi1 - lμvi - αμ*fdt)^2/(2.0*fdt)
      nλ   += 1.0
      ll   -= fdt*(exp(0.5*(lλvi + lλvi1))) + 
              fdt*(exp(0.5*(lμvi + lμvi1)))
    end

      #if speciation
    if λev
      ll += lλvi1
    # if extinction
    elseif μev
      ll += lμvi1
    # if fossilization
    elseif ψev
      ll += log(ψi)
    end
  end

  return ll, ix, ddλ, ddμ, ssλ, ssμ, nλ
end




"""
    llr_gbm_b_sep(lλp ::Array{Float64,1},
                  lμp ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  lμc ::Array{Float64,1},
                  αλ  ::Float64,
                  αμ  ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool,
                  μev ::Bool)

Returns the log-likelihood for a branch according to `fbdd`
separately (for gbm and bd).
"""
function llr_gbm_b_sep(lλp ::Array{Float64,1},
                       lμp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       lμc ::Array{Float64,1},
                       αλ  ::Float64,
                       αμ  ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       srδt::Float64,
                       λev ::Bool,
                       μev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλc)-2

    llrbmλ = llrbmμ = llrbdλ = llrbdμ = 0.0
    @turbo for i in Base.OneTo(nI)
      lλpi    = lλp[i]
      lλci    = lλc[i]
      lμpi    = lμp[i]
      lμci    = lμc[i]
      lλpi1   = lλp[i+1]
      lλci1   = lλc[i+1]
      lμpi1   = lμp[i+1]
      lμci1   = lμc[i+1]
      llrbmλ += (lλpi1 - lλpi - αλ*δt)^2 - (lλci1 - lλci - αλ*δt)^2
      llrbmμ += (lμpi1 - lμpi - αμ*δt)^2 - (lμci1 - lμci - αμ*δt)^2
      llrbdλ += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
      llrbdμ += exp(0.5*(lμpi + lμpi1)) - exp(0.5*(lμci + lμci1))
    end

    # standardized sum of squares
    ssrλ = llrbmλ/(2.0*δt)
    ssrμ = llrbmμ/(2.0*δt)

    # overall
    llrbmλ *= (-0.5/((σλ*srδt)^2))
    llrbmμ *= (-0.5/((σμ*srδt)^2))
    llrbm   = llrbmλ + llrbmμ
    llrbdλ *= (-δt)
    llrbdμ *= (-δt)

    lλpi1 = lλp[nI+2]
    lμpi1 = lμp[nI+2]
    lλci1 = lλc[nI+2]
    lμci1 = lμc[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lλpi    = lλp[nI+1]
      lλci    = lλc[nI+1]
      lμpi    = lμp[nI+1]
      lμci    = lμc[nI+1]
      ssrλ   += ((lλpi1 - lλpi - αλ*fdt)^2 - 
                 (lλci1 - lλci - αλ*fdt)^2)/(2.0*fdt)
      ssrμ   += ((lμpi1 - lμpi - αμ*fdt)^2 - 
                 (lμci1 - lμci - αμ*fdt)^2)/(2.0*fdt)
      srfdt   = sqrt(fdt)
      llrbm  += lrdnorm_bm_x(lλpi1, lλpi + αλ*fdt,
                             lλci1, lλci + αλ*fdt, srfdt*σλ) +
                lrdnorm_bm_x(lμpi1, lμpi + αμ*fdt, 
                             lμci1, lμci + αμ*fdt, srfdt*σμ)
      llrbdλ -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
      llrbdμ -= fdt*(exp(0.5*(lμpi + lμpi1)) - exp(0.5*(lμci + lμci1)))
    end
    llrbd = llrbdλ + llrbdμ
    if λev
      llrbd += lλpi1 - lλci1
    elseif μev
      llrbd += lμpi1 - lμci1
    end
  end

  return llrbm, llrbd, ssrλ, ssrμ
end




"""
    _ss_dd(tree::T,
           αλ  ::Float64,
           αμ  ::Float64,
           ddλ ::Float64,
           ddμ ::Float64,
           ssλ ::Float64,
           ssμ ::Float64,
           n   ::Float64) where {T <: iTfbd}

Returns the standardized sum of squares for rate `v`, the path number `n` 
and the delta drift `dd`.
"""
function _ss_dd(tree::T,
                αλ  ::Float64,
                αμ  ::Float64,
                ddλ ::Float64,
                ddμ ::Float64,
                ssλ ::Float64,
                ssμ ::Float64,
                n   ::Float64) where {T <: iTfbd}

  ddλ0, ddμ0, ssλ0, ssμ0, n0 = 
    _ss_dd_b(lλ(tree), lμ(tree), αλ, αμ, dt(tree), fdt(tree))

  ddλ += ddλ0
  ddμ += ddμ0
  ssλ += ssλ0
  ssμ += ssμ0
  n   += n0

  if def1(tree)
    ddλ, ddμ, ssλ, ssμ, n = 
      _ss_dd(tree.d1, αλ, αμ, ddλ, ddμ, ssλ, ssμ, n)
    if def2(tree)
      ddλ, ddμ, ssλ, ssμ, n = 
        _ss_dd(tree.d2, αλ, αμ, ddλ, ddμ, ssλ, ssμ, n)
    end
  end

  return ddλ, ddμ, ssλ, ssμ, n
end




"""
    _ss_dd_b(lλv::Array{Float64,1},
             lμv::Array{Float64,1},
             αλ  ::Float64,
             αμ  ::Float64,
             δt ::Float64,
             fdt::Float64)

Returns the standardized sum of squares for rate `v`, the path number `n` 
and the delta drift `dd`.
"""
function _ss_dd_b(lλv::Array{Float64,1},
                  lμv::Array{Float64,1},
                  αλ  ::Float64,
                  αμ  ::Float64,
                  δt ::Float64,
                  fdt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ = ssμ = 0.0
    @turbo for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      ssλ  += (lλvi1 - lλvi - αλ*δt)^2
      ssμ  += (lμvi1 - lμvi - αμ*δt)^2
    end

    # standardize
    invt = 1.0/(2.0*δt)
    ssλ *= invt
    ssμ *= invt

    n = Float64(nI)
    # add final non-standard `δt`
    if fdt > 0.0
      invt = 1.0/(2.0*fdt)
      lλvi  = lλv[nI+1]
      lμvi  = lμv[nI+1]
      lλvi1 = lλv[nI+2]
      lμvi1 = lμv[nI+2]
      ssλ += invt * (lλvi1 - lλvi - αλ*fdt)^2
      ssμ += invt * (lμvi1 - lμvi - αμ*fdt)^2
      n += 1.0
    end
  end

  return (lλv[nI+2] - lλv[1]), (lμv[nI+2] - lμv[1]), ssλ, ssμ, n
end




"""
    _ss(tree::T, 
        αλ  ::Float64, 
        αμ  ::Float64,
        ssλ ::Float64, 
        ssμ ::Float64) where {T <: iTfbd}

Returns the standardized sum of squares for the gbm part of a branch
for `fbdd`.
"""
function _ss(tree::T, 
             αλ  ::Float64, 
             αμ  ::Float64,
             ssλ ::Float64, 
             ssμ ::Float64) where {T <: iTfbd}

  ssλ0, ssμ0 = _ss_b(lλ(tree), lμ(tree), αλ, αμ, dt(tree), fdt(tree))
  ssλ += ssλ0
  ssμ += ssμ0

  if def1(tree)
    ssλ, ssμ = _ss(tree.d1, αλ, αμ, ssλ, ssμ)
    if def2(tree)
      ssλ, ssμ = _ss(tree.d2, αλ, αμ, ssλ, ssμ)
    end
  end

  return ssλ, ssμ
end




"""
    _ss_b(lλv::Array{Float64,1},
          lμv::Array{Float64,1},
          αλ  ::Float64, 
          αμ  ::Float64,
          δt ::Float64,
          fdt::Float64)

Returns the standardized sum of squares for the gbm part of a branch
for `bdd`.
"""
function _ss_b(lλv::Array{Float64,1},
               lμv::Array{Float64,1},
               αλ  ::Float64, 
               αμ  ::Float64,
               δt ::Float64,
               fdt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    ssλ = ssμ = 0.0
    @turbo for i in Base.OneTo(nI)
      ssλ  += (lλv[i+1] - lλv[i] - αλ*δt)^2
      ssμ  += (lμv[i+1] - lμv[i] - αμ*δt)^2
    end

    # add to global likelihood
    invt = 1.0/(2.0*δt)
    ssλ *= invt
    ssμ *= invt

    # add final non-standard `δt`
    if fdt > 0.0
      invt = 1.0/(2.0*fdt)
      ssλ += invt * (lλv[nI+2] - lλv[nI+1] - αλ*fdt)^2
      ssμ += invt * (lμv[nI+2] - lμv[nI+1] - αμ*fdt)^2
    end
  end

  return ssλ, ssμ
end





