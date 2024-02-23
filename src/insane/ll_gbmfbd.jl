#=

`fbdd` likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 10 03 2022
=#





"""
    llik_gbm(Ξ   ::Vector{iTfbd},
             idf ::Vector{iBffs},
             α   ::Float64,
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
                  α   ::Float64,
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
      lli, ix = llik_gbm(Ξ[i], α, σλ, σμ, ψ, bst[i], ψts, eix[i], δt, srδt, nep)
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
             α   ::Float64,
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
                  α   ::Float64,
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
        ll_gbm_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, t, ψts, ix, nep,
          δt, fdt(tree), srδt, true, false, false)

      ll1, ix1 = llik_gbm(tree.d1, α, σλ, σμ, ψ, t - ei, ψts, ix, δt, srδt, nep)
      ll2, ix1 = llik_gbm(tree.d2, α, σλ, σμ, ψ, t - ei, ψts, ix, δt, srδt, nep)

      ll += ll1 + ll2
    else
      ll, ix = ll_gbm_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, t, ψts, ix, nep,
        δt, fdt(tree), srδt, false, false, true)
      ll1, ix = 
        llik_gbm(tree.d1, α, σλ, σμ, ψ, t - e(tree), ψts, ix, δt, srδt, nep)

      ll += ll1
    end
  else
    ll, ix = ll_gbm_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, t, ψts, ix, nep,
               δt, fdt(tree), srδt, false, isextinct(tree), isfossil(tree))
  end

  return ll, ix
end





"""
    ll_gbm_b(lλv ::Array{Float64,1},
             lμv ::Array{Float64,1},
             α   ::Float64,
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
                  α   ::Float64,
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
      llλ  += (lλvi1 - lλvi - α*δt)^2
      llμ  += (lμvi1 - lμvi)^2
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
      ll   += ldnorm_bm(lλvi1, lλvi + α*fdt, srfdt*σλ)                +
              ldnorm_bm(lμvi1, lμvi, srfdt*σμ)                        -
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
                α   ::Float64,
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
function llik_gbm_ss(tree::iTfbd,
                     α   ::Float64,
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
      ll, ix, ddλ, ssλ, ssμ, nλ, irλ, irμ =
        ll_gbm_ss_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, t, ψts, ix, nep,
          δt, fdt(tree), srδt, true, false, false)

      ll1, ix1, ddλ1, ssλ1, ssμ1, nλ1, irλ1, irμ1, ns, ne =
        llik_gbm_ss(tree.d1, α, σλ, σμ, ψ, t - ei, ψts, ix, 
          δt, srδt, nep, ns, ne)
      ll2, ix1, ddλ2, ssλ2, ssμ2, nλ2, irλ2, irμ2, ns, ne =
        llik_gbm_ss(tree.d2, α, σλ, σμ, ψ, t - ei, ψts, ix, 
          δt, srδt, nep, ns, ne)

      ll  += ll1  + ll2
      ddλ  += ddλ1  + ddλ2
      ssλ += ssλ1 + ssλ2
      ssμ += ssμ1 + ssμ2
      nλ  += nλ1  + nλ2
      irλ += irλ1 + irλ2
      irμ += irμ1 + irμ2
    else
      ll, ix, ddλ, ssλ, ssμ, nλ, irλ, irμ =
        ll_gbm_ss_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, t, ψts, ix, nep,
          δt, fdt(tree), srδt, false, false, true)

      ll1, ix, ddλ1, ssλ1, ssμ1, nλ1, irλ1, irμ1, ns, ne =
        llik_gbm_ss(tree.d1, α, σλ, σμ, ψ, t - e(tree), ψts, ix, 
          δt, srδt, nep, ns, ne)

      ll  += ll1
      ddλ  += ddλ1
      ssλ += ssλ1
      ssμ += ssμ1
      nλ  += nλ1
      irλ += irλ1
      irμ += irμ1
    end
  else
    ie  = isextinct(tree)
    ne += Float64(ie)

    ll, ix, ddλ, ssλ, ssμ, nλ, irλ, irμ =
      ll_gbm_ss_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, t, ψts, ix, nep,
        δt, fdt(tree), srδt, false, ie, isfossil(tree))
  end

  return ll, ix, ddλ, ssλ, ssμ, nλ, irλ, irμ, ns, ne
end




"""
    ll_gbm_ss_b(lλv ::Array{Float64,1},
                lμv ::Array{Float64,1},
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64,
                fdt ::Float64,
                srδt::Float64,
                λev ::Bool,
                μev ::Bool)

Returns the log-likelihood for a branch according to `fbdd`.
"""
function ll_gbm_ss_b(lλv ::Array{Float64,1},
                     lμv ::Array{Float64,1},
                     α   ::Float64,
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
      llλ  += (lλvi1 - lλvi - α*δt)^2
      llμ  += (lμvi1 - lμvi)^2
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
    # add to global likelihood
    irλ = llbdλ*δt
    irμ = llbdμ*δt
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

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi  = lλv[nI+1]
      lμvi  = lμv[nI+1]
      srfdt = sqrt(fdt)
      ll   += ldnorm_bm(lλvi1, lλvi + α*fdt, srfdt*σλ) +
             ldnorm_bm(lμvi1, lμvi, srfdt*σμ)
      ssλ  += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
      ssμ  += (lμvi1 - lμvi)^2/(2.0*fdt)
      nλ   += 1.0
      iriλ  = fdt*(exp(0.5*(lλvi + lλvi1)))
      iriμ  = fdt*(exp(0.5*(lμvi + lμvi1)))
      ll   -= (iriλ + iriμ)
      irλ  += iriλ 
      irμ  += iriμ
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

  return ll, ix, ddλ, ssλ, ssμ, nλ, irλ, irμ
end



