#=

`gbmfbd` likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 10 03 2022
=#




"""
    llik_gbm(tree::iTfbd,
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTfbd` according to `gbmfbd`.
"""
function llik_gbm(tree::iTfbd,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  ψ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  if def1(tree)
    if def2(tree)
      ll_gbm_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, δt, fdt(tree), srδt,
        true, false, false)                                         +
      llik_gbm(tree.d1, α, σλ, σμ, ψ, δt, srδt)                     +
      llik_gbm(tree.d2, α, σλ, σμ, ψ, δt, srδt)
    else
      ll_gbm_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, δt, fdt(tree), srδt,
        false, false, true)                                         +
      llik_gbm(tree.d1, α, σλ, σμ, ψ, δt, srδt)
    end
  else
    ll_gbm_b(lλ(tree), lμ(tree), α, σλ, σμ, ψ, δt, fdt(tree), srδt,
        false, isextinct(tree), false)
  end
end




"""
    llik_gbm(Ξ   ::Vector{iTfbd},
             idf ::Vector{iBffs},
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             ψ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTfbd` according to `gbm-bd`.
"""
function llik_gbm(Ξ   ::Vector{iTfbd},
                  idf ::Vector{iBffs},
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  ψ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      bi  = idf[i]
      ll += llik_gbm(Ξ[i], α, σλ, σμ, ψ, δt, srδt)
      if !it(bi) && !isfossil(bi)
        ll += λt(bi)
      end
    end
  end

  return ll
end




"""
    ll_gbm_b(lλv ::Array{Float64,1},
             lμv ::Array{Float64,1},
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             ψ   ::Float64,
             δt  ::Float64,
             fdt ::Float64,
             srδt::Float64,
             λev ::Bool,
             μev ::Bool,
             ψev ::Bool)

Returns the log-likelihood for a branch according to `gbmfbd`.
"""
function ll_gbm_b(lλv ::Array{Float64,1},
                  lμv ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  ψ   ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool,
                  μev ::Bool,
                  ψev ::Bool)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ  = 0.0
    llμ  = 0.0
    llbd = 0.0
    @avx for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi - α*δt)^2
      llμ  += (lμvi1 - lμvi)^2
      llbd += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1))
    end

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) +
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    # add to global likelihood
    ll -= llbd*δt

    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi = lλv[nI+1]
      lμvi = lμv[nI+1]
      srfdt = sqrt(fdt)
      ll += ldnorm_bm(lλvi1, lλvi + α*fdt, srfdt*σλ)
      ll += ldnorm_bm(lμvi1, lμvi, srfdt*σμ)
      ll -= fdt*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)))
    end

    #if speciation
    if λev
      ll += lλvi1
    # if extinction
    elseif μev
      ll += lμvi1
    # if fossilization
    elseif ψev
      ll += log(ψ)
    end
  end

  return ll
end




"""
    llik_gbm_ss(tree::iTfbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

Returns the log-likelihood for a `iTfbd` according to `gbmfbd`.
"""
function llik_gbm_ss(tree::iTfbd,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     ψ   ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  if def1(tree)
    ll, dλ, ssλ, ssμ, nλ =
      ll_gbm_b_ss(lλ(tree), lμ(tree), α, σλ, σμ, ψ, δt, fdt(tree), srδt,
        false, false, true)
    ll1, dλ1, ssλ1, ssμ1, nλ1 =
      llik_gbm_ss(tree.d1, α, σλ, σμ, ψ, δt, srδt)
    ll  += ll1
    dλ  += dλ1
    ssλ += ssλ1
    ssμ += ssμ1
    nλ  += nλ1

    if def2(tree)
      ll2, dλ2, ssλ2, ssμ2, nλ2 =
        llik_gbm_ss(tree.d2, α, σλ, σμ, ψ, δt, srδt)
      ll  += ll2
      dλ  += dλ2
      ssλ += ssλ2
      ssμ += ssμ2
      nλ  += nλ2
    end
  else
    ll, dλ, ssλ, ssμ, nλ =
      ll_gbm_b_ss(lλ(tree), lμ(tree), α, σλ, σμ, ψ, δt, fdt(tree), srδt,
          false, isextinct(tree), false)
  end

  return ll, dλ, ssλ, ssμ, nλ
end




"""
    ll_gbm_b_ss(lλv ::Array{Float64,1},
                lμv ::Array{Float64,1},
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64,
                fdt ::Float64,
                srδt::Float64,
                λev ::Bool,
                μev ::Bool)

Returns the log-likelihood for a branch according to `gbmfbd`.
"""
function ll_gbm_b_ss(lλv ::Array{Float64,1},
                     lμv ::Array{Float64,1},
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     ψ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64,
                     λev ::Bool,
                     μev ::Bool,
                     ψev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2

    llλ  = 0.0
    llμ  = 0.0
    llbd = 0.0
    @avx for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lμvi  = lμv[i]
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi - α*δt)^2
      llμ  += (lμvi1 - lμvi)^2
      llbd += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1))
    end

    # standardized sum of squares
    ssλ = llλ/(2.0*δt)
    ssμ = llμ/(2.0*δt)
    nλ  = Float64(nI)

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) +
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    # add to global likelihood
    ll -= llbd*δt

    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]

    dλ = lλvi1 - lλv[1]

    # add final non-standard `δt`
    if fdt > 0.0
      lλvi  = lλv[nI+1]
      lμvi  = lμv[nI+1]
      srfdt = sqrt(fdt)
      ll  += ldnorm_bm(lλvi1, lλvi + α*fdt, srfdt*σλ)
      ll  += ldnorm_bm(lμvi1, lμvi, srfdt*σμ)
      ll  -= fdt*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)))
      ssλ += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
      ssμ += (lμvi1 - lμvi)^2/(2.0*fdt)
      nλ  += 1.0
    end

      #if speciation
    if λev
      ll += lλvi1
    # if extinction
    elseif μev
      ll += lμvi1
    # if fossilization
    elseif ψev
      ll += log(ψ)
    end
  end

  return ll, dλ, ssλ, ssμ, nλ
end




"""
    _deltaλ(tree::iTfbd)

Returns the log-likelihood ratio for a `iTfbd` according
to GBM birth-death for a `α` proposal.
"""
function _deltaλ(tree::iTfbd)

  lλv = lλ(tree)

  if def1(tree)
    lλv[end] - lλv[1] + _deltaλ(tree.d1) + 
    (def2(tree) ? _deltaλ(tree.d2) : 0.0)
  else
    lλv[end] - lλv[1]
  end
end








