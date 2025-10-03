#=

clads fossilized birth-death MH proposals for internal updates

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    _stem_update!(ξi      ::cTfbd,
                  eds     ::Float64,
                  λ1      ::Float64,
                  λ2      ::Float64,
                  μ1      ::Float64,
                  μ2      ::Float64,
                  αλ      ::Float64,
                  αμ      ::Float64,
                  σλ      ::Float64,
                  σμ      ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ddλ     ::Float64,
                  ddμ     ::Float64,
                  ssλ     ::Float64,
                  ssμ     ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  λ0_prior::NTuple{2,Float64},
                  μ0_prior::NTuple{2,Float64},
                  surv    ::Int64)

Do `clads` update for stem root.
"""
function _stem_update!(ξi      ::cTfbd,
                       eds     ::Float64,
                       λ1      ::Float64,
                       λ2      ::Float64,
                       μ1      ::Float64,
                       μ2      ::Float64,
                       αλ      ::Float64,
                       αμ      ::Float64,
                       σλ      ::Float64,
                       σμ      ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ddλ     ::Float64,
                       ddμ     ::Float64,
                       ssλ     ::Float64,
                       ssμ     ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       λ0_prior::NTuple{2,Float64},
                       μ0_prior::NTuple{2,Float64},
                       surv    ::Int64)

  @inbounds begin
    λi, μi = lλ(ξi), lλ(ξi)
    ei = e(ξi)

    ## node proposal
    # speciation
    λr = trioprop(λ1 - αλ, λ2 - αλ, λ0_prior[1], σλ^2, σλ^2, λ0_prior[2])
    # extinction
    μr = trioprop(μ1 - αμ, μ2 - αμ, μ0_prior[1], σμ^2, σμ^2, μ0_prior[2])

    llrbm = llrdnorm2_μ(λ1, λ2, λr + αλ, λi + αλ, σλ) + 
            llrdnorm2_μ(μ1, μ2, μr + αμ, μi + αμ, σμ)
    llrbd = λr - λi + (ei + eds)*(exp(λi) - exp(λr) + exp(μi) - exp(μr))

    if lU < llrbd + log(1000.0/mc)

      mp     = m_surv_cladsfbd(th, λr, μr, αλ, αμ, σλ, σμ, 1_000, surv)
      llrbd += log(mp/mc)

      if -randexp() < llrbd
        llc += llrbm + llrbd
        prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
               llrdnorm_x(μr, μi, μ0_prior[1], μ0_prior[2])
        ddλ += 2.0*(λi - λr)
        ddμ += 2.0*(μi - μr)
        ssλ += 0.5*(
                (λ1 - λr - αλ)^2 + (λ2 - λr - αλ)^2 - 
                (λ1 - λi - αλ)^2 - (λ2 - λi - αλ)^2)
        ssμ += 0.5*(
                (μ1 - μr - αμ)^2 + (μ2 - μr - αμ)^2 - 
                (μ1 - μi - αμ)^2 - (μ2 - μi - αμ)^2)
        mc  = mp
        λi, μi  = λr, μr
        setlλ!(ξi, λi)
        setlμ!(ξi, μi)
      end
    end
  end

  return llc, prc, ddλ, ddμ, ssλ, ssμ, mc, λi, μi
end




"""
    _crown_update!(ξi      ::cTfbd,
                   ξ1      ::cTfbd,
                   ξ2      ::cTfbd,
                   αλ      ::Float64,
                   αμ      ::Float64,
                   σλ      ::Float64,
                   σμ      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   ddλ     ::Float64,
                   ddμ     ::Float64,
                   ssλ     ::Float64,
                   ssμ     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   λ0_prior::NTuple{2,Float64},
                   μ0_prior::NTuple{2,Float64},
                   surv    ::Int64)

Do `clads` update for crown root.
"""
function _crown_update!(ξi      ::cTfbd,
                        ξ1      ::cTfbd,
                        ξ2      ::cTfbd,
                        αλ      ::Float64,
                        αμ      ::Float64,
                        σλ      ::Float64,
                        σμ      ::Float64,
                        llc     ::Float64,
                        prc     ::Float64,
                        ddλ     ::Float64,
                        ddμ     ::Float64,
                        ssλ     ::Float64,
                        ssμ     ::Float64,
                        mc      ::Float64,
                        th      ::Float64,
                        λ0_prior::NTuple{2,Float64},
                        μ0_prior::NTuple{2,Float64},
                        surv    ::Int64)

  @inbounds begin
    λi, λ1, λ2 = lλ(ξi), lλ(ξ1), lλ(ξ2)
    μi, μ1, μ2 = lμ(ξi), lμ(ξ1), lμ(ξ2)

    ## node proposal
    # speciation
    λr = trioprop(λ1 - αλ, λ2 - αλ, λ0_prior[1], σλ^2, σλ^2, λ0_prior[2])
    # extinction
    μr = trioprop(μ1 - αμ, μ2 - αμ, μ0_prior[1], σμ^2, σμ^2, μ0_prior[2])

    # survival ratio
    mp  = m_surv_cladsfbd(th, λr, μr, αλ, αμ, σλ, σμ, 1_000, surv)
    llr = log(mp/mc)

    if -randexp() < llr
      llc += llrdnorm2_μ(λ1, λ2, λr + αλ, λi + αλ, σλ) +
             llrdnorm2_μ(μ1, μ2, μr + αμ, μi + αμ, σμ) + llr
      prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2]) + 
             llrdnorm_x(μr, μi, μ0_prior[1], μ0_prior[2])
      ddλ += 2.0*(λi - λr)
      ddμ += 2.0*(μi - μr)
      ssλ += 0.5*((λ1 - λr - αλ)^2 + (λ2 - λr - αλ)^2 - 
                  (λ1 - λi - αλ)^2 - (λ2 - λi - αλ)^2)
      ssμ += 0.5*((μ1 - μr - αμ)^2 + (μ2 - μr - αμ)^2 - 
                  (μ1 - μi - αμ)^2 - (μ2 - μi - αμ)^2)
      mc  = mp
      setlλ!(ξi, λr)
      setlμ!(ξi, μr)
    end
  end

  return llc, prc, ddλ, ddμ, ssλ, ssμ, mc
end




"""
    _update_internal!(tree::T,
                      bi  ::iBffs,
                      eas ::Float64,
                      λa  ::Float64,
                      μa  ::Float64,
                      αλ  ::Float64,
                      αμ  ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      eds ::Float64,
                      λ1  ::Float64,
                      λ2  ::Float64,
                      μ1  ::Float64,
                      μ2  ::Float64,
                      llc ::Float64,
                      ddλ ::Float64,
                      ddμ ::Float64,
                      ssλ ::Float64,
                      ssμ ::Float64,
                      ter ::Bool) where {T <: cT}

Do `clads` internal rate updates on a decoupled tree recursively.
"""
function _update_internal!(tree::T,
                           bi  ::iBffs,
                           eas ::Float64,
                           λa  ::Float64,
                           μa  ::Float64,
                           αλ  ::Float64,
                           αμ  ::Float64,
                           σλ  ::Float64,
                           σμ  ::Float64,
                           eds ::Float64,
                           λ1  ::Float64,
                           λ2  ::Float64,
                           μ1  ::Float64,
                           μ2  ::Float64,
                           llc ::Float64,
                           ddλ ::Float64,
                           ddμ ::Float64,
                           ssλ ::Float64,
                           ssμ ::Float64,
                           ter ::Bool) where {T <: cT}

  if def1(tree)
    if def2(tree)
      llc, ddλ, ddμ, ssλ, ssμ, λa, μa = 
        update_triad!(tree, eas, λa, μa, αλ, αμ, σλ, σμ, 
          llc, ddλ, ddμ, ssλ, ssμ)

      llc, ddλ, ddμ, ssλ, ssμ, λx, μx =
        _update_internal!(tree.d1, bi, 0.0, λa, μa, αλ, αμ, σλ, σμ, eds, λ1, λ2, 
          μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ, ter)
      llc, ddλ, ddμ, ssλ, ssμ, λx, μx =
        _update_internal!(tree.d2, bi, 0.0, λa, μa, αλ, αμ, σλ, σμ, eds, λ1, λ2, 
          μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ, ter)
    else
      llc, ddλ, ddμ, ssλ, ssμ, λx, μx =
        _update_internal!(tree.d1, bi, eas + e(tree), λa, μa, αλ, αμ, σλ, σμ, 
          eds, λ1, λ2, μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ, ter)
    end
  else 
    # if leads to non-speciation
    if !isfix(tree) || ter || isnan(λ1)
      llc, ddλ, ddμ, ssλ, ssμ = 
        update_tip!(tree, eas, λa, μa, eds, αλ, αμ, σλ, σμ, 
          llc, ddλ, ddμ, ssλ, ssμ)
    # if leads to eventual speciation
    else
      llc, ddλ, ddμ, ssλ, ssμ = 
        update_faketip!(tree, bi, eas, λa, μa, eds, λ1, λ2, μ1, μ2, 
          αλ, αμ, σλ, σμ, llc, ddλ, ddμ, ssλ, ssμ)
    end
  end

  return llc, ddλ, ddμ, ssλ, ssμ, λa, μa
end




"""
    update_triad!(tree::T,
                  eas ::Float64,
                  λa  ::Float64,
                  μa  ::Float64,
                  αλ  ::Float64,
                  αμ  ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ddμ ::Float64,
                  ssλ ::Float64,
                  ssμ ::Float64) where {T <: cT}

Make a trio proposal for clads.
"""
function update_triad!(tree::T,
                       eas ::Float64,
                       λa  ::Float64,
                       μa  ::Float64,
                       αλ  ::Float64,
                       αμ  ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ddμ ::Float64,
                       ssλ ::Float64,
                       ssμ ::Float64) where {T <: cT}

  @inbounds begin

    ei = e(tree)
    λi, λ1, λ2 = lλ(tree), lλ(tree.d1), lλ(tree.d2)
    μi, μ1, μ2 = lμ(tree), lμ(tree.d1), lμ(tree.d2)

    # node proposal
    λn = trioprop(λa + αλ, λ1 - αλ, λ2 - αλ, σλ)
    μn = trioprop(μa + αμ, μ1 - αμ, μ2 - αμ, σμ)

    # likelihood ratios
    llrbm = llrdnorm3(λa + αλ, λ1 - αλ, λ2 - αλ, λn, λi, σλ) + 
            llrdnorm3(μa + αμ, μ1 - αμ, μ2 - αμ, μn, μi, σμ)
    llrbd = λn - λi + (ei + eas)*(exp(λi) - exp(λn) + exp(μi) - exp(μn))

    if -randexp() < llrbd
      llc += llrbm + llrbd
      ddλ += (λi - λn)
      ddμ += (μi - μn)
      ssλ += 0.5*(
              (λn - λa - αλ)^2 + (λ1 - λn - αλ)^2 + (λ2 - λn - αλ)^2 -
              (λi - λa - αλ)^2 - (λ1 - λi - αλ)^2 - (λ2 - λi - αλ)^2)
      ssμ += 0.5*(
              (μn - μa - αμ)^2 + (μ1 - μn - αμ)^2 + (μ2 - μn - αμ)^2 -
              (μi - μa - αμ)^2 - (μ1 - μi - αμ)^2 - (μ2 - μi - αμ)^2)
      λi   = λn
      μi   = μn
      setlλ!(tree, λi)
      setlμ!(tree, μi)
    end
  end

  return llc, ddλ, ddμ, ssλ, ssμ, λi, μi
end




"""
    update_tip!(tree::cTfbd,
                eas ::Float64,
                λa  ::Float64,
                μa  ::Float64,
                eds ::Float64,
                αλ  ::Float64,
                αμ  ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                llc ::Float64,
                ddλ ::Float64,
                ddμ ::Float64,
                ssλ ::Float64,
                ssμ ::Float64)

Make a `clads` tip proposal.
"""
function update_tip!(tree::cTfbd,
                     eas ::Float64,
                     λa  ::Float64,
                     μa  ::Float64,
                     eds ::Float64,
                     αλ  ::Float64,
                     αμ  ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     llc ::Float64,
                     ddλ ::Float64,
                     ddμ ::Float64,
                     ssλ ::Float64,
                     ssμ ::Float64)

  @inbounds begin

    ei = e(tree)
    λi, μi = lλ(tree), lμ(tree)

    # node proposal
    λn = rnorm(λa + αλ, σλ)
    μn = rnorm(μa + αμ, σμ)

    # likelihood ratios
    llrbm = llrdnorm_x(λn, λi, λa + αλ, σλ^2) + 
            llrdnorm_x(μn, μi, μa + αμ, σμ^2)
    llrbd = (eas + ei + eds) * (exp(λi) - exp(λn) + exp(μi) - exp(μn))

    if isextinct(tree)
      llrbd += μn - μi
    end

    if -randexp() < llrbd
      llc += llrbm + llrbd
      ddλ += λn - λi
      ddμ += μn - μi
      ssλ += 0.5*((λn - λa - αλ)^2 - (λi - λa - αλ)^2)
      ssμ += 0.5*((μn - μa - αμ)^2 - (μi - μa - αμ)^2)
      setlλ!(tree, λn)
      setlμ!(tree, μn)
    end
  end

  return llc, ddλ, ddμ, ssλ, ssμ
end




"""
    update_faketip!(tree::T,
                    bi  ::iBffs,
                    eas ::Float64,
                    λa  ::Float64,
                    μa  ::Float64,
                    eds ::Float64,
                    λ1  ::Float64,
                    λ2  ::Float64,
                    μ1  ::Float64,
                    μ2  ::Float64,
                    αλ  ::Float64,
                    αμ  ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    llc ::Float64,
                    ddλ ::Float64,
                    ddμ ::Float64,
                    ssλ ::Float64,
                    ssμ ::Float64) where {T <: cT}

Make a `clads` tip proposal.
"""
function update_faketip!(tree::T,
                         bi  ::iBffs,
                         eas ::Float64,
                         λa  ::Float64,
                         μa  ::Float64,
                         eds ::Float64,
                         λ1  ::Float64,
                         λ2  ::Float64,
                         μ1  ::Float64,
                         μ2  ::Float64,
                         αλ  ::Float64,
                         αμ  ::Float64,
                         σλ  ::Float64,
                         σμ  ::Float64,
                         llc ::Float64,
                         ddλ ::Float64,
                         ddμ ::Float64,
                         ssλ ::Float64,
                         ssμ ::Float64) where {T <: cT}
  @inbounds begin

    ei = e(tree)
    λi, μi = lλ(tree), lμ(tree)

    # node proposal
    λn = trioprop(λa + αλ, λ1 - αλ, λ2 - αλ, σλ)
    μn = trioprop(μa + αμ, μ1 - αμ, μ2 - αμ, σμ)

    # likelihood ratios
    llrbm = llrdnorm3(λa + αλ, λ1 - αλ, λ2 - αλ, λn, λi, σλ) + 
            llrdnorm3(μa + αμ, μ1 - αμ, μ2 - αμ, μn, μi, σμ)
    llrbd = λn - λi + (eas + ei + eds)*(exp(λi) - exp(λn) + exp(μi) - exp(μn))

    if -randexp() < llrbd
      llc += llrbm + llrbd
      ddλ += (λi - λn)
      ddμ += (μi - μn)
      ssλ += 0.5*(
              (λn - λa - αλ)^2 + (λ1 - λn - αλ)^2 + (λ2 - λn - αλ)^2 -
              (λi - λa - αλ)^2 - (λ1 - λi - αλ)^2 - (λ2 - λi - αλ)^2)
      ssμ += 0.5*(
              (μn - μa - αμ)^2 + (μ1 - μn - αμ)^2 + (μ2 - μn - αμ)^2 -
              (μi - μa - αμ)^2 - (μ1 - μi - αμ)^2 - (μ2 - μi - αμ)^2)
      λi   = λn
      μi   = μn
      setlλ!(tree, λi)
      setlμ!(tree, μi)
      setλt!(bi, λi)
      setμt!(bi, μi)
    end
  end

  return llc, ddλ, ddμ, ssλ, ssμ 
end



