#=

clads birth-death MH proposals for internal updates

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#





"""
    _stem_update!(ξi      ::cTbd,
                  eds     ::Float64,
                  λ1      ::Float64,
                  λ2      ::Float64,
                  μ1      ::Float64,
                  μ2      ::Float64,
                  α       ::Float64,
                  σλ      ::Float64,
                  σμ      ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ddλ     ::Float64,
                  ssλ     ::Float64,
                  ssμ     ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  λ0_prior::NTuple{2,Float64},
                  μ0_prior::NTuple{2,Float64},
                  surv    ::Int64)

Do `clads` update for stem root.
"""
function _stem_update!(ξi      ::cTbd,
                       eds     ::Float64,
                       λ1      ::Float64,
                       λ2      ::Float64,
                       μ1      ::Float64,
                       μ2      ::Float64,
                       α       ::Float64,
                       σλ      ::Float64,
                       σμ      ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ddλ     ::Float64,
                       ssλ     ::Float64,
                       ssμ     ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       λ0_prior::NTuple{2,Float64},
                       μ0_prior::NTuple{2,Float64},
                       surv    ::Int64)

  @inbounds begin
    λi, μi = lλ(ξi), lμ(ξi)
    ei = e(ξi)

    if def1(ξi)
      eds, λ1, λ2, μ1, μ2 = 0.0, lλ(ξi.d1), lλ(ξi.d2), lμ(ξi.d1), lμ(ξi.d2)
    end

    ## node proposal
    # speciation
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], σλ^2, σλ^2, λ0_prior[2])
   # extinction
    μr = trioprop(μ1,     μ2,     μ0_prior[1], σμ^2, σμ^2, μ0_prior[2])

    llrbd = λr - λi + (ei + eds)*(exp(λi) - exp(λr) + exp(μi) - exp(μr))

    lU = -randexp()

    if lU < llrbd + log(1000.0/mc)

      mp     = m_surv_cladsbd(th, λr, μr, α, σλ, σμ, 1_000, surv)
      llrbd += log(mp/mc)

      if lU < llrbd
        llc += llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ) + 
               llrdnorm2_μ(μ1, μ2,     μr,     μi, σμ) + llrbd
        prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
               llrdnorm_x(μr, μi, μ0_prior[1], μ0_prior[2])
        ddλ += 2.0*(λi - λr)
        ssλ += 0.5*(
                (λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
                (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
        ssμ += 0.5*((μ1 - μr)^2 + (μ2 - μr)^2 - (μ1 - μi)^2 - (μ2 - μi)^2)
        mc  = mp
        λi, μi  = λr, μr
        setlλ!(ξi, λi)
        setlμ!(ξi, μi)
      end
    end
  end

  return llc, prc, ddλ, ssλ, ssμ, mc, λi, μi
end




"""
    _crown_update!(ξi      ::cTbd,
                   ξ1      ::cTbd,
                   ξ2      ::cTbd,
                   α       ::Float64,
                   σλ      ::Float64,
                   σμ      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   ddλ     ::Float64,
                   ssλ     ::Float64,
                   ssμ     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   λ0_prior::NTuple{2,Float64},
                   μ0_prior::NTuple{2,Float64},
                   surv    ::Int64)

Do `clads` update for crown root.
"""
function _crown_update!(ξi      ::cTbd,
                        ξ1      ::cTbd,
                        ξ2      ::cTbd,
                        α       ::Float64,
                        σλ      ::Float64,
                        σμ      ::Float64,
                        llc     ::Float64,
                        prc     ::Float64,
                        ddλ     ::Float64,
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
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], σλ^2, σλ^2, λ0_prior[2])
    # extinction
    μr = trioprop(μ1,     μ2,     μ0_prior[1], σμ^2, σμ^2, μ0_prior[2])

    # survival ratio
    mp  = m_surv_cladsbd(th, λr, μr, α, σλ, σμ, 1_000, surv)
    llr = log(mp/mc)

    if -randexp() < llr
      llc += llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ) +
             llrdnorm2_μ(μ1, μ2,     μr,     μi, σμ) + llr
      prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2]) + 
             llrdnorm_x(μr, μi, μ0_prior[1], μ0_prior[2])
      ddλ += 2.0*(λi - λr)
      ssλ += 0.5*((λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
                  (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      ssμ += 0.5*((μ1 - μr)^2 + (μ2 - μr)^2 - (μ1 - μi)^2 - (μ2 - μi)^2)
      mc  = mp
      setlλ!(ξi, λr)
      setlμ!(ξi, μr)
    end
  end

  return llc, prc, ddλ, ssλ, ssμ, mc
end




"""
    _update_internal!(tree::T,
                      bi  ::iBffs,
                      eas ::Float64,
                      λa  ::Float64,
                      μa  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      eds ::Float64,
                      λ1  ::Float64,
                      λ2  ::Float64,
                      μ1  ::Float64,
                      μ2  ::Float64,
                      llc ::Float64,
                      ddλ ::Float64,
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
                           α   ::Float64,
                           σλ  ::Float64,
                           σμ  ::Float64,
                           eds ::Float64,
                           λ1  ::Float64,
                           λ2  ::Float64,
                           μ1  ::Float64,
                           μ2  ::Float64,
                           llc ::Float64,
                           ddλ ::Float64,
                           ssλ ::Float64,
                           ssμ ::Float64,
                           ter ::Bool) where {T <: cT}

  if def1(tree)
    if isfinite(λa)
      llc, ddλ, ssλ, ssμ, λa, μa = 
        update_triad!(tree, eas, λa, μa, α, σλ, σμ, llc, ddλ, ssλ, ssμ)
    else
      λa, μa = lλ(tree), lμ(tree)
    end

    llc, ddλ, ssλ, ssμ, λx, μx =
      _update_internal!(tree.d1, bi, 0.0, λa, μa, α, σλ, σμ, eds, λ1, λ2, 
        μ1, μ2, llc, ddλ, ssλ, ssμ, ter)
    llc, ddλ, ssλ, ssμ, λx, μx =
      _update_internal!(tree.d2, bi, 0.0, λa, μa, α, σλ, σμ, eds, λ1, λ2, 
        μ1, μ2, llc, ddλ, ssλ, ssμ, ter)
  else 
    if isfix(tree)
      # if leads to eventual speciation
      if isfinite(λ1)
        llc, ddλ, ssλ, ssμ = 
          update_faketip!(tree, bi, eas, λa, μa, eds, λ1, λ2, μ1, μ2, α, σλ, σμ,
            llc, ddλ, ssλ, ssμ)
      # if leads to non-speciation
      else
      llc, ddλ, ssλ, ssμ = 
        update_tip!(tree, eas, λa, μa, eds, α, σλ, σμ, llc, ddλ, ssλ, ssμ)
      end
    # if DA tip
    else
      llc, ddλ, ssλ, ssμ = 
        update_tip!(tree, eas, λa, μa, 0.0, α, σλ, σμ, llc, ddλ, ssλ, ssμ)
    end
  end

  return llc, ddλ, ssλ, ssμ, λa, μa
end




"""
    update_triad!(tree::T,
                  eas ::Float64,
                  λa  ::Float64,
                  μa  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64,
                  ssμ ::Float64) where {T <: cT}

Make a trio proposal for clads.
"""
function update_triad!(tree::T,
                       eas ::Float64,
                       λa  ::Float64,
                       μa  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       ssμ ::Float64) where {T <: cT}

  @inbounds begin

    ei = e(tree)
    λi, λ1, λ2 = lλ(tree), lλ(tree.d1), lλ(tree.d2)
    μi, μ1, μ2 = lμ(tree), lμ(tree.d1), lμ(tree.d2)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)
    μn = trioprop(μa,     μ1,     μ2,     σμ)

    # likelihood ratio
    llrbd = λn - λi + (ei + eas)*(exp(λi) - exp(λn) + exp(μi) - exp(μn))

    if -randexp() < llrbd
      llc += llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ) + 
             llrdnorm3(μa,     μ1,     μ2,     μn, μi, σμ) + llrbd
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      ssμ += 0.5*(
              (μn - μa)^2 + (μ1 - μn)^2 + (μ2 - μn)^2 -
              (μi - μa)^2 - (μ1 - μi)^2 - (μ2 - μi)^2)
      λi   = λn
      μi   = μn
      setlλ!(tree, λi)
      setlμ!(tree, μi)
    end
  end

  return llc, ddλ, ssλ, ssμ, λi, μi
end





"""
    update_tip!(tree::cTbd,
                eas ::Float64,
                λa  ::Float64,
                eds ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                ϵ   ::Float64,
                llc ::Float64,
                ddλ ::Float64,
                ssλ ::Float64,
                seλ ::Float64)

Make a `clads` tip proposal.
"""
function update_tip!(tree::cTbd,
                     eas ::Float64,
                     λa  ::Float64,
                     μa  ::Float64,
                     eds ::Float64,
                     α   ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     llc ::Float64,
                     ddλ ::Float64,
                     ssλ ::Float64,
                     ssμ ::Float64)

  @inbounds begin

    ei = e(tree)
    λi = lλ(tree)
    μi = lμ(tree)

    # node proposal
    λn = rnorm(λa + α, σλ)
    μn = rnorm(μa,     σμ)

    # likelihood ratios
    llrbd = (eas + ei + eds) * (exp(λi) - exp(λn) + exp(μi) - exp(μn))

    if isextinct(tree)
      llrbd += μn - μi
    end

    if -randexp() < llrbd
      llc += llrdnorm_x(λn, λi, λa + α, σλ^2) + 
             llrdnorm_x(μn, μi, μa,     σμ^2) + llrbd
      ddλ += λn - λi
      ssλ += 0.5*((λn - λa - α)^2 - (λi - λa - α)^2)
      ssμ += 0.5*((μn - μa)^2 - (μi - μa)^2)
      setlλ!(tree, λn)
      setlμ!(tree, μn)
    end
  end

  return llc, ddλ, ssλ, ssμ
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
                    α   ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    llc ::Float64,
                    ddλ ::Float64,
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
                         α   ::Float64,
                         σλ  ::Float64,
                         σμ  ::Float64,
                         llc ::Float64,
                         ddλ ::Float64,
                         ssλ ::Float64,
                         ssμ ::Float64) where {T <: cT}
  @inbounds begin

    ei = e(tree)
    λi = lλ(tree)
    μi = lμ(tree)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)
    μn = trioprop(μa,     μ1,     μ2,     σμ)

    # likelihood ratios
    llrbd = λn - λi + (eas + ei + eds)*(exp(λi) - exp(λn) + exp(μi) - exp(μn))

    if -randexp() < llrbd
      llc += llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ) + 
             llrdnorm3(μa,     μ1,     μ2,     μn, μi, σμ) + llrbd
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      ssμ += 0.5*(
              (μn - μa)^2 + (μ1 - μn)^2 + (μ2 - μn)^2 -
              (μi - μa)^2 - (μ1 - μi)^2 - (μ2 - μi)^2)
      λi   = λn
      μi   = μn
      setlλ!(tree, λi)
      setlμ!(tree, μi)
      setλt!(bi, λi)
      setμt!(bi, μi)
    end
  end

  return llc, ddλ, ssλ, ssμ 
end



