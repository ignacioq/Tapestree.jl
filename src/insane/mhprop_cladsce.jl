#=

clads constant extinction MH proposals for internal updates

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    _stem_update!(ξi      ::cTce,
                  eds     ::Float64,
                  λ1      ::Float64,
                  λ2      ::Float64,
                  α       ::Float64,
                  σλ      ::Float64,
                  μ       ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ddλ     ::Float64,
                  ssλ     ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  λ0_prior::NTuple{2,Float64},
                  surv    ::Int64)

Do `clads` update for crown root.
"""
function _stem_update!(ξi      ::cTce,
                       eds     ::Float64,
                       λ1      ::Float64,
                       λ2      ::Float64,
                       α       ::Float64,
                       σλ      ::Float64,
                       μ       ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ddλ     ::Float64,
                       ssλ     ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       λ0_prior::NTuple{2,Float64},
                       surv    ::Int64)

  @inbounds begin
    λi = lλ(ξi)
    ei = e(ξi)

    if def1(ξi)
      eds, λ1, λ2 = 0.0, lλ(ξi.d1), lλ(ξi.d2)
    end

    # node proposal
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], 
                  σλ^2,     σλ^2, λ0_prior[2])

    llrce = λr - λi + (ei + eds)*(exp(λi) - exp(λr))

    lU = -randexp()

    if lU < llrce + log(1000.0/mc)

      mp     = m_surv_cladsce(th, λr, α, σλ, μ, 1_000, surv)
      llrce += log(mp/mc)

      if lU < llrce
        llc += llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ) + llrce
        prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
        ddλ += 2.0*(λi - λr)
        ssλ += 0.5*(
                (λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
                (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
        mc  = mp
        λi  = λr
        setlλ!(ξi, λi)
      end
    end
  end

  return llc, prc, ddλ, ssλ, mc, λi
end




"""
    _crown_update!(ξi      ::cTce,
                   ξ1      ::cTce,
                   ξ2      ::cTce,
                   α       ::Float64,
                   σλ      ::Float64,
                   μ       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   ddλ     ::Float64,
                   ssλ     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   λ0_prior::NTuple{2,Float64},
                   surv    ::Int64)

Do `clads` update for crown root.
"""
function _crown_update!(ξi      ::cTce,
                        ξ1      ::cTce,
                        ξ2      ::cTce,
                        α       ::Float64,
                        σλ      ::Float64,
                        μ       ::Float64,
                        llc     ::Float64,
                        prc     ::Float64,
                        ddλ     ::Float64,
                        ssλ     ::Float64,
                        mc      ::Float64,
                        th      ::Float64,
                        λ0_prior::NTuple{2,Float64},
                        surv    ::Int64)

  @inbounds begin
    λi = lλ(ξi)
    λ1 = lλ(ξ1)
    λ2 = lλ(ξ2)

    # node proposal
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], 
                  σλ^2,   σλ^2,   λ0_prior[2])

    # survival ratio
    mp  = m_surv_cladsce(th, λr, α, σλ, μ, 1_000, surv)
    llr = log(mp/mc)

    if -randexp() < llr
      llc += llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ) + llr
      prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
      ddλ += 2.0*(λi - λr)
      ssλ += 0.5*((λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
                  (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      mc  = mp
      setlλ!(ξi, λr)
    end
  end

  return llc, prc, ddλ, ssλ, mc
end




"""
    _update_internal!(tree::T,
                      bi  ::iBffs,
                      eas ::Float64,
                      λa  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      eds ::Float64,
                      λ1  ::Float64,
                      λ2  ::Float64,
                      llc ::Float64,
                      ddλ ::Float64,
                      ssλ ::Float64,
                      ter ::Bool) where {T <: cT}

Do `clads` internal rate updates on a decoupled tree recursively.
"""
function _update_internal!(tree::T,
                           bi  ::iBffs,
                           eas ::Float64,
                           λa  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           eds ::Float64,
                           λ1  ::Float64,
                           λ2  ::Float64,
                           llc ::Float64,
                           ddλ ::Float64,
                           ssλ ::Float64,
                           ter ::Bool) where {T <: cT}

  if def1(tree)
    if isfinite(λa)
      llc, ddλ, ssλ, λa = 
        update_triad!(tree, eas, λa, α, σλ, llc, ddλ, ssλ)
    else
      λa = lλ(tree)
    end

    llc, ddλ, ssλ, λx =
      _update_internal!(tree.d1, bi, 0.0, λa, α, σλ, eds, λ1, λ2, 
        llc, ddλ, ssλ, ter)
    llc, ddλ, ssλ, λx =
      _update_internal!(tree.d2, bi, 0.0, λa, α, σλ, eds, λ1, λ2, 
        llc, ddλ, ssλ, ter)
  else 
    if isfix(tree)
      # if leads to eventual speciation
      if isfinite(λ1)
        llc, ddλ, ssλ = 
          update_faketip!(tree, bi, eas, λa, eds, λ1, λ2, α, σλ, llc, ddλ, ssλ)
      # if leads to non-speciation
      else
        llc, ddλ, ssλ = 
          update_tip!(tree, eas, λa, eds, α, σλ, llc, ddλ, ssλ)
      end
    # if DA tip
    else
      llc, ddλ, ssλ = 
        update_tip!(tree, eas, λa, 0.0, α, σλ, llc, ddλ, ssλ)
    end
  end

  return llc, ddλ, ssλ, λa
end




"""
    update_triad!(tree::T,
                  eas ::Float64,
                  λa  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64) where {T <: cT}

Make a trio proposal for clads.
"""
function update_triad!(tree::T,
                       eas ::Float64,
                       λa  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64) where {T <: cT}

  @inbounds begin

    λi = lλ(tree)
    λ1 = lλ(tree.d1)
    λ2 = lλ(tree.d2)
    ei = e(tree)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)

    # likelihood ratios
    llrce = λn - λi + (ei + eas)*(exp(λi) - exp(λn))

    if -randexp() < llrce
      llc += llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ) + llrce
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      λi   = λn
      setlλ!(tree, λn)
    end
  end

  return llc, ddλ, ssλ, λi
end





"""
    update_faketip!(tree::T,
                    bi  ::iBffs,
                    eas ::Float64,
                    λa  ::Float64,
                    eds ::Float64,
                    λ1  ::Float64,
                    λ2  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    llc ::Float64,
                    ddλ ::Float64,
                    ssλ ::Float64) where {T <: cT}

Make a `clads` tip proposal.
"""
function update_faketip!(tree::T,
                         bi  ::iBffs,
                         eas ::Float64,
                         λa  ::Float64,
                         eds ::Float64,
                         λ1  ::Float64,
                         λ2  ::Float64,
                         α   ::Float64,
                         σλ  ::Float64,
                         llc ::Float64,
                         ddλ ::Float64,
                         ssλ ::Float64) where {T <: cT}

  @inbounds begin

    λi = lλ(tree)
    ei = e(tree)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)

    # likelihood ratios
    llrce = λn - λi + (eas + ei + eds)*(exp(λi) - exp(λn))

    if -randexp() < llrce
      llc += llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ) + llrce
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      λi   = λn
      setlλ!(tree, λi)
      setλt!(bi, λi)
    end
  end

  return llc, ddλ, ssλ
end



