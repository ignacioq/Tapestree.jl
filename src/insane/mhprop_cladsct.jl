#=

clads constant extinction MH proposals for internal updates

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    _stem_update!(ξi      ::cTct,
                  eds     ::Float64,
                  λ1      ::Float64,
                  λ2      ::Float64,
                  α       ::Float64,
                  σλ      ::Float64,
                  ϵ       ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ddλ     ::Float64,
                  ssλ     ::Float64,
                  seλ     ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  λ0_prior::NTuple{2,Float64},
                  surv    ::Int64)

Do `clads` update for crown root.
"""
function _stem_update!(ξi      ::cTct,
                       eds     ::Float64,
                       λ1      ::Float64,
                       λ2      ::Float64,
                       α       ::Float64,
                       σλ      ::Float64,
                       ϵ       ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ddλ     ::Float64,
                       ssλ     ::Float64,
                       seλ     ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       λ0_prior::NTuple{2,Float64},
                       surv    ::Int64)

  @inbounds begin
    λi = lλ(ξi)
    ei = e(ξi)

    # node proposal
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], 
                  σλ^2,     σλ^2, λ0_prior[2])

    llrbm = llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ)
    eλr   = (ei + eds) * (exp(λi) - exp(λr))
    llrct = λr - λi + eλr * (1.0 + ϵ)

    if lU < llrct + log(1000.0/mc)

      mp     = m_surv_cladsct(th, λr, α, σλ, ϵ, 1_000, surv)
      llrct += log(mp/mc)

      if -randexp() < llrct
        llc += llrbm + llrct
        prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
        ddλ += 2.0*(λi - λr)
        ssλ += 0.5*(
                (λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
                (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
        seλ -= eλr
        mc   = mp
        λi   = λr
        setlλ!(ξi, λi)
      end
    end
  end

  return llc, prc, ddλ, ssλ, seλ, mc, λi
end




"""
    _crown_update!(ξi      ::cTct,
                   ξ1      ::cTct,
                   ξ2      ::cTct,
                   α       ::Float64,
                   σλ      ::Float64,
                   ϵ       ::Float64,
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
function _crown_update!(ξi      ::cTct,
                        ξ1      ::cTct,
                        ξ2      ::cTct,
                        α       ::Float64,
                        σλ      ::Float64,
                        ϵ       ::Float64,
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
    mp  = m_surv_cladsct(th, λr, α, σλ, ϵ, 1_000, surv)
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
    _update_internal!(tree::cTct,
                      bi  ::iBffs,
                      eas ::Float64,
                      λa  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      ϵ   ::Float64,
                      eds ::Float64,
                      λ1  ::Float64,
                      λ2  ::Float64,
                      llc ::Float64,
                      ddλ ::Float64,
                      ssλ ::Float64,
                      seλ ::Float64,
                      ter ::Bool)


Do `clads` internal rate updates on a decoupled tree recursively.
"""
function _update_internal!(tree::cTct,
                           bi  ::iBffs,
                           eas ::Float64,
                           λa  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           ϵ   ::Float64,
                           eds ::Float64,
                           λ1  ::Float64,
                           λ2  ::Float64,
                           llc ::Float64,
                           ddλ ::Float64,
                           ssλ ::Float64,
                           seλ ::Float64,
                           ter ::Bool)

  if def1(tree)
    llc, ddλ, ssλ, seλ, λa = 
      update_triad!(tree, eas, λa, α, σλ, ϵ, llc, ddλ, ssλ, seλ)

    llc, ddλ, ssλ, seλ, λx =
      _update_internal!(tree.d1, bi, 0.0, λa, α, σλ, ϵ, eds, λ1, λ2, 
        llc, ddλ, ssλ, seλ, ter)
    llc, ddλ, ssλ, seλ, λx =
      _update_internal!(tree.d2, bi, 0.0, λa, α, σλ, ϵ, eds, λ1, λ2, 
        llc, ddλ, ssλ, seλ, ter)
  else 
    # if real tip
    if !isfix(tree) || ter
      llc, ddλ, ssλ, seλ = 
        update_tip!(tree, eas, λa, 0.0, α, σλ, ϵ, llc, ddλ, ssλ, seλ)
    # if leads to non-speciation
    elseif isnan(λ1)
      llc, ddλ, ssλ, seλ = 
        update_tip!(tree, eas, λa, eds, α, σλ, ϵ, llc, ddλ, ssλ, seλ)
    # if leads to eventual speciation
    else
      llc, ddλ, ssλ, seλ = 
        update_faketip!(tree, bi, eas, λa, eds, λ1, λ2, α, σλ, ϵ, 
          llc, ddλ, ssλ, seλ)
    end
  end

  return llc, ddλ, ssλ, seλ, λa
end




"""
    update_triad!(tree::cTct,
                  eas ::Float64,
                  λa  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64,
                  seλ ::Float64)

Make a trio proposal for clads.
"""
function update_triad!(tree::cTct,
                       eas ::Float64,
                       λa  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       ϵ   ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       seλ ::Float64)

  @inbounds begin

    λi = lλ(tree)
    λ1 = lλ(tree.d1)
    λ2 = lλ(tree.d2)
    ei = e(tree)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)

    # likelihood ratios
    llrbm = llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ)
    eλr   = (ei + eas) * (exp(λi) - exp(λn))
    llrct = λn - λi + eλr * (1.0 + ϵ)

    if -randexp() < llrct
      llc += llrbm + llrct
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      seλ -= eλr
      λi   = λn
      setlλ!(tree, λn)
    end
  end

  return llc, ddλ, ssλ, seλ, λi
end




"""
    update_tip!(tree::cTct,
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
function update_tip!(tree::cTct,
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

  @inbounds begin

    λi = lλ(tree)
    ei = e(tree)

    # node proposal
    λn = rnorm(λa + α, σλ)

    # likelihood ratios
    llrbm = llrdnorm_x(λn, λi, λa + α, σλ^2)
    eλr   = (eas + ei + eds) * (exp(λi) - exp(λn))
    llrct = eλr * (1.0 + ϵ)

    if isextinct(tree)
      llrct += λn - λi
    end

    if -randexp() < llrct
      llc += llrbm + llrct
      ddλ += λn - λi
      ssλ += 0.5*((λn - λa - α)^2 - (λi - λa - α)^2)
      seλ -= eλr
      setlλ!(tree, λn)
    end
  end

  return llc, ddλ, ssλ, seλ
end




"""
    update_faketip!(tree::cTct,
                    bi  ::iBffs,
                    eas ::Float64,
                    λa  ::Float64,
                    eds ::Float64,
                    λ1  ::Float64,
                    λ2  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    llc ::Float64,
                    ddλ ::Float64,
                    ssλ ::Float64,
                    seλ ::Float64)

Make a `clads` final internal branch proposal.
"""
function update_faketip!(tree::cTct,
                         bi  ::iBffs,
                         eas ::Float64,
                         λa  ::Float64,
                         eds ::Float64,
                         λ1  ::Float64,
                         λ2  ::Float64,
                         α   ::Float64,
                         σλ  ::Float64,
                         ϵ   ::Float64,
                         llc ::Float64,
                         ddλ ::Float64,
                         ssλ ::Float64,
                         seλ ::Float64)

  @inbounds begin

    λi = lλ(tree)
    ei = e(tree)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)

    # likelihood ratios
    llrbm = llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ)
    eλr   = (eas + ei + eds) * (exp(λi) - exp(λn))
    llrct = λn - λi + eλr * (1.0 + ϵ)

    if -randexp() < llrct
      llc += llrbm + llrct
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      seλ -= eλr
      setlλ!(tree, λn)
      setλt!(bi, λn)
    end
  end

  return llc, ddλ, ssλ, seλ
end



