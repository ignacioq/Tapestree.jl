#=

Anagenetic GBM pure-birth MCMC MH proposals

Ignacio Quintero Mächler

t(-_-t)

Created 14 11 2021
=#




"""
    _stem_update!(ξi      ::cTpb,
                  ξ1      ::cTpb,
                  ξ2      ::cTpb,
                  α       ::Float64,
                  σλ      ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ddλ     ::Float64,
                  ssλ     ::Float64,
                  irλ     ::Float64,
                  λ0_prior::NTuple{2,Float64})

Do gbm update for crown root.
"""
function _stem_update!(ξi      ::cTpb,
                       λ1      ::Float64,
                       λ2      ::Float64,
                       α       ::Float64,
                       σλ      ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ddλ     ::Float64,
                       ssλ     ::Float64,
                       irλ     ::Float64,
                       λ0_prior::NTuple{2,Float64})

  @inbounds begin
    λi = lλ(ξi)
    ei = e(ξi)

    # node proposal
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], 
                  σλ^2,     σλ^2, λ0_prior[2])

    llrbm = llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ)
    irrλ  = ei*(exp(λi) - exp(λr))
    llrpb = λr - λi + irrλ

    if -randexp() < llrpb
      llc += llrbm + llrpb
      prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
      ddλ += 2.0*(λi - λr)
      ssλ += 0.5*(
              (λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
              (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      irλ -= irrλ
      setlλ!(ξi, λr)
    end
  end

  return llc, prc, ddλ, ssλ, irλ
end




"""
    _crown_update!(ξi      ::cTpb,
                   ξ1      ::cTpb,
                   ξ2      ::cTpb,
                   α       ::Float64,
                   σλ      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   ddλ     ::Float64,
                   ssλ     ::Float64,
                   λ0_prior::NTuple{2,Float64})

Do gbm update for crown root.
"""
function _crown_update!(ξi      ::cTpb,
                        ξ1      ::cTpb,
                        ξ2      ::cTpb,
                        α       ::Float64,
                        σλ      ::Float64,
                        llc     ::Float64,
                        prc     ::Float64,
                        ddλ     ::Float64,
                        ssλ     ::Float64,
                        λ0_prior::NTuple{2,Float64})

  @inbounds begin
    λi = lλ(ξi)
    λ1 = lλ(ξ1)
    λ2 = lλ(ξ2)

    # node proposal
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], 
                  σλ^2,   σλ^2,   λ0_prior[2])

    # updated automatically
    llc += llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ)
    prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
    ddλ += 2.0*(λi - λr)
    ssλ += 0.5*((λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
                (λ1 - λi - α)^2 - (λ2 - λi - α)^2)

    setlλ!(ξi, λr)
  end

  return llc, prc, ddλ, ssλ
end




"""
    _update_internal!(tree::cTpb,
                      λa  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      llc ::Float64,
                      ddλ ::Float64,
                      ssλ ::Float64,
                      irλ ::Float64,
                      ter ::Bool)

Do clads internal rate updates on a decoupled tree recursively.
"""
function _update_internal!(tree::cTpb,
                           λa  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           llc ::Float64,
                           ddλ ::Float64,
                           ssλ ::Float64,
                           irλ ::Float64,
                           ter ::Bool)

  if def1(tree)

    llc, ddλ, ssλ, irλ, λa = 
      update_triad!(tree, λa, α, σλ, llc, ddλ, ssλ, irλ)

    llc, ddλ, ssλ, irλ, λa =
      _update_internal!(tree.d1, λa, α, σλ, llc, ddλ, ssλ, irλ, ter)
    llc, ddλ, ssλ, irλ, λa =
      _update_internal!(tree.d2, λa, α, σλ, llc, ddλ, ssλ, irλ, ter)
  elseif !isfix(tree) || ter
    llc, ddλ, ssλ, irλ = 
      update_tip!(tree, λa, α, σλ, llc, ddλ, ssλ, irλ)
  end

  return llc, ddλ, ssλ, irλ
end




"""
    update_tip!(tree::cTpb,
                λa  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                llc ::Float64,
                ddλ ::Float64,
                ssλ ::Float64,
                irλ ::Float64)

Make a `clads` tip proposal.
"""
function update_tip!(tree::cTpb,
                     λa  ::Float64,
                     α   ::Float64,
                     σλ  ::Float64,
                     llc ::Float64,
                     ddλ ::Float64,
                     ssλ ::Float64,
                     irλ ::Float64)

  @inbounds begin

    λi = lλ(tree)
    ei = e(tree)

    # node proposal
    λn = rnorm(λa + α, σλ)

    # likelihood ratios
    llrbm = llrdnorm_x(λn, λi, λa + α, σλ^2)
    irrλ  = ei*(exp(λi) - exp(λn))
    llrpb = irrλ

    if -randexp() < llrpb
      llc += llrbm + llrpb
      ddλ += λn - λi
      ssλ += 0.5*((λn - λa - α)^2 - (λi - λa - α)^2)
      irλ -= irrλ
      setlλ!(tree, λn)
    end
  end

  return llc, ddλ, ssλ, irλ
end




"""
    update_triad!(ξi  ::cTpb,
                  ξ1  ::cTpb,
                  ξ2  ::cTpb,
                  λa  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64,
                  irλ ::Float64)

Make a `gbm` trio proposal.
"""
function update_triad!(ξi  ::cTpb,
                       ξ1  ::cTpb,
                       ξ2  ::cTpb,
                       λa  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       irλ ::Float64)

  @inbounds begin
    λi = lλ(ξi)
    λ1 = lλ(ξ1)
    λ2 = lλ(ξ2)
    ei = e(ξi)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)

    # likelihood ratios
    llrbm = llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ)
    irrλ  = ei*(exp(λi) - exp(λn))
    llrpb = λn - λi + irrλ

    if -randexp() < llrpb
      llc += llrbm + llrpb
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      irλ -= irrλ
      λi   = λn
      setlλ!(ξi, λn)
    end
  end

  return llc, ddλ, ssλ, irλ, λi
end




"""
    update_triad!(tree::cTpb,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64,
                  irλ ::Float64)

Make a trio proposal for clads.
"""
function update_triad!(tree::cTpb,
                       λa  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64,
                       irλ ::Float64)

  @inbounds begin

    λi = lλ(tree)
    λ1 = lλ(tree.d1)
    λ2 = lλ(tree.d2)
    ei = e(tree)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)

    # likelihood ratios
    llrbm = llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ)
    irrλ  = ei*(exp(λi) - exp(λn))
    llrpb = λn - λi + irrλ

    if -randexp() < llrpb
      llc += llrbm + llrpb
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      irλ -= irrλ
      λi   = λn
      setlλ!(ξi, λn)
    end
  end

  return llc, ddλ, ssλ, irλ, λi
end



