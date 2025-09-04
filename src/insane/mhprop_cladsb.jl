#=

clads pure-birth MH proposals for internal updates

Ignacio Quintero Mächler

t(-_-t)

Created 15 07 2025
=#




"""
    _stem_update!(ξi      ::cTb,
                  ξ1      ::cTb,
                  ξ2      ::cTb,
                  α       ::Float64,
                  σλ      ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ddλ     ::Float64,
                  ssλ     ::Float64,
                  λ0_prior::NTuple{2,Float64})

Do `clads` update for crown root.
"""
function _stem_update!(ξi      ::cTb,
                       λ1      ::Float64,
                       λ2      ::Float64,
                       α       ::Float64,
                       σλ      ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ddλ     ::Float64,
                       ssλ     ::Float64,
                       λ0_prior::NTuple{2,Float64})

  @inbounds begin
    λi = lλ(ξi)
    ei = e(ξi)

    # node proposal
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], 
                  σλ^2,     σλ^2, λ0_prior[2])

    llrbm = llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ)
    llrpb = λr - λi + ei*(exp(λi) - exp(λr))

    if -randexp() < llrpb
      llc += llrbm + llrpb
      prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
      ddλ += 2.0*(λi - λr)
      ssλ += 0.5*(
              (λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
              (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      setlλ!(ξi, λr)
    end
  end

  return llc, prc, ddλ, ssλ
end




"""
    _crown_update!(ξi      ::cTb,
                   ξ1      ::cTb,
                   ξ2      ::cTb,
                   α       ::Float64,
                   σλ      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   ddλ     ::Float64,
                   ssλ     ::Float64,
                   λ0_prior::NTuple{2,Float64})

Do `clads` update for crown root.
"""
function _crown_update!(ξi      ::cTb,
                        ξ1      ::cTb,
                        ξ2      ::cTb,
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
    _update_internal!(tree::cTb,
                      λa  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      llc ::Float64,
                      ddλ ::Float64,
                      ssλ ::Float64,
                      ter ::Bool)

Do `clads` internal rate updates on a decoupled tree recursively.
"""
function _update_internal!(tree::T,
                           λa  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           llc ::Float64,
                           ddλ ::Float64,
                           ssλ ::Float64,
                           ter ::Bool) where {T <: cT}

  if def1(tree)
    llc, ddλ, ssλ, λa = 
      update_triad!(tree, λa, α, σλ, llc, ddλ, ssλ)

    llc, ddλ, ssλ, λx =
      _update_internal!(tree.d1, λa, α, σλ, llc, ddλ, ssλ, ter)
    llc, ddλ, ssλ, λx =
      _update_internal!(tree.d2, λa, α, σλ, llc, ddλ, ssλ, ter)
  elseif !isfix(tree) || ter
    if !isnan(λa)
      llc, ddλ, ssλ = 
        update_tip!(tree, 0.0, λa, α, σλ, llc, ddλ, ssλ)
    end
  end

  return llc, ddλ, ssλ, λa
end




"""
    update_tip!(tree::T,
                eas ::Float64,
                λa  ::Float64,
                eds ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                llc ::Float64,
                ddλ ::Float64,
                ssλ ::Float64) where {T <: cT}

Make a `clads` tip proposal.
"""
function update_tip!(tree::T,
                     eas ::Float64,
                     λa  ::Float64,
                     eds ::Float64,
                     α   ::Float64,
                     σλ  ::Float64,
                     llc ::Float64,
                     ddλ ::Float64,
                     ssλ ::Float64) where {T <: cT}

  @inbounds begin

    λi = lλ(tree)
    ei = e(tree)

    # node proposal
    λn = rnorm(λa + α, σλ)

    # likelihood ratios
    llrbm = llrdnorm_x(λn, λi, λa + α, σλ^2)
    llrpb = (eas + ei + eds)*(exp(λi) - exp(λn))

    if -randexp() < llrpb
      llc += llrbm + llrpb
      ddλ += λn - λi
      ssλ += 0.5*((λn - λa - α)^2 - (λi - λa - α)^2)
      setlλ!(tree, λn)
    end
  end

  return llc, ddλ, ssλ
end




"""
    update_triad!(ξi  ::T,
                  ξ1  ::T,
                  ξ2  ::T,
                  λa  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64) where {T <: cT}

Make a `clads` trio proposal.
"""
function update_triad!(ξi  ::T,
                       ξ1  ::T,
                       ξ2  ::T,
                       λa  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64) where {T <: cT}

  @inbounds begin
    λi = lλ(ξi)
    λ1 = lλ(ξ1)
    λ2 = lλ(ξ2)
    ei = e(ξi)

    # node proposal
    λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)

    # likelihood ratios
    llrbm = llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ)
    llrpb = λn - λi + ei*(exp(λi) - exp(λn))

    if -randexp() < llrpb
      llc += llrbm + llrpb
      ddλ += (λi - λn)
      ssλ += 0.5*(
              (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
              (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
      λi   = λn
      setlλ!(ξi, λn)
    end
  end

  return llc, ddλ, ssλ, λi
end




"""
    update_triad!(tree::T,
                  λa  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ddλ ::Float64,
                  ssλ ::Float64) where {T <: cT}

Make a trio proposal for clads.
"""
function update_triad!(tree::T,
                       λa  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       ddλ ::Float64,
                       ssλ ::Float64) where {T <: cT}

  @inbounds begin

    λi = lλ(tree)

    if !isnan(λa)
      λ1 = lλ(tree.d1)
      λ2 = lλ(tree.d2)
      ei = e(tree)

      # node proposal
      λn = trioprop(λa + α, λ1 - α, λ2 - α, σλ)

      # likelihood ratios
      llrbm = llrdnorm3(λa + α, λ1 - α, λ2 - α, λn, λi, σλ)
      llrpb = λn - λi + ei*(exp(λi) - exp(λn))

      if -randexp() < llrpb
        llc += llrbm + llrpb
        ddλ += (λi - λn)
        ssλ += 0.5*(
                (λn - λa - α)^2 + (λ1 - λn - α)^2 + (λ2 - λn - α)^2 -
                (λi - λa - α)^2 - (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
        λi   = λn
        setlλ!(tree, λn)
      end
    end
  end

  return llc, ddλ, ssλ, λi
end



