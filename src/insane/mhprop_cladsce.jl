#=

clads constant extinction MH proposals for internal updates

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    _stem_update!(ξi      ::cTce,
                  ξ1      ::cTce,
                  ξ2      ::cTce,
                  α       ::Float64,
                  σλ      ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ddλ     ::Float64,
                  ssλ     ::Float64,
                  λ0_prior::NTuple{2,Float64})

Do `clads` update for crown root.
"""
function _stem_update!(ξi      ::cTce,
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

    # node proposal
    λr = trioprop(λ1 - α, λ2 - α, λ0_prior[1], 
                  σλ^2,     σλ^2, λ0_prior[2])

    llrbm = llrdnorm2_μ(λ1, λ2, λr + α, λi + α, σλ)
    llrce = λr - λi + ei*(exp(λi) - exp(λr))

    if lU < llrce + log(1000.0/mc)

      mp     = m_surv_cladsce(th, λr, α, σλ, μ, 1_000, surv)
      llrce += log(mp/mc)

      if -randexp() < llrce
        llc += llrbm + llrce
        prc += llrdnorm_x(λr, λi, λ0_prior[1], λ0_prior[2])
        ddλ += 2.0*(λi - λr)
        ssλ += 0.5*(
                (λ1 - λr - α)^2 + (λ2 - λr - α)^2 - 
                (λ1 - λi - α)^2 - (λ2 - λi - α)^2)
        mc  = mp
        setlλ!(ξi, λr)
      end
    end
  end

  return llc, prc, ddλ, ssλ, mc
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




