#=

clads fossilized birth-death MH proposals for internal updates

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    _stem_update!(bix     ::Int64,
                  Ξ       ::Vector{acTfbd},
                  idf     ::Vector{iBffs},
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
function _stem_update!(bix     ::Int64,
                       Ξ       ::Vector{acTfbd},
                       idf     ::Vector{iBffs},
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
    ξi = Ξ[bix]

    """
    here: proposal comes from relative likelihood for sh
    """


    # current rates
    λi, μi = lλ(ξi), lμ(ξi)

    ### node proposal
    # find bud number and sum
    ei, nb, sλ, sμ, iμ = sumλμbuds(bix, Ξ, idf, 0.0, 0.0, 0.0, 0.0, false)

    ## speciation
    # prior
    λ0m, λ0s = λ0_prior
    # conditional normal for proposal
    λs2 = 1.0/(1.0/λ0s + nb/σλ^2)
    λm  = λs2 * (λ0m/λ0s + (sλ - nb*αλ)/σλ^2)
    λr = rnorm(λm, sqrt(λs2))

    ## extinction
    # prior
    μ0m, μ0s = μ0_prior
    # conditional normal for proposal
    μs2 = 1.0/(1.0/μ0s + nb/σμ^2)
    μm  = μs2 * (μ0m/μ0s + (sμ - nb*αμ)/σμ^2)
    μr = rnorm(μm, sqrt(μs2))

    # likelihood ratio
    llr = nb*(λr - λi) + ei * (exp(λi) - exp(λr) + exp(μi) - exp(μr)) + 
          iμ ? 0.0 : (μr - μi)

    lU = -randexp()

    if lU < llr + log(1000.0/mc)

      # survival ratio
      mp   = m_surv_acladsfbd(th, λr, μr, αλ, αμ, σλ, σμ, 1_000, surv)
      llr += log(mp/mc)

      if lU < llr

        λi, μi = lλ(Ξ[bix]), lμ(Ξ[bix])
        λs2 = 1.0/(nb/σλ^2)
        μs2 = 1.0/(nb/σμ^2)

        llc += llrdnorm_x(λr, λi, λs2 * (sλ - nb*αλ)/σλ^2, λs2) + 
               llrdnorm_x(μr, μi, μs2 * (sμ - nb*αμ)/σμ^2, μs2) + 
               llr
        prc += llrdnorm_x(λr, λi, λ0m, λ0s) + 
               llrdnorm_x(μr, μi, μ0m, μ0s)
        ddλ += nb*(λi - λr)
        ddμ += nb*(μi - μr)
        ssλ += nb*(0.5*(λr^2 + (λi - λr)*sλ - λi^2) + αλ*(λr - λi))
        ssμ += nb*(0.5*(μr^2 + (μi - μr)*sμ - μi^2) + αμ*(μr - μi))
        mc   = mp

        setlλ!(ξi, λr)
        setlμ!(ξi, μr)
      end
    end
  end


  return llc, prc, ddλ, ddμ, ssλ, ssμ, mc
end




"""
    _crown_update!(bix     ::Int64,
                   Ξ       ::Vector{acTfbd},
                   idf     ::Vector{iBffs},
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
function _crown_update!(bix     ::Int64,
                        Ξ       ::Vector{acTfbd},
                        idf     ::Vector{iBffs},
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
    ξi = Ξ[bix]

    # current rates
    λi, μi = lλ(ξi), lμ(ξi)




    """
    here: proposal comes from relative likelihood for sh
    """




    ### node proposal
    # find bud number and sum
    ei, nb, sλ, sμ, iμ = sumλμbuds(bix, Ξ, idf, 0.0, 0.0, 0.0, 0.0, false)

    ## speciation
    # prior
    λ0m, λ0s = λ0_prior
    # conditional normal for proposal
    λs2 = 1.0/(1.0/λ0s + nb/σλ^2)
    λm  = λs2 * (λ0m/λ0s + (sλ - nb*αλ)/σλ^2)
    λr = rnorm(λm, sqrt(λs2))

    ## extinction
    # prior
    μ0m, μ0s = μ0_prior
    # conditional normal for proposal
    μs2 = 1.0/(1.0/μ0s + nb/σμ^2)
    μm  = μs2 * (μ0m/μ0s + (sμ - nb*αμ)/σμ^2)
    μr = rnorm(μm, sqrt(μs2))

    # likelihood ratio
    llr = (nb - 1.0)*(λr - λi) + ei * (exp(λi) - exp(λr) + exp(μi) - exp(μr)) + 
           iμ ? 0.0 : (μr - μi)

    lU = -randexp()

    if lU < llr + log(1000.0/mc)

      # survival ratio
      mp   = m_surv_acladsfbd(th, λr, μr, αλ, αμ, σλ, σμ, 1_000, surv)
      llr += log(mp/mc)

      if lU < llr

        λi, μi = lλ(Ξ[bix]), lμ(Ξ[bix])
        λs2 = 1.0/(nb/σλ^2)
        μs2 = 1.0/(nb/σμ^2)

        llc += llrdnorm_x(λr, λi, λs2 * (sλ - nb*αλ)/σλ^2, λs2) + 
               llrdnorm_x(μr, μi, μs2 * (sμ - nb*αμ)/σμ^2, μs2) + 
               llr
        prc += llrdnorm_x(λr, λi, λ0m, λ0s) + 
               llrdnorm_x(μr, μi, μ0m, μ0s)
        ddλ += nb*(λi - λr)
        ddμ += nb*(μi - μr)
        ssλ += nb*(0.5*(λr^2 + (λi - λr)*sλ - λi^2) + αλ*(λr - λi))
        ssμ += nb*(0.5*(μr^2 + (μi - μr)*sμ - μi^2) + αμ*(μr - μi))
        mc   = mp


        """
        here: create setters
        """


        setlλ!(ξi, λr)
        setlμ!(ξi, μr)
      end
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
                           ssμ ::Float64) where {T <: cT}

  if def1(tree)
    if def2(tree)
      if isfinite(λa)
        llc, ddλ, ddμ, ssλ, ssμ, λa, μa = 
          update_triad!(tree, eas, λa, μa, αλ, αμ, σλ, σμ, 
            llc, ddλ, ddμ, ssλ, ssμ)
      else
        λa, μa = lλ(tree), lμ(tree)
      end

      llc, ddλ, ddμ, ssλ, ssμ, λx, μx =
        _update_internal!(tree.d1, bi, 0.0, λa, μa, αλ, αμ, σλ, σμ, eds, λ1, λ2, 
          μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ)
      llc, ddλ, ddμ, ssλ, ssμ, λx, μx =
        _update_internal!(tree.d2, bi, 0.0, λa, μa, αλ, αμ, σλ, σμ, eds, λ1, λ2, 
          μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ)
    else
      llc, ddλ, ddμ, ssλ, ssμ, λx, μx =
        _update_internal!(tree.d1, bi, eas + e(tree), λa, μa, αλ, αμ, σλ, σμ, 
          eds, λ1, λ2, μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ)
    end
  else 
    if isfix(tree) 
      # if leads to eventual speciation
      if isfinite(λ1)
        llc, ddλ, ddμ, ssλ, ssμ = 
            update_faketip!(tree, bi, eas, λa, μa, eds, λ1, λ2, μ1, μ2, 
              αλ, αμ, σλ, σμ, llc, ddλ, ddμ, ssλ, ssμ)
      # if leads to non-speciation or eventual extinction
      else 
        llc, ddλ, ddμ, ssλ, ssμ = 
          update_tip!(tree, eas, λa, μa, eds, isfinite(μ1), αλ, αμ, σλ, σμ, 
            llc, ddλ, ddμ, ssλ, ssμ)
      end
    else
      llc, ddλ, ddμ, ssλ, ssμ = 
        update_tip!(tree, eas, λa, μa, 0.0, false, αλ, αμ, σλ, σμ, 
          llc, ddλ, ddμ, ssλ, ssμ)
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
    update_tip!(tree::acTfbd,
                eas ::Float64,
                λa  ::Float64,
                μa  ::Float64,
                eds ::Float64,
                eμ  ::Bool,
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
function update_tip!(tree::acTfbd,
                     eas ::Float64,
                     λa  ::Float64,
                     μa  ::Float64,
                     eds ::Float64,
                     eμ  ::Bool,
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

    if isextinct(tree) || eμ
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



