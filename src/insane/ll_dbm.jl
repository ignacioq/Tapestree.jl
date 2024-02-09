#=

Diffused Brownian motion likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#



"""
    llik_dbm(Ξ   ::Vector{sTxs},
             γ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for `sTxs` according to diffused Brownian motion.
"""
function llik_dbm(Ξ   ::Vector{sTxs},
                  γ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)
  ll = 0.0
  for ξ in Ξ
    ll += llik_dbm(ξ, γ, δt, srδt)
  end

  return ll
end




"""
    llik_dbm_v(Ξ   ::Vector{sTxs},
               γ   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Returns the log-likelihood for `sTxs` according to diffused Brownian motion.
"""
function llik_dbm_v(Ξ   ::Vector{sTxs},
                    γ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  @inbounds begin
    xn = lastindex(Ξ)
    ll = zeros(xn)
    for i in Base.OneTo(xn)
      ll[i] = llik_dbm(Ξ[i], γ, δt, srδt)
    end
  end

  return ll
end




"""
    llik_dbm(tree::sTxs,
             γ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `sTxs` according to diffused Brownian motion.
"""
function llik_dbm(tree::sTxs,
                  γ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  if def1(tree)
    if def2(tree)
      llik_dbm(tree.d1, γ, δt, srδt)                        +
      llik_dbm(tree.d2, γ, δt, srδt)                        +
      ll_dbm_b(xv(tree), lσ(tree), γ, δt, fdt(tree), srδt)
    else
      llik_dbm(tree.d1, γ, δt, srδt)                        +
      ll_dbm_b(xv(tree), lσ(tree), γ, δt, fdt(tree), srδt)
    end
  else
    ll_dbm_b(xv(tree), lσ(tree), γ, δt, fdt(tree), srδt)
  end
end



"""
    ll_dbm_b(x   ::Array{Float64,1},
             σ   ::Array{Float64,1},
             γ   ::Float64,
             δt  ::Float64,
             fdt ::Float64,
             srδt::Float64)

Returns the log-likelihood for a branch according to diffused Brownian motion.
"""
function ll_dbm_b(x   ::Array{Float64,1},
                  σ   ::Array{Float64,1},
                  γ   ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2

    llx = 0.0
    llσ = 0.0
    @turbo for i in Base.OneTo(nI)
      σi   = σ[i]
      σi1  = σ[i+1]
      llσ += (σi1 - σi)^2
      σa   = exp(0.5*(σi1 + σi))
      llx += -0.5*((x[i+1] - x[i])/(σa*srδt))^2 - log(σa*srδt)
    end

    # estimate global likelihood
    nct  = Float64(nI)*(0.5*log(2.0π))
    llx -= nct
    ll   = llx +
           llσ*(-0.5/((γ*srδt)^2)) - Float64(nI)*(log(γ*srδt)) - nct

    # add final non-standard `δt`
    if fdt > 0.0
      σi    = σ[nI+1]
      σi1   = σ[nI+2]
      σa    = exp(0.5*(σi1 + σi))
      srfdt = sqrt(fdt)
      ll   += ldnorm_bm(x[nI+2], x[nI+1], srfdt*σa) +
              ldnorm_bm(σi1,          σi, srfdt*γ)
    end
  end

  return ll
end




"""
    ll_dbm_ss_b(x   ::Array{Float64,1},
                σ   ::Array{Float64,1},
                γ   ::Float64,
                δt  ::Float64,
                fdt ::Float64,
                srδt::Float64)

Returns the log-likelihood for a branch according to diffused Brownian motion.
"""
function ll_dbm_ss_b(x   ::Array{Float64,1},
                     σ   ::Array{Float64,1},
                     γ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2

    llx = 0.0
    llσ = 0.0
    @turbo for i in Base.OneTo(nI)
      σi   = σ[i]
      σi1  = σ[i+1]
      llσ += (σi1 - σi)^2
      σa   = exp(0.5*(σi1 + σi))
      llx += -0.5*((x[i+1] - x[i])/(σa*srδt))^2 - log(σa*srδt)
    end

    # estimate standard squares
    ss = llσ/(2.0*δt)
    
    # estimate global likelihood
    nct  = Float64(nI)*(0.5*log(2.0π))
    llx -= nct
    ll   = llx +
           llσ*(-0.5/((γ*srδt)^2)) - Float64(nI)*(log(γ*srδt)) - nct

    # add final non-standard `δt`
    if fdt > 0.0
      σi    = σ[nI+1]
      σi1   = σ[nI+2]
      σa    = exp(0.5*(σi1 + σi))
      srfdt = sqrt(fdt)
      ll   += ldnorm_bm(x[nI+2], x[nI+1], srfdt*σa) +
              ldnorm_bm(σi1,          σi, srfdt*γ)
      ss   += (σi1 - σi)^2/(2.0*fdt)
    end
  end

  return ll, ss
end




"""
    llr_dbm(x   ::Array{Float64,1},
            σp  ::Array{Float64,1},
            σc  ::Array{Float64,1},
            γ   ::Float64,
            δt  ::Float64,
            fdt ::Float64,
            srδt::Float64)

Returns the log-likelihood ratio for a `σ(t)` path proposal.
"""
function llr_dbm(x   ::Array{Float64,1},
                 σp  ::Array{Float64,1},
                 σc  ::Array{Float64,1},
                 γ   ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2

    llr = 0.0
    @turbo for i in Base.OneTo(nI)
      σpa  = exp(0.5*(σp[i+1] + σp[i]))
      σca  = exp(0.5*(σc[i+1] + σc[i]))
      llr += -(0.5*(x[i+1] - x[i])^2/δt)*(1.0/σpa^2 - 1.0/σca^2) + log(σca/σpa)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      σpa   = exp(0.5*(σp[nI+2] + σp[nI+1]))
      σca   = exp(0.5*(σc[nI+2] + σc[nI+1]))
      llr  += -(0.5*(x[nI+2] - x[nI+1])^2/fdt)*(1.0/σpa^2 - 1.0/σca^2) +
               log(σca/σpa)
    end
  end

  return llr
end


