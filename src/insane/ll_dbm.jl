#=

Diffused Brownian motion likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#



"""
    llik_dbm(Ξ   ::Vector{sTxs},
             γ   ::Float64,
             δt  ::Float64)

Returns the log-likelihood for `sTxs` according to diffused Brownian motion.
"""
function llik_dbm(Ξ   ::Vector{sTxs},
                  γ   ::Float64,
                  δt  ::Float64)
  ll = 0.0
  for ξ in Ξ
    ll += llik_dbm(ξ, γ, δt)
  end

  return ll
end




"""
    llik_dbm_v!(ll  ::Vector{Float64},
                 Ξ   ::Vector{sTxs},
                 γ   ::Float64,
                 δt  ::Float64)

Returns the log-likelihood for `sTxs` according to diffused Brownian motion.
"""
function llik_dbm_v!(ll  ::Vector{Float64},
                     Ξ   ::Vector{sTxs},
                     γ   ::Float64,
                     δt  ::Float64)

  @inbounds begin
    xn = lastindex(Ξ)
    for i in Base.OneTo(xn)
      ll[i] = llik_dbm(Ξ[i], γ, δt)
    end
  end

  return nothing
end




"""
    llik_dbm(tree::sTxs,
             γ   ::Float64,
             δt  ::Float64)

Returns the log-likelihood for a `sTxs` according to diffused Brownian motion.
"""
function llik_dbm(tree::sTxs,
                  γ   ::Float64,
                  δt  ::Float64)

  if def1(tree)
    if def2(tree)
      llik_dbm(tree.d1, γ, δt)                        +
      llik_dbm(tree.d2, γ, δt)                        +
      ll_dbm_b(xv(tree), lσ(tree), γ, δt, fdt(tree))
    else
      llik_dbm(tree.d1, γ, δt)                        +
      ll_dbm_b(xv(tree), lσ(tree), γ, δt, fdt(tree))
    end
  else
    ll_dbm_b(xv(tree), lσ(tree), γ, δt, fdt(tree))
  end
end



"""
    ll_dbm_b(x   ::Array{Float64,1},
             σ   ::Array{Float64,1},
             γ   ::Float64,
             δt  ::Float64,
             fdt ::Float64)

Returns the log-likelihood for a branch according to diffused Brownian motion.
"""
function ll_dbm_b(x   ::Array{Float64,1},
                  σ   ::Array{Float64,1},
                  γ   ::Float64,
                  δt  ::Float64,
                  fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2
    n  = Float64(nI)

    llx = llσ = 0.0
    @turbo for i in Base.OneTo(nI)
      σi   = σ[i]
      σi1  = σ[i+1]
      llσ += (σi1 - σi)^2
      llx += -0.5*(x[i+1] - x[i])^2/(exp(σi1 + σi)*δt) - 0.5*(σi1 + σi)
    end

    # estimate global likelihood
    ll = llx - 0.5*n*log(δt) + 
         llσ*(-0.5/(γ^2*δt)) - 0.5*n*log(γ^2*δt) - n*log(2.0π)

    # add final non-standard `δt`
    if fdt > 0.0
      σi  = σ[nI+1]
      σi1 = σ[nI+2]
      ll += -0.5*(x[nI+2] - x[nI+1])^2/(exp(σi1 + σi)*fdt) - 0.5*(σi1 + σi) - 
             0.5*log(fdt) + 
            (σi1 - σi)^2*(-0.5/(γ^2*fdt)) - 0.5*log(γ^2*fdt) - log(2.0π)
    end
  end

  return ll
end




"""
    ll_dbm_ss_b(x   ::Array{Float64,1},
                σ   ::Array{Float64,1},
                γ   ::Float64,
                δt  ::Float64,
                fdt ::Float64)

Returns the log-likelihood for a branch according to diffused Brownian motion.
"""
function ll_dbm_ss_b(x   ::Array{Float64,1},
                     σ   ::Array{Float64,1},
                     γ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2
    n  = Float64(nI)

    llx = llσ = 0.0
    @turbo for i in Base.OneTo(nI)
      σi   = σ[i]
      σi1  = σ[i+1]
      llσ += (σi1 - σi)^2
      llx += -0.5*(x[i+1] - x[i])^2/(exp(σi1 + σi)*δt) - 0.5*(σi1 + σi)
    end

    # estimate standard squares
    ss = llσ/(2.0*δt)

    # estimate global likelihood
    ll = llx - 0.5*n*log(δt) + 
         llσ*(-0.5/(γ^2*δt)) - 0.5*n*log(γ^2*δt) - n*log(2.0π)

    # add final non-standard `δt`
    if fdt > 0.0
      σi  = σ[nI+1]
      σi1 = σ[nI+2]
      dσ2 = (σi1 - σi)^2
      ll += -0.5*(x[nI+2] - x[nI+1])^2/(exp(σi1 + σi)*fdt) - 0.5*(σi1 + σi) - 
             0.5*log(fdt) + 
            dσ2*(-0.5/(γ^2*fdt)) - 0.5*log(γ^2*fdt) - log(2.0π)
      ss += dσ2/(2.0*fdt)
    end
  end

  return ll, ss
end




"""
    llr_dbm(x   ::Array{Float64,1},
            σp  ::Array{Float64,1},
            σc  ::Array{Float64,1},
            δt  ::Float64,
            fdt ::Float64)

Returns the log-likelihood ratio for a `σ(t)` path proposal.
"""
function llr_dbm(x   ::Array{Float64,1},
                 σp  ::Array{Float64,1},
                 σc  ::Array{Float64,1},
                 δt  ::Float64,
                 fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2

    acr = 0.0
    @turbo for i in Base.OneTo(nI)
      σci  = σc[i]
      σci1 = σc[i+1]
      σpi  = σp[i]
      σpi1 = σp[i+1]
      acr += -(0.5*(x[i+1] - x[i])^2/δt)*
              (1.0/exp(σpi + σpi1) - 1.0/exp(σci + σci1)) + 
               0.5*(σci + σci1 - σpi - σpi1)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      σci  = σc[nI+1]
      σci1 = σc[nI+2]
      σpi  = σp[nI+1]
      σpi1 = σp[nI+2]
      acr += -(0.5*(x[nI+2] - x[nI+1])^2/fdt) *
              (1.0/exp(σpi + σpi1) - 1.0/exp(σci + σci1)) + 
               0.5*(σci + σci1 - σpi - σpi1)
    end
  end

  return acr
end




"""
  llr_scale(Ξ   ::Vector{sTxs},
            s   ::Vector{Float64},
            δt  ::Float64)

Returns the log-likelihood ratio for a `σ(t)` scaled proposal.
"""
function llr_scale(Ξ   ::Vector{sTxs},
                   s   ::Float64,
                   δt  ::Float64)

  @inbounds begin
    xn  = lastindex(Ξ)
    llr = zeros(xn)
    for i in Base.OneTo(xn)
      ξ = Ξ[i]
      llr[i] = llr_scale(xv(ξ), lσ(ξ), s, δt, fdt(ξ))
    end
  end

  return llr
end




"""
    llr_scale(x   ::Array{Float64,1},
              σc  ::Array{Float64,1},
              s   ::Float64,
              δt  ::Float64,
              fdt ::Float64)

Returns the log-likelihood ratio for a `σ(t)` scaled proposal.
"""
function llr_scale(x   ::Array{Float64,1},
                   σc  ::Array{Float64,1},
                   s   ::Float64,
                   δt  ::Float64,
                   fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI  = lastindex(x)-2
    e2s = exp(2.0*s)

    llr = 0.0
    @turbo for i in Base.OneTo(nI)
      σca  = exp(σc[i+1] + σc[i])
      σpa  = σca * e2s
      llr -= 0.5*(x[i+1] - x[i])^2/δt * (1.0/σpa - 1.0/σca)
    end
    llr -= Float64(nI)*s

    # add final non-standard `δt`
    if fdt > 0.0
      σca  = exp(σc[nI+2] + σc[nI+1])
      σpa  = σca * e2s
      llr += -(0.5*(x[nI+2] - x[nI+1])^2/fdt)*(1.0/σpa - 1.0/σca) - s
    end
  end

  return llr
end


