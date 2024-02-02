#=

Diffused Brownian motion likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#




"""
    llik_gbm(Ξ   ::Vector{sTxs},
             γ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for `sTxs` according to diffused Brownian motion.
"""
function llik_gbm(Ξ   ::Vector{sTxs},
                  γ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      ll += llik_dbm(Ξ[i], γ, δt, srδt)
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




#=
will need to make likelihood that estimates ssσ at the same time
=#



