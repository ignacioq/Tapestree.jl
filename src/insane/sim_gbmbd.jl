#=

Anagenetic GBM birth-death Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sim_gbm(t   ::Float64,
            λt  ::Float64,
            μt  ::Float64,
            σ_λ ::Float64,
            σ_μ ::Float64,
            dt  ::Float64,
            srdt::Float64)

Simulate `iTgbmbd` according to a geometric Brownian motion.
"""
function sim_gbm(t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
                 σ_λ ::Float64,
                 σ_μ ::Float64,
                 dt  ::Float64,
                 srdt::Float64)

  λv = Float64[λt]
  μv = Float64[μt]
  bt = 0.0
  ts = Float64[bt]

  while true

    t  -= dt
    bt += dt

    λt1 = rnorm(λt, srdt*σ_λ)
    μt1 = rnorm(μt, srdt*σ_μ)

    push!(λv, λt1)
    push!(μv, μt1)
    push!(ts, bt)

    if 0.0 > t - dt
      return iTgbmbd(nothing, nothing, bt, false, false, ts, λv, μv)
    end

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, dt)
      # if speciation
      if λorμ(λm, μm)
        return iTgbmbd(sim_gbm(t, λt1, μt1, σ_λ, σ_μ, dt, srdt), 
                       sim_gbm(t, λt1, μt1, σ_λ, σ_μ, dt, srdt), 
                bt, false, false, ts, λv, μv)
      # if extinction
      else
        return iTgbmbd(nothing, nothing, bt, true, false, ts, λv, μv)
      end
    end

    λt = λt1
    μt = μt1
  end
end




"""
    divev(λ::Float64, μ::Float64, dt::Float64)

Return true if diversification event.
"""
divev(λ::Float64, μ::Float64, dt::Float64) = @fastmath rand() < (λ + μ)*dt 




"""
    rnorm(μ::Float64, σ::Float64)

Generate a normal variable with mean `μ` and variance `σ`.
"""
rnorm(μ::Float64, σ::Float64) = @fastmath randn()*σ + μ




