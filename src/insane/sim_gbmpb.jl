#=

Anagenetic GBM pure-birth Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sim_gbm(t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
                 σλ ::Float64,
                 σ_μ ::Float64,
                 dt  ::Float64,
                 srdt::Float64)

Simulate `iTree` according to a geometric Brownian motion.
"""
function sim_gbm(t   ::Float64,
                 λt  ::Float64,
                 σλ  ::Float64,
                 dt  ::Float64,
                 srdt::Float64)

  λv = Float64[λt]
  bt = 0.0
  ts = Float64[bt]

  while true

    t  -= dt
    bt += dt

    λt1 = rnorm(λt, srdt*σλ)

    push!(λv, λt1)
    push!(ts, bt)

    if 0.0 > t - dt
      return iTgbmpb(nothing, nothing, bt, ts, λv)
    end

    λm = exp(0.5*(λt+λt1))

    if divev(λm, dt)
      return iTgbmpb(sim_gbm(t, λt1, σλ, dt, srdt), 
                     sim_gbm(t, λt1, σλ, dt, srdt), 
              bt, ts, λv)
    end

    λt = λt1
  end
end




"""
    divev(λ::Float64, dt::Float64)

Return true if diversification event.
"""
divev(λ::Float64, dt::Float64) = @fastmath rand() < λ*dt 

