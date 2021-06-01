#=

Anagenetic GBM pure-birth Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sim_gbm(t   ::Float64,
            λt  ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Simulate `iTgbmbd` according to a pure-birth geometric Brownian motion.
"""
function sim_gbm(t   ::Float64,
                 λt  ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

  λv = Float64[λt]
  bt = 0.0

  while true

    t  -= δt
    bt += δt

    λt1 = rnorm(λt, srδt*σλ)

    push!(λv, λt1)
    push!(ts, bt)

    if t <= δt 
      return iTgbmpb(nothing, nothing, bt, δt, nsdt, λv)
    end

    λm = exp(0.5*(λt + λt1))

    if divev(λm, δt)
      return iTgbmpb(sim_gbm(t, λt1, σλ, δt, srδt), 
                     sim_gbm(t, λt1, σλ, δt, srδt), 
              bt, δt, nsdt, λv)
    end

    λt = λt1
  end
end




"""
    divev(λ::Float64, δt::Float64)

Return true if diversification event.
"""
divev(λ::Float64, δt::Float64) = @fastmath rand() < λ*δt 

