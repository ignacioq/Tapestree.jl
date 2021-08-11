#=

Anagenetic GBM pure-birth Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sim_gbm(t   ::Float64,
            λt  ::Float64,
            α   ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Simulate `iTgbmpb` according to a pure-birth geometric Brownian motion.
"""
function sim_gbm(t   ::Float64,
                 λt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

  λv = Float64[λt]
  bt = 0.0

  while true

    if t <= δt 
      bt  += t

      t   = max(0.0, t)
      λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, t)
        return iTgbmpb(sim_gbm(0.0, λt1, α, σλ, δt, srδt), 
                       sim_gbm(0.0, λt1, α, σλ, δt, srδt), 
                 bt, δt, t, λv)
      end

      return iTgbmpb(bt, δt, t, λv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt + α*δt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divev(λm, δt)
      return iTgbmpb(sim_gbm(t, λt1, α, σλ, δt, srδt), 
                     sim_gbm(t, λt1, α, σλ, δt, srδt), 
              bt, δt, δt, λv)
    end

    λt = λt1
  end
end




"""
    divev(λ::Float64, δt::Float64)

Return true if diversification event.
"""
divev(λ::Float64, δt::Float64) = @fastmath rand() < λ*δt 

