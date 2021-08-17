#=

Anagenetic GBM pure-birth Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sim_gbmpb(t   ::Float64,
              λt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64)

Simulate `iTgbmpb` according to a pure-birth geometric Brownian motion.
"""
function sim_gbmpb(t   ::Float64,
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
        return iTgbmpb(iTgbmpb(0.0, δt, 0.0, λv),
                       iTgbmpb(0.0, δt, 0.0, λv),
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
      return iTgbmpb(sim_gbmpb(t, λt1, α, σλ, δt, srδt), 
                     sim_gbmpb(t, λt1, α, σλ, δt, srδt), 
              bt, δt, δt, λv)
    end

    λt = λt1
  end
end






"""
    sim_gbmpb(t   ::Float64,
              λt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              nsp ::Int64,
              nlim::Int64)

Simulate `iTgbmpb` according to a pure-birth geometric Brownian motion.
"""
function sim_gbmpb(t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   nsp ::Int64,
                   nlim::Int64)

  if nsp < nlim

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
          nsp += 1
          return iTgbmpb(iTgbmpb(0.0, δt, 0.0, λv),
                         iTgbmpb(0.0, δt, 0.0, λv),
                         bt, δt, t, λv), nsp
        end

        return iTgbmpb(bt, δt, t, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nsp += 1
        td1, nsp = sim_gbmpb(t, λt1, α, σλ, δt, srδt, nsp, nlim)
        td2, nsp = sim_gbmpb(t, λt1, α, σλ, δt, srδt, nsp, nlim)
        
        return iTgbmpb(td1, td2, bt, δt, δt, λv), nsp
      end

      λt = λt1
    end
  else
    return iTgbmpb(), nsp
  end

end





"""
    divev(λ::Float64, δt::Float64)

Return true if diversification event.
"""
divev(λ::Float64, δt::Float64) = @fastmath rand() < λ*δt 

