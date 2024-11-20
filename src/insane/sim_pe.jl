#=

constant punkeek simulation

Ignacio Quintero Mächler

t(-_-t)

Created 20 11 2024
=#


"""
    sim_pe(t::Float64, λ::Float64, μ::Float64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ`.
"""
function sim_pe(t ::Float64,
                λ ::Float64,
                μ ::Float64,
                x0::Float64,
                σa::Float64,
                σk::Float64)

  tw = cbd_wait(λ, μ)


  if tw > t
    x1 = rnorm(x0, sqrt(t) * σa)
    return peT(t, false, x0, x1, false)
  end

  x1 = rnorm(x0, sqrt(tw) * σa)

  if λorμ(λ, μ)
    xk = rnorm(x1, σk)
    if rand() < 0.5
      return peT(sim_pe(t - tw, λ, μ, xk, σa, σk), 
                 sim_pe(t - tw, λ, μ, x1, σa, σk), 
                 tw, false, x0, x1, true, false)
    else
      return peT(sim_pe(t - tw, λ, μ, x1, σa, σk), 
                 sim_pe(t - tw, λ, μ, xk, σa, σk), 
                 tw, false, x0, x1, false, false)
    end
  else
    return peT(tw, true, x0, x1, false)
  end
end


















