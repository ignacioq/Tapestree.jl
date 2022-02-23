#=

constant birth-death simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    sim_cbd(t::Float64, λ::Float64, μ::Float64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ`.
"""
function sim_cbd(t ::Float64, 
                 λ ::Float64, 
                 μ ::Float64,
                 x0::Float64, 
                 σx::Float64)

  tw = cbd_wait(λ, μ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    return sTbdX(t, false, false, x0, x1)
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  if λorμ(λ, μ)
    return sTbdX(sim_cbd(t - tw, λ, μ, x1, σx), 
                 sim_cbd(t - tw, λ, μ, x1, σx), 
                 tw, false, false, x0, x1)
  else
    return sTbdX(tw, true, false, x0, x1)
  end
end





"""
    sim_cbd(t ::Float64, 
            λ ::Float64, 
            μ ::Float64,
            x0::Float64, 
            σx::Float64
            na::Int64)

Simulate a constant birth-death `iTree` with traits of height `t` 
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cbd(t ::Float64, 
                 λ ::Float64, 
                 μ ::Float64,
                 x0::Float64, 
                 σx::Float64,
                 na::Int64)

  tw = cbd_wait(λ, μ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    na += 1
    return sTbdX(t, false, false, x0, x1), na
  end

  x1 = rnorm(x0, sqrt(tw) * σx)


  if λorμ(λ, μ)
    d1, na = sim_cbd(t - tw, λ, μ, x1, σx, na)
    d2, na = sim_cbd(t - tw, λ, μ, x1, σx, na)

    return sTbdX(d1, d2, tw, false, false, x0, x1), na
  else
    return sTbdX(tw, true, false, x0, x1), na
  end
end




