#=

constant pure-birth simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#



"""
  sim_cpb(t::Float64, λ::Float64, x0::Float64, σx::Float64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cpb(t::Float64, λ::Float64, x0::Float64, σx::Float64)

  tw = cpb_wait(λ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    return sTpbX(t, false, x0, x1)
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  sTpbX(sim_cpb(t - tw, λ, x1, σx),
        sim_cpb(t - tw, λ, x1, σx),
        tw, false, x0, x1)
end




"""
    sim_cpb(t::Float64, λ::Float64, x0::Float64, σx::Float64, na::Int64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cpb(t::Float64, λ::Float64, x0::Float64, σx::Float64, na::Int64)

  tw = cpb_wait(λ)

  if tw > t
    na += 1
    x1 = rnorm(x0, sqrt(t) * σx)
    return sTpbX(t, false, x0, x1), na
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  d1, na = sim_cpb(t - tw, λ, x1, σx, na)
  d2, na = sim_cpb(t - tw, λ, x1, σx, na)

  sTpbX(d1, d2, tw, false, x0, x1), na
end





