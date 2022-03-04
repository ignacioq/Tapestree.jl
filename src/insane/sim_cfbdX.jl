#=

constant birth-death simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    sim_cfbd(t ::Float64,
             λ ::Float64,
             μ ::Float64,
             ψ ::Float64,
             x0::Float64,
             σx::Float64)

Simulate a constant fossilized birth-death `iTree` of height `t` 
with speciation rate `λ`, extinction rate `μ` and fossil sampling rate `ψ`.
"""
function sim_cfbd(t ::Float64,
                  λ ::Float64,
                  μ ::Float64,
                  ψ ::Float64,
                  x0::Float64,
                  σx::Float64)

  tw = cfbd_wait(λ, μ, ψ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    return sTfbdX(t, false, false, false, x0, x1)
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  if λevent(λ, μ, ψ)
    # speciation
    return sTfbdX(sim_cfbd(t - tw, λ, μ, ψ, x1, σx), 
                  sim_cfbd(t - tw, λ, μ, ψ, x1, σx), 
                  tw, false, false, false, x0, x1)
  elseif μevent(μ, ψ)
    # extinction
    return sTfbdX(tw, true, false, false, x0, x1)
  else
    # fossil sampling
    return sTfbdX(sim_cfbd(t - tw, λ, μ, ψ, x1, σx), 
                  tw, false, true, false, x0, x1)
  end
end





"""
    sim_cfbd(t ::Float64,
             λ ::Float64,
             μ ::Float64,
             ψ ::Float64,
             x0::Float64,
             σx::Float64,
             na::Int64)

Simulate a constant  fossilized birth-death `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cfbd(t ::Float64,
                  λ ::Float64,
                  μ ::Float64,
                  ψ ::Float64,
                  x0::Float64,
                  σx::Float64,
                  na::Int64)

  tw = cfbd_wait(λ, μ, ψ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    na += 1
    return sTfbdX(t, false, false, false, x0, x1), na, nf
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  # speciation
  if λevent(λ, μ, ψ)
    d1, na, nf = sim_cfbd(t - tw, λ, μ, ψ, x1, σx)
    d2, na, nf = sim_cfbd(t - tw, λ, μ, ψ, x1, σx)

    return sTfbdX(d1, d2, tw, false, false, false, x0, x1), na, nf
  # extinction
  elseif μevent(μ, ψ)
    return sTfbdX(tw, true, false, false, x0, x1), na, nf
  # fossil sampling
  else
    nf += 1
    d1, na, nf = sim_cfbd(t - tw, λ, μ, ψ, x1, σx)
    return sTfbdX(d1, tw, false, true, false, x0, x1), na, nf
  end
end




