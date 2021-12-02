#=

constant pure-birth simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    cpb_wait(λ::Float64)

Sample a per-lineage waiting time for pure-birth species 
with speciation rate `λ`.
"""
cpb_wait(λ::Float64) = rexp(λ)




"""
  sim_cpb(t::Float64, λ::Float64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`.
"""
function sim_cpb(t::Float64, λ::Float64)

  tw = cpb_wait(λ)

  if tw > t
    return sTpb(t)
  end

  sTpb(sim_cpb(t - tw, λ), sim_cpb(t - tw, λ), tw, false)
end




"""
  sim_cpb(t::Float64, λ::Float64, na::Int64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`.
"""
function sim_cpb(t::Float64, λ::Float64, na::Int64)

  tw = cpb_wait(λ)

  if tw > t
    na += 1
    return sTpb(t), na
  end
  d1, na = sim_cpb(t - tw, λ, na)
  d2, na = sim_cpb(t - tw, λ, na)

  sTpb(d1, d2, tw, false), na
end





