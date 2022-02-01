#=

Survival nodes for conditioning

Ignacio Quintero Mächler

t(-_-t)

Created 16 11 2021
=#

"""
    surv_m(t::Float64, λ::Float64, μ::Float64, ntry::Int64, stem ::Bool)

Sample the total number of `m` trials until both simulations survive 
for constant birth-death.
"""
function surv_m(t::Float64, λ::Float64, μ::Float64, ntry::Int64, stem ::Bool)

  ntries = 0
  m      = 0.0

  # if stem conditioning
  if stem

    while true
      m      += 1.0
      ntries += 1
      t, s, n = sim_cbd_surv(t, λ, μ, false, 1)

      s && break
      ntries == ntry && break
    end

  # if crown conditioning
  else

    while true
      m      += 1.0
      ntries += 1

      t1, s1, n1 = sim_cbd_surv(t, λ, μ, false, 1)
      t2, s2, n2 = sim_cbd_surv(t, λ, μ, false, 1)

      s1 && s2 && break
      ntries == ntry && break
    end
  end

  return m
end



