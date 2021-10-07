#=

constant fossilized birth-death simulation

Jérémy Andréoletti
Adapted from constant birth-death simulation by Ignacio Quintero Mächler

v(°_°v)

Created 07 10 2021
=#




"""
    cfbd_wait(n::Float64, λ::Float64, μ::Float64, ψ::Float64)

Sample a waiting time for constant fossilized birth-death when `n` species 
are alive with speciation rate `λ` and extinction rate `μ`.
"""
cfbd_wait(n::Float64, λ::Float64, μ::Float64, ψ::Float64) = rexp(n*(λ + μ + ψ))




"""
    cfbd_wait(λ::Float64, μ::Float64, ψ::Float64)

Sample a per-lineage waiting time for constant fossilized birth-death 
with speciation rate `λ` and extinction rate `μ`.
"""
cfbd_wait(λ::Float64, μ::Float64, ψ::Float64) = rexp(λ + μ + ψ)





"""
    rexp(r::Float64)

Generate an exponential sample with rate `r`.
"""
rexp(r::Float64) = @fastmath randexp()/r




"""
    λevent(λ::Float64, μ::Float64, ψ::Float64)

Return `true` if speciation event
"""
λevent(λ::Float64, μ::Float64, ψ::Float64) = (λ/(λ + μ + ψ)) > rand()




"""
    λevent(μ::Float64, ψ::Float64)

Return `true` if extinction event, conditionned on "not a speciation event"
"""
μevent(μ::Float64, ψ::Float64) = (μ/(μ + ψ)) > rand()




"""
    sim_cfbd(t::Float64, λ::Float64, μ::Float64, ψ::Float64)

Simulate a constant fossilized birth-death `iTree` of height `t` with speciation 
rate `λ` and extinction rate `μ`.
"""
function sim_cfbd(t::Float64, 
                  λ::Float64, 
                  μ::Float64, 
                  ψ::Float64)

  tw = cfbd_wait(λ, μ, ψ)

  if tw > t
    return sTfbd(t)
  end

  if λevent(λ, μ, ψ)
    # speciation
    return sTfbd(sim_cfbd(t - tw, λ, μ, ψ), sim_cfbd(t - tw, λ, μ, ψ), tw)
  elseif μevent(μ, ψ)
    # extinction
    return sTfbd(tw, true)
  else
    # fossil sampling
    return sTfbd(sim_cfbd(t - tw, λ, μ, ψ), tw, false, true, false)
  end
end




"""
   sim_cfbd_b(n::Int64, λ::Float64, μ::Float64, ψ::Float64)

Simulate constant fossilized birth-death in backward time.
"""
function sim_cfbd_b(n::Int64, 
                    λ::Float64, 
                    μ::Float64, 
                    ψ::Float64)

  nF = Float64(n)
  nI = n

  # disjoint trees vector 
  tv = sTfbd[]
  for i in Base.OneTo(nI)
    push!(tv, sTfbd(0.0))
  end

  # start simulation
  while true
    w = cfbd_wait(nF, λ, μ, ψ)

    for t in tv
      adde!(t, w)
    end

    # if speciation
    if λevent(λ, μ, ψ)
      if isone(nI)
        return tv[nI]
      else
        j, k = samp2(Base.OneTo(nI))
        tv[j] = sTfbd(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    # if extinction
    else
      nI += 1
      nF += 1.0
      push!(tv, sTfbd(0.0, true))
    end
  end
end




"""
    sim_cfbd_b(λ::Float64, 
               μ::Float64, 
               ψ::Float64, 
               mxth::Float64, 
               maxn::Int64)

Simulate constant fossilized birth-death in backward time conditioned on 
1 survival and not having a greater tree height than `mxth`.
"""
function sim_cfbd_b(λ::Float64, 
                    μ::Float64, 
                    ψ::Float64, 
                    mxth::Float64, 
                    maxn::Int64)

  nF = 1.0
  nI = 1

  # disjoint trees vector 
  tv = [sTfbd(0.0, false)]

  th = 0.0

  # start simulation
  while true
    w   = cfbd_wait(nF, λ, μ, ψ)

    # track backward time
    th += w

    if nI > maxn
      return tv[nI], (mxth + 0.1)
    end
    
    if th > mxth 
     return tv[nI], th
    end

    for t in tv
      adde!(t, w)
    end

    # if speciation
    if λevent(λ, μ, ψ)
      if isone(nI)
        return tv[nI], th
      else
        j, k = samp2(Base.OneTo(nI))
        tv[j] = sTfbd(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    # if extinction
    else
      nI += 1
      nF += 1.0
      push!(tv, sTfbd(0.0, true))
    end
  end
end




"""
    samp2(o::Base.OneTo{Int64})

Sample `2` without replacement from `o`.
"""
function samp2(o::Base.OneTo{Int64})
  j = rand(o)
  k = rand(o)
  while k == j
    k = rand(o)
  end
  return j, k
end






