#=

birth-death simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    cbd_wait(n::Float64, λ::Float64, μ::Float64)

Sample a waiting time for constant birth-death when `n` species 
are alive with speciation rate `λ` and extinction rate `μ`.
"""
cbd_wait(n::Float64, λ::Float64, μ::Float64) = rexp(n*(λ + μ))




"""
    cbd_wait(λ::Float64, μ::Float64)

Sample a per-lineage waiting time for constant birth-death species 
with speciation rate `λ` and extinction rate `μ`.
"""
cbd_wait(λ::Float64, μ::Float64) = rexp(λ + μ)




"""
    cpb_wait(λ::Float64)

Sample a per-lineage waiting time for pure-birth species 
with speciation rate `λ`.
"""
cpb_wait(λ::Float64) = rexp(λ)




"""
    rexp(r::Float64)

Generate an exponential sample with rate `r`.
"""
rexp(r::Float64) = @fastmath randexp()/r




"""
  sim_cpb(t::Float64, λ::Float64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`.
"""
function sim_cpb(t::Float64, λ::Float64)

  tw = cpb_wait(λ)

  if tw > t
    return iTree(t)
  end

  iTree(sim_cpb(t - tw, λ), sim_cpb(t - tw, λ), tw)
end




"""
    λorμ(λ::Float64, μ::Float64)

Return `true` if speciation event
"""
λorμ(λ::Float64, μ::Float64) = λ/(λ + μ) > rand() ? true : false




"""
    sim_cbd(t::Float64, λ::Float64, μ::Float64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ`.
"""
function sim_cbd(t::Float64, λ::Float64, μ::Float64)

  tw = cbd_wait(λ, μ)

  if tw > t
    return iTree(t)
  end

  if λorμ(λ, μ)
    return iTree(sim_cbd(t - tw, λ, μ), sim_cbd(t - tw, λ, μ), tw)
  else
    return iTree(tw, true)
  end
end




"""
   sim_cbd_b(n::Int64, λ::Float64, μ::Float64)

Simulate constant birth-death in backward time.
"""
function sim_cbd_b(n::Int64, λ::Float64, μ::Float64)

  nF = Float64(n)
  nI = n

  # disjoint trees vector 
  tv = iTree[]
  for i in Base.OneTo(nI)
    push!(tv, iTree(0.0))
  end

  # start simulation
  while true
    w = cbd_wait(nF, λ, μ)

    for t in tv
      addpe!(t, w)
    end

    # if speciation
    if λorμ(λ, μ)
      if isone(nI)
        return tv[nI] 
      else
        j, k = samp2(Base.OneTo(nI))
        tv[j] = iTree(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    # if extinction
    else
      nI += 1
      nF += 1.0
      push!(tv, iTree(0.0, true))
    end
  end
end




"""
    sim_cbd_b(λ::Float64, μ::Float64, mxth::Float64)

Simulate constant birth-death in backward time conditioned on extinction 
and not having a greater tree height than `mxth`.
"""
function sim_cbd_b(λ::Float64, μ::Float64, mxth::Float64)

  nF = 1.0
  nI = 1

  # disjoint trees vector 
  tv = [iTree(0.0, true)]

  th = 0.0

  # start simulation
  while true
    w   = cbd_wait(nF, λ, μ)

    # track backward time
    th += w

    if th > mxth 
     return tv[nI], th
    end

    for t in tv
      addpe!(t, w)
    end

    # if speciation
    if λorμ(λ, μ)
      if isone(nI)
        return tv[nI], th
      else
        j, k = samp2(Base.OneTo(nI))
        tv[j] = iTree(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    # if extinction
    else
      nI += 1
      nF += 1.0
      push!(tv, iTree(0.0, true))
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






