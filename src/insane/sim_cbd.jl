#=

constant birth-death simulation

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
    rexp(r::Float64)

Generate an exponential sample with rate `r`.
"""
rexp(r::Float64) = @fastmath randexp()/r




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
function sim_cbd(t::Float64, 
                 λ::Float64, 
                 μ::Float64, 
                 ::Type{T}) where {T <: iTree}

  tw = cbd_wait(λ, μ)

  if tw > t
    return T(t)
  end

  if λorμ(λ, μ)
    return T(sim_cbd(t - tw, λ, μ, T), sim_cbd(t - tw, λ, μ, T), tw)
  else
    return T(tw, true)
  end
end




"""
   sim_cbd_b(n::Int64, λ::Float64, μ::Float64)

Simulate constant birth-death in backward time.
"""
function sim_cbd_b(n::Int64, 
                   λ::Float64, 
                   μ::Float64, 
                   ::Type{T}) where {T <: iTree}

  nF = Float64(n)
  nI = n

  # disjoint trees vector 
  tv = T[]
  for i in Base.OneTo(nI)
    push!(tv, iT(0.0))
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
        tv[j] = T(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    # if extinction
    else
      nI += 1
      nF += 1.0
      push!(tv, T(0.0, true))
    end
  end
end




"""
    sim_cbd_b(λ::Float64, μ::Float64, mxth::Float64)

Simulate constant birth-death in backward time conditioned on 1 survival 
and not having a greater tree height than `mxth`.
"""
function sim_cbd_b(λ::Float64, 
                   μ::Float64, 
                   mxth::Float64, 
                   maxn::Int64, 
                   ::Type{T}) where {T <: iTree}

  nF = 1.0
  nI = 1

  # disjoint trees vector 
  tv = [T(0.0, true)]

  th = 0.0

  # start simulation
  while true
    w   = cbd_wait(nF, λ, μ)

    # track backward time
    th += w

    if nI > maxn
      return tv[nI], (mxth + 0.1)
    end
    
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
        tv[j] = T(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    # if extinction
    else
      nI += 1
      nF += 1.0
      push!(tv, T(0.0, true))
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






