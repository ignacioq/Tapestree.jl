#=

constant occurrence birth-death simulation

Jérémy Andréoletti
Adapted from constant birth-death simulation by Ignacio Quintero Mächler

v(°-°v)

Created 07 10 2021
=#




"""
    cobd_wait(n::Float64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

Sample a waiting time for constant occurrence birth-death when `n` species 
are alive with speciation rate `λ`, extinction rate `μ`, fossil sampling rates 
`ψ` (fossils included in the tree) and `ω`(fossil occurrences).
"""
cobd_wait(n::Float64, λ::Float64, μ::Float64, 
          ψ::Float64, ω::Float64) = rexp(n*(λ + μ + ψ + ω))




"""
    cobd_wait(λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

Sample a per-lineage waiting time for constant occurrence birth-death 
with speciation rate `λ`, extinction rate `μ`, fossil sampling rates 
`ψ` (fossils included in the tree) and `ω`(fossil occurrences).
"""
cobd_wait(λ::Float64, μ::Float64, ψ::Float64, ω::Float64) = rexp(λ + μ + ψ + ω)




"""
    λevent(λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

Return `true` if speciation event.
"""
λevent(λ::Float64, μ::Float64, ψ::Float64, ω::Float64) = (λ/(λ+μ+ψ+ω)) > rand()




"""
    λevent(μ::Float64, ψ::Float64, ω::Float64)

Return `true` if extinction event, conditionned on "not a speciation event".
"""
μevent(μ::Float64, ψ::Float64, ω::Float64) = (μ/(μ + ψ + ω)) > rand()




"""
    ψevent(ψ::Float64, ω::Float64)

Return `true` if the sampled fossil is included in the tree.
"""
ψevent(ψ::Float64, ω::Float64) = (ψ/(ψ + ω)) > rand()




"""
    sim_cobd(t::Float64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

Simulate a constant occurrence birth-death `iTree` of height `t` with speciation 
rate `λ`, extinction rate `μ`, fossil sampling rates `ψ` (fossils included in 
the tree) and `ω`(fossil occurrences).
"""
function sim_cobd(t::Float64, 
                  λ::Float64, 
                  μ::Float64, 
                  ψ::Float64,
                  ω::Float64)

  tw = cobd_wait(λ, μ, ψ, ω)
  ωtimes = Float64[]

  if tw > t
    return sTfbd(t), ωtimes
  end

  if λevent(λ, μ, ψ, ω)
    # speciation
    d1, ωtimes1 = sim_cobd(t - tw, λ, μ, ψ, ω)
    d2, ωtimes2 = sim_cobd(t - tw, λ, μ, ψ, ω)
    append!(ωtimes1, ωtimes2)
    return sTfbd(d1, d2, tw), ωtimes1
  elseif μevent(μ, ψ, ω)
    # extinction
    return sTfbd(tw, true), ωtimes
  
  elseif ψevent(ψ, ω)
    # fossil sampling (included in the tree)
    d, ωtimes = sim_cobd(t - tw, λ, μ, ψ, ω)
    return sTfbd(d, tw, false, true, false), ωtimes
  
  else
    # fossil occurrence sampling (not included in the tree)
    d, ωtimes = sim_cobd(t - tw, λ, μ, ψ, ω)
    push!(ωtimes, t - tw)
    adde!(d, tw)
    return d, ωtimes
  end
end




"""
    sim_cobd(t::Float64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64, 
             na::Int64, nfos::Int64)

Simulate a constant occurrence birth-death `iTree` of height `t` with speciation 
rate `λ`, extinction rate `μ`, fossil sampling rates `ψ` (fossils included in the 
tree) and `ω`(fossil occurrences). `na` initializes the number of alived tips.
"""
function sim_cobd(t::Float64, 
                  λ::Float64, 
                  μ::Float64, 
                  ψ::Float64,
                  ω::Float64,
                  na::Int64,
                  nfos::Int64)

  tw = cobd_wait(λ, μ, ψ, ω)
  ωtimes = Float64[]

  if tw > t
    na += 1
    return sTfbd(t), ωtimes, na, nfos
  end

  if λevent(λ, μ, ψ, ω)
    # speciation
    d1, ωtimes1, na, nfos = sim_cobd(t - tw, λ, μ, ψ, ω, na, nfos)
    d2, ωtimes2, na, nfos = sim_cobd(t - tw, λ, μ, ψ, ω, na, nfos)
    append!(ωtimes1, ωtimes2)
    return sTfbd(d1, d2, tw), ωtimes1, na, nfos

  elseif μevent(μ, ψ, ω)
    # extinction
    return sTfbd(tw, true), ωtimes, na, nfos

  elseif ψevent(ψ, ω)
    # fossil sampling (included in the tree)
    nfos += 1
    @show nfos
    d, ωtimes, na, nfos = sim_cobd(t - tw, λ, μ, ψ, ω, na, nfos)
    return sTfbd(d, tw, false, true, false), ωtimes, na, nfos
  
  else
    # fossil occurrence sampling (not included in the tree)
    d, ωtimes, na, nfos = sim_cobd(t - tw, λ, μ, ψ, ω, na, nfos)
    push!(ωtimes, t - tw)
    adde!(d, tw)
    return d, ωtimes, na, nfos
  end
end




"""
   sim_cobd_b(n::Int64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

Simulate constant occurrence birth-death in backward time.
"""
function sim_cobd_b(n::Int64, 
                    λ::Float64, 
                    μ::Float64, 
                    ψ::Float64,
                    ω::Float64)

  nF = Float64(n)
  nI = n
  ωtimes = Float64[]

  # disjoint trees vector 
  tv = sTfbd[]
  for i in Base.OneTo(nI)
    push!(tv, sTfbd(0.0))
  end

  th = 0.0

  # start simulation
  while true
    w = cobd_wait(nF, λ, μ, ψ, ω)

    # track backward time
    th += w

    for t in tv
      adde!(t, w)
    end

    # if speciation
    if λevent(λ, μ, ψ, ω)
      if isone(nI)
        return tv[nI], ωtimes
      else
        j, k = samp2(Base.OneTo(nI))
        tv[j] = sTfbd(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    
    # if extinction
    elseif μevent(μ, ψ, ω)
      nI += 1
      nF += 1.0
      push!(tv, sTfbd(0.0, true))
    
    # if fossil sampling (included in the tree)
    elseif ψevent(ψ, ω)
      j = rand(Base.OneTo(nI))
      tv[j] = sTfbd(tv[j], 0.0, false, true, false)
    
    # if fossil occurrence sampling (not included in the tree)
    else
      push!(ωtimes, th)
    end
  end
end




"""
    sim_cobd_b(λ::Float64, 
               μ::Float64, 
               ψ::Float64,
               ω::Float64, 
               mxth::Float64, 
               maxn::Int64)

Simulate constant occurrence birth-death in backward time conditioned on 
1 survival and not having a greater tree height than `mxth`.
"""
function sim_cobd_b(λ::Float64, 
                    μ::Float64, 
                    ψ::Float64,
                    ω::Float64, 
                    mxth::Float64, 
                    maxn::Int64)

  nF = 1.0
  nI = 1
  ωtimes = Float64[]

  # disjoint trees vector 
  tv = [sTfbd(0.0, false)]

  th = 0.0

  # start simulation
  while true
    w   = cobd_wait(nF, λ, μ, ψ, ω)

    # track backward time
    th += w

    if nI > maxn
      return tv[nI], ωtimes, (mxth + 0.1)
    end
    
    if th > mxth 
     return tv[nI], ωtimes, th
    end

    for t in tv
      adde!(t, w)
    end

    # if speciation
    if λevent(λ, μ, ψ, ω)
      if isone(nI)
        return tv[nI], ωtimes, th
      else
        j, k = samp2(Base.OneTo(nI))
        tv[j] = sTfbd(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    
    # if extinction
    elseif μevent(μ, ψ, ω)
      nI += 1
      nF += 1.0
      push!(tv, sTfbd(0.0, true))
    
    # if fossil sampling (included in the tree)
    elseif ψevent(ψ, ω)
      j = rand(Base.OneTo(nI))
      tv[j] = sTfbd(tv[j], 0.0, false, true, false)
    
    # if fossil occurrence sampling (not included in the tree)
    else
      push!(ωtimes, th)
    end
  end
end





