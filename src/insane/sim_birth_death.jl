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
    rexp(r::Float64)

Generate an exponential sample with rate `r`.
"""
rexp(r::Float64) = randexp()/r




"""
    λorμ(λ::Float64, μ::Float64)

Return `true` if speciation event
"""
λorμ(λ::Float64, μ::Float64) = λ/(λ + μ) > rand() ? true : false




"""
    addpe!(tree::itree, bs::Array{Bool,1}, pe::Float64)

Add `pe` to pendant edge of node given by binary search vector `bs`.
Binary search vector follows (left -> false, right -> true).
"""
function addpe!(tree::itree, bs::Array{Bool,1}, pe::Float64)
  for i in bs
    tree = i ? tree.d2 : tree.d1
  end
  tree.pe += pe
end




"""
    newλ!(tree::itree, bs::Array{Bool,1})

Add new daughters to node given by binary search vector `bs`.
Binary search vector follows (left -> false, right -> true).
"""
function newλ!(tree::itree, bs::Array{Bool,1})
  for i in bs
    tree = i ? tree.d2 : tree.d1
  end
  tree.d1 = itree()
  tree.d2 = itree()
end




"""
    sim_cbd(λ::Float64, μ::Float64, tf::Float64)

Simulate complete birth-death trees in `itree` format until time `tf`.
"""
function sim_cbd(λ::Float64, μ::Float64, tf::Float64)

  tc = 0.0   # current time
  nF = 1.0   # number of alive lineages (Float)
  nI = 1     # number of alive lineages (Int)

  nλ = 0     # number of speciation events
  nμ = 0     # number of extinction events

  tree  = itree()

  wa = Array{Bool,1}[]
  push!(wa, Bool[])

  while true

    # waiting time
    tw = cbd_wait(nF, λ, μ)

    tc += tw

    if tc > tf
      for bs in wa
        addpe!(tree, bs, tf - tc + tw)
      end
      break
    else
      for bs in wa
        addpe!(tree, bs, tw)
      end
    end

    # if speciation
    if λorμ(λ, μ)
      nλ += 1

      # which speciates
      wλ = rand(1:nI)
      bs = wa[wλ]

      # add new nodes
      newλ!(tree, bs)

      # add new lineage
      nI += 1
      nF += 1.0

      # change which alive
      push!(bs, false)
      push!(wa, push!(copy(bs[1:(end-1)]),true))

    # if extinction
    else
      nμ += 1

      # which goes extinct
      wλ  = rand(1:nI)
      deleteat!(wa, wλ)

      nI -= 1
      nF -= 1.0

      if iszero(nI) 
        @warn "tree has no survivors at time $tf"
        break
      end
    end
  end

  return tree, nλ, nμ
end




"""
    sim_cbd(λ::Float64, μ::Float64, tf::Float64)

Simulate complete birth-death trees in graph format until time `tf`.
"""
function sim_cbd(λ::Float64, μ::Float64, tf::Float64, graph::Bool)

  tc = 0.0   # current time
  nF = 1.0   # number of alive lineages (Float)
  nI = 1     # number of alive lineages (Int)
  pn = [1]   # parent nodes
  dn = [2]   # daughter nodes
  el = [0.0] # branch lengths
  wa = [1]   # which alive
  nc = 2     # node count
  nλ = 0     # number of speciation events
  nμ = 0     # number of extinction events

  while true

    # waiting time
    tw = cbd_wait(nF, λ, μ)

    tc += tw

    if tc > tf 
      for i in wa
        el[i] += (tf - tc + tw)
      end
      break
    else
      for i in wa
        el[i] += tw
      end
    end

    # if speciation
    if λorμ(λ, μ)

      nλ += 1
      # which speciates
      wλ  = rand(1:nI)
      wwλ = wa[wλ]

      # add new nodes
      push!(pn, dn[wwλ], dn[wwλ])
      nc += 1
      push!(dn, nc)
      nc += 1
      push!(dn, nc)

      # add new edges
      push!(el, 0.0, 0.0)

      # add new lineage
      nI += 1
      nF += 1.0

      # change which alive
      li = lastindex(dn)
      push!(wa, li-1, li)
      deleteat!(wa, wλ)

    # if extinction
    else
      nμ += 1

      # which goes extinct
      wλ  = rand(1:nI)
      wwλ = wa[wλ]

      deleteat!(wa, wλ)

      nI -= 1
      nF -= 1.0

      if iszero(nI) 
        @warn "tree has no survivors at time $tf"
        break
      end
    end
  end

  return pn, dn, el, nλ, nμ
end







