
"""

Utilities for Data Augmented Competition model

Ignacio Quintero Mächler

t(-_-t)

May 01 2017

"""




"""
  log-density of exponential 
"""
logdexp(x::Float64, λ::Float64) = @fastmath log(λ) - λ * x




"""
  log-density of normal
"""
logdnorm(x::Float64, μ::Float64, σ²::Float64) = 
 @fastmath -(0.5*log(2.0π) + 0.5*log(σ²) + abs2(x - μ)/(2.0 * σ²))




"""
  log-density of half-cauchy
"""
logdhcau(x::Float64, scl::Float64) = 
  @fastmath log(2 * scl/(π *(x * x + scl * scl)))




"""
  log-density of half-cauchy with scale 1
"""
logdhcau1(x::Float64) = 
  @fastmath log(2/(π * (x * x + 1)))




"""
  random exponential generator
"""
rexp(λ::Float64) = @fastmath log(rand()) * -(1/λ)





"""
  function for sampling a coin flip 
  with non equal probilities
"""
coinsamp(p0::Float64) = rand() < p0 ? 0 : 1





"""
  normalize probabilities to 1
"""
normlize(pt1::Float64, pt2::Float64) = pt1/(pt1 + pt2)





"""
  Markov chain probabilities 
  through fast analytical solution
"""
function Ptrfast(λ1::Float64, λ0::Float64, t::Float64)
  
  @fastmath begin

    sumλ ::Float64 = λ1 + λ0
    ex   ::Float64 = exp(-sumλ*t)
    sd1  ::Float64 = 1/sumλ
    λ1ex ::Float64 = λ1*ex
    λ0ex ::Float64 = λ0*ex

    ((sd1*(λ0 + λ1ex), sd1*(λ1 - λ1ex)), 
     (sd1*(λ0 - λ0ex), sd1*(λ1 + λ0ex)))
  end
end






"""
  Markov chain probabilities 
  through fast analytical solution
  conditional on starting value
"""
function Ptrfast_start(λ1::Float64, λ0::Float64, t::Float64, state::Int64)
  
  @fastmath begin

    sumλ ::Float64 = λ1 + λ0
    ex   ::Float64 = exp(-sumλ*t)
    sd1  ::Float64 = 1/sumλ
    λ1ex ::Float64 = λ1*ex
    λ0ex ::Float64 = λ0*ex

    if state == 0
      sd1*(λ0 + λ1ex), sd1*(λ1 - λ1ex) 
    else
      sd1*(λ0 - λ0ex), sd1*(λ1 + λ0ex)
    end
  end
end






"""
  Markov chain probabilities 
  through fast analytical solution
  conditional on starting value
"""
function Ptrfast_end(λ1::Float64, λ0::Float64, t::Float64, state::Int64)
  
  @fastmath begin

    sumλ ::Float64 = λ1 + λ0
    ex   ::Float64 = exp(-sumλ*t)
    sd1  ::Float64 = 1/sumλ
    λ1ex ::Float64 = λ1*ex
    λ0ex ::Float64 = λ0*ex

    if state == 0
      (sd1*(λ0 + λ1ex), sd1*(λ0 - λ0ex))
    else
      (sd1*(λ1 - λ1ex), sd1*(λ1 + λ0ex))
    end
  end
end







"""
  make ragged array with index for each edge in Yc
"""
function make_edgeind(childs::Array{Int64,1}, B::Array{Float64,2})

  bridx = Array{Int64,1}[]
  for b in childs
    push!(bridx, find(b .== B))
  end

  bridx
end







"""
  make ragged array of the cumulative delta times for each branch
"""
function make_edgeδt(bridx::Array{Array{Int64,1},1}, 
                     δt   ::Array{Float64,1}, 
                     m    ::Int64)
  
  brδt = Array{Float64,1}[]
  
  for j in 1:(length(bridx)-1)
    bi = copy(bridx[j] .- 1)
    for i in eachindex(bi)
      bi[i] = rowind(bi[i], m)
    end
    push!(brδt, cumsum(δt[bi]))
  end
  
  brδt
end









"""
  return index in vector "x" corresponding to a value 
  that is closest but smaller than "val" in sorted arrays 
  using a sort of uniroot algorithm
"""
function idxlessthan(x::Array{Float64,1}, val::Float64) 
  
  @inbounds begin

    a  ::Int64 = 1
    b  ::Int64 = endof(x)
  
    if x[b] < val
      return b
    end

    mid::Int64 = div(b,2)  

    while b-a > 1
      val < x[mid] ? b = mid : a = mid
      mid = div(b + a, 2)
    end

  end

  return a
end 








"""
  return index for closest value in sorted arrays 
  using a sort of uniroot algorithm
  FLOATS
"""
function indmindif_sorted(x::Array{Float64,1}, val::Float64) 
  a::Int64   = 1
  b::Int64   = endof(x)
  mid::Int64 = div(b,2)  

  while b-a > 1
    val < x[mid] ? b = mid : a = mid
    mid = div(b + a, 2)
  end

  abs(x[a] - val) < abs(x[b] - val) ? a : b
end 






"""
return index for closest value in sorted arrays 
using a sort of uniroot algorithm
INTEGERS
"""
function indmindif_sorted(x::Array{Int64,1}, val::Int64) 
  a::Int64   = 1
  b::Int64   = endof(x)
  mid::Int64 = div(b,2)  

  while b-a > 1
    val < x[mid] ? b = mid : a = mid
    mid = div(b + a, 2)
  end

  abs(x[a] - val) < abs(x[b] - val) ? a : b
end 





"""
  make branch triads:
  first number is the parent branch
  second and third numbers the daughters
"""
function maketriads(edges::Array{Int64,2})

  # internal nodes
  ins::Array{Int64,1} = unique(edges[:,1])[1:(end-1)]
  lins = length(ins)

  trios = Array{Int64,1}[]

  # for all internal nodes
  for i = Base.OneTo(lins)
    ndi  = ins[i]
    daus = find(edges[:,1] .== ndi)
    unshift!(daus, find(edges[:,2] .== ndi)[1])
    push!(trios, daus)
  end

  trios
end





