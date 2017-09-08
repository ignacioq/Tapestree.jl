#=

Utilities for Continuous Data Augmentation in 
Biogeographic Histories.

Ignacio Quintero Mächler

t(-_-t)

May 01 2017

=#




"""
  branch sampling for multiple areas,
  resampling if extinction is present
"""
function br_samp(ssii ::Array{Int64,1}, 
                 ssff ::Array{Int64,1},
                 λc   ::Array{Float64,1},
                 t    ::Float64,
                 narea::Int64)
  #time history
  t_hist = mult_rejsam(ssii, ssff, λc, t, narea)
  
  while ifext(t_hist, ssii, narea)
      t_hist = mult_rejsam(ssii, ssff, λc, t, narea)
  end

  t_hist
end





"""
  check if extinct
"""
function ifext(t_hist::Array{Array{Float64,1},1},
               ssii  ::Array{Int64,1}, 
               narea ::Int64)

  cs_hist = Float64[]
  hist_l  = Int64[]

  for i in eachindex(t_hist)
    append!(cs_hist, cumsum(t_hist[i]))
    append!(hist_l, fill(i,length(t_hist[i])))
  end

  # organize order of events
  sp  = sortperm(cs_hist)[1:(end - narea)]
  lhs = length(sp) + narea + 1

  # reconstruct state history
  s_hist = zeros(Int64, lhs,narea)
  for i in eachindex(ssii)
    s_hist[:,i] = ssii[i]
  end

  for i = eachindex(sp)
    setindex!(s_hist, 1 - s_hist[i,hist_l[sp[i]]], (i+1):lhs, hist_l[sp[i]])
  end

  # check if extinct
  for j = eachindex(sp)
    ss = 0
    for i = 1:narea
      ss += s_hist[j,i]
    end
    if ss == 0 
      return true
    end
  end

  return false
end




"""
  multistate branch sampling
"""
function mult_rejsam(ssii ::Array{Int64,1}, 
                     ssff ::Array{Int64,1},
                     λc   ::Array{Float64,1},
                     t    ::Float64,
                     narea::Int64)

  all_times = Array{Float64,1}[]

  for i in Base.OneTo(narea)
    push!(all_times, rejsam(ssii[i], ssff[i], λc[1], λc[2], t))
  end

  return all_times
end




"""
  rejection sampling for each branch
  condition on start and end point
"""
function rejsam(si::Int64, sf::Int64, λ1::Float64, λ0::Float64, t::Float64)
  
  sam::Tuple{Array{Float64,1},Int64} = brprop(si, λ1, λ0, t)
  
  while sam[2] != sf 
    sam = brprop(si, λ1, λ0, t)
  end

  return sam[1]
end




"""
  propose events for a branch
"""
function brprop(si::Int64, λ1::Float64, λ0::Float64, t::Float64)

  c_st   ::Int64            = si
  c_time ::Float64          = 0.0
  times  ::Array{Float64,1} = zeros(0)
  endtime::Float64          = t

  re::Float64 = c_st == 0 ? rexp(λ1) : rexp(λ0)

  c_time += re 

  while c_time < t
    push!(times, re)
    endtime  = t - c_time 
    c_st     = 1 - c_st
    re       = c_st == 0 ? rexp(λ1) : rexp(λ0)
    c_time  += re
  end

  push!(times, endtime)

  return times, c_st
end
