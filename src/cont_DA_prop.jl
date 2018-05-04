#=

Utilities for Continuous Data Augmentation in 
Biogeographic Histories.

Ignacio Quintero Mächler

t(-_-t)

May 01 2017

=#




"""
    rejsam!(times::Array{Float64,1}, 
            si   ::Int64, 
            sf   ::Int64, 
            λ1   ::Float64, 
            λ0   ::Float64, 
            t    ::Float64)

Rejection sampling for an area conditioned on start and end point.
"""
function rejsam!(times::Array{Float64,1}, 
                 si   ::Int64, 
                 sf   ::Int64, 
                 λ1   ::Float64, 
                 λ0   ::Float64, 
                 t    ::Float64)
  
  samf = brprop!(times, si, λ1, λ0, t)
  
  while samf != sf 
    samf = brprop!(times, si, λ1, λ0, t)
  end

  return nothing
end




"""
    brprop!(times::Array{Float64,1}, 
            si   ::Int64, 
            λ1   ::Float64, 
            λ0   ::Float64, 
            t    ::Float64)

Two state DA proposal for one area.
"""
function brprop!(times::Array{Float64,1}, 
                 si   ::Int64, 
                 λ1   ::Float64, 
                 λ0   ::Float64, 
                 t    ::Float64)
  empty!(times)

  c_st   ::Int64   = si
  c_time ::Float64 = 0.0
  endtime::Float64 = t

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

  return c_st
end





