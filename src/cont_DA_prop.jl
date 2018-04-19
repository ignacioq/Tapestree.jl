#=

Utilities for Continuous Data Augmentation in 
Biogeographic Histories.

Ignacio Quintero Mächler

t(-_-t)

May 01 2017

=#



"""
br_samp!(evs  ::Array{Array{Float64,1},1},
         ssii  ::Array{Int64,1}, 
         ssff  ::Array{Int64,1},
         λ1    ::Float64,
         λ0    ::Float64,
         t     ::Float64,
         narea ::Int64)

Sample branch biogeographic histories according to an independent
model that does not allow extinction.
"""
function br_samp!(evs   ::Array{Array{Float64,1},1},
                  ssii  ::Array{Int64,1}, 
                  ssff  ::Array{Int64,1},
                  λ1    ::Float64,
                  λ0    ::Float64,
                  t     ::Float64,
                  narea ::Int64)
  #time history
  mult_rejsam!(evs, ssii, ssff, λ1, λ0, t, narea)

  while ifext(evs, ssii, narea, t)
    mult_rejsam!(evs, ssii, ssff, λ1, λ0, t, narea)
  end

  return nothing
end





"""
    ifext(t_hist::Array{Array{Float64,1},1},
          ssii  ::Array{Int64,1}, 
          narea ::Int64,
          t     ::Float64)

Return true if lineage goes extinct.
"""
function ifext(t_hist::Array{Array{Float64,1},1},
               ssii  ::Array{Int64,1}, 
               narea ::Int64,
               t     ::Float64)

  # initial occupancy time
  ioc  = findfirst(ssii)::Int64
  ioct = t_hist[ioc][1]::Float64

  ntries = 0
  while ioct < t

    if ioc == narea
      ioc = 1
    else 
      ioc += 1
    end

    tc = 0.0
    cs = ssii[ioc]
    for ts in t_hist[ioc]::Array{Float64,1}
      tc += ts
      if ioct < tc 
        if cs == 1
          ioct   = tc 
          ntries = 0
          break
        else
          ntries += 1
          if ntries > narea
            return true
          end
          break
        end
      end
      cs = 1 - cs
    end

  end

  return false
end





"""
  mult_rejsam!(evs   ::Array{Array{Float64,1},1},
               ssii  ::Array{Int64,1}, 
               ssff  ::Array{Int64,1},
               λ1    ::Float64,
               λ0    ::Float64,
               t     ::Float64,
               narea ::Int64)

  Multi-area branch rejection independent model sampling.
"""
function mult_rejsam!(evs  ::Array{Array{Float64,1},1},
                      ssii ::Array{Int64,1}, 
                      ssff ::Array{Int64,1},
                      λ1   ::Float64,
                      λ0   ::Float64,
                      t    ::Float64,
                      narea::Int64)

  for i = Base.OneTo(narea)
    rejsam!(evs[i], ssii[i], ssff[i], λ1, λ0, t)
  end

  return nothing
end




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

  c_st   ::Int64            = si
  c_time ::Float64          = 0.0
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

  return c_st
end





