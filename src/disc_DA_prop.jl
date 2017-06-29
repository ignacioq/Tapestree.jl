"""

Utilities for Discrete Data Augmentation in 
Biogeographic Histories.

Ignacio Quintero Mächler

t(-_-t)

May 01 2017

"""





"""
  rejection sampling for each branch
  condition on start and end point
  for bit sampling
"""
function bit_rejsam!(bitv ::Array{Int64,1},
                     sf   ::Int64,
                     λ1   ::Float64, 
                     λ0   ::Float64, 
                     cumts::Array{Float64,1})
  
  bit_prop_hist!(bitv, λ1, λ0, cumts)

  while bitv[end] != sf 
    bit_prop_hist!(bitv, λ1, λ0, cumts)
  end

end





"""
  function for proposing bit histories 
  according to cumulative δtimes and assigning
  to bitv
  * Ugly code but slightly faster *
"""
function bit_prop_hist!(bitv ::Array{Int64,1},
                        λ1   ::Float64, 
                        λ0   ::Float64, 
                        cumts::Array{Float64,1})

  @fastmath begin

    lbitv = endof(bitv)::Int64
    cur_s = bitv[1]::Int64
    cur_t = 0.0
    s     = 2

    if cur_s == 0

      while true

        cur_t += rexp(λ1)::Float64
        f      = idxlessthan(cumts, cur_t)::Int64

        bitv[s:f] = cur_s
        
        if f == lbitv
          break
        end

        cur_s = 1 - cur_s
        s     = f + 1

        # same but with loss rate
        cur_t += rexp(λ0)::Float64
        f      = idxlessthan(cumts, cur_t)::Int64

        bitv[s:f] = cur_s

        if f == lbitv
          break
        end

        cur_s = 1 - cur_s
        s     = f + 1

      end

    else

      while true

        cur_t += rexp(λ0)::Float64
        f      = idxlessthan(cumts, cur_t)::Int64

        bitv[s:f] = cur_s
        
        if f == lbitv
          break
        end

        cur_s = 1 - cur_s
        s     = f + 1

        # same but with loss rate
        cur_t += rexp(λ1)::Float64
        f      = idxlessthan(cumts, cur_t)::Int64

        bitv[s:f] = cur_s

        if f == lbitv
          break
        end

        cur_s = 1 - cur_s
        s     = f + 1

      end
    end

  end
end





"""
  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  Deprecated functions for discrete biogeographic sampling
  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
"""



"""
  rejection sampling for each branch
  condition on start and end point
"""
function rejsam_cumsum(si::Int64, 
                       sf::Int64, 
                       λ1::Float64, 
                       λ0::Float64, 
                       t::Float64)
  
  evs, ssf = brprop_cumsum(si, λ1, λ0, t)

  while ssf != sf 
    evs, ssf = brprop_cumsum(si, λ1, λ0, t)
  end

  return evs
end





"""
  propose events for a branch
  with equal rates for gain and loss
  return the cumsum of times
  * Ugly code but slightly faster *
"""
function brprop_cumsum(si::Int64, λ1::Float64, λ0::Float64, t::Float64)

  cur_s::Int64            = si
  cur_t::Float64          = 0.
  cum_t::Array{Float64,1} = Float64[]

  if cur_s == 0
    cur_t += rexp(λ1)

    while cur_t < t
      push!(cum_t, cur_t)
      cur_s  = 1 - cur_s
      cur_t += rexp(λ0)
    
      if cur_t > t
        break
      end

      push!(cum_t, cur_t)
      cur_s  = 1 - cur_s
      cur_t += rexp(λ1)
    end

  else
    cur_t += rexp(λ0)

    while cur_t < t
      push!(cum_t, cur_t)
      cur_s  = 1 - cur_s
      cur_t += rexp(λ1)
    
      if cur_t > t
        break
      end

      push!(cum_t, cur_t)
      cur_s  = 1 - cur_s
      cur_t += rexp(λ0)
    end

  end

  push!(cum_t, t)

  return cum_t, cur_s
end

