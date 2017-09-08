#=

Bit rejection-sampling using discrete 
Data Augmentation for Biogeographic Histories.

Ignacio Quintero Mächler

t(-_-t)

May 01 2017

=#



"""
    bit_rejsam!(Y::Array{Int64,3}, idx::Array{Int64,1}, sf::Int64, λ1::Float64,  λ0::Float64,  cumts::Array{Float64,1})

Bit rejection-sample a single branch
given the start and end states and a
vector of cumulative δtimes. It assigns to `Y`
in place, avoiding extra memory allocation.
"""
# reprehensible code... but faster 
function bit_rejsam!(Y    ::Array{Int64,3},
                     idx  ::Array{Int64,1},
                     sf   ::Int64,
                     λ1   ::Float64, 
                     λ0   ::Float64, 
                     cumts::Array{Float64,1})

  idx_end = idx[end]::Int64

  @inbounds @fastmath begin

    cur_s = Y[idx[1]]::Int64
    cur_t = 0.0
    s     = idx[2]::Int64
    idx_1 = idx[1] - 1

    if cur_s == 0

      while true

        cur_t += rexp(λ1)::Float64
        f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

        Y[s:f] = cur_s
        
        if f == idx_end
          break
        end

        cur_s = 1 - cur_s
        s     = f + 1

        # same but with loss rate
        cur_t += rexp(λ0)::Float64
        f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

        Y[s:f] = cur_s

        if f == idx_end
          break
        end

        cur_s = 1 - cur_s
        s     = f + 1

      end

    else

      while true

        cur_t += rexp(λ0)::Float64
        f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

        Y[s:f] = cur_s
        
        if f == idx_end
          break
        end

        cur_s = 1 - cur_s
        s     = f + 1

        # same but with loss rate
        cur_t += rexp(λ1)::Float64
        f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

        Y[s:f] = cur_s

        if f == idx_end
          break
        end

        cur_s = 1 - cur_s
        s     = f + 1

      end
    end

  end

  # rejection sampling if end simulation state do not match observed state
  while Y[idx_end] != sf 

    @fastmath begin

      cur_s = Y[idx[1]]::Int64
      cur_t = 0.0
      s     = idx[2]

      if cur_s == 0

        while true

          cur_t += rexp(λ1)::Float64
          f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

          Y[s:f] = cur_s
          
          if f == idx_end
            break
          end

          cur_s = 1 - cur_s
          s     = f + 1

          # same but with loss rate
          cur_t += rexp(λ0)::Float64
          f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

          Y[s:f] = cur_s

          if f == idx_end
            break
          end

          cur_s = 1 - cur_s
          s     = f + 1

        end

      else

        while true

          cur_t += rexp(λ0)::Float64
          f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

          Y[s:f] = cur_s
          
          if f == idx_end
            break
          end

          cur_s = 1 - cur_s
          s     = f + 1

          # same but with loss rate
          cur_t += rexp(λ1)::Float64
          f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

          Y[s:f] = cur_s

          if f == idx_end
            break
          end

          cur_s = 1 - cur_s
          s     = f + 1

        end
      end
    end

  end

end







#=
  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  Deprecated functions for discrete biogeographic sampling
  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#


# """
#   assigns discrete values according to 
#   the continuous sampling to Yc
# """
# function assigndisceve!(si     ::Int64, 
#                         Y      ::Array{Int64,3}, 
#                         contsam::Array{Float64,1}, 
#                         bridx  ::Array{Int64,1}, 
#                         δtvec  ::Array{Float64,1})
#   s    ::Int64 = 1
#   cur_s::Int64 = si
 
#   @inbounds begin

#     lbr = endof(bridx)
    
#     for i=eachindex(contsam)
#       f = indmindif_sorted(δtvec, contsam[i])
#       setindex!(Y, cur_s, bridx[s:f]) 
#       cur_s = 1 - cur_s
#       s     = f == lbr ? f : (f + 1)
#     end

#   end
# end






# """
#   rejection sampling for each branch
#   condition on start and end point
# """
# function rejsam_cumsum(si::Int64, 
#                        sf::Int64, 
#                        λ1::Float64, 
#                        λ0::Float64, 
#                        t::Float64)
  
#   evs, ssf = brprop_cumsum(si, λ1, λ0, t)

#   while ssf != sf 
#     evs, ssf = brprop_cumsum(si, λ1, λ0, t)
#   end

#   return evs
# end





# """
#   propose events for a branch
#   with equal rates for gain and loss
#   return the cumsum of times
#   * Ugly code but slightly faster *
# """
# function brprop_cumsum(si::Int64, λ1::Float64, λ0::Float64, t::Float64)

#   cur_s::Int64            = si
#   cur_t::Float64          = 0.
#   cum_t::Array{Float64,1} = Float64[]

#   if cur_s == 0
#     cur_t += rexp(λ1)

#     while cur_t < t
#       push!(cum_t, cur_t)
#       cur_s  = 1 - cur_s
#       cur_t += rexp(λ0)
    
#       if cur_t > t
#         break
#       end

#       push!(cum_t, cur_t)
#       cur_s  = 1 - cur_s
#       cur_t += rexp(λ1)
#     end

#   else
#     cur_t += rexp(λ0)

#     while cur_t < t
#       push!(cum_t, cur_t)
#       cur_s  = 1 - cur_s
#       cur_t += rexp(λ1)
    
#       if cur_t > t
#         break
#       end

#       push!(cum_t, cur_t)
#       cur_s  = 1 - cur_s
#       cur_t += rexp(λ0)
#     end

#   end

#   push!(cum_t, t)

#   return cum_t, cur_s
# end

