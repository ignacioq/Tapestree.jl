#=

Bit rejection-sampling using discrete 
Data Augmentation for Biogeographic Histories.

Ignacio Quintero Mächler

t(-_-t)

May 01 2017

=#




# reprehensible code... but faster 
"""
    bit_rejsam!(Y    ::Array{Int64,3},
                idx  ::UnitRange{Int64},
                sf   ::Int64,
                λ1   ::Float64, 
                λ0   ::Float64, 
                cumts::Array{Float64,1})

Bit rejection-sample a single branch
given the start and end states and a
vector of cumulative δtimes. It assigns to `Y`
in place, avoiding extra memory allocation.
"""
function bit_rejsam!(Y    ::Array{Int64,3},
                     idx  ::UnitRange{Int64},
                     sf   ::Int64,
                     λ1   ::Float64, 
                     λ0   ::Float64, 
                     cumts::Array{Float64,1})

  idx_end = idx[end]::Int64

  @inbounds begin

    cur_s = Y[idx[1]]::Int64
    cur_t = 0.0
    s     = idx[2]::Int64
    idx_1 = (idx[1] - 1)::Int64

    if cur_s == 0

      while true

        cur_t += rexp(λ1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] .= cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λ0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] .= cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

      end

    else

      while true

        cur_t += rexp(λ0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] .= cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λ1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] .= cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

      end
    end

  end

  # rejection sampling if end simulation state do not match observed state
  while Y[idx_end] != sf 

    cur_s = Y[idx[1]]::Int64
    cur_t = 0.0
    s     = idx[2]::Int64

    if cur_s == 0

      while true

        cur_t += rexp(λ1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] .= cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λ0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] .= cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

      end

    else

      while true

        cur_t += rexp(λ0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] .= cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λ1)::Float64
        f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

        Y[s:f] .= cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

      end
    end

  end

end


