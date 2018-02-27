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
                     idx  ::UnitRange{Int64},
                     sf   ::Int64,
                     λ1   ::Float64, 
                     λ0   ::Float64, 
                     cumts::Array{Float64,1})

  idx_end = idx[end]::Int64

  @inbounds @fastmath begin

    cur_s = Y[idx[1]]::Int64
    cur_t = 0.0
    s     = idx[2]::Int64
    idx_1 = (idx[1] - 1)::Int64

    if cur_s == 0

      while true

        cur_t += rexp(λ1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λ0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64

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

        Y[s:f] = cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λ1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64

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

    @fastmath begin

      cur_s = Y[idx[1]]::Int64
      cur_t = 0.0
      s     = idx[2]::Int64

      if cur_s == 0

        while true

          cur_t += rexp(λ1)::Float64
          f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

          Y[s:f] = cur_s::Int64
          
          if f == idx_end
            break
          end

          cur_s = (1 - cur_s)::Int64
          s     = (f + 1)::Int64

          # same but with loss rate
          cur_t += rexp(λ0)::Float64
          f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

          Y[s:f] = cur_s::Int64

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

          Y[s:f] = cur_s::Int64
          
          if f == idx_end
            break
          end

          cur_s = (1 - cur_s)::Int64
          s     = (f + 1)::Int64

          # same but with loss rate
          cur_t += rexp(λ1)::Float64
          f      = idx_1 + idxlessthan(cumts, cur_t)::Int64

          Y[s:f] = cur_s::Int64

          if f == idx_end
            break
          end

          cur_s = (1 - cur_s)::Int64
          s     = (f + 1)::Int64

        end
      end
    end

  end

end




"""
    bit_rejsam!(Y::Array{Int64,3}, idx::Array{Int64,1}, sf::Int64, λ1::Float64, ω1::Float64, λ0::Float64, ω0::Float64, avg_Δx::Float64,  cumts::Array{Float64,1})

Bit rejection-sample a single branch
given the start and end states and a
vector of cumulative δtimes incorporating information in 
`Δx` and `ω1` & `ω0`. It assigns to `Y`
in place, avoiding extra memory allocation.
"""
# reprehensible code... but faster 
function bit_rejsam!(Y     ::Array{Int64,3},
                     idx   ::UnitRange{Int64},
                     sf    ::Int64,
                     λ1    ::Float64,
                     λ0    ::Float64,
                     ω1    ::Float64,
                     ω0    ::Float64,
                     avg_Δx::Float64, 
                     cumts ::Array{Float64,1})


  @inbounds @fastmath begin

    idx_end = idx[end]::Int64
    λt1::Float64 = f_λ(λ1,ω1,avg_Δx)
    λt0::Float64 = f_λ(λ0,ω0,avg_Δx)

    cur_s = Y[idx[1]]::Int64
    cur_t = 0.0
    s     = idx[2]::Int64
    idx_1 = (idx[1] - 1)::Int64

    if cur_s == 0

      while true

        cur_t += rexp(λt1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λt0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

      end

    else

      while true

        cur_t += rexp(λt0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λt1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

      end
    end

  end

  const ntries = 0
  # rejection sampling if end simulation state do not match observed state
  while Y[idx_end] != sf 

    ntries += 1
    if ntries > 3_000_000
      # warn("bitrejsam inefficient, sampling uniformly")
      Y[rand(idx):idx_end] = sf
    end

    cur_s = Y[idx[1]]::Int64
    cur_t = 0.0
    s     = idx[2]::Int64

    if cur_s == 0

      while true

        cur_t += rexp(λt1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64
        
        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = f + 1

        # same but with loss rate
        cur_t += rexp(λt0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

      end

    else

      while true

        cur_t += rexp(λt0)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

        # same but with loss rate
        cur_t += rexp(λt1)::Float64
        f      = (idx_1 + idxlessthan(cumts, cur_t))::Int64

        Y[s:f] = cur_s::Int64

        if f == idx_end
          break
        end

        cur_s = (1 - cur_s)::Int64
        s     = (f + 1)::Int64

      end
    end

  end

end


