#=

Markov coupled slice-sampling

Ignacio Quintero Mächler

t(-_-t)

21 07 2020

=#


"""
  make_temperature(T::Float64, nchains::Int64)

Set chain temperature vectors.
"""
make_temperature(T::Float64, nchains::Int64) = 
  [1.0/(1.0 + (T*(i-1))) for i in 1.0:1.0:Float64(nchains)],
  [1:nchains...]




"""
    wchains(Os::Array{Int64,1})

Sample two chains to swap.
"""
function wchains(Os::Array{Int64,1})
  j = rand(Os)
  k = rand(Os)
  while k == j
    k = rand(Os)
  end
  return j, k
end




"""
    swap(j  ::Int64, 
         k  ::Int64, 
         lhc::Array{Float64,1},
         p  ::Array{Array{Float64,1},1},
         fp ::Array{Array{Float64,1},1},
         Ts ::Array{Float64,1}
         lhf::Function)

Returns `true` if chains should be swapped.
"""
function swap(j  ::Int64,
              k  ::Int64,
              lhc::Array{Float64,1},
              p  ::Array{Array{Float64,1},1},
              fp ::Array{Array{Float64,1},1},
              Ts ::Array{Float64,1},
              lhf::Function)
  @inbounds begin
    -randexp() < 
       (lhf(p[k], fp[k], T[j]) - lhc[j] + 
        lhf(p[j], fp[j], T[k]) - lhc[k])
  end
end




"""
    swap_chains(Os ::Array{Int64,1},
                Ts ::Array{Float64,1},
                lhc::Array{Float64,1},
                p  ::Array{Array{Float64,1},1},
                fp ::Array{Array{Float64,1},1},
                lhf::Function)

Return new temperatures `Ts` and order of temperatures `Os`.
"""
function swap_chains(Os ::Array{Int64,1},
                     Ts ::Array{Float64,1},
                     lhc::Array{Float64,1},
                     p  ::Array{Array{Float64,1},1},
                     fp ::Array{Array{Float64,1},1},
                     lhf::Function)

  j, k =  wchains(Os)

  # if swap
  if swap(j, k, lhc, p, fp, Ts, lhf)

    Os[j] = k
    Os[k] = j

    Tj    = Ts[j]
    Ts[j] = Ts[k]
    Ts[k] = Tj
  end

  return Os, Ts
end





