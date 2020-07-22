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
  [1.0/(1.0 + (T*(i-1.0))) for i in 1.0:1.0:Float64(nchains)],
  [1:nchains...]




"""
    temperature!(Ts::Array{Float64,1}, 
                 Os::Array{Int64,1},
                 T ::Float64)

Set chain temperature vectors.
"""
function temperature!(Ts::Array{Float64,1}, 
                      Os::Array{Int64,1},
                      T ::Float64)
  for (i,v) in enumerate(Os)
    Ts[i] = 1.0/(1.0 + (T*(Float64(v)-1.0)))
  end
end




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
         Ts ::Array{Float64,1})

Returns `true` if chains should be swapped.
"""
function swap(j  ::Int64,
              k  ::Int64,
              lhc::Array{Float64,1},
              Ts ::Array{Float64,1})
  @inbounds begin
    -randexp() < (Ts[j] - Ts[k]) * (lhc[k]/Ts[k] - lhc[j]/Ts[j])
  end
end





"""
    swap_chains(Os ::Array{Int64,1},
                Ts ::Array{Float64,1},
                lhc::Array{Float64,1})

Return new temperatures `Ts` and order `Os`.
"""
function swap_chains(Os ::Array{Int64,1},
                     Ts ::Array{Float64,1},
                     lhc::Array{Float64,1})

  j, k = wchains(Os)

  # if swap
  if swap(j, k, lhc, Ts)

    Oj    = Os[j]
    Os[j] = Os[k]
    Os[k] = Oj

    Tj    = Ts[j]
    Ts[j] = Ts[k]
    Ts[k] = Tj

    lhc[j] = lhc[j]/Ts[k]*Ts[j]
    lhc[k] = lhc[k]/Ts[j]*Ts[k]

    return 1.0
  end

  return 0.0
end





