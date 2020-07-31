#=

Markov coupled slice-sampling

Ignacio Quintero Mächler

t(-_-t)

21 07 2020

=#


"""
<<<<<<< HEAD
  make_temperature(T::Float64, nchains::Int64)

Set chain temperature vectors.
"""
make_temperature(T::Float64, nchains::Int64) = 
  [1.0/(1.0 + (T*(i-1.0))) for i in 1.0:1.0:Float64(nchains)],
  [1:nchains...]
=======
  make_temperature(dt::Float64, ncch::Int64)

Set chain temperature vectors.
"""
make_temperature(dt::Float64, ncch::Int64) = 
  [1.0/(1.0 + (dt*(i-1.0))) for i in 1.0:1.0:Float64(ncch)],
  [1:ncch...]
>>>>>>> mc3




"""
<<<<<<< HEAD
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
=======
    temperature!(t::Array{Float64,1}, 
                 o ::Array{Int64,1},
                 dt ::Float64)

Set chain temperature vectors.
"""
function temperature!(t::Array{Float64,1}, 
                      o::Array{Int64,1},
                      dt ::Float64)
  for (i,v) in enumerate(o)
    t[i] = 1.0/(1.0 + (dt*(Float64(v)-1.0)))
>>>>>>> mc3
  end
end




"""
<<<<<<< HEAD
    wchains(Os::Array{Int64,1})

Sample two chains to swap.
"""
function wchains(Os::Array{Int64,1})
  j = rand(Os)
  k = rand(Os)
  while k == j
    k = rand(Os)
=======
    wchains(o::Array{Int64,1})

Sample two chains to swap.
"""
function wchains(o::Array{Int64,1})
  j = rand(o)
  k = rand(o)
  while k == j
    k = rand(o)
>>>>>>> mc3
  end
  return j, k
end




"""
<<<<<<< HEAD
    swap(j  ::Int64,
         k  ::Int64,
         lhc::Array{Float64,1},
         Ts ::Array{Float64,1})
=======
    wchains(o::SharedArray{Int64,1})

Sample two chains to swap.
"""
function wchains(o::SharedArray{Int64,1})
  j = rand(o)
  k = rand(o)
  while k == j
    k = rand(o)
  end
  return j, k
end



"""
    swap(j  ::Int64,
         k  ::Int64,
         lhc::Array{Float64,1},
         t ::Array{Float64,1})
>>>>>>> mc3

Returns `true` if chains should be swapped.
"""
function swap(j  ::Int64,
              k  ::Int64,
              lhc::Array{Float64,1},
<<<<<<< HEAD
              Ts ::Array{Float64,1})
  @inbounds begin
    -randexp() < (Ts[j] - Ts[k]) * (lhc[k]/Ts[k] - lhc[j]/Ts[j])
=======
              o ::Array{Int64,1},
              t ::Array{Float64,1})
  @inbounds begin
    -randexp() < ((t[o[j]] - t[o[k]]) * (lhc[k]/t[o[k]] - lhc[j]/t[o[j]]))
>>>>>>> mc3
  end
end




<<<<<<< HEAD

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
=======
"""
    swap(j  ::Int64,
         k  ::Int64,
         lhc::SharedArray{Float64,1},
         o::SharedArray{Int64,1},
         t ::Array{Float64,1})

Returns `true` if chains should be swapped.
"""
function swap(j  ::Int64,
              k  ::Int64,
              lhc::SharedArray{Float64,1},
              o ::SharedArray{Int64,1},
              t ::Array{Float64,1})
  @inbounds begin
    -randexp() < ((t[o[j]] - t[o[k]]) * (lhc[k]/t[o[k]] - lhc[j]/t[o[j]]))
  end
end




"""
    swap_chains!(j  ::Int64,
                 k  ::Int64,
                 o ::Array{Int64,1},
                 t ::Array{Float64,1},
                 lhc::Array{Float64,1})

Return new temperature order `o`.
"""
function swap_chains!(o  ::Array{Int64,1},
                      t  ::Array{Float64,1},
                      lhc::Array{Float64,1})
  j, k = wchains(o)

  if swap(j, k, lhc, o, t)
    oj   = o[j]
    o[j] = o[k]
    o[k] = oj

    lhc[j] = lhc[j]/t[o[k]]*t[o[j]]
    lhc[k] = lhc[k]/t[o[j]]*t[o[k]]
  end

  return nothing
>>>>>>> mc3
end




<<<<<<< HEAD
=======
"""
    swap_chains!(j  ::Int64,
                 k  ::Int64,
                 o  ::SharedArray{Int64,1},
                 t  ::Array{Float64,1},
                 lhc::SharedArray{Float64,1})

Return new temperature order `o`.
"""
function swap_chains!(j  ::Int64,
                      k  ::Int64,
                      o  ::SharedArray{Int64,1},
                      t  ::Array{Float64,1},
                      lhc::SharedArray{Float64,1})

  if swap(j, k, lhc, o, t)
    oj   = o[j]
    o[j] = o[k]
    o[k] = oj

    lhc[j] = lhc[j]/t[o[k]]*t[o[j]]
    lhc[k] = lhc[k]/t[o[j]]*t[o[k]]
  end

  return true
end

>>>>>>> mc3

