#=

Markov coupled slice-sampling

Ignacio Quintero Mächler

t(-_-t)

21 07 2020

=#


"""
  make_temperature(dt::Float64, ncch::Int64)

Set chain temperature vectors.
"""
make_temperature(dt::Float64, ncch::Int64) = 
  [1.0/(1.0 + (dt*(i-1.0))) for i in 1.0:1.0:Float64(ncch)],
  [1:ncch...]




"""
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
  end
end




"""
    wchains(o::Array{Int64,1})

Sample two chains to swap.
"""
function wchains(o::Array{Int64,1})
  j = rand(o)
  k = rand(o)
  while k == j
    k = rand(o)
  end
  return j, k
end




"""
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

Returns `true` if chains should be swapped.
"""
function swap(j  ::Int64,
              k  ::Int64,
              lhc::Array{Float64,1},
              o  ::Array{Int64,1},
              t  ::Array{Float64,1})
  @inbounds begin
    -randexp() < ((t[o[j]] - t[o[k]]) * (lhc[k]/t[o[k]] - lhc[j]/t[o[j]]))
  end
end




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
              o  ::SharedArray{Int64,1},
              t  ::Array{Float64,1})
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
end




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


