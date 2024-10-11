#=

insane lineage through time structure `Ltt`

Ignacio Quintero Mächler

t(-_-t)

Created 17 08 2021
=#




"""
    Ltt

A Composite type representing lineage through time:

  `n`: number of lineages at times `t`.
  `t`: times (right piece-wise).
"""
struct Ltt
  n::Array{Int64,1}
  t::Array{Float64,1}
end

# pretty-printing
Base.show(io::IO, x::Ltt) =
  print(io, "LTT starting at ", round(maximum(x.t), digits = 3), 
            " with ", x.n[end], " extant tips")




"""
    times_n(n::Int64, nt::Ltt

Return tuple(s) of start and end of times for which there were `n` lineages.
"""
@inline function times_n(n::Int64, nt::Ltt)

  ix = findall(isequal(n), nt.n)

  t = NTuple{2,Float64}[]
  @simd for i in ix
    push!(t, (nt.t[i], nt.t[i+1]))
  end

  return t
end




"""
    usample(t::Vector{Tuple{Float64, Float64}}, p::Float64)

Sample uniformly from a vector of start and end times given probability `p` of
sampling over the sum of those times.
"""
function usample(t::Vector{Tuple{Float64, Float64}}, p::Float64)

  tt = 0.0
  for (tii, tff) in t
    tt += tii - tff
  end

  # sample
  if rand() < p*tt

    s = rand()*tt

    # find true t
    tt = 0.0
    for (tii, tff) in t
      tt += tii - tff

      if s < tt
        return tii - s
      end
    end
  end

  return 0.0
end




"""
    nspt(x::Ltt, t::Float64) = x.n[searchsortedlast(x.t, t, rev = true)]

Return the number of species at time t
"""
function nspt(x::Ltt, t::Float64)
  ix = searchsortedlast(x.t, t, rev = true)
  if iszero(ix)
    return 0
  else
    return x.n[searchsortedlast(x.t, t, rev = true)]
  end
end

