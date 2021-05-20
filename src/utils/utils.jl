#=

Utility functions for Tapestree

Ignacio Quintero Mächler

February 6 2017

t(-_-t)

=#


"""
    rowind(x::Int64, nrow::Int64)

Get row indexing from matrix indexing.
"""
rowind(x::Int64, nrow::Int64) = mod1(x,nrow)





"""
    colind(x::Int64, nrow::Int64)

Get column indexing from matrix indexing
"""
colind(x::Int64, nrow::Int64) = cld(x, nrow)





"""
    vecind(row::Int64, col::Int64, nrow::Int64)

Get vector indexing from column and row.
"""
vecind(row::Int64, col::Int64, nrow::Int64) = row + nrow*(col - 1)





"""
    uniroot(f, approx = 1e-8, a = 0.0, b = 0.1)

Find the root of function between `0.0` and `b`.
"""
function uniroot(f; approx = 1e-8, a = 0.0, b = 0.1) 
  # choose b
  while sign(f(a)::Float64)::Float64 == sign(f(b)::Float64)::Float64
    b += 0.1
  end
  m::Float64 = (a + b)/2.0::Float64

  while abs(f(m)::Float64)::Float64 > approx
    if sign(f(a)::Float64)::Float64 == sign(f(m)::Float64)::Float64
      a = m::Float64
    else 
      b = m::Float64
    end
    m = (a + b)/2.0::Float64
  end
  return m::Float64
end 






"""
    indmindif_sorted(x::Array{Float64,1}, val::Float64)

Get index in sorted vector `x` corresponding to the value 
that is closest to `val` in sorted arrays 
using a sort of uniroot algorithm.
"""
function indmindif_sorted(x::Array{Float64,1}, val::Float64) 
  a::Int64   = 1
  b::Int64   = lastindex(x)
  mid::Int64 = div(b,2)  

  while b-a > 1
    val < x[mid] ? b = mid : a = mid
    mid = div(b + a, 2)
  end

  abs(x[a] - val) < abs(x[b] - val) ? a : b
end 




"""
    idxlessthan(x::Array{Float64,1}, val::Float64)

Get index in sorted vector `x` corresponding to the value 
that is closest to but less than `val` in sorted arrays 
using a sort of uniroot algorithm.
"""
function idxlessthan(x::Array{Float64,1}, val::Float64) 
  
  @inbounds begin

    a::Int64 = 1
    b::Int64 = lastindex(x)
  
    if x[b] < val
      return b
    end

    mid::Int64 = div(b,2)  

    while b-a > 1
      val < x[mid] ? b = mid : a = mid
      mid = div(b + a, 2)
    end

  end

  return a
end 





"""
    fIrand(s::Int64)

Fast sampling of an Integer from `0:(s-1)` based on 
`https://discourse.julialang.org/t/rand-1-10-vs-int-round-10-rand/14339/9`
"""
@inline function fIrand(s::Int64)
    __nearlydivisionless(s, UInt64, UInt128, UInt128(1) << 64, 64)
end




"""
    __nearlydivisionless(s::Integer, T, T2, maxval, nbits)

Fast sampling of an Integer from `1:s` based on 
`https://discourse.julialang.org/t/rand-1-10-vs-int-round-10-rand/14339/9`
"""
@inline function __nearlydivisionless(s::Integer, T, T2, maxval, nbits)
    l, m = drawl(s, maxval, T, T2)
    if (l < s)
        t = maxval % s
        while l < t
            l, m = drawl(s, maxval, T, T2)
        end
    end
    return Int64(m >> nbits)
end




"""
    drawl(s::Integer, maxval, T, T2)

Fast sampling of an Integer from `1:s` based on 
`https://discourse.julialang.org/t/rand-1-10-vs-int-round-10-rand/14339/9`
"""
@inline function drawl(s::Integer, maxval, T, T2)
    x = rand(T)
    m = convert(T2, x) * s
    l = m % maxval
    return (l, m)
end












