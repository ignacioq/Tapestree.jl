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

















