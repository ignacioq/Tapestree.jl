#=
Create efficient approximation function

Ignacio Quintero Mächler

t(-_-t)

Created 16 09 2017
Updated 17 01 2019
=#





"""
    make_approxf(x::Array{Float64,1}, y::Array{Float64,1})

Returns a linear approximation function given values of x and y.
"""
function make_approxf(x::Array{Float64,1}, 
                      y::Array{Float64,1})

  function f(val::Float64)
    @inbounds begin

      a, lp = idxrange(x, val)

      if lp == 1
        return linpred(val, x[a], x[a+1], y[a], y[a+1])::Float64
      else
        return y[a]::Float64
      end

    end
  end

  return f::Function
end





"""
    make_approxf(x::Array{Float64,1}, y::Array{Float64,2})

Returns linear approximation functions, one for each y column.
"""
function make_approxf(x::Array{Float64,1}, 
                      y::Array{Float64,2})

  nc = size(y, 2)::Int64

  function f(val::Float64, r::Array{Float64,1})

    @inbounds begin

      for i in Base.OneTo(nc)
          a, lp = idxrange(x, val)

          if lp
            r[i] = linpred(val, x[a], x[a+1], y[a,i], y[a+1, i])::Float64
          else
            r[i] = y[a,i]::Float64
          end
      end

      return nothing
    end
  end

  return f::Function
end





"""
    linpred(val::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64)

Estimate val according to linear interpolation for a range.
"""
linpred(val::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64) = 
  @fastmath (y1 + (val - x1)*(y2 - y1)/(x2 - x1))::Float64




"""
    idxrange(x::Array{Float64,1}, val::Float64)

Get indexes in sorted vector `x` corresponding to the range in which 
`val` is in using a sort of uniroot algorithm.
"""
function idxrange(x::Array{Float64,1}, val::Float64)
  
  @inbounds begin

    a::Int64 = 1

    if x[a] > val
      return a, false
    end

    b::Int64 = lastindex(x)
  
    if x[b] < val
      return b, false
    end

    mid::Int64 = div(b,2)

    while b-a > 1
      val < x[mid] ? b = mid : a = mid
      mid = div(b + a, 2)
    end

    if x[a] == val 
      return a, false
    elseif x[b] == val
      return b, false
    else
      return a, true
    end

  end
end 














