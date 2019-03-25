#=
Create efficient approximation function

Ignacio Quintero Mächler

t(-_-t)

Created 16 09 2017
Updated 17 01 2019
=#




"""
    make_approxf(N::Int64)

Returns the closure for making an approximation function given 
`y` Array dimension.
"""
function make_approxf(x::Array{Float64,1}, 
                      y::Array{Float64,N}) where {N}

  if N == 1
    nc = 1
  else
    nc = size(y, 2)::Int64
  end

  lex1 = quote end
  pop!(lex1.args)

  if N == 1
    push!(lex1.args, :(r[1] = linpred(val, x[a], x[a+1], y[a], y[a+1])::Float64))
  else
    # unroll loop
    for i = Base.OneTo(nc)
      push!(lex1.args, :(r[$i] = linpred(val, xa, xap1, y[a,$i], y[a+1, $i])::Float64))
    end

    # add one assignment
    pushfirst!(lex1.args, :(xap1 = x[a+1]::Float64))
    pushfirst!(lex1.args, :(xa   = x[a]::Float64))
  end

  lex2 = quote end
  pop!(lex2.args)

  if N == 1
    push!(lex2.args, :(r[1] = y[a]::Float64))
  else
    # unroll loop
    for i = Base.OneTo(nc)
      push!(lex2.args, :(r[$i] = y[a,$i]::Float64))
    end
  end

  lex = quote
    a, lp = idxrange(x, val)::Tuple{Int64, Bool}
    if lp 
      $lex1 
    else 
      $lex2 
    end
  end

  # aesthetic cleaning
  deleteat!(lex.args,[1,3])

  popfirst!(lex.args[2].args[2].args)
  lex.args[2].args[2] = lex.args[2].args[2].args[1]


  popfirst!(lex.args[2].args[3].args)
  lex.args[2].args[3] = lex.args[2].args[3].args[1]

  ex = quote
    function af_gen(x::Array{Float64,1}, 
                    y::Array{Float64,$N})
      function af!(val::Float64, r::Array{Float64,1}) 
        @inbounds begin
          $lex
        end
        return nothing
      end
      return af!
    end
    af! = let x0 = $x, y0 = $y
      af_gen(x0, y0)
    end
  end

  return eval(ex)
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














