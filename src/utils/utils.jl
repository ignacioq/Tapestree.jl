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




"""
    update_jacobian!(μ    ::Float64,
                     σ    ::Float64, 
                     xmean::Float64, 
                     xvar ::Float64,
                     f    ::Vector{Float64},
                     j    ::Matrix{Float64},
                     jf   ::Vector{Float64})

Update jacobian for method of moments to estimating parameters of the
truncated **Gaussian** at `0.0` for lower tail.
"""
function update_jacobian!(μ    ::Float64,
                          σ    ::Float64, 
                          xmean::Float64, 
                          xvar ::Float64,
                          f    ::Vector{Float64},
                          j    ::Matrix{Float64})


  ϕ0 = stdnorm(-μ/σ)
  Φ0 = 1.0 - stpnorm(-μ/σ)

  f[1] = μ + σ * ϕ0 / Φ0 - xmean
  f[2] = σ^2 - μ*σ*ϕ0 / Φ0 - (σ*ϕ0  / Φ0)^2 - xvar

  j[1,1] = 1.0 + ϕ0*(-μ/σ * Φ0 - ϕ0)/Φ0^2
  j[2,1] = -σ*ϕ0 * (((1.0 - (μ/σ)^2) * Φ0 - μ/σ*ϕ0 )/Φ0^2 +
            2.0 * (-μ/σ * ϕ0 * Φ0 -  ϕ0^2)/Φ0^3)
  j[1,2] = ϕ0*(((μ/σ)^2 + 1.0)*Φ0 + μ^2/σ*ϕ0)/Φ0^2
  j[2,2] = 2.0*σ - μ*ϕ0*( ((μ/σ)^2 + 1.0) * Φ0 + μ^2/σ * ϕ0 )/ Φ0^2 - 
           2.0*σ*ϕ0^2 * ( ((μ/σ)^2 + 1.0) * Φ0 + μ^2/σ * ϕ0 )/ Φ0^3

end




"""
    run_newton(μ0::Float64, σ0::Float64, xmean::Float64, xvar::Float64)

Run Newton's method to estimate the roots (parameter estimates) for 
truncated **Gaussian** with lower bound at `0.0`. 
"""
function run_newton(μ0::Float64, σ0::Float64, xmean::Float64, xvar::Float64)

  @inbounds begin

    f  = Vector{Float64}(undef,2)
    j  = Matrix{Float64}(undef,2,2)
    jf = Vector{Float64}(undef,2)
    x0 = Vector{Float64}(undef,2)
    x1 = Vector{Float64}(undef,2)

    x0[1] = μ0
    x0[2] = σ0

    while true
   
      update_jacobian!(x0[1], x0[2], xmean, xvar, f, j)

      mul!(jf, inv(j), f)

      x1[1] = x0[1] - jf[1]
      x1[2] = x0[2] - jf[2]

      if (x0[1] - x1[1]) < 1e-2 && (x0[2] - x1[2]) < 1e-2
        break
      end

      x0[1] = x1[1]
      x0[2] = x1[2]
    end
  end

  return x1
end



"""
  skewness(x::Vector{Float64}, mean::Float64, sd::Float64)

Return sample skewness.
"""
@inline function skewness(x::Vector{Float64}, mean::Float64, sd::Float64)
  isd = 1.0/sd
  ss = 0.0
  @simd for xi in x
    ss += ((xi - mean)*isd)^3
  end
  ss /= Float64(lastindex(x))

  return sign(ss)*min(0.99, abs(ss))
end




"""
    mom_skewnormal(sk::Float64, mean::Float64, sd::Float64)

Return parameters for the skew normal distribution given the 
method of moments.
"""
function mom_skewnormal(sk::Float64, mean::Float64, sd::Float64)

  d = sqrt(π*0.5 * abs(sk)^(2/3) / (abs(sk)^(2/3) + ((4.0-π)*0.5)^(2/3)))
  d = sign(sk)*d
  a = d/sqrt(1.0 - d^2)
  o = sd/sqrt(1.0 - 2.0*d^2/π)
  m = mean - o*d*sqrt(2.0/π)

  return m, o, a
end


