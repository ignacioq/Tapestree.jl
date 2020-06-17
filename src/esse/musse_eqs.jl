#=

MuSSE equations

Ignacio Quintero Mächler

t(-_-t)

September 16 2017

=#





"""
    make_musse(t,u,p,du)

MuSSE equation function for integration.
"""
function make_musse(k::Int64)

  Qf! = makeQ(k)

  # preallocate matrices & vectors
  Q ::Array{Float64,2} = zeros(k, k)
  Qs::Array{Float64,1} = zeros(k)
  Qe::Array{Float64,1} = zeros(k)

  # estimate indices previously
  qsi::UnitRange{Int64} = 1:k
  qei::UnitRange{Int64} = (k+1):2k

  function f(du, u, p, t)

    Qf!(Q, p)
    BLAS.gemv!('N', 1.0, Q, u[qsi], 0.0, Qs)
    BLAS.gemv!('N', 1.0, Q, u[qei], 0.0, Qe)

    for i in Base.OneTo(k)
      du[i]   =        - (p[i] + p[k+i])*u[i]   + 2.0*p[i]*u[i+k]*u[i]   + Qs[i]
      du[i+k] = p[k+i] - (p[i] + p[k+i])*u[i+k] +     p[i]*u[i+k]*u[i+k] + Qe[i]
    end

  end

  return f
end





"""
    makeQ(k::Int64)

Transform instantaneous transition rates into the instantaneous 
transition matrix `Q`.
"""
function makeQ(k::Int64)

  pidx::UnitRange{Int64}  = (2k+1):(k + k*k)
  Didx::StepRange{Int64}  = 1:(k+1):(k*k)
  Qidx::Array{Int64,1}    = setdiff(1:k*k, Didx)  
  pitr::Base.OneTo{Int64} = Base.OneTo(lastindex(pidx))
  kitr::Base.OneTo{Int64} = Base.OneTo(k)

  function f(Q::Array{Float64,2}, p::Array{Float64,1})

    @inbounds begin
      # assign parameters
      @simd for i in pitr
        Q[Qidx[i]] = p[pidx[i]]
      end

      #assign 0s
      for i in Didx
        Q[i] = 0.0
      end

      #row sums
      for j = kitr, i = kitr
        if i == j 
          continue
        end
        Q[Didx[i]] -= Q[i,j]
      end
    end

  return nothing
  end
end


