#=

ESSE equations

Ignacio Quintero Mächler

t(-_-t)

September 16 2017

=#



"""
    make_esse_s(k::Int64, x::Array{Float64,1}, y::Array{Float64,1})

Speciation ESSE equation function for integration, one z(t).
"""
function make_esse_s(k ::Int64, 
                     af::Function)

  # make Qf function
  Qf! = makeQ(k)

  # preallocate matrices & vectors
  bk::Int64 = k + (k*k)
  rt::Int64 = bk + k + 1
  Q ::Array{Float64,2} = zeros(k, k)
  Qs::Array{Float64,1} = zeros(k)
  Qe::Array{Float64,1} = zeros(k)

  # estimate indices previously
  qsi::UnitRange{Int64} = 1:k
  qei::UnitRange{Int64} = (k+1):2k

  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)

    Qf!(Q, p)
    @views BLAS.gemv!('N', 1.0, Q, u[qsi], 0.0, Qs)
    @views BLAS.gemv!('N', 1.0, Q, u[qei], 0.0, Qe)

    aft = af(t)::Float64

    for i in Base.OneTo(k)
      st      = p[i]*exp(p[bk+i]*aft)
      du[i]   =        - (st + p[k+i])*u[i]   + 2.0*st*u[i+k]*u[i]   + Qs[i]
      du[i+k] = p[k+i] - (st + p[k+i])*u[i+k] +     st*u[i+k]*u[i+k] + Qe[i]
    end

    return nothing
  end

  return f
end





"""
    make_esse_s(t,u,p,du)

Speciation ESSE equation function for integration, for `k` `z(t)` functions.
"""
function make_esse_s(k  ::Int64, 
                     af!::Function,
                     md ::Bool)

  # make Qf function
  Qf! = makeQ(k)

  # preallocate matrices & vectors
  bk::Int64 = k+k*k
  Q ::Array{Float64,2} = zeros(k, k)
  Qs::Array{Float64,1} = zeros(k)
  Qe::Array{Float64,1} = zeros(k)
  r ::Array{Float64,1} = Array{Float64}(undef,k)

  # estimate indices previously
  qsi::UnitRange{Int64} = 1:k
  qei::UnitRange{Int64} = (k+1):2k

  function f(du, u, p, t)

    Qf!(Q, p)
    @views BLAS.gemv!('N', 1.0, Q, u[qsi], 0.0, Qs)
    @views BLAS.gemv!('N', 1.0, Q, u[qei], 0.0, Qe)

    af!(t, r)

    for i in Base.OneTo(k)
      st      = p[i]*exp(p[bk+i]*r[i])::Float64
      du[i]   =        - (st + p[k+i])*u[i]   + 2.0*st*u[i+k]*u[i]   + Qs[i]
      du[i+k] = p[k+i] - (st + p[k+i])*u[i+k] +     st*u[i+k]*u[i+k] + Qe[i]
    end

    return nothing
  end

  return f
end





"""
    make_esse_e(t,u,p,du)

Extinction ESSE equation function for integration, one `z(t)`.
"""
function make_esse_e(k ::Int64, 
                     af::Function)

  # make Qf function
  Qf! = makeQ(k)

  # preallocate matrices & vectors
  bk::Int64 = k+k*k
  Q ::Array{Float64,2} = zeros(k, k)
  Qs::Array{Float64,1} = zeros(k)
  Qe::Array{Float64,1} = zeros(k)

  # estimate indices previously
  qsi::UnitRange{Int64} = 1:k
  qei::UnitRange{Int64} = (k+1):2k
  
  function f(du, u, p, t)

    Qf!(Q, p)
    @views BLAS.gemv!('N', 1.0, Q, u[qsi], 0.0, Qs)
    @views BLAS.gemv!('N', 1.0, Q, u[qei], 0.0, Qe)

    aft = af(t)::Float64

    for i in Base.OneTo(k)
      et      = p[k+i]*exp(p[bk+i]*aft)::Float64
      du[i]   =    - (p[i] + et)*u[i]   + 2.0*p[i]*u[i+k]*u[i]   + Qs[i]
      du[i+k] = et - (p[i] + et)*u[i+k] +     p[i]*u[i+k]*u[i+k] + Qe[i]
    end

    return nothing
  end

  return f
end





"""
    make_esse_e(t,u,p,du)

Extinction ESSE equation function for integration, for `k` `z(t)` functions.

"""
function make_esse_e(k  ::Int64, 
                     af!::Function,
                     md ::Bool)

  # make Qf function
  Qf! = makeQ(k)

  # preallocate matrices & vectors
  bk::Int64 = k+k*k
  Q ::Array{Float64,2} = zeros(k, k)
  Qs::Array{Float64,1} = zeros(k)
  Qe::Array{Float64,1} = zeros(k)
  r ::Array{Float64,1} = Array{Float64}(undef, k)

  # estimate indices previously
  qsi::UnitRange{Int64} = 1:k
  qei::UnitRange{Int64} = (k+1):2k

  function f(du, u, p, t)

    Qf!(Q, p)
    @views BLAS.gemv!('N', 1.0, Q, u[qsi], 0.0, Qs)
    @views BLAS.gemv!('N', 1.0, Q, u[qei], 0.0, Qe)

    af!(t, r)

    for i in Base.OneTo(k)
      et      = p[k+i]*exp(p[bk+i]*r[i])::Float64
      du[i]   =    - (p[i] + et)*u[i]   + 2.0*p[i]*u[i+k]*u[i]   + Qs[i]
      du[i+k] = et - (p[i] + et)*u[i+k] +     p[i]*u[i+k]*u[i+k] + Qe[i]
    end

    return nothing
  end

  return f
end





"""
    make_esse_q(t,u,p,du)

Rates ESSE equation function for integration for one `z(t)` function.
"""
function make_esse_q(k ::Int64, 
                     af::Function)

  # make Qt function
  Qft! = makeQt(k)

  Q ::Array{Float64,2} = zeros(k, k)
  Qs::Array{Float64,1} = zeros(k)
  Qe::Array{Float64,1} = zeros(k)

  # estimate indices previously
  qsi::UnitRange{Int64} = 1:k
  qei::UnitRange{Int64} = (k+1):2k

  function f(du, u, p, t)

    Qft!(Q, p, af(t))
    @views BLAS.gemv!('N', 1.0, Q, u[qsi], 0.0, Qs)
    @views BLAS.gemv!('N', 1.0, Q, u[qei], 0.0, Qe)

    for i in Base.OneTo(k)
      du[i]   =        - (p[i] + p[k+i])*u[i]   + 2.0*p[i]*u[i+k]*u[i]   + Qs[i]
      du[i+k] = p[k+i] - (p[i] + p[k+i])*u[i+k] +     p[i]*u[i+k]*u[i+k] + Qe[i]
    end

    return nothing
  end

  return f
end





"""
    make_esse_q(t,u,p,du)

Rates ESSE equation function for integration for multiple time functions.
"""
function make_esse_q(k  ::Int64, 
                     af!::Function,
                     md ::Bool)

  # make Qt function
  Qft! = makeQt(k, af!)

  Q ::Array{Float64,2} = zeros(k, k)
  Qs::Array{Float64,1} = zeros(k)
  Qe::Array{Float64,1} = zeros(k)

  # estimate indices previously
  qsi::UnitRange{Int64} = 1:k
  qei::UnitRange{Int64} = (k+1):2k

  function f(du, u, p, t)

    Qft!(Q, p, t) 
    @views BLAS.gemv!('N', 1.0, Q, u[qsi], 0.0, Qs)
    @views BLAS.gemv!('N', 1.0, Q, u[qei], 0.0, Qe)

    for i in Base.OneTo(k)
      du[i]   =        - (p[i] + p[k+i])*u[i]   + 2.0*p[i]*u[i+k]*u[i]   + Qs[i]
      du[i+k] = p[k+i] - (p[i] + p[k+i])*u[i+k] +     p[i]*u[i+k]*u[i+k] + Qe[i]
    end

    return nothing
  end

  return f
end





"""
    makeQt(k)

Transform instantaneous transition rates into 
the instantaneous transition matrix `Q` for **one**
`z(t)` functions.
"""
function makeQt(k::Int64)

  pidx::UnitRange{Int64} = (2k+1):(k+k*k)
  Didx::StepRange{Int64} = 1:(k+1):(k*k)
  lidx::Int64            = lastindex(pidx)
  bk  ::Int64            = k+k*k
  Qidx::Array{Int64,1}   = setdiff(1:k*k, Didx)

  function f(Q  ::Array{Float64,2}, 
             p  ::Array{Float64,1}, 
             aft::Float64)

    @inbounds begin
      # assign parameters
      for i in Base.OneTo(lidx)
        Q[Qidx[i]] = p[pidx[i]]*exp(p[bk + i]*aft)::Float64
      end

      #assign 0s
      for i in Didx
        Q[i] = 0.0
      end

      #row sums
      for j = Base.OneTo(k), i = Base.OneTo(k)
        if i == j 
          continue
        end
        Q[Didx[i]] -= Q[i,j]::Float64
      end

    end

    return nothing
  end
end





"""
    makeQt(k::Int64, af!::Function)

Transform instantaneous transition rates into 
the instantaneous transition matrix `Q` for **multiple**
`z(t)` functions.
"""
function makeQt(k::Int64, af!::Function)

  pidx::UnitRange{Int64} = (2k+1):(k+k*k)
  Didx::StepRange{Int64} = 1:(k+1):(k*k)
  lidx::Int64            = lastindex(pidx)
  bk  ::Int64            = k+k*k
  Qidx::Array{Int64,1}   = setdiff(1:k*k, Didx)
  r   ::Array{Float64,1} = Array{Float64}(undef, k*k - k)

  function f(Q  ::Array{Float64,2}, 
             p  ::Array{Float64,1},
             t  ::Float64)

    @inbounds begin
      # estimate z(t)s
      af!(t, r)

      for i in Base.OneTo(lidx)
        Q[Qidx[i]] = p[pidx[i]]*exp(p[bk+i]*r[i])::Float64
      end

      #assign 0s
      for i in Didx
        Q[i] = 0.0
      end

      #row sums
      for j = Base.OneTo(k), i = Base.OneTo(k)
        if i == j 
          continue
        end
        Q[Didx[i]] -= Q[i,j]::Float64
      end
    end

    return nothing
  end

end


