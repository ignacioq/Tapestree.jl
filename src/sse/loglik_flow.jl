#=

Log-likelihood for ODEs, integrating along the tree
using flow algorithm

Ignacio Quintero Mächler

t(-_-t)

Created 23 10 2019

=#





"""
    make_Et(Et::Function,
            p0::Array{Float64,1}, 
            u0::Array{Float64,1},
            ts::Array{Float64,1},
            ti::Float64,

Make Extinction through time, `E(t)`.
"""
function make_Et(Et::Function,
                 p0::Array{Float64,1}, 
                 u0::Array{Float64,1},
                 ts::Array{Float64,1},
                 ti::Float64,
                 tf::Float64)

  prob = ODEProblem(Et, u0, (ti,tf), p0)

  int = init(prob,
             Tsit5(),
             saveat          = ts,
             save_everystep  = false, 
             calck           = false,
             force_dtmin     = true,
             save_start      = true,
             initialize_save = false,
             maxiters        = 100_000_000,
             verbose         = false)

  function f(p::Array{Float64,1})
    int.p = p
    reinit!(int)
    solve!(int).u::Array{Array{Float64,1},1}
  end

  return f
end





"""
    make_Gt(At::Function, 
            p0::Array{Array{Float64,1},1}, 
            u0::Array{Float64,2},
            ts::Array{Float64,1},
            ti::Float64,
            tf::Float64)

Make flow equation for likelihoods through time, `G(t)`.
"""
function make_Gt(At::Function, 
                 p0::Array{Array{Float64,1},1}, 
                 u0::Array{Float64,2},
                 ts::Array{Float64,1},
                 ti::Float64,
                 tf::Float64)

  prob = ODEProblem(At, u0, (ti,tf), p0)

  int = init(prob,
             Tsit5(),
             saveat          = ts,
             save_everystep  = false, 
             calck           = false,
             force_dtmin     = true,
             save_start      = true,
             initialize_save = false,
             maxiters        = 100_000_000,
             verbose         = false)

  function f(p::Array{Array{Float64,1},1})
    int.p = p
    reinit!(int)
    solve!(int).u::Array{Array{Float64,2},1}
  end

  return f
end





"""
    make_loglik(Gt     ::Function, 
                Et     ::Function,
                X      ::Array{Array{Float64,1},1}, 
                triads ::Array{Array{Int64,1},1},
                lbts   ::Array{Int64,2},
                bts    ::Array{Float64,1},
                ns     ::Int64, 
                ned    ::Int64, 
                nets   ::Int64,
                ntip   ::Int64,
                λevent!::Function,
                rootll ::Function)


Make log-likelihood function using the flow algorithm.
"""
function make_loglik(Gt     ::Function, 
                     Et     ::Function,
                     X      ::Array{Array{Float64,1},1}, 
                     triads ::Array{Array{Int64,1},1},
                     lbts   ::Array{Int64,2},
                     bts    ::Array{Float64,1},
                     ns     ::Int64, 
                     ned    ::Int64, 
                     nets   ::Int64,
                     λevent!::Function,
                     rootll ::Function)

  # preallocate arrays
  X0  = Array{Float64,1}(undef,ns)
  Xp1 = Array{Float64,1}(undef,ns)
  Xp2 = Array{Float64,1}(undef,ns)
  Xr  = Array{Float64,1}(undef,ns)
  wg  = Array{Float64,1}(undef,ns)
  Er  = Array{Float64,1}(undef,ns)

  rtime  = bts[end]
  netsm1 = nets - 1
  nbts   = length(bts)

  # start ll function for parameters
  function f(p::Array{Float64,1})

    @inbounds begin

      Ets = Et(p)

      # push parameters as the last vector in Ets
      push!(Ets, p)

      # estimate `Gts` according to `Ets` and `p`
      Gts = Gt(Ets)

      ll = 0.0

      if rank(Gts[nbts]) == ns

        for i in Base.OneTo(nbts) 
          if rank(Gts[i]) != ns
            @show Gts[nbts]
            @show Gts[i]
          end

        end

        for trio in triads

          pr, d1, d2 = trio

          ## first daughter
          X0 = Gts[lbts[d1,2]]\X[d1]
          mul!(Xp1, Gts[lbts[d1,1]], X0)

          check_negs(Xp1, ns) && return -Inf

          ## second daughter
          X0 = Gts[lbts[d2,2]]\X[d2]
          mul!(Xp2, Gts[lbts[d2,1]], X0)

          check_negs(Xp2, ns) && return -Inf

          λtime = bts[lbts[d1,1]]

          λevent!(λtime, X, Xp1, Xp2, p, pr)

          if pr != ned
            mdlus  = sum(X[pr])
            X[pr] /= mdlus
            ll += log(mdlus)
          end
        end
      else

        for trio in triads

          pr, d1, d2 = trio

          ## first daughter
          ldiv!(X0, qr!(Gts[lbts[d1,2]], Val(true)), X[d1])
          mul!(Xp1, Gts[lbts[d1,1]], X0)

          check_negs(Xp1, ns) && return -Inf

          ## second daughter
          ldiv!(X0, qr!(Gts[lbts[d2,2]], Val(true)), X[d2])
          mul!(Xp2, Gts[lbts[d2,1]], X0)

          check_negs(Xp2, ns) && return -Inf

          λtime = bts[lbts[d1,1]]

          λevent!(λtime, X, Xp1, Xp2, p, pr)

          if pr != ned
            mdlus  = sum(X[pr])
            X[pr] /= mdlus
            ll += log(mdlus)
          end
        end

      end

      # assign root likelihoods
      @simd for i in Base.OneTo(ns)
        Xr[i] = X[ned][i]
      end

      # assign root extinction
      @simd for i in Base.OneTo(ns)
        Er[i] = Ets[netsm1][i]
      end

      # estimate weights
      normbysum!(Xr, wg, ns)

      ll += log(rootll(rtime, Xr, Er, wg, p)::Float64)
    end

    return ll
  end

  return f
end






