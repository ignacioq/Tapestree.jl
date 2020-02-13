#=

Log-likelihood for ODEs, integrating along the tree
using flow algorithm

Ignacio Quintero Mächler

t(-_-t)

Created 23 10 2019

=#

# λ0, λ1, μ0, μ1, q01, q10

using DifferentialEquations
using LinearAlgebra

function fbisse(du::Array{Float64,1}, 
                u::Array{Float64,1}, 
                p::Array{Float64,1}, 
                t::Float64)
  @inbounds begin
    du[1]= -1.0 * (p[1] + p[3] + p[5]) * u[1] + p[5] * u[2] + 2.0 * p[1] * u[3] * u[1]
    du[2]= -1.0 * (p[2] + p[4] + p[6]) * u[2] + p[6] * u[1] + 2.0 * p[2] * u[4] * u[2]
    du[3]= -1.0 * (p[1] + p[3] + p[5]) * u[3] + p[3] + p[5] * u[4] + p[1] * u[3]^2
    du[4]= -1.0 * (p[2] + p[4] + p[6]) * u[4] + p[4] + p[6] * u[3] + p[2] * u[4]^2
  end
  
  return nothing
end

function fbisseE(du::Array{Float64,1}, 
                 u::Array{Float64,1}, 
                 p::Array{Float64,1}, 
                 t::Float64)
  @inbounds begin
    du[1] = -1.0 * (p[1] + p[3] + p[5]) * u[1] + p[3] + p[5] * u[2] + p[1] * u[1]^2
    du[2] = -1.0 * (p[2] + p[4] + p[6]) * u[2] + p[4] + p[6] * u[1] + p[2] * u[2]^2
  end
  
  return nothing
end


function make_fbisseM(afE!, idxl::Int64)

  rE  = Array{Float64,1}(undef,2)

  function f(du::Array{Float64,2}, 
             u::Array{Float64,2}, 
             p::Array{Array{Float64,1},1}, 
             t::Float64)
    @inbounds begin
      par = p[idxl]
      afE!(t, rE, p)
      du[1,1] = -1.0 * (par[1] + par[3] + par[5]) * u[1,1] + par[5] * u[2,1] + 2.0 * par[1] * rE[1] * u[1,1]
      du[2,1] = -1.0 * (par[2] + par[4] + par[6]) * u[2,1] + par[6] * u[1,1] + 2.0 * par[2] * rE[2] * u[2,1]
      du[1,2] = -1.0 * (par[1] + par[3] + par[5]) * u[1,2] + par[5] * u[2,2] + 2.0 * par[1] * rE[1] * u[1,2]
      du[2,2] = -1.0 * (par[2] + par[4] + par[6]) * u[2,2] + par[6] * u[1,2] + 2.0 * par[2] * rE[2] * u[2,2]
    end

    return nothing
  end
end


# make E function
function make_Et(Et,
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


# make E or G functions
function make_Gt(At, 
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
make_loglik(ode_make_ext,
            ode_make_lik,
            tv  ::Dict{Int64,Array{Float64,1}},
            ed  ::Array{Int64,2},
            el  ::Array{Float64,1},
            bts ::Array{Float64,1},
            p0  ::Array{Float64,1},
            E0  ::Array{Float64,1},
            k   ::Int64,
            h   ::Int64,
            ntip::Int64;
            Eδt::Float64 = 0.01,
            ti ::Float64 = 0.0)

Make log-likelihood function using the flow algorithm.
"""
function make_loglik(ode_make_lik,
                     ode_ext,
                     tv  ::Dict{Int64,Array{Float64,1}},
                     ed  ::Array{Int64,2},
                     el  ::Array{Float64,1},
                     bts ::Array{Float64,1},
                     p0  ::Array{Float64,1},
                     E0  ::Array{Float64,1},
                     k   ::Int64,
                     h   ::Int64,
                     ntip::Int64;
                     Eδt::Float64 = 0.01,
                     ti ::Float64 = 0.0)





"""
*** all the MAKE functions and just pass the MADE functions ***

I should just put this out ----->>>>
"""


  # number of states
  ns = k*h

  # Estimate extinction at `ts` times
  tf = maximum(bts)
  ts = [ti:Eδt:tf...]

  # Make extinction integral
  Et = make_Et(ode_ext, p0, E0, ts, ti, tf)

  # estimate first Extinction probabilities at times `ts`
  Ets  = Et(p0)
  nets = length(Ets) + 1

  # Make extinction approximated function
  # ** this make order is crucial
  afE! = make_af(ts,  Ets, Val(ns))

  # make likelihood integral
  ode_intf = ode_make_lik(afE!, nets)

  # push parameters as the last vector in Ets
  push!(Ets, p0)

  # make Gt function
  Gt = make_Gt(ode_intf, Ets, Matrix{Float64}(I, ns, ns), bts, ti, tf)


"""
<<<<----- until here (or so)
"""




  # estimate `Gts` according to `Ets` and `p0`
  Gts = Gt(Ets)

  # sort branching times
  sort!(bts)

  # add the present `0.0`
  pushfirst!(bts, 0.0)

  nbts = lastindex(bts)

  ## link edges with initial and end branching times 
  abts = abs_time_branches(el, ed, ntip)
  lbts = Array{Int64,2}(undef,size(abts))

  # preallocate absolute minimum matrix of differences
  minM = Array{Float64,1}(undef, nbts)
  # find links
  for i in Base.OneTo(ned*2)
    for j in Base.OneTo(nbts)
      minM[j] = abs(abts[i] - bts[j])
    end
    lbts[i] = argmin(minM)
  end

  # make internal node triads
  triads = maketriads(ed)

  # initialize likelihood vectors
  X = [Array{Float64,1}(undef,ns) for i in Base.OneTo(ned)]

  # intialize root likelihood
  Xr = Array{Float64,1}(undef,ns)

  # initialize X for tips 
  for i in Base.OneTo(ned)
    edi = ed[i,2]
    edi > ntip && continue
    X[i] = tv[edi]
  end

  # preallocate arrays
  Xp1 = Array{Float64,1}(undef,ns)
  Xp2 = Array{Float64,1}(undef,ns)
  wg  = Array{Float64,1}(undef,ns)
 


  # start ll function for parameters
  function f(p::Array{Float64,1})

    @inbounds begin

      Ets = Et(p)

      # push parameters as the last vector in Ets
      push!(Ets, p)

      # estimate `Gts` according to `Ets` and `p`
      Gts = Gt(Ets)

      # initialize loglik
      ll = 0.0

      for trio in triads
        pr, d1, d2 = trio

        ## first daughter
        X0 = Gts[lbts[d1,2]]\X[d1]
        mul!(Xp1, Gts[lbts[d1,1]], X0)

        ## second daughter
        X0 = Gts[lbts[d2,2]]\X[d2]
        mul!(Xp2, Gts[lbts[d2,1]], X0)

        # combine with speciation rates
        for k in Base.OneTo(ns)
          X[pr][k] = Xp1[k] * Xp2[k] * p[k]
        end

        mdlus = sum(X[pr])
        X[pr] /= mdlus
        ll += log(mdlus)
      end

      # Root conditioning
      for k in Base.OneTo(ns)
        Xr[k] = X[ned-1][k] * X[ned][k] * p[k]
      end

      #estimate weights
      ss = sum(Xr)
      for k in Base.OneTo(ns)
        wg[k] = Xr[k]/ss
      end

      # condition on no extinction
      α = 0.0
      for k in Base.OneTo(ns)
        α += wg[k] * p[k] * (1.0 - Ets[nets-1][k])^2
      end

      Xr /= α

      ll += log(sum(Xr .* pr))
    end

    return ll
  end

  return f
end









# make tree in R
using RCall

ntip = 5

tr, bts = make_ape_tree(ntip, 0.5, 0.1, order= "postorder")

ed  = copy(tr.ed)
el  = copy(tr.el)
ned = size(ed,1)


# make tip values
tv = Dict{Int64,Array{Float64,1}}()
for i in 1:ntip
  st1 = rand(0:1)
  tv[i] = Float64[st1, 1 - st1]
end

p0 = rand(6)
E0 = [0.0,0.0]
k = 2
h = 1

llf = make_loglik(make_fbisseM, fbisseE,
            tv, ed, el, copy(bts), p0, E0, k, h,
            ntip)





p = rand(6) 
llf(p)


















