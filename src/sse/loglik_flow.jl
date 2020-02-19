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
    prepare_ll(ode_make_lik,
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

Prepare EGeoSSE likelihoods for flow algorithm given input model
"""
function prepare_ll(ode_make_lik,
                    ode_ext,
                    cov_mod::NTuple{M,String},
                    tv  ::Dict{Int64,Array{Float64,1}},
                    ed  ::Array{Int64,2},
                    el  ::Array{Float64,1},
                    bts ::Array{Float64,1},
                    p0  ::Array{Float64,1},
                    E0  ::Array{Float64,1},
                    h   ::Int64,
                    ntip::Int64;
                    Eδt::Float64 = 0.01,
                    ti ::Float64 = 0.0) where{M}

  # k areas
  k = length(tv[1])::Int64

  # number of covariates
  ny = size(y,2)

  # number of states
  ns = (2^k - 1)*h

  # make z(t) approximation from discrete data `af!()`
  af! = make_af(x, y, Val(ny))

  # define model
  model = define_mod(cov_mod, k, h, ny)

  # make dictionary with relevant parameters
  pardic = build_par_names(k, h, ny, model)

  # get number of parameters
  npars = length(pardic)

  # find hidden factors for hidden states 
  phid = Int64[] 
  if h > 1
    re = Regex(".*_[1-"*string(h-1)*"]\$")
    for (k,v) in pardic 
      occursin(re, k) && push!(phid, v)
    end
    sort!(phid)
  end

  # create factor parameter vector
  fp = zeros(npars)

  # generate initial parameter values
  p  = fill(0.1,npars)
  βs = h*(k^2 + 2k + h) + 1
  δ  = Float64(length(tv)-1)/sum(el)
  p[βs:end]             .= 0.0                  # set βs
  p[1:(k+1)*h]          .= δ + rand()*δ         # set λs
  p[(k+1)*h+1:h*(2k+1)] .= p[1] - δ             # set μs

  # parameter update
  pupd = 1:npars

  # parameters constraint and fixed to 0
  dcp, dcfp, zp, zfp = 
    set_constraints(constraints, pardic, k, h, ny, model)

  # remove hidden factors from being updated from `p`
  pupd = setdiff(pupd, phid)

  # force pars in zerp to 0
  for i in zp
    p[i] = 0.0
  end

  # remove contraints from being updated
  pupd = setdiff(pupd, values(dcp)) 
  phid = setdiff(phid, values(dcfp)) 

  # remove fixed to zero parameters from being updated
  pupd = setdiff(pupd, zp)
  phid = setdiff(phid, zfp)

  # check if there are hidden factors hold to 0 that also have a forced equality
  for (k,v) in dcfp
    if in(k, zfp) || in(v, zfp) 
      filter!(x -> x ≠ k, phid)
      filter!(x -> x ≠ v, phid)
    end
  end

  # divide between non-negative and negative values
  nnps = filter(x -> βs >  x, pupd)
  nps  = filter(x -> βs <= x, pupd)

  # make hidden factors assigning 
  assign_hidfacs! = 
    make_assign_hidfacs(Val(k), Val(h), Val(ny), Val(model))

  # force same parameter values for constraints
  for wp in keys(dcp)
    while haskey(dcp, wp)
      tp = dcp[wp]
      p[tp] = p[wp]
      wp = tp
    end
  end

  for wp in keys(dcfp)
    while haskey(dcfp, wp)
      tp = dcfp[wp]
      fp[tp] = fp[wp]
      wp = tp
    end
  end

  # assign hidden factors
  assign_hidfacs!(p, fp)

  ntip = length(tv)

  # Estimate extinction at `ts` times
  tf = maximum(bts)
  ts = [ti:Eδt:tf...]

  egeohisse_E = make_egeohisse_E(Val(k), Val(h), Val(ny), Val(model), af!)

  # Make extinction integral
  Et = make_Et(egeohisse_E, p0, E0, ts, ti, tf)

  # estimate first Extinction probabilities at times `ts`
  Ets  = Et(p0)
  nets = length(Ets) + 1

  # Make extinction approximated function
  # ** this make order is crucial **
  afE! = make_af(ts,  Ets, Val(ns))

  # make likelihood integral
  ode_intf = 
    make_egeohisse_M(Val(k), Val(h), Val(ny), Val(model), af!, afE!, nets)

  # push parameters as the last vector in Ets
  push!(Ets, p0)

  # make Gt function
  Gt = make_Gt(ode_intf, Ets, Matrix{Float64}(I, ns, ns), bts, ti, tf)

  # sort branching times
  sort!(bts)

  # add the present `0.0`
  pushfirst!(bts, 0.0)

  nbts = lastindex(bts)

  # estimate `Gts` according to `Ets` and `p0`
  Gts = Gt(Ets)

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
  X = [zeros(ns) for i in Base.OneTo(ned)]

  # create states 
  S = create_states(k, h)
  
  ne    = size(ed,1)
  child = ed[:,2]
  wtp   = findall(child .<= ntip)

  # assign states to terminal branches
  for wi in wtp 
    wig = Set(findall(map(x -> isone(x), tv[child[wi]])))
    X[wi][findall(map(x -> isequal(x.g, wig), S))] .= 1.0
  end

  return Gt, Et, X, triads, lbts, ns, ned, nets
end





"""
    prepare_ll(ode_make_lik,
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

Prepare MuSSE likelihoods for flow algorithm given input model
"""
function prepare_ll(ode_make_lik,
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
  # ** this make order is crucial **
  afE! = make_af(ts,  Ets, Val(ns))

  # make likelihood integral
  ode_intf = ode_make_lik(afE!, nets)

  # push parameters as the last vector in Ets
  push!(Ets, p0)

  # make Gt function
  Gt = make_Gt(ode_intf, Ets, Matrix{Float64}(I, ns, ns), bts, ti, tf)

  # sort branching times
  sort!(bts)

  # add the present `0.0`
  pushfirst!(bts, 0.0)

  nbts = lastindex(bts)

  # estimate `Gts` according to `Ets` and `p0`
  Gts = Gt(Ets)

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

  # initialize X for tips for GeoHiSSE
  X = [Array{Float64,1}(undef,ns) for i in Base.OneTo(ned)]

  # initialize X for tips 
  for i in Base.OneTo(ned)
    edi = ed[i,2]
    edi > ntip && continue
    X[i] = tv[edi]
  end

  return Gt, Et, X, triads, lbts, ns, ned, nets
end





"""
    make_loglik(Gt    ::Function, 
                Et    ::Function, 
                triads::Array{Array{Int64,1},1}, 
                ns    ::Int64, 
                ned   ::Int64, 
                nets  ::Int64)


Make log-likelihood function using the flow algorithm.
"""
function make_loglik(Gt    ::Function, 
                     Et    ::Function,
                     X     ::Array{Array{Float64,1},1}, 
                     triads::Array{Array{Int64,1},1},
                     lbts  ::Array{Int64,2},
                     ns    ::Int64, 
                     ned   ::Int64, 
                     nets  ::Int64)

  # preallocate arrays
  Xp1 = Array{Float64,1}(undef,ns)
  Xp2 = Array{Float64,1}(undef,ns)
  Xr  = Array{Float64,1}(undef,ns)
  wg  = Array{Float64,1}(undef,ns)
  XR  = Array{Float64,1}(undef,ns)

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


"""
HERE --> apply proper node combination
"""

        # combine with speciation rates
        for k in Base.OneTo(ns)
          X[pr][k] = Xp1[k] * Xp2[k] * p[k]
        end

        mdlus  = sum(X[pr])
        X[pr] /= mdlus
        global ll    += log(mdlus)
      end



"""
HERE --> apply proper root conditioning
"""


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

      rmul!(Xr, 1.0/α)

      XR .= Xr .* wg

      ll += log(sum(XR))
    end

    return ll
  end

  return f
end







# # make tree in R
using RCall

ntip = 50

tr, bts = make_ape_tree(ntip, 0.5, 0.0, order = "postorder")

ed  = copy(tr.ed)
el  = copy(tr.el)
ned = size(ed,1)


# make tip values GeoSSE
function sdat() 
  r = rand(0.0:1.0,k)
  while sum(r) < 1.0
    r = rand(0.0:1.0,k)
  end
  return r
end

# create dictionaries
tv = Dict(i => sdat() for i = 1:ntip)



E0 = zeros(6)
k = 2
h = 2





# Gt, Et, X, triads, lbts, ns, ned, nets = 
#   prepare_ll(make_fbisseM, fbisseE, tv, ed, el, copy(bts), p0, E0, k, h, ntip)


# llf = make_loglik(Gt, Et, X, triads, lbts, ns, ned, nets)


# p = rand(6) 

# llf(p)

# @benchmark llf(p)











