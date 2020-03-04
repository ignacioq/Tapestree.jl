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

Prepare EGeoSSE likelihoods for flow algorithm given input model.
"""
function prepare_ll(cov_mod::NTuple{M,String},
                    tv     ::Dict{Int64,Array{Float64,1}},
                    ed     ::Array{Int64,2},
                    el     ::Array{Float64,1},
                    bts    ::Array{Float64,1},
                    E0     ::Array{Float64,1},
                    h      ::Int64,
                    Eδt    ::Float64 = 0.01,
                    ti     ::Float64 = 0.0) where{M}

  # k areas
  k = length(tv[1])::Int64

  # number of covariates
  ny = size(y,2)

  # number of states
  ns = (2^k - 1)*h

  # number of tips
  ntip = length(tv)

  # add root of length 0
  ed = cat(ed, [2*ntip ntip + 1], dims = 1)
  push!(el, 0.0)

  # number of edges
  ned = size(ed,1)

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

  # Estimate extinction at `ts` times
  tf = maximum(bts)
  ts = [ti:Eδt:tf...]

  egeohisse_E = make_egeohisse_E(Val(k), Val(h), Val(ny), Val(model), af!)

  # Make extinction integral
  Et = make_Et(egeohisse_E, p, E0, ts, ti, tf)

  # estimate first Extinction probabilities at times `ts`
  Ets  = Et(p)
  nets = length(Ets) + 1

  # Make extinction approximated function
  # ** this make order is crucial **
  afE! = make_af(ts,  Ets, Val(ns))

  # make likelihood integral
  ode_intf = 
    make_egeohisse_M(Val(k), Val(h), Val(ny), Val(model), af!, afE!, nets)

  # push parameters as the last vector in Ets
  push!(Ets, p)

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

  child = ed[:,2]
  wtp   = findall(child .<= ntip)

  # assign states to terminal branches
  for wi in wtp 
    wig = Set(findall(map(x -> isone(x), tv[child[wi]])))
    X[wi][findall(map(x -> isequal(x.g, wig), S))] .= 1.0
  end

  # make λevent!
  λevent! = make_λevent(h, k, ny, true, model, af!)

  # make root likelihood conditioning
  rootll = make_rootll(h, k, ny, model, af!)

  return Gt, Et, X, triads, lbts, ns, ned, nets, λevent!, rootll
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
    make_λevent(h    ::Int64, 
                k    ::Int64, 
                ny   ::Int64, 
                flow ::Bool,
                model::NTuple{3,Bool},
                af!  ::Function)

Make function for λevent.
"""
function make_λevent(h    ::Int64, 
                     k    ::Int64, 
                     ny   ::Int64, 
                     flow ::Bool,
                     model::NTuple{3,Bool},
                     af!  ::Function)

  λts = Array{Float64,1}(undef,h*k)
  r   = Array{Float64,1}(undef,ny)

  if flow
    λevent! = (t   ::Float64, 
               llik::Array{Array{Float64,1},1},
               ud1 ::Array{Float64,1},
               ud2 ::Array{Float64,1},
               p   ::Array{Float64,1},
               pr  ::Int64) ->
    begin
      λevent_full(t, llik, ud1, ud2, p, pr, λts, r, af!,
        Val(k), Val(h), Val(ny), Val(model))
      return nothing
    end
  else
       # make speciation events and closure
    λevent! = (t   ::Float64, 
               llik::Array{Float64,1},
               ud1 ::Array{Float64,1},
               ud2 ::Array{Float64,1},
               p   ::Array{Float64,1}) ->
    begin
      λevent_full(t, llik, ud1, ud2, p, λts, r, af!,
        Val(k), Val(h), Val(ny), Val(model))
      return nothing
    end
  end

  return λevent!
end





"""
    make_rootll(h    ::Int64, 
                k    ::Int64, 
                ny   ::Int64, 
                model::NTuple{3,Bool},
                af!  ::Function)

Make root conditioning likelihood function.
"""
function make_rootll(h    ::Int64, 
                     k    ::Int64, 
                     ny   ::Int64, 
                     model::NTuple{3,Bool},
                     af!  ::Function)

  λts = Array{Float64,1}(undef,h*k)
  r   = Array{Float64,1}(undef,ny)

  rootll = (t   ::Float64,
            llik::Array{Float64,1},
            extp::Array{Float64,1},
            w   ::Array{Float64,1},
            p   ::Array{Float64,1}) -> 
    begin
      rootll_full(t, llik, extp, w, p, λts, r, af!,
        Val(k), Val(h), Val(ny), Val(model))::Float64
    end


  return rootll
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
  Xp1 = Array{Float64,1}(undef,ns)
  Xp2 = Array{Float64,1}(undef,ns)
  Xr  = Array{Float64,1}(undef,ns)
  wg  = Array{Float64,1}(undef,ns)
  Er  = Array{Float64,1}(undef,ns)

  rtime  = bts[end]
  netsm1 = nets - 1

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

        λtime = bts[lbts[d1,1]]

        λevent!(λtime, X, Xp1, Xp2, p, pr)

        if pr != ned
          mdlus  = sum(X[pr])
          X[pr] /= mdlus
          ll += log(mdlus)
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






