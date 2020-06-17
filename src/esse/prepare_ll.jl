#=

Log-likelihood preparation given the data

Ignacio Quintero Mächler

t(-_-t)

Created 05 03 2020

=#





"""
    prepare_ll(X    ::Array{Array{Float64,1},1},
               p    ::Array{Float64,1},
               E0   ::Array{Float64,1},
               ns   ::Int64,
               k    ::Int64,
               h    ::Int64,
               ny   ::Int64,
               model::NTuple{3, Bool},
               abts ::Array{Float64,2},
               af!  ::Function)

Prepare **EGeoHiSSE** likelihoods using the **pruning** algorithm 
given input data.
"""
function prepare_ll(X    ::Array{Array{Float64,1},1},
                    p    ::Array{Float64,1},
                    E0   ::Array{Float64,1},
                    ns   ::Int64,
                    k    ::Int64,
                    h    ::Int64,
                    ny   ::Int64,
                    model::NTuple{3, Bool},
                    abts ::Array{Float64,2},
                    af!  ::Function)

  # add extinction to likelihood vector
  for i in Base.OneTo(lastindex(X))
    append!(X[i], E0)
  end

  # make them vectors for indexing efficiency
  abts1 = abts[:,1]
  abts2 = abts[:,2]

  # make speciation events and closure
  λevent! = make_λevent(h, k, ny, false, model, af!)

  # make root likelihood conditioning
  rootll = make_rootll(h, k, ny, model, af!)

  # make ODE function
  ode_fun = make_egeohisse(Val(k), Val(h), Val(ny), Val(model), af!)

  # make integral problem
  prob = ODEProblem(ode_fun, zeros(2*h*(2^k-1)), (0.0,1.0), p)

  int = init(prob, Tsit5(),
             save_everystep  = false, 
             calck           = false,
             force_dtmin     = true,
             save_start      = false,
             initialize_save = false,
             maxiters        = 100_000_000,
             verbose         = false)

  return X, int, λevent!, rootll, abts1, abts2
end





"""
    prepare_ll(cov_mod::NTuple{M,String},
               tv     ::Dict{Int64,Array{Float64,1}},
               ed     ::Array{Int64,2},
               el     ::Array{Float64,1},
               bts    ::Array{Float64,1},
               E0     ::Array{Float64,1},
               h      ::Int64;
               Eδt    ::Float64 = 0.01,
               ti     ::Float64 = 0.0) where{M}

Prepare **EGeoHiSSE** likelihoods using the **flow** algorithm given input data.
"""
function prepare_ll(p    ::Array{Float64,1},
                    bts  ::Array{Float64,1},
                    E0   ::Array{Float64,1},
                    k    ::Int64,
                    h    ::Int64,
                    ny   ::Int64,
                    ns   ::Int64,
                    ned  ::Int64,
                    model::NTuple{3, Bool},
                    Eδt  ::Float64,
                    ti   ::Float64,
                    abts ::Array{Float64,2},
                    af!  ::Function)

  # Estimate extinction at `ts` times
  tf = maximum(bts)
  ts = [ti:Eδt:tf...]

  # Make extinction integral
  egeohisse_E = make_egeohisse_E(Val(k), Val(h), Val(ny), Val(model), af!)

  # Make extinction function at times `ts`
  Et = make_Et(egeohisse_E, p, E0, ts, ti, tf)

  # estimate first Extinction probabilities at times `ts`
  Ets  = Et(p)
  nets = length(Ets) + 1

  # Make extinction approximated function
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

  # connect branches with branching times
  lbts = Array{Int64,2}(undef,size(abts))

  # preallocate absolute minimum matrix of differences
  minM = Array{Float64,1}(undef, nbts)
  # find links
  for i in Base.OneTo(2ned)
    for j in Base.OneTo(nbts)
      minM[j] = abs(abts[i] - bts[j])
    end
    lbts[i] = argmin(minM)
  end

  # make λevent!
  λevent! = make_λevent(h, k, ny, true, model, af!)

  # make root likelihood conditioning
  rootll = make_rootll(h, k, ny, model, af!)

  return Gt, Et, lbts, nets, λevent!, rootll
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

Prepare **MuSSE** likelihoods for flow algorithm given input model
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




