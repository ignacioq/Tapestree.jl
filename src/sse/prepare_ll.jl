#=

Log-likelihood preparation given the data

Ignacio Quintero Mächler

t(-_-t)

Created 05 03 2020

=#





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
function prepare_ll(cov_mod    ::NTuple{M,String},
                    tv         ::Dict{Int64,Array{Float64,1}},
                    x          ::Array{Float64,1},
                    y          ::Array{Float64,N},
                    ed         ::Array{Int64,2},
                    el         ::Array{Float64,1},
                    bts        ::Array{Float64,1},
                    E0         ::Array{Float64,1},
                    h          ::Int64;
                    constraints::NTuple{O,String} = (" ",),
                    Eδt        ::Float64          = 0.01,
                    ti         ::Float64          = 0.0) where {N,M,O}

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
    for i in Base.OneTo(2ned)
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

  return Gt, Et, X, p, fp, triads, lbts, ns, ned, nets, pupd, phid, nnps, nps,
    dcp, dcfp, pardic, k, h, ny, model, λevent!, rootll, assign_hidfacs!
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

Prepare **EGeoHiSSE** likelihoods using the **pruning** algorithm 
given input data.
"""
function prepare_ll(cov_mod    ::NTuple{M,String},
                    tv         ::Dict{Int64,Array{Float64,1}},
                    x          ::Array{Float64,1},
                    y          ::Array{Float64,N},
                    ed         ::Array{Int64,2},
                    el         ::Array{Float64,1},
                    bts        ::Array{Float64,1},
                    E0         ::Array{Float64,1},
                    h          ::Int64;
                    constraints::NTuple{O,String} = (" ",))







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




