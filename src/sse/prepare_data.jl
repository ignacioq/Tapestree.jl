#=

Preparation of the data

Ignacio Quintero Mächler

t(-_-t)

Created 16 03 2020

=#



"""
    prepare_data(cov_mod    ::NTuple{M,String},
                 tv         ::Dict{Int64,Array{Float64,1}},
                 x          ::Array{Float64,1},
                 y          ::Array{Float64,N},
                 ed         ::Array{Int64,2},
                 el         ::Array{Float64,1},
                 ρ          ::Array{Float64,1},
                 h          ::Int64,
                 constraints::NTuple{O,String},
                 mvpars     ::NTuple{P,String})

Prepare data for **EGeoHiSSE** likelihood calculations.
"""
function prepare_data(cov_mod    ::NTuple{M,String},
                      tv         ::Dict{Int64,Array{Float64,1}},
                      x          ::Array{Float64,1},
                      y          ::Array{Float64,N},
                      ed         ::Array{Int64,2},
                      el         ::Array{Float64,1},
                      ρ          ::Array{Float64,1},
                      h          ::Int64,
                      constraints::NTuple{O,String},
                      mvpars     ::NTuple{P,String}) where {N,M,O,P}

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
    re   = Regex(".*_[1-"*string(h-1)*"]\$")
    for (k,v) in pardic 
      occursin(re, k) && push!(phid, v)
    end
    sort!(phid)
  end

  # create factor parameter vector
  fp = zeros(npars)

  # generate initial parameter values
  p  = fill(1e-1, npars)
  βs = h*(k^2 + 2k + h) + 1
  δ  = Float64(ntip-1)/sum(el)
  p[βs:end]             .= 0.0                  # set βs
  p[1:(k+1)*h]          .= δ + rand()*δ         # set λs
  p[(k+1)*h+1:h*(2k+1)] .= p[1] - δ             # set μs

  # parameter update
  pupd = 1:npars

  # parameters constraints and fixed to 0
  dcp, dcfp, zp, zfp = 
    set_constraints(constraints, pardic, k, h, ny, model)

  # set multivariate sampled parameters
  mvps0 = set_multivariate(mvpars, pardic)

  # remove hidden factors from being updated from `p`
  pupd = setdiff(pupd, phid)

  # force pars in zerp to 0
  for i in zp
    p[i] = 0.0
  end

  # remove constraints and fixed to zero parameters from being updated
  setdiff!(pupd, values(dcp), zp)
  setdiff!(phid, values(dcfp), zfp)

  # check if there are hidden factors forced to 0 that also have a forced equality
  for (k,v) in dcfp
    if in(k, zfp) || in(v, zfp) 
      filter!(x -> x ≠ k, phid)
      filter!(x -> x ≠ v, phid)
    end
  end

  # divide between non-negative and negative values
  nnps = filter(x -> βs >  x, pupd)
  nps  = filter(x -> βs <= x, pupd)

  # divide multivariate updates into hidden and non-hidden updates
  mvps  = Array{Int64,1}[]
  nngps = Array{Bool,1}[]
  mvhfs = Array{Int64,1}[]
  for (p1, p2) in mvps0

    p1nn = in(p1, nnps)
    p1n  = in(p1, nps)
    p2nn = in(p2, nnps)
    p2n  = in(p2, nps)

    if p1nn && p2nn
      setdiff!(nnps, p1, p2)
      push!(mvps, [p1,p2])
      push!(nngps, [p1nn,p2nn])
    end

    if p1n && p2nn
      setdiff!(nps,  p1)
      setdiff!(nnps, p2)
      push!(mvps, [p1,p2])
      push!(nngps, [p1nn, p2nn])
    end

    if p1nn && p2n
      setdiff!(nnps, p1)
      setdiff!(nps,  p2)
      push!(mvps, [p1,p2])
      push!(nngps, [p1nn, p2nn])
    end

    if p1n && p2n
      setdiff!(nps, p1, p2)
      push!(mvps, [p1,p2])
      push!(nngps, [p1nn, p2nn])
    end

    if in(p1, phid) && in(p2, phid)
      setdiff!(phid, p1, p2)
      push!(mvhfs, [p1,p2])
    end
  end

  # make hidden factors assigning 
  assign_hidfacs! = make_assign_hidfacs(Val(k), Val(h), Val(ny), Val(model))

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

  # extinction at time 0 with sampling fraction `ρ_i`
  if isone(lastindex(ρ))
    E0 = fill(1.0 - ρ[1], ns)
  elseif lastindex(ρ) == ns
    E0 = ones(Float64, ns)
    for i in Base.OneTo(ns)
      E0[i] -= ρ[i] 
    end
  else
    @error "Length of sampling fraction vector $(lastindex(ρ)) does not match length of the number of states $ns"
  end

  # create geographic & hidden states 
  S = create_states(k, h)

  # get absolute times of branches as related to z(t)
  abts = abs_time_branches(el, ed, ntip)

  # get branching times
  bts = map(x -> round(x; digits = 9), abts[:,1])
  bts = unique!(bts)

  sort!(bts)

  # make trios
  trios = maketriads(ed, rev = true)

  # preallocate tip likelihoods
  X = [zeros(ns) for i in Base.OneTo(ned)]

  child = ed[:,2]
  wtp   = findall(child .<= ntip)

  # assign states to terminal branches
  for wi in wtp 
    wig = Set(findall(map(x -> isone(x), tv[child[wi]])))
    X[wi][findall(map(x -> isequal(x.g, wig), S))] .= 1.0
  end

  return X, p, fp, trios, ns, ned, pupd, phid, nnps, nps, mvps, nngps, mvhfs,
    dcp, dcfp, pardic, k, h, ny, model, af!, assign_hidfacs!, abts, bts, E0

end





