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
                 bts        ::Array{Float64,1},
                 E0         ::Array{Float64,1},
                 h          ::Int64,
                 constraints::NTuple{O,String}) 

Prepare data for **EGeoHiSSE** likelihood calculations.
"""
function prepare_data(cov_mod    ::NTuple{M,String},
                      tv         ::Dict{Int64,Array{Float64,1}},
                      x          ::Array{Float64,1},
                      y          ::Array{Float64,N},
                      ed         ::Array{Int64,2},
                      el         ::Array{Float64,1},
                      E0         ::Array{Float64,1},
                      h          ::Int64,
                      constraints::NTuple{O,String}) where {N,M,O}

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

  # extinction at time 0
  if E0 == [0.0,0.0]
    E0 = zeros(ns)
  end

  # create geographic & hidden states 
  S = create_states(k, h)

  # get absolute times of branches as related to z(t)
  abts = abs_time_branches(el, ed, ntip)

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

  return X, p, fp, trios, ns, ned, pupd, phid, nnps, nps, dcp, dcfp, 
    pardic, k, h, ny, model, af!, assign_hidfacs!, abts, E0

end





