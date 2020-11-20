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
                 ncch       ::Int64,
                 constraints::NTuple{O,String},
                 mvpars     ::NTuple{P,String}) where {N,M,O,P}

Prepare data for **ESSE.g** likelihood calculations.
"""
function prepare_data(cov_mod    ::NTuple{M,String},
                      tv         ::Dict{Int64,Array{Float64,1}},
                      x          ::Array{Float64,1},
                      y          ::Array{Float64,N},
                      ed         ::Array{Int64,2},
                      el         ::Array{Float64,1},
                      ρ          ::Array{Float64,1},
                      h          ::Int64,
                      ncch       ::Int64,
                      constraints::NTuple{O,String},
                      mvpars     ::NTuple{P,String},
                      parallel   ::Bool) where {N,M,O,P}

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
  af! = make_af(x, y, Val{ny})

  # define model
  model = define_mod(cov_mod, k, h, ny)

  # make dictionary with relevant parameters
  pardic = build_par_names(k, h, ny, model)

  # get number of parameters
  npars = length(pardic)

  # find λ hidden factors for hidden states 
  phid = Int64[] 
  if h > 1
    re   = Regex("^lambda_.*_[1-"*string(h-1)*"]\$")
    for (k,v) in pardic 
      occursin(re, k) && push!(phid, v)
    end
    sort!(phid)
  end

  # create factor parameter vector
  fp = zeros(k > 1 ? (k+1)*h : h)

  # generate initial parameter values
  p  = fill(1e-1, npars)

  # Pure-birth MLE
  δ  = Float64(ntip-1)/sum(el)

  if isone(k)
    βs = (2*h + h*(h-1)) + 1
    p[βs:end]      .= 0.0           # set βs
    p[1:h]         .= δ + 0.1*rand()*δ  # set λs
    p[(h+1):(2*h)] .= p[1] - δ      # set μs
  else
    βs = h*(k^2 + 2k + h) + 1
    p[βs:end]             .= 0.0                  # set βs
    p[1:(k+1)*h]          .= δ + 0.1*rand()*δ         # set λs
    p[(k+1)*h+1:h*(2k+1)] .= p[1] - δ             # set μs
  end

  # make vector for ncch
  p  = [copy(p)  for i in Base.OneTo(ncch)]
  fp = [copy(fp) for i in Base.OneTo(ncch)]

  # parameter update
  pupd = [1:npars...]

  # parameters constraints and fixed to 0
  dcp, zp, zfp = 
    set_constraints(constraints, pardic, k, h, ny, model)

  # set multivariate sampled parameters
  mvps0 = set_multivariate(mvpars, pardic)

  # remove hidden factors from being updated from `p`
  setdiff!(pupd, phid)

  # force pars in zerp to 0
  for c in Base.OneTo(ncch), i in zp
    p[c][i] = 0.0
  end

  # remove constraints and fixed to zero parameters from being updated
  setdiff!(pupd, values(dcp), zp)
  setdiff!(phid, zfp)

  # divide between non-negative and negative values
  nnps = filter(x -> βs >  x, pupd)
  nps  = filter(x -> βs <= x, pupd)

  # divide multivariate updates for λ into hidden and non-hidden updates
  mvps  = Array{Int64,1}[]
  mvhfs = Array{Int64,1}[]
  nngps = Array{Bool,1}[]
  hfgps = Array{Bool,1}[]
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

    if in(p1, phid)
      setdiff!(phid, p1)
      push!(mvhfs, [p1,p2])
      if in(p2, phid)
        setdiff!(phid, p2)
        push!(hfgps, [true, true])
      else
        p2nn ? setdiff!(nnps, p2) : setdiff!(nps, p2)
        push!(hfgps, [true, false])
      end
      setdiff!(phid, p1, p2)
    elseif in(p2, phid)
      p1nn ? setdiff!(nnps, p1) : setdiff!(nps, p1)
      setdiff!(phid, p2)
      push!(mvhfs, [p1,p2])
      push!(hfgps, [false, true])
    end
  end

  @debug " nps   = $nps 
          nnps  = $nnps 
          mvps  = $mvps 
          nngps = $nngps 
          mvhfs = $mvhfs 
          hfgps = $hfgps
          phid  = $phid"

  # make hidden factors assigning 
  assign_hidfacs! = make_assign_hidfacs(Val{k}, Val{h})

  # force same parameter values for constraints
  for c in Base.OneTo(ncch)
    for wp in keys(dcp)
      while haskey(dcp, wp)
        tp = dcp[wp]
        p[c][tp] = p[c][wp]
        wp = tp
      end
    end

    # assign hidden factors
    assign_hidfacs!(p[c], fp[c])
  end

  # if parallel
  if parallel
    p  = distribute(p)
    fp = distribute(fp)
  end

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

  return X, p, fp, trios, ns, ned, pupd, phid, nnps, nps, 
    mvps, nngps, mvhfs, hfgps, dcp, pardic, k, h, ny, model, 
    af!, assign_hidfacs!, abts, bts, E0

end





