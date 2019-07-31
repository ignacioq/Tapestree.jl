#=

Slice sampling for EGeoHiSSE

Ignacio Quintero Mächler

t(-_-t)

November 20 2017

=#




"""
    slice_sampler(tip_val     ::Dict{Int64,Array{Float64,1}},
                  ed          ::Array{Int64,2},
                  el          ::Array{Float64,1},
                  x           ::Array{Float64,1},
                  y           ::Array{Float64},
                  cov_mod     ::String,
                  out_file    ::String,
                  h           ::Int64;
                  constraints ::NTuple{N,String}  = (" ",),
                  niter       ::Int64             = 10_000,
                  nthin       ::Int64             = 10,
                  λpriors     ::Float64           = .1,
                  μpriors     ::Float64           = .1,
                  gpriors     ::Float64           = .1,
                  lpriors     ::Float64           = .1,
                  qpriors     ::Float64           = .1,
                  βpriors     ::NTuple{2,Float64} = (0.0, 5.0),
                  optimal_w   ::Float64           = 0.8,
                  screen_print::Int64             = 5) where {N}

Run slice-sampling Markov Chain for EGeoSSE model.
"""
function slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
                       ed         ::Array{Int64,2},
                       el         ::Array{Float64,1},
                       x          ::Array{Float64,1},
                       y          ::Array{Float64},
                       cov_mod    ::NTuple{M,String},
                       out_file   ::String,
                       h          ::Int64;
                       constraints::NTuple{N,String}  = (" ",),
                       niter      ::Int64             = 10_000,
                       nthin      ::Int64             = 10,
                       λpriors    ::Float64           = .1,
                       μpriors    ::Float64           = .1,
                       gpriors    ::Float64           = .1,
                       lpriors    ::Float64           = .1,
                       qpriors    ::Float64           = .1,
                       βpriors    ::NTuple{2,Float64} = (0.0, 5.0),
                       hpriors    ::Float64           = .1,
                       optimal_w  ::Float64           = 0.8,
                       screen_print::Int64            = 5) where {M,N}

  # k areas
  k = length(tip_val[1])::Int64

  # number of covariates
  ny = size(y,2)

  # make z(t) approximation from discrete data `af!()`
  af! = make_af(x, y, Val(ny))

  # define model
  model = define_mod(cov_mod, k, h, ny)

  # make dictionary with relevant parameters
  pardic = build_par_names(k, h, ny, model)

  # get number of parameters
  npars = length(pardic)

  # find hidden states for hidden states 
  # (Warning: identifies hidden state 10)
  phid = Int64[] 
  for (k,v) in pardic 
    !occursin("0", k) && push!(phid, v)
  end

  # create factor parameter vector
  fp = zeros(npars)

  # generate initial parameter values
  p  = fill(0.1,npars)
  βs = h*(k^2 + 2k + h) + 1
  δ  = log(Float64(length(tip_val)-1))/sum(el)
  p[βs:end]             .= 0.0             # set βs
  p[1:(k+1)*h]          .= δ + rand()*δ  # set λs
  p[(k+1)*h+1:h*(2k+1)] .= p[1] - δ        # set μs

  # parameter update
  pupd = 1:npars

  # parameters constraint and fixed to 0
  conp, zerp = set_constraints(constraints, pardic)

  # force pars in zerp to 0
  p[zerp] .= 0.0

  # remove contraints from being updated
  pupd = setdiff(pupd, values(conp)) 
  phid = setdiff(phid, values(conp)) 
  # remove fixed to zero parameters from being updated
  pupd = setdiff(pupd, zerp)
  phid = setdiff(phid, zerp) 
  # remove hidden factors fro being updated from `p`
  pupd = setdiff(pupd, phid)

  nnps = filter(x -> βs >  x, pupd)
  nps  = filter(x -> βs <= x, pupd)

  # force same parameter values for constraints
  for wp in keys(conp)
    while haskey(conp, wp)
      tp = conp[wp]
      p[tp] = p[wp]
      wp = tp
    end
  end

  # make ODE function
  ode_fun = make_egeohisse(Val(k), Val(h), Val(ny), Val(model), af!)

  # make ODE solver
  ode_solve = make_solver(ode_fun, p, zeros(2*h*(2^k-1)))

  # create likelihood function
  llf = make_llf(tip_val, ed, el, ode_solve, af!, 
    Val(k), Val(h), Val(ny), Val(model))

  # create prior function
  lpf = make_lpf(pupd, phid, 
    λpriors, μpriors, gpriors, lpriors, qpriors, βpriors, hpriors, k, h, model)

  # create posterior functions
  lhf = make_lhf(llf, lpf, conp, Val(k), Val(h), Val(ny), Val(model))

  #=
  run slice sampling
  =#

  # estimate optimal w
  p, fp, w = w_sampler(lhf, p, fp, nnps, nps, phid, npars, optimal_w, screen_print)

  # slice-sampler
  its, hlog, ps = 
    loop_slice_sampler(lhf, p, fp, nnps, nps, phid, w, npars, niter, nthin, screen_print)

  # save samples
  R = hcat(its, hlog, ps)

  # column names
  col_nam = ["Iteration", "Posterior"]

  for (k,v) in sort!(collect(pardic), by = x -> x[2])
    push!(col_nam, k)
  end

  R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

  writedlm(out_file*".log", R)

  return R
end





"""
    slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
                  ed         ::Array{Int64,2},
                  el         ::Array{Float64,1},
                  x          ::Array{Float64,1},
                  y          ::Array{Float64},
                  esse_mod   ::String,
                  out_file   ::String;
                  constraints::NTuple{N,String}  = (" ",),
                  niter      ::Int64             = 10_000,
                  nthin      ::Int64             = 10,
                  λpriors    ::Float64           = .1,
                  μpriors    ::Float64           = .1,
                  qpriors    ::Float64           = .1,
                  βpriors    ::NTuple{2,Float64} = (-1.0,1.0),
                  optimal_w  ::Float64           = 0.8) where {N}

Run slice-sampling Markov Chain for ESSE model.
"""
function slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
                       ed         ::Array{Int64,2},
                       el         ::Array{Float64,1},
                       x          ::Array{Float64,1},
                       y          ::Array{Float64},
                       esse_mod   ::String,
                       out_file   ::String;
                       constraints::NTuple{N,String}  = (" ",),
                       niter      ::Int64             = 10_000,
                       nthin      ::Int64             = 10,
                       λpriors    ::Float64           = .1,
                       μpriors    ::Float64           = .1,
                       qpriors    ::Float64           = .1,
                       βpriors    ::NTuple{2,Float64} = (-1.0,1.0),
                       optimal_w  ::Float64           = 0.8) where {N}

  # k areas
  k = length(tip_val[1])::Int64

  # make z(t) approximation from discrete data
  af = make_approxf(x, y)

  # make specific ode
  mod_ode, npars, pardic, md, ws = define_mod(esse_mod, k, af, size(y,2) > 1)

  # make initial p
  #= 
    TODO -> make smarter initial pars
  =#

  p  = fill(0.2,npars)
  βs = k + k*k + 1
  δ  = log(Float64(length(tip_val)-1))/sum(el)
  p[βs:end]   .= 0.0         # set βs
  p[1:k]      .= δ + rand()  # set λs
  p[(k+1):2k] .= p[1] - δ    # set μs

  # parameter update
  pupd = Base.OneTo(npars)::Base.OneTo{Int64}

  #constraints
  conp = set_constraints(constraints, pardic)::Dict{Int64,Int64}

  # remove contraints for being updated
  pupd = setdiff(pupd, keys(conp)) ::Array{Int64,1}
  nnps = filter(x -> βs >  x, pupd)::Array{Int64,1}
  nps  = filter(x -> βs <= x, pupd)::Array{Int64,1}

  # create likelihood, prior and posterior functions
  llf = make_llf(tip_val, ed, el, mod_ode, af, p, md, ws, sbrlen = sum(el)/10.)
  lpf = make_lpf(λpriors, μpriors, qpriors, βpriors, k, (npars != 2k*k))
  lhf = make_lhf(llf, lpf, conp)

  # set up slice-sampling
  nlogs = fld(niter,nthin)
  its   = zeros(Float64,nlogs)
  h     = zeros(Float64,nlogs)
  ps    = zeros(Float64,nlogs,npars)

  lthin, lit = 0, 0

  # estimate w
  p, w = w_sampler(lhf, p, nnps, nps, npars, optimal_w)::NTuple{2,Array{Float64,1}}

  # start iterations
  prog = Progress(niter, 5, "running slice-sampler....", 20)

  hc::Float64 = lhf(p)

  for it in Base.OneTo(niter) 

    for j in nnps
      S     = (hc - randexp())::Float64
      L, R  = find_nonneg_int(p, j, S, lhf, w[j])::NTuple{2,Float64}
      p, hc = sample_int(p, j, L, R, S, lhf)::Tuple{Array{Float64,1},Float64}
    end

    for j in nps
      S     = (hc - randexp())::Float64
      L, R  = find_real_int(p, j, S, lhf, w[j])::NTuple{2,Float64}
      p, hc = sample_int(p, j, L, R, S, lhf)::Tuple{Array{Float64,1},Float64}
    end

    # log samples
    lthin += 1
    if lthin == nthin
      @inbounds begin
        lit += 1
        setindex!(its, it, lit)
        setindex!(h,   hc, lit)
        setindex!(ps,   p, lit, :)
      end
      lthin = 0
    end

    next!(prog)
  end

  # save samples
  R = hcat(its, h, ps)

  # column names
  col_nam = ["Iteration", "Posterior"]

  # add λ names
  for i in 0:(k-1)
    push!(col_nam, "lamdba$i")
  end

  # add μ names
  for i in 0:(k-1)
    push!(col_nam, "mu$i")
  end

  # add q names
  for j in 0:(k-1), i in 0:(k-1)
    if i == j 
      continue
    end
    push!(col_nam, "q$i$j")
  end

  # add β names
  if in("beta01", keys(pardic))
    for j in 0:(k-1), i in 0:(k-1)
      if i == j 
        continue
      end
      push!(col_nam, "beta$i$j")
    end
  else
    for i in 0:(k-1)
      push!(col_nam, "beta$i")
    end
  end

  R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

  writedlm(out_file*".log", R)

  return R
end





"""
    slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
                  edges      ::Array{Int64,2},
                  edlen      ::Array{Float64,1},
                  out_file   ::String;
                  constraints::String  = NaN,
                  niter      ::Int64   = 10_000,
                  nthin      ::Int64   = 10,
                  model      ::String  = "musse",
                  λpriors    ::Float64 = .1,
                  μpriors    ::Float64 = .1,
                  qpriors    ::Float64 = .1)

Run slice-sampling Markov Chain for MuSSE model.
"""
function slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
                       edges      ::Array{Int64,2},
                       edlen      ::Array{Float64,1},
                       out_file   ::String;
                       constraints::NTuple{N,String} = (" ",),
                       niter      ::Int64            = 10_000,
                       nthin      ::Int64            = 10,
                       λpriors    ::Float64          = .1,
                       μpriors    ::Float64          = .1,
                       qpriors    ::Float64          = .1) where {N}

  # k areas
  k = length(tip_val[1])

  # make ode and define parameters
  mod_ode = make_musse(k)

  npars = k + k*k 
  p = fill(0.2,npars)
  δ = length(tip_val)/sum(edlen)
  p[1:k]      .= δ + rand()    # set λs
  p[(k+1):2k] .= p[1] - δ      # set μs

  pardic = build_par_names(k)

  # parameter update
  pupd = Base.OneTo(npars)

  #constraints
  conp = set_constraints(constraints, pardic)

  # remove contraints for being updated
  pupd = setdiff(pupd, keys(conp))

  # create likelihood, prior and posterior functions
  llf = make_llf(tip_val, edges, edlen, mod_ode, p, sbrlen = sum(edlen))
  lpf = make_lpf(λpriors, μpriors, qpriors, k)
  lhf = make_lhf(llf, lpf, conp)

  # set up slice-sampling
  nlogs = fld(niter,nthin)
  its   = zeros(Float64,nlogs)
  h     = zeros(Float64,nlogs)
  ps    = zeros(Float64,nlogs,npars)

  lthin = 0::Int64
  lit   = 0::Int64

  # cat
  printstyled("running MuSSE model \n", color=:green)

  # estimate w
  p, w = w_sampler(lhf, p, pupd, npars)

  # start iterations
  prog = Progress(niter, 5, "running slice-sampler....", 20)

  hc = lhf(p)::Float64

  for it in Base.OneTo(niter) 

    for j in pupd
      S     = hc - Random.randexp()
      L, R  = find_nonneg_int(p, j, S, lhf, w[j])
      p, hc = sample_int(p, j, L, R, S, lhf)
    end

    # log samples
    lthin += 1
    if lthin == nthin
      @inbounds begin
        lit += 1
        setindex!(its, it, lit)
        setindex!(h,   hc, lit)
        setindex!(ps,   p, lit, :)
      end
      lthin = 0
    end

    next!(prog)
  end

  # save samples
  R = hcat(its, h, ps)

  # column names
  col_nam = ["Iteration", "Posterior"]

  # add λ names
  for i in 0:(k-1)
    push!(col_nam, "lamdba$i")
  end

  # add μ names
  for i in 0:(k-1)
    push!(col_nam, "mu$i")
  end

  # add q names
  for j in 0:(k-1), i in 0:(k-1)
    if i == j 
      continue
    end
    push!(col_nam, "q$i$j")
  end

  R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

  writedlm(out_file*".log", R)

  return R
end




