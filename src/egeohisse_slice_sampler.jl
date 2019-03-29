#=

Slice sampling for EGeoHiSSE

Ignacio Quintero Mächler

t(-_-t)

November 20 2017

=#




"""
    slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
                  ed         ::Array{Int64,2},
                  el         ::Array{Float64,1},
                  x          ::Array{Float64,1},
                  y          ::Array{Float64},
                  cov_mod    ::String,
                  out_file   ::String;
                  h          ::Int64             = 2,
                  constraints::NTuple{N,String}  = (" ",),
                  niter      ::Int64             = 10_000,
                  nthin      ::Int64             = 10,
                  λpriors    ::Float64           = .1,
                  μpriors    ::Float64           = .1,
                  gpriors    ::Float64           = .1,
                  lpriors    ::Float64           = .1,
                  qpriors    ::Float64           = .1,
                  βpriors    ::NTuple{2,Float64} = (0.0, 5.0),
                  optimal_w  ::Float64           = 0.8) where {N}

Run slice-sampling Markov Chain for EGeoSSE model.
"""
function slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
                       ed         ::Array{Int64,2},
                       el         ::Array{Float64,1},
                       x          ::Array{Float64,1},
                       y          ::Array{Float64},
                       cov_mod    ::String,
                       out_file   ::String;
                       h          ::Int64             = 2,
                       constraints::NTuple{N,String}  = (" ",),
                       niter      ::Int64             = 10_000,
                       nthin      ::Int64             = 10,
                       λpriors    ::Float64           = .1,
                       μpriors    ::Float64           = .1,
                       gpriors    ::Float64           = .1,
                       lpriors    ::Float64           = .1,
                       qpriors    ::Float64           = .1,
                       βpriors    ::NTuple{2,Float64} = (0.0, 5.0),
                       optimal_w  ::Float64           = 0.8) where {N}

  # k areas
  k  = length(tip_val[1])::Int64

  # number of covariates
  ny = size(y,2)

  # make z(t) approximation from discrete data `af!()`
  make_approxf(x, y)

  # define model
  model = define_mod(cov_mod, k, h, ny)

  # make specific ode
  make_egeohisse(k, h, ny, af!, model, :ode_fun)

  # make dictionary with relevant parameters
  pardic = build_par_names(k, h, ny, model)

  # get number of parameters
  npars = length(pardic)

  # make initial p
  #= 
    TODO -> make smarter initial pars
  =#

  p  = fill(0.1,npars)
  βs = h*(k^2 + 2k + h) + 1
  δ  = log(Float64(length(tip_val)-1))/sum(el)
  p[βs:end]             .= 0.0             # set βs
  p[1:(k+1)*h]          .= δ + rand()*0.1  # set λs
  p[(k+1)*h+1:h*(2k+1)] .= p[1] - δ        # set μs

  # parameter update
  pupd = Base.OneTo(npars)

  # parameters constraints and fixed to 0
  conp, zerp = set_constraints(constraints, pardic)

  # force pars in zerp to 0
  p[zerp] .= 0.0

  # remove contraints from being updated
  pupd = setdiff(pupd, keys(conp)) 
  # remove fixed to zero parameters from being updated
  pupd = setdiff(pupd, zerp)
  nnps = filter(x -> βs >  x, pupd)
  nps  = filter(x -> βs <= x, pupd)


  # create likelihood, prior and posterior functions
  llf   = make_llf(tip_val, ed, el, ode_fun, af!, p, h, model)
  make_lpf(λpriors, μpriors, lpriors, gpriors, qpriors, βpriors, k, h, model[3])  
  lhf   = make_lhf(llf, lpf, conp)

  lpf(p)
  println("could evaluate lpf")

  llf(p)
  println("could evaluate llf")

  lhf(p)
  println("could evaluate lhf")

  # estimate optimal w
  p, w = w_sampler(lhf, p, nnps, nps, npars, optimal_w)

  println(" did w sampling")

  # set up slice-sampling
  nlogs = fld(niter,nthin)
  its   = zeros(Float64,nlogs)
  h     = zeros(Float64,nlogs)
  ps    = zeros(Float64,nlogs,npars)

  lthin, lit = 0, 0

  # preallocate pp
  pp = copy(p)

  # start iterations
  prog = Progress(niter, 5, "running slice-sampler....", 20)

  # intiial posterior
  hc = lhf(p)

  for it in Base.OneTo(niter) 

    for j in nnps
      S     = (hc - Random.randexp())
      L, R  = find_nonneg_int(p, pp, j, S, lhf, w[j])
      p, hc = sample_int(p, pp, j, L, R, S, lhf)
    end

    for j in nps
      S     = (hc - Random.randexp())
      L, R  = find_real_int(p, pp, j, S, lhf, w[j])
      p, hc = sample_int(p, pp, j, L, R, S, lhf)
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

  for (k,v) in sort!(collect(pardic), by = x -> x[2])
    push!(col_nam, k)
  end

  R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

  writedlm(out_file*".log", R)

  return R
end






"""
    define_mod(egeohisse_mod::String,
               k            ::Int64,
               h            ::Int64,
               ny           ::Int64)

Defines EGeoHiSSE model for `k` areas, `h` hidden states and `ny` covariates.
"""
function define_mod(egeohisse_mod::String,
                    k            ::Int64,
                    h            ::Int64,
                    ny           ::Int64)

  if occursin(r"^[s|S][A-za-z]*", egeohisse_mod)         # if speciation
    model   = (true,false,false)
    printstyled("running speciation EGeoHiSSE model with:
  $k single areas 
  $h hidden states 
  $ny covariates \n", 
                 color=:green)
  elseif occursin(r"^[e|E][A-za-z]*", egeohisse_mod)     # if extinction
    model   = (false,true,false)
    printstyled("running extinction EGeoHiSSE model with:
  $k single areas 
  $h hidden states 
  $ny covariates \n",
                 color=:green)
  elseif occursin(r"^[t|T|r|R|q|Q][A-za-z]*", egeohisse_mod) # if transition
    model   = (false,false,true)
    printstyled("running transition EGeoHiSSE model with:
  $k single areas
  $h hidden states
  $ny covariates \n",
                 color=:green)
  else 
    printstyled("running GeoHiSSE model with:
  $k single areas
  $h hidden states 
  0 covariates \n\n",
                color=:green)
  end

  return  model
end



