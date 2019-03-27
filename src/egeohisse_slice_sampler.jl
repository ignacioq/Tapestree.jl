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
                  esse_mod   ::String,
                  out_file   ::String;
                  constraints::NTuple{N,String} = (" ",),
                  niter      ::Int64            = 10_000,
                  nthin      ::Int64            = 10,
                  λpriors    ::Float64          = .1,
                  μpriors    ::Float64          = .1,
                  qpriors    ::Float64          = .1)

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
                       qpriors    ::Float64           = .1,
                       βpriors    ::NTuple{2,Float64} = (-1.0,1.0),
                       optimal_w  ::Float64           = 0.8) where {N}

  # k areas
  k  = length(tip_val[1])::Int64

  # number of covariates
  ny = size(y,2)

  # make z(t) approximation from discrete data `af!()`
  make_approxf(x, y)

  # make specific ode
  ode_fun, npars, pardic, model = 
    define_mod(cov_mod, k, h, ny, af!)

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
  llf = make_llf(tip_val, ed, el, ode_fun, af!, p, h, model)


  # define priors

  


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
    define_mod(egeohisse_mod::String,
               k            ::Int64,
               h            ::Int64,
               ny           ::Int64,
               af!          ::Function)

Defines EGeoHiSSE model for `k` areas, `h` hidden states and `ny` covariates.
"""
function define_mod(egeohisse_mod::String,
                    k            ::Int64,
                    h            ::Int64,
                    ny           ::Int64,
                    af!          ::Function)

  if occursin(r"^[s|S][A-za-z]*", egeohisse_mod)         # if speciation
    model   = (true,false,false)
    printstyled("running speciation EGeoHiSSE model with $k single areas, $h hidden states and $ny covariates \n", 
                 color=:green)
  elseif occursin(r"^[e|E][A-za-z]*", egeohisse_mod)     # if extinction
    model   = (false,true,false)
    printstyled("running extinction EGeoHiSSE model with $k single areas, $h hidden states and $ny covariates \n", 
                 color=:green)
  elseif occursin(r"^[t|T|r|R|q|Q][A-za-z]*", egeohisse_mod) # if transition
    model   = (false,false,true)
    printstyled("running transition EGeoHiSSE model with $k single areas, $h hidden states and $ny covariates \n", 
                 color=:green)
  else 
    printstyled("running GeoHiSSE model with $k single areas and 
                $h hidden states but no covariates \n", 
                color=:green)
  end

  mod_ode = make_egeohisse(k, h, ny, af!, model, :mod_ode)
  pardic  = build_par_names(k, h, ny, model)
  npars   = length(pardic)

  return mod_ode, npars, pardic, model
end



