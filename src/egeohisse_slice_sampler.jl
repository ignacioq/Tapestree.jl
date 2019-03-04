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

Run slice-sampling Markov Chain for MuSSE model.
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
    define_mod(egeohisse_mod::String,
               k            ::Int64,
               h            ::Int64,
               af           ::Function,
               md           ::Bool)

Defines EGeoHiSSE model for `k` areas and `h` hidden states.
"""
function define_mod(egeohisse_mod::String,
                    k            ::Int64,
                    h            ::Int64,
                    af           ::Function,
                    md           ::Bool)

  if occursin(r"^[s|S][A-za-z]*", esse_mod)         # if speciation
    mod_ode = md ? make_esse_s(k, af, md) : make_esse_s(k, af)
    npars   = 2k + k*k
    pardic  = build_par_names(k, h, (true, false, false))
    ws      = true
    printstyled("running speciation Geographical ESSE model with $k 
                 single areas and $h hidden states \n", color=:green)

  elseif occursin(r"^[e|E][A-za-z]*", esse_mod)     # if extinction
    mod_ode = md ? make_esse_e(k, af, md) : make_esse_e(k, af)
    npars   = (2k + k*k)::Int64
    pardic  = build_par_names(k, h, (false, true, false))
    ws      = false
    printstyled("running extinction Geographical ESSE model with $k 
                 single areas and $h hidden states \n", color=:green)

  elseif occursin(r"^[t|T|r|R|q|Q][A-za-z]*", esse_mod) # if transition
    mod_ode = md ? make_esse_q(k, af, md) : make_esse_q(k, af)
    npars   = 2k*k::Int64
    pardic  = build_par_names(k, h, (false, false, true))
    ws      = false
    printstyled("running transition Geographical ESSE model with $k 
                 single areas and $h hidden states \n", color=:green)

  else 
    error("esse_mod does not match any of the alternatives: 
          speciation, extinction or transition")
  end

  return mod_ode, npars, pardic, md, ws

end



