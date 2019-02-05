#=

Slice sampling for MuSSE

Ignacio Quintero Mächler

t(-_-t)

September 26 2017

=#





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





