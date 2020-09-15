#=

pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    insane_cpb(tree    ::iTpb, 
               out_file::String;
               λprior  ::Float64 = 0.1,
               niter   ::Int64   = 1_000,
               nthin   ::Int64   = 10,
               nburn   ::Int64   = 200,
               tune_int::Int64   = 100,
               λtni    ::Float64 = 1.0,
               obj_ar  ::Float64 = 0.4)

Run insane for constant pure-birth.
"""
function insane_cpb(tree    ::iTpb, 
                    out_file::String;
                    λprior  ::Float64 = 0.1,
                    niter   ::Int64   = 1_000,
                    nthin   ::Int64   = 10,
                    nburn   ::Int64   = 200,
                    tune_int::Int64   = 100,
                    λtni    ::Float64 = 1.0,
                    obj_ar  ::Float64 = 0.4)

  # tree characters
  tl = treelength(tree)
  nt = sntn(tree)

  scalef = makescalef(obj_ar)

  # adaptive phase
  llc, prc, λc, λtn = 
    mcmc_burn_cpb(tree, tl, nt, tune_int, λprior, nburn, λtni, scalef)

  # mcmc
  R =  mcmc_cpb(tree, llc, prc, λc, λprior, niter, nthin, λtn)

  pardic = Dict(("lambda" => 1))

  write_ssr(R, pardic, out_file)

  return R
end




"""
    mcmc_burn_cpb(tree    ::iTpb, 
                  tl      ::Float64,
                  nt      ::Int64,
                  tune_int::Int64,
                  λprior  ::Float64,
                  nburn   ::Int64,
                  λtni    ::Float64, 
                  scalef  ::Function)

MCMC chain for constant pure-birth.
"""
function mcmc_burn_cpb(tree    ::iTpb, 
                       tl      ::Float64,
                       nt      ::Int64,
                       tune_int::Int64,
                       λprior  ::Float64,
                       nburn   ::Int64,
                       λtni    ::Float64, 
                       scalef  ::Function)

  # initialize acceptance log
  ltn = 0
  lup = 0.0
  lac = 0.0
  λtn = λtni

  # starting parameters
  λc  = Float64(nt-1)/tl
  llc = llik_cpb(tree, λc)
  prc = logdexp(λc, λprior)

  for it in Base.OneTo(nburn)

    # parameter proposals
    λp = mulupt(λc, λtn)::Float64

    # one could make a ratio likelihood function
    llp = llik_cpb(tree, λp)
    prr = llrdexp_x(λp, λc, λprior)

    if -randexp() < (llp - llc + prr + log(λp/λc))
      llc  = llp::Float64
      prc += prr::Float64
      λc   = λp::Float64
      lac += 1.0
    end

    # log tuning parameters
    ltn += 1
    lup += 1.0
    if ltn == tune_int
      λtn = scalef(λtn,lac/lup)
      ltn = 0
    end

  end

  return llc, prc, λc, λtn
end




"""
    mcmc_cpb(tree  ::iTpb,
             llc   ::Float64,
             prc   ::Float64,
             λc    ::Float64,
             λprior::Float64,
             niter ::Int64,
             nthin ::Int64,
             λtn   ::Float64)

MCMC chain for constant pure-birth.
"""
function mcmc_cpb(tree  ::iTpb,
                  llc   ::Float64,
                  prc   ::Float64,
                  λc    ::Float64,
                  λprior::Float64,
                  niter ::Int64,
                  nthin ::Int64,
                  λtn   ::Float64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  R = Array{Float64,2}(undef, nlogs, 4)

  for it in Base.OneTo(niter)

    # parameter proposals
    λp = mulupt(λc, rand() < 0.3 ? λtn : 4.0*λtn)::Float64

    # one could make a ratio likelihood function
    llp = llik_cpb(tree, λp)
    prr = llrdexp_x(λp, λc, λprior)

    if -randexp() < (llp - llc + prr + log(λp/λc))
      llc  = llp::Float64
      prc += prr::Float64
      λc   = λp::Float64
    end

    # log parameters
    lthin += 1
    if lthin == nthin
      lit += 1
      @inbounds begin
        R[lit,1] = Float64(lit)
        R[lit,2] = llc
        R[lit,3] = prc
        R[lit,4] = λc
      end
      lthin = 0
    end

  end

  return R
end






"""
  write_ssr(R       ::Array{Float64,2}, 
            pardic  ::Dict{String,Int64},
            out_file::String)

Write the samples from an MC sampler data frame 
given a Dictionary of parameters.
"""
function write_ssr(R       ::Array{Float64,2}, 
                   pardic  ::Dict{String,Int64},
                   out_file::String)

  # column names
  col_nam = ["Iteration", "Likelihood", "Prior"]

  for (k,v) in sort!(collect(pardic), by = x -> x[2])
    push!(col_nam, k)
  end

  R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

  writedlm(out_file*".log", R)
end

