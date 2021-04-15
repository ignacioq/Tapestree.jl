#=

constant birth-death MCMC wrapper

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cbd(tree    ::sTbd, 
               out_file::String;
               λprior  ::Float64           = 0.1,
               μprior  ::Float64           = 0.1,
               niter   ::Int64             = 1_000,
               nthin   ::Int64             = 10,
               nburn   ::Int64             = 200,
               tune_int::Int64             = 100,
               ϵi      ::Float64           = 0.4,
               λi      ::Float64           = NaN,
               μi      ::Float64           = NaN,
               λtni    ::Float64           = 1.0,
               μtni    ::Float64           = 1.0,
               obj_ar  ::Float64           = 0.4,
               pupdp   ::NTuple{4,Float64} = (0.4,0.4,0.1,0.1),
               prints  ::Int64              = 5)

Run insane for constant pure-birth.
"""
function insane_cbd(tree        ::sTbd, 
                    out_file    ::String;
                    augmentation::String            = "fs",
                    λprior      ::Float64           = 0.1,
                    μprior      ::Float64           = 0.1,
                    niter       ::Int64             = 1_000,
                    nthin       ::Int64             = 10,
                    nburn       ::Int64             = 200,
                    tune_int    ::Int64             = 100,
                    ϵi          ::Float64           = 0.4,
                    λi          ::Float64           = NaN,
                    μi          ::Float64           = NaN,
                    λtni        ::Float64           = 1.0,
                    μtni        ::Float64           = 1.0,
                    ntry        ::Int64             = 1,
                    obj_ar      ::Float64           = 0.4,
                    pupdp       ::Array{Float64,1}  = [0.4,0.4,0.1,0.1],
                    prints      ::Int64             = 5)

  # forward simulation
  if occursin(r"^[f|F][A-za-z]*", augmentation)
    R, tree = insane_cbd_fs(tree, out_file, 
      λprior = λprior, 
      μprior = μprior, 
      niter = niter, 
      nthin = nthin, 
      nburn = nburn, 
      tune_int = tune_int, 
      ϵi = ϵi, 
      λi = λi, 
      μi = μi, 
      ntry = ntry, 
      λtni = λtni, 
      μtni = μtni, 
      obj_ar = obj_ar, 
      pupdp = pupdp, 
      prints = prints)

  # graft/prune
  elseif occursin(r"^[g|G][A-za-z]*", augmentation)
    R, tree = insane_cbd_gp(tree, out_file, 
      λprior = λprior, 
      μprior = μprior, 
      niter = niter, 
      nthin = nthin, 
      nburn = nburn, 
      tune_int = tune_int, 
      ϵi = ϵi, 
      λi = λi, 
      μi = μi, 
      λtni = λtni, 
      μtni = μtni, 
      obj_ar = obj_ar, 
      pupdp = pupdp, 
      prints = prints)
  else
    @error string(augmentation," does not match `fs` or `gp`")
  end

  return R, tree
end





