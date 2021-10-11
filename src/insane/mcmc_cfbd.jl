#=

constant fossilized birth-death MCMC wrapper

Jérémy Andréoletti
Adapted from fossilized birth-death MCMC wrapper by Ignacio Quintero Mächler

v(°-°v)

Created 07 10 2021
=#




"""
    insane_cfbd(tree        ::sTbd, 
                out_file    ::String;
                augmentation::String            = "fs",
                λprior      ::Float64           = 0.1,
                μprior      ::Float64           = 0.1,
                ψprior      ::Float64           = 0.1,
                niter       ::Int64             = 1_000,
                nthin       ::Int64             = 10,
                nburn       ::Int64             = 200,
                tune_int    ::Int64             = 100,
                ϵi          ::Float64           = 0.4,
                λi          ::Float64           = NaN,
                μi          ::Float64           = NaN,
                ψi          ::Float64           = NaN,
                λtni        ::Float64           = 1.0,
                μtni        ::Float64           = 1.0,
                ψtni        ::Float64           = 1.0,
                obj_ar      ::Float64           = 0.4,
                pupdp       ::NTuple{3,Float64} = (0.2,0.2,0.5),
                prints      ::Int64             = 5)

Run insane for constant fossilized pure-birth.
"""
function insane_cfbd(tree        ::sTbd, 
                     out_file    ::String;
                     augmentation::String            = "fs",
                     λprior      ::Float64           = 0.1,
                     μprior      ::Float64           = 0.1,
                     ψprior      ::Float64           = 0.1,
                     niter       ::Int64             = 1_000,
                     nthin       ::Int64             = 10,
                     nburn       ::Int64             = 200,
                     tune_int    ::Int64             = 100,
                     ϵi          ::Float64           = 0.4,
                     λi          ::Float64           = NaN,
                     μi          ::Float64           = NaN,
                     ψi          ::Float64           = NaN,
                     λtni        ::Float64           = 1.0,
                     μtni        ::Float64           = 1.0,
                     ψtni        ::Float64           = 1.0,
                     obj_ar      ::Float64           = 0.4,
                     pupdp       ::NTuple{3,Float64} = (0.2,0.2,0.2),
                     prints      ::Int64             = 5)

  # forward simulation
  if occursin(r"^[f|F][A-za-z]*", augmentation)
    R, tree = insane_cfbd_fs(tree, out_file, λprior, μprior, ψprior, 
      niter, nthin, nburn, tune_int, ϵi, λi, μi, ψi, λtni, μtni, ψtni, 
      obj_ar, pupdp, prints)

  # graft/prune
  elseif occursin(r"^[g|G][A-za-z]*", augmentation)
    R, tree = insane_cfbd_gp(tree, out_file, λprior, μprior, ψprior, 
      niter, nthin, nburn, tune_int, ϵi, λi, μi, ψi, λtni, μtni, ψtni, 
      obj_ar, pupdp, prints)
  else
    @error string(augmentation," does not match `fs` or `gp`")
  end

  return R, tree
end





