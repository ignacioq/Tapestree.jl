#=

Biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#



"""
    compete(tree_file::String, data_file::String)

Run Compete. Wrapper around all functions.
"""
function compete(tree_file::String,
                 data_file::String,
                 out_file ::String;
                 min_dt   ::Float64           = 0.01,
                 niter    ::Int64             = 50_000,
                 nburn    ::Int64             = 50_000,
                 nthin    ::Int64             = 500,
                 saveXY   ::Tuple{Bool,Int64} = (false, 1_000),
                 ωxprior  ::NTuple{2,Float64} = (0.,10.),
                 ω1prior  ::NTuple{2,Float64} = (0.,10.),
                 ω0prior  ::NTuple{2,Float64} = (0.,10.),
                 σ²prior  ::Float64           = 1e-1,
                 λprior   ::Float64           = 1e-1,
                 weight   ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),
                 λ1i      ::Float64           = 1.0,
                 λ0i      ::Float64           = 0.4,
                 ωxi      ::Float64           = 0.0,
                 ω1i      ::Float64           = 0.0,
                 ω0i      ::Float64           = 0.0,
                 fix_ωx   ::Bool              = false,
                 fix_ω1   ::Bool              = false,
                 fix_ω0   ::Bool              = false,
                 bbprop   ::Bool              = true,
                 delim    ::Char              = '\t',
                 eol      ::Char              = '\r')

  tip_values, tip_areas, tree, bts = 
    read_data(tree_file, data_file)

  X, Y, B, ncoup, δt, tree, si = 
    initialize_data(tip_values, tip_areas, min_dt, tree, bts)

  R = compete_mcmc(X, Y, ncoup, δt, tree.ed, tree.el, B,
                   niter    = niter,
                   nthin    = nthin,
                   nburn    = nburn,
                   saveXY   = saveXY,
                   ωxprior  = ωxprior,
                   ω1prior  = ω1prior,
                   ω0prior  = ω0prior,
                   σ²prior  = σ²prior,
                   λprior   = λprior,
                   out_file = out_file,
                   weight   = weight,
                   λ1i      = λ1i,
                   λ0i      = λ0i,
                   ωxi      = ωxi,
                   ω1i      = ω1i,
                   ω0i      = ω0i,
                   σ²i      = si,
                   stbrl    = maximum(bts),
                   fix_ωx   = fix_ωx,
                   fix_ω1   = fix_ω1,
                   fix_ω0   = fix_ω0)

  return R

end




"""
    compete_for_sims(tree_file::String, data_file::String)

Run Compete from simulations. Wrapper around all functions.
"""
function compete_for_sims(tip_values::Dict{Int64,Float64}, 
                          tip_areas ::Dict{Int64,Array{Int64,1}},
                          tree      ::rtree, 
                          bts       ::Array{Float64,1},
                          out_file  ::String;
                          min_dt    ::Float64           = 0.01,
                          niter     ::Int64             = 500_000,
                          nburn     ::Int64             = 500_000,
                          nthin     ::Int64             = 1_000,
                          saveXY    ::Tuple{Bool,Int64} = (false, 1_000),
                          ωxprior   ::NTuple{2,Float64} = (0.,10.),
                          ω1prior   ::NTuple{2,Float64} = (0.,10.),
                          ω0prior   ::NTuple{2,Float64} = (0.,10.),
                          σ²prior   ::Float64           = 1e-1,
                          λprior    ::Float64           = 1e-1,
                          weight    ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),
                          λ1i       ::Float64           = 1.0,
                          λ0i       ::Float64           = 0.4,                          ωxi       ::Float64           = 0.,
                          ω1i       ::Float64           = 0.,
                          ω0i       ::Float64           = 0.,
                          fix_ωx    ::Bool              = false,
                          fix_ω1    ::Bool              = false,
                          fix_ω0    ::Bool              = false,
                          bbprop    ::Bool              = true,
                          delim     ::Char              = '\t',
                          eol       ::Char              = '\r')

  X, Y, B, ncoup, δt, tree, si = 
    initialize_data(tip_values, tip_areas, min_dt, tree, bts)

  R = compete_mcmc(X, Y, ncoup, δt, 
                   deepcopy(tree.ed), 
                   deepcopy(tree.el), 
                   B,
                   niter     = niter,
                   nthin     = nthin,
                   nburn     = nburn,
                   saveXY    = saveXY,
                   ωxprior   = ωxprior,
                   ω1prior   = ω1prior,
                   ω0prior   = ω0prior,
                   σ²prior   = σ²prior,
                   λprior    = λprior,
                   out_file  = out_file,
                   weight    = weight,
                   λ1i      = λ1i,
                   λ0i      = λ0i,
                   ωxi      = ωxi,
                   ω1i      = ω1i,
                   ω0i      = ω0i,
                   σ²i      = si,
                   stbrl    = maximum(bts),
                   fix_ω1   = fix_ω1,
                   fix_ω0   = fix_ω0,
                   fix_ωx   = fix_ωx)

  return R

end

