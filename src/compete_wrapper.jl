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
                 min_dt   ::Float64                = 0.01,
                 niter    ::Int64                  = 500_000,
                 nburn    ::Int64                  = 500_000,
                 nthin    ::Int64                  = 1_000,
                 ωxprior  ::Tuple{Float64,Float64} = (0.,10.),
                 ωλprior  ::Tuple{Float64,Float64} = (0.,10.),
                 ωμprior  ::Tuple{Float64,Float64} = (0.,10.),
                 σ²prior  ::Float64                = 1e-1,
                 λprior   ::Float64                = 1e-1,
                 λi       ::Float64                = 1.,
                 ωxi      ::Float64                = 0.,
                 ωλi      ::Float64                = 0.,
                 ωμi      ::Float64                = 0.,
                 stbrl    ::Float64                = 1.,
                 fix_ωλ_ωμ::Bool                   = true,
                 delim    ::Char                   = '\t',
                 eol      ::Char                   = '\r')

  tip_values, tip_areas, tree, bts = 
    read_data(tree_file, data_file)

  X, Y, B, ncoup, δt, tree, si = 
    initialize_data(tip_values, tip_areas, min_dt, tree, bts)

  R = compete_mcmc(X, Y, ncoup, δt, tree.ed, tree.el, B,
                   niter     = niter,
                   nthin     = nthin,
                   nburn     = nburn,
                   ωxprior   = ωxprior,
                   ωλprior   = ωλprior,
                   ωμprior   = ωμprior,
                   σ²prior   = σ²prior,
                   λprior    = λprior,
                   out_file  = out_file,
                   λi        = λi,
                   ωxi       = ωxi,
                   ωλi       = ωλi,
                   ωμi       = ωμi,
                   σ²i       = si,
                   stbrl     = stbrl,
                   fix_ωλ_ωμ = fix_ωλ_ωμ)

  return R

end





