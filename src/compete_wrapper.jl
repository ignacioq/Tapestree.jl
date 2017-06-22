"""
Wrapper for running competition model

Ignacio Quintero

t(-_-t)

June 20 2017
"""

function compete(tree_file::String,
                 data_file::String;
                 m        ::Int64                  = 100,
                 niter    ::Int64                  = 500_000,
                 nthin    ::Int64                  = 1_000,
                 nburn    ::Int64                  = 500_000,
                 ωxprior  ::Tuple{Float64,Float64} = (0.,10.),
                 ωλprior  ::Tuple{Float64,Float64} = (0.,10.),
                 ωμprior  ::Tuple{Float64,Float64} = (0.,10.),
                 σ²prior  ::Float64                = 1e-1,
                 λprior   ::Float64                = 1e-1,
                 dir_out   ::String                 = "/data/turnover/model/",
                 out_file ::String                 = "compete",
                 λi       ::Float64                = 1.,
                 ωxi      ::Float64                = 0.,
                 ωλi      ::Float64                = 0.,
                 ωμi      ::Float64                = 0.,
                 stbrl    ::Float64                = 1.,
                 fix_ωλ_ωμ::Bool                   = true,
                 delim    ::Char                   = '\t',
                 eol      ::Char                   = '\r')



  tip_values, tip_areas, tree, bts = 
    read_data(tree_file, data_file, delim = delim, eol = eol)


  X, Y, B, ncoup, δt, tree, si = 
    initialize_data(tip_values, tip_areas, m, tree, bts)


  R = compete_mcmc(X, Y, ncoup, δt, tree.ed, tree.el, B,
                   niter     = niter,
                   nthin     = nthin,
                   nburn     = nburn,
                   ωxprior   = ωxprior,
                   ωλprior   = ωλprior,
                   ωμprior   = ωμprior,
                   σ²prior   = σ²prior,
                   λprior    = λprior,
                   dir_out   = dir_out,
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





