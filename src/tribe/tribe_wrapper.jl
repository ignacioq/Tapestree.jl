#=

tribe wrapper.

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#



"""
    tribe(tree_file   ::String,
          data_file   ::String,
          out_file    ::String;
          min_dt      ::Float64           = 0.01,
          niter       ::Int64             = 50_000,
          nburn       ::Int64             = 50_000,
          nthin       ::Int64             = 500,
          saveXY      ::Tuple{Bool,Int64} = (false, 1_000),
          saveDM      ::Tuple{Bool,Int64} = (false, 1_000),
          ωxprior     ::NTuple{2,Float64} = (0.,10.),
          ω1prior     ::NTuple{2,Float64} = (0.,10.),
          ω0prior     ::NTuple{2,Float64} = (0.,10.),
          σ²prior     ::Float64           = 1e-1,
          λprior      ::Float64           = 1e-1,
          weight      ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),
          λ1i         ::Float64           = 1.0,
          λ0i         ::Float64           = 0.5,
          ωxi         ::Float64           = 0.0,
          ω1i         ::Float64           = 0.0,
          ω0i         ::Float64           = 0.0,
          fix_ωx      ::Bool              = false,
          fix_ω1      ::Bool              = false,
          fix_ω0      ::Bool              = false,
          delim       ::Char              = '\t',
          eol         ::Char              = '\r',
          screen_print::Int64             = 5)

Run **tribe** model. 

...
# Arguments
- `tree_file::String`: full path to tree file.
- `data_file::String`: full path to data file.
- `out_file::String`: full path to write MCMC output.
- `min_dt::Float64 = 0.01`: the percentage of tree height allowed for 
discretization (lower values are more precise but take longer).
- `niter::Int64 = 50_000`: the number of iterations.
- `nburn::Int64 = 50_000`: the number of iterations in the adaptive burn-in phase.
- `nthin::Int64 = 500`: the iteration sampling frequency.
- `saveXY::Tuple{Bool,Int64} = (false, 1_000)`: first index to 
save (or not) data augmented histories, second index for sampling frequency.
- `saveDM::Tuple{Bool,Int64} = (false, 1_000)`: a tuple of length 2: first is a boolean to save (or not) data augmented deterministic effects, second an integer for sampling frequency.
- `ωxprior::NTuple{2,Float64} = (0.,10.)`: a tuple of length 2 for the normal prior of ωx, first the mean, second the variance.
- `ω1prior::NTuple{2,Float64} = (0.,10.)`: a tuple of length 2 for the normal prior of ω1, first the mean, second the variance.
- `ω0prior::NTuple{2,Float64} = (0.,10.)`: a tuple of length 2 for the normal prior of ω0, first the mean, second the variance.
- `σ²prior::Float64 = 1e-1`: a float for the mean of the exponential prior for σ².
- `λprior::Float64 = 1e-1`: a float for the mean of the exponential prior for both λs.
- `weight::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3)`: a tuple of length 5 specifying the probabilities to update σ², ωx, ω1 & ω0, and λ1 & λ0 respectively.
- `λ1i::Float64 = 1.0`: a float for the starting value for λ1.
- `λ0i::Float64 = 0.5`: a float for the starting value for λ0.
- `ωxi::Float64 = 0.0`: a float for the starting value for ωx.
- `ω1i::Float64 = 0.0`: a float for the starting value for ω1.
- `ω0i::Float64 = 0.0`: a float for the starting value for ω0.
- `fix_ωx::Bool = false`: a boolean to make inference without ωx.
- `fix_ω1::Bool = false`: a boolean to make inference without ω1.
- `fix_ω0::Bool = false`: a boolean to make inference without ω0.
- `delim::Char= '\t'`: for ddlm.
- `eol::Char= '\r'`: for ddlm.
- `screen_print::Int64 = 5`: seconds to wait to update screen log.
...

...
# Returned values
  - Array of the mcmc parameters.
...
"""
function tribe(tree_file   ::String,
               data_file   ::String,
               out_file    ::String;
               min_dt      ::Float64           = 0.01,
               niter       ::Int64             = 50_000,
               nburn       ::Int64             = 50_000,
               nthin       ::Int64             = 500,
               saveXY      ::Tuple{Bool,Int64} = (false, 1_000),
               saveDM      ::Tuple{Bool,Int64} = (false, 1_000),
               ωxprior     ::NTuple{2,Float64} = (0.,10.),
               ω1prior     ::NTuple{2,Float64} = (0.,10.),
               ω0prior     ::NTuple{2,Float64} = (0.,10.),
               σ²prior     ::Float64           = 1e-1,
               λprior      ::Float64           = 1e-1,
               weight      ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),
               λ1i         ::Float64           = 1.0,
               λ0i         ::Float64           = 0.5,
               ωxi         ::Float64           = 0.0,
               ω1i         ::Float64           = 0.0,
               ω0i         ::Float64           = 0.0,
               fix_ωx      ::Bool              = false,
               fix_ω1      ::Bool              = false,
               fix_ω0      ::Bool              = false,
               delim       ::Char              = '\t',
               eol         ::Char              = '\r',
               screen_print::Int64             = 5)

  tip_values, tip_areas, tree, bts = 
    read_data_tribe(tree_file, data_file)

  X, Y, B, ncoup, δt, tree, si = 
    initialize_data(tip_values, tip_areas, min_dt, tree, bts)

  R = tribe_mcmc(X, Y, ncoup, δt,
                 deepcopy(tree.ed), 
                 deepcopy(tree.el), 
                 B,
                 niter    = niter,
                 nthin    = nthin,
                 nburn    = nburn,
                 saveXY   = saveXY,
                 saveDM   = saveDM,
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
                 fix_ωx   = fix_ωx,
                 fix_ω1   = fix_ω1,
                 fix_ω0   = fix_ω0,
                 screen_print = screen_print)

  return R

end






"""
    tribe(tip_values::Dict{Int64,Float64}, 
          tip_areas ::Dict{Int64,Array{Int64,1}},
          tree      ::rtree, 
          bts       ::Array{Float64,1},
          out_file  ::String;
          min_dt    ::Float64           = 0.01,
          niter     ::Int64             = 500_000,
          nburn     ::Int64             = 500_000,
          nthin     ::Int64             = 1_000,
          saveXY    ::Tuple{Bool,Int64} = (false, 1_000),
          saveDM    ::Tuple{Bool,Int64} = (false, 1_000),
          ωxprior   ::NTuple{2,Float64} = (0.,10.),
          ω1prior   ::NTuple{2,Float64} = (0.,10.),
          ω0prior   ::NTuple{2,Float64} = (0.,10.),
          σ²prior   ::Float64           = 1e-1,
          λprior    ::Float64           = 1e-1,
          weight    ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),
          λ1i       ::Float64           = 1.0,
          λ0i       ::Float64           = 0.4,
          ωxi       ::Float64           = 0.,
          ω1i       ::Float64           = 0.,
          ω0i       ::Float64           = 0.,
          fix_ωx    ::Bool              = false,
          fix_ω1    ::Bool              = false,
          fix_ω0    ::Bool              = false,
          delim     ::Char              = '\t',
          eol       ::Char              = '\r')

Run tribe for simulations. Wrapper for all functions.
"""
function tribe(tip_values::Dict{Int64,Float64}, 
               tip_areas ::Dict{Int64,Array{Int64,1}},
               tree      ::rtree, 
               bts       ::Array{Float64,1},
               out_file  ::String;
               min_dt    ::Float64           = 0.01,
               niter     ::Int64             = 500_000,
               nburn     ::Int64             = 500_000,
               nthin     ::Int64             = 1_000,
               saveXY    ::Tuple{Bool,Int64} = (false, 1_000),
               saveDM    ::Tuple{Bool,Int64} = (false, 1_000),
               ωxprior   ::NTuple{2,Float64} = (0.,10.),
               ω1prior   ::NTuple{2,Float64} = (0.,10.),
               ω0prior   ::NTuple{2,Float64} = (0.,10.),
               σ²prior   ::Float64           = 1e-1,
               λprior    ::Float64           = 1e-1,
               weight    ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),
               λ1i       ::Float64           = 1.0,
               λ0i       ::Float64           = 0.4,
               ωxi       ::Float64           = 0.,
               ω1i       ::Float64           = 0.,
               ω0i       ::Float64           = 0.,
               fix_ωx    ::Bool              = false,
               fix_ω1    ::Bool              = false,
               fix_ω0    ::Bool              = false,
               delim     ::Char              = '\t',
               eol       ::Char              = '\r',
               screen_print::Int64           = 5)

  X, Y, B, ncoup, δt, tree, si = 
    initialize_data(tip_values, tip_areas, min_dt, tree, bts)

  R = tribe_mcmc(X, Y, ncoup, δt,
                 deepcopy(tree.ed), 
                 deepcopy(tree.el), 
                 B,
                 niter    = niter,
                 nthin    = nthin,
                 nburn    = nburn,
                 saveXY   = saveXY,
                 saveDM   = saveDM,
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
                 fix_ωx   = fix_ωx,
                 fix_ω1   = fix_ω1,
                 fix_ω0   = fix_ω0,
                 screen_print = screen_print)

  return R

end






"""
    tribe(out_file::String;
          niter   ::Int64             = 500_000,
          nburn   ::Int64             = 500_000,
          nthin   ::Int64             = 1_000,
          ωxprior ::NTuple{2,Float64} = (0.,10.),
          ω1prior ::NTuple{2,Float64} = (0.,10.),
          ω0prior ::NTuple{2,Float64} = (0.,10.),
          σ²prior ::Float64           = 1e-1,
          λprior  ::Float64           = 1e-1,
          weight  ::NTuple{4,Float64} = (0.15,0.05,0.02,0.02),
          σ²i     ::Float64           = 1.,
          ωxi     ::Float64           = 0.,
          ω1i     ::Float64           = 0.01,
          ω0i     ::Float64           = 0.01,
          λ1i     ::Float64           = 1.0,
          λ0i     ::Float64           = 0.2,
          fix_ωx  ::Bool              = false,
          fix_ω1  ::Bool              = false,
          fix_ω0  ::Bool              = false)

Run tribe **under the prior**. Wrapper for all functions.
"""
function tribe(out_file::String;
               niter   ::Int64             = 500_000,
               nburn   ::Int64             = 500_000,
               nthin   ::Int64             = 1_000,
               ωxprior ::NTuple{2,Float64} = (0.,10.),
               ω1prior ::NTuple{2,Float64} = (0.,10.),
               ω0prior ::NTuple{2,Float64} = (0.,10.),
               σ²prior ::Float64           = 1e-1,
               λprior  ::Float64           = 1e-1,
               weight  ::NTuple{4,Float64} = (0.15,0.05,0.02,0.02),
               σ²i     ::Float64           = 1.,
               ωxi     ::Float64           = 0.,
               ω1i     ::Float64           = 0.01,
               ω0i     ::Float64           = 0.01,
               λ1i     ::Float64           = 1.0,
               λ0i     ::Float64           = 0.2,
               fix_ωx  ::Bool              = false,
               fix_ω1  ::Bool              = false,
               fix_ω0  ::Bool              = false,
               screen_print::Int64         = 5)

  R = tribe_mcmc(out_file,
                 niter    = niter,
                 nthin    = nthin,
                 nburn    = nburn,
                 ωxprior  = ωxprior,
                 ω1prior  = ω1prior,
                 ω0prior  = ω0prior,
                 σ²prior  = σ²prior,
                 λprior   = λprior,
                 weight   = weight,
                 λ1i      = λ1i,
                 λ0i      = λ0i,
                 ωxi      = ωxi,
                 ω1i      = ω1i,
                 ω0i      = ω0i,
                 σ²i      = σ²i,
                 fix_ωx   = fix_ωx,
                 fix_ω1   = fix_ω1,
                 fix_ω0   = fix_ω0,
                 screen_print = screen_print)

  return R

end



