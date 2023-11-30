#=

**esse** wrapper

Ignacio Quintero Mächler

t(-_-t)

September 26 2017

=#




"""
    esse(tree_file   ::String,
         out_file    ::String,
         h           ::Int64;
         states_file ::String            = "NaN",
         envdata_file::String            = "NaN",
         cov_mod     ::NTuple{M,String}  = ("",),
         node_ps     ::Tuple{Bool,Int64} = (true, 10),
         out_states  ::String            = "",
         constraints ::NTuple{N,String}  = (" ",),
         mvpars      ::NTuple{O,String}  = (" ",),
         niter       ::Int64             = 10_000,
         nthin       ::Int64             = 10,
         nburn       ::Int64             = 200,
         tune_int    ::Int64             = 100,
         nswap       ::Int64             = 10,
         ncch        ::Int64             = 1,
         parallel    ::Bool              = ncch > 1,
         dt           ::Float64          = 0.2,
         ntakew      ::Int64             = 100,
         winit       ::Float64           = 2.0,
         scale_y     ::NTuple{2,Bool}    = (true, false),
         algorithm   ::String            = "pruning",
         mc          ::String            = "slice",
         λpriors     ::Float64           = .1,
         μpriors     ::Float64           = .1,
         gpriors     ::Float64           = .1,
         lpriors     ::Float64           = .1,
         qpriors     ::Float64           = .1,
         βpriors     ::NTuple{2,Float64} = (0.0, 10.0),
         hpriors     ::Float64           = .1,
         optimal_w   ::Float64           = 0.8,
         tni         ::Float64           = 1.0,
         obj_ar      ::Float64           = 0.6,
         screen_print::Int64             = 5,
         Eδt         ::Float64           = 1e-3,
         ti          ::Float64           = 0.0,
         ρ           ::Array{Float64,1}  = [1.0]) where {M,N,O}

Run geographic **esse**. See tutorial for how these files should be specified.

# Arguments
- `tree_file::String`: full path to tree file.
- `out_file::String`: full path to write MCMC output.
- `h::Int64`: number of hidden states.
- `states_file ::String = "NaN"`: full path to states file. If `"NaN"`, no 
observed states are used (only hidden states).
- `envdata_file::String = "NaN"`: full path to covariates file. If `"NaN"`, no 
covariates are used (i.e., constant rates).
- `cov_mod::NTuple{M,String} = ("",)`: specifies which rates are affected by 
covariates: `s` for speciation, `e` for extinction, and `g` for colonization.
More than 1 is possible.
- `node_ps::Tuple{Bool,Int64} = (true, 10)`: first index specifies if posterior
marginal probabilities for nodes should be computed, second index the number of
iterations to be computed. 
- `out_states::String = ""`: full path to write node probabilities output.
- `constraints::NTuple{N,String} = (" ",)`: constraints for the model 
parameters.
- `mvpars::NTuple{O,String} = (" ",)`: which parameters should be multivariate 
when using slice sampling for better convergence.
- `niter::Int64 = 10_000`: number of iterations.
- `nthin::Int64 = 10`: frequency at which to record MCMC state.
- `nburn::Int64 = 200`: number of iterations to discard as burn-in.
- `tune_int::Int64 = 100`: number of iterations during `nburn` to tune proposal
window for MH.
- `nswap::Int64 = 10`: every iteration to try to swap chain likelihoods in MC3.
- `ncch::Int64 = 1`: number of chains.
- `parallel::Bool = false`: if parallel run.
- `dt::Float64 = 0.2`: temperature for MC3.
- `ntakew::Int64 = 100`: number of iterations from `nburn` to tune the window
for slice sampling.
- `winit::Float64 = 2.0`: initial window for slice sampling.
- `scale_y::NTuple{2,Bool} = (true, false)`: first index if scale covariates `y`
to [0,1], second, if scale covariates `y` all together between [0,1].
- `algorithm::String = "pruning"`: likelihood algorithm between `pruning` 
(recommended) or `flow`.
- `mc::String = "slice"`: which sampling `slice` (slice-sampling) or 
`mh` (metropolis-hasting).
- `λpriors::Float64 = 0.1`: rate of Exponential prior for speciation.
- `μpriors::Float64 = 0.1`: rate of Exponential prior for global extinction.
- `gpriors::Float64 = 0.1`: rate of Exponential prior for colonization.
- `lpriors::Float64 = 0.1`: rate of Exponential prior for local extinction.
- `qpriors::Float64 = 0.1`: rate of Exponential prior for hidden state 
transitions.
- `βpriors::NTuple{2,Float64} = (0.0, 10.0)`: mean and variance of Normal 
prior for effect of covariates.
- `hpriors::Float64 = 0.1`: rate of Exponential prior for differences between
hidden states.
- `optimal_w::Float64 = 0.8`: optimal window.
- `tni::Float64 = 1.0`: initial tuning for rates.
- `obj_ar::Float64 = 0.23`: objective acceptance rate.
- `screen_print::Int64 = 5`: seconds to wait to update screen log.
- `Eδt::Float64 = 1e-3`: for flow algorithm.
- `ti::Float64 = 0.0`: for flow algorithm.
- `ρ::Array{Float64,1} = [1.0]`: sampling fraction for each state (each area and
widespread).


# Returned values
  - Array of the mcmc parameters.
"""
function esse(tree_file   ::String,
              out_file    ::String,
              h           ::Int64;
              states_file ::String            = "NaN",
              envdata_file::String            = "NaN",
              cov_mod     ::NTuple{M,String}  = ("",),
              node_ps     ::Tuple{Bool,Int64} = (true, 10),
              out_states  ::String            = "",
              constraints ::NTuple{N,String}  = (" ",),
              mvpars      ::NTuple{O,String}  = (" ",),
              niter       ::Int64             = 10_000,
              nthin       ::Int64             = 10,
              nburn       ::Int64             = 200,
              tune_int    ::Int64             = 100,
              nswap       ::Int64             = 10,
              ncch        ::Int64             = 1,
              parallel    ::Bool              = ncch > 1,
              dt           ::Float64          = 0.2,
              ntakew      ::Int64             = 100,
              winit       ::Float64           = 2.0,
              scale_y     ::NTuple{2,Bool}    = (true, false),
              algorithm   ::String            = "pruning",
              mc          ::String            = "mh",
              λpriors     ::Float64           = 0.1,
              μpriors     ::Float64           = 0.1,
              gpriors     ::Float64           = 0.1,
              lpriors     ::Float64           = 0.1,
              qpriors     ::Float64           = 0.1,
              βpriors     ::NTuple{2,Float64} = (0.0, 10.0),
              hpriors     ::Float64           = 0.1,
              optimal_w   ::Float64           = 0.8,
              tni         ::Float64           = 1.0,
              obj_ar      ::Float64           = 0.23,
              screen_print::Int64             = 5,
              Eδt         ::Float64           = 1e-3,
              ti          ::Float64           = 0.0,
              ρ           ::Array{Float64,1}  = [1.0]) where {M,N,O}


  states = !occursin(r"^NaN$", states_file)
  enviro = !occursin(r"^NaN$", envdata_file)

  # read data 
  if states
    if enviro
      tv, ed, el, bts, x, y = 
        read_data_esse(states_file, tree_file, envdata_file)
    else
      tv, ed, el, bts = 
        read_data_sse(states_file, tree_file)
    end
  elseif enviro
    tv, ed, el, bts, x, y = 
      read_data_e(tree_file, envdata_file)
  else
    tv, ed, el, bts = 
      read_data_esse(tree_file)
  end

  if enviro
    # scale y
    if scale_y[1]
      # if scale each function separately or together
      if scale_y[2]
        ymin = minimum(y)
        ymax = maximum(y)
        for j in axes(y,2), i in axes(y,1)
          y[i,j] = (y[i,j] - ymin)/(ymax - ymin)
        end
      else
        ymin = minimum(y, dims = 1)
        ymax = maximum(y, dims = 1)
        for j in axes(y,2), i in axes(y,1)
          y[i,j] = (y[i,j] - ymin[j])/(ymax[j] - ymin[j])
        end
      end
    end
  end

  @info "Data for $(length(tv)) species successfully read"

  # prepare data
  if enviro
    X, p, fp, ed, trios, tdic, ns, ned, pupd, phid, nnps, nps, 
    mvps, nngps, mvhfs, hfgps, dcp, pardic, k, h, ny, model, 
    af!, assign_hidfacs!, abts, bts, E0 = 
      prepare_data(cov_mod, tv, x, y, ed, el, ρ, h, ncch, constraints, mvpars,
        parallel)
  else
    X, p, fp, ed, trios, tdic, ns, ned, pupd, phid, nnps, nps, 
    mvps, nngps, mvhfs, hfgps, dcp, pardic, k, h, ny, model, 
    af!, assign_hidfacs!, abts, bts, E0 = 
      prepare_data(tv, ed, el, ρ, h, ncch, constraints, mvpars, parallel)
  end

  @info "Data successfully prepared"

  @debug sort!(collect(pardic), by = x -> x[2])

  ## make likelihood function
  # flow algorithm
  if occursin(r"^[f|F][A-za-z]*", algorithm) 

    # prepare likelihood
    Gt, Et, lbts, nets, λevent!, rootll = 
      prepare_ll(p[1], bts, E0, k, h, ny, ns, ned, model, Eδt, ti, abts, af!)

    # make likelihood function
    llf = make_loglik(Gt, Et, X, trios, lbts, bts, ns, ned, nets, 
                      λevent!, rootll)

  # pruning algorithm
  elseif occursin(r"^[p|P][A-za-z]*", algorithm)

    # prepare likelihood
    X, U, A, int, λevent!, rootll, rootll_nj!, abts1, abts2 = 
      prepare_ll(X, p[1], E0, k, h, ny, ns, model, abts, af!)

    # make likelihood function
    llf = make_loglik(X, U, abts1, abts2, trios, int, 
      λevent!, rootll, ns, ned)

    llfnj = make_loglik_nj(X, tdic, abts1, abts2, trios, int, 
      λevent!, rootll_nj!, ns, ned)

    @info "Likelihood based on pruning algorithm prepared"

  else
    @error "no matching likelihood for algorithm: $algorithm"
  end

  λupds, μupds, lupds, gupds, qupds, βupds, hfps, βp_m, βp_v = 
    make_prior_updates(pupd, phid, mvhfs, hfgps, βpriors, k, h, ny, model)

  # create prior function
  lpf = make_lpf(λupds, μupds, lupds, gupds, qupds, βupds, hfps, 
          λpriors, μpriors, gpriors, lpriors, qpriors, βp_m, βp_v, hpriors)

  # create posterior functions
  lhf = make_lhf(llf, lpf, assign_hidfacs!, dcp)

  spf = make_state_posteriors(llf, lpf, llfnj, X, U, A, ns, ned, k, h)

  # number of parameters
  npars = length(pardic)

  if occursin(r"^[m|M][A-za-z]*[h|H][A-za-z]*", mc)

    @info "Running Metropolis-Hastings Markov chain"

    R = mcmcmh(lhf, p, fp, nnps, nps, phid, npars, niter, nthin, nburn, nswap, 
      ncch, tni, tune_int, dt, out_file, pardic, screen_print, obj_ar)

  elseif occursin(r"^[s|S][A-za-z]*", mc)

    @info "Running Slice-Sampler Markov chain"

    # run slice-sampler
    R = slice_sampler(lhf, p, fp, nnps, nps, phid, mvps, 
        nngps, mvhfs, hfgps, npars, niter, nthin, nburn, ntakew, nswap, 
        ncch, winit, optimal_w, dt, out_file, pardic, screen_print)
  else

    @error "No matching Markov chain: $mc"
  end

  # write chain output
  write_ssr(R, pardic, out_file)

  # run ancestral state marginal probabilities
  if node_ps[1]
    @info "Estimating node marginal states probabilities..."
    S = sample_node_ps(R, A, spf, node_ps[2], ns, ned)

    # make dictionary for names
    nodic = Dict{String, Int64}()

    ed2 = ed[:,2]
    for j in Base.OneTo(ned), i in Base.OneTo(ns)
      push!(nodic, string("node_", ed2[j],"_state_",i) => i + (j-1)*(ns))
    end

    # write ancestral node reconstruction
    write_ssr(S, nodic, out_states)
  end

  @info "Finished"

  return R
end




"""
    esse(tv          ::Dict{Int64,Array{Float64,1}},
         ed          ::Array{Int64,2}, 
         el          ::Array{Float64,1}, 
         x           ::Array{Float64,1},
         y           ::Array{Float64,L}, 
         cov_mod     ::NTuple{M,String},
         out_file    ::String,
         h           ::Int64;
         constraints ::NTuple{N,String}  = (" ",),
         mvpars      ::NTuple{O,String}  = ("lambda = beta",),
         niter       ::Int64             = 10_000,
         nthin       ::Int64             = 10,
         nburn       ::Int64             = 200,
         ncch        ::Int64             = 1,
         ntakew      ::Int64             = 100,
         winit       ::Float64             = 2.0,
         scale_y     ::NTuple{2,Bool}    = (true, false),
         algorithm   ::String            = "pruning",
         λpriors     ::Float64           = .1,
         μpriors     ::Float64           = .1,
         gpriors     ::Float64           = .1,
         lpriors     ::Float64           = .1,
         qpriors     ::Float64           = .1,
         βpriors     ::NTuple{2,Float64} = (0.0, 10.0),
         hpriors     ::Float64           = .1,
         optimal_w   ::Float64           = 0.8,
         screen_print::Int64             = 5,
         Eδt         ::Float64           = 1e-3,
         ti          ::Float64           = 0.0,
         ρ           ::Array{Float64,1}  = [1.0]) where {L,M,N,O}

Wrapper for running a **esse** model from simulations.
"""
function esse(tv          ::Dict{Int64,Array{Float64,1}},
              ed          ::Array{Int64,2}, 
              el          ::Array{Float64,1}, 
              x           ::Array{Float64,1},
              y           ::Array{Float64,L}, 
              cov_mod     ::NTuple{M,String},
              out_file    ::String,
              h           ::Int64;
              node_ps     ::Tuple{Bool,Int64} = (true, 10),
              out_states  ::String            = "",
              constraints ::NTuple{N,String}  = (" ",),
              mvpars      ::NTuple{O,String}  = (" ",),
              niter       ::Int64             = 10_000,
              nthin       ::Int64             = 10,
              nburn       ::Int64             = 200,
              tune_int    ::Int64             = 100,
              nswap       ::Int64             = 10,
              ncch        ::Int64             = 1,
              parallel    ::Bool              = false,
              dt           ::Float64          = 0.2,
              ntakew      ::Int64             = 100,
              winit       ::Float64           = 2.0,
              scale_y     ::NTuple{2,Bool}    = (true, false),
              algorithm   ::String            = "pruning",
              mc          ::String            = "slice",
              λpriors     ::Float64           = .1,
              μpriors     ::Float64           = .1,
              gpriors     ::Float64           = .1,
              lpriors     ::Float64           = .1,
              qpriors     ::Float64           = .1,
              βpriors     ::NTuple{2,Float64} = (0.0, 10.0),
              hpriors     ::Float64           = .1,
              optimal_w   ::Float64           = 0.8,
              tni         ::Float64           = 1.0,
              obj_ar      ::Float64           = 0.6,
              screen_print::Int64             = 5,
              Eδt         ::Float64           = 1e-3,
              ti          ::Float64           = 0.0,
              ρ           ::Array{Float64,1}  = [1.0]) where {L,M,N,O}


  # scale y
  if scale_y[1]
    # if scale each function separately or together
    if scale_y[2]
      ymin = minimum(y)
      ymax = maximum(y)
      for j in axes(y,2), i in axes(y,1)
        y[i,j] = (y[i,j] - ymin)/(ymax - ymin)
      end
    else
      ymin = minimum(y, dims = 1)
      ymax = maximum(y, dims = 1)
      for j in axes(y,2), i in axes(y,1)
        y[i,j] = (y[i,j] - ymin[j])/(ymax[j] - ymin[j])
      end
    end
  end

  # prepare data
  X, p, fp, ed, trios, tdic, ns, ned, pupd, phid, nnps, nps, 
  mvps, nngps, mvhfs, hfgps, dcp, pardic, k, h, ny, model, 
  af!, assign_hidfacs!, abts, bts, E0 = 
    prepare_data(cov_mod, tv, x, y, ed, el, ρ, h, ncch, constraints, mvpars,
      parallel)

  @info "Data successfully prepared"

  @debug sort!(collect(pardic), by = x -> x[2])

  ## make likelihood function
  # flow algorithm
  if occursin(r"^[f|F][A-za-z]*", algorithm) 

    # prepare likelihood
    Gt, Et, lbts, nets, λevent!, rootll = 
      prepare_ll(p[1], bts, E0, k, h, ny, ns, ned, model, Eδt, ti, abts, af!)

    # make likelihood function
    llf = make_loglik(Gt, Et, X, trios, lbts, bts, ns, ned, nets, 
                      λevent!, rootll)

  # pruning algorithm
  elseif occursin(r"^[p|P][A-za-z]*", algorithm)

    # prepare likelihood
    X, U, A, int, λevent!, rootll, rootll_nj!, abts1, abts2 = 
      prepare_ll(X, p[1], E0, k, h, ny, ns, model, abts, af!)

    # make likelihood function
    llf = make_loglik(X, U, abts1, abts2, trios, int, 
      λevent!, rootll, ns, ned)

    llfnj = make_loglik_nj(X, tdic, abts1, abts2, trios, int, 
      λevent!, rootll_nj!, ns, ned)

    @info "Likelihood based on pruning algorithm prepared"

  else
    @error "no matching likelihood for algorithm: $algorithm"
  end

  λupds, μupds, lupds, gupds, qupds, βupds, hfps, βp_m, βp_v = 
    make_prior_updates(pupd, phid, mvhfs, hfgps, βpriors, k, h, ny, model)

  # create prior function
  lpf = make_lpf(λupds, μupds, lupds, gupds, qupds, βupds, hfps, 
          λpriors, μpriors, gpriors, lpriors, qpriors, βp_m, βp_v, hpriors)

  # create posterior functions
  lhf = make_lhf(llf, lpf, assign_hidfacs!, dcp)

  spf = make_state_posteriors(llf, lpf, llfnj, X, U, A, ns, ned, k, h)

  # number of parameters
  npars = length(pardic)

  if occursin(r"^[m|M][A-za-z]*[h|H][A-za-z]*", mc)

    @info "Running Metropolis-Hastings Markov chain"

    R = mcmcmh(lhf, p, fp, nnps, nps, phid, npars, niter, nthin, nburn, nswap, 
      ncch, tni, tune_int, dt, out_file, pardic, screen_print, obj_ar)

  elseif occursin(r"^[s|S][A-za-z]*", mc)

    @info "Running Slice-Sampler Markov chain"

    # run slice-sampler
    R = slice_sampler(lhf, p, fp, nnps, nps, phid, mvps, 
        nngps, mvhfs, hfgps, npars, niter, nthin, nburn, ntakew, nswap, 
        ncch, winit, optimal_w, dt, out_file, pardic, screen_print)
  else

    @error "No matching Markov chain: $mc"
  end

  # write chain output
  write_ssr(R, pardic, out_file)

  # run ancestral state marginal probabilities
  if node_ps[1]
    @info "Estimating node marginal states probabilities..."
    S = sample_node_ps(R, A, spf, node_ps[2], ns, ned)

    # make dictionary for names
    nodic = Dict{String, Int64}()

    ed2 = ed[:,2]
    for j in Base.OneTo(ned), i in Base.OneTo(ns)
      push!(nodic, string("node_", ed2[j],"_state_",i) => i + (j-1)*(ns))
    end

    # write ancestral node reconstruction
    write_ssr(S, nodic, out_states)
  end

  @info "Finished"

  return R
end




"""
    read_data_esse(states_file ::String, 
                   tree_file   ::String, 
                   envdata_file::String)

Process tree and state and environmental data file to run **esse**.
"""
function read_data_esse(states_file ::String, 
                        tree_file   ::String, 
                        envdata_file::String)

  # read tree in postorder and assign to objects
  tree, bts = read_tree(tree_file, order = "postorder", branching_times = true)
  ntip = tree.nnod + 1
  ed   = tree.ed
  el   = tree.el
  tlab = tree.tlab

  # assign tip labels to edge numbers
  tip_labels = Dict{String,Integer}()
  for i in Base.OneTo(lastindex(tlab))
    tip_labels[tlab[i]] = i
  end

  # read states text file
  data = readdlm(states_file)

  if size(data,1) != ntip
    data = readdlm(states_file, '\t', '\r')
  end

  if size(data,1) != ntip
    data = readdlm(states_file, '\t', '\n')
  end

  if size(data,1) != ntip 
    error("Data file cannot be made of the right dimensions.\n Make sure the data file has the same number of rows as tips in the tree")
  end

  data_tlab    = convert(Array{String,1}, data[:,1])
  data_states  = convert(Array{Float64,2},  data[:,2:end])

  # create dictionary
  tip_states = Dict(tip_labels[val] => data_states[i,:] 
                   for (i,val) = enumerate(data_tlab))

  # process environmental data file
  envdata = readdlm(envdata_file)

  x = envdata[:,1]
  y = envdata[:,2:end]

  return tip_states, ed, el, bts, x, y
end




"""
    read_data_sse(states_file ::String, 
                  tree_file   ::String)

Process tree and state file to run **esse**.
"""
function read_data_sse(states_file ::String, 
                       tree_file   ::String)

  # read tree in postorder and assign to objects
  tree, bts = read_tree(tree_file, order = "postorder", branching_times = true)
  ntip = tree.nnod + 1
  ed   = tree.ed
  el   = tree.el
  tlab = tree.tlab

  # assign tip labels to edge numbers
  tip_labels = Dict{String,Integer}()
  for i in Base.OneTo(lastindex(tlab))
    tip_labels[tlab[i]] = i
  end

  # read states text file
  data = readdlm(states_file)

  if size(data,1) != ntip
    data = readdlm(states_file, '\t', '\r')
  end

  if size(data,1) != ntip
    data = readdlm(states_file, '\t', '\n')
  end

  if size(data,1) != ntip 
    error("Data file cannot be made of the right dimensions.\n Make sure the data file has the same number of rows as tips in the tree")
  end

  data_tlab    = convert(Array{String,1}, data[:,1])
  data_states  = convert(Array{Float64,2},  data[:,2:end])

  # create dictionary
  tip_states = Dict(tip_labels[val] => data_states[i,:] 
                   for (i,val) = enumerate(data_tlab))

  return tip_states, ed, el, bts
end



"""
    read_data_esse(tree_file   ::String, 
                   envdata_file::String)

Process tree and state and environmental data file to run **esse**.
"""
function read_data_e(tree_file   ::String, 
                     envdata_file::String)

  # read tree in postorder and assign to objects
  tree, bts = read_tree(tree_file, order = "postorder", branching_times = true)
  ntip = tree.nnod + 1
  ed   = tree.ed
  el   = tree.el
  tlab = tree.tlab

  # create dictionary
  tip_states = Dict(i => [1.0] for i = Base.OneTo(ntip))

  # process environmental data file
  envdata = readdlm(envdata_file)

  x = envdata[:,1]
  y = envdata[:,2:end]

  return tip_states, ed, el, bts, x, y
end




"""
    read_data_esse(states_file ::String, 
                   tree_file   ::String, 
                   envdata_file::String)

Process tree and state and environmental data file to run **esse**.
"""
function read_data_esse(tree_file   ::String)

  # read tree in postorder and assign to objects
  tree, bts = read_tree(tree_file, order = "postorder", branching_times = true)
  ntip = tree.nnod + 1
  ed   = tree.ed
  el   = tree.el
  tlab = tree.tlab

  # assign tip labels to edge numbers
  tip_labels = Dict{String,Integer}()
  for i in Base.OneTo(lastindex(tlab))
    tip_labels[tlab[i]] = i
  end

  # create dictionary
  tip_states = Dict(i => [1.0] for i = Base.OneTo(ntip))

  return tip_states, ed, el, bts
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
  col_nam = ["Iteration", "Posterior"]

  for (k,v) in sort!(collect(pardic), by = x -> x[2])
    push!(col_nam, k)
  end

  R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

  writedlm(out_file*".log", R)
end





"""
  write_ssr(R       ::SharedArray{Float64,2}, 
            pardic  ::Dict{String,Int64},
            out_file::String,
            cits    ::Array{UnitRange{Int64},1})

Write the samples from multiple chains of MC sampler data frame 
given a Dictionary of parameters.
"""
function write_ssr(R       ::SharedArray{Float64,2}, 
                   pardic  ::Dict{String,Int64},
                   out_file::String,
                   cits    ::Array{UnitRange{Int64},1},
                   ci      ::Int64)

  # column names
  col_nam = ["Iteration", "Posterior"]

  for (k,v) in sort!(collect(pardic), by = x -> x[2])
    push!(col_nam, k)
  end

  ri = vcat(reshape(col_nam, 1, lastindex(col_nam)), R[cits[ci],:])

  writedlm(out_file*"_chain_$ci.log", ri)
end







