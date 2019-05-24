#=

Tapestree Package module

Ignacio Quintero Mächler

t(-_-t)

22 June 2017

=#


module Tapestree

using RCall
using Optim
using ProgressMeter
using Random
using DelimitedFiles

export tribe, simulate_tribe

include("utils.jl")
include("tribe_utils.jl")
include("cont_DA_prop.jl")
include("disc_DA_prop.jl")
include("data_initializer.jl")
include("area_lineage_averages.jl")
include("loglik_functions.jl")
include("mcmc.jl")
include("parameter_updates.jl")
include("mcmc_burn.jl")
include("under_prior.jl")
include("proposal_functions.jl")
include("tribe_wrapper.jl")
include("sim_utils.jl")

end # module Tapestree
