#=

Tapestree Package module

Ignacio Quintero Mächler

t(-_-t)

22 June 2017

=#


module Tapestree

# from standard library
using Random
using DelimitedFiles
using Statistics
using LinearAlgebra

# `true' dependencies
using RCall
using Optim
using ProgressMeter
using DifferentialEquations

export tribe, 
       simulate_tribe, 
       slice_sampler, 
       simesse,
       simulate_sse

# General utilities
include("utils.jl")
include("tree_utils.jl")
include("density_functions.jl")

#TRIBE
include("tribe/tribe_utils.jl")
include("tribe/cont_DA_prop.jl")
include("tribe/disc_DA_prop.jl")
include("tribe/tribe_data_initializer.jl")
include("tribe/area_lineage_averages.jl")
include("tribe/tribe_loglik_functions.jl")
include("tribe/mcmc_utils.jl")
include("tribe/mcmc.jl")
include("tribe/parameter_updates.jl")
include("tribe/mcmc_burn.jl")
include("tribe/under_prior.jl")
include("tribe/proposal_functions.jl")
include("tribe/tribe_wrapper.jl")
include("tribe/tribe_simulation_utils.jl")

#*SSE
include("sse/musse_eqs.jl")
include("sse/esse_eqs.jl")
include("sse/states_handling.jl")
include("sse/loglik_mu.jl")
include("sse/simulate_esse.jl")
include("sse/approx_fun.jl")
include("sse/model_definitions.jl")
include("sse/simulate_geoesse.jl")
include("sse/sse_wrapper.jl")
include("sse/slice_sampler_utils.jl")
include("sse/slice_sampler.jl")
include("sse/ode_solve.jl")
include("sse/egeohisse_eqs.jl")
include("sse/loglik_geo.jl")

end # module Tapestree
