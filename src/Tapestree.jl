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
include("denisty_functions.jl")

#TRIBE
include("tribe_utils.jl")
include("cont_DA_prop.jl")
include("disc_DA_prop.jl")
include("data_initializer.jl")
include("area_lineage_averages.jl")
include("tribe_loglik_functions.jl")
include("mcmc.jl")
include("parameter_updates.jl")
include("mcmc_burn.jl")
include("under_prior.jl")
include("proposal_functions.jl")
include("tribe_wrapper.jl")
include("tribe_simulation_utils.jl")

#ESSE
include("esse_eqs.jl")
include("simulate_esse.jl")
include("approx_fun.jl")
include("model definitions.jl")
include("musse_eqs.jl")
include("simulate_geoesse.jl")
include("sse_wrapper.jl")
include("slice_sampler_utils.jl")
include("slice_sampler.jl")
include("ode_solve.jl")
include("egeohisse_eqs.jl")
include("loglik_geo.jl")

end # module Tapestree
