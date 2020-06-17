#=

  "TRIBE: Trait and Range Interspecific Biogeographic Evolution" submodule package

=#

module TRIBE

using Random: randexp
using DelimitedFiles: readdlm, writedlm
using ProgressMeter: Progress, next!
using Optim: minimizer, optimize, Options
using LinearAlgebra: BLAS.axpy!, BLAS.gemm!, eigvecs, eigvals, diagm
using RCall: reval, @rput

# other submodules dependencies
using ..Utils

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

end # module TRIBE
