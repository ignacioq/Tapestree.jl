#=

  "ESSE: Environmental and State Dependent Diversification" submodule package

=#

module ESSE

using Random: randexp, shuffle!
using DelimitedFiles: readdlm, writedlm
using ProgressMeter: Progress, next!
using OrdinaryDiffEq: ODEProblem, init, reinit!, solve!, Tsit5, DiffEqBase
using LinearAlgebra: BLAS.gemv!, rank, mul!, ldiv!, qr!
using Distributed: @distributed, @sync
using SharedArrays: SharedArray
using DistributedArrays: DArray, d_closeall, distribute
using RCall: @rput, reval

# other submodules dependencies
using ..Utils

# files
include("esse/musse_eqs.jl")
include("esse/esse_eqs.jl")
include("esse/states_handling.jl")
include("esse/loglik_mu.jl")
include("esse/simulate_esse.jl")
include("esse/approx_fun.jl")
include("esse/model_definitions.jl")
include("esse/simulate_sse_g.jl")
include("esse/simulate_esse_g.jl")
include("esse/sse_wrapper.jl")
include("esse/slice_sampler_utils.jl")
include("esse/slice_sampler.jl")
include("esse/mcmcmh.jl")
include("esse/mcmcmh_utils.jl")
include("esse/ode_solve.jl")
include("esse/egeohisse_eqs.jl")
include("esse/prepare_data.jl")
include("esse/prepare_ll.jl")
include("esse/loglik_geo.jl")
include("esse/marginal_ancestral_states.jl")
include("esse/loglik_flow.jl")
include("esse/mc3.jl")

end # module ESSE