#=
ESSE Package module

Ignacio Quintero Mächler

t(-_-t)

Created 05 02 2019
=#

__precompile__()


module ESSE

  using ProgressMeter
  using DifferentialEquations
  using LinearAlgebra
  using RCall
  using Statistics
  using Random
  using DelimitedFiles

  export slice_sampler, 
         simesse

  include("approx_fun.jl")
  include("esse_eqs.jl")
  include("musse_eqs.jl")
  include("tree_utils.jl")
  include("loglik.jl")
  include("ode_solve.jl")
  include("slice_sampler_utils.jl")
  include("esse_slice_sampler.jl")
  include("musse_slice_sampler.jl")
  include("wrapper.jl")
  include("utils.jl")
  include("simulate_esse.jl")

end # module ESSE





