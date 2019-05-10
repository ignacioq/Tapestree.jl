#=
ESSE Package module

Ignacio Quintero Mächler

t(-_-t)

Created 05 02 2019
=#

module ESSE

  using ProgressMeter
  using DifferentialEquations
  using LinearAlgebra
  using RCall
  using Statistics
  using Random
  using DelimitedFiles

  export slice_sampler, 
         simesse,
         simulate_sse

  include("esse_eqs.jl")
  include("simulate_esse.jl")
  include("utils.jl")
  include("approx_fun.jl")
  include("esse_slice_sampler.jl")
  include("musse_eqs.jl")
  include("simulate_geoesse.jl")
  include("wrapper.jl")
  include("geohisse_eqs.jl")
  include("musse_slice_sampler.jl")
  include("slice_sampler_utils.jl")
  include("egeohisse_slice_sampler.jl")
  include("ode_solve.jl")
  include("tree_utils.jl")
  include("egeohisse_eqs.jl")
  include("loglik_geo.jl")

end # module ESSE
