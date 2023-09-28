#=

  "INSANE: joINt Speciation And Niche Evolution"

=#

module INSANE

  using Random: randexp, randn!, shuffle!
  using SpecialFunctions: erf
  using DelimitedFiles: writedlm
  using ProgressMeter: Progress, next!
  using Statistics: quantile, mean
  using LoopVectorization: @avx
  using PlotUtils: cgrad
  using RecipesBase
  using Parsers: parse as Pparse
  using Distributions: Poisson, Uniform

  # other submodules dependencies
  using ..Utils

  const accerr = √eps()
  const epochs = Float64[720.0, 635.0, 538.8, 521.0, 509.0, 497.0, 485.4, 470.0, 458.4, 443.8, 443.8, 419.2, 393.3, 382.7, 358.9, 323.2, 298.9, 273.01, 259.51, 251.902, 247.2, 237, 201.3, 174.1, 163.5, 145.0, 100.5, 66.0, 56.0, 33.9, 23.03, 5.333, 2.58]

  # files
  include("insane/iTree.jl")
  include("insane/iTreeX.jl")
  include("insane/iTree_summary.jl")
  include("insane/iB.jl")
  include("insane/iB_manipulation.jl")
  include("insane/Ltt.jl")
  include("insane/iTree_data.jl")
  include("insane/iTree_manipulation.jl")
  include("insane/iTree_plot.jl")
  include("insane/iTree_IO.jl")
  include("insane/sim_cpb.jl")
  include("insane/sim_cbd.jl")
  include("insane/sim_cfbd.jl")
  include("insane/sim_gbmpb.jl")
  include("insane/sim_gbmbd.jl")
  include("insane/sim_gbmbd_efx.jl")
  include("insane/sim_gbmce.jl")
  include("insane/sim_gbmct.jl")
  include("insane/sim_gbmfbd.jl")
  include("insane/sim_shift.jl")
  include("insane/sim_cpbX.jl")
  include("insane/sim_cbdX.jl")
  include("insane/sim_cfbdX.jl")
  include("insane/ll_cpb.jl")
  include("insane/ll_cbd.jl")
  include("insane/ll_cfbd.jl")
  include("insane/ll_gbmpb.jl")
  include("insane/ll_gbmbd.jl")
  include("insane/ll_gbmce.jl")
  include("insane/ll_gbmct.jl")
  include("insane/ll_gbmfbd.jl")
  include("insane/ll_cpbX.jl")
  include("insane/ll_cbdX.jl")
  include("insane/ll_cfbdX.jl")
  include("insane/survival.jl")
  include("insane/mcmc_cpb.jl")
  include("insane/mcmc_cbd_gp.jl")
  include("insane/mcmc_cbd.jl")
  include("insane/mcmc_cfbd.jl")
  include("insane/mcmc_gbmpb.jl")
  include("insane/mcmc_gbmbd.jl")
  include("insane/mcmc_gbmbd_efx.jl")
  include("insane/mcmc_gbmce.jl")
  include("insane/mcmc_gbmct.jl")
  include("insane/mcmc_gbmfbd.jl")
  include("insane/mcmc_cpbX.jl")
  include("insane/mcmc_cbdX.jl")
  include("insane/mcmc_cfbdX.jl")
  include("insane/mhprop_gbmpb.jl")
  include("insane/mhprop_gbmbd.jl")
  include("insane/mhprop_gbmbd_efx.jl")
  include("insane/mhprop_gbmce.jl")
  include("insane/mhprop_gbmct.jl")
  include("insane/mhprop_gbmfbd.jl")
  include("insane/bm_utils.jl")
  include("insane/decoupled_trees.jl")
  include("insane/marginal_likelihood.jl")

  const iTd = Dict{String, DataType}("sTpb"  => sTpb,
                                     "sTbd"  => sTbd,
                                     "sTfbd" => sTfbd,
                                     "iTpb"  => iTpb,
                                     "iTce"  => iTce,
                                     "iTct"  => iTct,
                                     "iTbd"  => iTbd,
                                     "iTfbd" => iTfbd)
end # module INSANE