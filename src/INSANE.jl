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

  # other submodules dependencies
  using ..Utils

  const err = √eps()

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

end # module INSANE