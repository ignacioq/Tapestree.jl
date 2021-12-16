#=

  "Tapestree" package

=#

module Tapestree

__precompile__(true)

#=
 Submodules
=#

# Utilities
include("Utils.jl")

# ESSE: Environmental and State Dependent Diversification
include("ESSE.jl")

# TRIBE: Trait and Range Interspecific Biogeographic Evolution
include("TRIBE.jl")

# INSANE: joINt Speciation And Niche Evolution
include("INSANE.jl")


#=
 Exported functions 
=#

using .ESSE: esse, simulate_sse
export esse, simulate_sse

using .TRIBE: tribe, simulate_tribe
export tribe, simulate_tribe

using .INSANE: read_newick, write_newick, 
  sTpb, sTbd, iTgbmpb, iTgbmce, iTgbmct, iTgbmbd,
  sim_cpb, sim_cbd, sim_gbmpb, sim_gbmce, sim_gbmct, sim_gbmbd,
  insane_cpb, insane_cbd, iscrowntree, rm_stem!,
  insane_gbmpb, insane_gbmce, insane_gbmct, insane_gbmbd, 
  iquantile, extract_vector!, mcmc_array, lλ, lμ, 
  remove_extinct!, remove_unsampled!, ntipsalive, ntips
export read_newick, write_newick, 
  sTpb, sTbd, iTgbmpb, iTgbmce, iTgbmct, iTgbmbd,
  sim_cpb, sim_cbd, sim_gbmpb, sim_gbmce, sim_gbmct, sim_gbmbd,
  insane_cpb, insane_cbd, iscrowntree, rm_stem!,
  insane_gbmpb, insane_gbmce, insane_gbmct, insane_gbmbd, 
  iquantile, extract_vector!, mcmc_array, lλ, lμ, 
  remove_extinct!, remove_unsampled!, ntipsalive, ntips

end # module Tapestree
