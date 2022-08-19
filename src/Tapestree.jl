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
  sTpb, sTbd, sTfbd, iTpb, iTce, iTct, iTbd, iTbdX, iTfbd, iTfbdX,
  iTgbmpb, iTgbmce, iTgbmct, iTgbmbd,
  sim_cpb, sim_cbd, sim_gbmpb, sim_gbmce, sim_gbmct, sim_gbmbd,
  iscrowntree, rm_stem!,
  insane_cpb, insane_cbd, insane_cfbd, 
  insane_gbmpb, insane_gbmce, insane_gbmct, insane_gbmbd, insane_gbmfbd,
  iquantile, extract_vector!, mcmc_array, lλ, lμ, 
  remove_extinct, remove_unsampled, remove_fossils, 
  b, d, ld, lb, nd, t,
  ntipsalive, ntips, treeheight, treelength, _ctl, ltt, subclade, tiplabels
export read_newick, write_newick, 
  sTpb, sTbd, sTfbd, iTpb, iTce, iTct, iTbd, iTbdX, iTfbd, iTfbdX,
  iTgbmpb, iTgbmce, iTgbmct, iTgbmbd,
  sim_cpb, sim_cbd, sim_gbmpb, sim_gbmce, sim_gbmct, sim_gbmbd,
  iscrowntree, rm_stem!,
  insane_cpb, insane_cbd, insane_cfbd, 
  insane_gbmpb, insane_gbmce, insane_gbmct, insane_gbmbd, insane_gbmfbd,
  iquantile, extract_vector!, mcmc_array, lλ, lμ, 
  remove_extinct, remove_unsampled, remove_fossils, 
  b, d, ld, lb, nd, t,
  ntipsalive, ntips, treeheight, treelength, _ctl, ltt, subclade, tiplabels

end # module Tapestree
