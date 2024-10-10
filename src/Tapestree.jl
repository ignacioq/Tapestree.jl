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
using .ESSE: esse, simulate_sse, save_esse_sim
export esse, simulate_sse, save_esse_sim

using .TRIBE: tribe, simulate_tribe
export tribe, simulate_tribe

using .INSANE: read_newick, write_newick, iread, iwrite,
  sT_label, sTf_label, sTb, sTbd, sTfbd, iTb, iTce, iTct, iTbd, iTpbd, iTbdX, iTfbd, iTfbdX,
  iTgbmb, iTgbmce, iTgbmct, iTgbmbd,
  sim_cb, sim_cbd, sim_cfbd, sim_gbmb, sim_gbmce, sim_gbmct, sim_gbmbd, 
  sim_gbmfbd, sim_gbmpbd, sim_shift,
  iscrowntree, rm_stem!, fixtree!,
  insane_cb, insane_cbd, insane_cfbd, 
  insane_gbmb, insane_gbmce, insane_gbmct, insane_gbmbd, insane_gbmfbd,
  iquantile, imean, irange, sample,
  lb, lλ, lμ, remove_extinct, remove_unsampled, remove_fossils, fixedpos, fossilize!,
  e, b, d, ld, lb, nd, t, lt, dλ, dμ, dλc, dμc, epochs, l,
  ntips, ntipsalive, ntipsextinct, nfossils, ntipfossils, fixtree!, 
  trextract, reorder!, treeheight, treelength, _ctl, ltt, subclade, tiplabels, 
  time_rate, make_idf
export read_newick, write_newick, iread, iwrite,
  sT_label, sTf_label, sTb, sTbd, sTfbd, iTb, iTce, iTct, iTbd, iTpbd, iTbdX, iTfbd, iTfbdX,
  iTgbmb, iTgbmce, iTgbmct, iTgbmbd,
  sim_cb, sim_cbd, sim_cfbd, sim_gbmb, sim_gbmce, sim_gbmct, sim_gbmbd, 
  sim_gbmfbd, sim_gbmpbd, sim_shift,
  iscrowntree, rm_stem!, fixtree!,
  insane_cb, insane_cbd, insane_cfbd, 
  insane_gbmb, insane_gbmce, insane_gbmct, insane_gbmbd, insane_gbmfbd,
  iquantile, imean, irange, sample,
  lb, lλ, lμ, remove_extinct, remove_unsampled, remove_fossils, fixedpos, fossilize!,
  e, b, d, ld, lb, nd, t, lt, dλ, dμ, dλc, dμc, epochs, l,
  ntips, ntipsalive, ntipsextinct, nfossils, ntipfossils, fixtree!, 
  trextract, reorder!, treeheight, treelength, _ctl, ltt, subclade, tiplabels, 
  time_rate, make_idf

end # module Tapestree
