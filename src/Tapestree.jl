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

# INSANE: joINt Speciation And Niche Evolution
include("INSANE.jl")

# TRIBE: Trait and Range Interspecific Biogeographic Evolution
include("TRIBE.jl")

# ESSE: Environmental and State Dependent Diversification
include("ESSE.jl")

#=
 Exported functions 
=#
using .ESSE: esse, simulate_sse, save_esse_sim
export esse, simulate_sse, save_esse_sim

using .TRIBE: tribe, simulate_tribe
export tribe, simulate_tribe

using .INSANE: read_newick, write_newick, write_nexus, iread, iwrite,
  sT_label, sTf_label, sTpb, sTbd, sTfbd, iTpb, iTce, iTct, iTbd, iTfbd, sTxs,
  Ltt,
  sim_cpb, sim_cbd, sim_cfbd, sim_gbmpb, sim_gbmce, sim_gbmct, sim_gbmbd, 
  sim_gbmfbd, sim_shift, sim_dbm,
  iscrowntree, rm_stem!, fixtree!,
  insane_cpb, insane_cbd, insane_cfbd, 
  insane_gbmpb, insane_gbmce, insane_gbmct, insane_gbmbd, insane_gbmfbd,
  insane_dbm,
  iquantile, imean, irange, sample,
  e, lλ, lμ, xv, lσ2,
  birth, death, logbirth, logdeath, diversification, turnover, epochs,
  remove_extinct, remove_unsampled, remove_fossils, fixedpos, fossilize!,
  ntips, ntipsalive, ntipsextinct, nfossils, ntipfossils, fixtree!, 
  trextract, reorder!, treeheight, treelength, _ctl, ltt, subclade, 
  time_rate, make_idf,
  tiplabels, labels, label
export read_newick, write_newick, write_nexus, iread, iwrite,
  sT_label, sTf_label, sTpb, sTbd, sTfbd, iTpb, iTce, iTct, iTbd, iTfbd, sTxs,
  Ltt,
  sim_cpb, sim_cbd, sim_cfbd, sim_gbmpb, sim_gbmce, sim_gbmct, sim_gbmbd, 
  sim_gbmfbd, sim_shift, sim_dbm,
  iscrowntree, rm_stem!, fixtree!,
  insane_cpb, insane_cbd, insane_cfbd, 
  insane_gbmpb, insane_gbmce, insane_gbmct, insane_gbmbd, insane_gbmfbd,
  insane_dbm,
  iquantile, imean, irange, sample,
  e, lλ, lμ, xv, lσ2,
  birth, death, logbirth, logdeath, diversification, turnover, epochs,
  remove_extinct, remove_unsampled, remove_fossils, fixedpos, fossilize!,
  ntips, ntipsalive, ntipsextinct, nfossils, ntipfossils, fixtree!, 
  trextract, reorder!, treeheight, treelength, _ctl, ltt, subclade, 
  time_rate, make_idf,
  tiplabels, labels, label

end # module Tapestree
