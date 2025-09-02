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

using .INSANE: sT_label, sTf_label, sTb, sTbd, sTpe, sTfbd, sTfpe, 
  iTb, iTce, iTct, iTbd, iTpbd, iTfbd, sTxs, iTd, Ltt,
  read_newick, write_newick, write_nexus, iread, iwrite,
  sim_cb, sim_cbd, sim_cfbd, sim_cobd, sim_gbmb, sim_gbmce, sim_gbmct, sim_gbmbd, 
  sim_gbmfbd, sim_gbmobd, sim_gbmpbd, sim_shift, sim_occurrences, sim_dbm,
  iscrowntree, rm_stem!, fixtree!,
  insane_cb, insane_cbd, insane_cfbd, insane_cobd,
  insane_gbmb, insane_gbmce, insane_gbmct, insane_gbmbd, insane_gbmfbd,
  insane_gbmobd, insane_dbm,
  iquantile, imean, irange, sample,
  e, lλ, lμ, xv, lσ2,
  birth, death, logbirth, logdeath, diversification, turnover, epochs,
  remove_extinct, remove_unsampled, remove_fossils, fixedpos, fossilize!,
  ntips, ntipsalive, ntipsextinct, nfossils, ntipfossils, fixtree!, 
  trextract, reorder!, treeheight, treelength, _ctl, ltt, subclade, 
  time_rate, make_idf, plotω,
  tiplabels, labels, label
export sT_label, sTf_label, sTb, sTbd, sTpe, sTfbd, sTfpe, 
  iTb, iTce, iTct, iTbd, iTpbd, iTfbd, sTxs, iTd, Ltt,
  read_newick, write_newick, write_nexus, iread, iwrite,
  sim_cb, sim_cbd, sim_cfbd, sim_cobd, sim_gbmb, sim_gbmce, sim_gbmct, sim_gbmbd, 
  sim_gbmfbd, sim_gbmobd, sim_gbmpbd, sim_shift, sim_occurrences, sim_dbm,
  iscrowntree, rm_stem!, fixtree!,
  insane_cb, insane_cbd, insane_cfbd, insane_cobd,
  insane_gbmb, insane_gbmce, insane_gbmct, insane_gbmbd, insane_gbmfbd,
  insane_gbmobd, insane_dbm,
  iquantile, imean, irange, sample,
  e, lλ, lμ, xv, lσ2,
  birth, death, logbirth, logdeath, diversification, turnover, epochs,
  remove_extinct, remove_unsampled, remove_fossils, fixedpos, fossilize!,
  ntips, ntipsalive, ntipsextinct, nfossils, ntipfossils, fixtree!, 
  trextract, reorder!, treeheight, treelength, _ctl, ltt, subclade, 
  time_rate, make_idf, plotω,
  tiplabels, labels, label

end # module Tapestree
