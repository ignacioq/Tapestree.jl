#=

  Utilities submodule package

=#

module Utils

# package dependencies
using RCall: reval, rcopy

include("utils/utils.jl")
include("utils/tree_utils.jl")
include("utils/density_functions.jl")

export rtree, read_tree, make_ape_tree, maketriads, abs_time_branches,
  brts, tree_height, postorderedges, logdexp, logdunifU, logdunif, llrdexp_x,
  logdbeta, llrdbeta_x, logdnorm, logdnorm_tc, llrdnorm_ωxm, llrdnorm_σ², 
  llrdnorm_μ, llrdnorm_x, llrdnorm_xμ, logdtnorm, llrdtnorm_x, erf_custom,
  logdhcau, logdhcau1

end # module Utils
