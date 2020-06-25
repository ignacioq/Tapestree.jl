#=

  "INSANE: joINt Speciation And Niche Evolution"

=#

module INSANE

using Random: randexp
using DelimitedFiles: readdlm, writedlm
using ProgressMeter: Progress, next!
using Distributed: @sync, @distributed
using SharedArrays: SharedArray

# other submodules dependencies
using ..Utils

# files
include("insane/itree.jl")

end # module INSANE