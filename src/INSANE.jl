#=

  "INSANE: joINt Speciation And Niche Evolution"

=#

module INSANE

using RecipesBase: @recipe
using Random: randexp
using DelimitedFiles: writelm

# other submodules dependencies
using ..Utils

# files
include("insane/itree.jl")

end # module INSANE