#=

  "INSANE: joINt Speciation And Niche Evolution"

=#

module INSANE

using LightGraphs: SimpleDiGraph, SimpleGraphs.SimpleEdge, add_edge!, rem_edge!, add_vertex!, add_vertices!, rem_vertex!, inneighbors, outneighbors 
using Random: randexp

# other submodules dependencies
using ..Utils

# files
include("insane/itree.jl")

end # module INSANE