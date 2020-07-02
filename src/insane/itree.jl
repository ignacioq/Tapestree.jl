#=

insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#





"""
    itree

An immutable Composite type representing a phylogenetic tree using
graph theory for `insane` use, with the following fields:

dg: Directed Graph (nodes and edges configuration)
bl: branch lengths
th: tree height
nλ: number of speciation events
nμ: number of extinction events
nt: number of extant tips


    itree()

Constructs an empty `itree` object.
"""
struct itree
  dg::SimpleDiGraph{Int64}
  bl::Dict{LightGraphs.Edge{Int64}, Float64}
  th::Base.RefValue{Float64}
  nλ::Base.RefValue{Int64}
  nμ::Base.RefValue{Int64}
  nt::Base.RefValue{Int64}

  itree() = 
    new(SimpleDiGraph{Int64}(),
        Dict{LightGraphs.Edge{Int64}, Float64}(),
        Ref(0.0), Ref(0), Ref(0), Ref(0))

  Base.show(io::IO, tree::itree) = print(io, 
    "insane phylogenetic tree: ", 
    tree.nt[], " extant + ", 
    tree.nμ[], " extinct tips, ", 
    tree.nλ[], " internal nodes")
end






