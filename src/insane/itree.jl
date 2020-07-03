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
# struct itree
#   dg::SimpleDiGraph{Int64}
#   bl::Dict{LightGraphs.Edge{Int64}, Float64}
#   th::Base.RefValue{Float64}
#   nλ::Base.RefValue{Int64}
#   nμ::Base.RefValue{Int64}
#   nt::Base.RefValue{Int64}

#   # bare bones constructor
#   itree() = 
#     new(SimpleDiGraph{Int64}(),
#         Dict{LightGraphs.Edge{Int64}, Float64}(),
#         Ref(0.0), Ref(0), Ref(0), Ref(0))
# end




mutable struct itree
  d1::itree
  d2::itree
  l ::Float64

  # bare bones constructor
  itree() = (x = new(); x.d1 = x; x.d2 = x)
end



t = itree(itree(), itree(), 10.)
t.l = 1.0

t.d1 = 


t = itree(itree(), itree(), 0.5)




mutable struct tree
  d1::Union{tree, Nothing}
  d2::Union{tree, Nothing}
  l ::Float64

  # constructors
  tree() = new(nothing, nothing, 0.0)
  tree(l::Float64) = new(nothing, nothing, l)
  tree(d1::tree, d2::tree, l::Float64) = new(d1, d2, l)
end



#make tree
t = tree(tree(tree(0.6),tree(0.6), 0.4), tree(1.0), 0.0)

st = tree(tree(0.2),tree(0.2), 0.1)



