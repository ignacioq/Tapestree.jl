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
el: edge lengths
th: tree height
nλ: number of speciation events
nμ: number of extinction events
nt: number of extant tips


    itree()

Constructs an empty `itree` object.
"""
struct itree
  dg::SimpleDiGraph{Int64}
  el::Dict{LightGraphs.SimpleGraphs.SimpleEdge{Int64}, Float64}
  th::Base.RefValue{Float64}
  nλ::Base.RefValue{Int64}
  nμ::Base.RefValue{Int64}
  nt::Base.RefValue{Int64}

  itree() = 
    new(SimpleDiGraph{Int64}(),
        Dict{LightGraphs.SimpleGraphs.SimpleEdge{Int64}, Float64}(),
        Ref(0.0), Ref(0), Ref(0), Ref(0))

  Base.show(io::IO, tree::itree) = print(io, 
    "insane phylogenetic tree with ", tree.nt[] + tree.nμ[], " tips and ", tree.nλ[], " internal nodes")
end





"""
    nλ(tree::itree)

Return number of speciation events for an `itree` tree.
"""
nλ(tree::itree) = getproperty(tree, :nλ)[]





"""
    nμ(tree::itree)
Return number of extinction events for an `itree` tree.
"""
nμ(tree::itree) = getproperty(tree, :nμ)[]





"""
    nt(tree::itree)
Return number of extant tips for an `itree` tree.
"""
nt(tree::itree) = getproperty(tree, :ne)[]





"""
    th(tree::itree)

Return the tree height for an `itree` tree.
"""
th(tree::itree) = getproperty(tree, :th)[]





"""
    nλ(tree::itree)
"""
set_nλ!(tree::itree, n::Int64)   = getproperty(tree,:nλ)[] = n





"""
    nμ(tree::itree)
"""
set_nμ!(tree::itree, n::Int64)   = getproperty(tree,:nμ)[] = n





"""
    nt(tree::itree)
"""
set_nt!(tree::itree, n::Int64)   = getproperty(tree,:nt)[] = n





"""
    th(tree::itree)
"""
set_th!(tree::itree, h::Float64) = getproperty(tree,:th)[] = h





"""
    add_edge!(tree::itree, p::Int64, d::Int64, el::Float64)

Creates an edge from parent node `p` to daughter node `d` with length `el`.
"""
function add_edge!(tree::itree, p::Int64, d::Int64, el::Float64)
  LightGraphs.add_edge!(tree.dg, p, d)
  push!(tree.el, LightGraphs.SimpleGraphs.SimpleEdge(p, d) => el)
  return nothing
end





"""
    rm_edge!(tree::itree, p::Int64, d::Int64)

Removes the edge from parent node `p` to daughter node `d`.
"""
function rm_edge!(tree::itree, p::Int64, d::Int64)
  LightGraphs.rem_edge!(tree.dg, p, d)
  delete!(tree.el, LightGraphs.SimpleGraphs.SimpleEdge(p, d))
  return nothing
end





"""
    add_node!(tree::itree)

Add a node.
"""
function add_node!(tree::itree)
  LightGraphs.add_vertex!(tree.dg)
  return nothing
end





"""
    add_nodes!(tree::itree, nn::Int64)

Add `nn` nodes.
"""
function add_nodes!(tree::itree, nn::Int64)
  LightGraphs.add_vertices!(tree.dg, nn)
  return nothing
end





"""
    rm_node!(tree::itree, n::Int64)

Add a node.
"""
function rm_node!(tree::itree, n::Int64)
  LightGraphs.rem_vertex!(tree.dg, n)
  return nothing
end





"""
    pnode(tree::itree, n::Int64)

Retrieve parent node.
"""
pnode(tree::itree, n::Int64) = LightGraphs.inneighbors(tree.dg, n)





"""
    dnodes(tree::itree, n::Int64)

Retrieve daughter nodes.
"""
dnodes(tree::itree, n::Int64) = LightGraphs.outneighbors(tree.dg, n)





t = itree()
add_nodes!(t, 5)

add_edge!(t, 1, 2, 1.0)
add_edge!(t, 1, 3, 0.5)
add_edge!(t, 3, 4, 0.5)
add_edge!(t, 3, 5, 0.5)


pnode(t, 2)

dnodes(t, 3)

