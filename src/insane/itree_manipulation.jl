#=

insane tree manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#



"""
    graft!(tree ::itree, 
           stree::itree, 
           p    ::Int64,
           d    ::Int64,
           btime::Float64)

Graft stree into tree at time `btime` on `Edge(p,d)`.
"""
function graft!(tree ::itree, 
                stree::itree, 
                p    ::Int64,
                d    ::Int64,
                btime::Float64)

  # add node and subdivide branch lengths
  add_node!(tree)
  nnt = nn(tree)
  add_branch!(tree, p, nnt, btime)
  add_branch!(tree, nnt, d, blength(tree, p, d) - btime)

  # remove branch
  rm_branch!(tree, p, d)

  # merge into tree
  merge_itree!(tree, stree)

  # set numerics
  sum_nλ!(tree, nλ(stree)) 
  sum_nμ!(tree, nμ(stree)) 

  return nothing
end




"""
    merge_itree!(tree ::itree, stree::itree)

Merge stree into tree from the last node added to `tree`, which is the
one to merge the graphs with.
"""
function merge_itree!(tree ::itree,
                      stree::itree)
  nnt = nn(tree) - 1

  add_nodes!(tree, nn(stree))
  for b in branches(stree)
    p = b.src
    d = b.dst
    add_branch!(tree, p + nnt, d + nnt, blength(stree, p, d))
  end

  return nothing
end





"""
    set_nλ!(tree::itree, n::Int64)

Set number of speciation events.
"""
set_nλ!(tree::itree, n::Int64) = getproperty(tree,:nλ)[]  = n

"""
    set_nλ!(tree::itree, n::Int64)

Add to the number of speciation events.
"""
sum_nλ!(tree::itree, n::Int64) = getproperty(tree,:nλ)[] += n




"""
    set_nμ!(tree::itree, n::Int64)

Set number of extinction events.
"""
set_nμ!(tree::itree, n::Int64) = getproperty(tree,:nμ)[] = n

"""
    sum_nμ!(tree::itree, n::Int64)

Set number of extinction events.
"""
sum_nμ!(tree::itree, n::Int64) = getproperty(tree,:nμ)[] += n



"""
    nt(tree::itree)

Set number of extant tips.
"""
set_nt!(tree::itree, n::Int64) = getproperty(tree,:nt)[] = n




"""
    th(tree::itree)

Set tree height.
"""
set_th!(tree::itree, h::Float64) = getproperty(tree,:th)[] = h






"""
    add_branch!(tree::itree, p::Int64, d::Int64, el::Float64)

Creates a branch from parent node `p` to daughter node `d` with length `bl`.
"""
function add_branch!(tree::itree, p::Int64, d::Int64, bl::Float64)
  LightGraphs.add_edge!(tree.dg, p, d)
  push!(tree.bl, LightGraphs.Edge(p, d) => bl)
  return nothing
end





"""
    rm_branch!(tree::itree, p::Int64, d::Int64)

Removes the edge from parent node `p` to daughter node `d`.
"""
function rm_branch!(tree::itree, p::Int64, d::Int64)
  LightGraphs.rem_edge!(tree.dg, p, d)
  delete!(tree.bl, LightGraphs.Edge(p, d))
  return nothing
end




"""
    add_node!(tree::itree)

Add one node.
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





