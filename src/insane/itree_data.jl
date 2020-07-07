#=

insane tree data gatherers

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    istip(tree::itree)

Return if is either an extant or extinct tip node.
"""
istip(tree::itree) = isnothing(tree.d1) && isnothing(tree.d2)




"""
    isλ(tree::itree)

Return if is an extinction node.
"""
isextinct(tree::itree) = getproperty(tree,:isμ)




"""
    pe(tree::itree)

Return pendant edge.
"""
pe(tree::itree) = getproperty(tree,:pe)

"""
    pe(tree::itree)

Return pendant edge.
"""
pe(::Nothing) = 0.0




"""
    treelength(tree::itree)

Return the branch length sum of `tree`.
"""
treelength(tree::itree) = treelength(tree.d1) + treelength(tree.d2) + pe(tree)

"""
    treelength(::Nothing)

Return the branch length sum of `tree`.
"""
treelength(::Nothing) = 0.0




"""
    treeheight(tree::itree)

Return the tree height of `tree`.
"""
function treeheight(tree::itree)
  th1 = treeheight(tree.d1)
  th2 = treeheight(tree.d2)
  (th1 > th2 ? th1 : th2) + pe(tree)
end

"""
    treeheight(tree::itree)

Return the tree height of `tree`.
"""
treeheight(::Nothing) = 0.0




"""
    snn(tree::itree)

Return the number of descendant nodes for `tree`.
"""
snn(tree::itree) = snn(tree.d1) + snn(tree.d2) + 1

"""
    snn(::Nothing)

Return the number of descendant nodes for `tree`.
"""
snn(::Nothing) = 0




"""
    snin(tree::itree)

Return the number of internal nodes for `tree`.
"""
function snin(tree::itree)
    if istip(tree)
      return 0
    else
      return snin(tree.d1) + snin(tree.d2) + 1
    end
end

"""
    snin(::Nothing)

Return the number of internal nodes for `tree`.
"""
snin(::Nothing) = 0





"""
    sntn(tree::itree)

Return the number of tip nodes for `tree`.
"""
function sntn(tree::itree)
    if istip(tree)
      return 1
    else
      return sntn(tree.d1) + sntn(tree.d2)
    end
end


"""
    sntn(::Nothing)

Return the number of tip nodes for `tree`.
"""
sntn(::Nothing) = 0








# """
#     nλ(tree::itree)

# Return number of speciation events for an `itree` tree.
# """
# nλ(tree::itree) = getproperty(tree, :nλ)[]




# """
#     nμ(tree::itree)
# Return number of extinction events for an `itree` tree.
# """
# nμ(tree::itree) = getproperty(tree, :nμ)[]




# """
#     nt(tree::itree)
# Return number of extant tips for an `itree` tree.
# """
# nt(tree::itree) = getproperty(tree, :ne)[]




# """
#     th(tree::itree)

# Return the tree height for an `itree` tree.
# """
# th(tree::itree) = getproperty(tree, :th)[]





# """
#     nn(tree::itree)

# Return the number of nodes for an `itree` tree.
# """
# nn(tree::itree) = nv(getproperty(tree, :dg))





# """
#     blength(tree::itree, p::Int64, d::Int64)

# Return the branch length of a branch an `itree` tree.
# """
# blength(tree::itree, p::Int64, d::Int64) = tree.bl[Edge(p,d)]




# """
#     pnode(tree::itree, n::Int64)

# Retrieve parent node.
# """
# pnode(tree::itree, n::Int64) = LightGraphs.inneighbors(tree.dg, n)




# """
#     dnodes(tree::itree, n::Int64)

# Retrieve daughter nodes.
# """
# dnodes(tree::itree, n::Int64) = LightGraphs.outneighbors(tree.dg, n)





# """
#     branches(tree::itree)

# Basic iterator over branches.
# """
# branches(tree::itree) = LightGraphs.edges(tree.dg)





# """
#     nodes(tree::itree)

# Basic iterator over nodes.
# """
# nodes(tree::itree) = LightGraphs.vertices(tree.dg)



