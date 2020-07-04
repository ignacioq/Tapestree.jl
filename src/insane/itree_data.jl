#=

insane tree data gatherers

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#



"""
    treelength(tree::itree)

Return the branch length sum of `tree`.
"""
treelength(tree::itree) = treelength(tree.d1) + treelength(tree.d2) + tree.pe

"""
    treelength(::Nothing)

Return the branch length sum of `tree`.
"""
treelength(::Nothing) = 0.0




"""
    subnn(tree::itree)

Return the number of descendant nodes of `tree`.
"""
subnn(tree::itree) = subnn(tree.d1) + subnn(tree.d2) + 1

"""
    subnn(::Nothing)

Return the number of descendant nodes of `tree`.
"""
subnn(::Nothing) = 0





"""
    treeheight(tree::itree)

Return the number of descendant nodes of `tree`.
"""
treeheight(tree::itree) = 


"""
    subnn(::Nothing)

Return the number of descendant nodes of `tree`.
"""
treeheight(::Nothing) = 0.0


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



