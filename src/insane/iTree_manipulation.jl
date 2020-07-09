#=

insane tree manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#



setpe!(tree::iTree, pe::Float64) = tree.pe = pe

# remove extinct tips



t.d1

isextinct(t.d2)

remove_extinct(::Nothing) = nothing

function remove_extinct(tree::iTree)


  if isextinct(tree.d1)
    npe  = pe(tree) + pe(tree.d2)
    tree = tree.d2
    setpe!(tree, npe)
  elseif isextinct(tree.d2)
    npe  = pe(tree) + pe(tree.d1)
    tree = tree.d1
    setpe!(tree, npe)
  end

  tree.d1 = remove_extinct(tree.d1)
  tree.d2 = remove_extinct(tree.d2)

  return tree
end

"""
deletes in place 
"""

tree = deepcopy(t)

te = remove_extinct(deepcopy(tree))

plot(t)
plot(te)

snn(tree)
treeheight(tree)

snn(t)
treeheight(t)
snen(t)

plot(tree)

function sim_cpb(t::Float64, λ::Float64)

  tw = cpb_wait(λ)

  if tw > t
    return iTree(t)
  end

  iTree(sim_cpb(t - tw, λ), sim_cpb(t - tw, λ), tw)
end

