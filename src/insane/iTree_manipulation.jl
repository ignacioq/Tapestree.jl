#=

insane tree manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#


"""
  setpe!(tree::iTree, pe::Float64)

Set pendant edge for `tree`.
"""
setpe!(tree::iTree, pe::Float64) = setproperty!(tree, :pe, pe)



"""
  addpe!(tree::iTree, pe::Float64)

Add `pe` to pendant edge of `tree`.
"""
addpe!(tree::iTree, pe::Float64) = tree.pe += pe




"""
    remove_extinct(tree::iTree)

Remove extinct tips from `iTree`.
WARNING: it changes the `tree` object.
"""
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
    remove_extinct(::Nothing)

Remove extinct tips from `iTree`.
WARNING: it changes the `tree` object.
"""
remove_extinct(::Nothing) = nothing

