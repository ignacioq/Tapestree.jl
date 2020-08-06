#=

insane tree manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    randbranch(th   ::Float64,
               sth  ::Float64,
               idv  ::Array{iDir,1},
               wbr  ::BitArray{1})

Randomly select branch to graft a subtree with height `sth` onto tree 
with height `th`.
"""
function randbranch(th   ::Float64,
                    sth  ::Float64,
                    idv  ::Array{iDir,1},
                    wbr  ::BitArray{1})

  # uniformly set the height at which to insert `stree`
  h = sth + rand()*(th - sth)

  # get which branches are cut by `h` and sample one at random
  nb  = branchescut!(wbr, h, idv)
  rn  = rand(Base.OneTo(nb))
  br  = getbranch(wbr, rn, idv)
  dri = dr(br)
  ldr = lastindex(dri)

  return dri, ldr
end




"""
    graftree!(tree ::iTree,
              stree::iTree,
              dri  ::BitArray{1},
              h    ::Float64,
              ldr  ::Int64,
              thc  ::Float64;
              ix   ::Int64 = 0)

Graft `stree` into `tree` given the address `idr`.
"""
function graftree!(tree ::iTree,
                   stree::iTree,
                   dri  ::BitArray{1},
                   h    ::Float64,
                   ldr  ::Int64,
                   thc  ::Float64,
                   ix   ::Int64)

  if ix == ldr 
    if thc > h > (thc - pe(tree))
      npe = thc - h
      addpe!(tree, -npe)
      tree = rand() <= 0.5 ? iTree(tree, stree, npe, false, true) :
                             iTree(stree, tree, npe, false, true)
    else
      if isfix(tree.d1)
        tree.d1 = 
          graftree!(tree.d1, stree, dri, h, ldr, thc - pe(tree), ix)
      else
        tree.d2 =
          graftree!(tree.d2, stree, dri, h, ldr, thc - pe(tree), ix)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      ix += 1
      if dri[ix]
        tree.d1 = 
          graftree!(tree.d1, stree, dri, h, ldr, thc - pe(tree), ix)
      else
        tree.d2 = 
          graftree!(tree.d2, stree, dri, h, ldr, thc - pe(tree), ix)
      end
    elseif ifx1
      tree.d1 = 
        graftree!(tree.d1, stree, dri, h, ldr, thc - pe(tree), ix)
    else
      tree.d2 =
        graftree!(tree.d2, stree, dri, h, ldr, thc - pe(tree), ix)
    end
  end

  return tree
end





"""
    remove_extinct(tree::iTree)

Remove extinct tips from `iTree`.
"""
function remove_extinct(tree::iTree)

  tree.d1 = remove_extinct(tree.d1)
  tree.d2 = remove_extinct(tree.d2)

  if isextinct(tree.d1)
    npe  = pe(tree) + pe(tree.d2)
    tree = tree.d2
    setpe!(tree, npe)
  elseif isextinct(tree.d2)
    npe  = pe(tree) + pe(tree.d1)
    tree = tree.d1
    setpe!(tree, npe)
  end

  return tree
end

"""
    remove_extinct(::Nothing)

Remove extinct tips from `iTree`.
"""
remove_extinct(::Nothing) = nothing



"""
    fixtree!(tree::iTree)

Fix all `tree`.
"""
function fixtree!(tree::iTree)
  fix!(tree)
  fixtree!(tree.d1)
  fixtree!(tree.d2)
end

"""
    fixtree!(::Nothing)

Fix all `tree`.
"""
fixtree!(::Nothing) = nothing




"""
  fix!(tree::iTree)

Fix `tree`.
"""
fix!(tree::iTree) = setproperty!(tree, :fx, true)

"""
  fix!(::Nothing)

Fix `tree`.
"""
fix!(::Nothing) = nothing




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
  setd1!(tree::iTree, stree::iTree)

Set `d1` to `stree` in `tree`.
"""
setd1!(tree::iTree,  stree::iTree) = setproperty!(tree, :d1, stree)




"""
  setd2!(tree::iTree, stree::iTree)

Set `d2` to `stree` in `tree`.
"""
setd2!(tree::iTree,  stree::iTree) = setproperty!(tree, :d2, stree)



