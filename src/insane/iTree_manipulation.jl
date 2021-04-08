#=

insane tree manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    swapbranch!(tree::iTgbmbd,
                nbtr::iTgbmbd,
                dri ::BitArray{1}, 
                ldr ::Int64,
                it  ::Bool,
                ix  ::Int64)

Swap branch given by `dri` by `nbtr` and return the tree.
"""
function swapbranch!(tree::T,
                     nbtr::T,
                     dri ::BitArray{1}, 
                     ldr ::Int64,
                     it  ::Bool,
                     ix  ::Int64) where {T <: iTree}

  if ix === ldr
    if !it
      nbtr = addtree(nbtr, tree) 
    end
    return nbtr
  elseif ix < ldr
    ifx1 = isfix(tree.d1::T)
    if ifx1 && isfix(tree.d2::T)
      ix += 1
      if dri[ix]
        tree.d1 = 
          swapbranch!(tree.d1::T, nbtr::T, dri, ldr, it, ix)
      else
        tree.d2 = 
          swapbranch!(tree.d2::T, nbtr::T, dri, ldr, it, ix)
      end
    elseif ifx1
      tree.d1 = 
          swapbranch!(tree.d1::T, nbtr::T, dri, ldr, it, ix)
    else
      tree.d2 = 
          swapbranch!(tree.d2::T, nbtr::T, dri, ldr, it, ix)
    end
  end

  return tree
end




"""
    addtree(tree::sTbd, dtree::sTbd) 

Add `dtree` to not-extinct tip in `tree` as speciation event, making
sure that the daughters of `dtree` are fixed.
"""
function addtree(tree::T, dtree::T) where {T <: iTree}

  if istip(tree::T) && !isextinct(tree::T)

    dtree = fixds(dtree)

    tree.d1 = dtree.d1
    tree.d2 = dtree.d2

    return tree
  end

  if !isnothing(tree.d1)
    tree.d1 = addtree(tree.d1::T, dtree::T)
  end
  if !isnothing(tree.d2)
    tree.d2 = addtree(tree.d2::T, dtree::T)
  end

  return tree
end




"""
    fixds(tree::T)

Returns the first tree with both daughters fixed.
"""
function fixds(tree::T) where {T <: iTree}

  ifx1 = isfix(tree.d1::T)
  if ifx1 && isfix(tree.d2::T)
    return tree
  elseif ifx1
    tree = fixds(tree.d1::T)
  else
    tree = fixds(tree.d2::T)
  end

  return tree
end




"""
    gbm_copy!(treec::iTgbmbd,
              treep::iTgbmbd)

Copies the gbm birth-death in place for a fixed branch.
"""
function gbm_copy_f!(treec::iTgbmbd,
                     treep::iTgbmbd)

  copyto!(lλ(treec), lλ(treep))
  copyto!(lμ(treec), lμ(treep))

  if !istip(treec)
    ifx1 = isfix(treec.d1::iTgbmbd)
    if ifx1 && isfix(treec.d2::iTgbmbd)
      return nothing
    elseif ifx1
      gbm_copy_f!(treec.d1::iTgbmbd, treep.d1::iTgbmbd)
      gbm_copy!(  treec.d2::iTgbmbd, treep.d2::iTgbmbd)
    else
      gbm_copy!(  treec.d1::iTgbmbd, treep.d1::iTgbmbd)
      gbm_copy_f!(treec.d2::iTgbmbd, treep.d2::iTgbmbd)
    end
  end

  return nothing
end




"""
    gbm_copy!(treec::iTgbmbd,
              treep::iTgbmbd)

Copies the gbm birth-death in place.
"""
function gbm_copy!(treec::iTgbmbd,
                   treep::iTgbmbd)

  copyto!(lλ(treec), lλ(treep))
  copyto!(lμ(treec), lμ(treep))

  if !istip(treec)
    gbm_copy!(treec.d1::iTgbmbd, treep.d1::iTgbmbd)
    gbm_copy!(treec.d2::iTgbmbd, treep.d2::iTgbmbd)
  end

  return nothing
end





"""
    graftree!(tree ::T,
              stree::T,
              dri  ::BitArray{1},
              h    ::Float64,
              ldr  ::Int64,
              thc  ::Float64,
              ix   ::Int64) where {T <: iTree}

Graft `stree` into `tree` given the address `idr`.
"""
function graftree!(tree ::T,
                   stree::T,
                   dri  ::BitArray{1},
                   h    ::Float64,
                   ldr  ::Int64,
                   thc  ::Float64,
                   ix   ::Int64) where {T <: iTree}

  if ix == ldr 
    if thc > h > (thc - pe(tree))
      npe = thc - h
      addpe!(tree, -npe)
      tree = rand() <= 0.5 ? T(tree, stree, npe, false, true) :
                             T(stree, tree, npe, false, true)
    else
      if isfix(tree.d1)
        tree.d1 = 
          graftree!(tree.d1::T, stree, dri, h, ldr, thc - pe(tree), ix)
      else
        tree.d2 =
          graftree!(tree.d2::T, stree, dri, h, ldr, thc - pe(tree), ix)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      ix += 1
      if dri[ix]
        tree.d1 = 
          graftree!(tree.d1::T, stree, dri, h, ldr, thc - pe(tree), ix)
      else
        tree.d2 = 
          graftree!(tree.d2::T, stree, dri, h, ldr, thc - pe(tree), ix)
      end
    elseif ifx1
      tree.d1 = 
        graftree!(tree.d1::T, stree, dri, h, ldr, thc - pe(tree), ix)
    else
      tree.d2 =
        graftree!(tree.d2::T, stree, dri, h, ldr, thc - pe(tree), ix)
    end
  end

  return tree
end




"""
    prunetree!(tree::T, 
               dri ::BitArray{1}, 
               ldr ::Int64,
               wpr ::Int64,
               ix  ::Int64, 
               px  ::Int64) where {T <: iTree}

Prune tree at branch given by `dri` and grafted `wpr`.
"""
function prunetree!(tree::T, 
                    dri ::BitArray{1}, 
                    ldr ::Int64,
                    wpr ::Int64,
                    ix  ::Int64, 
                    px  ::Int64) where {T <: iTree}

  if ix == ldr
    if px == wpr
      if isfix(tree.d1::T)
        npe  = pe(tree) + pe(tree.d1)
        setpe!(tree.d1, npe)
        tree = tree.d1
      elseif isfix(tree.d2::T)
        npe  = pe(tree) + pe(tree.d2)
        setpe!(tree.d2, npe)
        tree = tree.d2
      end
    else
      px += 1
      if isfix(tree.d1::T)
        tree.d1 = 
          prunetree!(tree.d1::T, dri, ldr, wpr, ix, px)
      else
        tree.d2 =
          prunetree!(tree.d2::T, dri, ldr, wpr, ix, px)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::T)
    if ifx1 && isfix(tree.d2::T)
      ix += 1
      if dri[ix]
        tree.d1 = 
          prunetree!(tree.d1::T, dri, ldr, wpr, ix, px)
      else
        tree.d2 = 
          prunetree!(tree.d2::T, dri, ldr, wpr, ix, px)
      end
    elseif ifx1
      tree.d1 = 
        prunetree!(tree.d1::T, dri, ldr, wpr, ix, px)
    else
      tree.d2 =
        prunetree!(tree.d2::T, dri, ldr, wpr, ix, px)
    end
  end

  return tree
end




"""
    remove_extinct(tree::T) where {T <: iTree}

Remove extinct tips from `iTree`.
"""
function remove_extinct(tree::T) where {T <: iTree}

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
    fixtree!(tree::T) where {T <: iTree}

Fix all `tree`.
"""
function fixtree!(tree::T) where {T <: iTree}
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
  fix!(tree::T) where {T <: iTree}

Fix `tree`.
"""
fix!(tree::T) where {T <: iTree} = setproperty!(tree, :fx, true)

"""
  fix!(::Nothing)

Fix `tree`.
"""
fix!(::Nothing) = nothing




"""
  setpe!(tree::T, pe::Float64) where {T <: iTree}

Set pendant edge for `tree`.
"""
setpe!(tree::T, pe::Float64) where {T <: iTree} = setproperty!(tree, :pe, pe)




"""
  addpe!(tree::T, pe::Float64) where {T <: iTree}

Add `pe` to pendant edge of `tree`.
"""
addpe!(tree::T, pe::Float64) where {T <: iTree} = tree.pe += pe




"""
  setd1!(tree::T,  stree::T) where {T <: iTree}

Set `d1` to `stree` in `tree`.
"""
setd1!(tree::T,  stree::T) where {T <: iTree} = setproperty!(tree, :d1, stree)

setd1!(tree::T,  ::Nothing) where {T <: iTree} = 
  setproperty!(tree, :d1, nothing)




"""
  setd2!(tree::T,  stree::T) where {T <: iTree}

Set `d2` to `stree` in `tree`.
"""
setd2!(tree::T,  stree::T) where {T <: iTree} = setproperty!(tree, :d2, stree)

setd2!(tree::T,  ::Nothing) where {T <: iTree} = 
  setproperty!(tree, :d2, nothing)


