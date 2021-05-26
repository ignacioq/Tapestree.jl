#=

insane tree manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#





"""
    swapbranch!(treep::T,
                treec::T,
                nbtr ::T,
                dri  ::BitArray{1}, 
                ldr  ::Int64,
                it   ::Bool,
                ix   ::Int64) where {T <: iTree}

Swap branch given by `dri` by `nbtr` and return the tree.
"""
function swapbranch!(treep::T,
                     treec::T,
                     nbtr ::T,
                     dri  ::BitArray{1}, 
                     ldr  ::Int64,
                     it   ::Bool,
                     ix   ::Int64) where {T <: iTree}

  if ix === ldr
    nbtrp = deepcopy(nbtr)
    if !it
      nbtrp = addtree(nbtrp, treep) 
      nbtr  = addtree(nbtr,  treec) 
    end
    return nbtrp, nbtr
  elseif ix < ldr
    ifx1 = isfix(treec.d1::T)
    if ifx1 && isfix(treec.d2::T)
      ix += 1
      if dri[ix]
        treep.d1, treec.d1 = 
          swapbranch!(treep.d1::T, treec.d1::T, nbtr::T, dri, ldr, it, ix)
      else
        treep.d2, treec.d2 = 
          swapbranch!(treep.d2::T, treec.d2::T, nbtr::T, dri, ldr, it, ix)
      end
    elseif ifx1
      treep.d1, treec.d1 = 
          swapbranch!(treep.d1::T, treec.d1::T, nbtr::T, dri, ldr, it, ix)
    else
      treep.d2, treec.d2 = 
          swapbranch!(treep.d2::T, treec.d2::T, nbtr::T, dri, ldr, it, ix)
    end
  end

  return treep, treec
end





"""
    swapbranch!(tree::T,
                nbtr::T,
                dri ::BitArray{1}, 
                ldr ::Int64,
                it  ::Bool,
                ix  ::Int64) where {T <: iTree}

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

    dtree = fixdstree(dtree)

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
    gbm_copy_f!(tree::iTgbmbd,
                bbλ ::Array{Float64,1},
                bbμ ::Array{Float64,1},
                ii  ::Int64)

Copies the gbm birth-death in place for a fixed branch into the 
help arrays `bbλ` and `bbμ`.
"""
function gbm_copy_f!(tree::iTgbmbd,
                     bbλ ::Array{Float64,1},
                     bbμ ::Array{Float64,1},
                     ii  ::Int64)

  lλv = lλ(tree)
  lμv = lμ(tree)
  lt  = lastindex(lλv)

  @simd for i in Base.OneTo(lt)
    ii     += 1
    bbλ[ii] = lλv[i]
    bbμ[ii] = lμv[i]
  end

  if !istip(tree)
    ifx1 = isfix(tree.d1::iTgbmbd)
    if ifx1 && isfix(tree.d2::iTgbmbd)
      return nothing
    elseif ifx1
      gbm_copy_f!(tree.d1::iTgbmbd, bbλ, bbμ, ii-1)
    else
      gbm_copy_f!(tree.d2::iTgbmbd, bbλ, bbμ, ii-1)
    end
  end

  return nothing
end





"""
    gbm_copy_f!(treec::iTgbmbd,
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
    fixalive!(tree::T) where T <: iTree

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::T) where T <: iTree

  if istip(tree::T) && !isextinct(tree::T)
    fix!(tree::T)
    return true
  end

  if !isnothing(tree.d2)
    f = fixalive!(tree.d2::T)
    if f 
      fix!(tree)
      return true
    end
  end

  if !isnothing(tree.d1)
    f = fixalive!(tree.d1::T)
    if f 
      fix!(tree)
      return true
    end
  end

  return false
end





"""
    fixrtip!(tree::T, na::Int64) where T <: iTree

Fixes the the path for a random non extinct tip.
"""
function fixrtip!(tree::T, na::Int64) where T <: iTree

  fix!(tree)

  if !istip(tree)
    if isextinct(tree.d1::T)
      fixrtip!(tree.d2::T, na)
    elseif isextinct(tree.d2::T)
      fixrtip!(tree.d1::T, na)
    else
      na1 = snan(tree.d1::T)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        fixrtip!(tree.d2::T, na - na1)
      else
        fixrtip!(tree.d1::T, na1)
      end
    end
  end
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
    remove_extinct(tree::iTgbmbd)

Remove extinct tips from `iTgbmbd`.
"""
function remove_extinct(tree::iTgbmbd)

  tree.d1 = remove_extinct(tree.d1)
  tree.d2 = remove_extinct(tree.d2)

  if isextinct(tree.d1)
    if isextinct(tree.d2)
      tree.d1 = nothing
      tree.d2 = nothing
      setproperty!(tree, :iμ, true)
    else
      ppr = pe(tree)
      npe = ppr + pe(tree.d2)
      ts2 = ts(tree.d2)
      ls  = lastindex(ts2)

      @simd for i in Base.OneTo(ls) 
        ts2[i] += ppr
      end

      ts0 = ts(tree)
      lλ0 = lλ(tree)
      lμ0 = lμ(tree)

      pop!(ts0)
      pop!(lλ0)
      pop!(lμ0)

      prepend!(ts2, ts0) 
      prepend!(lλ(tree.d2), lλ0) 
      prepend!(lμ(tree.d2), lμ0)

      tree = tree.d2
      setpe!(tree, npe)

    end
  elseif isextinct(tree.d2)
    ppr = pe(tree)
    npe = ppr + pe(tree.d1)

    ts1 = ts(tree.d1)

    @simd for i in Base.OneTo(lastindex(ts1)) 
      ts1[i] += ppr
    end

    ts0 = ts(tree)
    lλ0 = lλ(tree)
    lμ0 = lμ(tree)

    pop!(ts0)
    pop!(lλ0)
    pop!(lμ0)

    prepend!(ts1, ts0) 
    prepend!(lλ(tree.d1), lλ0) 
    prepend!(lμ(tree.d1), lμ0)

    tree = tree.d1
    setpe!(tree, npe)
  end

  return tree
end




"""
    remove_extinct(treev::Array{T,1}) where {T <: iTree}

Remove extinct taxa for a vector of trees.
"""
function remove_extinct(treev::Array{T,1}) where {T <: iTree}

  treevne = T[]
  for t in treev
    push!(treevne, remove_extinct(deepcopy(t)))
  end

  return treevne
end




"""
    remove_extinct(tree::T) where {T <: sT}

Remove extinct tips from `sT`.
"""
function remove_extinct(tree::sTbd)

  tree.d1 = remove_extinct(tree.d1)
  tree.d2 = remove_extinct(tree.d2)

  if isextinct(tree.d1)
    if isextinct(tree.d2)
      tree.d1 = nothing
      tree.d2 = nothing
      setproperty!(tree, :iμ, true)
    else
      npe  = pe(tree) + pe(tree.d2)
      tree = tree.d2
      setpe!(tree, npe)
    end
  elseif isextinct(tree.d2)
    npe  = pe(tree) + pe(tree.d1)
    tree = tree.d1
    setpe!(tree, npe)
  end

  return tree
end

"""
    remove_extinct(::Nothing)

Remove extinct tips from `sT`.
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


