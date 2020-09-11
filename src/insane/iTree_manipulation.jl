#=

insane tree manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




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


