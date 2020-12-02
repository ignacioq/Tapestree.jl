#=

insane tree data gatherers

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#



"""
    isfix(tree::T) where {T <: iTree} 

Return if is either an extant or extinct tip node.
"""
isfix(tree::T) where {T <: iTree} = getproperty(tree,:fx)

isfix(::Nothing) = true




"""
    istip(tree::T) where {T <: iTree}

Return if is either an extant or extinct tip node.
"""
istip(tree::T) where {T <: iTree} = isnothing(tree.d1) && isnothing(tree.d2)

istip(::Nothing) = false




"""
    isextinct(tree::T) where {T <: iTree}

Return if is an extinction node.
"""
isextinct(tree::T) where {T <: iTree} = getproperty(tree,:iμ)

"""
    isextinct(::Nothing)

Return if is an extinction node.
"""
isextinct(::Nothing) = false




"""
    pe(tree::T) where {T <: iTree}

Return pendant edge.
"""
pe(tree::T) where {T <: iTree} = getproperty(tree,:pe)

"""
    pe(::Nothing)

Return pendant edge.
"""
pe(::Nothing) = 0.0




"""
    treelength(tree::T) where {T <: iTree}

Return the branch length sum of `tree`.
"""
treelength(tree::T) where {T <: iTree} = treelength(tree.d1) + treelength(tree.d2) + pe(tree)

"""
    treelength(::Nothing)

Return the branch length sum of `tree`.
"""
treelength(::Nothing) = 0.0




"""
    treeheight(tree::T) where {T <: iTree}

Return the tree height of `tree`.
"""
function treeheight(tree::T) where {T <: iTree}
  th1 = treeheight(tree.d1)
  th2 = treeheight(tree.d2)
  (th1 > th2 ? th1 : th2) + pe(tree)
end

"""
    treeheight(::Nothing)

Return the tree height of `tree`.
"""
treeheight(::Nothing) = 0.0




"""
    snn(tree::T) where {T <: iTree}

Return the number of descendant nodes for `tree`.
"""
snn(tree::T) where {T <: iTree} = snn(tree.d1) + snn(tree.d2) + 1

"""
    snn(::Nothing)

Return the number of descendant nodes for `tree`.
"""
snn(::Nothing) = 0




"""
    snin(tree::T) where {T <: iTree}

Return the number of internal nodes for `tree`.
"""
function snin(tree::T) where {T <: iTree}
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
    sntn(tree::T) where {T <: iTree}

Return the number of tip nodes for `tree`.
"""
function sntn(tree::T) where {T <: iTree}
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




"""
    snen(tree::T) where {T <: iTree}

Return the number of extinct tip nodes for `tree`.
"""
function snen(tree::T) where {T <: iTree}
    if isextinct(tree)
      return 1
    else
      return snen(tree.d1) + snen(tree.d2)
    end
end


"""
    snen(::Nothing)

Return the number of extinct tip nodes for `tree`.
"""
snen(::Nothing) = 0





"""
    streeheight(tree::T,
                h   ::Float64, 
                th  ::Float64,
                dri ::BitArray{1}, 
                ldr ::Int64,
                wpr ::Int64,
                ix  ::Int64, 
                px  ::Int64) where {T <: iTree}

Return the height of a grafted subtree in `tree`. 
"""
function streeheight(tree::T,
                     h   ::Float64, 
                     th  ::Float64,
                     dri ::BitArray{1}, 
                     ldr ::Int64,
                     wpr ::Int64,
                     ix  ::Int64, 
                     px  ::Int64) where {T <: iTree}

  if ix == ldr
    if px == wpr
      if isfix(tree.d1::T)
        return h - pe(tree), treeheight(tree.d2)
      elseif isfix(tree.d2::T)
        return h - pe(tree), treeheight(tree.d1)
      end
    else
      px += 1
      if isfix(tree.d1::T)
        h, th = 
          streeheight(tree.d1::T, h - pe(tree), th, dri, ldr, wpr, ix, px)
      else
        h, th =
          streeheight(tree.d2::T, h - pe(tree), th, dri, ldr, wpr, ix, px)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::T)
    if ifx1 && isfix(tree.d2::T)
      ix += 1
      if dri[ix]
        h, th = 
          streeheight(tree.d1::T, h - pe(tree), th, dri, ldr, wpr, ix, px)
      else
        h, th = 
          streeheight(tree.d2::T, h - pe(tree), th, dri, ldr, wpr, ix, px)
      end
    elseif ifx1
      h, th = 
        streeheight(tree.d1::T, h - pe(tree), th, dri, ldr, wpr, ix, px)
    else
      h, th =
        streeheight(tree.d2::T, h - pe(tree), th, dri, ldr, wpr, ix, px)
    end
  end
end




"""
    ts(tree::T) where {T <: iTgbm} 

Return pendant edge.
"""
ts(tree::T) where {T <: iTgbm} = getproperty(tree,:ts)

"""
    pe(tree::iTree)

Return pendant edge.
"""
ts(::Nothing) = 0.0




"""
    lλ(tree::T) where {T <: iTgbm}

Return pendant edge.
"""
lλ(tree::T) where {T <: iTgbm} = getproperty(tree,:lλ)

"""
    pe(tree::iTree)

Return pendant edge.
"""
lλ(::Nothing) = 0.0




"""
    lμ(tree::iTgbm)

Return pendant edge.
"""
lμ(tree::iTgbmbd) = getproperty(tree,:lμ)

"""
    pe(tree::iTree)

Return pendant edge.
"""
lμ(::Nothing) = 0.0
