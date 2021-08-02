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
isfix(tree::T) where {T <: iTree} = getproperty(tree, :fx)




"""
    istip(tree::T) where {T <: iTree}

Return if is either an extant or extinct tip node.
"""
istip(tree::T) where {T <: iTree} = !isdefined(tree, :d1)




"""
    isextinct(tree::T) where {T <: iTree}

Return if is an extinction node.
"""
isextinct(tree::T) where {T <: iTree} = getproperty(tree, :iμ)




"""
    isalive(tree::T) where {T <: iTree}

Return if is an extinction node.
"""
isalive(tree::T) where {T <: iTree} = !getproperty(tree, :iμ)




"""
    e(tree::T) where {T <: iTree}

Return edge length.
"""
e(tree::T) where {T <: iTree} = getproperty(tree, :e)




"""
    dt(tree::T) where {T <: iTree}

Return `δt`.
"""
dt(tree::T) where {T <: iTree} = getproperty(tree, :dt)




"""
    fdt(tree::T) where {T <: iTree}

Return final `δt`.
"""
fdt(tree::T) where {T <: iTree} = getproperty(tree, :fdt)





"""
    treelength(tree::T) where {T <: iTree}

Return the branch length sum of `tree`.
"""
function treelength(tree::T, l::Float64) where {T <: iTree}
  l += e(tree)
  if isdefined(tree, :d1)
    l = treelength(tree.d1, l)::Float64
    l = treelength(tree.d2, l)::Float64
  end
  return l
end




"""
    treeheight(tree::T, th1::Float64, th2::Float64)

Return the tree height of `tree`.
"""
function treeheight(tree::T) where {T <: iTree}
  if isdefined(tree, :d1)
    th1 = treeheight(tree.d1)
    th2 = treeheight(tree.d2)
    return (th1 > th2 ? th1 : th2) + e(tree)
  end
  return e(tree)
end




"""
    snn(tree::T, n::Int64) where {T <: iTree}

Return the number of descendant nodes for `tree`.
"""
function snn(tree::T, n::Int64) where {T <: iTree} 
  n += 1
  if isdefined(tree, :d1)
    n = snn(tree.d1, n)
    n = snn(tree.d2, n)
  end
  
  return n
end




"""
    snin(tree::T, n::Int64) where {T <: iTree}

Return the number of internal nodes for `tree`.
"""
function snin(tree::T, n::Int64) where {T <: iTree}
  if isdefined(tree, :d1)
    n += 1
    n = snin(tree.d1, n)
    n = snin(tree.d2, n)
  end

  return n
end




"""
    snen(tree::T) where {T <: iTree}

Return the number of extinct nodes for `tree`.
"""
function snen(tree::T, n::Int64) where {T <: iTree}
  if isdefined(tree, :d1)
    n = snen(tree.d1, n)
    n = snen(tree.d2, n)
  elseif isextinct(tree)
    n += 1
  end

  return n
end




"""
    snenF(tree::T) where {T <: iTree}

Return the number of extinct nodes for `tree` as Float64.
"""
function snenF(tree::T, n::Int64) where {T <: iTree}
  if isdefined(tree, :d1)
    n = snen(tree.d1, n)
    n = snen(tree.d2, n)
  elseif isextinct(tree)
    n += 1
  end

  return n
end



"""
    snenF(tree::T, n::Float64) where {T <: iTree}

Return the number of extinct nodes for `tree` as Float64.
"""
function snenF(tree::T, n::Float64) where {T <: iTree}

  if isdefined(tree, :d1)
    n = snenF(tree.d1, n)
    n = snenF(tree.d2, n)
  elseif isextinct(tree)
    n += 1.0
  end

  return n
end




"""
    snan(tree::T, n::Int64) where {T <: iTree}

Return the number of alive nodes for `tree`.
"""
function snan(tree::T, n::Int64) where {T <: iTree}

  if isdefined(tree, :d1)
    n = snan(tree.d1, n)
    n = snan(tree.d2, n)
  elseif isalive(tree)
    n += 1
  end

  return n
end




"""
    sntn(tree::T, n::Int64) where {T <: iTree}

Return the number of tip nodes for `tree`.
"""
function sntn(tree::T, n::Int64) where {T <: iTree}

  if isdefined(tree, :d1)
    n = sntn(tree.d1, n)
    n = sntn(tree.d2, n)
  else
    n += 1
  end

  return n
end




"""
    treelength_ne(tree::T, l::Float64, n::Float64)

Return the tree length and Float number of extinct tip nodes for `tree`.
"""
function treelength_ne(tree::T, l::Float64, n::Float64) where {T <: iTree}

  l += e(tree)
  if isdefined(tree, :d1)
    l, n = treelength_ne(tree.d1, l, n)
    l, n = treelength_ne(tree.d2, l, n)
  else
    if isextinct(tree)
      n += 1.0
    end
  end

  return l, n
end




"""
    lλ(tree::T) where {T <: iTgbm}

Return pendant edge.
"""
lλ(tree::T) where {T <: iTgbm} = getproperty(tree, :lλ)




"""
    lμ(tree::iTgbm)

Return pendant edge.
"""
lμ(tree::iTgbmbd) = getproperty(tree,:lμ)




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

  if ix === ldr
    if px === wpr
      if isfix(tree.d1)
        return h - e(tree), treeheight(tree.d2, 0.0, 0.0)
      elseif isfix(tree.d2)
        return h - e(tree), treeheight(tree.d1, 0.0, 0.0)
      end
    else
      px += 1
      if isfix(tree.d1)
        h, th = 
          streeheight(tree.d1, h - e(tree), th, dri, ldr, wpr, ix, px)
      else
        h, th =
          streeheight(tree.d2, h - e(tree), th, dri, ldr, wpr, ix, px)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      ix += 1
      if dri[ix]
        h, th = 
          streeheight(tree.d1, h - e(tree), th, dri, ldr, wpr, ix, px)
      else
        h, th = 
          streeheight(tree.d2, h - e(tree), th, dri, ldr, wpr, ix, px)
      end
    elseif ifx1
      h, th = 
        streeheight(tree.d1, h - e(tree), th, dri, ldr, wpr, ix, px)
    else
      h, th =
        streeheight(tree.d2, h - e(tree), th, dri, ldr, wpr, ix, px)
    end
  end
end





"""
    λμ01(tree::iTgbmbd,
         dri ::BitArray{1}, 
         ldr ::Int64,
         ix  ::Int64,
         λ0  ::Float64,
         μ0  ::Float64)

Return the `λ`,`μ` at the start `0` and end `1` of a fixed branch. 
"""
function λμ01(tree::iTgbmbd,
              dri ::BitArray{1}, 
              ldr ::Int64,
              ix  ::Int64,
              λ0  ::Float64,
              μ0  ::Float64)

  if ix === ldr
    if isnan(λ0)
      λ0 = lλ(tree)[1]
      μ0 = lμ(tree)[1]
    end

    if istip(tree) && islive(tree)
      λ1 = lλ(tree)[end]
      μ1 = lμ(tree)[end]

      return λ0, μ0, λ1, μ1
    else
      ifx1 = isfix(tree.d1::iTgbmbd)
      if ifx1 && isfix(tree.d2::iTgbmbd)
        λ1 = lλ(tree)[end]
        μ1 = lμ(tree)[end]

        return λ0, μ0, λ1, μ1
      elseif ifx1
        λμ01(tree.d1::iTgbmbd, dri, ldr, ix, λ0, μ0)
      else
        λμ01(tree.d2::iTgbmbd, dri, ldr, ix, λ0, μ0)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTgbmbd)
    if ifx1 && isfix(tree.d2::iTgbmbd)
      ix += 1
      if dri[ix]
        λμ01(tree.d1::iTgbmbd, dri, ldr, ix, λ0, μ0)
      else
        λμ01(tree.d2::iTgbmbd, dri, ldr, ix, λ0, μ0)
      end
    elseif ifx1
      λμ01(tree.d1::iTgbmbd, dri, ldr, ix, λ0, μ0)
    else
      λμ01(tree.d2::iTgbmbd, dri, ldr, ix, λ0, μ0)
    end
  end
end





"""
    λμath(tree::iTgbmbd,
          h    ::Float64, 
          th   ::Float64,
          dri  ::BitArray{1}, 
          ldr  ::Int64,
          wpr  ::Int64,
          ix   ::Int64, 
          px   ::Int64)

Return the `λ`,`μ` and time of `tree` at closest to `h`. 
"""
function λμath(tree::iTgbmbd,
               h   ::Float64, 
               th  ::Float64,
               dri ::BitArray{1}, 
               ldr ::Int64,
               ix  ::Int64)

  if ix === ldr
    pei = e(tree)
    if th > h > (th - pei)

      bh  = th - h
      tsi = ts(tree)
      hi  = indmindif_sorted(tsi, bh)
      λh  = lλ(tree)[hi]
      μh  = lμ(tree)[hi]
      nh  = th - tsi[hi]

      return λh, μh, nh, hi
    else
      if isfix(tree.d1::iTgbmbd)
        λμath(tree.d1::iTgbmbd, h, th - pei, dri, ldr, ix)
      else
        λμath(tree.d2::iTgbmbd, h, th - pei, dri, ldr, ix)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTgbmbd)
    if ifx1 && isfix(tree.d2::iTgbmbd)
      ix += 1
      if dri[ix]
        λμath(tree.d1::iTgbmbd, h, th - e(tree), dri, ldr, ix)
      else
        λμath(tree.d2::iTgbmbd, h, th - e(tree), dri, ldr, ix)
      end
    elseif ifx1
      λμath(tree.d1::iTgbmbd, h, th - e(tree), dri, ldr, ix)
    else
      λμath(tree.d2::iTgbmbd, h, th - e(tree), dri, ldr, ix)
    end
  end
end





"""
    ifxe(tree::iTgbmbd)

Return `true` if fix branch goes extinct for GBM forward simulation. 
"""
function ifxe(tree::T) where T <: iTree

  if istip(tree)
    if isextinct(tree)
      return true
    else
      return false
    end
  else
    if isfix(tree.d1)
      ifxe(tree.d1)
    else
      ifxe(tree.d2)
    end
  end
end





"""
    fixds(tree::T) where {T <: iTgbm}

Return the first fixed daughters `d1` and `d2`. Works only for 
internal branches.
"""
function fixds(tree::T) where {T <: iTgbm}
  ifx1 = isfix(tree.d1::T)
  if ifx1 && isfix(tree.d2::T)
    return tree.d1::T, tree.d2::T
  elseif ifx1
    fixds(tree.d1::T)
  else
    fixds(tree.d2::T)
  end
end





"""
    fixdstree(tree::T) where {T <: iTree}

Returns the first tree with both daughters fixed.
"""
function fixdstree(tree::T) where {T <: iTree}

  ifx1 = isfix(tree.d1::T)
  if ifx1 && isfix(tree.d2::T)
    return tree
  elseif ifx1
    tree = fixdstree(tree.d1::T)
  else
    tree = fixdstree(tree.d2::T)
  end

  return tree
end





"""
    fixd1(tree::iTgbmbd)

Return the first fixed daughter `d1`. Works only for internal branches
"""
function fixd1(tree::iTgbmbd)
  ifx1 = isfix(tree.d1::iTgbmbd)
  if ifx1 && isfix(tree.d2::iTgbmbd)
    return tree.d1
  elseif ifx1
    fixd1(tree.d1::iTgbmbd)
  else
    fixd1(tree.d2::iTgbmbd)
  end
end




"""
    fixd2(tree::iTgbmbd)

Return the first fixed daughter `d2`. Works only for internal branches
"""
function fixd2(tree::iTgbmbd)
  ifx1 = isfix(tree.d1::iTgbmbd)
  if ifx1 && isfix(tree.d2::iTgbmbd)
    return tree.d2
  elseif ifx1
    fixd2(tree.d1::iTgbmbd)
  else
    fixd2(tree.d2::iTgbmbd)
  end
end




"""
    drtree(tree::T, dri::BitArray{1}, ldr::Int64, ix::Int64)

Returns the tree given by `dr`. Assuming no unfix branches.
"""
function drtree(tree::T, dri::BitArray{1}, ldr::Int64, ix::Int64) where T <: iTree
  if ldr === ix
    return tree
  elseif ldr > ix
    ix += 1
    if dri[ix]
      drtree(tree.d1::iTgbmbd, dri, ldr, ix)
    else
      drtree(tree.d2::iTgbmbd, dri, ldr, ix)
    end
  end
end



"""
    makebbv!(tree::T, 
             bbλ ::Array{Array{Float64,1},1}, 
             tsv ::Array{Array{Float64,1},1}) where {T <: iTgbm}

Make `bbv` vector with allocated `bb` (brownian bridges) and 
with `tsv` vector of branches times `ts`.
"""
function makebbv!(tree::T, 
                  bbλ ::Array{Array{Float64,1},1}, 
                  tsv ::Array{Array{Float64,1},1}) where {T <: iTgbm}

  push!(tsv, [e(tree), fdt(tree)])
  push!(bbλ, lλ(tree))

  if isdefined(tree, :d1)
    makebbv!(tree.d1, bbλ, tsv)
    makebbv!(tree.d2, bbλ, tsv)
  end

  return nothing
end




"""
    makebbv!(tree::iTgbmbd, 
             bbλ ::Array{Array{Float64,1},1}, 
             bbμ ::Array{Array{Float64,1},1}, 
             tsv ::Array{Array{Float64,1},1})

Make `bbv` vector with allocated `bb` (brownian bridges) and 
with `tsv` vector of branches times `ts`.
"""
function makebbv!(tree::iTgbmbd, 
                  bbλ ::Array{Array{Float64,1},1}, 
                  bbμ ::Array{Array{Float64,1},1}, 
                  tsv ::Array{Array{Float64,1},1})

  push!(tsv, [e(tree), fdt(tree)])
  push!(bbλ, lλ(tree))
  push!(bbμ, lμ(tree))

  if isdefined(tree, :d1)
    makebbv!(tree.d1, bbλ, bbμ, tsv)
    makebbv!(tree.d2, bbλ, bbμ, tsv)
  end

  return nothing
end




