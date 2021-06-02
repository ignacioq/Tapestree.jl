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

isfix(::Nothing) = false




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
pe(tree::T) where {T <: iTree} = getproperty(tree, :pe)

"""
    pe(::Nothing)

Return pendant edge.
"""
pe(::Nothing) = nothing




"""
    dt(tree::T) where {T <: iTree}

Return `δt`.
"""
dt(tree::T) where {T <: iTree} = getproperty(tree, :dt)

"""
    dt(::Nothing)

Return `δt`.
"""
dt(::Nothing) = nothing





"""
    fdt(tree::T) where {T <: iTree}

Return final `δt`.
"""
fdt(tree::T) where {T <: iTree} = getproperty(tree, :fdt)

"""
    fdt(::Nothing)

Return final `δt`.
"""
fdt(::Nothing) = nothing






"""
    treelength(tree::T) where {T <: iTree}

Return the branch length sum of `tree`.
"""
treelength(tree::T) where {T <: iTree} = 
  treelength(tree.d1) + treelength(tree.d2) + pe(tree)

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
    snan(tree::T) where {T <: iTree}

Return the number of alive nodes for `tree`.
"""
function snan(tree::T) where {T <: iTree}
    if istip(tree) && !isextinct(tree)
      return 1
    else
      return snan(tree.d1) + snan(tree.d2)
    end
end


"""
    snan(::Nothing)

Return the number of alive nodes for `tree`.
"""
snan(::Nothing) = 0





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
    lλ(tree::T) where {T <: iTgbm}

Return pendant edge.
"""
lλ(tree::T) where {T <: iTgbm} = getproperty(tree,:lλ)

"""
    pe(tree::iTree)

Return pendant edge.
"""
lλ(::Nothing) = nothing




"""
    lμ(tree::iTgbm)

Return pendant edge.
"""
lμ(tree::iTgbmbd) = getproperty(tree,:lμ)

"""
    pe(tree::iTree)

Return pendant edge.
"""
lμ(::Nothing) = nothing




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

    if istip(tree) && !isextinct(tree)
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

  if ix == ldr
    pei = pe(tree)
    if th > h > (th - pei)

      bh  = th - h
      tsi = ts(tree)
      hi = indmindif_sorted(tsi, bh)
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
        λμath(tree.d1::iTgbmbd, h, th - pe(tree), dri, ldr, ix)
      else
        λμath(tree.d2::iTgbmbd, h, th - pe(tree), dri, ldr, ix)
      end
    elseif ifx1
      λμath(tree.d1::iTgbmbd, h, th - pe(tree), dri, ldr, ix)
    else
      λμath(tree.d2::iTgbmbd, h, th - pe(tree), dri, ldr, ix)
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
    if isfix(tree.d1::T)
      ifxe(tree.d1::T)
    else
      ifxe(tree.d2::T)
    end
  end
end





"""
    fixds(tree::iTgbmbd)

Return the first fixed daughters `d1` and `d2`. Works only for 
internal branches.
"""
function fixds(tree::iTgbmbd)
  ifx1 = isfix(tree.d1::iTgbmbd)
  if ifx1 && isfix(tree.d2::iTgbmbd)
    return tree.d1::iTgbmbd, tree.d2::iTgbmbd
  elseif ifx1
    fixds(tree.d1::iTgbmbd)
  else
    fixds(tree.d2::iTgbmbd)
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

  push!(tsv, [pe(tree), fdt(tree)])
  push!(bbλ, lλ(tree))
  push!(bbμ, lμ(tree))

  makebbv!(tree.d1, bbλ, bbμ, tsv)
  makebbv!(tree.d2, bbλ, bbμ, tsv)

  return nothing
end




"""
    makebbv!(tree::Nothing, 
             bbλ ::Array{Array{Float64,1},1}, 
             bbμ ::Array{Array{Float64,1},1}, 
             tsv ::Array{Array{Float64,1},1})

Make `bbv` vector with allocated `bb` (brownian bridges) and 
with `tsv` vector of branches times `ts`.
"""
makebbv!(tree::Nothing, 
         bbλ ::Array{Array{Float64,1},1}, 
         bbμ ::Array{Array{Float64,1},1}, 
         tsv ::Array{Array{Float64,1},1}) = nothing


