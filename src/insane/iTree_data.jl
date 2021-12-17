#=

insane tree data gatherers

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#



"""
    iscrowntree(tree::T) where {T <: iTree}

Return if the tree is a crown tree.
"""
iscrowntree(tree::T) where {T <: iTree} = iszero(e(tree))




"""
    isfix(tree::T) where {T <: iTree} 

Return if is a fixed (i.e. observed) node.
"""
isfix(tree::T) where {T <: iTree} = getproperty(tree, :fx)




"""
    istip(tree::T) where {T <: iTree} = !isdefined(tree, :d1)
    istip(tree::sTfbd) = !isdefined(tree, :d1) && !isdefined(tree, :d2)

Return if is either an extant or extinct tip node.
"""
istip(tree::T) where {T <: iTree} = !isdefined(tree, :d1)
istip(tree::sTfbd) = !isdefined(tree, :d1) && !isdefined(tree, :d2)




"""
    isextinct(tree::T) where {T <: iTree}

Return if is an extinction node.
"""
isextinct(tree::T) where {T <: iTree} = getproperty(tree, :iμ)




"""
    isextinct(tree::sTpb)
    isextinct(tree::iTgbmpb)

Return if is an extinction node.
"""
isextinct(tree::sTpb) = false
isextinct(tree::iTgbmpb) = false




"""
    isalive(tree::T) where {T <: iTree}

Return if is an alive node.
"""
isalive(tree::T) where {T <: iTree} = !isextinct(tree)




"""
    isfossil(tree::T) where {T <: iTree}

Return if is a fossil tip node.
"""
isfossil(tree::T) where {T <: iTree} = getproperty(tree, :iψ)




"""
    isfossil(tree::sTpb)
    isfossil(tree::sTbd)
    isfossil(tree::iTgbmpb)
    isfossil(tree::iTgbmce)
    isfossil(tree::iTgbmct)
    isfossil(tree::iTgbmbd)

Return if is a fossil tip node : false because not allowed for those tree types.
"""
isfossil(tree::sTpb) = false
isfossil(tree::sTbd) = false
isfossil(tree::iTgbmpb) = false
isfossil(tree::iTgbmce) = false
isfossil(tree::iTgbmct) = false
isfossil(tree::iTgbmbd) = false




"""
    issampledancestor(tree::T) where {T <: iTree}

Return if is a sampled ancestor, i.e. a fossil internal node.
"""
issampledancestor(tree::T) where {T <: iTree} = isfossil(tree) && !istip(tree)




"""
    e(tree::T) where {T <: iTree}

Return edge length.
"""
e(tree::T) where {T <: iTree} = getproperty(tree, :e)





"""
    l(tree::sT_label)
    l(tree::sTf_label)

Return label.
"""
l(tree::sT_label) = getproperty(tree, :l)
l(tree::sTf_label) = getproperty(tree, :l)




"""
    tiplabels!(tree::sT_label)
    tiplabels!(tree::sTf_label)

Returns tip labels for `sT_label` and `sTf_label`.
"""
tiplabels(tree::sT_label) = _tiplabels!(tree, String[])
tiplabels(tree::sTf_label) = _tiplabels!(tree, String[])




"""
    _tiplabels!(tree::sT_label, labels::Array{String,1})

Returns tip labels for `sT_label`.
"""
function _tiplabels!(tree::sT_label, labels::Array{String,1})

  if !isdefined(tree, :d1)
    push!(labels, l(tree))
  else
    _tiplabels!(tree.d1, labels)
    _tiplabels!(tree.d2, labels)
  end
end




"""
    _tiplabels!(tree::sTf_label, labels::Array{String,1})

Returns tip labels for `sTf_label`.
"""
function _tiplabels!(tree::sTf_label, labels::Array{String,1})
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  if !defd1 && !defd2
    push!(labels, l(tree))
  else
    if defd1 _tiplabels!(tree.d1, labels) end
    if defd2 _tiplabels!(tree.d2, labels) end
  end
  return labels
end




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
treelength(tree::T) where {T <: iTree} = _treelength(tree, 0.0)




"""
    _treelength(tree::T, l::Float64) where {T <: iTree}

Return the branch length sum of `tree`, initialized at `l`.
"""
function _treelength(tree::T, l::Float64) where {T <: iTree}
  l += e(tree)
  if isdefined(tree, :d1)
    l = _treelength(tree.d1, l)::Float64
    l = _treelength(tree.d2, l)::Float64
  end
  return l
end




"""
    _treelength(tree::sTfbd, l::Float64)

Return the branch length sum of `tree`, initialized at `l`.
"""
function _treelength(tree::sTfbd, l::Float64)
  l += e(tree)
  
  if isdefined(tree, :d1) l = _treelength(tree.d1, l)::Float64 end
  if isdefined(tree, :d2) l = _treelength(tree.d2, l)::Float64 end
  
  return l
end




"""
    _ctl(tree::T, l::Float64) where {T <: iTree}
Return the branch length sum of `tree` based on `δt` and `fδt` 
for debugging purposes.
"""
function _ctl(tree::T, l::Float64) where {T <: iTgbm}

  l += Float64(lastindex(lλ(tree)) - 2)*dt(tree) + fdt(tree)

  if isdefined(tree, :d1)
    l = _ctl(tree.d1, l)
    l = _ctl(tree.d2, l)
  end

  return l
end





"""
    _ctl(tree::sTfbd, l::Float64)

Return the branch length sum of `tree` based on `δt` and `fδt` 
for debugging purposes.
"""
function _ctl(tree::sTfbd, l::Float64)

  l += Float64(lastindex(lλ(tree)) - 2)*dt(tree) + fdt(tree)

  if isdefined(tree, :d1) l = _ctl(tree.d1, l)  end
  if isdefined(tree, :d2) l = _ctl(tree.d2, l)  end

  return l
end




"""
    treeheight(tree::T) where {T <: iTree}

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
    treeheight(tree::T) where {T <: Union{sTfbd, sTf_label}}

Return the tree height of `tree`.
"""
function treeheight(tree::T) where {T <: Union{sTfbd, sTf_label}}
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  if defd1 && defd2
    th1 = treeheight(tree.d1)
    th2 = treeheight(tree.d2)
    return max(th1,th2) + e(tree)
  end

  if defd1 return treeheight(tree.d1) + e(tree) end
  if defd2 return treeheight(tree.d2) + e(tree) end

  return e(tree)
end




"""
    nnodes(tree::T) where {T <: iTree}

Return the number of descendant nodes for `tree`.
"""
nnodes(tree::T) where {T <: iTree} = _nnodes(tree, 0)




"""
    _nnodes(tree::T, n::Int64) where {T <: iTree}

Return the number of descendant nodes for `tree`, initialized at `n`.
"""
function _nnodes(tree::T, n::Int64) where {T <: iTree} 
  n += 1
  if isdefined(tree, :d1)
    n = _nnodes(tree.d1, n)
    n = _nnodes(tree.d2, n)
  end

  return n
end




"""
    _nnodes(tree::T, n::Int64) where {T <: Union{sTfbd, sTf_label}}

Return the number of descendant nodes for `tree`, initialized at `n`.
"""
function _nnodes(tree::T, n::Int64) where {T <: Union{sTfbd, sTf_label}}
  n += 1
  
  if isdefined(tree, :d1) n = _nnodes(tree.d1, n) end
  if isdefined(tree, :d2) n = _nnodes(tree.d2, n) end
  
  return n
end




"""
    nnodesinternal(tree::T) where {T <: iTree}

Return the number of internal nodes for `tree`.
"""
nnodesinternal(tree::T) where {T <: iTree} = _nnodesinternal(tree, 0)




"""
    _nnodesinternal(tree::T, n::Int64) where {T <: iTree}

Return the number of internal nodes for `tree`, initialized at `n`.
"""
function _nnodesinternal(tree::T, n::Int64) where {T <: iTree}
  if isdefined(tree, :d1)
    n += 1
    n = _nnodesinternal(tree.d1, n)
    n = _nnodesinternal(tree.d2, n)
  end

  return n
end




"""
    _nnodesinternal(tree::sTfbd, n::Int64)

Return the number of internal nodes for `tree`, initialized at `n`.
"""
function _nnodesinternal(tree::sTfbd, n::Int64)
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  if defd1 n = _nnodesinternal(tree.d1, n) end
  if defd2 n = _nnodesinternal(tree.d2, n) end
  if defd1 || defd2 n += 1 end   # considers sampled ancestors as internal nodes

  return n
end




"""
    ntips(tree::T) where {T <: iTree}

Return the number of tip nodes for `tree`.
"""
ntips(tree::T) where {T <: iTree} = _ntips(tree, 0)




"""
    _ntips(tree::T, n::Int64) where {T <: iTree}

Return the number of tip nodes for `tree`, initialized at `n`.
"""
function _ntips(tree::T, n::Int64) where {T <: iTree}

  if isdefined(tree, :d1) 
    n = _ntips(tree.d1, n)
    n = _ntips(tree.d2, n)
  else
    n += 1
  end

  return n
end



"""
    _ntips(tree::T, n::Int64) where {T <: Union{sTfbd, sTf_label}}

Return the number of tip nodes for `tree`, initialized at `n`.
"""
function _ntips(tree::T, n::Int64) where {T <: Union{sTfbd, sTf_label}}
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if defd1 n = _ntips(tree.d1, n) end
  if defd2 n = _ntips(tree.d2, n) end
  if !defd1 && !defd2 n += 1 end

  return n
end




"""
    ntipsalive(tree::T) where {T <: iTree}

Return the number of alive nodes for `tree`.
"""
ntipsalive(tree::T) where {T <: iTree} = _ntipsalive(tree, 0)



"""
    _ntipsalive(tree::T, n::Int64) where {T <: iTree}
Return the number of alive nodes for `tree`, initialized at `n`.
"""
function _ntipsalive(tree::T, n::Int64) where {T <: iTree}

  if isdefined(tree, :d1)
    n = _ntipsalive(tree.d1, n)
    n = _ntipsalive(tree.d2, n)
  elseif isalive(tree)
    n += 1
  end

  return n
end




"""
    _ntipsalive(tree::sTfbd, n::Int64)

Return the number of alive nodes for `tree`, initialized at `n`.
"""
function _ntipsalive(tree::sTfbd, n::Int64)
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  if defd1 n = _ntipsalive(tree.d1, n) end
  if defd2 n = _ntipsalive(tree.d2, n) end
  if !defd1 && !defd2 && isalive(tree) n += 1 end

  return n
end




"""
    ntipsextinct(tree::T) where {T <: iTree}

Return the number of extinct nodes for `tree`.
"""
ntipsextinct(tree::T) where {T <: iTree} = _ntipsextinct(tree, 0)




"""
    _ntipsextinct(tree::T, n::Int64) where {T <: iTree}

Return the number of extinct nodes for `tree`, initialized at `n`.
"""
function _ntipsextinct(tree::T, n::Int64) where {T <: iTree}

  if isdefined(tree, :d1) 
    n = _ntipsextinct(tree.d1, n) 
    n = _ntipsextinct(tree.d2, n) 
  else
    if isextinct(tree) 
      n += 1
    end
  end

  return n
end




"""
    _ntipsextinct(tree::sTfbd, n::Int64)

Return the number of extinct nodes for `tree`, initialized at `n`.
"""
function _ntipsextinct(tree::sTfbd, n::Int64)
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  if defd1 n = _ntipsextinct(tree.d1, n) end
  if defd2 n = _ntipsextinct(tree.d2, n) end
  if isextinct(tree) n += 1 end

  return n
end




"""
    ntipsextinctF(tree::T) where {T <: iTree}

Return the number of extinct nodes for `tree` as Float64.
"""
ntipsextinctF(tree::T) where {T <: iTree} = _ntipsextinctF(tree, 0.0)




"""
    _ntipsextinctF(tree::T, n::Float64) where {T <: iTree}

Return the number of extinct nodes for `tree` as Float64, initialized at `n`.
"""
function _ntipsextinctF(tree::T, n::Float64) where {T <: iTree}
  
  if isdefined(tree, :d1) 
    n = _ntipsextinctF(tree.d1, n) 
    n = _ntipsextinctF(tree.d2, n) 
  else
    if isextinct(tree) 
      n += 1.0
    end
  end

  return n
end




"""
    _ntipsextinctF(tree::sTfbd, n::Float64)

Return the number of extinct nodes for `tree` as Float64, initialized at `n`.
"""
function _ntipsextinctF(tree::sTfbd, n::Float64)
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  if defd1 n = _ntipsextinctF(tree.d1, n) end
  if defd2 n = _ntipsextinctF(tree.d2, n) end
  if isextinct(tree) n += 1.0 end

  return n
end




"""
    nfossils(tree::T) where {T <: iTree}

Return the number of fossil nodes for `tree`.
"""
nfossils(tree::T) where {T <: iTree} = _nfossils(tree, 0)




"""
    _nfossils(tree::T, n::Int64) where {T <: iTree}

Return the number of fossil nodes for `tree`, initialized at `n`.
"""
function _nfossils(tree::T, n::Int64) where {T <: iTree}
  if isfossil(tree)
    n += 1
  end
  
  if isdefined(tree, :d1) n = _nfossils(tree.d1, n) end
  if isdefined(tree, :d2) n = _nfossils(tree.d2, n) end

  return n
end




"""
    nsampledancestors(tree::T) where {T <: iTree}

Return the number of fossil nodes for `tree`.
"""
nsampledancestors(tree::T) where {T <: iTree} = _nsampledancestors(tree, 0)




"""
    _nsampledancestors(tree::T, n::Int64) where {T <: iTree}

Return the number of fossil nodes for `tree`, initialized at `n`.
"""
function _nsampledancestors(tree::T, n::Int64) where {T <: iTree}
  if issampledancestor(tree)
    n += 1
  end
  
  if isdefined(tree, :d1) n = _nsampledancestors(tree.d1, n) end
  if isdefined(tree, :d2) n = _nsampledancestors(tree.d2, n) end

  return n
end




"""
    treelength_ns(tree::T, 
                  l   ::Float64, 
                  n   ::Float64) where {T <: iTree}

Return the tree length and speciation events events.
"""
function treelength_ns(tree::T, 
                       l   ::Float64, 
                       n   ::Float64) where {T <: iTree}

  l += e(tree)
  if isdefined(tree, :d1)
    n += 1.0
    l, n = treelength_ns(tree.d1, l, n)
    l, n = treelength_ns(tree.d2, l, n)
  end

  return l, n
end




"""
    survives(tree::T) where {T <: iTree}

Return if the tree survives until present.
"""
function survives(tree::T) where {T <: iTree}

  if istip(tree) && isalive(tree) && !isfossil(tree)
    return true
  else
    return (isdefined(tree, :d1) && survives(tree.d1)) ||
           (isdefined(tree, :d2) && survives(tree.d2))
  end
end




"""
    treelength_ne(tree::T, 
                  l   ::Float64, 
                  n   ::Float64) where {T <: iTree}

Return the tree length and extinction events.
"""
function treelength_ne(tree::T, 
                       l   ::Float64, 
                       n   ::Float64) where {T <: iTree}

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
    treelength_ne(tree::sTfbd, 
                  l   ::Float64, 
                  n   ::Float64)

Return the tree length and extinction events.
"""
function treelength_ne(tree::sTfbd, 
                       l   ::Float64, 
                       n   ::Float64)

  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  l += e(tree)
  
  if defd1  l, n = treelength_ne(tree.d1, l, n) end
  if defd2  l, n = treelength_ne(tree.d2, l, n) end

  if !defd1 && !defd2 && isextinct(tree)
    n += 1.0
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
    fixtip(tree::T) where {T <: iTree}

Return the first fixed tip.
"""
function fixtip(tree::T) where {T <: iTree}
  if istip(tree)
    return tree
  elseif isfix(tree.d1::T)
    fixtip(tree.d1::T)
  else
    fixtip(tree.d2::T)
  end
end




"""
    fixds(tree::T) where {T <: iTgbm}

Return the first fixed daughters `d1` and `d2`. Works only for 
internal branches.
"""
function fixds(tree::T) where {T <: iTree}
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
    fixds(tree::sTfbd)

Return the first fixed daughters `d1` and `d2`. Works only for 
internal branches.
"""
function fixds(tree::sTfbd)
  ifx1 = isdefined(tree, :d1) && isfix(tree.d1::sTfbd)
  if ifx1 && isdefined(tree, :d2) && isfix(tree.d2::sTfbd)
    return tree.d1::sTfbd, tree.d2::sTfbd
  elseif ifx1
    fixds(tree.d1::sTfbd)
  else
    fixds(tree.d2::sTfbd)
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
    fixdstree(tree::sTfbd)

Returns the first tree with every daughter fixed.
"""
function fixdstree(tree::sTfbd)

  ifx1 = isdefined(tree, :d1) && isfix(tree.d1::sTfbd)
  ifx2 = isdefined(tree, :d2) && isfix(tree.d2::sTfbd)

  if (ifx1 && ifx2) || isfossil(tree)
    # every defined daughter is fixed, or sampled ancestor
    return tree
  elseif ifx1
    tree = fixdstree(tree.d1::sTfbd)
  else
    tree = fixdstree(tree.d2::sTfbd)
  end
end




"""
    drtree(tree::T, dri::BitArray{1}, ldr::Int64, ix::Int64)

Returns the tree given by `dr`, circumventing unfix branches.
"""
function drtree(tree::T, 
                dri ::BitArray{1}, 
                ldr ::Int64, 
                ix  ::Int64) where T <: iTgbm
  if ldr === ix
    return tree
  elseif ldr > ix
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      ix += 1
      if dri[ix]
        drtree(tree.d1::T, dri, ldr, ix)
      else
        drtree(tree.d2::T, dri, ldr, ix)
      end
    elseif ifx1
      drtree(tree.d1::T, dri, ldr, ix)
    else
      drtree(tree.d2::T, dri, ldr, ix)
    end
  end
end




"""
    drtree(tree::sTfbd, dri::BitArray{1}, ldr::Int64, ix::Int64)

Returns the tree given by `dr`, circumventing unfix branches.
"""
function drtree(tree::sTfbd, 
                dri ::BitArray{1}, 
                ldr ::Int64, 
                ix  ::Int64)
  if ldr === ix
    return tree
  elseif ldr > ix
    ifx1 = isdefined(tree, :d1) && isfix(tree.d1)
    if ifx1 && isdefined(tree, :d2) && isfix(tree.d2)
      ix += 1
      if dri[ix]
        drtree(tree.d1::sTfbd, dri, ldr, ix)
      else
        drtree(tree.d2::sTfbd, dri, ldr, ix)
      end
    elseif ifx1
      drtree(tree.d1::sTfbd, dri, ldr, ix)
    else
      drtree(tree.d2::sTfbd, dri, ldr, ix)
    end
  end
end




"""
    drtree(treec::T, treep::T, dri::BitArray{1}, ldr::Int64, ix::Int64)

Returns the trees given by `dr`, circumventing unfix branches.
"""
function drtree(treec::T, 
                treep::T, 
                dri  ::BitArray{1}, 
                ldr  ::Int64, 
                ix   ::Int64) where T <: iTree

  if ldr === ix
    return treec, treep
  elseif ldr > ix
    ifx1 = isfix(treec.d1)
    if ifx1 && isfix(treec.d2)
      ix += 1
      if dri[ix]
        drtree(treec.d1::T, treep.d1::T, dri, ldr, ix)
      else
        drtree(treec.d2::T, treep.d2::T, dri, ldr, ix)
      end
    elseif ifx1
      drtree(treec.d1::T, treep.d1::T, dri, ldr, ix)
    else
      drtree(treec.d2::T, treep.d2::T, dri, ldr, ix)
    end
  end
end




"""
    drtree(treec::sTfbd, treep::sTfbd, dri::BitArray{1}, ldr::Int64, ix::Int64)

Returns the trees given by `dr`, circumventing unfix branches.
"""
function drtree(treec::sTfbd, 
                treep::sTfbd, 
                dri  ::BitArray{1}, 
                ldr  ::Int64, 
                ix   ::Int64)

  if ldr === ix
    return treec, treep
  elseif ldr > ix
    ifx1 = isdefined(tree, :d1) && isfix(tree.d1)
    if ifx1 && isdefined(tree, :d2) && isfix(tree.d2)
      ix += 1
      if dri[ix]
        drtree(treec.d1::sTfbd, treep.d1::sTfbd, dri, ldr, ix)
      else
        drtree(treec.d2::sTfbd, treep.d2::sTfbd, dri, ldr, ix)
      end
    elseif ifx1
      drtree(treec.d1::sTfbd, treep.d1::sTfbd, dri, ldr, ix)
    else
      drtree(treec.d2::sTfbd, treep.d2::sTfbd, dri, ldr, ix)
    end
  end
end




"""
    survivaldr(tree::sTfbd)

Returns the directory of the first node with 2 surviving daughters.
"""
survivaldr(tree::sTfbd) = _survivaldr(tree, BitArray{1}())




"""
    _survivaldr(tree::sTfbd, bit::BitArray{1})

Returns the directory of the first node with 2 surviving daughters, initialized
at `bit`.
"""
function _survivaldr(tree::sTfbd, bit::BitArray{1})

  if isdefined(tree, :d1) && survives(tree.d1)
    if isdefined(tree, :d2) && survives(tree.d2)
      # both daughters survive
      return bit
    else
      # sampled ancestor or only d1 surviving
      push!(bit,true)
      _survivaldr(tree.d1, bit)
    end
  elseif isdefined(tree, :d2) && survives(tree.d2)
    # sampled ancestor or only d2 surviving
    push!(bit,false)
    _survivaldr(tree.d2, bit)
  else
    # the input tree had only one extant species, return this extant tip
    return bit
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
  push!(bbλ, copy(lλ(tree)))

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
  push!(bbλ, copy(lλ(tree)))
  push!(bbμ, copy(lμ(tree)))

  if isdefined(tree, :d1)
    makebbv!(tree.d1, bbλ, bbμ, tsv)
    makebbv!(tree.d2, bbλ, bbμ, tsv)
  end

  return nothing
end




"""
    eventimes(tree::T) where {T <: iTree}

Return speciation and extinction event times.
"""
function eventimes(tree::T) where {T <: iTree}
  se = Float64[]
  ee = Float64[]

  _eventimes!(tree, 0.0, se, ee)

  return se, ee
end





"""
    _eventimes!(tree::T, 
                t   ::Float64, 
                se  ::Array{Float64,1}, 
                ee  ::Array{Float64,1}) where {T <: iTree}
Recursive structure that returns speciation and extinction event times.
"""
function _eventimes!(tree::T, 
                     t   ::Float64, 
                     se  ::Array{Float64,1}, 
                     ee  ::Array{Float64,1}) where {T <: iTree}

  et = e(tree)
  if isextinct(tree)
    push!(ee, t + et)
  elseif isdefined(tree, :d1)
    push!(se, t + et)

    _eventimes!(tree.d1, t + et, se, ee)
    _eventimes!(tree.d2, t + et, se, ee)
  end

  return nothing
end




"""
    _eventimes!(tree::sTfbd, 
                t   ::Float64, 
                se  ::Array{Float64,1}, 
                ee  ::Array{Float64,1})

Recursive structure that returns speciation and extinction event times.
"""
function _eventimes!(tree::sTfbd, 
                     t   ::Float64, 
                     se  ::Array{Float64,1}, 
                     ee  ::Array{Float64,1})

  et = e(tree)
  if isextinct(tree)
    push!(ee, t + et)
  else
    defd1 = isdefined(tree, :d1)
    defd2 = isdefined(tree, :d2)
    
    if defd1 && defd2 push!(se, t + et) end
    if defd1  _eventimes!(tree.d1, t + et, se, ee)  end
    if defd2  _eventimes!(tree.d2, t + et, se, ee)  end
  end

  return nothing
end




"""
    ltt(tree::T) where {T <: iTree}

Returns number of species through time.
"""
@inline function ltt(tree::T) where {T <: iTree}

  # speciation and extinction events
  se, ee = eventimes(tree)

  # which ones are extinctions when appended
  ii = lastindex(se)

  append!(se, ee)
  l = lastindex(se)

  sp = sortperm(se)
  n  = ones(Int64, l+1)

  @inbounds begin
    @simd for i in Base.OneTo(l)
      if sp[i] > ii
        n[i+1] = n[i] - 1
      else
        n[i+1] = n[i] + 1
      end
    end
  end

  push!(se, 0.0)
  sort!(se)

  return Ltt(n, se)
end




"""
    ltt2(tree::T) where {T <: iTree}

Returns number of species through time.
"""
function ltt2(tree::T) where {T <: iTree}
  # speciation and extinction events
  se, ee = eventimes(tree)
  # start with 1 lineage
  push!(se, 0.0)
  lse = lastindex(se)
  lee = lastindex(ee)

  events = append!(se, ee)
  jumps_order = sortperm(events)
  
  jumps = ones(Int64,lse)
  append!(jumps,  fill(-1,lee))
  jumps = jumps[jumps_order]

  return Ltt(cumsum!(jumps,jumps), events[jumps_order])
end




"""
    treeapply(tree::T, FUN::Function) where {T <: iTree}

Returns a recursive vector structure with requested data for all tree nodes.
"""
function treeapply(tree::T, FUN::Function) where {T <: iTree}
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  if defd1 && defd2
    return [FUN(tree),treeapply(tree.d1,FUN),treeapply(tree.d2,FUN)]
  end
  
  if defd1
    return [FUN(tree),treeapply(tree.d1,FUN)]
  end
  
  if defd2
    return [FUN(tree),treeapply(tree.d2,FUN)]
  end

  return FUN(tree)
end




"""
    nlin_t(tree::T, t::Float64, tc::Float64) where {T <: iTree}

Number of lineages at time t
"""
function nlin_t(tree::T, t::Float64, tc::Float64) where {T <: iTree}

  el = e(tree)
  if tc < t
    if t - 1e-12 <= tc + el
      return 1
    elseif isdefined(tree, :d1)
      return nlin_t(tree.d1, t, tc + el) +
             nlin_t(tree.d2, t, tc + el)
    else
      return 0
    end
  end
end




"""
    _λat!(tree::T, 
          c   ::Float64,
          λs  ::Vector{Float64},
          t   ::Float64) where {T <: iTgbm}

Return speciation rates, `λs`, at time `c` for `tree`.
"""
function _λat!(tree::T, 
               c   ::Float64,
               λs  ::Vector{Float64},
               t   ::Float64) where {T <: iTgbm}

  et = e(tree)

  if (t + et) >= c 
    if !isfix(tree)

      lλv = lλ(tree)
      δt  = dt(tree)
      fδt = fdt(tree)

      # find final lλ
      if isapprox(c - t, et)
        Ix = lastindex(lλv) - 1
        ix = Float64(Ix) - 1.0
      else
        ix  = fld(c - t, δt)
        Ix  = Int64(ix) + 1
      end
      tii = ix*δt
      tff = tii + δt
      if tff > et
        tff = tii + fδt
      end
      eλ = linpred(c - t, tii, tff, lλv[Ix], lλv[Ix+1])

      push!(λs, eλ)

      return nothing
    end
  elseif isdefined(tree, :d1)
      _λat!(tree.d1, c, λs, t + et)
      _λat!(tree.d2, c, λs, t + et)
  end

  return nothing
end




