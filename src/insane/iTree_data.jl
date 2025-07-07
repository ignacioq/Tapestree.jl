#=

insane tree data gatherers

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    def1(tree::T) where {T <: iTree}

Wrapper for if defined daughter 1.
"""
def1(tree::T) where {T <: iTree} = isdefined(tree, :d1)




"""
    def2(tree::T) where {T <: iTree}

Wrapper for if defined daughter 2.
"""
def2(tree::T) where {T <: iTree} = isdefined(tree, :d2)




"""
    xi(tree::T) where {T <: Tx}

Return initial trait value.
"""
xi(tree::T) where {T <: Tpe} = getproperty(tree, :xi)
xi(tree::T) where {T <: Tx} = getproperty(tree, :xi)



"""
    xf(tree::T) where {T <: Tx}

Return final trait value.
"""
xf(tree::T) where {T <: Tx} = getproperty(tree, :xf)




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
isfix(tree::Tlabel) = true
isfix(tree::sTxs)   = true



"""
    sh(tree::T) where {T <: iTree}

Return `true` if punkeek shift is in `d1`, `false` if in `d2`
"""
sh(tree::T) where {T <: Tpe} = getproperty(tree, :sh)



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
isextinct(tree::Tlabel) = false
isextinct(tree::sTpb)   = false
isextinct(tree::cTpb)   = false
isextinct(tree::iTpb)   = false
isextinct(tree::sTxs)   = false




"""
    isalive(tree::T) where {T <: iTree}

Return if is an alive node.
"""
isalive(tree::T) where {T <: iTree} = !isextinct(tree) && !isfossil(tree)




"""
    anyalive(tree::T) where {T <: iTree}

Return if at least one tip is alive.
"""
function anyalive(tree::T) where {T <: iTree}
  return _anyalive(tree, false)
end




"""
    _anyalive(tree::T, a::Bool) where {T <: iTree}

Return if at least one tip is alive.
"""
function _anyalive(tree::T, a::Bool) where {T <: iTree}

  if !a
    if def1(tree)
      a = _anyalive(tree.d1, a)
      if def2(tree)
        a = _anyalive(tree.d2, a)
      end
    else
      a = isalive(tree)
    end
  end

  return a
end



"""
    isfossil(tree::T) where {T <: iTree}

Return if is a fossil tip node.
"""
isfossil(tree::T) where {T <: iTree} = getproperty(tree, :iψ)
isfossil(tree::sT_label) = false
isfossil(tree::sTpb)     = false
isfossil(tree::sTbd)     = false
isfossil(tree::iTpb)     = false
isfossil(tree::cTpb)     = false
isfossil(tree::iTce)     = false
isfossil(tree::iTct)     = false
isfossil(tree::iTbd)     = false
isfossil(tree::sTpe)     = false
isfossil(tree::sTxs)     = false




"""
    isinternalfossiltip(tree::T, f::Bool) where {T <: iTree}

Return if the fix has a sampled ancestor, i.e. a fossil internal node.
"""
isinternalfossil(tree::T) where {T <: iTree} =
  _isinternalfossil(tree, false)




"""
    _isinternalfossil(tree::T, f::Bool) where {T <: iTree}

Return if the fix has a sampled ancestor, i.e. a fossil internal node.
"""
function _isinternalfossil(tree::T, f::Bool) where {T <: iTree}

  if isfix(tree) && isfossil(tree) && istip(tree)
    return true
  elseif !istip(tree)
    if isfix(tree.d1) && def1(tree)
      f = _isinternalfossil(tree.d1, f)
    elseif def2(tree)
      f = _isinternalfossil(tree.d2, f)
    end
  end

  return f
end




"""
    e(tree::T) where {T <: iTree}

Return edge length.
"""
e(tree::T) where {T <: iTree} = getproperty(tree, :e)




"""
    label(tree::T) where {T <: iTree}
    label(tree::sT_label)
    label(tree::sTf_label)

Return label.
"""
label(tree::T) where {T <: iTree} = nothing
label(tree::sT_label)  = getproperty(tree, :l)
label(tree::sTf_label) = getproperty(tree, :l)




"""
    tiplabels(tree::T) where {T <: Tlabel}

Returns tip labels for `sT_label` and `sTf_label`.
"""
tiplabels(tree::T) where {T <: Tlabel} = _tiplabels!(tree, String[])



"""
    _tiplabels!(tree::sTf_label, labels::Array{String,1})

Returns tip labels for `sTf_label`.
"""
function _tiplabels!(tree::T, labels::Array{String,1}) where {T <: sT}

  if istip(tree)
    push!(labels, label(tree))
  else
    _tiplabels!(tree.d1, labels)
    if def2(tree)
      _tiplabels!(tree.d2, labels)
    end
  end
  return labels
end




"""
    alivetiplabels(tree::T) where {T <: Tlabel}

Returns tip labels for `sT_label` and `sTf_label`.
"""
alivetiplabels(tree::T) where {T <: Tlabel} = _alivetiplabels!(tree, String[])




"""
    _alivetiplabels!(tree::sTf_label, labels::Array{String,1})

Returns tip labels for `sTf_label`.
"""
function _alivetiplabels!(tree::T, labels::Array{String,1}) where {T <: sT}

  if istip(tree)
    if !isfossil(tree)
      push!(labels, label(tree))
    end
  else
    _alivetiplabels!(tree.d1, labels)
    if def2(tree)
      _alivetiplabels!(tree.d2, labels)
    end
  end
  return labels
end



"""
    fossillabels(tree::T) where {T <: Tlabel}

Returns tip labels for `sT_label` and `sTf_label`.
"""
fossillabels(tree::T) where {T <: Tlabel} = _fossillabels!(tree, String[])



"""
    _fossillabels!(tree::sTf_label, labels::Array{String,1})

Returns tip labels for `sTf_label`.
"""
function _fossillabels!(tree::T, labels::Array{String,1}) where {T <: sT}

  if isfossil(tree) 
    push!(labels, label(tree))
  end

  if def1(tree)
    _fossillabels!(tree.d1, labels)
    if def2(tree)
      _fossillabels!(tree.d2, labels)
    end
  end


  return labels
end




"""
    fossiltiplabels(tree::T) where {T <: Tlabel}

Returns tip labels for `sT_label` and `sTf_label`.
"""
fossiltiplabels(tree::T) where {T <: Tlabel} = _fossiltiplabels!(tree, String[])



"""
    _fossiltiplabels!(tree::sTf_label, labels::Array{String,1})

Returns tip labels for `sTf_label`.
"""
function _fossiltiplabels!(tree::T, labels::Array{String,1}) where {T <: sT}

  if istip(tree)
    if isfossil(tree)
      push!(labels, label(tree))
    end
  else
    _fossiltiplabels!(tree.d1, labels)
    if def2(tree)
      _fossiltiplabels!(tree.d2, labels)
    end
  end
  return labels
end




"""
    subclade(tree::iTree, ix::Int64)

Return the minimum stem subclade according to recursive position `ix`.
"""
subclade(tree::T, ix::Int64) where {T <: iTree} = 
  _subclade(tree, 0, ix, T())[2]

"""
    _subclade(tree::iTree, i::Int64,  ix::Int64)

Return the minimum stem subclade according to recursive position `ix`.
"""
function _subclade(tree::T, i::Int64,  ix::Int64, t0::T) where {T <: iTree}

  i += 1

  if i <= ix
    if i === ix
      return i + 1, T(tree)
    else
      if def1(tree)
        i, t0 = _subclade(tree.d1, i, ix, t0)
        if def2(tree)
          i, t0 = _subclade(tree.d2, i, ix, t0)
        end
      end
    end
  end

  return i, t0
end



"""
    subclade_fx(tree::iTree, ix::Int64)

Return the minimum stem subclade according to recursive position `ix`.
"""
subclade_fx(tree::T, ix::Int64) where {T <: iTree} = 
  _subclade_fx(tree, 0, ix, T(), true)[2]

"""
    _subclade_fx(tree::T, i::Int64, ix::Int64, t0::T, sis::Bool) where {T <: iTree}

Return the minimum stem subclade according to recursive position `ix`.
"""
function _subclade_fx(tree::T, i::Int64, ix::Int64, t0::T, sis::Bool) where {T <: iTree}

  if isfix(tree) && sis
    i += 1
  end

  if i <= ix
    if i === ix
      return i + 1, T(tree)
    else
      if def1(tree)
        if def2(tree)
          i, t0 = _subclade_fx(tree.d1, i, ix, t0, isfix(tree.d2))
          i, t0 = _subclade_fx(tree.d2, i, ix, t0, isfix(tree.d1))
        else
          i, t0 = _subclade_fx(tree.d1, i, ix, t0, true)
        end
      end
    end
  end

  return i, t0
end



"""
    subclade(trees::Vector{T}, 
                  ltree::sT_label, 
                  tips ::Vector{String},
                  stem ::Bool) where {T <: iTree}

Return the minimum subclade that includes tip labels in `tips`.
"""
function subclade(trees::Vector{T}, 
                  ltree::sT_label, 
                  tips ::Vector{String},
                  stem ::Bool) where {T <: iTree}

  vT = T[]
  for t in trees
    push!(vT, subclade(t, ltree, tips, stem))
  end

  return vT
end




"""
    subclade(tree::sT_label, tips::Vector{String})

Return the minimum subclade that includes tip labels in `tips`.
"""
function subclade(tree::sT_label, tips::Vector{String}, stem::Bool)

  ls, n2v = labels(tree)
  if stem
    return _subclade_stem(tree,  1:lastindex(ls), tips, ls, n2v)
  else
    return _subclade_crown(tree, 1:lastindex(ls), tips, ls, n2v)
  end
end




"""
    subclade(tree::iTree, 
             ltree::sT_label, 
             tips ::Vector{String}, 
             stem ::Bool)

Return the minimum subclade that includes tip labels in `tips`.
"""
function subclade(tree::iTree, 
                  ltree::sT_label, 
                  tips ::Vector{String}, 
                  stem ::Bool)

  ls, n2v = labels(tree, ltree)
  if stem
    return _subclade_stem(tree,  1:lastindex(ls), tips, ls, n2v)
  else
    return _subclade_crown(tree, 1:lastindex(ls), tips, ls, n2v)
  end
end





"""
    _subclade_stem(tree::iTree,
                   ixs ::UnitRange,
                   tips::Vector{String},
                   ls  ::Vector{String},
                   n2v ::Vector{Int64})

Return the minimum subclade that includes tip labels
in `tips` (recursive function).
"""
function _subclade_stem(tree::iTree,
                        ixs ::UnitRange,
                        tips::Vector{String},
                        ls  ::Vector{String},
                        n2v ::Vector{Int64})

  @inbounds begin

    n2i = n2v[ixs[1]]
    ix1 = ixs[n2i*2 + 1:lastindex(ixs)]
    ix2 = ixs[2:n2i*2]

    it1 = false
    od1 = true
    for i in ix1
      lsi = ls[i]
      if in(lsi, tips)
        it1 = true
      elseif !(lsi === "" || lsi === "tda")
        od1 = false
      end
      it1 && !od1 && break
    end

    it2 = false
    od2 = true
    for i in ix2
      lsi = ls[i]
      if in(lsi, tips)
        it2 = true
      elseif !(lsi === "" || lsi === "tda")
        od2 = false
      end
      it2 && !od2 && break
    end

    if it1
      if !it2
        if od1
          return tree.d1
        end 
        tree = _subclade_stem(tree.d1, ix1, tips, ls, n2v)
      end
    elseif it2
      if od2
        return tree.d2
      else
        tree = _subclade_stem(tree.d2, ix2, tips, ls, n2v)
      end
    end
  end

  return tree
end





"""
    _subclade(tree::iTree,
              ixs ::UnitRange,
              tips::Vector{String},
              ls  ::Vector{String},
              n2v ::Vector{Int64})

Return the minimum subclade that includes tip labels
in `tips` (recursive function).
"""
function _subclade_crown(tree::iTree,
                         ixs ::UnitRange,
                         tips::Vector{String},
                         ls  ::Vector{String},
                         n2v ::Vector{Int64})

  @inbounds begin

    n2i = n2v[ixs[1]]
    ix1 = ixs[n2i*2 + 1:lastindex(ixs)]
    ix2 = ixs[2:n2i*2]

    it1 = false
    for si in tips
      for i in ix1
        if si === ls[i]
          it1 = true
          break
        end
        it1 && break
      end
      it1 && break
    end

    it2 = false
    for si in tips
      for i in ix2
        if si === ls[i]
          it2 = true
          break
        end
        it2 && break
      end
      it2 && break
    end

    if it1
      if !it2
        tree = _subclade_crown(tree.d1, ix1, tips, ls, n2v)
      end
    elseif it2
      tree = _subclade_crown(tree.d2, ix2, tips, ls, n2v)
    end
  end

  return tree
end




"""
    subclades_n(tree::T, 
                n_min::Int64,
                n_max::Int64) where {T <: iTree}

Return all the subclades with less than `n_max` and more than `n_min` species.
"""
function subclades_n(tree::T, 
                     n_min::Int64,
                     n_max::Int64) where {T <: iTree}

  strees = T[]
  _subclades_n!(strees,  tree, ntips(tree), n_min, n_max)

  return strees
end




"""
    _subclades_n!(strees::Vector{T},
                  tree  ::T,
                  n_min ::Int64,
                  n_max ::Int64) where {T <: iTree}

Return the minimum subclades_n that includes tip labels
in `tips` (recursive function).
"""
function _subclades_n!(strees::Vector{T},
                       tree  ::T,
                       n     ::Int64,
                       n_min ::Int64,
                       n_max ::Int64) where {T <: iTree}
  if def1(tree)
    n1 = ntips(tree.d1)

    if n_min <= n1 <= n_max
      push!(strees, T(tree.d1))
    else
      _subclades_n!(strees, tree.d1, n1, n_min, n_max)
    end

    if def2(tree)
      n2 = n - n1
      if n_min <= n2 <= n_max
        push!(strees, T(tree.d2))
      else
        _subclades_n!(strees, tree.d2, n2, n_min, n_max)
      end
    end
  end

  return nothing
end


"""
    labels(tree::Tlabel)

Return labels and left node order.
"""
function labels(tree::Tlabel)
  ls = String[]
  _make_ls!(tree, ls)
  return ls
end



"""
    make_ls!(tree::Tlabel, ls::Vector{String})

Return labels.
"""
function _make_ls!(tree::Tlabel, ls::Vector{String})

  if def1(tree)
    _make_ls!(tree.d1, ls)
    if def2(tree)
      _make_ls!(tree.d2, ls)
    elseif label(tree) != ""
        push!(ls, label(tree))
    end
  elseif label(tree) != ""
    push!(ls, label(tree))
  end
end




"""
    labels(tree::T, ltree::sT_label) where {T <: iTree}

Return labels and left node order.
"""
function labels(tree::T, ltree::sT_label) where {T <: iTree}

  ls  = String[]
  n2v = Int64[]
  make_ls!(tree, ltree, ls, n2v)

  reverse!(ls)
  reverse!(n2v)

  return ls, n2v
end




"""
    make_ls!(tree ::T,
             ltree::sT_label,
             ls   ::Array{String,1},
             n2v  ::Array{Int64,1}) where {T <: iTree}

Return labels and left node order (recursive function) for data augmented
tree.
"""
function make_ls!(tree ::T,
                  ltree::sT_label,
                  ls   ::Array{String,1},
                  n2v  ::Array{Int64,1}) where {T <: iTree}

  if istip(tree)
    if isfix(tree)
      push!(ls, label(ltree))
    else
      push!(ls, "tda")
    end

    push!(n2v, 0)
    return 1
  end

  if isfix(tree.d1) && isfix(tree.d2)
    n1 = make_ls!(tree.d1, ltree.d1, ls, n2v)
    n2 = make_ls!(tree.d2, ltree.d2, ls, n2v)
    push!(ls, label(ltree))

  else
    n1 = make_ls!(tree.d1, ltree, ls, n2v)
    n2 = make_ls!(tree.d2, ltree, ls, n2v)

    push!(ls, "")
  end

  push!(n2v, n2)

  return n1 + n2
end




"""
    fixedpos(tree::T)

Returns vector of the position of fixed tips across all tips 
of a data augmented tree.
"""
function fixedpos(tree::T) where {T <: iTree}

  fp = Int64[]
  _fixedpos!(tree, 0, fp)

  return fp
end



"""
    _fixedpos!(tree::T, i::Int64, fp::Vector{Int64}) where {T <: iTree}

Returns vector of the position of fixed tips across all tips 
of a data augmented tree.
"""

function _fixedpos!(tree::T, i::Int64, fp::Vector{Int64}) where {T <: iTree}

  if istip(tree) 
    i += 1
    if isfix(tree)
      push!(fp, i)
    end
  else
    i = _fixedpos!(tree.d1, i, fp)
    if def2(tree)
      i = _fixedpos!(tree.d2, i, fp)
    end
  end

  return i
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
  if def1(tree)
    l = _treelength(tree.d1, l)::Float64
    l = _treelength(tree.d2, l)::Float64
  end
  return l
end




"""
    _treelength(tree::T, l::Float64) where {T <: iTf}

Return the branch length sum of `tree`, initialized at `l`.
"""
function _treelength(tree::T, l::Float64) where {T <: iTf}
  l += e(tree)

  if def1(tree)
    l = _treelength(tree.d1, l)::Float64
    if def2(tree)
      l = _treelength(tree.d2, l)::Float64
    end
  end

  return l
end




"""
    treelength(tree::T, ets::Vector{Float64})  where {T <: iTf}

Return the branch length sum of `tree` at different epochs, initialized at `l`.
"""
function treelength(tree::T, ets::Vector{Float64}) where {T <: iTf}
  nep = lastindex(ets) + 1
  ls  = zeros(nep)
  _treelength!(tree, treeheight(tree), ls, ets, 1, nep)

  return ls
end




"""
    _treelength(tree::T,
                t   ::Float64,
                ls  ::Vector{Float64},
                ets ::Vector{Float64},
                ix  ::Int64,
                nep ::Int64) where {T <: iTf}

Return the branch length sum of `tree` at different epochs recursively.
"""
function _treelength!(tree::T,
                      t   ::Float64,
                      ls  ::Vector{Float64},
                      ets ::Vector{Float64},
                      ix  ::Int64,
                      nep ::Int64) where {T <: iTf}
  @inbounds begin

    ei  = e(tree)

    # if epoch change
    while ix < nep && t - ei < ets[ix]
      li      = t - ets[ix]
      ls[ix] += li
      ei     -= li
      t       = ets[ix]
      ix     += 1
    end

    ls[ix] += ei
    t      -= ei

    if def1(tree)
      _treelength!(tree.d1, t, ls, ets, ix, nep)
      if def2(tree)
        _treelength!(tree.d2, t, ls, ets, ix, nep)
      end
    end
  end

  return nothing
end




"""
    irange(tree::T, f::Function) where {T <: iTree}

Return the extrema of the output of function `f` on `tree`.
"""
function irange(tree::T, f::Function) where {T <: iTree}

  mn, mx = extrema(f(tree))

  if def1(tree)
    mnd, mxd = irange(tree.d1, f)
    mn = min(mn, mnd)
    mx = max(mx, mxd)

    if def2(tree)
      mnd, mxd = irange(tree.d2, f)
      mn = min(mn, mnd)
      mx = max(mx, mxd)
    end
  end

  return mn, mx
end




"""
    _ctl(tree::T, f::Function, l::Float64) where {T <: iTf}

Return the branch length sum of `tree` based on `δt` and `fδt`
for debugging purposes.
"""
function _ctl(tree::T, f::Function, l::Float64) where {T <: iTree}

  l += max(0.0, Float64(lastindex(f(tree)) - 2)*dt(tree)) + fdt(tree)

  if def1(tree)
    l = _ctl(tree.d1, f, l)
    if def2(tree)
      l = _ctl(tree.d2, f, l)
    end
  end

  return l
end




"""
    treeheight(tree::T) where {T <: iTree}

Return the tree height of `tree`.
"""
function treeheight(tree::T) where {T <: iTree}
  if def1(tree)
    th1 = treeheight(tree.d1)
    th2 = treeheight(tree.d2)
    return (th1 > th2 ? th1 : th2) + e(tree)
  end
  return e(tree)
end




"""
    treeheight(tree::T) where {T <: iTf}

Return the tree height of `tree`.
"""
function treeheight(tree::T) where {T <: iTf}

  if def2(tree)
    th1 = treeheight(tree.d1)
    th2 = treeheight(tree.d2)
    return max(th1,th2) + e(tree)
  elseif def1(tree)
    return treeheight(tree.d1) + e(tree)
  end

  return e(tree)
end





"""
    treeheight(tree::T, nd::Int64) where {T <: iTree}

Return the tree height of `tree`.
"""
function treeheight(tree::T, nd::Int64) where {T <: iTree}
  if def1(tree)
    th1 = treeheight(tree.d1)
    th2 = treeheight(tree.d2)
    return round(max(th1,th2) + e(tree), digits = nd)
  end
  return e(tree)
end




"""
    treeheight(tree::T, nd::Int64) where {T <: iTf}

Return the tree height of `tree`.
"""
function treeheight(tree::T, nd::Int64) where {T <: iTf}

  if def2(tree)
    th1 = treeheight(tree.d1)
    th2 = treeheight(tree.d2)
    return round(max(th1,th2) + e(tree), digits = nd)
  elseif def1(tree)
    return round(treeheight(tree.d1) + e(tree), digits = nd)
  end

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
  if def1(tree)
    n = _nnodes(tree.d1, n)
    n = _nnodes(tree.d2, n)
  end

  return n
end




"""
    _nnodes(tree::T, n::Int64) where {T <: iTf}

Return the number of descendant nodes for `tree`, initialized at `n`.
"""
function _nnodes(tree::T, n::Int64) where {T <: iTf}
  n += 1

  if def1(tree)
    n = _nnodes(tree.d1, n)
    if def2(tree)
      n = _nnodes(tree.d2, n)
    end
  end

  return n
end





"""
    nedgesF(treev::Vector{T}) where {T <: iTree}

Return the number of descendant nodes for `tree`.
"""
nedgesF(tree::T) where {T <: iTree} =  _nedgesF(tree, 0.0)



"""
    nedgesF(treev::Vector{T}) where {T <: iTree}

Return the number of descendant nodes for `tree`.
"""
function nedgesF(treev::Vector{T}) where {T <: iTree} 

  n = 0.0
  for t in treev
    n = _nedgesF(t, n)
  end
  return n
end



"""
    _nedgesF(tree::T, n::Float64) where {T <: iTree}

Return the number of descendant nodes for `tree`, initialized at `n`.
"""
function _nedgesF(tree::T, n::Float64) where {T <: iTree}
  n += 1.0

  if def1(tree)
    n = _nedgesF(tree.d1, n)
    if def2(tree)
      n = _nedgesF(tree.d2, n)
    end
  end

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
  if def1(tree)
    n += 1
    n = _nnodesinternal(tree.d1, n)
    n = _nnodesinternal(tree.d2, n)
  end

  return n
end




"""
    _nnodesinternal(tree::T, n::Int64) where {T <: iTf}

Return the number of internal nodes for `tree`, initialized at `n`.
"""
function _nnodesinternal(tree::T, n::Int64) where {T <: iTf}

  if def1(tree)
    n += 1
    n = _nnodesinternal(tree.d1, n)
    if def2(tree)
      n = _nnodesinternal(tree.d2, n)
    end
  end

  return n
end



"""
    nnodesbifurcation(idf::Vector{iBffs})

Return the number of internal nodes for `tree`, initialized at `n`.
"""
function nnodesbifurcation(idf::Vector{iBffs})

  n = 0.0
  for bi in idf
    if d2(bi) > 0
      n += 1.0
    end
  end

  return n
end




"""
    nnodesbifurcation(tree::T) where {T <: iTree}
    nnodesbifurcation(tree::T) where {T <: iTf}

Return the number of bifurcation nodes for `tree`.
"""
nnodesbifurcation(tree::T) where {T <: iTree} = _nnodesinternal(tree, 0)
nnodesbifurcation(tree::T) where {T <: iTf}   = _nnodesbifurcation(tree, 0)




"""
    _nnodesbifurcation(tree::T, n::Int64) where {T <: iTf}

Return the number of internal nodes for `tree`, initialized at `n`.
"""
function _nnodesbifurcation(tree::T, n::Int64) where {T <: iTf}

  if def1(tree)
    n = _nnodesbifurcation(tree.d1, n)
    if def2(tree)
      n += 1
      n = _nnodesbifurcation(tree.d2, n)
    end
  end

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

  if def1(tree)
    n = _ntips(tree.d1, n)
    n = _ntips(tree.d2, n)
  else
    n += 1
  end

  return n
end




"""
    _ntips(tree::T, n::Int64) where {T <: iTf}

Return the number of tip nodes for `tree`, initialized at `n`.
"""
function _ntips(tree::T, n::Int64) where {T <: iTf}

  if def1(tree)
    n = _ntips(tree.d1, n)
    if def2(tree)
      n = _ntips(tree.d2, n)
    end
  else
    n += 1
  end

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

  if def1(tree)
    n = _ntipsalive(tree.d1, n)
    n = _ntipsalive(tree.d2, n)
  elseif isalive(tree)
    n += 1
  end

  return n
end




"""
    _ntipsalive(tree::T, n::Int64) where {T <: iTf}

Return the number of alive nodes for `tree`, initialized at `n`.
"""
function _ntipsalive(tree::T, n::Int64) where {T <: iTf}

  if def1(tree)
    n = _ntipsalive(tree.d1, n)
    if def2(tree)
      n = _ntipsalive(tree.d2, n)
    end
  elseif isalive(tree) && !isfossil(tree)
    n += 1
  end

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

  if def1(tree)
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
    _ntipsextinct(tree::T, n::Int64) where {T <: iTf}

Return the number of extinct nodes for `tree`, initialized at `n`.
"""
function _ntipsextinct(tree::T, n::Int64) where {T <: iTf}

  if istip(tree)
    if isextinct(tree)
      n += 1
    end
  else
    n = _ntipsextinct(tree.d1, n)
    if def2(tree)
      n = _ntipsextinct(tree.d2, n)
    end
  end

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

  if def1(tree)
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
    _ntipsextinctF(tree::T, n::Float64) where {T <: iTf}

Return the number of extinct nodes for `tree` as Float64, initialized at `n`.
"""
function _ntipsextinctF(tree::T, n::Float64) where {T <: iTf}

  if istip(tree)
    if isextinct(tree)
      n += 1.0
    end
  else
    n = _ntipsextinctF(tree.d1, n)
    if def2(tree)
      n = _ntipsextinctF(tree.d2, n)
    end
  end

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

  if def1(tree)
    n = _nfossils(tree.d1, n)
    if def2(tree)
      n = _nfossils(tree.d2, n)
    end
  end

  return n
end




"""
    ntipfossils(tree::T) where {T <: iTree}

Return the number of fossil nodes for `tree`.
"""
ntipfossils(tree::T) where {T <: iTree} = _ntipfossils(tree, 0)




"""
    _ntipfossils(tree::T, n::Int64) where {T <: iTree}

Return the number of fossil nodes for `tree`, initialized at `n`.
"""
function _ntipfossils(tree::T, n::Int64) where {T <: iTree}

  if def1(tree)
    n = _ntipfossils(tree.d1, n)
    if def2(tree)
      n = _ntipfossils(tree.d2, n)
    end
  elseif isfossil(tree)
    n += 1
  end

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
  if def1(tree)
    n += 1.0
    l, n = treelength_ns(tree.d1, l, n)
    l, n = treelength_ns(tree.d2, l, n)
  end

  return l, n
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
  if def1(tree)
    l, n = treelength_ne(tree.d1, l, n)
    l, n = treelength_ne(tree.d2, l, n)
  elseif isextinct(tree)
    n += 1.0
  end

  return l, n
end




"""
    treelength_ne(tree::T,
                  l   ::Float64,
                  n   ::Float64) where {T <: iTf}

Return the tree length and extinction events.
"""
function treelength_ne(tree::T,
                       l   ::Float64,
                       n   ::Float64) where {T <: iTf}

  l += e(tree)
  if def1(tree)
    l, n = treelength_ne(tree.d1, l, n)
    if def2(tree)
      l, n = treelength_ne(tree.d2, l, n)
    end
  elseif isextinct(tree)
    n += 1.0
  end

  return l, n
end




"""
    lλ(tree::T) where {T <: iT}

Return pendant edge.
"""
lλ(tree::T) where {T <: iT} = getproperty(tree, :lλ)
lλ(tree::T) where {T <: cT} = getproperty(tree, :lλ)




"""
    lμ(tree::iTbdU)

Return pendant edge.
"""
lμ(tree::iTbdU) = getproperty(tree,:lμ)




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
    λμ01(tree::iTbd,
         dri ::BitArray{1},
         ldr ::Int64,
         ix  ::Int64,
         λ0  ::Float64,
         μ0  ::Float64)

Return the `λ`,`μ` at the start `0` and end `1` of a fixed branch.
"""
function λμ01(tree::iTbd,
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
      ifx1 = isfix(tree.d1::iTbd)
      if ifx1 && isfix(tree.d2::iTbd)
        λ1 = lλ(tree)[end]
        μ1 = lμ(tree)[end]

        return λ0, μ0, λ1, μ1
      elseif ifx1
        λμ01(tree.d1::iTbd, dri, ldr, ix, λ0, μ0)
      else
        λμ01(tree.d2::iTbd, dri, ldr, ix, λ0, μ0)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTbd)
    if ifx1 && isfix(tree.d2::iTbd)
      ix += 1
      if dri[ix]
        λμ01(tree.d1::iTbd, dri, ldr, ix, λ0, μ0)
      else
        λμ01(tree.d2::iTbd, dri, ldr, ix, λ0, μ0)
      end
    elseif ifx1
      λμ01(tree.d1::iTbd, dri, ldr, ix, λ0, μ0)
    else
      λμ01(tree.d2::iTbd, dri, ldr, ix, λ0, μ0)
    end
  end
end





"""
    λμath(tree::iTbd,
          h    ::Float64,
          th   ::Float64,
          dri  ::BitArray{1},
          ldr  ::Int64,
          wpr  ::Int64,
          ix   ::Int64,
          px   ::Int64)

Return the `λ`,`μ` and time of `tree` at closest to `h`.
"""
function λμath(tree::iTbd,
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
      if isfix(tree.d1::iTbd)
        λμath(tree.d1::iTbd, h, th - pei, dri, ldr, ix)
      else
        λμath(tree.d2::iTbd, h, th - pei, dri, ldr, ix)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTbd)
    if ifx1 && isfix(tree.d2::iTbd)
      ix += 1
      if dri[ix]
        λμath(tree.d1::iTbd, h, th - e(tree), dri, ldr, ix)
      else
        λμath(tree.d2::iTbd, h, th - e(tree), dri, ldr, ix)
      end
    elseif ifx1
      λμath(tree.d1::iTbd, h, th - e(tree), dri, ldr, ix)
    else
      λμath(tree.d2::iTbd, h, th - e(tree), dri, ldr, ix)
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
    fixtrait(tree::T) where {T <: iTree}

Return the trait of the final fixed tip or fossil.
"""
function fixtrait(tree::T) where {T <: iTree}
  if isfossil(tree) || istip(tree)
    return xf(tree)
  elseif isfix(tree.d1::T)
    fixtrait(tree.d1::T)
  else
    fixtrait(tree.d2::T)
  end
end




"""
    anyfossil(tree::T) where {T <: iTree}

Return true first fixed tip.
"""
function anyfossil(tree::T) where {T <: iTree}

  if isfossil(tree)
    return true
  elseif def1(tree)
    anyfossil(tree.d1)
    if def2(tree.d2)
      anyfossil(tree.d1)
    end
  end

  return false
end




"""
    makebbv!(tree::T,
             bbλ ::Array{Array{Float64,1},1},
             tsv ::Array{Array{Float64,1},1}) where {T <: iT}

Make `bbv` vector with allocated `bb` (brownian bridges) and
with `tsv` vector of branches times `ts`.
"""
function makebbv!(tree::T,
                  bbλ ::Array{Array{Float64,1},1},
                  tsv ::Array{Array{Float64,1},1}) where {T <: iT}

  push!(tsv, [e(tree), fdt(tree)])
  push!(bbλ, copy(lλ(tree)))

  if def1(tree)
    makebbv!(tree.d1, bbλ, tsv)
    makebbv!(tree.d2, bbλ, tsv)
  end

  return nothing
end




"""
    makebbv!(tree::iTbd,
             bbλ ::Array{Array{Float64,1},1},
             bbμ ::Array{Array{Float64,1},1},
             tsv ::Array{Array{Float64,1},1})
Make `bbv` vector with allocated `bb` (brownian bridges) and
with `tsv` vector of branches times `ts`.
"""
function makebbv!(tree::iTbd,
                  bbλ ::Array{Array{Float64,1},1},
                  bbμ ::Array{Array{Float64,1},1},
                  tsv ::Array{Array{Float64,1},1})

  push!(tsv, [e(tree), fdt(tree)])
  push!(bbλ, copy(lλ(tree)))
  push!(bbμ, copy(lμ(tree)))

  if def1(tree)
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

  _eventimes!(tree, treeheight(tree), se, ee)

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
    push!(ee, t - et)
  elseif def1(tree)
    push!(se, t - et)

    _eventimes!(tree.d1, t - et, se, ee)
    _eventimes!(tree.d2, t - et, se, ee)
  end

  return nothing
end




"""
    _eventimes!(tree::T,
                t   ::Float64,
                se  ::Array{Float64,1},
                ee  ::Array{Float64,1}) where {T <: iTf}

Recursive structure that returns speciation and extinction event times.
"""
function _eventimes!(tree::T,
                     t   ::Float64,
                     se  ::Array{Float64,1},
                     ee  ::Array{Float64,1}) where {T <: iTf}

  et = e(tree)
  if isextinct(tree)
    push!(ee, t - et)
  elseif def1(tree)
    if def2(tree)
      push!(se, t - et)
      _eventimes!(tree.d1, t - et, se, ee)
      _eventimes!(tree.d2, t - et, se, ee)
    else
      _eventimes!(tree.d1, t - et, se, ee)
    end
  elseif isfossil(tree)
    push!(ee, t - et)
  end

  return nothing
end



"""
    _time1_λ!(tree::T,
              t   ::Float64) where {T <: iTree}

Returns time to first speciation
"""
function _time1_λ!(tree ::T,
                   tstem::Float64) where {T <: iTree}

  tstem += e(tree)

  if def1(tree)
    if def2(tree)
      return tstem
    else
      tstem = _time1_λ!(tree.d1, tstem)
    end
  end

  return tstem
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
  lse = lastindex(se)

  sp = sortperm(se, rev = true)
  n  = ones(Int64, lse+1)

  @inbounds begin
    @simd for i in Base.OneTo(lse)
      if sp[i] > ii
        n[i+1] = n[i] - 1
      else
        n[i+1] = n[i] + 1
      end
    end
  end

  sort!(se, rev = true)
  # if crown or stem
  pushfirst!(se, se[1] + _time1_λ!(tree, 0.0))

  # last no events if alive
  push!(n,  n[end])
  push!(se, 0.0)

  return Ltt(n, se)
end




"""
    ltt(tree::Vector{T}) where {T <: iTree}

Returns number of species through time for a tree vector.
"""
function ltt(trees::Vector{T}) where {T <: iTree}
  ltv = Ltt[]
  for t in trees
    push!(ltv, ltt(t))
  end

  return ltv
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
    elseif def1(tree)
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
          t   ::Float64) where {T <: iT}

Return speciation rates, `λs`, at time `c` for `tree`.
"""
function _λat!(tree::T,
               c   ::Float64,
               λs  ::Vector{Float64},
               t   ::Float64) where {T <: iT}

  et = e(tree)

  if (t + et) > c
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
  elseif def1(tree)
      _λat!(tree.d1, c, λs, t + et)
      _λat!(tree.d2, c, λs, t + et)
  end

  return nothing
end




"""
    _xatt!(tree::T,
           c   ::Float64,
           xs  ::Vector{Float64},
           t   ::Float64,
           x   ::Float64,
           s  ::Bool) where {T <: Tpe}

Return trait `x` at time `c` for `tree`.
"""
function _xatt!(tree::T,
                c   ::Float64,
                xs  ::Vector{Float64},
                t   ::Float64,
                x   ::Float64,
                s  ::Bool) where {T <: Tpe}

  et = e(tree)

  if (t + et) >= c - accerr

    xfi = linpred(c, t, t + et, xi(tree), xf(tree))

    if isfix(tree)
      x = xfi
      s = sh(tree)
    end

    push!(xs, xfi)

    return x, s
  elseif def1(tree)
    x, s = _xatt!(tree.d1, c, xs, t + et, x, s)

    if def2(tree)
      x, s = _xatt!(tree.d2, c, xs, t + et, x, s)
    end
  end

  return x, s
end




"""
    _xisatt!(tree::T,
             c   ::Float64,
             xis ::Vector{Float64},
             es  ::Vector{Float64},
             t   ::Float64,
             na  ::Int64,
             xic ::Float64) where {T <: Tpe}

Return initial traits and edge lengths for those alive at time `c` for `tree`.
"""
function _xisatt!(tree::T,
                  c   ::Float64,
                  xis ::Vector{Float64},
                  es  ::Vector{Float64},
                  t   ::Float64,
                  na  ::Int64,
                  xic ::Float64) where {T <: Tpe}

  et = e(tree)

  if (t + et) >= c - accerr
    na += 1

    if isfix(tree)
      xic = xi(tree)
    end

    push!(xis, xi(tree))
    push!(es, c - t)

    return na, xic
  elseif def1(tree)
    na, xic = _xisatt!(tree.d1, c, xis, es, t + et, na, xic)
    if def2(tree)
      na, xic = _xisatt!(tree.d2, c, xis, es, t + et, na, xic)
    end
  end

  return na, xic
end




"""
   _xisatt!(tree::T,
            c   ::Float64,
            xis ::Vector{Float64},
            es  ::Vector{Float64},
            t   ::Float64,
            na  ::Int64,
            xic ::Float64,
            xc  ::Float64) where {T <: Tpe}

Return initial traits and edge lengths for those alive at time `c` for `tree`.
"""
function _xisatt!(tree::T,
                  c   ::Float64,
                  xis ::Vector{Float64},
                  es  ::Vector{Float64},
                  t   ::Float64,
                  na  ::Int64,
                  xic ::Float64,
                  xc  ::Float64) where {T <: Tpe}

  et = e(tree)

  if (t + et) >= c - accerr
    na += 1

    if isfix(tree)
      xic = xi(tree)
      xc  = xf(tree)
    end

    push!(xis, xi(tree))
    push!(es, c - t)

    return na, xic, xc
  elseif def1(tree)
    na, xic, xc = _xisatt!(tree.d1, c, xis, es, t + et, na, xic, xc)
    if def2(tree)
      na, xic, xc = _xisatt!(tree.d2, c, xis, es, t + et, na, xic, xc)
    end
  end

  return na, xic, xc
end




"""
   _xifsatt!(tree::T,
            c   ::Float64,
            xis ::Vector{Float64},
            es  ::Vector{Float64},
            t   ::Float64,
            na  ::Int64,
            xic ::Float64,
            xc  ::Float64) where {T <: Tpe}

Return initial traits and edge lengths for those alive at time `c` for `tree`.
"""
function _xifsatt!(tree::T,
                   c   ::Float64,
                   xis ::Vector{Float64},
                   xfs ::Vector{Float64},
                   es  ::Vector{Float64},
                   t   ::Float64,
                   na  ::Int64,
                   xic ::Float64,
                   xc  ::Float64) where {T <: Tpe}

  et = e(tree)

  if (t + et) >= c - accerr
    na += 1

    if isfix(tree)
      xic = xi(tree)
      xc  = xf(tree)
    end

    push!(xis, xi(tree))
    push!(xfs, xf(tree))
    push!(es, c - t)

    return na, xic, xc
  elseif def1(tree)
    na, xic, xc = _xifsatt!(tree.d1, c, xis, xfs, es, t + et, na, xic, xc)
    if def2(tree)
      na, xic, xc = _xifsatt!(tree.d2, c, xis, xfs, es, t + et, na, xic, xc)
    end
  end

  return na, xic, xc
end




"""
    fixed_xt(tree::T)  where {T <: Tx}

Make joint proposal to match simulation with tip fixed `x` value.
"""
function fixed_xt(tree::T)  where {T <: Tx}

  if istip(tree)
    return xf(tree)
  else
    if isfix(tree.d1)
      xt = fixed_xt(tree.d1)
    else
      xt = fixed_xt(tree.d2)
    end
  end

  return xt
end




"""
    fossil_xt(tree::T)

Make joint proposal to match simulation with tip fixed `x` value.
"""
function fossil_xt(tree::sTfbdx)

  if isfossil(tree)
    return xf(tree)
  elseif isfix(tree.d1)
    xt = fossil_xt(tree.d1)
  else
    xt = fossil_xt(tree.d2)
  end

  return xt
end




"""
    fossil_xt(tree::T)

Make joint proposal to match simulation with tip fixed `x` value.
"""
function fossil_xt(tree::iTfbdx)

  if isfossil(tree)
    return xv(tree)[end]
  elseif isfix(tree.d1)
    xt = fossil_xt(tree.d1)
  else
    xt = fossil_xt(tree.d2)
  end

  return xt
end




"""
    xv(tree::T) where {T <: Tx}

Return trait vector.
"""
xv(tree::T) where {T <: Tx} = getproperty(tree, :xv)




"""
    lσ2(tree::T) where {T <: Tx}

Return trait rate vector.
"""
lσ2(tree::T) where {T <: Txs} = getproperty(tree, :lσ2)




"""
    branchlengths(tree::T)

Extract all branch lengths.
"""
function branchlengths(tree::iTree)
  bls = Float64[]
  _branchlengths!(bls, tree)
  return bls
end




"""
    _branchlengths!(bls::Vector{Float64}, tree::iTree)

Extract all branch lengths.
"""
function _branchlengths!(bls::Vector{Float64}, tree::iTree)

  push!(bls, e(tree))

  if def1(tree)
    _branchlengths!(bls, tree.d1)
    if def2(tree)
      _branchlengths!(bls, tree.d2)
    end
  end
end

 


"""
    trextract(tree::iTree, f::Function)

Perform function `f` in each recursive tree in `tree`.
"""
function trextract(tree::iTree, f::Function)
  tvs = typeof(f(tree))[]
  _trextract!(tvs, tree, f)
  return tvs
end


"""
    _trextract!(tvs::Vector{T}, tree::iTree, f::Function) where {T}

Perform function `f` recursively in each recursive tree in `tree`.
"""
function _trextract!(tvs::Vector{T}, tree::iTree, f::Function) where {T}

  push!(tvs, f(tree))

  if def1(tree)
    _trextract!(tvs, tree.d1, f)
    if def2(tree)
      _trextract!(tvs, tree.d2, f)
    end
  end
end




"""
    tiptraits(tree::T) where {T <: Tpe}

Extract tip rates from a `Tpe` tree.
"""
function tiptraits(tree::T) where {T <: Tpe}
  x = Float64[]
  _tiptraits!(x, tree)
  return x
end


"""
    _tiptraits!(x::Vector{Float64}, tree::T) where {T <: Tpe}

Extract tip rates from a `Tpe` tree.
"""
function _tiptraits!(x::Vector{Float64}, tree::T) where {T <: Tpe}

  if def1(tree)
    _tiptraits!(x, tree.d1)
    if def2(tree)
      _tiptraits!(x, tree.d2)
    end
  else
    push!(x, xf(tree))
  end
end




"""
    traits(tree::T) where {T <: Tpe}

Extract tip rates from a `Tpe` tree.
"""
function traits(tree::T) where {T <: Tpe}
  x = Float64[]
  _traits!(x, tree)
  return x
end


"""
    _traits!(x::Vector{Float64}, tree::T) where {T <: Tpe}

Extract tip rates from a `Tpe` tree.
"""
function _traits!(x::Vector{Float64}, tree::T) where {T <: Tpe}

  if def1(tree)
    _traits!(x, tree.d1)
    if def2(tree)
      _traits!(x, tree.d2)
    else
      push!(x, xf(tree))
    end
  else
    push!(x, xf(tree))
  end
end



"""
    check_cpe!(tree::sTpe)

Extract tip rates from a sTpe tree
"""
function _check_cpe(tree::sTpe)

  if def1(tree)
    if def2(tree)

      ch = if sh(tree)
        xf(tree) === xi(tree.d2)
      else
        xf(tree) === xi(tree.d1)
      end

      if !ch 
        @show tree, sh(tree), xf(tree), xi(tree.d2), xi(tree.d1)
      end

      return _check_cpe(tree.d1) && _check_cpe(tree.d2) && ch 
    else
      return _check_cpe(tree.d1)
    end
  end

  return true
end




