#=

insane tree manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    scale_rate!(tree::T, f::Function, s::Float64) where {T <: cT}

Add `s` to vector retrieved using function `f`.
"""
function scale_rate!(tree::T, f::Function, s::Float64) where {T <: cT}

  f(tree, s)

  if def1(tree)
    scale_rate!(tree.d1, f, s)
    if def2(tree)
      scale_rate!(tree.d2, f, s)
    end
  end

  return nothing
end




"""
    scale_rate!(tree::iTree, f::Function, s::Float64)

Add `s` to vector retrieved using function `f`.
"""
function scale_rate!(tree::iTree, f::Function, s::Float64)

  v = f(tree)
  @turbo for i in Base.OneTo(lastindex(v))
    v[i] += s
  end

  if def1(tree)
    scale_rate!(tree.d1, f, s)
    if def2(tree)
      scale_rate!(tree.d2, f, s)
    end
  end

  return nothing
end




"""
    reorder(tree::T) where {T <: iTree}

Reorder order of daughter branches according to number of tips, with daughter
1 always having more than daughter 2.
"""
function reorder!(tree::T) where {T <: iTree}
  tree, n = _reorder!(tree)
  return tree
end




"""
    _reorder!(tree::T) where {T <: iTree}

Reorder order of daughter branches according to number of tips, with daughter
1 always having more than daughter 2.
"""
function _reorder!(tree::T) where {T <: iTree}

  if istip(tree)
    n = 1
  else
    t1, n1 = _reorder!(tree.d1)

    if def2(tree)
      t2, n2 = _reorder!(tree.d2)
      if n1 < n2
        tree.d1 = t2
        tree.d2 = t1
      end
      n = n1 + n2
    else
      n = n1
    end
  end
  return tree, n
end




"""
    reorder!(tree::T, treeda::D) where {T <: iTree, D <: iTree}

Reorder a data augmented `treeda` tree with `tree`.
"""
function reorder!(tree::T, treeda::D) where {T <: iTree, D <: iTree}

  if def1(tree)
    if def2(tree)

      t1  = tree.d1
      td1 = treeda.d1
      e1  = e(tree.d1)
      ed1 = fixedge(treeda.d1)

      if !isapprox(e1, ed1)
        treeda.d1 = treeda.d2
        treeda.d2 = td1
      end

      reorder!(tree.d1, fixtree(treeda.d1))
      reorder!(tree.d2, fixtree(treeda.d2))
    else
      reorder!(tree.d1, fixtree(treeda.d1))
    end  
  end
end




"""
    rm_stem(tree::T) where {T <: iTree}

Removes stem branch.
"""
rm_stem!(tree::T) where {T <: iTree} =
  _rm_stem(tree)




"""
    _rm_stem(tree::T) where {T <: sT}

Remove stem branch.
"""
_rm_stem(tree::T) where {T <: sT} = sete!(tree, 0.0)




"""
    _rm_stem(tree::T) where {T <: iT}

Remove stem branch.
"""
function _rm_stem(tree::T) where {T <: iT}

  sete!(  tree, 0.0)
  setfdt!(tree, 0.0)
  lλe = lλ(tree)[end]
  setlλ!( tree, [lλe,lλe])

  return tree
end




"""
    _rm_stem(tree::iTbd)

Remove stem branch.
"""
function _rm_stem(tree::iTbdU)

  sete!(  tree, 0.0)
  setfdt!(tree, 0.0)
  lλe = lλ(tree)[end]
  lμe = lμ(tree)[end]
  setlλ!( tree, [lλe, lλe])
  setlμ!( tree, [lμe, lμe])

  return tree
end





"""
    prune_tips(tree::T, tips::Vector{String}) where {T <: Tlabel}

Prune `tips` from tree.
"""
function prune_tips(tree::T, tips::Vector{String}) where {T <: Tlabel}
  tips = Set(tips)

  return _prune_tips!(T(tree::T), tips)
end




"""
    prune_tips(treev::Vector{T}, tips::Vector{String}) where {T <: Tlabel}

Prune `tips` from tree.
"""
function prune_tips(treev::Vector{T}, tips::Vector{String}) where {T <: Tlabel}
  treevne = T[]
  for t in treev
    push!(treevne, prune_tips(t, tips))
  end
  return treevne
end




"""
    _prune_tips!(tree::T, tips::Vector{String}) where {T <: Tlabel}

Prune `tips` from tree.
"""
function _prune_tips!(tree::T, tips::Set{String}) where {T <: iTf}

  if def1(tree)
    tree.d1 = _prune_tips!(tree.d1, tips)
    if def2(tree)
      tree.d2 = _prune_tips!(tree.d2, tips)
      if in(label(tree.d1), tips)
        setdiff!(tips, label(tree.d1))
        if in(label(tree.d2), tips)
          return T(e(tree), label(tree.d2))
        else
          ne  = e(tree) + e(tree.d2)
          tree = tree.d2
          sete!(tree, ne)
        end
      elseif in(label(tree.d2), tips)
        setdiff!(tips, label(tree.d2))
        ne  = e(tree) + e(tree.d1)
        tree = tree.d1
        sete!(tree, ne)
      end
      return tree
    else
      if in(label(tree.d1), tips)
        setdiff!(tips, label(tree.d1))
        return T(e(tree), label(tree))
      end
    end
  end

  return tree
end





"""
    cutbottom(tree::T, c::Float64) where {T <: iTree}

Cut the bottom part of the tree after `c`.
"""
cutbottom(tree::T, c::Float64) where {T <: iTree} = _cutbottom(tree, c, 0.0)




"""
    _cutbottom(tree::sTb,
               c   ::Float64,
               t   ::Float64)

Cut the bottom part of the tree after `c`, starting at time `t`.
"""
function _cutbottom(tree::sTb,
                    c   ::Float64,
                    t   ::Float64)

  et = e(tree)

  if (t + et) > c
    tree = sTb(c - t)
  else
    if def1(tree)
      tree.d1 = _cutbottom(tree.d1, c, t + et)
      tree.d2 = _cutbottom(tree.d2, c, t + et)
    end
  end

  return tree
end




"""
    _cutbottom(tree::sTbd,
               c   ::Float64,
               t   ::Float64)

Cut the bottom part of the tree after `c`, starting at time `t`.
"""
function _cutbottom(tree::sTbd,
                    c   ::Float64,
                    t   ::Float64)

  et = e(tree)

  if (t + et) > c
    tree = sTbd(c - t, false, isfix(tree))
  else
    if def1(tree)
      tree.d1 = _cutbottom(tree.d1, c, t + et)
      tree.d2 = _cutbottom(tree.d2, c, t + et)
    end
  end

  return tree
end




"""
    _cutbottom(tree::sTfbd,
               c   ::Float64,
               t   ::Float64)

Cut the bottom part of the tree after `c`, starting at time `t`.
"""
function _cutbottom(tree::sTfbd,
                    c   ::Float64,
                    t   ::Float64)

  et = e(tree)

  if (t + et) > c
    tree = sTfbd(c - t, false, false, isfix(tree))
  else
    if def1(tree) tree.d1 = _cutbottom(tree.d1, c, t + et) end
    if def2(tree) tree.d2 = _cutbottom(tree.d2, c, t + et) end
  end

  return tree
end




"""
    _cutbottom(tree::iTb,
               c   ::Float64,
               t   ::Float64)

Cut the bottom part of the tree after `c`, starting at time `t`.
"""
function _cutbottom(tree::iTb,
                    c   ::Float64,
                    t   ::Float64)

  et = e(tree)

  if (t + et) > c

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

    lλv = lλv[1:Ix]

    push!(lλv, eλ)

    tree = iTb(c - t, δt, c - t - tii, true, lλv)

  else
    if def1(tree)
      tree.d1 = _cutbottom(tree.d1, c, t + et)
      tree.d2 = _cutbottom(tree.d2, c, t + et)
    end
  end

  return tree
end




"""
    _cutbottom(tree::T,
               c   ::Float64,
               t   ::Float64) where {T <: iT}

Cut the bottom part of the tree after `c`, starting at time `t`.
"""
function _cutbottom(tree::T,
                    c   ::Float64,
                    t   ::Float64) where {T <: iT}

  et = e(tree)

  if (t + et) > c

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

    lλv = lλv[1:Ix]

    push!(lλv, eλ)

    tree = T(c - t, δt, c - t - tii, false, isfix(tree), lλv)

  else
    if def1(tree)
      tree.d1 = _cutbottom(tree.d1, c, t + et)
      tree.d2 = _cutbottom(tree.d2, c, t + et)
    end
  end

  return tree
end




"""
    _cutbottom(tree::iTbd,
               c   ::Float64,
               t   ::Float64)

Cut the bottom part of the tree after `c`, starting at time `t`.
"""
function _cutbottom(tree::iTbd,
                    c   ::Float64,
                    t   ::Float64)

  et = e(tree)

  if (t + et) > c

    lλv = lλ(tree)
    lμv = lμ(tree)
    δt  = dt(tree)
    fδt = fdt(tree)

    # find final lλ & lμ
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
    eμ = linpred(c - t, tii, tff, lμv[Ix], lμv[Ix+1])

    lλv = lλv[1:Ix]
    lμv = lμv[1:Ix]

    push!(lλv, eλ)
    push!(lμv, eμ)

    tree = iTbd(c - t, δt, c - t - tii, false, isfix(tree), lλv, lμv)

  else
    if def1(tree)
      tree.d1 = _cutbottom(tree.d1, c, t + et)
      tree.d2 = _cutbottom(tree.d2, c, t + et)
    end
  end

  return tree
end




"""
    _cutbottom(tree::iTfbd,
               c   ::Float64,
               t   ::Float64)

Cut the bottom part of the tree after `c`, starting at time `t`.
"""
function _cutbottom(tree::iTfbd,
                    c   ::Float64,
                    t   ::Float64)

  et = e(tree)

  if (t + et) > c

    lλv = lλ(tree)
    lμv = lμ(tree)
    δt  = dt(tree)
    fδt = fdt(tree)

    # find final lλ & lμ
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
    eμ = linpred(c - t, tii, tff, lμv[Ix], lμv[Ix+1])

    lλv = lλv[1:Ix]
    lμv = lμv[1:Ix]

    push!(lλv, eλ)
    push!(lμv, eμ)

    tree = iTfbd(c - t, δt, c - t - tii, false, false, isfix(tree), lλv, lμv)

  else
    if def1(tree)
      tree.d1 = _cutbottom(tree.d1, c, t + et)
      if def2(tree)
        tree.d2 = _cutbottom(tree.d2, c, t + et)
      end
    end
  end

  return tree
end



"""
    fossilizefixedtip!(tree::T) where {T <: iTf}

Change all alive tips to fossil tips.
"""
function fossilizefixedtip!(tree::T) where {T <: iTf}

  if istip(tree)
    fossilize!(tree)
  elseif isfix(tree.d1)
    fossilizefixedtip!(tree.d1::T)
  else
    fossilizefixedtip!(tree.d2::T)
  end
end




"""
    fossilizepasttips!(tree::T) where {T <: iTf}

Change all past tips to fossil tips.
"""
fossilizepasttips!(tree::T, ne::Float64) where {T <: iTf} = 
  _fossilizepasttips!(tree::T, treeheight(tree), ne)




"""
    _fossilizepasttips!(tree::T,
                        t   ::Float64,
                        ne  ::Float64) where {T <: iTf}

Change all past tips to fossil tips, initialized at tree height `t`.
"""
function _fossilizepasttips!(tree::T,
                             t   ::Float64,
                             ne  ::Float64) where {T <: iTf}

  t -= e(tree)

  if def1(tree)
    _fossilizepasttips!(tree.d1::T, t, ne)
    if def2(tree)
      _fossilizepasttips!(tree.d2::T, t, ne)
    end
  else
    if t > ne
      fossilize!(tree)
    end
  end
end




"""
    fossilize!(tree::T, sp_label::String) where {T <: iTf}

Fossilize a given tree given its name.
"""
function fossilize!(tree::T, sp_label::String) where {T <: sTf_label}

  if def1(tree)
    fossilize!(tree.d1::T, sp_label)
    if def2(tree)
      fossilize!(tree.d2::T, sp_label)
    end
  elseif label(tree) == sp_label 
    fossilize!(tree)
  end
end




"""
    defossilize!(tree::T, sp_label::String) where {T <: iTf}

Fossilize a given tree given its name.
"""
function defossilize!(tree::T, sp_label::String) where {T <: sTf_label}

  if def1(tree)
    defossilize!(tree.d1::T, sp_label)
    if def2(tree)
      defossilize!(tree.d2::T, sp_label)
    end
  elseif label(tree) == sp_label 
    defossilize!(tree)
  end
end




"""
    extinguishpasttips!(tree::T) where {T <: iTf}

Change all past tips to extinct tips.
"""
extinguishpasttips!(tree::T, ne::Float64) where {T <: iTf} = 
  _extinguishpasttips!(tree::T, treeheight(tree), ne)




"""
    _extinguishpasttips!(tree::T,
                         t   ::Float64,
                         ne  ::Float64) where {T <: iTf}

Change all past tips to extinct tips, initialized at tree height `t`.
"""
function _extinguishpasttips!(tree::T,
                              t   ::Float64,
                              ne  ::Float64) where {T <: iTf}

  t -= e(tree)

  if def1(tree)
    _extinguishpasttips!(tree.d1::T, t, ne)
    if def2(tree)
      _extinguishpasttips!(tree.d2::T, t, ne)
    end
  else
    if t > ne
      extinguish!(tree)
    end
  end
end




"""
    extinguish!(tree::T, label::String) where {T <: iTf}

Extinguish a given tree given its name.
"""
function extinguish!(tree::T, label::String) where {T <: sTf_label}

  if def1(tree)
    extinguish!(tree.d1::T, label)
    if def2(tree)
      extinguish!(tree.d2::T, label)
    end
  elseif l(tree) == label 
    extinguish!(tree)
  end
end





"""
    fixrtip!(tree::T, na::Int64, λf::Float64) where {T <: iT}

Fixes the the path for a random non extinct tip and returns final `λ(t)`.
"""
function fixrtip!(tree::T, na::Int64, λf::Float64) where {T <: iT}

  fix!(tree)

  if def1(tree)
    if isextinct(tree.d1)
      λf = fixrtip!(tree.d2, na, λf)
    elseif isextinct(tree.d2)
      λf = fixrtip!(tree.d1, na, λf)
    else
      na1 = ntipsalive(tree.d1)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        λf = fixrtip!(tree.d2, na - na1, λf)
      else
        λf = fixrtip!(tree.d1, na1,      λf)
      end
    end
  else
    λf = lλ(tree)[end]
  end

  return λf
end




"""
    fixrtip!(tree::T, na::Int64, xt::Float64) where T <: Tx

Fixes the the path for a random non extinct tip and returns final `λ(t)`.
"""
function fixrtip!(tree::T, na::Int64, xt::Float64) where T <: Tx

  fix!(tree)

  if def1(tree)
    if isextinct(tree.d1)
      xt = fixrtip!(tree.d2, na, xt)
    elseif isextinct(tree.d2)
      xt = fixrtip!(tree.d1, na, xt)
    else
      na1 = ntipsalive(tree.d1)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        xt = fixrtip!(tree.d2, na - na1, xt)
      else
        xt = fixrtip!(tree.d1, na1,      xt)
      end
    end
  else
    xt = xf(tree)
  end

  return xt
end




"""
    fixrtip!(tree::T, na::Int64) where T <: Tx

Fixes the the path for a random non extinct tip and returns final tip.
"""
function fixrtip!(tree::T, na::Int64) where T <: Tx

  fix!(tree)

  if def1(tree)
    if isextinct(tree.d1)
      fixrtip!(tree.d2, na)
    elseif isextinct(tree.d2)
      fixrtip!(tree.d1, na)
    else
      na1 = ntipsalive(tree.d1)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        fixrtip!(tree.d2, na - na1)
      else
        fixrtip!(tree.d1, na1)
      end
    end
  else
    return tree
  end
end





"""
    fixrtip!(tree::iTbd,
             na  ::Int64,
             λf  ::Float64,
             μf  ::Float64)

Fixes the the path for a random non extinct tip.
"""
function fixrtip!(tree::T,
                  na  ::Int64,
                  λf  ::Float64,
                  μf  ::Float64) where {T <: iTbdU}

  fix!(tree)

  if def1(tree)
    if isextinct(tree.d1)
      λf, μf =
        fixrtip!(tree.d2, na, λf, μf)
    elseif isextinct(tree.d2)
      λf, μf =
        fixrtip!(tree.d1, na, λf, μf)
    else
      na1 = ntipsalive(tree.d1)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        λf, μf =
          fixrtip!(tree.d2, na - na1, λf, μf)
      else
        λf, μf =
          fixrtip!(tree.d1, na1, λf, μf)
      end
    end
  else

    λv   = lλ(tree)
    l    = lastindex(λv)
    λf   = λv[l]
    μf   = lμ(tree)[l]
  end

  return λf, μf
end




"""
    fixalive!(tree::T, λf::Float64) where {T <:iT}

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::T, λf::Float64) where {T <:iT}

  if istip(tree)
    if isalive(tree)
      fix!(tree)
      λf   = lλ(tree)[end]
      return true, λf
    end
  else
    f, λf = fixalive!(tree.d2, λf)
    if f
      fix!(tree)
      return true, λf
    end
    f, λf = fixalive!(tree.d1, λf)
    if f
      fix!(tree)
      return true, λf
    end
  end

  return false, λf
end




"""
    fixalive!(tree::T) where T <: iTree
Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::T) where T <: iTree

  if istip(tree::T)
    if isalive(tree::T)
      fix!(tree::T)
      return true
    end
  else
    f = fixalive!(tree.d2::T)
    if f
      fix!(tree)
      return true
    end
    f = fixalive!(tree.d1::T)
    if f
      fix!(tree)
      return true
    end
  end

  return false
end




"""
    fixalive!(tree::T, xt::Float64) where T <: iTree

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::T, xt::Float64) where T <: Tx

  if istip(tree::T)
    if isalive(tree::T)
      fix!(tree::T)
      xt = xf(tree)
      return true, xt
    end
  else
    f, xt = fixalive!(tree.d2::T, xt)
    if f
      fix!(tree)
      return true, xt
    end
    f, xt = fixalive!(tree.d1::T, xt)
    if f
      fix!(tree)
      return true, xt
    end
  end

  return false, xt
end




"""
    fixalive!(tree::iTbd,
              λf  ::Float64,
              μf  ::Float64)

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::T,
                   λf  ::Float64,
                   μf  ::Float64) where {T <: iTbdU}

  if istip(tree)
    if isalive(tree)
      fix!(tree)
      λv   = lλ(tree)
      μv   = lμ(tree)
      l    = lastindex(λv)
      λf   = λv[l]
      μf   = lμ(tree)[l]

      return true, λf, μf
    end
  else
    f, λf, μf =
      fixalive!(tree.d2, λf, μf)
    if f
      fix!(tree)
      return true, λf, μf
    end
    f, λf, μf =
      fixalive!(tree.d1, λf, μf)
    if f
      fix!(tree)
      return true, λf, μf
    end
  end

  return false, λf, μf
end



"""
    fixalive!(tree::T) where {T <: iTf}

Fixes the path from root to the only species alive.
"""
function fixalive!(tree::T) where {T <: iTf}

  if istip(tree::T) && isalive(tree::T)
    fix!(tree::T)
    return true
  end

  if def2(tree)
    f = fixalive!(tree.d2::T)
    if f
      fix!(tree)
      return true
    end
  end

  if def1(tree)
    f = fixalive!(tree.d1::T)
    if f
      fix!(tree)
      return true
    end
  end

  return false
end




"""
    fixrtip!(tree::T) where T <: iTree

Fixes the the path for a random non extinct tip.
"""
fixrtip!(tree::T) where T <: iTree = _fixrtip!(tree, ntipsalive(tree))



"""
    _fixrtip!(tree::T, na::Int64) where {T <: iT}

Fixes the the path for a random non extinct tip among `na`.
"""
function _fixrtip!(tree::T, na::Int64) where T <: iTree

  fix!(tree)

  if def2(tree)
    if isextinct(tree.d1)
      _fixrtip!(tree.d2, na)
    elseif isextinct(tree.d2)
      _fixrtip!(tree.d1, na)
    else
      na1 = ntipsalive(tree.d1)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        _fixrtip!(tree.d2, na - na1)
      else
        _fixrtip!(tree.d1, na1)
      end
    end
  end
end



"""
    _fixrtip!(tree::T, na::Int64) where {T <: iT}

Fixes the the path for a random non extinct tip among `na`.
"""
function _fixfossrtip!(tree::T, na::Int64) where T <: iTf

  fix!(tree)

  if def2(tree)
    if isextinct(tree.d1)
      _fixfossrtip!(tree.d2, na)
    elseif isextinct(tree.d2)
      _fixfossrtip!(tree.d1, na)
    else
      na1 = ntipsalive(tree.d1)
      # probability proportional to number of lineages
      if (fIrand(na) + 1) > na1
        _fixfossrtip!(tree.d2, na - na1)
      else
        _fixfossrtip!(tree.d1, na1)
      end
    end
  else
    tree.iψ = true
  end
end




"""
    fixtipsρ!(tree::T, ρ::Float64) where {T <: iTree}
Fixes the the path to each extant tip with probability `ρ`.
"""
function fixtipsρ!(tree::T, ρ::Float64) where {T <: iTree}

  if istip(tree)
    if isalive(tree)
      if rand()<ρ
        fix!(tree)
        return true
      end
    end
  
  else
    f1 = fixtipsρ!(tree.d1, ρ)
    f2 = fixtipsρ!(tree.d2, ρ)
    if f1 || f2
      fix!(tree)
      return true
    end
  end

  return false
end




"""
    fixtip1!(tree::T, wi::Int64, ix::Int64) where {T <: iTree}
Fixes the the path to tip `wi` in d1 order.
"""
function fixtip1!(tree::T, wi::Int64, ix::Int64) where {T <: iTree}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        return true, ix
      end
    end
  else
    f, ix = fixtip1!(tree.d1, wi, ix)
    if f
      fix!(tree)
      return true, ix
    end
    f, ix = fixtip1!(tree.d2, wi, ix)
    if f
      fix!(tree)
      return true, ix
    end
  end

  return false, ix
end




"""
    fixtip2!(tree::T, wi::Int64, ix::Int64) where {T <: iTree}

Fixes the the path to tip `wi` in d2 order.
"""
function fixtip2!(tree::T, wi::Int64, ix::Int64) where {T <: iTree}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        return true, ix
      end
    end
  else
    f, ix = fixtip2!(tree.d2, wi, ix)
    if f
      fix!(tree)
      return true, ix
    end
    f, ix = fixtip2!(tree.d1, wi, ix)
    if f
      fix!(tree)
      return true, ix
    end
  end

  return false, ix
end




"""
    fixtip1!(tree::T, wi::Int64, ix::Int64) where {T <: iTree}

Fixes the the path to tip `wi` in d1 order.
"""
function fixtip1!(tree::T, wi::Int64, ix::Int64, xc::Float64) where {T <: iTree}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        setxf!(tree, xc)
        return true, ix
      end
    end
  else
    f, ix = fixtip1!(tree.d1, wi, ix, xc)
    if f
      fix!(tree)
      return true, ix
    end
    f, ix = fixtip1!(tree.d2, wi, ix, xc)
    if f
      fix!(tree)
      return true, ix
    end
  end

  return false, ix
end




"""
    fixtip2!(tree::T, wi::Int64, ix::Int64) where {T <: iTree}

Fixes the the path to tip `wi` in d2 order.
"""
function fixtip2!(tree::T, wi::Int64, ix::Int64, xc::Float64) where {T <: iTree}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        setxf!(tree, xc)
        return true, ix
      end
    end
  else
    f, ix = fixtip2!(tree.d2, wi, ix, xc)
    if f
      fix!(tree)
      return true, ix
    end
    f, ix = fixtip2!(tree.d1, wi, ix, xc)
    if f
      fix!(tree)
      return true, ix
    end
  end

  return false, ix
end





"""
    fixtip1!(tree::T, wi::Int64, ix::Int64, s::Bool) where {T <: iTree}

Fixes the the path to tip `wi` in d1 order.
"""
function fixtip1!(tree::T, wi::Int64, ix::Int64, s::Bool) where {T <: iTree}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        setsh!(tree, s)
        return true, ix
      end
    end
  else
    f, ix = fixtip1!(tree.d1, wi, ix, s)
    if f
      fix!(tree)
      return true, ix
    end
    f, ix = fixtip1!(tree.d2, wi, ix, s)
    if f
      fix!(tree)
      return true, ix
    end
  end

  return false, ix
end




"""
    fixtip2!(tree::T, wi::Int64, ix::Int64, s::Bool) where {T <: iTree}

Fixes the the path to tip `wi` in d2 order.
"""
function fixtip2!(tree::T, wi::Int64, ix::Int64, s::Bool) where {T <: iTree}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        setsh!(tree, s)
        return true, ix
      end
    end
  else
    f, ix = fixtip2!(tree.d2, wi, ix, s)
    if f
      fix!(tree)
      return true, ix
    end
    f, ix = fixtip2!(tree.d1, wi, ix, s)
    if f
      fix!(tree)
      return true, ix
    end
  end

  return false, ix
end




"""
    fixtip1!(tree::T, 
             wi  ::Int64, 
             ix  ::Int64, 
             λa  ::Float64, 
             ei  ::Float64, 
             λi  ::Float64) where {T <: cT}

Fixes the the path to tip `wi` in d1 order.
"""
function fixtip1!(tree::T, 
                  wi  ::Int64, 
                  ix  ::Int64, 
                  λa  ::Float64, 
                  ei  ::Float64, 
                  λi  ::Float64) where {T <: cT}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        setlλ!(tree, λi)
        return true, ix, λa, e(tree)
      end
    end
  else
    f, ix, λa, ei = fixtip1!(tree.d1, wi, ix, lλ(tree), ei, λi)
    if f
      fix!(tree)
      return true, ix, λa, ei
    end
    f, ix, λa, ei = fixtip1!(tree.d2, wi, ix, lλ(tree), ei, λi)
    if f
      fix!(tree)
      return true, ix, λa, ei
    end
  end

  return false, ix, λa, ei
end




"""
    fixtip2!(tree::T, 
             wi  ::Int64, 
             ix  ::Int64, 
             λa  ::Float64, 
             ei  ::Float64, 
             λi  ::Float64) where {T <: cT}

Fixes the the path to tip `wi` in d2 order.
"""
function fixtip2!(tree::T, 
                  wi  ::Int64, 
                  ix  ::Int64, 
                  λa  ::Float64, 
                  ei  ::Float64, 
                  λi  ::Float64) where {T <: cT}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        setlλ!(tree, λi)
        return true, ix, λa, e(tree)
      end
    end
  else
    f, ix, λa, ei = fixtip2!(tree.d2, wi, ix, lλ(tree), ei, λi)
    if f
      fix!(tree)
      return true, ix, λa, ei
    end
    f, ix, λa, ei = fixtip2!(tree.d1, wi, ix, lλ(tree), ei, λi)
    if f
      fix!(tree)
      return true, ix, λa, ei
    end
  end

  return false, ix, λa, ei
end




"""
    fixtip1!(tree::T,
             wi  ::Int64,
             ix  ::Int64,
             xc  ::Float64,
             σx  ::Float64,
             δt  ::Float64,
             srδt::Float64) where {T <: Tx}

Fixes the the path to tip `wi` in d1 order.
"""
function fixtip1!(tree::T,
                  wi  ::Int64,
                  ix  ::Int64,
                  xc  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  srδt::Float64) where {T <: Tx}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        xvi = xv(tree)
        bb!(xvi, xvi[1], xc, σx, δt, fdt(tree), srδt)

        return true, ix
      end
    end
  else
    f, ix = fixtip1!(tree.d1, wi, ix, xc, σx, δt, srδt)
    if f
      fix!(tree)
      return true, ix
    end
    f, ix = fixtip1!(tree.d2, wi, ix, xc, σx, δt, srδt)
    if f
      fix!(tree)
      return true, ix
    end
  end

  return false, ix
end




"""
    fixtip2!(tree::T,
             wi  ::Int64,
             ix  ::Int64,
             xc  ::Float64,
             σx  ::Float64,
             δt  ::Float64,
             srδt::Float64,
             llr ::Float64,
             ssrx::Float64) where {T <: Tx}

Fixes the the path to tip `wi` in d2 order.
"""
function fixtip2!(tree::T,
                  wi  ::Int64,
                  ix  ::Int64,
                  xc  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  srδt::Float64) where {T <: Tx}

  if istip(tree)
    if isalive(tree)
      ix += 1
      if ix === wi
        fix!(tree)
        xvi = xv(tree)
        bb!(xvi, xvi[1], xc, σx, δt, fdt(tree), srδt)

        return true, ix
      end
    end
  else
    f, ix = fixtip2!(tree.d2, wi, ix, xc, σx, δt, srδt)
    if f
      fix!(tree)
      return true, ix
    end
    f, ix = fixtip2!(tree.d1, wi, ix, xc, σx, δt, srδt)
    if f
      fix!(tree)
      return true, ix
    end
  end

  return false, ix
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

  if ix === ldr
    if thc > h > (thc - e(tree))
      ne = thc - h
      adde!(tree, -ne)
      tree = rand() <= 0.5 ? T(tree, stree, ne, false, true) :
                             T(stree, tree, ne, false, true)
    else
      if isfix(tree.d1)
        tree.d1 =
          graftree!(tree.d1::T, stree, dri, h, ldr, thc - e(tree), ix)
      else
        tree.d2 =
          graftree!(tree.d2::T, stree, dri, h, ldr, thc - e(tree), ix)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      ix += 1
      if dri[ix]
        tree.d1 =
          graftree!(tree.d1::T, stree, dri, h, ldr, thc - e(tree), ix)
      else
        tree.d2 =
          graftree!(tree.d2::T, stree, dri, h, ldr, thc - e(tree), ix)
      end
    elseif ifx1
      tree.d1 =
        graftree!(tree.d1::T, stree, dri, h, ldr, thc - e(tree), ix)
    else
      tree.d2 =
        graftree!(tree.d2::T, stree, dri, h, ldr, thc - e(tree), ix)
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

  if ix === ldr
    if px === wpr
      if isfix(tree.d1::T)
        ne  = e(tree) + e(tree.d1)
        sete!(tree.d1, ne)
        tree = tree.d1
      elseif isfix(tree.d2::T)
        ne  = e(tree) + e(tree.d2)
        sete!(tree.d2, ne)
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
    remove_unsampled(tree::T) where {T <: iTree}

Remove unsampled tips from `iTree`.
"""
function remove_unsampled(tree::T) where {T <: iTree}
  return _remove_unsampled!(T(tree::T))
end




"""
    remove_unsampled(treev::Vector{T}) where {T <: iTree}

Remove unsampled taxa for a vector of trees.
"""
function remove_unsampled(treev::Vector{T}) where {T <: iTree}

  treevne = T[]
  for t in treev
    push!(treevne, remove_unsampled(t))
  end

  return treevne
end




"""
    _remove_unsampled!(tree::sTbd)

Remove unsampled tips (extinct and extant not sampled).
"""
function _remove_unsampled!(tree::sTbd)

  if def1(tree)
    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return sTbd(e(tree), isextinct(tree), isfix(tree))
      else
        ne  = e(tree) + e(tree.d2)
        tree = tree.d2
        sete!(tree, ne)
      end
    elseif !isfix(tree.d2)
      ne  = e(tree) + e(tree.d1)
      tree = tree.d1
      sete!(tree, ne)
    end
  end

  return tree
end




"""
    _remove_unsampled!(tree::sTpe)

Remove unsampled tips (extinct and extant not sampled).
"""
function _remove_unsampled!(tree::sTpe)

  if def1(tree)
    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return sTpe(e(tree), isextinct(tree), xi(tree), xf(tree), sh(tree), isfix(tree))
      else
        ne = e(tree) + e(tree.d2)
        nx = xi(tree)
        tree = tree.d2
        sete!(tree, ne)
        setxi!(tree, nx)
      end
    elseif !isfix(tree.d2)
      ne = e(tree) + e(tree.d1)
      nx = xi(tree)
      tree = tree.d1
      sete!(tree, ne)
      setxi!(tree, nx)
    end
  end

  return tree
end



"""
    _remove_unsampled!(tree::sTfpe)

Remove unsampled tips (extinct and extant not sampled).
"""
function _remove_unsampled!(tree::sTfpe)

  if def1(tree)
    tree.d1 = _remove_unsampled!(tree.d1)
    if def2(tree)
      tree.d2 = _remove_unsampled!(tree.d2)
      if !isfix(tree.d1)
        if !isfix(tree.d2)
          return sTfpe(e(tree), isextinct(tree), isfossil(tree), xi(tree), xf(tree), sh(tree), isfix(tree))
        else
          ne = e(tree) + e(tree.d2)
          nx = xi(tree)
          tree = tree.d2
          sete!(tree, ne)
          setxi!(tree, nx)
        end
      elseif !isfix(tree.d2)
        ne = e(tree) + e(tree.d1)
        nx = xi(tree)
        tree = tree.d1
        sete!(tree, ne)
        setxi!(tree, nx)
      end
    else
      if !isfix(tree.d1)
        return sTfpe(e(tree), isextinct(tree), isfossil(tree), xi(tree), xf(tree), sh(tree), isfix(tree))
      end
    end
  end

  return tree
end




"""
    _remove_unsampled!(tree::T) where {T <: iTf}

Remove unsampled tips.
"""
function _remove_unsampled!(tree::T) where {T <: iTf}

  if def1(tree)
    tree.d1 = _remove_unsampled!(tree.d1)
    if def2(tree)
      tree.d2 = _remove_unsampled!(tree.d2)

      if !isfix(tree.d1)
        if !isfix(tree.d2)
          return T(e(tree), isextinct(tree), isfossil(tree), isfix(tree))
        else
          ne  = e(tree) + e(tree.d2)
          tree = tree.d2
          sete!(tree, ne)
        end
      elseif !isfix(tree.d2)
        ne  = e(tree) + e(tree.d1)
        tree = tree.d1
        sete!(tree, ne)
      end
      return tree
    end

    if !isfix(tree.d1)
      return T(e(tree), isextinct(tree), isfossil(tree), isfix(tree))
    end
  end

  return tree
end




"""
    _remove_unsampled!(tree::cTb)

Remove extinct tips from `cTb`.
"""
function _remove_unsampled!(tree::cTb)

  if def1(tree)

    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return cTb(e(tree), isfix(tree), lλ(tree))
      else
        e0   = e(tree)
        e2   = e(tree.d2)
        λ0 = lλ(tree)
        λ2 = lλ(tree.d2)
        tree = tree.d2
        sete!(tree, e0 + e2)
        setlλ!(tree, ((e0*λ0) + (e2*λ2)) /(e0 + e2))
      end
    elseif !isfix(tree.d2)
      e0   = e(tree)
      e1   = e(tree.d1)
      λ0 = lλ(tree)
      λ1 = lλ(tree.d1)
      tree = tree.d1
      sete!(tree, e0 + e1)
      setlλ!(tree, ((e0*λ0) + (e1*λ1)) /(e0 + e1))
    end
  end
  return tree
end




"""
    _remove_unsampled!(tree::T) where {T <: cT}

Remove extinct tips from `cTce`.
"""
function _remove_unsampled!(tree::T) where {T <: cT}

  if def1(tree)

    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return T(e(tree), isextinct(tree), isfix(tree), lλ(tree))
      else
        e0   = e(tree)
        e2   = e(tree.d2)
        λ0 = lλ(tree)
        λ2 = lλ(tree.d2)
        tree = tree.d2
        sete!(tree, e0 + e2)
        setlλ!(tree, ((e0*λ0) + (e2*λ2)) /(e0 + e2))
      end
    elseif !isfix(tree.d2)
      e0   = e(tree)
      e1   = e(tree.d1)
      λ0 = lλ(tree)
      λ1 = lλ(tree.d1)
      tree = tree.d1
      sete!(tree, e0 + e1)
      setlλ!(tree, ((e0*λ0) + (e1*λ1)) /(e0 + e1))
    end
  end
  return tree
end




"""
    _remove_unsampled!(tree::cTbd)

Remove extinct tips from `cTbd`.
"""
function _remove_unsampled!(tree::cTbd)

  if def1(tree)

    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return cTbd(e(tree), isextinct(tree), isfix(tree), lλ(tree), lμ(tree))
      else
        e0   = e(tree)
        e2   = e(tree.d2)
        λ0 = lλ(tree)
        λ2 = lλ(tree.d2)
        μ0 = lμ(tree)
        μ2 = lμ(tree.d2)
        tree = tree.d2
        sete!(tree, e0 + e2)
        setlλ!(tree, ((e0*λ0) + (e2*λ2)) /(e0 + e2))
        setlμ!(tree, ((e0*μ0) + (e2*μ2)) /(e0 + e2))
      end
    elseif !isfix(tree.d2)
      e0   = e(tree)
      e1   = e(tree.d1)
      λ0 = lλ(tree)
      λ1 = lλ(tree.d1)
      μ0 = lμ(tree)
      μ1 = lμ(tree.d1)
      tree = tree.d1
      sete!(tree, e0 + e1)
      setlλ!(tree, ((e0*λ0) + (e1*λ1)) /(e0 + e1))
      setlμ!(tree, ((e0*μ0) + (e1*μ1)) /(e0 + e1))
    end
  end
  return tree
end




"""
    _remove_unsampled!(tree::iTb)

Remove extinct tips from `iTb`.
"""
function _remove_unsampled!(tree::iTb)

  if def1(tree)

    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return iTb(e(tree), dt(tree), fdt(tree), isfix(tree), lλ(tree))
      else
        ne  = e(tree) + e(tree.d2)
        lλ0 = lλ(tree)
        lλ2 = lλ(tree.d2)

        fdt2 = fdt(tree.d2)
        pop!(lλ0)
        iszero(fdt2) && popfirst!(lλ2)
        prepend!(lλ2, lλ0)
        fdt0 = fdt(tree) + fdt2
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d2
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    elseif !isfix(tree.d2)
      ne  = e(tree) + e(tree.d1)
      lλ0 = lλ(tree)
      lλ1 = lλ(tree.d1)

      fdt1 = fdt(tree.d1)
      pop!(lλ0)
      iszero(fdt1) && popfirst!(lλ1)
      prepend!(lλ1, lλ0)
      fdt0 = fdt(tree) + fdt1
      if fdt0 > dt(tree)
        fdt0 -= dt(tree)
      end
      tree = tree.d1
      sete!(tree, ne)
      setfdt!(tree, fdt0)
    end
  end

  return tree
end




"""
    _remove_unsampled!(tree::iTce)

Remove unsampled tips from `iTce`.
"""
function _remove_unsampled!(tree::iTce)

  if def1(tree)

    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return iTce(e(tree), dt(tree), fdt(tree),
          isextinct(tree), isfix(tree), lλ(tree))
      else
        ne  = e(tree) + e(tree.d2)
        lλ0 = lλ(tree)
        lλ2 = lλ(tree.d2)

        fdt2 = fdt(tree.d2)
        pop!(lλ0)
        iszero(fdt2) && popfirst!(lλ2)
        prepend!(lλ2, lλ0)
        fdt0 = fdt(tree) + fdt2
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d2
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    elseif !isfix(tree.d2)
      ne  = e(tree) + e(tree.d1)
      lλ0 = lλ(tree)
      lλ1 = lλ(tree.d1)

      fdt1 = fdt(tree.d1)
      pop!(lλ0)
      iszero(fdt1) && popfirst!(lλ1)
      prepend!(lλ1, lλ0)
      fdt0 = fdt(tree) + fdt1
      if fdt0 > dt(tree)
        fdt0 -= dt(tree)
      end
      tree = tree.d1
      sete!(tree, ne)
      setfdt!(tree, fdt0)
    end
  end

  return tree
end




"""
    _remove_unsampled!(tree::iTct)
Remove extinct tips from `iTct`.
"""
function _remove_unsampled!(tree::iTct)

  if def1(tree)

    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return iTct(e(tree), dt(tree), fdt(tree),
          isextinct(tree), isfix(tree), lλ(tree))
      else
        ne  = e(tree) + e(tree.d2)
        lλ0 = lλ(tree)
        lλ2 = lλ(tree.d2)

        fdt2 = fdt(tree.d2)
        pop!(lλ0)
        iszero(fdt2) && popfirst!(lλ2)
        prepend!(lλ2, lλ0)
        fdt0 = fdt(tree) + fdt2
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d2
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    elseif !isfix(tree.d2)
      ne  = e(tree) + e(tree.d1)
      lλ0 = lλ(tree)
      lλ1 = lλ(tree.d1)

      fdt1 = fdt(tree.d1)
      pop!(lλ0)
      iszero(fdt1) && popfirst!(lλ1)
      prepend!(lλ1, lλ0)
      fdt0 = fdt(tree) + fdt1
      if fdt0 > dt(tree)
        fdt0 -= dt(tree)
      end
      tree = tree.d1
      sete!(tree, ne)
      setfdt!(tree, fdt0)
    end
  end

  return tree
end




"""
    _remove_unsampled!(tree::iTbd)

Remove extinct tips from `iTbd`.
"""
function _remove_unsampled!(tree::iTbd)

  if def1(tree)

    tree.d1 = _remove_unsampled!(tree.d1)
    tree.d2 = _remove_unsampled!(tree.d2)

    if !isfix(tree.d1)
      if !isfix(tree.d2)
        return iTbd(e(tree), dt(tree), fdt(tree),
          isextinct(tree), isfix(tree), lλ(tree), lμ(tree))
      else
        ne = e(tree) + e(tree.d2)
        lλ0 = lλ(tree)
        lμ0 = lμ(tree)
        lλ2 = lλ(tree.d2)
        lμ2 = lμ(tree.d2)
        fdt2 = fdt(tree.d2)
        pop!(lλ0)
        pop!(lμ0)
        if iszero(fdt2)
          popfirst!(lλ2)
          popfirst!(lμ2)
        end
        prepend!(lλ2, lλ0)
        prepend!(lμ2, lμ0)

        fdt0 = fdt(tree) + fdt(tree.d2)
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d2
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    elseif !isfix(tree.d2)
      ne = e(tree) + e(tree.d1)
      lλ0 = lλ(tree)
      lμ0 = lμ(tree)
      lλ1 = lλ(tree.d1)
      lμ1 = lμ(tree.d1)
      fdt1 = fdt(tree.d1)
      pop!(lλ0)
      pop!(lμ0)
      if iszero(fdt1)
        popfirst!(lλ1)
        popfirst!(lμ1)
      end
      prepend!(lλ1, lλ0)
      prepend!(lμ1, lμ0)

      fdt0 = fdt(tree) + fdt(tree.d1)
      if fdt0 > dt(tree)
        fdt0 -= dt(tree)
      end
      tree = tree.d1
      sete!(tree, ne)
      setfdt!(tree, fdt0)
    end
  end

  return tree
end





"""
    _remove_unsampled!(tree::iTfbd)

Remove extinct tips from `iTfbd`.
"""
function _remove_unsampled!(tree::iTfbd)

  if def1(tree)
    tree.d1 = _remove_unsampled!(tree.d1)
    if def2(tree)
      tree.d2 = _remove_unsampled!(tree.d2)
      if !isfix(tree.d1)
        if !isfix(tree.d2)
          return iTfbd(e(tree), dt(tree), fdt(tree),
            isextinct(tree), isfossil(tree), isfix(tree), lλ(tree), lμ(tree))
        else
          ne = e(tree) + e(tree.d2)
          lλ0 = lλ(tree)
          lμ0 = lμ(tree)
          lλ2 = lλ(tree.d2)
          lμ2 = lμ(tree.d2)
          fdt2 = fdt(tree.d2)
          pop!(lλ0)
          pop!(lμ0)
          if iszero(fdt2)
            popfirst!(lλ2)
            popfirst!(lμ2)
          end
          prepend!(lλ2, lλ0)
          prepend!(lμ2, lμ0)

          fdt0 = fdt(tree) + fdt(tree.d2)
          if fdt0 > dt(tree)
            fdt0 -= dt(tree)
          end
          tree = tree.d2
          sete!(tree, ne)
          setfdt!(tree, fdt0)
        end
      elseif !isfix(tree.d2)
        ne = e(tree) + e(tree.d1)
        lλ0 = lλ(tree)
        lμ0 = lμ(tree)
        lλ1 = lλ(tree.d1)
        lμ1 = lμ(tree.d1)
        fdt1 = fdt(tree.d1)
        pop!(lλ0)
        pop!(lμ0)
        if iszero(fdt1)
          popfirst!(lλ1)
          popfirst!(lμ1)
        end
        prepend!(lλ1, lλ0)
        prepend!(lμ1, lμ0)

        fdt0 = fdt(tree) + fdt(tree.d1)
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d1
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    else
      if !isfix(tree.d1)
          return iTfbd(e(tree), dt(tree), fdt(tree),
            isextinct(tree), isfossil(tree), isfix(tree), lλ(tree), lμ(tree))
      end
    end
  end

  return tree
end




"""
    remove_extinct(tree::T) where {T <: iTree}

Remove extinct tips from `iTce`.
"""
function remove_extinct(tree::T) where {T <: iTree}
  return _remove_extinct!(T(tree::T))
end




"""
    remove_extinct(treev::Vector{T}) where {T <: iTree}

Remove extinct taxa for a vector of trees.
"""
function remove_extinct(treev::Vector{T}) where {T <: iTree}

  treevne = T[]
  for t in treev
    push!(treevne, remove_extinct(t))
  end

  return treevne
end




"""
    _remove_extinct!(tree::sTbd)

Remove extinct tips (except fossil tips).
"""
function _remove_extinct!(tree::sTbd)

  if def1(tree)
    tree.d1 = _remove_extinct!(tree.d1)
    tree.d2 = _remove_extinct!(tree.d2)

    if isextinct(tree.d1)
      if isextinct(tree.d2)
        return sTbd(e(tree), true, isfix(tree))
      else
        ne  = e(tree) + e(tree.d2)
        tree = tree.d2
        sete!(tree, ne)
      end
    elseif isextinct(tree.d2)
      ne  = e(tree) + e(tree.d1)
      tree = tree.d1
      sete!(tree, ne)
    end
  end

  return tree
end




"""
    _remove_extinct!(tree::sTpe)

Remove extinct tips (except fossil tips).
"""
function _remove_extinct!(tree::sTpe)

  if def1(tree)
    tree.d1 = _remove_extinct!(tree.d1)
    tree.d2 = _remove_extinct!(tree.d2)

    if isextinct(tree.d1)
      if isextinct(tree.d2)
        return sTpe(e(tree), true, xi(tree), xf(tree), sh(tree), isfix(tree))
      else
        ne   = e(tree) + e(tree.d2)
        xii  = xi(tree)
        tree = tree.d2
        sete!(tree, ne)
        setxi!(tree, xii)
      end
    elseif isextinct(tree.d2)
      ne   = e(tree) + e(tree.d1)
      xii  = xi(tree)
      tree = tree.d1
      sete!(tree, ne)
      setxi!(tree, xii)
    end
  end

  return tree
end




"""
    _remove_extinct!(tree::T) where {T <: iTf}

Remove unsampled tips.
"""
function _remove_extinct!(tree::T) where {T <: iTf}

  if def1(tree)
    tree.d1 = _remove_extinct!(tree.d1)
    if def2(tree)
      tree.d2 = _remove_extinct!(tree.d2)

      if isextinct(tree.d1)
        if isextinct(tree.d2)
          return T(e(tree), true, isfossil(tree), isfix(tree))
        else
          ne  = e(tree) + e(tree.d2)
          tree = tree.d2
          sete!(tree, ne)
        end
      elseif isextinct(tree.d2)
        ne  = e(tree) + e(tree.d1)
        tree = tree.d1
        sete!(tree, ne)
      end
      return tree
    else
      if isextinct(tree.d1)
        return T(e(tree), false, true, isfix(tree))
      end
    end
  end

  return tree
end




"""
    _remove_extinct!(tree::sTfpe)

Remove extinct tips (except fossil tips).
"""
function _remove_extinct!(tree::sTfpe)

  if def1(tree)
    tree.d1 = _remove_extinct!(tree.d1)
    if def2(tree)
      tree.d2 = _remove_extinct!(tree.d2)

      if isextinct(tree.d1)
        if isextinct(tree.d2)
          return sTfpe(e(tree), true, isfossil(tree), 
                   xi(tree), xf(tree), sh(tree), isfix(tree))
        else
          ne   = e(tree) + e(tree.d2)
          xii  = xi(tree)
          tree = tree.d2
          sete!(tree, ne)
          setxi!(tree, xii)
        end
      elseif isextinct(tree.d2)
        ne   = e(tree) + e(tree.d1)
        xii  = xi(tree)
        tree = tree.d1
        sete!(tree, ne)
        setxi!(tree, xii)
      end
    else
      if isextinct(tree.d1)
        return sTfpe(e(tree), false, true, 
                xi(tree), xf(tree), sh(tree), isfix(tree))
      end
    end
  end

  return tree
end




"""
    _remove_extinct!(tree::iTce)

Remove extinct tips from `iTce`.
"""
function _remove_extinct!(tree::iTce)

  if def1(tree)

    tree.d1 = _remove_extinct!(tree.d1)
    tree.d2 = _remove_extinct!(tree.d2)

    if isextinct(tree.d1)
      if isextinct(tree.d2)
        return iTce(e(tree), dt(tree), fdt(tree),
          true, isfix(tree), lλ(tree))
      else
        ne  = e(tree) + e(tree.d2)
        lλ0 = lλ(tree)
        lλ2 = lλ(tree.d2)

        fdt2 = fdt(tree.d2)
        pop!(lλ0)
        iszero(fdt2) && popfirst!(lλ2)
        prepend!(lλ2, lλ0)
        fdt0 = fdt(tree) + fdt2
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d2
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    elseif isextinct(tree.d2)
      ne  = e(tree) + e(tree.d1)
      lλ0 = lλ(tree)
      lλ1 = lλ(tree.d1)

      fdt1 = fdt(tree.d1)
      pop!(lλ0)
      iszero(fdt1) && popfirst!(lλ1)
      prepend!(lλ1, lλ0)
      fdt0 = fdt(tree) + fdt1
      if fdt0 > dt(tree)
        fdt0 -= dt(tree)
      end
      tree = tree.d1
      sete!(tree, ne)
      setfdt!(tree, fdt0)
    end
  end

  return tree
end




"""
    _remove_extinct!(tree::iTct)

Remove extinct tips from `iTct`.
"""
function _remove_extinct!(tree::iTct)

  if def1(tree)

    tree.d1 = _remove_extinct!(tree.d1)
    tree.d2 = _remove_extinct!(tree.d2)

    if isextinct(tree.d1)
      if isextinct(tree.d2)
        return iTct(e(tree), dt(tree), fdt(tree),
          true, isfix(tree), lλ(tree))
      else
        ne  = e(tree) + e(tree.d2)
        lλ0 = lλ(tree)
        lλ2 = lλ(tree.d2)

        fdt2 = fdt(tree.d2)
        pop!(lλ0)
        iszero(fdt2) && popfirst!(lλ2)
        prepend!(lλ2, lλ0)
        fdt0 = fdt(tree) + fdt2
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d2
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    elseif isextinct(tree.d2)
      ne  = e(tree) + e(tree.d1)
      lλ0 = lλ(tree)
      lλ1 = lλ(tree.d1)

      fdt1 = fdt(tree.d1)
      pop!(lλ0)
      iszero(fdt1) && popfirst!(lλ1)
      prepend!(lλ1, lλ0)
      fdt0 = fdt(tree) + fdt1
      if fdt0 > dt(tree)
        fdt0 -= dt(tree)
      end
      tree = tree.d1
      sete!(tree, ne)
      setfdt!(tree, fdt0)
    end
  end

  return tree
end




"""
    _remove_extinct!(tree::iTbd)

Remove extinct tips from `iTbd`.
"""
function _remove_extinct!(tree::iTbd)

  if def1(tree)

    tree.d1 = _remove_extinct!(tree.d1)
    tree.d2 = _remove_extinct!(tree.d2)

    if isextinct(tree.d1)
      if isextinct(tree.d2)
        return iTbd(e(tree), dt(tree), fdt(tree),
          true, isfix(tree), lλ(tree), lμ(tree))
      else
        ne = e(tree) + e(tree.d2)
        lλ0 = lλ(tree)
        lμ0 = lμ(tree)
        lλ2 = lλ(tree.d2)
        lμ2 = lμ(tree.d2)
        fdt2 = fdt(tree.d2)
        pop!(lλ0)
        pop!(lμ0)
        if iszero(fdt2)
          popfirst!(lλ2)
          popfirst!(lμ2)
        end
        prepend!(lλ2, lλ0)
        prepend!(lμ2, lμ0)

        fdt0 = fdt(tree) + fdt(tree.d2)
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d2
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    elseif isextinct(tree.d2)
      ne = e(tree) + e(tree.d1)
      lλ0 = lλ(tree)
      lμ0 = lμ(tree)
      lλ1 = lλ(tree.d1)
      lμ1 = lμ(tree.d1)
      fdt1 = fdt(tree.d1)
      pop!(lλ0)
      pop!(lμ0)
      if iszero(fdt1)
        popfirst!(lλ1)
        popfirst!(lμ1)
      end
      prepend!(lλ1, lλ0)
      prepend!(lμ1, lμ0)

      fdt0 = fdt(tree) + fdt(tree.d1)
      if fdt0 > dt(tree)
        fdt0 -= dt(tree)
      end
      tree = tree.d1
      sete!(tree, ne)
      setfdt!(tree, fdt0)
    end
  end

  return tree
end





"""
    _remove_extinct!(tree::iTfbd)

Remove extinct tips from `iTfbd`.
"""
function _remove_extinct!(tree::iTfbd)

  if def1(tree)
    tree.d1 = _remove_extinct!(tree.d1)

    if def2(tree)
      tree.d2 = _remove_extinct!(tree.d2)

      if isextinct(tree.d1)
        if isextinct(tree.d2)
          return iTfbd(e(tree), dt(tree), fdt(tree),
            true, isfossil(tree), isfix(tree), lλ(tree), lμ(tree))
        else
          ne = e(tree) + e(tree.d2)
          lλ0 = lλ(tree)
          lμ0 = lμ(tree)
          lλ2 = lλ(tree.d2)
          lμ2 = lμ(tree.d2)
          fdt2 = fdt(tree.d2)
          pop!(lλ0)
          pop!(lμ0)
          if iszero(fdt2)
            popfirst!(lλ2)
            popfirst!(lμ2)
          end
          prepend!(lλ2, lλ0)
          prepend!(lμ2, lμ0)

          fdt0 = fdt(tree) + fdt(tree.d2)
          if fdt0 > dt(tree)
            fdt0 -= dt(tree)
          end
          tree = tree.d2
          sete!(tree, ne)
          setfdt!(tree, fdt0)
        end
      elseif isextinct(tree.d2)
        ne = e(tree) + e(tree.d1)
        lλ0 = lλ(tree)
        lμ0 = lμ(tree)
        lλ1 = lλ(tree.d1)
        lμ1 = lμ(tree.d1)
        fdt1 = fdt(tree.d1)
        pop!(lλ0)
        pop!(lμ0)
        if iszero(fdt1)
          popfirst!(lλ1)
          popfirst!(lμ1)
        end
        prepend!(lλ1, lλ0)
        prepend!(lμ1, lμ0)

        fdt0 = fdt(tree) + fdt(tree.d1)
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d1
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    else
      if isextinct(tree.d1)
        return iTfbd(e(tree), dt(tree), fdt(tree),
            isextinct(tree), isfossil(tree), isfix(tree), lλ(tree), lμ(tree))
      end
    end
  end

  return tree
end





"""
    _remove_extinct!(tree::iTpbd)

Remove extinct tips from `iTpbd`.
"""
function _remove_extinct!(tree::iTpbd)

  if def1(tree)
    tree.d1 = _remove_extinct!(tree.d1)

    if def2(tree)
      tree.d2 = _remove_extinct!(tree.d2)

      if isextinct(tree.d1)
        if isextinct(tree.d2)
          return iTpbd(e(tree), dt(tree), fdt(tree), true, iscomplete(tree), 
            isfix(tree), lb(tree), lλ(tree), lμ(tree))
        else
          ne = e(tree) + e(tree.d2)
          lb0 = lb(tree)
          lλ0 = lλ(tree)
          lμ0 = lμ(tree)
          lb2 = lb(tree.d2)
          lλ2 = lλ(tree.d2)
          lμ2 = lμ(tree.d2)
          fdt2 = fdt(tree.d2)
          pop!(lb0)
          pop!(lλ0)
          pop!(lμ0)
          if iszero(fdt2)
            popfirst!(lb2)
            popfirst!(lλ2)
            popfirst!(lμ2)
          end
          prepend!(lb2, lb0)
          prepend!(lλ2, lλ0)
          prepend!(lμ2, lμ0)

          fdt0 = fdt(tree) + fdt(tree.d2)
          if fdt0 > dt(tree)
            fdt0 -= dt(tree)
          end
          tree = tree.d2
          sete!(tree, ne)
          setfdt!(tree, fdt0)
        end
      elseif isextinct(tree.d2)
        ne = e(tree) + e(tree.d1)
        lb0 = lb(tree)
        lλ0 = lλ(tree)
        lμ0 = lμ(tree)
        lb1 = lb(tree.d1)
        lλ1 = lλ(tree.d1)
        lμ1 = lμ(tree.d1)
        fdt1 = fdt(tree.d1)
        pop!(lb0)
        pop!(lλ0)
        pop!(lμ0)
        if iszero(fdt1)
          popfirst!(lb1)
          popfirst!(lλ1)
          popfirst!(lμ1)
        end
        prepend!(lb1, lb0)
        prepend!(lλ1, lλ0)
        prepend!(lμ1, lμ0)

        fdt0 = fdt(tree) + fdt(tree.d1)
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d1
        sete!(tree, ne)
        setfdt!(tree, fdt0)
      end
    else
      if isextinct(tree.d1)
        return iTpbd(e(tree), dt(tree), fdt(tree), true, 
          iscomplete(tree), isfix(tree), lb(tree), lλ(tree), lμ(tree))
      end
    end
  end

  return tree
end





"""
    remove_fossils(tree::T) where {T <: iTf}

Remove fossils.
"""
function remove_fossils(tree::T) where {T <: iTf}
  return _remove_fossils!(T(tree::T))
end




"""
    _remove_fossils!(tree::T) where {T <: iTf}

Remove fossils.
"""
function _remove_fossils!(tree::T) where {T <: iTf}

  if def1(tree)
    tree.d1 = _remove_fossils!(tree.d1)
    if def2(tree)
      tree.d2 = _remove_fossils!(tree.d2)
      if isfossil(tree.d2)
        adde!(tree.d1, e(tree))
        tree = tree.d1
      elseif isfossil(tree.d1)
        adde!(tree.d2, e(tree))
        tree = tree.d2
      end
    else
      t1 = _remove_fossils!(tree.d1)
      #defossilize!(t1)
      adde!(t1, e(tree))
      return t1
    end
  else
    #isfossil(tree) && defossilize!(tree)
  end

  return tree
end



"""
    _remove_fossils!(tree::iTfbd)

Remove fossils.
"""
function _remove_fossils!(tree::iTfbd)

  if def1(tree)
    if def2(tree)
      tree.d1 = _remove_fossils!(tree.d1)
      tree.d2 = _remove_fossils!(tree.d2)
    else
      t1 = _remove_fossils!(tree.d1)
      defossilize!(t1)
      lλ1 = lλ(t1)
      lμ1 = lμ(t1)
      dti = dt(tree)
      e1  = e(t1)
      t0  = 0.0
      tn  = dti - fdt(tree)
      i   = 1
      while e1 > tn + dti + accerr
        lλ1[i] = linpred(tn, t0, t0 + dti, lλ1[i], lλ1[i+1])
        lμ1[i] = linpred(tn, t0, t0 + dti, lμ1[i], lμ1[i+1])
        tn += dti
        t0 += dti
        i  += 1
      end
      if fdt(tree) < dti
        lλ1[i] = lλ1[i+1]
        lμ1[i] = lμ1[i+1]
        if (e1 - t0) > dti + accerr || e1 < tn
          pop!(lλ1)
          pop!(lμ1)
        end
      end
      lλv = lλ(tree)
      lμv = lμ(tree)
      pop!(lλv)
      pop!(lμv)
      prepend!(lλ1, lλv)
      prepend!(lμ1, lμv)
      adde!(t1, e(tree))
      if tn < e1
        setfdt!(t1, e1 - tn)
      else
        setfdt!(t1, e1 + fdt(tree))
      end

      return t1
    end
  else
    isfossil(tree) && defossilize!(tree)
  end

  return tree
end



"""
    prune_fossils(treev::Vector{sTf_label})

Prune fossils.
"""
function prune_fossils(treev::Vector{sTf_label})

  treevne = sTf_label[]
  for t in treev
    push!(treevne, prune_fossils(t))
  end

  return treevne
end




"""
    prune_fossils(tree::sTf_label)

Prune fossils.
"""
function prune_fossils(tree::sTf_label)
  return _prune_fossils!(sTf_label(tree::sTf_label))
end



"""
    _prune_fossils!(tree::sTf_label)

Prune fossils.
"""
function _prune_fossils!(tree::sTf_label)

  if def1(tree)
    tree.d1 = _prune_fossils!(tree.d1)
    if def2(tree)
      tree.d2 = _prune_fossils!(tree.d2)

      if isfossil(tree.d1)
        if isfossil(tree.d2)
          return sTf_label(e(tree), isextinct(tree), true, label(tree))
        else
          ne  = e(tree) + e(tree.d2)
          tree = tree.d2
          sete!(tree, ne)
        end
      elseif isfossil(tree.d2)
        ne  = e(tree) + e(tree.d1)
        tree = tree.d1
        sete!(tree, ne)
      end
      return tree
    else
        ne  = e(tree) + e(tree.d1)
        tree = tree.d1
        sete!(tree, ne)
    end
  end

  return tree
end




"""
    remove_sampled_ancestors(tree::T, p::Float64) where {T <: iTf}

Remove sampled ancestor fossils with a certain probability `p`.
"""
function remove_sampled_ancestors(tree::T, p::Float64) where {T <: iTf}
  t, i = _remove_sampled_ancestors!(T(tree::T),  p, false)
  return t
end




"""
    _remove_sampled_ancestors!(tree::T) where {T <: iTf}
Remove sampled ancestor fossils with a certain probability `p`.
"""
function _remove_sampled_ancestors!(tree::T,  p::Float64, i::Bool) where {T <: iTf}

  if def1(tree)
    if def2(tree)
      tree.d1, i = _remove_sampled_ancestors!(tree.d1, p, i)
      tree.d2, i = _remove_sampled_ancestors!(tree.d2, p, i)
    else
      t1, i = _remove_sampled_ancestors!(tree.d1, p, i)

      if i || rand() < p
        adde!(t1, e(tree))
        return t1, false
      else
        tree.d1 = t1
        return tree, true
      end
    end
  end

  return tree, i
end




"""
    remove_incipient(tree::T) where {T <: iTpbd}

Remove incipient lineages.
"""
function remove_incipient(tree::T) where {T <: iTpbd}
  tree_resampled = _resample_representative!(T(tree))
  return _remove_incipient!(tree_resampled::T)
end





"""
    _remove_incipient!(tree::iTpbd, representative::Bool)

Randomly sample a representative lineage among `na` lineages alive.
"""
function _resample_representative!(tree::iTpbd; representative::Bool=iscomplete(tree))
  
  representative ? setgood!(tree) : setincipient!(tree)

  if def1(tree)
    if def2(tree)
      if representative
        # count the number of alive lineages of the same species for each daughter
        na1 = ntipsalivespecies(tree.d1)
        na2 = ntipsalivespecies(tree.d2)

        # if all descendants speciate or go extinct
        if iszero(na1) && iszero(na2)
          tree.d1 = _resample_representative!(tree.d1)
          tree.d2 = _resample_representative!(tree.d2)
        # otherwise uniformly sample a representative daughter
        else
          if rand() < na1/(na1 + na2)
            tree.d1 = _resample_representative!(tree.d1, representative=true)
            tree.d2 = _resample_representative!(tree.d2, representative=false)
          else
            tree.d1 = _resample_representative!(tree.d1, representative=false)
            tree.d2 = _resample_representative!(tree.d2, representative=true)
          end
        end
      else
        # daughters of a non-representative lineage are also non-representative
        tree.d1 = _resample_representative!(tree.d1, representative=false)
        tree.d2 = _resample_representative!(tree.d2, representative=false)
      end
    
    # completion event: creates a new representative lineage
    else
      tree.d1 = _resample_representative!(tree.d1, representative=true)
    end
  end

  return tree
end




"""
    _remove_incipient!(tree::iTpbd)

Remove extinct tips from `iTpbd`.
"""
function _remove_incipient!(tree::iTpbd)

  if def1(tree)
    tree.d1 = _remove_incipient!(tree.d1)
    if def2(tree)
      tree.d2 = _remove_incipient!(tree.d2)

      # the parent of two good lineages is also good
      if iscomplete(tree.d1)
        if iscomplete(tree.d2)
          setgood!(tree)
        end
      
      # if the second daughter only is good, remove the first one
      elseif iscomplete(tree.d2)
        ne = e(tree) + e(tree.d2)
        lb0 = lb(tree)
        lλ0 = lλ(tree)
        lμ0 = lμ(tree)
        lb2 = lb(tree.d2)
        lλ2 = lλ(tree.d2)
        lμ2 = lμ(tree.d2)
        fdt2 = fdt(tree.d2)
        pop!(lb0)
        pop!(lλ0)
        pop!(lμ0)
        if iszero(fdt2)
          popfirst!(lb2)
          popfirst!(lλ2)
          popfirst!(lμ2)
        end
        prepend!(lb2, lb0)
        prepend!(lλ2, lλ0)
        prepend!(lμ2, lμ0)

        fdt0 = fdt(tree) + fdt(tree.d2)
        if fdt0 > dt(tree)
          fdt0 -= dt(tree)
        end
        tree = tree.d2
        sete!(tree, ne)
        setfdt!(tree, fdt0)        
      end
    end
  end
  if def1(tree)
    # if the first daughter only is good
    if (def2(tree) && !iscomplete(tree.d2)) || !def2(tree)
      ne = e(tree) + e(tree.d1)
      lb0 = lb(tree)
      lλ0 = lλ(tree)
      lμ0 = lμ(tree)
      lb1 = lb(tree.d1)
      lλ1 = lλ(tree.d1)
      lμ1 = lμ(tree.d1)
      fdt1 = fdt(tree.d1)
      pop!(lb0)
      pop!(lλ0)
      pop!(lμ0)
      if iszero(fdt1)
        popfirst!(lb1)
        popfirst!(lλ1)
        popfirst!(lμ1)
      end
      prepend!(lb1, lb0)
      prepend!(lλ1, lλ0)
      prepend!(lμ1, lμ0)

      fdt0 = fdt(tree) + fdt(tree.d1)
      if fdt0 > dt(tree)
        fdt0 -= dt(tree)
      end
      tree = tree.d1
      sete!(tree, ne)
      setfdt!(tree, fdt0)
    end
  end

  return tree
end




"""
    reconstructed(tree::T) where {T <: iTree}
    reconstructed(tree::T) where {T <: iTf}
Returns the reconstructed tree, i.e. the observed tree from sampled extant
tips and fossils.
"""
# For all trees without fossils, it simply means removing extinct lineages
reconstructed(tree::T) where {T <: iTree} = remove_unsampled(tree::T)
reconstructed(tree::T) where {T <: iTf} = _reconstructed!(T(tree::T))




"""
    _reconstructed!(tree::T) where {T <: iTf}
Returns the reconstructed tree, i.e. the observed tree from sampled extant
tips and fossils.
"""
function _reconstructed!(tree::T) where {T <: iTf}
  defd1 = def1(tree)
  defd2 = def2(tree)

  if defd1 tree.d1 = _reconstructed!(tree.d1) end
  if defd2 tree.d2 = _reconstructed!(tree.d2) end

  if !defd1 && !defd2
    return tree
  end

  extd1 = defd1 && isextinct(tree.d1)
  extd2 = defd2 && isextinct(tree.d2)

  # 2 extinct daughters -> extinct tip
  if extd1 && extd2
    return T(e(tree), true, false, isfix(tree))
  end

  # sampled ancestor with extinct daughter -> fossil tip (labelled extinct)
  if extd1 && !defd2
    lab = extd1 ? l(tree.d1) : l(tree.d2)
    if isnothing(lab)
      return T(e(tree), false, true, isfix(tree))
    else
      return T(e(tree), false, true, isfix(tree), lab)
    end
  end

  # 1 extinct and 1 alive branch -> keep only the alive one
  if extd1 ne = e(tree)+e(tree.d2); tree = tree.d2; sete!(tree,ne) end
  if extd2 ne = e(tree)+e(tree.d1); tree = tree.d1; sete!(tree,ne) end

  return tree
end




"""
    fixtree!(tree::T) where {T <: iTree}

Fix all `tree`.
"""
function fixtree!(tree::T) where {T <: iTree}
  fix!(tree)
  if def1(tree)
    fixtree!(tree.d1)
    fixtree!(tree.d2)
  end
end




"""
    fixtree!(tree::T) where {T <: iTf}

Fix all `tree`.
"""
function fixtree!(tree::T) where {T <: iTf}
  fix!(tree)
  if def1(tree) fixtree!(tree.d1) end
  if def2(tree) fixtree!(tree.d2) end
end




"""
    fix!(tree::T) where {T <: iTree}

Fix `tree`.
"""
fix!(tree::T) where {T <: iTree} = setproperty!(tree, :fx, true)




"""
    setlλ!(tree::T, lλ::Array{Float64,1}) where {T <: iT}

Set speciation `lλ` for `tree`.
"""
setlλ!(tree::T, lλ::Array{Float64,1}) where {T <: iT} =
  setproperty!(tree, :lλ, lλ)
setlλ!(tree::T, lλ::Float64) where {T <: cT} = setproperty!(tree, :lλ, lλ)




"""
    setlμ!(tree::T, lμ::Array{Float64,1}) where {T <: iTbd}

Set number of `δt` for `tree`.
"""
setlμ!(tree::T, lμ::Array{Float64,1}) where {T <: iT} =
  setproperty!(tree, :lμ, lμ)
setlμ!(tree::T, lμ::Float64) where {T <: cT} = setproperty!(tree, :lμ, lμ)




"""
    addlλ!(tree::cTb, a::Float64)

Add to `lλ` for `tree`.
"""
addlλ!(tree::T, a::Float64) where {T <: cT} = setproperty!(tree, :lλ, lλ(tree) + a)




"""
    addlμ!(tree::cTb, a::Float64)

Add to `lμ` for `tree`.
"""
addlμ!(tree::T, a::Float64) where {T <: cT} = setproperty!(tree, :lμ, lμ(tree) + a)




"""
  setfdt!(tree::T, fdt::Float64) where {T <: iTree}

Set number of `δt` for `tree`.
"""
setfdt!(tree::T, fdt::Float64) where {T <: iTree} =
  setproperty!(tree, :fdt, fdt)




"""
  setdt!(tree::T, dt::Float64) where {T <: iTree}

Set number of `δt` for `tree`.
"""
setdt!(tree::T, dt::Float64) where {T <: iTree} =
  setproperty!(tree, :dt, dt)




"""
  sete!(tree::T, e::Float64) where {T <: iTree}

Set endant edge for `tree`.
"""
sete!(tree::T, e::Float64) where {T <: iTree} = setproperty!(tree, :e, e)




"""
  adde!(tree::T, e::Float64) where {T <: iTree}

Add `e` to edge of `tree`.
"""
adde!(tree::T, e::Float64) where {T <: iTree} = tree.e += e




"""
  setd1!(tree::T,  stree::T) where {T <: iTree}

Set `d1` to `stree` in `tree`.
"""
setd1!(tree::T,  stree::T) where {T <: iTree} = setproperty!(tree, :d1, stree)




"""
  setd2!(tree::T,  stree::T) where {T <: iTree}

Set `d2` to `stree` in `tree`.
"""
setd2!(tree::T,  stree::T) where {T <: iTree} = setproperty!(tree, :d2, stree)




"""
  fossilize!(tree::iTf)

Make tree a fossil.
"""
fossilize!(tree::iTf) = setproperty!(tree, :iψ, true)




"""
  defossilize!(tree::iTf)

Make tree a non fossil.
"""
defossilize!(tree::iTf) = setproperty!(tree, :iψ, false)



"""
  extinguish!(tree::iTf)

Make tree a non fossil.
"""
extinguish!(tree::iTf) = setproperty!(tree, :iμ, true)



"""
  setgood!(tree::iTpbd)

Make tree a good/representative lineage.
"""
setgood!(tree::iTpbd) = setproperty!(tree, :ig, true)



"""
  setincipient!(tree::iTpbd)

Make tree an incipient lineage.
"""
setincipient!(tree::iTpbd) = setproperty!(tree, :ig, false)



"""
  setxi!(tree::T, x::Float64) where {T <: Tx}

Set `x` as initial trait in tree.
"""
setxi!(tree::T, x::Float64) where {T <: Tx} =
  setproperty!(tree, :xi, x)




"""
  setxf!(tree::T, x::Float64) where {T <: Tx}

Set `x` as final trait in tree.
"""
setxf!(tree::T, x::Float64) where {T <: Tx} =
  setproperty!(tree, :xf, x)




"""
  setsh!(tree::T, x::Bool) where {T <: Tx} =

Set `x` as shift to d1 (true) or d2 (false).
"""
setsh!(tree::T, x::Bool) where {T <: Tx} =
  setproperty!(tree, :sh, x)




"""
    setupstreamλ!(λi ::Float64,
                  i  ::Int64,
                  Ξ  ::Vector{T},
                  idf::Vector{iBffs}) where {T <: cT}

Set the speciation rate of the upstream ancestors, if any, for 
middle branches to `λi`.
"""
function setupstreamλ!(λi ::Float64,
                       i  ::Int64,
                       Ξ  ::Vector{T},
                       idf::Vector{iBffs}) where {T <: cT}

  @inbounds begin
    bi = idf[i]

    # if branch is non-cladogenetic
    if iszero(d2(bi))
      ξi   = Ξ[i]

      if def2(ξi)
        lξi = fixtip(ξi)
        setlλ!(lξi, λi)
      else
        setlλ!(ξi, λi)
        setupstreamλ!(λi, pa(bi), Ξ, idf)
      end
    end
  end

  return nothing
end




"""
    setupstreamλμ!(λi ::Float64,
                   μi ::Float64,
                   i  ::Int64,
                   Ξ  ::Vector{cTbd},
                   idf::Vector{iBffs})

Set the speciation and extinction rate of the upstream ancestors, if any, for 
middle branches to `λi`.
"""
function setupstreamλμ!(λi ::Float64,
                        μi ::Float64,
                        i  ::Int64,
                        Ξ  ::Vector{cTbd},
                        idf::Vector{iBffs})

  @inbounds begin
    bi = idf[i]

    # if branch is non-cladogenetic
    if iszero(d2(bi))
      ξi   = Ξ[i]

      if def2(ξi)
        lξi = fixtip(ξi)
        setlλ!(lξi, λi)
        setlμ!(lξi, μi)
      else
        setlλ!(ξi, λi)
        setlμ!(ξi, μi)
        setupstreamλμ!(λi, μi, pa(bi), Ξ, idf)
      end
    end
  end

  return nothing
end




"""
    setdownstreamλ!(λi ::Float64,
                    i  ::Int64,
                    Ξ  ::Vector{T}, 
                    idf::Vector{iBffs}) where {T <: cT}

Set the speciation rate of the downstream daughters, if any, for 
middle branches to `λi`.
"""
function setdownstreamλ!(λi ::Float64,
                         i  ::Int64,
                         Ξ  ::Vector{T}, 
                         idf::Vector{iBffs}) where {T <: cT}

  @inbounds begin

    ξi = Ξ[i]
    setlλ!(ξi, λi)

    if istip(ξi)
      bi = idf[i]
      i1 = d1(bi)
      if i1 > 0
        if iszero(d2(bi))
          setdownstreamλ!(λi, i1, Ξ, idf)
        else
          setλt!(bi, λi)
        end
      end
    end
  end

  return nothing
end





"""
    setdownstreamλμ!(λi ::Float64,
                     μi ::Float64,
                     i  ::Int64,
                     Ξ  ::Vector{cTbd}, 
                     idf::Vector{iBffs})

Set the speciation rate of the downstream daughters, if any, for 
middle branches to `λi` and `μi`.
"""
function setdownstreamλμ!(λi ::Float64,
                          μi ::Float64,
                          i  ::Int64,
                          Ξ  ::Vector{cTbd}, 
                          idf::Vector{iBffs})

  @inbounds begin

    ξi = Ξ[i]
    setlλ!(ξi, λi)
    setlμ!(ξi, μi)

    if istip(ξi)
      bi = idf[i]
      i1 = d1(bi)
      if i1 > 0
        if iszero(d2(bi))
          setdownstreamλ!(λi, μi, i1, Ξ, idf)
        else
          setλt!(bi, λi)
          setμt!(bi, μi)
        end
      end
    end
  end

  return nothing
end



