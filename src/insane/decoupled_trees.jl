#=

decoupled tree utilities

Ignacio Quintero Mächler

t(-_-t)

Created 05 11 2020
=#



"""
    copy_Ξ(Ξ::Vector{T}) where {T <: iTree}

Make edge tree `Ξ` from the edge directory.
"""
function copy_Ξ(Ξ::Vector{T}) where {T <: iTree}
  Ξc = T[]
  for ξ in Ξ
    push!(Ξc, T(ξ))
  end
  return Ξc
end




"""
    sTpb!(Ξ::Vector{sTpb}, tree::sT_label)

Make edge tree `Ξ` from the recursive tree.
"""
function sTpb!(Ξ::Vector{sTpb}, tree::sT_label)

  push!(Ξ, sTpb(e(tree), true))
  if def1(tree)
    sTpb!(Ξ, tree.d2)
    sTpb!(Ξ, tree.d1)
  end
end




"""
    make_Ξ(idf::Vector{iBffs}, xr::Vector{Float64}, ::Type{sTpbX})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, xr::Vector{Float64}, ::Type{sTpbX})
  Ξ = sTpbX[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    paix = iszero(paix) ? 1 : paix
    xii  = xr[paix]
    xfi  = xr[i]
    ξ = sTpbX(e(idfi), true, xii, xfi)
    push!(Ξ, ξ)
  end

  return Ξ
end




"""
    make_Ξ(idf::Vector{iBffs}, xr::Vector{Float64}, ::Type{sTbdX})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, xr::Vector{Float64}, ::Type{sTbdX})
  Ξ = sTbdX[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    paix = iszero(paix) ? 1 : paix
    xii  = xr[paix]
    xfi  = xr[i]
    ξ = sTbdX(e(idfi), false, true, xii, xfi)
    push!(Ξ, ξ)
  end

  return Ξ
end




"""
    sTbd!(Ξ::Vector{sTbd}, tree::sT_label)

Make edge tree `Ξ` from the recursive tree.
"""
function sTbd!(Ξ::Vector{sTbd}, tree::sT_label)

  push!(Ξ, sTbd(e(tree), false, true))
  if def1(tree)
    sTbd!(Ξ, tree.d2)
    sTbd!(Ξ, tree.d1)
  end
end




"""
    sTfbd!(Ξ::Vector{sTfbd}, tree::sTf_label)

Make edge tree `Ξ` from the recursive tree.
"""
function sTfbd!(Ξ::Vector{sTfbd}, tree::sTf_label)

  # no fossil can be a `true` tip nor extinct
  if istip(tree) && isfossil(tree)
    # add first daughter tree for tip fossil that is extinct with a 0 edge.
    push!(Ξ, sTfbd(
               sTfbd(0.0, true, false, true),
               e(tree), false, true, true))
  else
    push!(Ξ, sTfbd(e(tree), false, isfossil(tree), true))
  end

  if def2(tree) sTfbd!(Ξ, tree.d2) end
  if def1(tree) sTfbd!(Ξ, tree.d1) end
end




"""
    make_Ξ(idf::Vector{iBffs}, ::Type{sTbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, ::Type{sTbd})
  Ξ = sTbd[]
  for i in Base.OneTo(lastindex(idf))
    ξ = sTbd(e(idf[i]), false, true)
    push!(Ξ, ξ)
  end
  return Ξ
end




"""
    make_Ξ(idf::Vector{iBffs}, ::Type{sTfbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, ::Type{sTfbd})
  Ξ = sTfbd[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    ξ = sTfbd(e(bi), ie(bi), ifos(bi), true)
    push!(Ξ, ξ)
  end
  return Ξ
end




"""
    iTgbmpb!(Ξ   ::Vector{iTgbmpb},
             tree::sT_label,
             δt  ::Float64,
             srδt::Float64,
             lλa ::Float64,
             α   ::Float64,
             σλ  ::Float64)

Make edge tree `Ξ` from the recursive tree.
"""
function iTgbmpb!(Ξ   ::Vector{iTgbmpb},
                  tree::sT_label,
                  δt  ::Float64,
                  srδt::Float64,
                  lλa ::Float64,
                  α   ::Float64,
                  σλ  ::Float64)

  et = e(tree)

  if iszero(et)
    lλv  = Float64[lλa, lλa]
    fdti = 0.0
    l    = 2
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    l   = lastindex(lλv)
  end

  push!(Ξ, iTgbmpb(et, true, δt, fdti, lλv))
  if def1(tree)
    iTgbmpb!(Ξ, tree.d2, δt, srδt, lλv[l], α, σλ)
    iTgbmpb!(Ξ, tree.d1, δt, srδt, lλv[l], α, σλ)
  end
end




"""
    iTgbmce!(Ξ   ::Vector{iTgbmce},
             tree::sT_label,
             δt  ::Float64,
             srδt::Float64,
             lλa ::Float64,
             α   ::Float64,
             σλ  ::Float64)

Make edge tree `Ξ` from the recursive tree.
"""
function iTgbmce!(Ξ   ::Vector{iTgbmce},
                  tree::sT_label,
                  δt  ::Float64,
                  srδt::Float64,
                  lλa ::Float64,
                  α   ::Float64,
                  σλ  ::Float64)

  et = e(tree)

  if iszero(et)
    lλv  = Float64[lλa, lλa]
    fdti = 0.0
    l    = 2
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    l   = lastindex(lλv)
  end

  push!(Ξ, iTgbmce(et, δt, fdti, false, true, lλv))
  if def1(tree)
    iTgbmce!(Ξ, tree.d2, δt, srδt, lλv[l], α, σλ)
    iTgbmce!(Ξ, tree.d1, δt, srδt, lλv[l], α, σλ)
  end
end




"""
    iTgbmct!(Ξ   ::Vector{iTgbmct},
             tree::sT_label,
             δt  ::Float64,
             srδt::Float64,
             lλa ::Float64,
             α   ::Float64,
             σλ  ::Float64)

Make edge tree `Ξ` from the recursive tree.
"""
function iTgbmct!(Ξ   ::Vector{iTgbmct},
                  tree::sT_label,
                  δt  ::Float64,
                  srδt::Float64,
                  lλa ::Float64,
                  α   ::Float64,
                  σλ  ::Float64)

  et = e(tree)

  if iszero(et)
    lλv  = Float64[lλa, lλa]
    fdti = 0.0
    l    = 2
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    l   = lastindex(lλv)
  end

  push!(Ξ, iTgbmct(et, δt, fdti, false, true, lλv))
  if def1(tree)
    iTgbmct!(Ξ, tree.d2, δt, srδt, lλv[l], α, σλ)
    iTgbmct!(Ξ, tree.d1, δt, srδt, lλv[l], α, σλ)
  end
end




"""
    iTgbmbd!(Ξ   ::Vector{iTgbmbd},
             tree::sT_label,
             δt  ::Float64,
             srδt::Float64,
             lλa ::Float64,
             lμa ::Float64,
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64)

Make edge tree `Ξ` from the recursive tree.
"""
function iTgbmbd!(Ξ   ::Vector{iTgbmbd},
                  tree::sT_label,
                  δt  ::Float64,
                  srδt::Float64,
                  lλa ::Float64,
                  lμa ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64)

  et = e(tree)

  if iszero(et)
    lλv  = Float64[lλa, lλa]
    lμv  = Float64[lμa, lμa]
    fdti = 0.0
    l    = 2
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    lμv = sim_bm(lμa, α, σμ, δt, fdti, srδt, nt)
    l   = nt + 2
  end

  push!(Ξ, iTgbmbd(et, δt, fdti, false, true, lλv, lμv))
  if def1(tree)
    iTgbmbd!(Ξ, tree.d2, δt, srδt, lλv[l], lμv[l], α, σλ, σμ)
    iTgbmbd!(Ξ, tree.d1, δt, srδt, lλv[l], lμv[l], α, σλ, σμ)
  end
end




"""
    couple(Ξ::Vector{T},
           idf::Vector{iBffs},
           ix ::Int64) where {T <: iTree}

Build tree from decoupled tree.
"""
function couple(Ξ::Vector{T},
                idf::Vector{iBffs},
                ix ::Int64) where {T <: iTree}

  bi = idf[ix]
  ξi = Ξ[ix]
  if !it(bi)
    ξit = fixtip(ξi)
    if !iszero(d1(bi)) ξit.d1 = couple(Ξ, idf, d1(bi)) end
    if !iszero(d2(bi)) ξit.d2 = couple(Ξ, idf, d2(bi)) end
  end

  return ξi
end




"""
    treelength(Ξ::Vector{T}) where {T<: iTree}

Return the branch length sum of `Ξ`.
"""
function treelength(Ξ::Vector{T}) where {T<: iTree}
  L = 0.0
  for ξ in Ξ
    L += _treelength(ξ, 0.0)
  end
  return L
end




"""
    _ctl(Ξ::Vector{T}) where {T <: iTgbm}

Return the branch length sum of `Ξ` based on `δt` and `fδt`
for debugging purposes.
"""
function _ctl(Ξ::Vector{T}) where {T <: iTgbm}
  L = 0.0
  for ξ in Ξ
    L += _ctl(ξ, 0.0)
  end
  return L
end




"""
    nnodesinternal(Ξ::Vector{T}) where {T<: iTree}

Return the number of internal nodes in `Ξ`.
"""
function nnodesinternal(Ξ::Vector{T}) where {T<: iTree}
  n = 0
  for ξ in Ξ
    n += _nnodesinternal(ξ, 0)
  end
  n += Float64(lastindex(Ξ) - 1)/2.0

  return n
end




"""
    nnodesbifurcation(Ξ::Vector{T}) where {T<: iTree}

Return the number of bifurcating nodes in `Ξ`.
"""
function nnodesbifurcation(Ξ::Vector{T}) where {T<: iTree}
  ns = 0
  nf = 0
  for ξ in Ξ
    ns += _nnodesbifurcation(ξ, 0)
    nf += isfixfossil(ξ)
  end
  ns += Float64(lastindex(Ξ) - nf - 1)*0.5

  return n
end




"""
    ntipsextinct(Ξ::Vector{T}) where {T<: iTree}

Return the number of extinct nodes in `Ξ`.
"""
function ntipsextinct(Ξ::Vector{T}) where {T<: iTree}
  n = 0
  for ξ in Ξ
    n += _ntipsextinct(ξ, 0)
  end
  return n
end




"""
    nfossils(Ξ::Vector{T}) where {T<: iTree}

Return the number of fossil nodes in `Ξ`.
"""
function nfossils(Ξ::Vector{T}) where {T<: iTree}
  n = 0
  for ξ in Ξ
    n += _nfossils(ξ, 0)
  end
  return n
end




"""
    sss_gbm(Ξ::Vector{T}, α::Float64) where {T <: iTgbm}

Returns the standardized sum of squares a `iTgbm` according
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(Ξ::Vector{T}, α::Float64) where {T <: iTgbm}

  n   = 0.0
  ssλ = 0.0
  for ξi in Ξ
    ssλ, n = _sss_gbm(ξi, α, ssλ, n)
  end

  return ssλ, n
end




"""
    sss_gbm(Ξ::Vector{iTgbmbd}, α::Float64)

Returns the standardized sum of squares a `iTgbm` according
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(Ξ::Vector{iTgbmbd}, α::Float64)

  n   = 0.0
  ssλ = 0.0
  ssμ = 0.0
  for ξi in Ξ
    ssλ, ssμ, n = _sss_gbm(ξi, α, ssλ, ssμ, n)
  end

  return ssλ, ssμ, n
end




"""
    Σλ_gbm(Ξ::Vector{T}) where {T<: iTgbm}

Return the internal nodes of `Ξ`.
"""
function Σλ_gbm(Ξ::Vector{T}) where {T <: iTgbm}
  Σλ = 0.0
  for ξ in Ξ
    Σλ += Σλ_gbm(ξ)
  end
  return Σλ
end

