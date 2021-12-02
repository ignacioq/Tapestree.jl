#=

decoupled tree utilities

Ignacio Quintero Mächler

t(-_-t)

Created 05 11 2020
=#




"""
    sTbd!(Ψ::Vector{sTpb}, tree::sT_label)

Make edge tree `Ψ` from the edge directory.
"""
function sTpb!(Ψ::Vector{sTpb}, tree::sT_label)

  push!(Ψ, sTpb(e(tree), true))
  if isdefined(tree, :d1)
    sTpb!(Ψ, tree.d2)
    sTpb!(Ψ, tree.d1)
  end
end



"""
    sTbd!(Ψ::Vector{sTbd}, tree::sT_label)

Make edge tree `Ψ` from the edge directory.
"""
function sTbd!(Ψ::Vector{sTbd}, tree::sT_label)

  push!(Ψ, sTbd(e(tree), false, true))
  if isdefined(tree, :d1)
    sTbd!(Ψ, tree.d2)
    sTbd!(Ψ, tree.d1)
  end
end




"""
    make_Ψ(idf::Vector{iBffs})

Make edge tree `Ψ` from the edge directory.
"""
function make_Ψ(idf::Vector{iBffs}, ::Type{sTbd})
  Ψ = sTbd[]
  for i in Base.OneTo(lastindex(idf))
    ψ = sTbd(e(idf[i]), false, true)
    push!(Ψ, ψ)
  end
  return Ψ
end




"""
    iTgbmpb!(Ψ   ::Vector{iTgbmpb},
             tree::sT_label,
             δt  ::Float64, 
             srδt::Float64, 
             lλa ::Float64,
             α   ::Float64,
             σλ  ::Float64)

Make edge tree `Ψ` from the edge directory.
"""
function iTgbmpb!(Ψ   ::Vector{iTgbmpb},
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
    end
    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    l   = lastindex(lλv)
  end

  push!(Ψ, iTgbmpb(et, true, δt, fdti, lλv))
  if isdefined(tree, :d1)
    iTgbmpb!(Ψ, tree.d2, δt, srδt, lλv[l], α, σλ) 
    iTgbmpb!(Ψ, tree.d1, δt, srδt, lλv[l], α, σλ)
  end
end





"""
    iTgbmce!(Ψ   ::Vector{iTgbmce},
             tree::sT_label,
             δt  ::Float64, 
             srδt::Float64, 
             lλa ::Float64,
             α   ::Float64,
             σλ  ::Float64)

Make edge tree `Ψ` from the edge directory.
"""
function iTgbmce!(Ψ   ::Vector{iTgbmce},
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
    end
    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    l   = lastindex(lλv)
  end

  push!(Ψ, iTgbmce(et, δt, fdti, false, true, lλv))
  if isdefined(tree, :d1)
    iTgbmce!(Ψ, tree.d2, δt, srδt, lλv[l], α, σλ) 
    iTgbmce!(Ψ, tree.d1, δt, srδt, lλv[l], α, σλ)
  end
end




"""
    iTgbmct!(Ψ   ::Vector{iTgbmct},
             tree::sT_label,
             δt  ::Float64, 
             srδt::Float64, 
             lλa ::Float64,
             α   ::Float64,
             σλ  ::Float64)

Make edge tree `Ψ` from the edge directory.
"""
function iTgbmct!(Ψ   ::Vector{iTgbmct},
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
    end
    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    l   = lastindex(lλv)
  end

  push!(Ψ, iTgbmct(et, δt, fdti, false, true, lλv))
  if isdefined(tree, :d1)
    iTgbmct!(Ψ, tree.d2, δt, srδt, lλv[l], α, σλ) 
    iTgbmct!(Ψ, tree.d1, δt, srδt, lλv[l], α, σλ)
  end
end




"""
    iTgbmbd!(Ψ   ::Vector{iTgbmbd},
             tree::sT_label,
             δt  ::Float64, 
             srδt::Float64, 
             lλa ::Float64,
             lμa ::Float64,
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64)

Make edge tree `Ψ` from the edge directory.
"""
function iTgbmbd!(Ψ   ::Vector{iTgbmbd},
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
    end
    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    lμv = sim_bm(lμa, α, σμ, δt, fdti, srδt, nt)
    l   = nt + 2
  end

  push!(Ψ, iTgbmbd(et, δt, fdti, false, true, lλv, lμv))
  if isdefined(tree, :d1)
    iTgbmbd!(Ψ, tree.d2, δt, srδt, lλv[l], lμv[l], α, σλ, σμ)
    iTgbmbd!(Ψ, tree.d1, δt, srδt, lλv[l], lμv[l], α, σλ, σμ)
  end
end




"""
    couple(psi::Vector{T},
           idf::Vector{iBffs},
           ix ::Int64) where {T <: iTree}

Build tree from decoupled tree.
"""
function couple(psi::Vector{T},
                idf::Vector{iBffs},
                ix ::Int64) where {T <: iTree}

  bi = idf[ix]
  ψi = psi[ix]
  if !it(bi)
    ψit = fixtip(ψi)
    ψit.d1 = couple(psi, idf, d1(bi))
    ψit.d2 = couple(psi, idf, d2(bi))
  end

  return ψi
end




"""
    treelength(psi::Vector{T}) where {T<: iTree}

Return the branch length sum of `Ψ`.
"""
function treelength(psi::Vector{T}) where {T<: iTree}
  L = 0.0
  for ψ in psi
    L += _treelength(ψ, 0.0)
  end
  return L
end





"""
    _ctl(tree::Vector{T}) where {T <: iTgbm}

Return the branch length sum of `tree` based on `δt` and `fδt` 
for debugging purposes.
"""
function _ctl(tree::Vector{T}) where {T <: iTgbm}

  L = 0.0
  for ψ in psi
    L += _ctl(ψ, 0.0)
  end
  return L
end




"""
    nnodesinternal(psi::Vector{T}) where {T<: iTree}

Return the internal nodes of `Ψ`.
"""
function nnodesinternal(psi::Vector{T}) where {T<: iTree}
  n = 0
  for ψ in psi
    n += _nnodesinternal(ψ, 0)
  end
  n += Float64(lastindex(psi) - 1)/2.0

  return n
end




"""
    ntipsextinct(psi::Vector{T}) where {T<: iTree}

Return the internal nodes of `Ψ`.
"""
function ntipsextinct(psi::Vector{T}) where {T<: iTree}
  n = 0
  for ψ in psi
    n += _ntipsextinct(ψ, 0)
  end
  return n
end





"""
    sss_gbm(psi::Vector{T}, α::Float64) where {T <: iTgbm}

Returns the standardized sum of squares a `iTgbm` according 
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(psi::Vector{T}, α::Float64) where {T <: iTgbm}

  n   = 0.0
  ssλ = 0.0
  for ψi in psi
    ssλ, n = _sss_gbm(ψi, α, ssλ, n)
  end

  return ssλ, n
end





"""
    sss_gbm(psi::Vector{iTgbmbd}, α::Float64)

Returns the standardized sum of squares a `iTgbm` according 
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(psi::Vector{iTgbmbd}, α::Float64)

  n   = 0.0
  ssλ = 0.0
  ssμ = 0.0
  for ψi in psi
    ssλ, ssμ, n = _sss_gbm(ψi, α, ssλ, ssμ, n)
  end

  return ssλ, ssμ, n
end





"""
    Σλ_gbm(psi::Vector{T}) where {T<: iTgbm}

Return the internal nodes of `Ψ`.
"""
function Σλ_gbm(psi::Vector{T}) where {T <: iTgbm}
  Σλ = 0.0
  for ψ in psi
    Σλ += Σλ_gbm(ψ)
  end
  return Σλ
end

