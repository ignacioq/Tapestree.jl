#=

decoupled tree utilities

Ignacio Quintero Mächler

t(-_-t)

Created 05 11 2020
=#




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

  push!(Ψ, iTgbmpb(et, δt, fdti, lλv))
  if isdefined(tree, :d1)
    iTgbmpb!(Ψ, tree.d2, δt, srδt, lλv[l], α, σλ) 
    iTgbmpb!(Ψ, tree.d1, δt, srδt, lλv[l], α, σλ)
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
  if !iszero(d1(bi))
    ψit = fixtip(ψi)
    ψit.d1 = couple(psi, idf, d1(bi))
    ψit.d2 = couple(psi, idf, d2(bi))
  end

  return ψi
end




"""
    treelength(psi::Vector{sTbd})

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
    nnodesinternal(psi::Vector{sTbd})

Return the internal nodes of `Ψ`.
"""
function nnodesinternal(psi::Vector{sTbd})
  n = 0
  for ψ in psi
    n += _nnodesinternal(ψ, 0)
  end
  n += Float64(lastindex(psi) - 1)/2.0

  return n
end




"""
    ntipsextinct(psi::Vector{sTbd})

Return the internal nodes of `Ψ`.
"""
function ntipsextinct(psi::Vector{sTbd})
  n = 0
  for ψ in psi
    n += _ntipsextinct(ψ, 0)
  end
  return n
end

