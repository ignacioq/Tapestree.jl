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
function make_Ψ(idf::Vector{iBffs})
  Ψ = sTbd[]
  for i in Base.OneTo(lastindex(idf))
    ψ = sTbd(e(idf[i]), false, true)
    push!(Ψ, ψ)
  end
  return Ψ
end




"""
    build(psi::Vector{T},
          idf::Vector{iBffs},
          ix ::Int64) where {T <: iTree}

Build tree from decoupled tree.
"""
function build(psi::Vector{T},
               idf::Vector{iBffs},
               ix ::Int64) where {T <: iTree}

  bi = idf[ix]
  ψi = psi[ix]
  if !iszero(d1(bi))
    ψit = fixtip(ψi)
    ψit.d1 = build(psi, idf, d1(bi))
    ψit.d2 = build(psi, idf, d2(bi))
  end

  return ψi
end




"""
    treelength(psi::Vector{sTbd})

Return the branch length sum of `Ψ`.
"""
function treelength(psi::Vector{sTbd})
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

