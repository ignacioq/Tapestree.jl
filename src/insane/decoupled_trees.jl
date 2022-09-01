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
    push!(Ξ, sTpbX(e(idfi), true, xii, xfi))
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
    push!(Ξ, sTbdX(e(idfi), false, true, xii, xfi))
  end

  return Ξ
end




"""
    make_Ξ(idf::Vector{iBffs}, xr::Vector{Float64}, ::Type{sTfbdX})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, xr::Vector{Float64}, ::Type{sTfbdX})
  Ξ = sTfbdX[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    paix = iszero(paix) ? 1 : paix
    xii  = xr[paix]
    xfi  = xr[i]
    iψ = isfossil(idfi)
    if iψ && it(idfi)
      push!(Ξ, sTfbdX(
                 sTfbdX(1.0e-10, true, false, false, xfi, xfi),
                 e(idfi), false, true, true, xii, xfi))
    else
      push!(Ξ, sTfbdX(e(idfi), false, iψ, true, xii, xfi))
    end
  end

  return Ξ
end




"""
    make_Ξ(idf ::Vector{iBffs},
           xr  ::Vector{Float64},
           lλa ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           σx  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTpbX})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                xr  ::Vector{Float64},
                lλa ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTpbX})

  lλi = lλa
  Ξ   = iTpbX[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    paix = iszero(paix) ? 1 : paix
    xii  = xr[paix]
    xfi  = xr[i]
    et   = e(idfi)
    if i > 1 
      lλi = λt(idf[paix])
    end

    if iszero(et)
      lλv = Float64[lλi, lλi]
      xv  = Float64[xii, xfi]
      fdti = 0.0
      l    = 2
    else
      nt, fdti = divrem(et, δt, RoundDown)
      nt = Int64(nt)

      if iszero(fdti)
        fdti = δt
        nt  -= 1
      end

      lλv = bm(lλi, α,   σλ, δt, fdti, srδt, nt)
      xv  = bb(xii, xfi, σx, δt, fdti, srδt, nt)
      l   = nt + 2
    end
    setλt!(idfi, lλv[l])
    push!(λst(idfi), lλv[l])
    push!(Ξ, iTpbX(et, true, δt, fdti, lλv, xv))
  end

  return Ξ
end





"""
     make_Ξ(idf ::Vector{iBffs},
            xr  ::Vector{Float64},
            lλa ::Float64,
            lμa ::Float64,
            α   ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            σx  ::Float64,
            δt  ::Float64,
            srδt::Float64,
            ::Type{iTbdX})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                xr  ::Vector{Float64},
                lλa ::Float64,
                lμa ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTbdX})

  lλi = lλa
  lμi = lμa
  Ξ   = iTbdX[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    paix = iszero(paix) ? 1 : paix
    xii  = xr[paix]
    xfi  = xr[i]
    et   = e(idfi)
    if i > 1 
      lλi = λt(idf[paix])
      lμi = μt(idf[paix])
    end

    if iszero(et)
      lλv  = Float64[lλi, lλi]
      lμv  = Float64[lμi, lμi]
      xv   = Float64[xii, xfi]
      fdti = 0.0
      l    = 2
    else
      nt, fdti = divrem(et, δt, RoundDown)
      nt = Int64(nt)

      if iszero(fdti)
        fdti = δt
        nt  -= 1
      end

      lλv = bm(lλi,   α, σλ, δt, fdti, srδt, nt)
      lμv = bm(lμi, 0.0, σμ, δt, fdti, srδt, nt)
      xv  = bb(xii, xfi, σx, δt, fdti, srδt, nt)
      l   = nt + 2
    end
    setλt!(idfi, lλv[l])
    setμt!(idfi, lμv[l])
    push!(λst(idfi), lλv[l])
    push!(μst(idfi), lμv[l])
    push!(Ξ, iTbdX(et, δt, fdti, false, true, lλv, lμv, xv))
  end

  return Ξ
end




"""
     make_Ξ(idf ::Vector{iBffs},
            xr  ::Vector{Float64},
            lλa ::Float64,
            lμa ::Float64,
            α   ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            σx  ::Float64,
            δt  ::Float64,
            srδt::Float64,
            ::Type{iTfbdX})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                xr  ::Vector{Float64},
                lλa ::Float64,
                lμa ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                σx  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTfbdX})

  lλi = lλa
  lμi = lμa
  Ξ   = iTfbdX[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    paix = iszero(paix) ? 1 : paix
    xii  = xr[paix]
    xfi  = xr[i]
    et   = e(idfi)
    if i > 1 
      lλi = λt(idf[paix])
      lμi = μt(idf[paix])
    end

    if iszero(et)
      lλv  = Float64[lλi, lλi]
      lμv  = Float64[lμi, lμi]
      xv   = Float64[xii, xfi]
      fdti = 0.0
      l    = 2
    else
      nt, fdti = divrem(et, δt, RoundDown)
      nt = Int64(nt)

      if iszero(fdti)
        fdti = δt
        nt  -= 1
      end

      lλv = bm(lλi,   α, σλ, δt, fdti, srδt, nt)
      lμv = bm(lμi, 0.0, σμ, δt, fdti, srδt, nt)
      xv  = bb(xii, xfi, σx, δt, fdti, srδt, nt)
      l   = nt + 2
    end

    if it(idfi) && isfossil(idfi)
      lλl = lλv[l]
      lμl = lμv[l]
      xvl = xv[l]
      push!(Ξ, iTfbdX(
                 iTfbdX(1.0e-10, δt, 1.0e-10, true, false, false, 
                   Float64[lλl, rnorm(lλl + α*1.0e-10, 1.0e-5*σλ)], 
                   Float64[lμl, rnorm(lμl,             1.0e-5*σμ)],
                   Float64[xvl, rnorm(xvl,             1.0e-5*σx)]),
                 et, δt, fdti, false, true, true, lλv, lμv, xv))
    else
      push!(Ξ, iTfbdX(et, δt, fdti, false, isfossil(idfi), true, lλv, lμv, xv))
    end

    setλt!(idfi, lλv[l])
    setμt!(idfi, lμv[l])
    push!(λst(idfi), lλv[l])
    push!(μst(idfi), lμv[l])
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
    make_Ξ(idf::Vector{iBffs}, ::Type{sTfbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, ::Type{sTfbd})
  Ξ = sTfbd[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    iψ = isfossil(idfi)
    if iψ && it(idfi)
      push!(Ξ, sTfbd(sTfbd(1.0e-10, true, false, false),
                     e(idfi), false, true, true))
    else
      push!(Ξ, sTfbd(e(idfi), false, iψ, true))
    end
  end

  return Ξ
end





"""
    make_Ξ(idf::Vector{iBffs}, ::Type{sTbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, ::Type{sTbd})
  Ξ = sTbd[]
  for i in Base.OneTo(lastindex(idf))
    push!(Ξ, sTbd(e(idf[i]), false, true))
  end
  return Ξ
end




"""
    make_Ξ(idf ::Vector{iBffs},
           lλa ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTpb})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                lλa ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTpb})

  lλi = lλa
  Ξ   = iTpb[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    et   = e(idfi)
    if i > 1 
      lλi = λt(idf[paix])
    end

    if iszero(et)
      lλv = Float64[lλi, lλi]
      fdti = 0.0
      l    = 2
    else
      nt, fdti = divrem(et, δt, RoundDown)
      nt = Int64(nt)

      if iszero(fdti)
        fdti = δt
        nt  -= 1
      end

      lλv = bm(lλi, α, σλ, δt, fdti, srδt, nt)
      l   = nt + 2
    end
    setλt!(idfi, lλv[l])
    push!(Ξ, iTpb(et, true, δt, fdti, lλv))
  end

  return Ξ
end




"""
    make_Ξ(idf ::Vector{iBffs},
           lλa ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTce})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                lλa ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTce})

  lλi = lλa
  Ξ   = iTce[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    et   = e(idfi)
    if i > 1 
      lλi = λt(idf[paix])
    end

    if iszero(et)
      lλv = Float64[lλi, lλi]
      fdti = 0.0
      l    = 2
    else
      nt, fdti = divrem(et, δt, RoundDown)
      nt = Int64(nt)

      if iszero(fdti)
        fdti = δt
        nt  -= 1
      end

      lλv = bm(lλi, α, σλ, δt, fdti, srδt, nt)
      l   = nt + 2
    end
    setλt!(idfi, lλv[l])
    push!(Ξ, iTce(et, δt, fdti, false, true, lλv))
  end

  return Ξ
end




"""
    make_Ξ(idf ::Vector{iBffs},
           lλa ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTct})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                lλa ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTct})

  lλi = lλa
  Ξ   = iTct[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    et   = e(idfi)
    if i > 1 
      lλi = λt(idf[paix])
    end

    if iszero(et)
      lλv = Float64[lλi, lλi]
      fdti = 0.0
      l    = 2
    else
      nt, fdti = divrem(et, δt, RoundDown)
      nt = Int64(nt)

      if iszero(fdti)
        fdti = δt
        nt  -= 1
      end

      lλv = bm(lλi, α, σλ, δt, fdti, srδt, nt)
      l   = nt + 2
    end
    setλt!(idfi, lλv[l])
    push!(Ξ, iTct(et, δt, fdti, false, true, lλv))
  end

  return Ξ
end





"""
    make_Ξ(idf ::Vector{iBffs},
           lλa ::Float64,
           lμa ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                λ   ::Float64,
                μ   ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTbd})

  Ξ   = iTbd[]
  _make_Ξ!(Ξ, 1, log(λ), log(μ), α, σλ, σμ, δt, srδt, idf, iTbd)

  return Ξ
end





"""
    _make_Ξ!(idf ::Vector{iBffs},
             i   ::Int64,
             lλ0 ::Float64,
             lμ0 ::Float64,
             α   ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64,
             ::Type{iTbd})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{iTbd},
                  i   ::Int64,
                  lλ0 ::Float64,
                  lμ0 ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  idf ::Vector{iBffs},
                  ::Type{iTbd})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)
  io = dn(bi)
  et = e(bi)

  if iszero(et)
    lλv  = [lλ0, lλ0]
    lμv  = [lμ0, lμ0]
    fdti = 0.0
    nts  = 0
  else
    ntF, fdti = divrem(et, δt, RoundDown)
    nts = Int64(ntF)

    if iszero(fdti) || isone(io)
      fdti  = δt
      nts  -= 1
    end

    lλv = bm(lλ0,   α, σλ, δt, fdti, srδt, nts)
    lμv = bm(lμ0, 0.0, σμ, δt, fdti, srδt, nts)
  end

  l = nts + 2

  setλt!(bi, lλv[l])
  push!(Ξ, iTbd(et, δt, fdti, false, true, lλv, lμv))

  if i1 > 0 
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλv[l], lμv[l], α, σλ, σμ, δt, srδt, idf, iTbd)
      _make_Ξ!(Ξ, i1, lλv[l], lμv[l], α, σλ, σμ, δt, srδt, idf, iTbd)
    else
      _make_Ξ!(Ξ, i1, lλv[l], lμv[l], α, σλ, σμ, δt, srδt, idf, iTbd)
    end
  end

  return nothing
end






"""
    make_Ξ(idf ::Vector{iBffs},
           lλa ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           tv  ::Vector{Float64},
           ev  ::Vector{Float64},
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTbd})

Make edge tree `Ξ` from the edge directory for fixed extinction.
"""
function make_Ξ(idf ::Vector{iBffs},
                lλa ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                tv  ::Vector{Float64},
                le  ::Vector{Float64},
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTbd})

  lλi = lλa
  Ξ   = iTbd[]
  ixiv = Int64[] # start point for fixed branches
  ixfv = Int64[] # end point for fixed branches
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    et   = e(idfi)
    if i > 1 
      lλi = λt(idf[paix])
      lμi = μt(idf[paix])
    end

    if iszero(et)
      lλv  = Float64[lλi, lλi]
      lμi  = linpred(ti(idfi), tv[1], tv[2], le[1], le[2])
      lμv  = Float64[lμi, lμi]
      fdti = 0.0
      l    = 2
      push!(ixiv, 1)
      push!(ixfv, 1)
    else
      nts, fdti = divrem(et, δt, RoundDown)
      nts = Int64(nts)

      if iszero(fdti)
        fdti = δt
        nts  -= 1
      end
      # speciation 
      lλv = bm(lλi, α, σλ, δt, fdti, srδt, nts)

      # extinction
      tii = ti(idfi)
      tif = tf(idfi)
      ix  = findfirst(x -> x < tii, tv) - 1
      push!(ixiv, ix)
      tc  = tii
      lμv = Float64[]
      push!(lμv, linpred(tii, tv[ix], tv[ix+1], le[ix], le[ix+1]))
      for i in Base.OneTo(nts)
        tc -= δt
        ix = findnext(x -> x < tc, tv, ix) - 1
        push!(lμv, linpred(tc, tv[ix], tv[ix+1], le[ix], le[ix+1]))
        ix += 1
      end
      ix = findnext(x -> x <= tif, tv, ix) - 1
      push!(lμv, linpred(tif, tv[ix], tv[ix+1], le[ix], le[ix+1]))
      l   = nts + 2
      push!(ixfv, ix)
    end

    setλt!(idfi, lλv[l])
    setμt!(idfi, lμv[l])
    push!(Ξ, iTbd(et, δt, fdti, false, true, lλv, lμv))
  end

  return Ξ, ixiv, ixfv
end




"""
    make_Ξ(idf ::Vector{iBffs},
           lλa ::Float64,
           lμa ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTfbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                lλa ::Float64,
                lμa ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTfbd})

  lλi = lλa
  lμi = lμa
  Ξ   = iTfbd[]
  for i in Base.OneTo(lastindex(idf))
    idfi = idf[i]
    paix = pa(idfi)
    et   = e(idfi)
    if i > 1 
      lλi = λt(idf[paix])
      lμi = μt(idf[paix])
    end

    if iszero(et)
      lλv  = Float64[lλi, lλi]
      lμv  = Float64[lμi, lμi]
      fdti = 0.0
      l    = 2
    else
      nt, fdti = divrem(et, δt, RoundDown)
      nt = Int64(nt)

      if iszero(fdti)
        fdti = δt
        nt  -= 1
      end

      lλv = bm(lλi, α,   σλ, δt, fdti, srδt, nt)
      lμv = bm(lμi, 0.0, σμ, δt, fdti, srδt, nt)
      l   = nt + 2
    end

    if it(idfi) && isfossil(idfi)
      lλl = lλv[l]
      lμl = lμv[l]
      push!(Ξ, iTfbd(
                 iTfbd(1.0e-10, δt, 1.0e-10, true, false, false, 
                   Float64[lλl, rnorm(lλl + α*1.0e-10, 1.0e-5*σλ)], 
                   Float64[lμl, rnorm(lμl,             1.0e-5*σμ)]),
                 et, δt, fdti, false, true, true, lλv, lμv))
    else
      push!(Ξ, iTfbd(et, δt, fdti, false, isfossil(idfi), true, lλv, lμv))
    end
    setλt!(idfi, lλv[l])
    setμt!(idfi, lμv[l])
  end

  return Ξ
end




"""
    couple(Ξ::Vector{T},
           idf::Vector{iBffs},
           ix ::Int64) where {T <: iTree}

Build tree from decoupled tree.
"""
function couple(Ξ  ::Vector{T},
                idf::Vector{iBffs},
                ix ::Int64) where {T <: iTree}

  bi  = idf[ix]
  ξi  = Ξ[ix]
  dni = dn(bi)
  if dni > 0
    ξit = fixtip(ξi)
    ξit.d1 = couple(Ξ, idf, d1(bi))
    if dni > 1 
      ξit.d2 = couple(Ξ, idf, d2(bi))
    end
  end

  return ξi
end




"""
    couple(Ξ  ::Vector{iTbd},
           idf::Vector{iBffs},
           ix ::Int64)

Build tree from decoupled tree.
"""
function couple(Ξ  ::Vector{iTbd},
                idf::Vector{iBffs},
                ix ::Int64)

  bi = idf[ix]
  ξi = iTbd(Ξ[ix])
  i1 = d1(bi)
  i2 = d2(bi)

  if i1 > 0
    ξit    = fixtip(ξi)
    if i2 > 0
      ξit.d1 = couple(Ξ, idf, i1)
      ξit.d2 = couple(Ξ, idf, i2)
    else
      ξd1 = couple(Ξ, idf, i1)
      sete!(ξit, e(ξit) + e(ξd1))
      setfdt!(ξit, fdt(ξd1))
      lλv = lλ(ξit)
      lμv = lμ(ξit)
      pop!(lλv)
      pop!(lμv)
      append!(lλv, lλ(ξd1))
      append!(lμv, lμ(ξd1))

      if def1(ξd1)
        ξit.d1 = ξd1.d1
        ξit.d2 = ξd1.d2
      end
    end
  end

  return ξi
end


T = couple(Ξ, idf, 1)





"""
    treelength(Ξ::Vector{T}) where {T <: iTree}

Return the branch length sum of `Ξ`.
"""
function treelength(Ξ::Vector{T}, ) where {T <: iTree}
  L = 0.0
  for ξ in Ξ
    L += _treelength(ξ, 0.0)
  end
  return L
end




"""
    treelength(Ξ  ::Vector{T},
               ets::Vector{Float64},
               bst::Vector{Float64},
               eix::Vector{Int64})  where {T <: iTf}

Return the branch length sum of `tree` at different epochs, initialized at `l`.
"""
function treelength(Ξ  ::Vector{T},
                    ets::Vector{Float64},
                    bst::Vector{Float64},
                    eix::Vector{Int64}) where {T <: iTf}

  nep = lastindex(ets) + 1
  ls  = zeros(nep)
  for i in Base.OneTo(lastindex(Ξ))
    _treelength!(Ξ[i], bst[i], ls, ets, eix[i], nep)
  end

  return ls
end




"""
    _ctl(Ξ::Vector{T}) where {T <: iT}

Return the branch length sum of `Ξ` based on `δt` and `fδt`
for debugging purposes.
"""
function _ctl(Ξ::Vector{T}) where {T <: iT}
  L = 0.0
  for ξ in Ξ
    L += _ctl(ξ, 0.0)
  end
  return L
end




"""
    nnodesinternal(Ξ::Vector{T}) where {T <: iTree}

Return the number of internal nodes in `Ξ`.
"""
function nnodesinternal(Ξ::Vector{T}) where {T <: iTree}
  n = 0
  for ξ in Ξ
    n += _nnodesinternal(ξ, 0)
  end
  n += Float64(lastindex(Ξ) - 1)/2.0

  return n
end




"""
    nnodesbifurcation(Ξ::Vector{T}) where {T <: iTree}

Return the number of bifurcating nodes in `Ξ`.
"""
function nnodesbifurcation(Ξ::Vector{T}) where {T <: iTf}
  ns = 0
  nf = 0
  for ξ in Ξ
    ns += _nnodesbifurcation(ξ, 0)
    nf += isinternalfossil(ξ)
  end
  ns += 0.5 * Float64(lastindex(Ξ) - nf - 1)

  return ns
end




"""
    ntipsextinct(Ξ::Vector{T}) where {T <: iTree}

Return the number of extinct nodes in `Ξ`.
"""
function ntipsextinct(Ξ::Vector{T}) where {T <: iTree}
  n = 0
  for ξ in Ξ
    n += _ntipsextinct(ξ, 0)
  end
  return n
end




"""
    nfossils(Ξ::Vector{T}) where {T <: iTree}

Return the number of fossil nodes in `Ξ`.
"""
function nfossils(Ξ::Vector{T}) where {T <: iTree}
  n = 0
  for ξ in Ξ
    n += _nfossils(ξ, 0)
  end
  return n
end




"""
    sss_gbm(Ξ::Vector{T}, α::Float64) where {T <: iT}

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(Ξ::Vector{T}, α::Float64) where {T <: iT}

  n   = 0.0
  ssλ = 0.0
  for ξi in Ξ
    ssλ, n = _sss_gbm(ξi, α, ssλ, n)
  end

  return ssλ, n
end




"""
    sss_gbm(Ξ::Vector{T}, α::Float64) where {T <: iTbdU}

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function sss_gbm(Ξ::Vector{T}, α::Float64) where {T <: iTbdU}

  n   = 0.0
  ssλ = 0.0
  ssμ = 0.0
  for ξi in Ξ
    ssλ, ssμ, n = _sss_gbm(ξi, α, ssλ, ssμ, n)
  end

  return ssλ, ssμ, n
end




"""
    sss_gbm(Ξ::Vector{T}, α::Float64) where {T <: iTbdU}

Returns the standardized sum of squares for a `iTX` according
to GBM lambda and X.
"""
function sss_gbm(Ξ::Vector{T}, α::Float64, βλ::Float64) where {T <: iTX}

  n   = 0.0
  ssλ = 0.0
  ssx = 0.0
  for ξi in Ξ
    ssλ, ssx, n = _sss_gbm(ξi, α, βλ, ssλ, ssx, n)
  end

  return ssλ, ssx, n
end




"""
    sss_gbm(Ξ::Vector{T}, α::Float64) where {T <: iTbdU}

Returns the standardized sum of squares for a `iTX` according
to GBM lambda and X.
"""
function sss_gbm(Ξ::Vector{T}, α::Float64, βλ::Float64) where {T <: iTbdUX}

  n   = 0.0
  ssλ = 0.0
  ssμ = 0.0
  ssx = 0.0
  for ξi in Ξ
    ssλ, ssμ, ssx, n = _sss_gbm(ξi, α, βλ, ssλ, ssμ, ssx, n)
  end

  return ssλ, ssμ, ssx, n
end




"""
    Σλ_gbm(Ξ::Vector{T}) where {T<: iT}

Return the sum over `λ` gbm.
"""
function Σλ_gbm(Ξ::Vector{T}) where {T <: iT}
  Σλ = 0.0
  for ξ in Ξ
    Σλ += Σλ_gbm(ξ)
  end
  return Σλ
end

