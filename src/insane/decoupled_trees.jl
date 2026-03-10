#=

decoupled tree utilities

Ignacio Quintero Mächler

t(-_-t)

Created 05 11 2020
=#




"""
    make_Ξ(idf::Vector{iBffs}, ::Type{sTb})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, ::Type{sTb})
  Ξ = sTb[]
  for bi in idf
    push!(Ξ, sTb(e(bi), true))
  end
  return Ξ
end




"""
    make_Ξ(idf::Vector{iBffs}, ::Type{sTbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, ::Type{sTbd})
  Ξ = sTbd[]
  for bi in idf
    push!(Ξ, sTbd(e(bi), false, true))
  end
  return Ξ
end




"""
    make_Ξ(idf::Vector{iBffs}, ::Type{sTfbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, ::Type{sTfbd})
  Ξ = sTfbd[]
  for bi in idf
    if isfossil(bi) && iszero(d1(bi))
      push!(Ξ, sTfbd(sTfbd(0.0, true, false, false),
                     e(bi), false, true, true))
    else
      push!(Ξ, sTfbd(e(bi), false, isfossil(bi), true))
    end
  end

  return Ξ
end




"""
    make_Ξ(idf ::Vector{iBffs},
           λ   ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{T}) where {T <: iT}

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                λ   ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{T}) where {T <: iT}

  Ξ = T[]
  _make_Ξ!(Ξ, 1, log(λ), α, σλ, δt, srδt, idf)

  return Ξ
end




"""
    _make_Ξ!(Ξ   ::Vector{T},
             i   ::Int64,
             lλ0 ::Float64,
             α   ::Float64,
             σλ  ::Float64,
             δt  ::Float64,
             srδt::Float64,
             idf ::Vector{iBffs},
             ::Type{T}) where {T <: iT}

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{T},
                  i   ::Int64,
                  lλ0 ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  idf ::Vector{iBffs}) where {T <: iT}

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)
  et = e(bi)

  if iszero(et)
    lλv  = [lλ0, lλ0]
    fdti = 0.0
    nts  = 0
  else
    ntF, fdti = divrem(et, δt, RoundDown)

    if isapprox(fdti, δt)
      ntF += 1.0
      fdti = δt
    end

    nts = Int64(ntF)

    if iszero(fdti) || (i1 > 0 && iszero(i2)) 
      fdti  = δt
      nts  -= 1
    end

    lλv = bm(lλ0, α, σλ, δt, fdti, srδt, nts)
  end

  l = nts + 2

  setλt!(bi, lλv[l])
  if T === iTb
    push!(Ξ, iTb(et, δt, fdti, true, lλv))
  else
    push!(Ξ, T(et, δt, fdti, false, true, lλv))
  end

  if i1 > 0 
    _make_Ξ!(Ξ, i1, lλv[l], α, σλ, δt, srδt, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλv[l], α, σλ, δt, srδt, idf)
    end
  end

  return nothing
end



"""
    make_Ξ(idf::Vector{iBffs}, λ::Float64, ::Type{T}) where {T <: cT}

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs}, λ::Float64, ::Type{T}) where {T <: cT}

  Ξ = T[]
  _make_Ξ!(Ξ, 1, log(λ), idf)

  return Ξ
end




"""
     _make_Ξ!(Ξ   ::Vector{cTb},
              i   ::Int64,
              lλ0 ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              idf ::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{cTb},
                  i   ::Int64,
                  lλ0 ::Float64,
                  idf ::Vector{iBffs})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)

  setλt!(bi, lλ0)
  push!(Ξ, cTb(e(bi), true, lλ0))

  if i1 > 0 
    _make_Ξ!(Ξ, i1, lλ0, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλ0, idf)
    end
  end

  return nothing
end




"""
     _make_Ξ!(Ξ   ::Vector{cTce},
              i   ::Int64,
              lλ0 ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              idf ::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{cTce},
                  i   ::Int64,
                  lλ0 ::Float64,
                  idf ::Vector{iBffs})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)

  setλt!(bi, lλ0)
  push!(Ξ, cTce(e(bi), false, true, lλ0))

  if i1 > 0 
    _make_Ξ!(Ξ, i1, lλ0, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλ0, idf)
    end
  end

  return nothing
end




"""
     _make_Ξ!(Ξ   ::Vector{cTct},
              i   ::Int64,
              lλ0 ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              idf ::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{cTct},
                  i   ::Int64,
                  lλ0 ::Float64,
                  idf ::Vector{iBffs})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)

  setλt!(bi, lλ0)
  push!(Ξ, cTct(e(bi), false, true, lλ0))

  if i1 > 0 
    _make_Ξ!(Ξ, i1, lλ0, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλ0, idf)
    end
  end

  return nothing
end




"""
    make_Ξ(idf ::Vector{iBffs},
                λ   ::Float64,
                μ   ::Float64,
                ::Type{T}) where {T <: Union{cTbd, cTfbd, acTfbd}}

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                λ   ::Float64,
                μ   ::Float64,
                ::Type{T}) where {T <: Union{cTbd, cTfbd, acTfbd}}
  Ξ = T[]
  _make_Ξ!(Ξ, 1, log(λ), log(μ), idf)

  return Ξ
end




"""
    _make_Ξ!(Ξ  ::Vector{cTbd},
             i  ::Int64,
             lλ0::Float64,
             lμ0::Float64,
             idf::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ  ::Vector{cTbd},
                  i  ::Int64,
                  lλ0::Float64,
                  lμ0::Float64,
                  idf::Vector{iBffs})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)

  setλt!(bi, lλ0)
  setμt!(bi, lμ0)
  push!(Ξ, cTbd(e(bi), false, true, lλ0, lμ0))

  if i1 > 0 
    _make_Ξ!(Ξ, i1, lλ0, lμ0, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλ0, lμ0, idf)
    end
  end

  return nothing
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

  Ξ = iTbd[]
  _make_Ξ!(Ξ, 1, log(λ), log(μ), α, σλ, σμ, δt, srδt, idf)

  return Ξ
end




"""
    _make_Ξ!(Ξ   ::Vector{iTbd},
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
                  idf ::Vector{iBffs})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)
  et = e(bi)

  if iszero(et)
    lλv  = [lλ0, lλ0]
    lμv  = [lμ0, lμ0]
    fdti = 0.0
    nts  = 0
  else
    ntF, fdti = divrem(et, δt, RoundDown)

    if isapprox(fdti, δt)
      ntF += 1.0
      fdti = δt
    end

    nts = Int64(ntF)

    if iszero(fdti) || (i1 > 0 && iszero(i2)) 
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
    _make_Ξ!(Ξ, i1, lλv[l], lμv[l], α, σλ, σμ, δt, srδt, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλv[l], lμv[l], α, σλ, σμ, δt, srδt, idf)
    end
  end

  return nothing
end




"""
    make_Ξ(idf ::Vector{iBffs},
           λ   ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           tv  ::Vector{Float64},
           le  ::Vector{Float64},
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                λ   ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                tv  ::Vector{Vector{Float64}},
                le  ::Vector{Vector{Float64}},
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTbd})

  Ξ    = iTbd[]
  ixiv = Int64[] # start point for fixed branches
  ixfv = Int64[] # end point for fixed branches
  _make_Ξ!(Ξ, ixiv, ixfv, 1, log(λ), α, σλ, tv, le, δt, srδt, idf)

  return Ξ, ixiv, ixfv
end




"""
    _make_Ξ!(Ξ   ::Vector{iTbd},
             ixiv::Vector{Int64},
             ixfv::Vector{Int64},
             i   ::Int64,
             lλ0 ::Float64,
             α   ::Float64,
             σλ  ::Float64,
             tv  ::Vector{Float64},
             le  ::Vector{Float64},
             δt  ::Float64,
             srδt::Float64,
             idf ::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{iTbd},
                  ixiv::Vector{Int64},
                  ixfv::Vector{Int64},
                  i   ::Int64,
                  lλ0 ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  tv  ::Vector{Vector{Float64}},
                  le  ::Vector{Vector{Float64}},
                  δt  ::Float64,
                  srδt::Float64,
                  idf ::Vector{iBffs})

  bi  = idf[i]
  tvi = tv[i]
  lei = le[i]
  i1  = d1(bi)
  i2  = d2(bi)
  et  = e(bi)

  if iszero(et)
    lλv  = [lλ0, lλ0]
    tii  = ti(bi)
    ix   = findfirst(x -> x < tii, tvi) - 1
    lμi  = linpred(tii, tvi[ix], tvi[ix+1], lei[ix], lei[ix+1])
    lμv  = Float64[lμi, lμi]
    fdti = 0.0
    nts  = 0
    push!(ixiv, ix)
    ix = findnext(x -> x <= tf(bi), tvi, ix) - 1
    push!(ixfv, ix)
  else

    ntF, fdti = divrem(et, δt, RoundDown)

    if isapprox(fdti, δt)
      ntF += 1.0
      fdti = δt
    end

    nts = Int64(ntF)

    if iszero(fdti) || (i1 > 0 && iszero(i2)) 
      fdti  = δt
      nts  -= 1
    end

    # speciation
    lλv = bm(lλ0, α, σλ, δt, fdti, srδt, nts)

    # extinction
    tii = ti(bi)
    tif = tf(bi)
    ix  = findfirst(x -> x < tii, tvi) - 1
    push!(ixiv, ix)
    tc  = tii
    lμv = Float64[]
    push!(lμv, linpred(tii, tvi[ix], tvi[ix+1], lei[ix], lei[ix+1]))
    for i in Base.OneTo(nts)
      tc -= δt
      ix = findnext(x -> x < tc, tvi, ix) - 1
      push!(lμv, linpred(tc, tvi[ix], tvi[ix+1], lei[ix], lei[ix+1]))
      ix += 1
    end
    ix = findnext(x -> x <= tif, tvi, ix) - 1
    push!(lμv, linpred(tif, tvi[ix], tvi[ix+1], lei[ix], lei[ix+1]))
    push!(ixfv, ix)
  end

  l = nts + 2

  setλt!(bi, lλv[l])
  push!(Ξ, iTbd(et, δt, fdti, false, true, lλv, lμv))

  if i1 > 0
    _make_Ξ!(Ξ, ixiv, ixfv, i1, lλv[l], α, σλ, tv, le, δt, srδt, idf)
    if i2 > 0
      _make_Ξ!(Ξ, ixiv, ixfv, i2, lλv[l], α, σλ, tv, le, δt, srδt, idf)
    end
  end

  return nothing
end




"""
    make_Ξ(idf ::Vector{iBffs},
           λ   ::Float64,
           μ   ::Float64,
           αλ  ::Float64,
           αμ  ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ::Type{cTfbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                λ   ::Float64,
                μ   ::Float64,
                ::Type{T}) where {T <: Union{cTfbd, acTfbd}}
  Ξ = T[]
  _make_Ξ!(Ξ, 1, log(λ), log(μ), idf)

  return Ξ
end




"""
    _make_Ξ!(Ξ   ::Vector{cTfbd},
             i   ::Int64,
             lλ0 ::Float64,
             lμ0 ::Float64,
             αλ  ::Float64,
             αμ  ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64,
             idf ::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{cTfbd},
                  i   ::Int64,
                  lλ0 ::Float64,
                  lμ0 ::Float64,
                  idf ::Vector{iBffs})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)

  setλt!(bi, lλ0)
  setμt!(bi, lμ0)


  if isfossil(bi) && iszero(i1)
    push!(Ξ, cTfbd(
               cTfbd(0.0, true, false, false, lλ0, lμ0),
               e(bi), false, true, true, lλ0, lμ0))
  else
    push!(Ξ, cTfbd(e(bi), false, isfossil(bi), true, lλ0, lμ0))
  end

  if i1 > 0 
    _make_Ξ!(Ξ, i1, lλ0, lμ0, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλ0, lμ0, idf)
    end
  end

  return nothing
end




"""
    _make_Ξ!(Ξ   ::Vector{acTfbd},
             i   ::Int64,
             lλ0 ::Float64,
             lμ0 ::Float64,
             αλ  ::Float64,
             αμ  ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64,
             idf ::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{acTfbd},
                  i   ::Int64,
                  lλ0 ::Float64,
                  lμ0 ::Float64,
                  idf ::Vector{iBffs})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)

  setλt!(bi, lλ0)
  setμt!(bi, lμ0)

  if isfossil(bi) && iszero(i1)
    push!(Ξ, acTfbd(
               acTfbd(0.0, true, false, false, false, lλ0, lμ0),
               e(bi), false, true, false, true, lλ0, lμ0))
  else
    push!(Ξ, acTfbd(e(bi), false, isfossil(bi), rand(Bool), true, lλ0, lμ0))
  end

  if i1 > 0 
    _make_Ξ!(Ξ, i1, lλ0, lμ0, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλ0, lμ0, idf)
    end
  end

  return nothing
end




"""
    make_Ξ(idf ::Vector{iBffs},
           λ   ::Float64,
           μ   ::Float64,
           αλ  ::Float64,
           αμ  ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{iTfbd})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                λ   ::Float64,
                μ   ::Float64,
                αλ  ::Float64,
                αμ  ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{iTfbd})

  Ξ = iTfbd[]
  _make_Ξ!(Ξ, 1, log(λ), log(μ), αλ, αμ, σλ, σμ, δt, srδt, idf)

  return Ξ
end




"""
    _make_Ξ!(Ξ   ::Vector{iTfbd},
             i   ::Int64,
             lλ0 ::Float64,
             lμ0 ::Float64,
             αλ  ::Float64,
             αμ  ::Float64,
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64,
             idf ::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{iTfbd},
                  i   ::Int64,
                  lλ0 ::Float64,
                  lμ0 ::Float64,
                  αλ  ::Float64,
                  αμ  ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  idf ::Vector{iBffs})

  bi = idf[i]
  i1 = d1(bi)
  i2 = d2(bi)
  et = e(bi)

  if iszero(et)
    lλv  = [lλ0, lλ0]
    lμv  = [lμ0, lμ0]
    fdti = 0.0
    nts  = 0
  else
    ntF, fdti = divrem(et, δt, RoundDown)

    if isapprox(fdti, δt)
      ntF += 1.0
      fdti = δt
    end

    nts = Int64(ntF)

    if iszero(fdti) || (i1 > 0 && iszero(i2) && !isfossil(bi))
      fdti  = δt
      nts  -= 1
    end

    lλv = bm(lλ0, αλ, σλ, δt, fdti, srδt, nts)
    lμv = bm(lμ0, αμ, σμ, δt, fdti, srδt, nts)
  end

  l = nts + 2

  if isfossil(bi) && iszero(i1)
    lλl = lλv[l]
    lμl = lμv[l]

    push!(Ξ, iTfbd(
               iTfbd(0.0, δt, 0.0, true, false, false, 
                 [lλl, rnorm(lλl + αλ*0.0, 0.0*σλ)], 
                 [lμl, rnorm(lμl + αμ*0.0, 0.0*σμ)]),
               et, δt, fdti, false, true, true, lλv, lμv))
  else
    push!(Ξ, iTfbd(et, δt, fdti, false, isfossil(bi), true, lλv, lμv))
  end
  setλt!(bi, lλv[l])

  if i1 > 0 
    _make_Ξ!(Ξ, i1, lλv[l], lμv[l], αλ, αμ, σλ, σμ, δt, srδt, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, lλv[l], lμv[l], αλ, αμ, σλ, σμ, δt, srδt, idf)
    end
  end

  return nothing
end




"""
    make_Ξ(idf ::Vector{iBffs},
           xr  ::Vector{Float64},
           σi  ::Float64,
           γi  ::Float64,
           δt  ::Float64,
           srδt::Float64,
           ::Type{sTxs})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf ::Vector{iBffs},
                xr  ::Vector{Float64},
                σi  ::Float64,
                γi  ::Float64,
                δt  ::Float64,
                srδt::Float64,
                ::Type{sTxs})

  Ξ = sTxs[]
  _make_Ξ!(Ξ, 1, xr, σi, γi, δt, srδt, idf)

  return Ξ
end




"""
    _make_Ξ!(Ξ   ::Vector{sTxs},
             i   ::Int64,
             xr  ::Vector{Float64},
             σi  ::Float64,
             γi  ::Float64,
             δt  ::Float64,
             srδt::Float64,
             idf ::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ   ::Vector{sTxs},
                  i   ::Int64,
                  xr  ::Vector{Float64},
                  σi  ::Float64,
                  γi  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  idf ::Vector{iBffs})

  bi  = idf[i]
  i1  = d1(bi)
  i2  = d2(bi)
  ip  = pa(bi)
  ip  = iszero(ip) ? 1 : ip
  et  = e(bi)
  xii = xr[ip]
  xfi = xr[i]

  if iszero(et)
    xv   = Float64[xii, xii]
    lσ2  = Float64[σi, σi]
    fdti = 0.0
  else
    ntF, fdti = divrem(et, δt, RoundDown)

    if isapprox(fdti, δt)
      ntF += 1.0
      fdti = δt
    end

    if iszero(fdti) || (i1 > 0 && iszero(i2) && !isfossil(bi))
      fdti  = δt
      ntF  -= 1.0
    end

    xv, lσ2 = dbb(xii, xfi, σi, σi, γi, δt, fdti, srδt, Int64(ntF))
  end

  push!(Ξ, sTxs(et, δt, fdti, xv, lσ2))

  if i1 > 0 
    _make_Ξ!(Ξ, i1, xr, σi, γi, δt, srδt, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, xr, σi, γi, δt, srδt, idf)
    end
  end

  return nothing
end






"""
    make_Ξ(idf::Vector{iBffs},
           xr ::Vector{Float64},
           σk ::Float64,
           ::Type{sTpe})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs},
                xr ::Vector{Float64},
                σk::Float64,
                ::Type{sTpe})

  Ξ = sTpe[]
  _make_Ξ!(Ξ, 1, false, xr, σk, idf)

  return Ξ
end




"""
    _make_Ξ!(Ξ  ::Vector{sTpe},
             i  ::Int64,
             xr ::Vector{Float64},
             σk::Float64,
             idf::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ  ::Vector{sTpe},
                  i  ::Int64,
                  clb::Bool,
                  xr ::Vector{Float64},
                  σk ::Float64,
                  idf::Vector{iBffs})

  bi  = idf[i]
  i1  = d1(bi)
  i2  = d2(bi)
  ip  = pa(bi)
  ip  = iszero(ip) ? 1 : ip
  et  = e(bi)
  xii = xr[ip]
  xfi = xr[i]

  if clb # if it is cladogenetic bud
    dx   = xfi - xii
    xii += sign(dx)*σk
  end

  shi = rand(Bool)
  push!(Ξ, sTpe(et, false, xii, xfi, shi, true))

  if i1 > 0 
    if i2 > 0 
      _make_Ξ!(Ξ, i1, shi, xr, σk, idf)
      _make_Ξ!(Ξ, i2, !shi, xr, σk, idf)
    else
      _make_Ξ!(Ξ, i1, false, xr, σk, idf)
    end

  end

  return nothing
end




"""
    make_Ξ(idf::Vector{iBffs},
           xr ::Vector{Float64},
           σai::Float64,
           σki::Float64,
           ::Type{sTfpe})

Make edge tree `Ξ` from the edge directory.
"""
function make_Ξ(idf::Vector{iBffs},
                xr ::Vector{Float64},
                σai::Float64,
                σki::Float64,
                ::Type{sTfpe})

  Ξ = sTfpe[]
  _make_Ξ!(Ξ, 1, xr, σai, σki, idf)

  return Ξ
end




"""
    _make_Ξ!(Ξ  ::Vector{sTfpe},
             i  ::Int64,
             xr ::Vector{Float64},
             σai::Float64,
             σki::Float64,
             idf::Vector{iBffs})

Make edge tree `Ξ` from the edge directory.
"""
function _make_Ξ!(Ξ  ::Vector{sTfpe},
                  i  ::Int64,
                  xr ::Vector{Float64},
                  σai::Float64,
                  σki::Float64,
                  idf::Vector{iBffs})

  bi  = idf[i]
  i1  = d1(bi)
  i2  = d2(bi)
  ip  = pa(bi)
  ip  = iszero(ip) ? 1 : ip
  et  = e(bi)
  xii = xr[ip]
  xfi = xr[i]
  shi = rand(Bool)

  if isfossil(bi) && iszero(d1(bi))
    push!(Ξ, sTfpe(sTfpe(1e-10, true, false, xfi, xfi, false, false),
                   e(bi), false, true, xii, xfi, shi, true))
  else
    push!(Ξ, sTfpe(e(bi), false, isfossil(bi), xii, xfi, shi, true))
  end

  if i1 > 0 
    _make_Ξ!(Ξ, i1, xr, σai, σki, idf)
    if i2 > 0 
      _make_Ξ!(Ξ, i2, xr, σai, σki, idf)
    end
  end

  return nothing
end




"""
    couple(Ξ::Vector{T},
           idf::Vector{iBffs},
           ix ::Int64) where {T <: sT}

Build tree from decoupled tree.
"""
function couple(Ξ  ::Vector{T},
                idf::Vector{iBffs},
                ix ::Int64) where {T <: sT}

  bi  = idf[ix]
  ξi  = T(Ξ[ix])
  i1  = d1(bi)
  i2  = d2(bi)

  if i1 > 0
    ξit = fixtip(ξi)
    if i2 > 0 
      ξit.d1 = couple(Ξ, idf, i1)
      ξit.d2 = couple(Ξ, idf, i2)
    elseif isfossil(bi)
      ξit.d1 = couple(Ξ, idf, i1)
    else
      ξd1 = couple(Ξ, idf, i1)
      adde!(ξit, e(ξd1))
      if isfossil(ξd1)
        fossilize!(ξit)
      end
      if def1(ξd1)
        ξit.d1 = ξd1.d1
        if def2(ξd1)
          ξit.d2 = ξd1.d2
        end
      end
    end
  end

  return ξi
end




"""
    couple(Ξ::Vector{T},
           idf::Vector{iBffs},
           ix ::Int64) where {T <: sT}

Build tree from decoupled tree.
"""
function couple(Ξ::Vector{T}, idf::Vector{iBffs}, ix::Int64) where {T <: Tpe}

  bi  = idf[ix]
  ξi  = T(Ξ[ix])
  i1  = d1(bi)
  i2  = d2(bi)

  if i1 > 0
    ξit = fixtip(ξi)
    if i2 > 0 
      ξit.d1 = couple(Ξ, idf, i1)
      ξit.d2 = couple(Ξ, idf, i2)
    elseif isfossil(bi)
      ξit.d1 = couple(Ξ, idf, i1)
    else
      ξd1 = couple(Ξ, idf, i1)

      adde!(ξit, e(ξd1))
      setxf!(ξit, xf(ξd1))
      setsh!(ξit, sh(ξd1))
      if isfossil(ξd1)
        fossilize!(ξit)
      end
      if def1(ξd1)
        ξit.d1 = ξd1.d1
        if def2(ξd1)
          ξit.d2 = ξd1.d2
        end
      end
    end
  end

  return ξi
end




"""
    couple(Ξ  ::Vector{T},
           idf::Vector{iBffs},
           ix ::Int64) where {T <: cT}

Build tree from decoupled tree.
"""
function couple(Ξ  ::Vector{T},
                idf::Vector{iBffs},
                ix ::Int64) where {T <: cT}

  bi  = idf[ix]
  ξi  = T(Ξ[ix])
  i1  = d1(bi)
  i2  = d2(bi)

  if i1 > 0
    ξit = fixtip(ξi)
    if i2 > 0 
      ξit.d1 = couple(Ξ, idf, i1)
      ξit.d2 = couple(Ξ, idf, i2)
    elseif isfossil(bi)
      ξit.d1 = couple(Ξ, idf, i1)
    else
      ξd1 = couple(Ξ, idf, i1)
      if isfossil(ξd1)
        fossilize!(ξit)
      end
      adde!(ξit, e(ξd1))
      if def1(ξd1)
        ξit.d1 = ξd1.d1
        if def2(ξd1)
          ξit.d2 = ξd1.d2
        end
      end
    end
  end

  return ξi
end




"""
    couple(Ξ::Vector{T},
           idf::Vector{iBffs},
           ix ::Int64) where {T <: iTree}

Build tree from decoupled tree.
"""
function couple(Ξ  ::Vector{T},
                idf::Vector{iBffs},
                ix ::Int64) where {T <: iT}

  bi  = idf[ix]
  ξi  = T(Ξ[ix])
  i1  = d1(bi)
  i2  = d2(bi)

  if i1 > 0
    ξit = fixtip(ξi)
    if i2 > 0 
      ξit.d1 = couple(Ξ, idf, i1)
      ξit.d2 = couple(Ξ, idf, i2)
    else
      ξd1 = couple(Ξ, idf, i1)
      lλv = lλ(ξit)
      if iszero(e(ξit))
        empty!(lλv)
      else
        pop!(lλv)
      end
      append!(lλv, lλ(ξd1))

      adde!(ξit, e(ξd1))
      setfdt!(ξit, fdt(ξd1))

      if def1(ξd1)
        ξit.d1 = ξd1.d1
        ξit.d2 = ξd1.d2
      end
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
function couple(Ξ  ::Vector{T},
                idf::Vector{iBffs},
                ix ::Int64) where {T <: iTbdU}

  bi = idf[ix]
  ξi = T(Ξ[ix])
  i1 = d1(bi)
  i2 = d2(bi)

  if i1 > 0
    ξit = fixtip(ξi)
    if i2 > 0
      ξit.d1 = couple(Ξ, idf, i1)
      ξit.d2 = couple(Ξ, idf, i2)
    elseif isfossil(bi)
      ξit.d1 = couple(Ξ, idf, i1)
    else
      ξd1 = couple(Ξ, idf, i1)
      if isfossil(ξd1)
        fossilize!(ξit)
      end
      lλv = lλ(ξit)
      lμv = lμ(ξit)
      if iszero(e(ξit))
        empty!(lλv)
        empty!(lμv) 
      else
        pop!(lλv)
        pop!(lμv)
      end
      append!(lλv, lλ(ξd1))
      append!(lμv, lμ(ξd1))

      adde!(ξit, e(ξd1))
      setfdt!(ξit, fdt(ξd1))
      if def1(ξd1)
        ξit.d1 = ξd1.d1
        if def2(ξd1)
          ξit.d2 = ξd1.d2
        end
      end
    end
  end

  return ξi
end




"""
    couple(Ξ  ::Vector{sTxs},
           idf::Vector{iBffs},
           ix ::Int64)

Build tree from decoupled tree.
"""
function couple(Ξ  ::Vector{sTxs},
                idf::Vector{iBffs},
                ix ::Int64)

  bi = idf[ix]
  ξi = sTxs(Ξ[ix])
  i1 = d1(bi)
  i2 = d2(bi)

  if i1 > 0
    ξit = fixtip(ξi)
    if i2 > 0
      ξit.d1 = couple(Ξ, idf, i1)
      ξit.d2 = couple(Ξ, idf, i2)
    elseif isfossil(bi)
      ξit.d1 = couple(Ξ, idf, i1)
    else
      ξd1 = couple(Ξ, idf, i1)
      xvi = xv(ξit)
      lσ2i = lσ2(ξit)
      if iszero(e(ξit))
        empty!(xvi)
        empty!(lσ2i) 
      else
        pop!(xvi)
        pop!(lσ2i)
      end
      append!(xvi, xv(ξd1))
      append!(lσ2i, lσ2(ξd1))

      adde!(ξit, e(ξd1))
      setfdt!(ξit, fdt(ξd1))
      if def1(ξd1)
        ξit.d1 = ξd1.d1
        if def2(ξd1)
          ξit.d2 = ξd1.d2
        end
      end
    end
  end

  return ξi
end




"""
    decouple!(Ωtimes::Vector{Vector{Float64}},
              idf   ::Vector{iBffs},
              ix    ::Int64)

Build a decoupled occurrence record from an occurrence record.
"""
function decouple!(Ωtimes::Vector{Vector{Float64}},
                   ωtimes::Vector{Float64},
                   idf   ::Vector{iBffs})

  for ix in Base.OneTo(lastindex(idf))
    bi = idf[ix]
    tii = ti(bi)
    tfi = tf(bi)
    push!(Ωtimes, filter(ωtime -> tii.>ωtime, ωtimes))
  end
end




"""
    couple(Ωtimes::Vector{Vector{Float64}})

Build an occurrence record from a decoupled occurrence record.
"""
couple(Ωtimes::Vector{Vector{Float64}}) = unique(vcat(Ωtimes...))




"""
    treelength(Ξ::Vector{T}) where {T <: iTree}

Return the branch length sum of `Ξ`.
"""
function treelength(Ξ::Vector{T}) where {T <: iTree}
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
    _ctl(Ξ::Vector{T}, f::Function) where {T <: iTree}

Return the branch length sum of `Ξ` based on `δt` and `fδt`
for debugging purposes.
"""
function _ctl(Ξ::Vector{T}, f::Function) where {T <: iTree}
  L = 0.0
  for ξ in Ξ
    L = _ctl(ξ, f, L)
  end

  return L
end




"""
    nnodesinternal(Ξ::Vector{T}) where {T <: iTree}

Return the number of internal nodes in `Ξ`. 

Warning: only works when branches are not subdivided.
"""
function nnodesinternal(Ξ::Vector{T}) where {T <: iTree}

  n = -0.5
  for ξ in Ξ
    n = _nnodesinternal(ξ, n) + 0.5
  end

  return n
end




"""
    nnodesbifurcation(Ξ::Vector{T}) where {T <: iTree}

Return the number of bifurcating nodes in `Ξ`.

Warning: only works when branches are not subdivided.
"""
function nnodesbifurcation(Ξ::Vector{T}) where {T <: iTree}

  n = -0.5
  for ξ in Ξ
    n = _nnodesinternal(ξ, n) + 0.5 - Int64(anyfossil(ξ))
  end

  return n
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
    _ss_ir_dd(Ξ::Vector{T}, f::Function, α::Float64) where {T <: iTree}

Returns the standardized sum of squares of a diffusion without drift `α`.
"""
function _ss_ir_dd(Ξ::Vector{T}, f::Function, α::Float64) where {T <: iTree}

  dd = ss = n = ir = 0.0
  for ξi in Ξ
    dd, ss, n, ir = _ss_ir_dd(ξi, f, α, dd, ss, n, ir)
  end

  return dd, ss, n, ir
end




"""
    _dd_ss(Ξ::Vector{T}, idf::Vector{iBffs}, α::Float64) where {T <: cT}

Returns the standardized sum of squares of a diffusion without drift `α`.
"""
function _dd_ss(Ξ::Vector{T}, idf::Vector{iBffs}, α::Float64) where {T <: cT}

  dd = ss = 0.0
  for i in Base.OneTo(lastindex(Ξ))

    dd, ss = _dd_ss(Ξ[i], α, dd, ss)

    bi  = idf[i]
    bi2 = d2(bi)

    if bi2 > 0
      lλi = λt(bi)
      lλ1 = lλ(Ξ[d1(bi)])
      lλ2 = lλ(Ξ[bi2])

      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
      dd += lλ1 + lλ2 - 2.0*lλi
    end
  end

  return dd, ss
end



"""
    _dd_ss_seλ(Ξ::Vector{cTct}, idf::Vector{iBffs}, α::Float64)

Returns the standardized sum of squares of a diffusion without drift `α`.
"""
function _dd_ss_seλ(Ξ::Vector{cTct}, idf::Vector{iBffs}, α::Float64)

  dd = ss = seλ = 0.0
  for i in Base.OneTo(lastindex(Ξ))

    dd, ss, seλ = _dd_ss_seλ(Ξ[i], α, dd, ss, seλ)

    bi  = idf[i]
    bi2 = d2(bi)

    if bi2 > 0
      lλi = λt(bi)
      lλ1 = lλ(Ξ[d1(bi)])
      lλ2 = lλ(Ξ[bi2])

      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
      dd += lλ1 + lλ2 - 2.0*lλi
    end
  end

  return dd, ss, seλ
end




"""
    _dd_ss(Ξ::Vector{T}, idf::Vector{iBffs}, α::Float64) where {T <: cT}

Returns the standardized sum of squares of a diffusion without drift `α`.
"""
function _dd_ss(Ξ::Vector{cTbd}, idf::Vector{iBffs}, α::Float64)

  dd = ssλ = ssμ = 0.0
  for i in Base.OneTo(lastindex(Ξ))

    dd, ssλ, ssμ = _dd_ss(Ξ[i], α, dd, ssλ, ssμ)

    bi  = idf[i]
    bi2 = d2(bi)

    if bi2 > 0
      lλi  = λt(bi)
      lμi  = μt(bi)
      ξ1  = Ξ[d1(bi)]
      ξ2  = Ξ[bi2]
      lλ1, lλ2 = lλ(ξ1), lλ(ξ2)

      dd  += lλ1 + lλ2 - 2.0*lλi
      ssλ += 0.5*((lλ1    - lλi - α)^2 + (lλ2    - lλi - α)^2)
      ssμ += 0.5*((lμ(ξ1) - lμi)^2  + (lμ(ξ2) - lμi)^2)
    end
  end

  return dd, ssλ, ssμ
end




"""
    _dd_ss(Ξ::Vector{cTfbd}, idf::Vector{iBffs}, α::Float64)

Returns the standardized sum of squares of a diffusion without drift `α`.
"""
function _dd_ss(Ξ::Vector{cTfbd}, idf::Vector{iBffs}, αλ::Float64, αμ::Float64)

  ddλ = ddμ = ssλ = ssμ = 0.0
  for i in Base.OneTo(lastindex(Ξ))

    ddλ, ddμ, ssλ, ssμ = _dd_ss(Ξ[i], αλ, αμ, ddλ, ddμ, ssλ, ssμ)

    bi  = idf[i]
    bi2 = d2(bi)

    if bi2 > 0
      lλi = λt(bi)
      lμi = μt(bi)
      ξ1  = Ξ[d1(bi)]
      ξ2  = Ξ[bi2]
      lλ1, lλ2 = lλ(ξ1), lλ(ξ2)
      lμ1, lμ2 = lμ(ξ1), lμ(ξ2)

      ddλ += lλ1 + lλ2 - 2.0*lλi
      ddμ += lμ1 + lμ2 - 2.0*lμi
      ssλ += 0.5*((lλ1 - lλi - αλ)^2 + (lλ2 - lλi - αλ)^2)
      ssμ += 0.5*((lμ1 - lμi - αμ)^2 + (lμ2 - lμi - αμ)^2)
    end
  end

  return ddλ, ddμ, ssλ, ssμ
end





"""
    _dd_ss(Ξ::Vector{acTfbd}, idf::Vector{iBffs}, α::Float64)

Returns the standardized sum of squares of a diffusion without drift `α`.
"""
function _dd_ss(Ξ::Vector{acTfbd}, idf::Vector{iBffs}, αλ::Float64, αμ::Float64)

  ddλ = ddμ = ssλ = ssμ = 0.0
  for i in Base.OneTo(lastindex(Ξ))
    ξi = Ξ[i]
    ddλ, ddμ, ssλ, ssμ = _dd_ss(ξi, αλ, αμ, ddλ, ddμ, ssλ, ssμ)

    bi  = idf[i]
    bi2 = d2(bi)
    if bi2 > 0
      lξi = fixtip(ξi)
      ξb  = if sh(lξi) Ξ[d1(bi)] else Ξ[bi2] end
      lλi, lλb = lλ(lξi), lλ(ξb)
      lμi, lμb = lμ(lξi), lμ(ξb)

      ddλ += lλb - lλi
      ddμ += lμb - lμi
      ssλ += 0.5*(lλb - lλi - αλ)^2
      ssμ += 0.5*(lμb - lμi - αμ)^2
    end
  end

  return ddλ, ddμ, ssλ, ssμ
end





"""
    _dd_ss(Ξ::Vector{T}, α::Float64) where {T <: iT}

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function _dd_ss(Ξ::Vector{T}, α::Float64) where {T <: iT}

  dd = ss = n = 0.0
  for ξi in Ξ
    dd, ss, n = _dd_ss(ξi, α, dd, ss, n)
  end

  return dd, ss, n
end




"""
    _dd_ss(Ξ::Vector{iTbd}, α::Float64)

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function _dd_ss(Ξ::Vector{iTbd}, α::Float64)

  dd = ssλ = ssμ = n = 0.0
  for ξi in Ξ
    dd, ssλ, ssμ, n = _dd_ss(ξi, α, dd, ssλ, ssμ, n)
  end

  return dd, ssλ, ssμ, n
end




"""
    _dd_ss(Ξ::Vector{iTfbd}, α::Float64)

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function _dd_ss(Ξ::Vector{iTfbd}, αλ::Float64, αμ::Float64)

  ddλ = ddμ = ssλ = ssμ = n = 0.0
  for ξi in Ξ
    ddλ, ddμ, ssλ, ssμ, n = _dd_ss(ξi, αλ, αμ, ddλ, ddμ, ssλ, ssμ, n)
  end

  return ddλ, ddμ, ssλ, ssμ, n
end



"""
    _dd_ss(Ξ::Vector{T}, α::Float64) where {T <: iTbd}

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function _dd_ss(Ξ::Vector{T}, α::Float64) where {T <: iTbd}

  dd = ssλ = ssμ = n = 0.0
  for ξi in Ξ
    dd, ssλ, ssμ, n = _dd_ss(ξi, α, dd, ssλ, ssμ, n)
  end

  return dd, ssλ, ssμ, n
end




"""
    sss_v(Ξ ::Vector{T}, 
          fx::Function, 
          fσ::Function, 
          αx::Float64) where {T <: iTree}

Returns the delta drift for `fx` and `fσ`, and the standardized sum of 
squares for `fσ` in vector form.
"""
function sss_v(Ξ ::Vector{T}, 
               fx::Function, 
               fσ::Function, 
               ασ::Float64) where {T <: iTree}

  nv  = lastindex(Ξ)
  Lσ  = zeros(nv)
  Δσ  = zeros(nv)
  ddσ = zeros(nv)
  ss  = zeros(nv)
  n   = zeros(nv)
  for i in Base.OneTo(nv)
    Lσ[i], Δσ[i], ddσ[i], ss[i], n[i] = 
      _dd_ss(Ξ[i], fx, fσ, ασ, 0.0, 0.0, 0.0, 0.0, 0.0)
  end

  return Lσ, Δσ, ddσ, ss, n
end




"""
    _dd_ss(Ξ::Vector{T}, α::Float64) where {T <: iTfbd}

Returns the standardized sum of squares a `iT` according
to GBM birth-death for a `σ` proposal.
"""
function _dd_ss(Ξ::Vector{T}, αλ::Float64, αμ::Float64) where {T <: iTfbd}

  ddλ = ddμ = ssλ = ssμ = n = 0.0
  for ξi in Ξ
    ddλ, ddμ, ssλ, ssμ, n = 
      _dd_ss(ξi, αλ, αμ, ddλ, ddμ, ssλ, ssμ, n)
  end

  return ddλ, ddμ, ssλ, ssμ, n
end



"""
    _ss(Ξ::Vector{T}, α::Float64) where {T <: iT}

Returns the standardized sum of squares a for rate `f` a `σ` proposal.
"""
function _ss(Ξ::Vector{T}, f::Function, α::Float64) where {T <: iTree}

  ss = 0.0
  for ξi in Ξ
    ss += _ss(ξi, f, α)
  end

  return ss
end



"""
    _ss!(ss::Vector{Float64}, 
         Ξ ::Vector{T}, 
         f ::Function, 
         α ::Float64) where {T <: iTree}

Returns the standardized sum of squares a for rate `f` a `σ` proposal.
"""
function _ss!(ss::Vector{Float64}, 
              Ξ ::Vector{T}, 
              f ::Function, 
              α ::Float64) where {T <: iTree}

  for i in Base.OneTo(lastindex(Ξ))
    ss[i] = _ss(Ξ[i], f, α)
  end

  return nothing
end




"""
    _ss(Ξ::Vector{T}, α::Float64) where {T <: iT}

Returns the standardized sum of squares a for rate `f` a `σ` proposal.
"""
function _ss(Ξ::Vector{T}, α::Float64) where {T <: iTree}

  ssλ = ssμ = 0.0
  for ξi in Ξ
    ssλ, ssμ = _ss(ξi, α, ssλ, ssμ)
  end

  return ssλ, ssμ
end



"""
    _ss(Ξ::Vector{T}, αλ::Float64, αμ::Float64) where {T <: iTfbd}

Returns the standardized sum of squares a for rate `f` a `σ` proposal.
"""
function _ss(Ξ::Vector{T}, αλ::Float64, αμ::Float64) where {T <: iTfbd}

  ssλ = ssμ = 0.0
  for ξi in Ξ
    ssλ, ssμ = _ss(ξi, αλ, αμ, ssλ, ssμ)
  end

  return ssλ, ssμ
end




"""
    _ss(Ξ::Vector{cTb}, idf::Vector{iBffs}, f::Function, α::Float64)

Returns the standardized sum of squares of a diffusion without drift `α`.
"""
function _ss(Ξ::Vector{T}, idf::Vector{iBffs}, α::Float64) where {T <: cT}

  ss = 0.0
  for i in Base.OneTo(lastindex(Ξ))

    ss = _ss(Ξ[i], α, ss)

    bi  = idf[i]
    bi2 = d2(bi)

    if bi2 > 0
      lλi = λt(bi)
      lλ1 = lλ(Ξ[d1(bi)])
      lλ2 = lλ(Ξ[bi2])
      ss += 0.5*((lλ1 - lλi - α)^2 + (lλ2 - lλi - α)^2)
    end
  end

  return ss
end




"""
    _ss(Ξ  ::Vector{cTfbd},
        idf::Vector{iBffs},
        α  ::Float64,
        fT ::Function,
        fB ::Function)

Returns the standardized sum of squares for drift `α` for function `f`.
"""
function _ss(Ξ  ::Vector{cTfbd},
             idf::Vector{iBffs},
             α  ::Float64,
             fT ::Function,
             fB ::Function)

  ss = 0.0
  for i in Base.OneTo(lastindex(Ξ))

    ss = _ss(Ξ[i], α, ss, fT)

    bi  = idf[i]
    bi2 = d2(bi)

    if bi2 > 0
      fi = fB(bi)
      f1 = fT(Ξ[d1(bi)])
      f2 = fT(Ξ[bi2])
      ss += 0.5*((f1 - fi - α)^2 + (f2 - fi - α)^2)
    end
  end

  return ss
end




"""
    _ss(Ξ  ::Vector{acTfbd},
        idf::Vector{iBffs},
        α  ::Float64,
        f  ::Function)

Returns the standardized sum of squares for drift `α` for function `f`.
"""
function _ss(Ξ  ::Vector{acTfbd},
             idf::Vector{iBffs},
             α  ::Float64,
             f  ::Function)

  ss = 0.0
  for i in Base.OneTo(lastindex(Ξ))
    ξi = Ξ[i]
    ss = _ss(ξi, α, ss, f)

    bi  = idf[i]
    bi2 = d2(bi)
    if bi2 > 0
      lξi    = fixtip(ξi)
      ξb     = if sh(lξi) Ξ[d1(bi)] else Ξ[bi2] end
      fi, fb = f(lξi), f(ξb)

      ss += 0.5*(fb - fi - α)^2
    end
  end

  return ss
end






"""
    _ir(Ξ::Vector{T}) where {T <: cT}

Returns the `integrated rate`.
"""
function _ir(Ξ::Vector{T}) where {T <: cT}

  ir = 0.0
  for ξi in Ξ
    ir = _ir(ξi, ir)
  end

  return ir
end




"""
    _irbd(Ξ::Vector{T}) where {T <: iTree}

Returns the standardized sum of squares a for rate `f` a `σ` proposal.
"""
function _irbd(Ξ::Vector{T}) where {T <: cTbdU}

  irλ = irμ = 0.0
  for ξi in Ξ
    irλ, irμ = _ir(ξi, irλ, irμ)
  end

  return irλ, irμ
end




"""
    _ir(Ξ::Vector{T}) where {T <: iTree}

Returns the standardized sum of squares a for rate `f` a `σ` proposal.
"""
function _ir(Ξ::Vector{T}) where {T <: iTree}

  irλ = irμ = 0.0
  for ξi in Ξ
    irλ, irμ = _ir(ξi, irλ, irμ)
  end

  return irλ, irμ
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




"""
    scale_rate!(Ξ::Vector{T}, f::Function, s::Float64)

Add `s` to vector retrieved using function `f`.
"""
function scale_rate!(Ξ::Vector{T}, f::Function, s::Float64) where {T <: iTree}

  for ξ in Ξ
    scale_rate!(ξ, f, s)
  end

  return nothing
end




"""
    scale_rate!(idf::Vector{iBffs}, f::Function, s::Float64)

Add `s` to vector retrieved using function `f`.
"""
function scale_rate!(idf::Vector{iBffs}, f::Function, s::Float64)

  for bi in idf
    if d2(bi) > 0
      f(bi, s)
    end
  end

  return nothing
end


