#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    iTree

An abstract type for all composite recursive types 
representing a binary phylogenetic tree for `insane` use
"""
abstract type iTree end




"""
    sT

An abstract type for all composite recursive types 
representing a simple binary phylogenetic tree for `insane` use
"""
abstract type sT <: iTree end




"""
    sTpb

The simplest composite recursive type of supertype `iTree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge

    sTpb()

Constructs an empty `sTpb` object.

    sTpb(e::Float64)

Constructs an empty `sTpb` object with edge `e`.

    sTpb(d1::sTpb, d2::sTpb, e::Float64)

Constructs an `sTpb` object with two `sTpb` daughters and edge `e`.
"""
mutable struct sTpb <: sT
  d1::sTpb
  d2::sTpb
  e ::Float64

  sTpb() = new()
  sTpb(e::Float64) = (x = new(); x.e = e; x)
  sTpb(d1::sTpb, d2::sTpb, e::Float64) = new(d1, d2, e)
end

# pretty-printing
Base.show(io::IO, t::sTpb) = 
  print(io, "insane simple pure-birth tree with ", sntn(t,0), " tips")




"""
    sTbd

The simplest composite recursive type of supertype `sTbdree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: is an extinction node
  fx: if it is fix

    sTbd()

Constructs an empty `sTbd` object.

    sTbd(e::Float64)

Constructs an empty `sTbd` object with edge `e`.

    sTbd(d1::sTbd, d2::sTbd, e::Float64)

Constructs an `sTbd` object with two `sTbd` daughters and edge `e`.
"""
mutable struct sTbd <: sT
  d1::sTbd
  d2::sTbd
  e ::Float64
  iμ::Bool
  fx::Bool

  sTbd() = new()
  sTbd(e::Float64) = 
    (x = new(); x.e = e; x.iμ = false; x.fx = false; x)
  sTbd(e::Float64, iμ::Bool) = 
    (x = new(); x.e = e; x.iμ = iμ; x)
  sTbd(e::Float64, iμ::Bool, fx::Bool) = 
    (x = new(); x.e = e; x.iμ = iμ; x.fx = fx; x)
  sTbd(d1::sTbd, d2::sTbd, e::Float64) = 
    new(d1, d2, e, false, false) 
  sTbd(d1::sTbd, d2::sTbd, e::Float64, iμ::Bool, fx::Bool) = 
    new(d1, d2, e, iμ, fx)
end

# pretty-printing
Base.show(io::IO, t::sTbd) = 
  print(io, "insane simple birth-death tree with ", sntn(t,0), " tips (", snen(t,0)," extinct)")




"""
    iTgbm

An abstract type for all composite recursive types 
representing a binary phylogenetic tree with Geometric
Brownian motion rates for `insane` use
"""
abstract type iTgbm <: iTree end




"""
    iTgbmpb

A composite recursive type of supertype `iTgbm` 
representing a binary phylogenetic tree with no extinction
and `λ` evolving as a Geometric Brownian motion  for `insane` use, 
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   edge
  dt:  choice of time lag
  fdt: final `δt`
  lλ:  array of a Brownian motion evolution of `log(λ)`

    iTgbmpb()

Constructs an empty `iTgbmpb` object.

    iTgbmpb(e::Float64)

Constructs an empty `iTgbmpb` object with pendant edge `pe`.

    iTgbmpb(d1::iTgbmpb, d2::iTgbmpb, e::Float64)

Constructs an `iTgbmpb` object with two `iTgbmpb` daughters and pendant edge `pe`.
"""
mutable struct iTgbmpb <: iTgbm
  d1 ::iTgbmpb
  d2 ::iTgbmpb
  e  ::Float64
  dt ::Float64
  fdt::Float64
  lλ ::Array{Float64,1}

  iTgbmpb() = new()

  iTgbmpb(e::Float64, dt::Float64, fdt::Float64, lλ::Array{Float64,1}) = 
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt; x.lλ = lλ; x)

  iTgbmpb(d1::iTgbmpb, d2::iTgbmpb, e::Float64, 
    dt::Float64, fdt::Float64, lλ::Array{Float64,1}) = 
      new(d1, d2, e, dt, fdt, lλ)
end


# pretty-printing
Base.show(io::IO, t::iTgbmpb) = 
  print(io, "insane pb-gbm tree with ", sntn(t,0), " tips")



"""
    iTgbmpb(tree::sTpb, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64,
            α   ::Float64,
            σλ  ::Float64)

Promotes an `sTpb` to `iTgbmpb` according to some values for `λ` diffusion.
"""
function iTgbmpb(tree::sTpb, 
                 δt  ::Float64, 
                 srδt::Float64, 
                 lλa ::Float64,
                 α   ::Float64,
                 σλ  ::Float64)

  et = e(tree)

  if iszero(et)
    if isdefined(tree, :d1)
      iTgbmpb(iTgbmpb(tree.d1, δt, srδt, lλa, α, σλ), 
              iTgbmpb(tree.d2, δt, srδt, lλa, α, σλ),
              et, δt, 0.0, Float64[lλa, lλa])
    else
      iTgbmpb(et, δt, 0.0, Float64[lλa, lλa])
    end
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
    end

    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)
    l   = lastindex(lλv)

    if isdefined(tree, :d1)
      iTgbmpb(iTgbmpb(tree.d1, δt, srδt, lλv[l], α, σλ), 
              iTgbmpb(tree.d2, δt, srδt, lλv[l], α, σλ),
              et, δt, fdti, lλv)
    else
      iTgbmpb(et, δt, fdti, lλv)
    end
  end
end





"""
    iTgbmpb(e0::Array{Int64,1}, 
            e1::Array{Int64,1}, 
            el::Array{Float64,1}, 
            λs::Array{Array{Float64,1},1}, 
            ea::Array{Int64,1}, 
            ni::Int64, 
            ei::Int64,
            δt::Float64)

Transform edge structure to `iTgbmpb`.
"""
function iTgbmpb(e0::Array{Int64,1}, 
                 e1::Array{Int64,1}, 
                 el::Array{Float64,1}, 
                 λs::Array{Array{Float64,1},1}, 
                 ea::Array{Int64,1}, 
                 ni::Int64, 
                 ei::Int64,
                 δt::Float64)

  # if tip
  if in(ei, ea)
    return iTgbmpb(el[ei], δt, δt, λs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTgbmpb(iTgbmpb(e0, e1, el, λs, ea, n1, ei1, δt), 
                   iTgbmpb(e0, e1, el, λs, ea, n2, ei2, δt), 
                   el[ei], δt, 
                   (el[ei] == 0.0 ? 0.0 : δt), 
                   λs[ei])
  end
end




"""
    sTpb(tree::iTgbmpb)

Demotes a tree of type `iTgbmce` to `sTpb`.
"""
function sTpb(tree::iTgbmpb)
  if isdefined(tree, :d1)
    sTpb(sTpb(tree.d1), sTpb(tree.d2), e(tree))
  else
    sTpb(e(tree))
  end
end




"""
    iTgbmce

A composite recursive type of supertype `iTgbm` 
representing a binary phylogenetic tree with  `λ` evolving as a 
Geometric Brownian motion and constant `μ` for `insane` use, 
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    edge
  iμ:   if extinct node
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`

  iTgbmce(d1 ::iTgbmce, 
          d2 ::iTgbmce, 
          e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1})
"""
mutable struct iTgbmce <: iTgbm
  d1 ::iTgbmce
  d2 ::iTgbmce
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}

  iTgbmce() = new()

  iTgbmce(e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1}) = 
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt; 
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x)

  iTgbmce(d1 ::iTgbmce, 
          d2 ::iTgbmce, 
          e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1}) = 
    new(d1, d2, e, dt, fdt, iμ, fx, lλ)
end

# pretty-printing
Base.show(io::IO, t::iTgbmce) = 
  print(io, "insane gbm-ce tree with ", sntn(t,0), " tips (", snen(t,0)," extinct)")




"""
    sTbd(tree::iTgbmce)

Demotes a tree of type `iTgbmce` to `sTbd`.
"""
function sTbd(tree::iTgbmce)
  if isdefined(tree, :d1)
    sTbd(sTbd(tree.d1), sTbd(tree.d2), e(tree), isextinct(tree), false)
  else
    sTbd(e(tree), isextinct(tree))
  end
end




"""
    iTgbmce(tree::sTbd, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            σλ  ::Float64,
            σμ  ::Float64)

Promotes an `sTbd` to `iTgbmce` according to some values for `λ` and `μ` 
diffusion.
"""
function iTgbmce(tree::sTbd, 
                 δt  ::Float64, 
                 srδt::Float64, 
                 lλa ::Float64,
                 α   ::Float64,
                 σλ  ::Float64)

  et = e(tree)

  if iszero(et)
    if isdefined(tree, :d1)
      iTgbmce(iTgbmce(tree.d1, δt, srδt, lλa, α, σλ), 
              iTgbmce(tree.d2, δt, srδt, lλa, α, σλ),
              et, δt, 0.0, isextinct(tree), isfix(tree), Float64[lλa, lλa])
    else
      iTgbmce(et, δt, 0.0, isextinct(tree), isfix(tree), Float64[lλa, lλa])
    end

  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)

    if iszero(fdti)
      fdti = δt
    end

    l = lastindex(lλv)

    if isdefined(tree, :d1)
      iTgbmce(iTgbmce(tree.d1, δt, srδt, lλv[l], α, σλ), 
              iTgbmce(tree.d2, δt, srδt, lλv[l], α, σλ),
              et, δt, fdti, isextinct(tree), isfix(tree), lλv)
    else
      iTgbmce(et, δt, fdti, isextinct(tree), isfix(tree), lλv)
    end
  end
end




"""
    iTgbmce(e0::Array{Int64,1}, 
            e1::Array{Int64,1}, 
            el::Array{Float64,1}, 
            λs::Array{Array{Float64,1},1}, 
            ea::Array{Int64,1}, 
            ee::Array{Int64,1},
            ni::Int64, 
            ei::Int64,
            δt::Float64)

Transform edge structure to `iTgbmce`.
"""
function iTgbmce(e0::Array{Int64,1}, 
                 e1::Array{Int64,1}, 
                 el::Array{Float64,1}, 
                 λs::Array{Array{Float64,1},1}, 
                 ea::Array{Int64,1}, 
                 ee::Array{Int64,1},
                 ni::Int64, 
                 ei::Int64,
                 δt::Float64)

  # if tip
  if in(ei, ea)
    return iTgbmce(el[ei], δt, δt, false, false, λs[ei])
  # if extinct
  elseif in(ei, ee)
    return iTgbmce(el[ei], δt, δt, true, false, λs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTgbmce(iTgbmce(e0, e1, el, λs, ea, ee, n1, ei1, δt), 
                   iTgbmce(e0, e1, el, λs, ea, ee, n2, ei2, δt), 
                   el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), false, false, λs[ei])
  end
end




"""
    iTgbmct

A composite recursive type of supertype `iTgbm` 
representing a binary phylogenetic tree with  `λ` evolving as a 
Geometric Brownian motion and constant `μ` for `insane` use, 
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    edge
  iμ:   if extinct node
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`

  iTgbmct(d1 ::iTgbmct, 
          d2 ::iTgbmct, 
          e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1})
"""
mutable struct iTgbmct <: iTgbm
  d1 ::iTgbmct
  d2 ::iTgbmct
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}

  iTgbmct() = new()

  iTgbmct(e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1}) = 
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt; 
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x)

  iTgbmct(d1 ::iTgbmct, 
          d2 ::iTgbmct, 
          e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1}) = 
    new(d1, d2, e, dt, fdt, iμ, fx, lλ)
end


# pretty-printing
Base.show(io::IO, t::iTgbmct) = 
  print(io, "insane gbm-ct tree with ", sntn(t,0), " tips (", snen(t,0)," extinct)")




"""
    sTbd(tree::iTgbmct)

Demotes a tree of type `iTgbmct` to `sTbd`.
"""
function sTbd(tree::iTgbmct)
  if isdefined(tree, :d1)
    sTbd(sTbd(tree.d1), sTbd(tree.d2), e(tree), isextinct(tree), false)
  else
    sTbd(e(tree), isextinct(tree))
  end
end




"""
    iTgbmct(tree::sTbd, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            σλ  ::Float64,
            σμ  ::Float64)

Promotes an `sTbd` to `iTgbmct` according to some values for `λ` and `μ` 
diffusion.
"""
function iTgbmct(tree::sTbd, 
                 δt  ::Float64, 
                 srδt::Float64, 
                 lλa ::Float64, 
                 α   ::Float64,
                 σλ  ::Float64)

  et = e(tree)

  # if crown root
  if iszero(et)
    if isdefined(tree, :d1)
      iTgbmct(iTgbmct(tree.d1, δt, srδt, lλa, α, σλ), 
              iTgbmct(tree.d2, δt, srδt, lλa, α, σλ),
              et, δt, 0.0, isextinct(tree), isfix(tree), Float64[lλa, lλa])
    else
      iTgbmct(et, δt, 0.0, isextinct(tree), isfix(tree), Float64[lλa, lλa])
    end
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    lλv = sim_bm(lλa, α, σλ, δt, fdti, srδt, nt)

    if iszero(fdti)
      fdti = δt
    end

    l = lastindex(lλv)

    if isdefined(tree, :d1)
      iTgbmct(iTgbmct(tree.d1, δt, srδt, lλv[l], α, σλ), 
              iTgbmct(tree.d2, δt, srδt, lλv[l], α, σλ),
              et, δt, fdti, isextinct(tree), isfix(tree), lλv)
    else
      iTgbmct(et, δt, fdti, isextinct(tree), isfix(tree), lλv)
    end
  end
end




"""
    iTgbmct(e0::Array{Int64,1}, 
            e1::Array{Int64,1}, 
            el::Array{Float64,1}, 
            λs::Array{Array{Float64,1},1}, 
            ea::Array{Int64,1}, 
            ee::Array{Int64,1},
            ni::Int64, 
            ei::Int64,
            δt::Float64)

Transform edge structure to `iTgbmct`.
"""
function iTgbmct(e0::Array{Int64,1}, 
                 e1::Array{Int64,1}, 
                 el::Array{Float64,1}, 
                 λs::Array{Array{Float64,1},1}, 
                 ea::Array{Int64,1}, 
                 ee::Array{Int64,1},
                 ni::Int64, 
                 ei::Int64,
                 δt::Float64)

  # if tip
  if in(ei, ea)
    return iTgbmct(el[ei], δt, δt, false, false, λs[ei])
  # if extinct
  elseif in(ei, ee)
    return iTgbmct(el[ei], δt, δt, true, false, λs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTgbmct(iTgbmct(e0, e1, el, λs, ea, ee, n1, ei1, δt), 
                   iTgbmct(e0, e1, el, λs, ea, ee, n2, ei2, δt), 
                   el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), false, false, λs[ei])
  end
end




"""
    iTgbmbd

A composite recursive type of supertype `iTgbm` 
representing a binary phylogenetic tree with  `λ` and `μ` 
evolving as a Geometric Brownian motion  for `insane` use, 
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    pendant edge
  iμ:   if extinct node
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`
  lμ:   array of a Brownian motion evolution of `log(μ)`

  iTgbmbd(d1 ::iTgbmbd, 
          d2 ::iTgbmbd, 
          e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1},
          lμ ::Array{Float64,1})
"""
mutable struct iTgbmbd <: iTgbm
  d1 ::iTgbmbd
  d2 ::iTgbmbd
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  lμ ::Array{Float64,1}

  iTgbmbd() = new()

  iTgbmbd(e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1},
          lμ ::Array{Float64,1}) = 
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt; 
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)

  iTgbmbd(d1 ::iTgbmbd, 
          d2 ::iTgbmbd, 
          e  ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1},
          lμ ::Array{Float64,1}) = 
    new(d1, d2, e, dt, fdt, iμ, fx, lλ, lμ)
end


# pretty-printing
Base.show(io::IO, t::iTgbmbd) = 
  print(io, "insane bd-gbm tree with ", sntn(t,0), " tips (", snen(t,0)," extinct)")




"""
    sTbd(tree::iTgbmbd)

Demotes a tree of type `iTgbmbd` to `sTbd`.
"""
function sTbd(tree::iTgbmbd)
  if isdefined(tree, :d1)
    sTbd(sTbd(tree.d1), sTbd(tree.d2), e(tree), isextinct(tree), false)
  else
    sTbd(e(tree), isextinct(tree), false)
  end
end




"""
    iTgbmbd(tree::sTbd, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            lμa ::Float64, 
            α   ::Float64,
            σλ  ::Float64,
            σμ  ::Float64)

Promotes an `sTbd` to `iTgbmbd` according to some values for `λ` and `μ` 
diffusion.
"""
function iTgbmbd(tree::sTbd, 
                 δt  ::Float64, 
                 srδt::Float64, 
                 lλa ::Float64, 
                 lμa ::Float64, 
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64)

  et  = e(tree)

  # if crown root
  if iszero(et)
    if isdefined(tree, :d1)
      iTgbmbd(iTgbmbd(tree.d1, δt, srδt, lλa, lμa, α, σλ, σμ), 
              iTgbmbd(tree.d2, δt, srδt, lλa, lμa, α, σλ, σμ),
              e(tree), δt, 0.0, isextinct(tree), isfix(tree), 
              Float64[lλa, lλa], Float64[lμa, lμa])
    else
      iTgbmbd(e(tree), δt, 0.0, isextinct(tree), isfix(tree), 
              Float64[lλa, lλa], Float64[lμa, lμa])
    end
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    lλv = sim_bm(lλa, α,   σλ, δt, fdti, srδt, nt)
    lμv = sim_bm(lμa, 0.0, σμ, δt, fdti, srδt, nt)

    if iszero(fdti)
      fdti = δt
    end

    l = lastindex(lμv)

    if isdefined(tree, :d1)
      iTgbmbd(iTgbmbd(tree.d1, δt, srδt, lλv[l], lμv[l], α, σλ, σμ), 
              iTgbmbd(tree.d2, δt, srδt, lλv[l], lμv[l], α, σλ, σμ),
              e(tree), δt, fdti, isextinct(tree), isfix(tree), lλv, lμv)
    else
      iTgbmbd(e(tree), δt, fdti, isextinct(tree), isfix(tree), lλv, lμv)
    end
  end
end




"""
    iTgbmbd(e0::Array{Int64,1}, 
            e1::Array{Int64,1}, 
            el::Array{Float64,1}, 
            λs::Array{Array{Float64,1},1}, 
            μs::Array{Array{Float64,1},1}, 
            ea::Array{Int64,1}, 
            ee::Array{Int64,1},
            ni::Int64, 
            ei::Int64,
            δt::Float64)

Transform edge structure to `iTgbmbd`.
"""
function iTgbmbd(e0::Array{Int64,1}, 
                 e1::Array{Int64,1}, 
                 el::Array{Float64,1}, 
                 λs::Array{Array{Float64,1},1}, 
                 μs::Array{Array{Float64,1},1}, 
                 ea::Array{Int64,1}, 
                 ee::Array{Int64,1},
                 ni::Int64, 
                 ei::Int64,
                 δt::Float64)

  # if tip
  if in(ei, ea)
    return iTgbmbd(el[ei], δt, δt, false, false, λs[ei], μs[ei])
  # if extinct
  elseif in(ei, ee)
    return iTgbmbd(el[ei], δt, δt, true, false, λs[ei], μs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTgbmbd(iTgbmbd(e0, e1, el, λs, μs, ea, ee, n1, ei1, δt), 
                   iTgbmbd(e0, e1, el, λs, μs, ea, ee, n2, ei2, δt), 
                   el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), 
                   false, false, λs[ei], μs[ei])
  end
end
