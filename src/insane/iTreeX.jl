#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    iTreeX

An abstract type for all composite recursive types representing 
a simple binary phylogenetic tree for `insane` use
"""
abstract type iTreeX <: iTree end



"""
    sT

An abstract type for all composite recursive types representing 
a simple binary phylogenetic tree for `insane` use
"""
abstract type sTX <: iTreeX end




"""
    sTf

An abstract type for all composite recursive types representing 
a simple binary phylogenetic tree with fossils for `insane` use
"""
abstract type sTfX <: sTX end




"""
    sTpbX

The simplest composite recursive type of supertype `sT` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  fx: if fix
  xi: initial trait value
  xf: final trait value

    sTpbX()

Constructs an empty `sTpbX` object.

    sTpbX(e::Float64, fx::Bool, x::Float64)

Constructs an `sTpbX` object with two `sTpbX` daughters and edge `e`, 
fix information `fx`, initial node trait `xi` and final `xf`.
"""
mutable struct sTpbX <: sTX
  d1::sTpbX
  d2::sTpbX
  e ::Float64
  fx::Bool
  xi::Float64
  xf::Float64

  sTpbX() = new()
  sTpbX(e::Float64, fx::Bool, xi::Float64, xf::Float64) = 
    (t = new(); t.e = e; t.fx = fx; t.xi = xi; t.xf = xf; t)
  sTpbX(d1::sTpb, d2::sTpb, e::Float64, fx::Bool, x::Float64) = 
    new(d1, d2, e, fx, xi, xf)
end

# pretty-printing
Base.show(io::IO, t::sTpbX) = 
  print(io, "insane pure-birth and trait tree with ", ntips(t), " tips")











"""
    sT_label

A composite recursive type of supertype `sT` 
representing a labelled binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  l:  label

    sT_label()

Constructs an empty `sT_label` object.

    sT_label(e::Float64)

Constructs an empty `sT_label` object with edge `e`.

    sT_label(d1::sT_label, d2::sT_label, e::Float64)

Constructs an `sT_label` object with two `sT_label` daughters and edge `e`.
"""
mutable struct sT_label <: sT
  d1::sT_label
  d2::sT_label
  e ::Float64
  l ::String

  sT_label() = new()
  sT_label(e::Float64, l::String) = (x = new(); x.e = e; x.l = l; x)
  sT_label(d1::sT_label, 
           d2::sT_label, 
           e ::Float64,
           l ::String) = new(d1, d2, e, l)
end

# pretty-printing
Base.show(io::IO, t::sT_label) = 
  print(io, "insane simple labelled tree with ", ntips(t), " tips")




"""
    sT_label(tree::T) where {T <: iTree}

Demotes a tree to `sT_label`.
"""
function sT_label(tree::T) where {T <: iTree}
  _sT_label(tree::T, 0)[1]
end




"""
    sT_label(tree::T) where {T <: iTree}

Demotes a tree to `sT_label`.
"""
function _sT_label(tree::T, i::Int64) where {T <: iTree}
  if isdefined(tree, :d1)
    t1, i = _sT_label(tree.d1, i)
    t2, i = _sT_label(tree.d2, i)
    tree = sT_label(t1, t2, e(tree), "")
  else
    i += 1
    tree = sT_label(e(tree), string("t",i))
  end
  return tree, i
end




"""
    sTf_label

A composite recursive type of supertype `sTf` representing a labelled binary 
phylogenetic tree with fossils for `insane` use, with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: if it is extinct
  iψ: if it is a fossil
  l:  label

    sTf_label()

Constructs an empty `sTf_label` object.

    sTf_label(e::Float64, (iμ::Bool), (iψ::Bool), (l ::String))

Constructs an empty `sTf_label` object with edge `e` and label `l`.

    sTf_label(d1::sTf_label e::Float64, l ::String)

Constructs an `sTf_label` object with one `sTf_label` daughter and edge `e`.

    sTf_label(d1::sTf_label, d2::sTf_label, e::Float64, l ::String)

Constructs an `sTf_label` object with two `sTf_label` daughters and edge `e`.
"""
mutable struct sTf_label <: sTf
  d1::sTf_label
  d2::sTf_label
  e ::Float64
  iμ::Bool
  iψ::Bool
  l ::String

  sTf_label() = new()
  sTf_label(e::Float64) = (x=new(); x.e=e; x.iμ=false; x.iψ=false; x.l=""; x)
  sTf_label(e::Float64, 
            l::String)  = (x=new(); x.e=e; x.iμ=false; x.iψ=false; x.l=l; x)
  sTf_label(e ::Float64, 
            iμ::Bool,
            iψ::Bool)   = (x=new(); x.e=e; x.iμ=iμ; x.iψ=iψ; x.l=""; x)
  sTf_label(e ::Float64, 
            iμ::Bool,
            iψ::Bool,
            l ::String) = (x=new(); x.e=e; x.iμ=iμ; x.iψ=iψ; x.l=l; x)
  sTf_label(d1::sTf_label, 
            e ::Float64,
            l ::String) = (x=new(); x.d1=d1; x.e=e; x.iμ=false; x.iψ=true; 
                           x.l=l; x)
  sTf_label(d1::sTf_label, 
            d2::sTf_label, 
            e ::Float64,
            l ::String) = (x=new(); x.d1=d1; x.d2=d2; x.e=e; x.iμ=false; 
                           x.iψ=false; x.l=l; x)
end

# pretty-printing
Base.show(io::IO, t::sTf_label) = 
  print(io, "insane simple labelled tree with ", ntips(t), " tips (", 
            ntipsextinct(t)," extinct) and ", nfossils(t)," fossils")




"""
    sTf_label(tree::T) where {T <: iTree}

Demotes a tree to `sTf_label`.
"""
function sTf_label(tree::T) where {T <: iTree}
  _sTf_label(tree::T, 0)[1]
end




"""
    _sTf_label(tree::T, i::Int64) where {T <: iTree}

Demotes a tree to `sTf_label`, initialized with label i.
"""
function _sTf_label(tree::T, i::Int64) where {T <: iTree}
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if defd1 && defd2
    t1, i = _sTf_label(tree.d1, i)
    t2, i = _sTf_label(tree.d2, i)
    tree = sTf_label(t1, t2, e(tree), l(tree))
  elseif defd1
    t1, i = _sTf_label(tree.d1, i)
    tree = sTf_label(t1, e(tree), l(tree))
  elseif defd2
    t2, i = _sTf_label(tree.d2, i)
    tree = sTf_label(t2, e(tree), l(tree))
  else
    i += 1
    lab = ifelse(isempty(l(tree)),string("t",i),l(tree))
    if isdefined(tree, :iμ)
      tree = sTf_label(e(tree), isextinct(tree), isfossil(tree), lab)
    else
      tree = sTf_label(e(tree), lab)
    end
  end
  return tree, i
end




"""
    sTpb(tree::sTpb)

Demotes a tree of type `sT_label` to `sTpb`.
"""
function sTpb(tree::sTpb)
  if isdefined(tree, :d1)
    sTpb(sTpb(tree.d1), sTpb(tree.d2), e(tree), isfix(tree))
  else
    sTpb(e(tree), isfix(tree))
  end
end



"""
    sTbd

The simplest composite recursive type of supertype `sT` 
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
    (x = new(); x.e = e; x.iμ = iμ; x.fx = false; x)
  sTbd(e::Float64, iμ::Bool, fx::Bool) = 
    (x = new(); x.e = e; x.iμ = iμ; x.fx = fx; x)
  sTbd(d1::sTbd, d2::sTbd, e::Float64) = 
    new(d1, d2, e, false, false) 
  sTbd(d1::sTbd, d2::sTbd, e::Float64, iμ::Bool, fx::Bool) = 
    new(d1, d2, e, iμ, fx)
end

# pretty-printing
Base.show(io::IO, t::sTbd) = 
  print(io, "insane simple birth-death tree with ", ntips(t), " tips (", 
    ntipsextinct(t)," extinct)")





"""
    sTbd(tree::sT_label)

Transforms a tree of type `sT_label` to `sTbd`.
"""
function sTbd(tree::sT_label)
  if isdefined(tree, :d1)
    sTbd(sTbd(tree.d1), sTbd(tree.d2), e(tree), false, false)
  else
    sTbd(e(tree), false)
  end
end




"""
    sTbd(tree::sTbd)

Produce a new copy of `sTbd`.
"""
function sTbd(tree::sTbd)
  if isdefined(tree, :d1)
    sTbd(sTbd(tree.d1), sTbd(tree.d2), e(tree), isextinct(tree), isfix(tree))
  else
    sTbd(e(tree), isextinct(tree), isfix(tree))
  end
end




"""
    sTfbd

The simplest composite recursive type of supertype `sTf` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: is an extinction node
  iψ: is a fossil node
  fx: if it is fix

    sTfbd()

Constructs an empty `sTfbd` object.

    sTfbd(e::Float64)

Constructs an empty `sTfbd` object with edge `e`.

    sTfbd(d1::sTfbd, d2::sTfbd, e::Float64)

Constructs an `sTfbd` object with two `sTfbd` daughters and edge `e`.

    sTfbd(d1::sTfbd, e::Float64)

Constructs an `sTfbd` object with one sampled ancestor, one `sTfbd` daughter and 
edge `e`.
"""
mutable struct sTfbd <: sTf
  d1::sTfbd
  d2::sTfbd
  e ::Float64
  iμ::Bool
  iψ::Bool
  fx::Bool

  sTfbd() = new()
  sTfbd(e::Float64) = 
    (x = new(); x.e = e; x.iμ = false; x.iψ = false; x.fx = false; x)
  sTfbd(e::Float64, iμ::Bool) = 
    (x = new(); x.e = e; x.iμ = iμ; x.iψ = false; x.fx = false; x)
  sTfbd(e::Float64, iμ::Bool, iψ::Bool) = 
    (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = false; x)
  sTfbd(e::Float64, iμ::Bool, iψ::Bool, fx::Bool) = 
    (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x)
  sTfbd(d1::sTfbd, d2::sTfbd, e::Float64) = 
    new(d1, d2, e, false, false, false)
  sTfbd(d1::sTfbd, e::Float64) = 
    (x = new(); x.d1 = d1; x.e = e; x.iμ = false; x.iψ = false; x.fx = false; x)
  sTfbd(d1::sTfbd, d2::sTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool) = 
    new(d1, d2, e, iμ, iψ, fx)
  sTfbd(d1::sTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool) = 
    (x = new(); x.d1 = d1; x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x)
end

# pretty-printing
Base.show(io::IO, t::sTfbd) = 
  print(io, "insane simple fossilized birth-death tree with ", ntips(t), 
    " tips (", ntipsextinct(t)," extinct) and ", nfossils(t)," fossils")



"""
    sTfbd(tree::sTfbd)

Transforms a tree of type `sT_label` to `sTfbd`.
"""
function sTfbd(tree::sTfbd)
  if isdefined(tree, :d1)
    sTfbd(sTfbd(tree.d1), sTfbd(tree.d2), e(tree), 
      isextinct(tree), isfossil(tree),  isfix(tree))
  else
    sTfbd(e(tree), isextinct(tree), isfossil(tree),  isfix(tree))
  end
end




"""
    sTfbd(tree::sTf_label)

Transforms a tree of type `sTf_label` to `sTfbd`.
"""
function sTfbd(tree::sTf_label)
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if defd1 && defd2
    sTfbd(sTfbd(tree.d1), sTfbd(tree.d2), e(tree), false, false, false)
  elseif defd1
    sTfbd(sTfbd(tree.d1), e(tree), false, true, false)
  elseif defd2
    sTfbd(sTfbd(tree.d2), e(tree), false, true, false)
  else
    sTfbd(e(tree), false)
  end
end




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
  fx:  if fix tree
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
  fx ::Bool
  dt ::Float64
  fdt::Float64
  lλ ::Array{Float64,1}

  iTgbmpb() = new()

  iTgbmpb(e::Float64, fx::Bool, dt::Float64, fdt::Float64, 
    lλ::Array{Float64,1}) = 
      (x = new(); x.e = e; x.fx = fx; x.dt = dt; x.fdt = fdt; x.lλ = lλ; x)

  iTgbmpb(d1::iTgbmpb, d2::iTgbmpb, e::Float64, fx::Bool,
    dt::Float64, fdt::Float64, lλ::Array{Float64,1}) = 
      new(d1, d2, e, fx, dt, fdt, lλ)
end


# pretty-printing
Base.show(io::IO, t::iTgbmpb) = 
  print(io, "insane pb-gbm tree with ", ntips(t), " tips")



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
              et, isfix(tree), δt, 0.0, Float64[lλa, lλa])
    else
      iTgbmpb(et, isfix(tree), δt, 0.0, Float64[lλa, lλa])
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
              et, isfix(tree), δt, fdti, lλv)
    else
      iTgbmpb(et, isfix(tree), δt, fdti, lλv)
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
    return iTgbmpb(el[ei], true, δt, δt, λs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTgbmpb(iTgbmpb(e0, e1, el, λs, ea, n1, ei1, δt), 
                   iTgbmpb(e0, e1, el, λs, ea, n2, ei2, δt), 
                   el[ei], true, δt, (el[ei] == 0.0 ? 0.0 : δt), λs[ei])
  end
end




"""
    iTgbmpb(tree::iTgbmpb)

Produce a new copy of `iTgbmpb`.
"""
function iTgbmpb(tree::iTgbmpb)
  if isdefined(tree, :d1)
    iTgbmpb(iTgbmpb(tree.d1), iTgbmpb(tree.d2), 
      e(tree), isfix(tree), dt(tree), fdt(tree), copy(lλ(tree)))
  else
    iTgbmpb(e(tree), isfix(tree), dt(tree), fdt(tree), copy(lλ(tree)))
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
  print(io, "insane gbm-ce tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    iTgbmce(tree::iTgbmce)

Produce a new copy of `iTgbmce`.
"""
function iTgbmce(tree::iTgbmce)
  if isdefined(tree, :d1)
    iTgbmce(iTgbmce(tree.d1), iTgbmce(tree.d2), 
      e(tree), dt(tree), fdt(tree), isextinct(tree), 
      isfix(tree), copy(lλ(tree)))
  else
    iTgbmce(e(tree), dt(tree), fdt(tree), isextinct(tree), 
      isfix(tree), copy(lλ(tree)))
  end
end




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
  print(io, "insane gbm-ct tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    iTgbmct(tree::iTgbmct)

Produce a new copy of `iTgbmct`.
"""
function iTgbmct(tree::iTgbmct)
  if isdefined(tree, :d1)
    iTgbmct(iTgbmct(tree.d1), iTgbmct(tree.d2), 
      e(tree), dt(tree), fdt(tree), isextinct(tree), 
      isfix(tree), copy(lλ(tree)))
  else
    iTgbmct(e(tree), dt(tree), fdt(tree), isextinct(tree), 
      isfix(tree), copy(lλ(tree)))
  end
end




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
  print(io, "insane bd-gbm tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    iTgbmbd(tree::iTgbmbd)

Produce a new copy of `iTgbmbd`.
"""
function iTgbmbd(tree::iTgbmbd)
  if isdefined(tree, :d1)
    iTgbmbd(iTgbmbd(tree.d1), iTgbmbd(tree.d2), 
      e(tree), dt(tree), fdt(tree), isextinct(tree), 
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
  else
    iTgbmbd(e(tree), dt(tree), fdt(tree), isextinct(tree), 
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
  end
end




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
