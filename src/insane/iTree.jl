#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    iTree

An abstract type for all composite recursive types representing
a binary phylogenetic tree for `insane` use
"""
abstract type iTree end




"""
    sT

An abstract type for all composite recursive types representing
a simple binary phylogenetic tree for `insane` use
"""
abstract type sT <: iTree end





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

  sT_label(e::Float64, l::String) =
    (x = new(); x.e = e; x.l = l; x)
  sT_label(d1::sT_label, d2::sT_label, e::Float64, l::String) =
    new(d1, d2, e, l)
end

# pretty-printing
function Base.show(io::IO, t::sT_label)
  nt = ntips(t)
  print(io, "insane simple labelled tree with ", nt, " tip",
    (isone(nt) ? "" : "s"))
end




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
  if def1(tree)
    t1, i = _sT_label(tree.d1, i)
    t2, i = _sT_label(tree.d2, i)
    tree  = sT_label(t1, t2, e(tree), "")
  else
    i += 1
    tree = sT_label(e(tree), string("t",i))
  end
  return tree, i
end




"""
    sTf_label

A composite recursive type of supertype `sT` representing a labelled binary
phylogenetic tree with fossils for `insane` use, with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: if it is extinct
  iψ: if it is a fossil
  l:  label

    sTf_label(e::Float64, iμ::Bool, iψ::Bool, l::String)

Constructs an empty `sTf_label` object with edge `e` and label `l`.

    sTf_label(d1::sTf_label e::Float64, l::String)

Constructs an `sTf_label` object with one `sTf_label` daughter and edge `e`.

    sTf_label(d1::sTf_label, d2::sTf_label, e::Float64, l ::String)

Constructs an `sTf_label` object with two `sTf_label` daughters and edge `e`.
"""
mutable struct sTf_label <: sT
  d1::sTf_label
  d2::sTf_label
  e ::Float64
  iμ::Bool
  iψ::Bool
  l ::String

  sTf_label(e::Float64, l::String) =
    (x=new(); x.e=e; x.iμ=false; x.iψ=false; x.l=l; x)
  sTf_label(e::Float64, iμ::Bool, iψ::Bool) =
    (x=new(); x.e=e; x.iμ=iμ; x.iψ=iψ; x.l=""; x)
  sTf_label(e::Float64, iμ::Bool, iψ::Bool, l::String) =
    (x=new(); x.e=e; x.iμ=iμ; x.iψ=iψ; x.l=l; x)
  sTf_label(d1::sTf_label, e::Float64, l::String) =
    (x=new(); x.d1=d1; x.e=e; x.iμ=false; x.iψ=true; x.l=l; x)
  sTf_label(d1::sTf_label, d2::sTf_label, e::Float64, l::String) =
    new(d1, d2, e, false, false, l)
end

# pretty-printing
Base.show(io::IO, t::sTf_label) =
  print(io, "insane simple fossil labelled tree with ", ntips(t),
    " tips (", ntipsextinct(t)," extinct) and ", nfossils(t)," fossils")




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

  if def1(tree)
    if def2(tree)
      t1, i = _sTf_label(tree.d1, i)
      t2, i = _sTf_label(tree.d2, i)
      tree  = sTf_label(t1, t2, e(tree), l(tree))
    else
      t1, i = _sTf_label(tree.d1, i)
      tree = sTf_label(t1, e(tree), l(tree))
    end
  elseif def2(tree)
    t2, i = _sTf_label(tree.d2, i)
    tree = sTf_label(t2, e(tree), l(tree))
  else
    i += 1
    lab = isempty(l(tree)) ? string("t",i) : l(tree)
    if isdefined(tree, :iμ)
      tree = sTf_label(e(tree), isextinct(tree), isfossil(tree), lab)
    else
      tree = sTf_label(e(tree), lab)
    end
  end

  return tree, i
end




"""
    sTpb

The simplest composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  fx: if fix

    sTpb(e::Float64)

Constructs an empty `sTpb` object with edge `e`.

    sTpb(d1::sTpb, d2::sTpb, e::Float64)

Constructs an `sTpb` object with two `sTpb` daughters and edge `e`.
"""
mutable struct sTpb <: sT
  d1::sTpb
  d2::sTpb
  e ::Float64
  fx::Bool

  sTpb() = new()
  sTpb(e::Float64) = (x = new(); x.e = e; x.fx = false; x)
  sTpb(e::Float64, fx::Bool) = (x = new(); x.e = e; x.fx = fx; x)
  sTpb(d1::sTpb, d2::sTpb, e::Float64, fx::Bool) = new(d1, d2, e, fx)
end

# pretty-printing
Base.show(io::IO, t::sTpb) =
  print(io, "insane simple pure-birth tree with ", ntips(t), " tips")




"""
    sTpb(tree::sTpb)

Creates a new `sTpb` copy.
"""
function sTpb(tree::sTpb)
  if def1(tree)
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
  if def1(tree)
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
  if def1(tree)
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

    sTfbd(e::Float64)

Constructs an empty `sTfbd` object with edge `e`.

    sTfbd(d1::sTfbd, d2::sTfbd, e::Float64)

Constructs an `sTfbd` object with two `sTfbd` daughters and edge `e`.

    sTfbd(d1::sTfbd, e::Float64)

Constructs an `sTfbd` object with one sampled ancestor, one `sTfbd` daughter and
edge `e`.
"""
mutable struct sTfbd <: sT
  d1::sTfbd
  d2::sTfbd
  e ::Float64
  iμ::Bool
  iψ::Bool
  fx::Bool

  sTfbd() = new()
  sTfbd(e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
    (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x)
  sTfbd(d1::sTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
    (x = new(); x.d1 = d1; x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x)
  sTfbd(d1::sTfbd, d2::sTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
    new(d1, d2, e, iμ, iψ, fx)
end

# pretty-printing
function Base.show(io::IO, t::sTfbd)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane simple fossil tree with ",
    nt , " tip",  (isone(nt) ? "" : "s" ),
    ", (", ntipsextinct(t)," extinct) and ",
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    sTfbd(tree::sTfbd)

Creates a copy of a `sTfbd` tree.
"""
function sTfbd(tree::sTfbd)
  if def1(tree)
    if def2(tree)
      sTfbd(sTfbd(tree.d1), sTfbd(tree.d2), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree))
    else
      sTfbd(sTfbd(tree.d1), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree))
    end
  else
    sTfbd(e(tree), isextinct(tree), isfossil(tree), isfix(tree))
  end
end




"""
    sTfbd_wofe(tree::sTfbd)

Creates a copy of a `sTfbd` tree without fossils extinct tips.
"""
function sTfbd_wofe(tree::sTfbd)
  if def1(tree) && def2(tree)
      sTfbd(sTfbd_wofe(tree.d1), sTfbd_wofe(tree.d2), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree))
  else
    sTfbd(e(tree), isextinct(tree), isfossil(tree), isfix(tree))
  end
end




"""
    iT

An abstract type for all composite recursive types
representing a binary phylogenetic tree with Geometric
Brownian motion rates for `insane` use
"""
abstract type iT <: iTree end




"""
    iTpb

A composite recursive type of supertype `iT`
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

    iTpb()

Constructs an empty `iTpb` object.

    iTpb(e::Float64)

Constructs an empty `iTpb` object with pendant edge `pe`.

    iTpb(d1::iTpb, d2::iTpb, e::Float64)

Constructs an `iTpb` object with two `iTpb` daughters and pendant edge `pe`.
"""
mutable struct iTpb <: iT
  d1 ::iTpb
  d2 ::iTpb
  e  ::Float64
  dt ::Float64
  fdt::Float64
  fx ::Bool
  lλ ::Array{Float64,1}

  iTpb() = new()
  iTpb(e::Float64, dt::Float64, fdt::Float64, fx::Bool, lλ::Array{Float64,1}) =
      (x = new(); x.e = e; x.dt = dt; x.fdt = fdt; x.fx = fx; x.lλ = lλ; x)
  iTpb(d1::iTpb, d2::iTpb, e::Float64, dt::Float64, fdt::Float64, 
    fx::Bool, lλ::Array{Float64,1}) =
      new(d1, d2, e, dt, fdt, fx, lλ)
end


# pretty-printing
Base.show(io::IO, t::iTpb) =
  print(io, "insane pb-gbm tree with ", ntips(t), " tips")




"""
    iTpb(e0::Array{Int64,1},
         e1::Array{Int64,1},
         el::Array{Float64,1},
         λs::Array{Array{Float64,1},1},
         ea::Array{Int64,1},
         ni::Int64,
         ei::Int64,
         δt::Float64)

Transform edge structure to `iTpb`.
"""
function iTpb(e0::Array{Int64,1},
              e1::Array{Int64,1},
              el::Array{Float64,1},
              λs::Array{Array{Float64,1},1},
              ea::Array{Int64,1},
              ni::Int64,
              ei::Int64,
              δt::Float64)

  # if tip
  if in(ei, ea)
    return iTpb(el[ei], δt, δt, true, λs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTpb(iTpb(e0, e1, el, λs, ea, n1, ei1, δt),
                iTpb(e0, e1, el, λs, ea, n2, ei2, δt),
                el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), true, λs[ei])
  end
end




"""
    iTpb(tree::iTpb)

Produce a new copy of `iTpb`.
"""
function iTpb(tree::iTpb)
  if def1(tree)
    iTpb(iTpb(tree.d1), iTpb(tree.d2),
      e(tree), dt(tree), fdt(tree), isfix(tree), copy(lλ(tree)))
  else
    iTpb(e(tree), dt(tree), fdt(tree), isfix(tree), copy(lλ(tree)))
  end
end




"""
    iTce

A composite recursive type of supertype `iT`
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

  iTce(d1 ::iTce,
          d2 ::iTce,
          e  ::Float64,
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1})
"""
mutable struct iTce <: iT
  d1 ::iTce
  d2 ::iTce
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}

  iTce() = new()
  iTce(e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       iμ ::Bool,
       fx ::Bool,
       lλ ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x)
  iTce(d1 ::iTce,
       d2 ::iTce,
       e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       iμ ::Bool,
       fx ::Bool,
       lλ ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, iμ, fx, lλ)
end

# pretty-printing
Base.show(io::IO, t::iTce) =
  print(io, "insane gbm-ce tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    iTce(tree::iTce)

Produce a new copy of `iTce`.
"""
function iTce(tree::iTce)
  if def1(tree)
    iTce(iTce(tree.d1), iTce(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)))
  else
    iTce(e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)))
  end
end




"""
    iTce(e0::Array{Int64,1},
            e1::Array{Int64,1},
            el::Array{Float64,1},
            λs::Array{Array{Float64,1},1},
            ea::Array{Int64,1},
            ee::Array{Int64,1},
            ni::Int64,
            ei::Int64,
            δt::Float64)

Transform edge structure to `iTce`.
"""
function iTce(e0::Array{Int64,1},
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
    return iTce(el[ei], δt, δt, false, false, λs[ei])
  # if extinct
  elseif in(ei, ee)
    return iTce(el[ei], δt, δt, true, false, λs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTce(iTce(e0, e1, el, λs, ea, ee, n1, ei1, δt),
                   iTce(e0, e1, el, λs, ea, ee, n2, ei2, δt),
                   el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), false, false, λs[ei])
  end
end




"""
    iTct

A composite recursive type of supertype `iT`
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

  iTct(d1 ::iTct,
          d2 ::iTct,
          e  ::Float64,
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1})
"""
mutable struct iTct <: iT
  d1 ::iTct
  d2 ::iTct
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}

  iTct() = new()
  iTct(e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       iμ ::Bool,
       fx ::Bool,
       lλ ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x)
  iTct(d1 ::iTct,
       d2 ::iTct,
       e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       iμ ::Bool,
       fx ::Bool,
       lλ ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, iμ, fx, lλ)
end


# pretty-printing
Base.show(io::IO, t::iTct) =
  print(io, "insane gbm-ct tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    iTct(tree::iTct)

Produce a new copy of `iTct`.
"""
function iTct(tree::iTct)
  if def1(tree)
    iTct(iTct(tree.d1), iTct(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)))
  else
    iTct(e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)))
  end
end




"""
    iTct(e0::Array{Int64,1},
         e1::Array{Int64,1},
         el::Array{Float64,1},
         λs::Array{Array{Float64,1},1},
         ea::Array{Int64,1},
         ee::Array{Int64,1},
         ni::Int64,
         ei::Int64,
         δt::Float64)

Transform edge structure to `iTct`.
"""
function iTct(e0::Array{Int64,1},
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
    return iTct(el[ei], δt, δt, false, false, λs[ei])
  # if extinct
  elseif in(ei, ee)
    return iTct(el[ei], δt, δt, true, false, λs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTct(iTct(e0, e1, el, λs, ea, ee, n1, ei1, δt),
                   iTct(e0, e1, el, λs, ea, ee, n2, ei2, δt),
                   el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), false, false, λs[ei])
  end
end




"""
    iTbd

A composite recursive type of supertype `iT`
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

  iTbd(d1 ::iTbd,
          d2 ::iTbd,
          e  ::Float64,
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1},
          lμ ::Array{Float64,1})
"""
mutable struct iTbd <: iT
  d1 ::iTbd
  d2 ::iTbd
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  lμ ::Array{Float64,1}

  iTbd() = new()
  iTbd(e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       iμ ::Bool,
       fx ::Bool,
       lλ ::Array{Float64,1},
       lμ ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
  iTbd(d1 ::iTbd,
       d2 ::iTbd,
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
Base.show(io::IO, t::iTbd) =
  print(io, "insane gbm-bd tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    iTbd(tree::iTbd)

Produce a new copy of `iTbd`.
"""
function iTbd(tree::iTbd)
  if def1(tree)
    iTbd(iTbd(tree.d1), iTbd(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
  else
    iTbd(e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
  end
end




"""
    iTbd(e0::Array{Int64,1},
         e1::Array{Int64,1},
         el::Array{Float64,1},
         λs::Array{Array{Float64,1},1},
         μs::Array{Array{Float64,1},1},
         ea::Array{Int64,1},
         ee::Array{Int64,1},
         ni::Int64,
         ei::Int64,
         δt::Float64)

Transform edge structure to `iTbd`.
"""
function iTbd(e0::Array{Int64,1},
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
    return iTbd(el[ei], δt, δt, false, false, λs[ei], μs[ei])
  # if extinct
  elseif in(ei, ee)
    return iTbd(el[ei], δt, δt, true, false, λs[ei], μs[ei])
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTbd(iTbd(e0, e1, el, λs, μs, ea, ee, n1, ei1, δt),
                iTbd(e0, e1, el, λs, μs, ea, ee, n2, ei2, δt),
              el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt),
              false, false, λs[ei], μs[ei])
  end
end




"""
    iTfbd

A composite recursive type of supertype `iT`
representing a binary phylogenetic tree with  `λ` and `μ`
evolving as a Geometric Brownian motion with fossils for `insane` use,
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    pendant edge
  iμ:   if extinct node
  iψ:   if is fossil
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`
  lμ:   array of a Brownian motion evolution of `log(μ)`

  iTfbd(d1 ::iTfbd,
        d2 ::iTfbd,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        iψ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1})
"""
mutable struct iTfbd <: iT
  d1 ::iTfbd
  d2 ::iTfbd
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  iψ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  lμ ::Array{Float64,1}


  iTfbd() = new()
  iTfbd(e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        iψ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
  iTfbd(d1 ::iTfbd,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        iψ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1}) =
    (x = new(); x.d1 = d1; x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
  iTfbd(d1 ::iTfbd,
        d2 ::iTfbd,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        iψ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, iμ, iψ, fx, lλ, lμ)
end


# pretty-printing
function Base.show(io::IO, t::iTfbd)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane gbm-bd fossil tree with ",
    nt , " tip",  (isone(nt) ? "" : "s" ),
    ", (", ntipsextinct(t)," extinct) and ",
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    iTfbd(tree::iTfbd)

Produce a new copy of `iTfbd`.
"""
function iTfbd(tree::iTfbd)
  if def1(tree)
    if def2(tree)
      iTfbd(iTfbd(tree.d1), iTfbd(tree.d2),
        e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
        isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
    else
      iTfbd(iTfbd(tree.d1),
        e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
        isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
    end
  else
    iTfbd(e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
  end
end




"""
    iTfbd_wofe(tree::iTfbd)

Creates a copy of a `iTfbd` tree without fossils extinct tips.
"""
function iTfbd_wofe(tree::iTfbd)
  if def1(tree) && def2(tree)
    iTfbd(iTfbd_wofe(tree.d1), iTfbd_wofe(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
  else
    iTfbd(e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
  end
end




"""
    iTfbd(e0::Array{Int64,1},
          e1::Array{Int64,1},
          el::Array{Float64,1},
          λs::Array{Array{Float64,1},1},
          μs::Array{Array{Float64,1},1},
          ea::Array{Int64,1},
          ee::Array{Int64,1},
          ef::Array{Int64,1},
          ni::Int64,
          ei::Int64,
          δt::Float64)

Transform edge structure to `iTfbd`.
"""
function iTfbd(e0::Array{Int64,1},
               e1::Array{Int64,1},
               el::Array{Float64,1},
               λs::Array{Array{Float64,1},1},
               μs::Array{Array{Float64,1},1},
               ea::Array{Int64,1},
               ee::Array{Int64,1},
               ef::Array{Int64,1},
               ni::Int64,
               ei::Int64,
               δt::Float64)

  # if tip
  if in(ei, ea)
    return iTfbd(el[ei], δt, δt, false, false, false, λs[ei], μs[ei])
  
  # if extinct
  elseif in(ei, ee)
    return iTfbd(el[ei], δt, δt, true, false, false, λs[ei], μs[ei])

  # if fossil
  elseif in(ei, ef)
    ei1 = findfirst(isequal(ni), e0)
    n1  = e1[ei1]
    return iTfbd(iTfbd(e0, e1, el, λs, μs, ea, ee, ef, n1, ei1, δt),
                 el[ei], δt, δt, false, true, false, λs[ei], μs[ei])

  # if internal
  else
    ei1, ei2 = findall(isequal(ni), e0)
    n1, n2   = e1[ei1:ei2]
    return iTfbd(iTfbd(e0, e1, el, λs, μs, ea, ee, ef, n1, ei1, δt),
                 iTfbd(e0, e1, el, λs, μs, ea, ee, ef, n2, ei2, δt),
                 el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt),
                 false, false, false, λs[ei], μs[ei])
  end
end



