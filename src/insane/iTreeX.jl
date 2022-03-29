#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




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

    sTpbX(e::Float64, fx::Bool, xi::Float64, xf::Float64)

Constructs an `sTpbX` object with two `sTpbX` daughters and edge `e`,
fix information `fx`, initial node trait `xi` and final `xf`.
"""
mutable struct sTpbX <: sT
  d1::sTpbX
  d2::sTpbX
  e ::Float64
  fx::Bool
  xi::Float64
  xf::Float64

  sTpbX() = new()
  sTpbX(e::Float64, fx::Bool, xi::Float64, xf::Float64) =
    (t = new(); t.e = e; t.fx = fx; t.xi = xi; t.xf = xf; t)
  sTpbX(d1::sTpbX, d2::sTpbX, e::Float64, fx::Bool, xi::Float64, xf::Float64) =
    new(d1, d2, e, fx, xi, xf)
end

# pretty-printing
Base.show(io::IO, t::sTpbX) =
  print(io, "insane pure-birth with trait tree with ", ntips(t), " tips")




"""
    sTpbX(tree::sTpbX)

Creates a copy of `sTpbX`.
"""
function sTpbX(tree::sTpbX)
  if def1(tree)
    sTpbX(sTpbX(tree.d1), sTpbX(tree.d2),
      e(tree), isfix(tree), xi(tree), xf(tree))
  else
    sTpbX(e(tree), isfix(tree), xi(tree), xf(tree))
  end
end




"""
    sTbdX

The simplest composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: if extinct
  fx: if fix
  xi: initial trait value
  xf: final trait value


    sTbdX(e::Float64, iμ::Bool, fx::Bool, xi::Float64, xf::Float64)

Constructs an `sTbdX` object with two `sTbdX` daughters and edge `e`,
fix information `fx`, initial node trait `xi` and final `xf`.
"""
mutable struct sTbdX <: sT
  d1::sTbdX
  d2::sTbdX
  e ::Float64
  iμ::Bool
  fx::Bool
  xi::Float64
  xf::Float64

  sTbdX() = new()
  sTbdX(e::Float64, iμ::Bool, fx::Bool,  xi::Float64, xf::Float64) =
    (t = new(); t.e = e; t.iμ = iμ;t.fx = fx; t.xi = xi; t.xf = xf; t)
  sTbdX(d1::sTbdX, d2::sTbdX, e::Float64, iμ::Bool, fx::Bool,
    xi::Float64, xf::Float64) =
      new(d1, d2, e, iμ, fx, xi, xf)
end

# pretty-printing
Base.show(io::IO, t::sTbdX) =
  print(io, "insane birth-death tree with trait with ", ntips(t), " tips (",
    ntipsextinct(t)," extinct)")




"""
    sTbdX(tree::sTbdX)

Creates a copy of `sTbdX`.
"""
function sTbdX(tree::sTbdX)
  if def1(tree)
    sTbdX(sTbdX(tree.d1), sTbdX(tree.d2),
      e(tree), isextinct(tree), isfix(tree), xi(tree), xf(tree))
  else
    sTbdX(e(tree), isextinct(tree), isfix(tree), xi(tree), xf(tree))
  end
end




"""
    sTfbdX

The simplest composite recursive type of supertype `sTf`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: is an extinction node
  iψ: is a fossil node
  fx: if it is fix
  xi: initial trait value
  xf: final trait value

    sTfbdXX(e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

    sTfbdXX(d1::sTfbdX, e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

    sTfbdXX(d1::sTfbdX, d2::sTfbdX, e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

Constructs an `sTfbdX` object with one sampled ancestor, one `sTfbdX` daughter and
edge `e`.
"""
mutable struct sTfbdX <: sT
  d1::sTfbdX
  d2::sTfbdX
  e ::Float64
  iμ::Bool
  iψ::Bool
  fx::Bool
  xi::Float64
  xf::Float64

  sTfbdX() = new()
  sTfbdX(e::Float64, iμ::Bool, iψ::Bool, fx::Bool, xi::Float64, xf::Float64) =
    (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.xi = xi; x.xf = xf;
      x)
  sTfbdX(d1::sTfbdX, e::Float64, iμ::Bool, iψ::Bool, fx::Bool,
    xi::Float64, xf::Float64) =
    (x = new(); x.d1 = d1; x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx;
      x.xi = xi; x.xf = xf; x)
  sTfbdX(d1::sTfbdX, d2::sTfbdX, e::Float64, iμ::Bool, iψ::Bool, fx::Bool,
    xi::Float64, xf::Float64) = new(d1, d2, e, iμ, iψ, fx, xi, xf)
end

# pretty-printing
function Base.show(io::IO, t::sTfbdX)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane simple fossil tree with traits with ",
    nt , " tip",  (isone(nt) ? "" : "s" ),
    ", (", ntipsextinct(t)," extinct) and ",
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    sTfbdX(tree::sTfbdX)

Produces a copy of `sTfbdX`.
"""
function sTfbdX(tree::sTfbdX)
  if def1(tree)
    d1 = sTfbdX(tree.d1)
    if def2(tree)
      sTfbdX(d1, sTfbdX(tree.d2), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree), xi(tree), xf(tree))
    else
      sTfbdX(d1, e(tree),
        isextinct(tree), isfossil(tree), isfix(tree), xi(tree), xf(tree))
    end
  else
    sTfbdX(e(tree), isextinct(tree), isfossil(tree), isfix(tree),
      xi(tree), xf(tree))
  end
end





"""
    iTpbX

A composite recursive type of supertype `iT`
representing a binary phylogenetic tree with  `λ` evolving as a
Geometric Brownian motion and no extinction,
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    edge
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`
  xv:   array of a Brownian motion evolution of `X`.

  iTpbX(d1 ::iTpbX,
          d2 ::iTpbX,
          e  ::Float64,
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1})
"""
mutable struct iTpbX <: iT
  d1 ::iTpbX
  d2 ::iTpbX
  e  ::Float64
  dt ::Float64
  fdt::Float64
  fx ::Bool
  lλ ::Array{Float64,1}
  xv ::Array{Float64,1}

  iTpbX() = new()
  iTpbX(e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        fx ::Bool,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.fx = fx; x.lλ = lλ; x.xv = xv; x)
  iTpbX(d1 ::iTpbX,
        d2 ::iTpbX,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        fx ::Bool,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, fx, lλ, xv)
end

# pretty-printing
Base.show(io::IO, t::iTpbX) =
  print(io, "insane gbm-pb tree with traits with", ntips(t), 
    " tips (", ntipsextinct(t)," extinct)")




"""
    iTpbX(tree::iTpbX)

Produce a new copy of `iTpbX`.
"""
function iTpbX(tree::iTpbX)
  if def1(tree)
    iTpbX(iTpbX(tree.d1), iTpbX(tree.d2),
      e(tree), dt(tree), fdt(tree),
      isfix(tree), copy(lλ(tree)), copy(xv(tree)))
  else
    iTpbX(e(tree), dt(tree), fdt(tree),
      isfix(tree), copy(lλ(tree)), copy(xv(tree)))
  end
end





"""
    iTceX

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
  xv:   array of a Brownian motion evolution of `X`.

  iTceX(d1 ::iTceX,
          d2 ::iTceX,
          e  ::Float64,
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1})
"""
mutable struct iTceX <: iT
  d1 ::iTceX
  d2 ::iTceX
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  xv ::Array{Float64,1}

  iTceX() = new()
  iTceX(e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.xv = xv; x)
  iTceX(d1 ::iTceX,
        d2 ::iTceX,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, iμ, fx, lλ, xv)
end

# pretty-printing
Base.show(io::IO, t::iTceX) =
  print(io, "insane gbm-ce tree with traits with ", ntips(t), 
    " tips (", ntipsextinct(t)," extinct)")




"""
    iTceX(tree::iTceX)

Produce a new copy of `iTceX`.
"""
function iTceX(tree::iTceX)
  if def1(tree)
    iTceX(iTceX(tree.d1), iTceX(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(xv(tree)))
  else
    iTceX(e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(xv(tree)))
  end
end





#=
Union Types
=#




"""
    Union type for label trees

Tlabel = Union{sT_label, sTf_label}
"""
Tlabel = Union{sT_label, sTf_label}




"""
    Union type for trait data

sTX = Union{sTpbX, sTbdX, sTfbdX}
"""
sTX = Union{sTpbX, sTbdX, sTfbdX}



"""
    Union type for trait data

iTX = Union{iTpbX, iTceX}
"""
iTX = Union{iTpbX, iTceX}




"""
    Union type for fossil data

iTf = Union{sTfbd, sTfbdX, iTfbd}
"""
iTf = Union{sTf_label, sTfbd, sTfbdX, iTfbd}



"""
    Union type for gbm-bd data

iTbdU = Union{iTbd, iTfbd}
"""
iTbdU = Union{iTbd, iTfbd}






#=
Type aliases (for compatibility with older versions)
=#


const iTgbmpb = iTpb

const iTgbmce = iTce

const iTgbmct = iTct

const iTgbmbd = iTbd
