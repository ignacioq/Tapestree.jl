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
mutable struct sTfbdX <: sTf
  d1::sTfbdX
  d2::sTfbdX
  e ::Float64
  iμ::Bool
  iψ::Bool
  fx::Bool
  xi::Float64
  xf::Float64

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
Union type for trait data
"""
sTX = Union{sTpbX, sTbdX, sTfbdX}

