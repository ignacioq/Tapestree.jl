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

    sTpbX(e::Float64, fx::Bool, xi::Float64, xf::Float64)

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
  if isdefined(tree, :d1)
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
mutable struct sTbdX <: sTX
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
  if isdefined(tree, :d1)
    sTbdX(sTbdX(tree.d1), sTbdX(tree.d2), 
      e(tree), isextinct(tree), isfix(tree), xi(tree), xf(tree))
  else
    sTbdX(e(tree), isextinct(tree), isfix(tree), xi(tree), xf(tree))
  end
end


