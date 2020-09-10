#=

insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    iT

The simplest composite recursive type of supertype `iTree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node

    iT()

Constructs an empty `iT` object.

    iT(pe::Float64)

Constructs an empty `iT` object with pendant edge `pe`.

    iT(d1::iT, d2::iT, pe::Float64)

Constructs an `iT` object with two `iT` daughters and pendant edge `pe`.
"""
mutable struct iT <: iTree
  d1::Union{iT, Nothing}
  d2::Union{iT, Nothing}
  pe::Float64
  iμ::Bool
  fx::Bool

  # inner constructor
  iT(d1::Union{iT, Nothing}, d2::Union{iT, Nothing}, 
    pe::Float64, iμ::Bool, fx::Bool) = new(d1, d2, pe, iμ, fx)
end

# outer constructors
iT() = iT(nothing, nothing, 0.0, false, false)

iT(pe::Float64) = iT(nothing, nothing, pe, false, false)

iT(pe::Float64, iμ::Bool) = iT(nothing, nothing, pe, iμ, false)

iT(d1::iT, d2::iT, pe::Float64) = iT(d1, d2, pe, false, false)

iT(d1::iT, d2::iT, pe::Float64, iμ::Bool) = 
  iT(d1, d2, pe, iμ, false)

# pretty-printing
Base.show(io::IO, t::iT) = 
  print(io, "insane simple tree with ", sntn(t), " tips (", snen(t)," extinct)")

  