#=

insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    iTsimple

The simplest composite recursive type of supertype `iTree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node

    iTsimple()

Constructs an empty `iTsimple` object.

    iTsimple(pe::Float64)

Constructs an empty `iTsimple` object with pendant edge `pe`.

    iTsimple(d1::iTsimple, d2::iTsimple, pe::Float64)

Constructs an `iTsimple` object with two `iTsimple` daughters and pendant edge `pe`.
"""
mutable struct iTsimple <: iTree
  d1::Union{iTsimple, Nothing}
  d2::Union{iTsimple, Nothing}
  pe::Float64
  iμ::Bool
  fx::Bool

  # inner constructor
  iTsimple(d1::Union{iTsimple, Nothing}, d2::Union{iTsimple, Nothing}, 
    pe::Float64, iμ::Bool, fx::Bool) = new(d1, d2, pe, iμ, fx)
end

# outer constructors
iTsimple() = iTsimple(nothing, nothing, 0.0, false, false)

iTsimple(pe::Float64) = iTsimple(nothing, nothing, pe, false, false)

iTsimple(pe::Float64, iμ::Bool) = iTsimple(nothing, nothing, pe, iμ, false)

iTsimple(d1::iTsimple, d2::iTsimple, pe::Float64) = iTsimple(d1, d2, pe, false, false)

iTsimple(d1::iTsimple, d2::iTsimple, pe::Float64, iμ::Bool) = 
  iTsimple(d1, d2, pe, iμ, false)

# pretty-printing
Base.show(io::IO, t::iTsimple) = 
  print(io, "insane tree with ", sntn(t), " tips (", snen(t)," extinct)")

  