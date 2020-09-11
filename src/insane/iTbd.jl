#=

insane tree structure for birth-death

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    iTbd

The simplest composite recursive type of supertype `iTbdree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node

    iTbd()

Constructs an empty `iTbd` object.

    iTbd(pe::Float64)

Constructs an empty `iTbd` object with pendant edge `pe`.

    iTbd(d1::iTbd, d2::iTbd, pe::Float64)

Constructs an `iTbd` object with two `iTbd` daughters and pendant edge `pe`.
"""
mutable struct iTbd <: iT
  d1::Union{iTbd, Nothing}
  d2::Union{iTbd, Nothing}
  pe::Float64
  iμ::Bool
  fx::Bool

  # inner constructor
  iTbd(d1::Union{iTbd, Nothing}, d2::Union{iTbd, Nothing}, 
    pe::Float64, iμ::Bool, fx::Bool) = new(d1, d2, pe, iμ, fx)
end

# outer constructors
iTbd() = iTbd(nothing, nothing, 0.0, false, false)

iTbd(pe::Float64) = iTbd(nothing, nothing, pe, false, false)

iTbd(pe::Float64, iμ::Bool) = iTbd(nothing, nothing, pe, iμ, false)

iTbd(d1::iTbd, d2::iTbd, pe::Float64) = iTbd(d1, d2, pe, false, false)

iTbd(d1::iTbd, d2::iTbd, pe::Float64, iμ::Bool) = 
  iTbd(d1, d2, pe, iμ, false)

# pretty-printing
Base.show(io::IO, t::iTbd) = 
  print(io, "insane simple birth-death tree with ", sntn(t), " tips (", snen(t)," extinct)")

  