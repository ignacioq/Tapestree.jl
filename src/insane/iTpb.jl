#=

insane tree structure for pure-birth

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#




"""
    iTpb

The simplest composite recursive type of supertype `iTpbree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node

    iTpb()

Constructs an empty `iTpb` object.

    iTpb(pe::Float64)

Constructs an empty `iTpb` object with pendant edge `pe`.

    iTpb(d1::iTpb, d2::iTpb, pe::Float64)

Constructs an `iTpb` object with two `iTpb` daughters and pendant edge `pe`.
"""
mutable struct iTpb <: iT
  d1::Union{iTpb, Nothing}
  d2::Union{iTpb, Nothing}
  pe::Float64

  # inner constructor
  iTpb(d1::Union{iTpb, Nothing}, d2::Union{iTpb, Nothing}, pe::Float64) = 
    new(d1, d2, pe)
end

# outer constructors
iTpb() = iTpb(nothing, nothing, 0.0)

iTpb(pe::Float64) = iTpb(nothing, nothing, pe)


# pretty-printing
Base.show(io::IO, t::iTpb) = 
  print(io, "insane simple pure-birth tree with ", sntn(t), " tips")

  