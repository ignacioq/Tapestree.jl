#=

Geometric Brownian motion insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    iTgbm

A simple composite recursive type of supertype `iTree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node

    iTgbm()

Constructs an empty `iTgbm` object.

    iTgbm(pe::Float64)

Constructs an empty `iTgbm` object with pendant edge `pe`.

    iTgbm(d1::iTgbm, d2::iTgbm, pe::Float64)

Constructs an `iTgbm` object with two `iTgbm` daughters and pendant edge `pe`.
"""
mutable struct iTgbm <: iTree
  d1::Union{iTgbm, Nothing}
  d2::Union{iTgbm, Nothing}
  pe::Float64
  iμ::Bool
  fx::Bool

  # inner constructor
  iTgbm(d1::Union{iTgbm, Nothing}, d2::Union{iTgbm, Nothing}, 
    pe::Float64, iμ::Bool, fx::Bool) = new(d1, d2, pe, iμ, fx)
end

# outer constructors
iTgbm() = iTgbm(nothing, nothing, 0.0, false, false)

iTgbm(pe::Float64) = iTgbm(nothing, nothing, pe, false, false)

iTgbm(pe::Float64, iμ::Bool) = iTgbm(nothing, nothing, pe, iμ, false)

iTgbm(d1::iTgbm, d2::iTgbm, pe::Float64) = iTgbm(d1, d2, pe, false, false)

iTgbm(d1::iTgbm, d2::iTgbm, pe::Float64, iμ::Bool) = 
  iTgbm(d1, d2, pe, iμ, false)

# pretty-printing
Base.show(io::IO, t::iTgbm) = 
  print(io, "insane tree with ", sntn(t), " tips (", snen(t)," extinct)")
