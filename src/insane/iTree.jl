#=

insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#



"""
    iTree

A Composite recursive type representing a binary phylogenetic tree 
for `insane` use, with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  at: absolute node time
  iμ: is an extinction node

    iTree()

Constructs an empty `iTree` object.

    iTree(pe::Float64)

Constructs an empty `iTree` object with pendant edge `pe`.

    iTree(d1::iTree, d2::iTree, pe::Float64)

Constructs an `iTree` object with two `iTree` daughters and pendant edge `pe`.

    iTree(d1::iTree, d2::iTree, pe::Float64, at::Float64)

Constructs an `iTree` object with two `iTree` daughters, pendant edge `pe` at absolute time `at`.
"""
mutable struct iTree
  d1::Union{iTree, Nothing}
  d2::Union{iTree, Nothing}
  pe::Float64
  at::Float64
  iμ::Bool

  # constructors
  iTree() = new(nothing, nothing, 0.0, 0.0, false)
  iTree(pe::Float64) = new(nothing, nothing, pe, 0.0, false)
  iTree(pe::Float64, iμ::Bool) = 
    new(nothing, nothing, pe, 0.0, iμ)
  iTree(d1::iTree, d2::iTree, pe::Float64) = 
    new(d1, d2, pe, 0.0, false)
  iTree(d1::iTree, d2::iTree, pe::Float64, iμ::Bool) = 
    new(d1, d2, pe, 0.0, iμ)
  iTree(d1::iTree, d2::iTree, pe::Float64, at::Float64) = 
    new(d1, d2, pe, at, false)
  iTree(d1::iTree, d2::iTree, pe::Float64, at::Float64, 
        isλ::Bool, iμ::Bool, ist::Bool) = 
    new(d1, d2, pe, at, iμ)
end



# pretty-printing
Base.show(io::IO, t::iTree) = 
  print(io, "insane tree with ", sntn(t), " tips (", snen(t)," extinct)")





