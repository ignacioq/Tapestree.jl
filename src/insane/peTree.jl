#=

punkeek tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 20 11 2024
=#


"""
    peT

The composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `punkeek` use,
with the following fields:

  d1::peT
  d2::peT
  e ::Float64
  iμ::Bool
  xi::Float64
  xf::Float64
  sh::Bool
  fx::Bool

Constructs a tip.

  peT(e::Float64, iμ::Bool, xi::Float64, xf::Float64, sh::Bool, fx::Bool) =
    (x = new(); x.e = e; x.iμ = iμ; x.xi = xi; x.xf = xf; x.sh = sh;  x.fx = fx; x)

Constructs an `peT` object with two `peT` daughters and edge `e`.

  peT(d1::peT, d2::peT, e::Float64, iμ::Bool, xi::Float64, xf::Float64, sh::Bool, fx::Bool) =
    new(d1, d2, e, iμ, xi, xf, sh, fx)
"""
mutable struct peT <: sT
  d1::peT
  d2::peT
  e ::Float64
  iμ::Bool
  xi::Float64
  xf::Float64
  sh::Bool
  fx::Bool

  peT() = new()
  peT(e::Float64, iμ::Bool, xi::Float64, xf::Float64, fx::Bool) =
    (x = new(); x.e = e; x.iμ = iμ; x.xi = xi; x.xf = xf; x.sh = false; x.fx = fx; x)
  peT(d1::peT, d2::peT, e::Float64, iμ::Bool, xi::Float64, xf::Float64, sh::Bool, fx::Bool) =
    new(d1, d2, e, iμ, xi, xf, sh, fx)
end

# pretty-printing
Base.show(io::IO, t::peT) =
  print(io, "insane simple punkeek tree with ", ntips(t), " tips (",
    ntipsextinct(t)," extinct)")



# """
#     peT(tree::sT_label)

# Transform a tree of type `sT_label` to `peT`.
# """
# function peT(tree::sT_label)
#   if def1(tree)
#     peT(peT(tree.d1), peT(tree.d2), e(tree), false, false)
#   else
#     peT(e(tree), false)
#   end
# end




# """
#     peT(tree::peT)

# Produce a new copy of `peT`.
# """
# function peT(tree::peT)
#   if def1(tree)
#     peT(peT(tree.d1), peT(tree.d2), e(tree), isextinct(tree), isfix(tree))
#   else
#     peT(e(tree), isextinct(tree), isfix(tree))
#   end
# end




# """
#     sTfbd

# The simplest composite recursive type of supertype `sTf`
# representing a binary phylogenetic tree for `insane` use,
# with the following fields:

#   d1: daughter tree 1
#   d2: daughter tree 2
#   e:  edge
#   iμ: is an extinction node
#   iψ: is a fossil node
#   fx: if it is fix

#     sTfbd(e::Float64)

# Constructs an empty `sTfbd` object with edge `e`.

#     sTfbd(d1::sTfbd, d2::sTfbd, e::Float64)

# Constructs an `sTfbd` object with two `sTfbd` daughters and edge `e`.

#     sTfbd(d1::sTfbd, e::Float64)

# Constructs an `sTfbd` object with one sampled ancestor, one `sTfbd` daughter and
# edge `e`.
# """
# mutable struct sTfbd <: sT
#   d1::sTfbd
#   d2::sTfbd
#   e ::Float64
#   iμ::Bool
#   iψ::Bool
#   fx::Bool

#   sTfbd() = new()
#   sTfbd(e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
#     (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x)
#   sTfbd(d1::sTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
#     (x = new(); x.d1 = d1; x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x)
#   sTfbd(d1::sTfbd, d2::sTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
#     new(d1, d2, e, iμ, iψ, fx)
# end

# # pretty-printing
# function Base.show(io::IO, t::sTfbd)
#   nt = ntips(t)
#   nf = nfossils(t)

#   print(io, "insane simple fossil tree with ",
#     nt , " tip",  (isone(nt) ? "" : "s" ),
#     ", (", ntipsextinct(t)," extinct) and ",
#     nf," fossil", (isone(nf) ? "" : "s" ))
# end




# """
#     sTfbd(tree::sTfbd)

# Creates a copy of a `sTfbd` tree.
# """
# function sTfbd(tree::sTfbd)
#   if def1(tree)
#     if def2(tree)
#       sTfbd(sTfbd(tree.d1), sTfbd(tree.d2), e(tree),
#         isextinct(tree), isfossil(tree), isfix(tree))
#     else
#       sTfbd(sTfbd(tree.d1), e(tree),
#         isextinct(tree), isfossil(tree), isfix(tree))
#     end
#   else
#     sTfbd(e(tree), isextinct(tree), isfossil(tree), isfix(tree))
#   end
# end
