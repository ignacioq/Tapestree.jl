#=

punkeek tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 20 11 2024
=#


"""
    sTpe

The composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `punkeek` use,
with the following fields:

  d1::sTpe
  d2::sTpe
  e ::Float64
  iμ::Bool
  xi::Float64
  xf::Float64
  sh::Bool
  fx::Bool

Constructs a tip.

  sTpe(e::Float64, iμ::Bool, xi::Float64, xf::Float64, sh::Bool, fx::Bool) =
    (x = new(); x.e = e; x.iμ = iμ; x.xi = xi; x.xf = xf; x.sh = sh;  x.fx = fx; x)

Constructs an `sTpe` object with two `sTpe` daughters and edge `e`.

  sTpe(d1::sTpe, d2::sTpe, e::Float64, iμ::Bool, xi::Float64, xf::Float64, sh::Bool, fx::Bool) =
    new(d1, d2, e, iμ, xi, xf, sh, fx)
"""
mutable struct sTpe <: sT
  d1::sTpe
  d2::sTpe
  e ::Float64
  iμ::Bool
  xi::Float64
  xf::Float64
  sh::Bool
  fx::Bool

  sTpe() = new()
  sTpe(e::Float64, iμ::Bool, xi::Float64, xf::Float64, sh::Bool, fx::Bool) =
    (x = new(); x.e = e; x.iμ = iμ; x.xi = xi; x.xf = xf; x.sh = sh; x.fx = fx; x)
  sTpe(d1::sTpe, d2::sTpe, e::Float64, iμ::Bool, xi::Float64, xf::Float64, sh::Bool, fx::Bool) =
    new(d1, d2, e, iμ, xi, xf, sh, fx)
end

# pretty-printing
Base.show(io::IO, t::sTpe) =
  print(io, "insane simple punkeek tree with ", ntips(t), " tips (",
    ntipsextinct(t)," extinct)")




"""
    sTpe(tree::sTpe)

Produce a new copy of `sTpe`.
"""
function sTpe(tree::sTpe)
  if def1(tree)
    sTpe(sTpe(tree.d1), sTpe(tree.d2), e(tree), isextinct(tree), 
      xi(tree), xf(tree), sh(tree), isfix(tree))
  else
    sTpe(e(tree), isextinct(tree), xi(tree), xf(tree), false, isfix(tree))
  end
end


# """
#     sTpe(tree::sT_label)

# Transform a tree of type `sT_label` to `sTpe`.
# """
# function sTpe(tree::sT_label)
#   if def1(tree)
#     sTpe(sTpe(tree.d1), sTpe(tree.d2), e(tree), false, false)
#   else
#     sTpe(e(tree), false)
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
