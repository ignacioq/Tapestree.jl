#=

Abstract insane cT structure

Ignacio Quintero Mächler

t(-_-t)

Created 16 07 2025
=#




"""
    cT

An abstract type for all composite recursive types
representing a binary phylogenetic tree with cladogenetic rate shifts 
for `insane` use.
"""
abstract type cT <: iTree end




"""
    cTb

A composite recursive type of supertype `cT`
representing a binary phylogenetic tree with no extinction
and `λ` under cladogenetic rate shifts  for `insane` use,
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   edge
  fx:  if fix tree
  lλ:  `log(λ)`

    cTb()

Constructs an empty `cTb` object.

    cTb(e::Float64, fx::Bool, lλ::Float64)

Constructs an empty `cTb` object with pendant edge `pe`.

    cTb(d1::cTb, d2::cTb, e::Float64, fx::Bool, lλ::Float64)

Constructs an `cTb` object with two `cTb` daughters and pendant edge `pe`.
"""
mutable struct cTb <: cT
  d1 ::cTb
  d2 ::cTb
  e  ::Float64
  fx ::Bool
  lλ ::Float64

  cTb() = new()
  cTb(e::Float64, fx::Bool, lλ::Float64) =
      (x = new(); x.e = e; x.fx = fx; x.lλ = lλ; x)
  cTb(d1::cTb, d2::cTb, e::Float64, fx::Bool, lλ::Float64) =
      new(d1, d2, e, fx, lλ)
end


# pretty-printing
Base.show(io::IO, t::cTb) =
  print(io, "insane pb-clads tree with ", ntips(t), " tips")



# """
#     cTb(e0::Array{Int64,1},
#          e1::Array{Int64,1},
#          el::Array{Float64,1},
#          λs::Array{Array{Float64,1},1},
#          ea::Array{Int64,1},
#          ni::Int64,
#          ei::Int64,
#          δt::Float64)

# Transform edge structure to `cTb`.
# """
# function cTb(e0::Array{Int64,1},
#               e1::Array{Int64,1},
#               el::Array{Float64,1},
#               λs::Array{Array{Float64,1},1},
#               ea::Array{Int64,1},
#               ni::Int64,
#               ei::Int64,
#               δt::Float64)

#   # if tip
#   if in(ei, ea)
#     return cTb(el[ei], δt, δt, true, λs[ei])
#   else
#     ei1, ei2 = findall(isequal(ni), e0)
#     n1, n2   = e1[ei1:ei2]
#     return cTb(cTb(e0, e1, el, λs, ea, n1, ei1, δt),
#                 cTb(e0, e1, el, λs, ea, n2, ei2, δt),
#                 el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), true, λs[ei])
#   end
# end




"""
    cTb(tree::cTb)

Produce a new copy of `cTb`.
"""
function cTb(tree::cTb)
  if def1(tree)
    cTb(cTb(tree.d1), cTb(tree.d2), e(tree), isfix(tree), lλ(tree))
  else
    cTb(e(tree), isfix(tree), lλ(tree))
  end
end




"""
    cTce

A composite recursive type of supertype `cT`
representing a binary phylogenetic tree with constant extinction
and `λ` under cladogenetic rate shifts for `insane` use,
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   edge
  iμ:  if extinct node
  fx:  if fix tree
  lλ:  `log(λ)`

    cTce()

Constructs an empty `cTce` object.

    cTce(e::Float64, fx::Bool, lλ::Float64)

Constructs an empty `cTce` object with pendant edge `pe`.

    cTce(d1::cTce, d2::cTce, e::Float64, fx::Bool, lλ::Float64)

Constructs an `cTce` object with two `cTce` daughters and pendant edge `pe`.
"""
mutable struct cTce <: cT
  d1 ::cTce
  d2 ::cTce
  e  ::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Float64

  cTce() = new()
  cTce(e::Float64, iμ::Bool, fx::Bool, lλ::Float64) =
      (x = new(); x.e = e; x.iμ = iμ; x.fx = fx; x.lλ = lλ; x)
  cTce(d1::cTce, d2::cTce, e::Float64, iμ::Bool, fx::Bool, lλ::Float64) =
      new(d1, d2, e, iμ, fx, lλ)
end


# pretty-printing
Base.show(io::IO, t::cTce) =
  print(io, "insane ce-clads tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    cTce(tree::cTce)

Produce a new copy of `cTce`.
"""
function cTce(tree::cTce)
  if def1(tree)
    cTce(cTce(tree.d1), cTce(tree.d2), e(tree), isextinct(tree), 
      isfix(tree), lλ(tree))
  else
    cTce(e(tree), isextinct(tree), isfix(tree), lλ(tree))
  end
end




"""
    cTct

A composite recursive type of supertype `cT`
representing a binary phylogenetic tree with constant turnover
and `λ` under cladogenetic rate shifts for `insane` use,
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   edge
  iμ:  if extinct node
  fx:  if fix tree
  lλ:  `log(λ)`

    cTct()

Constructs an empty `cTct` object.

    cTct(e::Float64, fx::Bool, lλ::Float64)

Constructs an empty `cTct` object with pendant edge `pe`.

    cTct(d1::cTct, d2::cTct, e::Float64, fx::Bool, lλ::Float64)

Constructs an `cTct` object with two `cTct` daughters and pendant edge `pe`.
"""
mutable struct cTct <: cT
  d1 ::cTct
  d2 ::cTct
  e  ::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Float64

  cTct() = new()
  cTct(e::Float64, iμ::Bool, fx::Bool, lλ::Float64) =
      (x = new(); x.e = e; x.iμ = iμ; x.fx = fx; x.lλ = lλ; x)
  cTct(d1::cTct, d2::cTct, e::Float64, iμ::Bool, fx::Bool, lλ::Float64) =
      new(d1, d2, e, iμ, fx, lλ)
end


# pretty-printing
Base.show(io::IO, t::cTct) =
  print(io, "insane ct-clads tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    cTct(tree::cTct)

Produce a new copy of `cTct`.
"""
function cTct(tree::cTct)
  if def1(tree)
    cTct(cTct(tree.d1), cTct(tree.d2), e(tree), isextinct(tree), 
      isfix(tree), lλ(tree))
  else
    cTct(e(tree), isextinct(tree), isfix(tree), lλ(tree))
  end
end




"""
    cTbd

A composite recursive type of supertype `cT`
representing a binary phylogenetic tree with constant turnover
and `λ` under cladogenetic rate shifts for `insane` use,
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   edge
  iμ:  if extinct node
  fx:  if fix tree
  lλ:  `log(λ)`
  lμ:  `log(μ)`

    cTbd()

Constructs an empty `cTbd` object.

    cTbd(e::Float64, iμ::Bool, fx::Bool, lλ::Float64, lμ::Float64)

Constructs an empty `cTbd` object with pendant edge `pe`.

    cTbd(d1::cTbd, d2::cTbd, e::Float64, iμ::Bool, fx::Bool, lλ::Float64, lμ::Float64)

Constructs an `cTbd` object with two `cTbd` daughters and pendant edge `pe`.
"""
mutable struct cTbd <: cT
  d1 ::cTbd
  d2 ::cTbd
  e  ::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Float64
  lμ ::Float64

  cTbd() = new()
  cTbd(e::Float64, iμ::Bool, fx::Bool, lλ::Float64, lμ::Float64) =
      (x = new(); x.e = e; x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
  cTbd(d1::cTbd, d2::cTbd, e::Float64, iμ::Bool, fx::Bool, lλ::Float64, lμ::Float64) =
      new(d1, d2, e, iμ, fx, lλ, lμ)
end


# pretty-printing
Base.show(io::IO, t::cTbd) =
  print(io, "insane bd-clads tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




"""
    cTbd(tree::cTbd)

Produce a new copy of `cTbd`.
"""
function cTbd(tree::cTbd)
  if def1(tree)
    cTbd(cTbd(tree.d1), cTbd(tree.d2), e(tree), isextinct(tree), 
      isfix(tree), lλ(tree), lμ(tree))
  else
    cTbd(e(tree), isextinct(tree), isfix(tree), lλ(tree), lμ(tree))
  end
end




"""
    cTfbd

A composite recursive type of supertype `cT`
representing a binary phylogenetic tree with constant turnover
and `λ` under cladogenetic rate shifts for `insane` use,
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   edge
  iψ:  if fossil
  iμ:  if extinct node
  fx:  if fix tree
  lλ:  `log(λ)`
  lμ:  `log(μ)`


Constructs an empty `cTfbd` object.

    cTfbd()

Constructs an empty `cTfbd` object with pendant edge `e`.

    cTfbd(e::Float64, iμ::Bool, iψ::Bool, fx::Bool, lλ::Float64, lμ::Float64)

Constructs an `cTfbd` object with one `cTfbd` daughter and pendant edge `e`.

    cTfbd(d1::cTfbd, e::Float64, iψ::Bool, iμ::Bool, fx::Bool, lλ::Float64, lμ::Float64)

Constructs an `cTfbd` object with two `cTfbd` daughters and pendant edge `e`.

    cTfbd(d1::cTfbd, e::Float64, iψ::Bool, iμ::Bool, fx::Bool, lλ::Float64, lμ::Float64)
"""
mutable struct cTfbd <: cT
  d1 ::cTfbd
  d2 ::cTfbd
  e  ::Float64
  iμ ::Bool
  iψ ::Bool
  fx ::Bool
  lλ ::Float64
  lμ ::Float64

  cTfbd() = new()
  cTfbd(e::Float64, iμ::Bool, iψ::Bool, fx::Bool, lλ::Float64, lμ::Float64) =
      (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
  cTfbd(d1::cTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool, lλ::Float64, lμ::Float64) =
      (x = new(); x.d1 = d1; x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
  cTfbd(d1::cTfbd, d2::cTfbd, e::Float64, iψ::Bool, iμ::Bool, fx::Bool, lλ::Float64, lμ::Float64) =
      new(d1, d2, e, iμ, iψ, fx, lλ, lμ)
end


# pretty-printing
function Base.show(io::IO, t::cTfbd)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane clads-bd fossil tree with ",
    nt , " tip",  (isone(nt) ? "" : "s" ),
    ", (", ntipsextinct(t)," extinct) and ",
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    cTfbd(tree::cTfbd)

Produce a new copy of `cTfbd`.
"""
function cTfbd(tree::cTfbd)
  if def1(tree)
    if def2(tree)
      cTfbd(cTfbd(tree.d1), 
            cTfbd(tree.d2), e(tree), isextinct(tree), isfossil(tree), 
            isfix(tree), lλ(tree), lμ(tree))
    else
      cTfbd(cTfbd(tree.d1), e(tree), isextinct(tree), isfossil(tree),
            isfix(tree), lλ(tree), lμ(tree))
    end
  else
    cTfbd(e(tree), isextinct(tree), isfossil(tree), 
      isfix(tree), lλ(tree), lμ(tree))
  end
end




"""
    cTfbd_wofe(tree::cTfbd)

Creates a copy of a `cTfbd` tree without fossils extinct tips.
"""
function cTfbd_wofe(tree::cTfbd)
  if def1(tree) && def2(tree)
    cTfbd(cTfbd_wofe(tree.d1), cTfbd_wofe(tree.d2),
      e(tree), isextinct(tree), isfossil(tree),
      isfix(tree), lλ(tree), lμ(tree))
  else
    cTfbd(e(tree), isextinct(tree), isfossil(tree), 
      isfix(tree), lλ(tree), lμ(tree))
  end
end




