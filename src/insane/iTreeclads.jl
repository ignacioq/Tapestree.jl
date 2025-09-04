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




# """
#     iTce

# A composite recursive type of supertype `iT`
# representing a binary phylogenetic tree with  `λ` evolving as a
# Geometric Brownian motion and constant `μ` for `insane` use,
# with the following fields:

#   d1:   daughter tree 1
#   d2:   daughter tree 2
#   e:    edge
#   iμ:   if extinct node
#   fx:   if fix (observed) node
#   dt:   choice of time lag
#   fdt:  final `dt`
#   lλ:   array of a Brownian motion evolution of `log(λ)`

#   iTce(d1 ::iTce,
#           d2 ::iTce,
#           e  ::Float64,
#           dt ::Float64,
#           fdt::Float64,
#           iμ ::Bool,
#           fx ::Bool,
#           lλ ::Array{Float64,1})
# """
# mutable struct iTce <: iT
#   d1 ::iTce
#   d2 ::iTce
#   e  ::Float64
#   dt ::Float64
#   fdt::Float64
#   iμ ::Bool
#   fx ::Bool
#   lλ ::Array{Float64,1}

#   iTce() = new()
#   iTce(e  ::Float64,
#        dt ::Float64,
#        fdt::Float64,
#        iμ ::Bool,
#        fx ::Bool,
#        lλ ::Array{Float64,1}) =
#     (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
#       x.iμ = iμ; x.fx = fx; x.lλ = lλ; x)
#   iTce(d1 ::iTce,
#        d2 ::iTce,
#        e  ::Float64,
#        dt ::Float64,
#        fdt::Float64,
#        iμ ::Bool,
#        fx ::Bool,
#        lλ ::Array{Float64,1}) =
#     new(d1, d2, e, dt, fdt, iμ, fx, lλ)
# end

# # pretty-printing
# Base.show(io::IO, t::iTce) =
#   print(io, "insane gbm-ce tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




# """
#     iTce(tree::iTce)

# Produce a new copy of `iTce`.
# """
# function iTce(tree::iTce)
#   if def1(tree)
#     iTce(iTce(tree.d1), iTce(tree.d2),
#       e(tree), dt(tree), fdt(tree), isextinct(tree),
#       isfix(tree), copy(lλ(tree)))
#   else
#     iTce(e(tree), dt(tree), fdt(tree), isextinct(tree),
#       isfix(tree), copy(lλ(tree)))
#   end
# end




# """
#     iTce(e0::Array{Int64,1},
#             e1::Array{Int64,1},
#             el::Array{Float64,1},
#             λs::Array{Array{Float64,1},1},
#             ea::Array{Int64,1},
#             ee::Array{Int64,1},
#             ni::Int64,
#             ei::Int64,
#             δt::Float64)

# Transform edge structure to `iTce`.
# """
# function iTce(e0::Array{Int64,1},
#               e1::Array{Int64,1},
#               el::Array{Float64,1},
#               λs::Array{Array{Float64,1},1},
#               ea::Array{Int64,1},
#               ee::Array{Int64,1},
#               ni::Int64,
#               ei::Int64,
#               δt::Float64)

#   # if tip
#   if in(ei, ea)
#     return iTce(el[ei], δt, δt, false, false, λs[ei])
#   # if extinct
#   elseif in(ei, ee)
#     return iTce(el[ei], δt, δt, true, false, λs[ei])
#   else
#     ei1, ei2 = findall(isequal(ni), e0)
#     n1, n2   = e1[ei1:ei2]
#     return iTce(iTce(e0, e1, el, λs, ea, ee, n1, ei1, δt),
#                    iTce(e0, e1, el, λs, ea, ee, n2, ei2, δt),
#                    el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), false, false, λs[ei])
#   end
# end




# """
#     iTct

# A composite recursive type of supertype `iT`
# representing a binary phylogenetic tree with  `λ` evolving as a
# Geometric Brownian motion and constant `μ` for `insane` use,
# with the following fields:

#   d1:   daughter tree 1
#   d2:   daughter tree 2
#   e:    edge
#   iμ:   if extinct node
#   fx:   if fix (observed) node
#   dt:   choice of time lag
#   fdt:  final `dt`
#   lλ:   array of a Brownian motion evolution of `log(λ)`

#   iTct(d1 ::iTct,
#           d2 ::iTct,
#           e  ::Float64,
#           dt ::Float64,
#           fdt::Float64,
#           iμ ::Bool,
#           fx ::Bool,
#           lλ ::Array{Float64,1})
# """
# mutable struct iTct <: iT
#   d1 ::iTct
#   d2 ::iTct
#   e  ::Float64
#   dt ::Float64
#   fdt::Float64
#   iμ ::Bool
#   fx ::Bool
#   lλ ::Array{Float64,1}

#   iTct() = new()
#   iTct(e  ::Float64,
#        dt ::Float64,
#        fdt::Float64,
#        iμ ::Bool,
#        fx ::Bool,
#        lλ ::Array{Float64,1}) =
#     (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
#       x.iμ = iμ; x.fx = fx; x.lλ = lλ; x)
#   iTct(d1 ::iTct,
#        d2 ::iTct,
#        e  ::Float64,
#        dt ::Float64,
#        fdt::Float64,
#        iμ ::Bool,
#        fx ::Bool,
#        lλ ::Array{Float64,1}) =
#     new(d1, d2, e, dt, fdt, iμ, fx, lλ)
# end


# # pretty-printing
# Base.show(io::IO, t::iTct) =
#   print(io, "insane gbm-ct tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




# """
#     iTct(tree::iTct)

# Produce a new copy of `iTct`.
# """
# function iTct(tree::iTct)
#   if def1(tree)
#     iTct(iTct(tree.d1), iTct(tree.d2),
#       e(tree), dt(tree), fdt(tree), isextinct(tree),
#       isfix(tree), copy(lλ(tree)))
#   else
#     iTct(e(tree), dt(tree), fdt(tree), isextinct(tree),
#       isfix(tree), copy(lλ(tree)))
#   end
# end




# """
#     iTct(e0::Array{Int64,1},
#          e1::Array{Int64,1},
#          el::Array{Float64,1},
#          λs::Array{Array{Float64,1},1},
#          ea::Array{Int64,1},
#          ee::Array{Int64,1},
#          ni::Int64,
#          ei::Int64,
#          δt::Float64)

# Transform edge structure to `iTct`.
# """
# function iTct(e0::Array{Int64,1},
#               e1::Array{Int64,1},
#               el::Array{Float64,1},
#               λs::Array{Array{Float64,1},1},
#               ea::Array{Int64,1},
#               ee::Array{Int64,1},
#               ni::Int64,
#               ei::Int64,
#               δt::Float64)

#   # if tip
#   if in(ei, ea)
#     return iTct(el[ei], δt, δt, false, false, λs[ei])
#   # if extinct
#   elseif in(ei, ee)
#     return iTct(el[ei], δt, δt, true, false, λs[ei])
#   else
#     ei1, ei2 = findall(isequal(ni), e0)
#     n1, n2   = e1[ei1:ei2]
#     return iTct(iTct(e0, e1, el, λs, ea, ee, n1, ei1, δt),
#                    iTct(e0, e1, el, λs, ea, ee, n2, ei2, δt),
#                    el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), false, false, λs[ei])
#   end
# end




# """
#     iTbd

# A composite recursive type of supertype `iT`
# representing a binary phylogenetic tree with  `λ` and `μ`
# evolving as a Geometric Brownian motion  for `insane` use,
# with the following fields:

#   d1:   daughter tree 1
#   d2:   daughter tree 2
#   e:    pendant edge
#   iμ:   if extinct node
#   fx:   if fix (observed) node
#   dt:   choice of time lag
#   fdt:  final `dt`
#   lλ:   array of a Brownian motion evolution of `log(λ)`
#   lμ:   array of a Brownian motion evolution of `log(μ)`

#   iTbd(d1 ::iTbd,
#           d2 ::iTbd,
#           e  ::Float64,
#           dt ::Float64,
#           fdt::Float64,
#           iμ ::Bool,
#           fx ::Bool,
#           lλ ::Array{Float64,1},
#           lμ ::Array{Float64,1})
# """
# mutable struct iTbd <: iT
#   d1 ::iTbd
#   d2 ::iTbd
#   e  ::Float64
#   dt ::Float64
#   fdt::Float64
#   iμ ::Bool
#   fx ::Bool
#   lλ ::Array{Float64,1}
#   lμ ::Array{Float64,1}

#   iTbd() = new()
#   iTbd(e  ::Float64,
#        dt ::Float64,
#        fdt::Float64,
#        iμ ::Bool,
#        fx ::Bool,
#        lλ ::Array{Float64,1},
#        lμ ::Array{Float64,1}) =
#     (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
#       x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
#   iTbd(d1 ::iTbd,
#        d2 ::iTbd,
#        e  ::Float64,
#        dt ::Float64,
#        fdt::Float64,
#        iμ ::Bool,
#        fx ::Bool,
#        lλ ::Array{Float64,1},
#        lμ ::Array{Float64,1}) =
#     new(d1, d2, e, dt, fdt, iμ, fx, lλ, lμ)
# end


# # pretty-printing
# Base.show(io::IO, t::iTbd) =
#   print(io, "insane gbm-bd tree with ", ntips(t), " tips (", ntipsextinct(t)," extinct)")




# """
#     iTbd(tree::iTbd)

# Produce a new copy of `iTbd`.
# """
# function iTbd(tree::iTbd)
#   if def1(tree)
#     iTbd(iTbd(tree.d1), iTbd(tree.d2),
#       e(tree), dt(tree), fdt(tree), isextinct(tree),
#       isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
#   else
#     iTbd(e(tree), dt(tree), fdt(tree), isextinct(tree),
#       isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
#   end
# end




# """
#     iTbd(e0::Array{Int64,1},
#          e1::Array{Int64,1},
#          el::Array{Float64,1},
#          λs::Array{Array{Float64,1},1},
#          μs::Array{Array{Float64,1},1},
#          ea::Array{Int64,1},
#          ee::Array{Int64,1},
#          ni::Int64,
#          ei::Int64,
#          δt::Float64)

# Transform edge structure to `iTbd`.
# """
# function iTbd(e0::Array{Int64,1},
#               e1::Array{Int64,1},
#               el::Array{Float64,1},
#               λs::Array{Array{Float64,1},1},
#               μs::Array{Array{Float64,1},1},
#               ea::Array{Int64,1},
#               ee::Array{Int64,1},
#               ni::Int64,
#               ei::Int64,
#               δt::Float64)

#   # if tip
#   if in(ei, ea)
#     return iTbd(el[ei], δt, δt, false, false, λs[ei], μs[ei])
#   # if extinct
#   elseif in(ei, ee)
#     return iTbd(el[ei], δt, δt, true, false, λs[ei], μs[ei])
#   else
#     ei1, ei2 = findall(isequal(ni), e0)
#     n1, n2   = e1[ei1:ei2]
#     return iTbd(iTbd(e0, e1, el, λs, μs, ea, ee, n1, ei1, δt),
#                 iTbd(e0, e1, el, λs, μs, ea, ee, n2, ei2, δt),
#               el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt),
#               false, false, λs[ei], μs[ei])
#   end
# end




# """
#     iTfbd

# A composite recursive type of supertype `iT`
# representing a binary phylogenetic tree with  `λ` and `μ`
# evolving as a Geometric Brownian motion with fossils for `insane` use,
# with the following fields:

#   d1:   daughter tree 1
#   d2:   daughter tree 2
#   e:    pendant edge
#   iμ:   if extinct node
#   iψ:   if is fossil
#   fx:   if fix (observed) node
#   dt:   choice of time lag
#   fdt:  final `dt`
#   lλ:   array of a Brownian motion evolution of `log(λ)`
#   lμ:   array of a Brownian motion evolution of `log(μ)`

#   iTfbd(d1 ::iTfbd,
#         d2 ::iTfbd,
#         e  ::Float64,
#         dt ::Float64,
#         fdt::Float64,
#         iμ ::Bool,
#         iψ ::Bool,
#         fx ::Bool,
#         lλ ::Array{Float64,1},
#         lμ ::Array{Float64,1})
# """
# mutable struct iTfbd <: iT
#   d1 ::iTfbd
#   d2 ::iTfbd
#   e  ::Float64
#   dt ::Float64
#   fdt::Float64
#   iμ ::Bool
#   iψ ::Bool
#   fx ::Bool
#   lλ ::Array{Float64,1}
#   lμ ::Array{Float64,1}


#   iTfbd() = new()
#   iTfbd(e  ::Float64,
#         dt ::Float64,
#         fdt::Float64,
#         iμ ::Bool,
#         iψ ::Bool,
#         fx ::Bool,
#         lλ ::Array{Float64,1},
#         lμ ::Array{Float64,1}) =
#     (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
#       x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
#   iTfbd(d1 ::iTfbd,
#         e  ::Float64,
#         dt ::Float64,
#         fdt::Float64,
#         iμ ::Bool,
#         iψ ::Bool,
#         fx ::Bool,
#         lλ ::Array{Float64,1},
#         lμ ::Array{Float64,1}) =
#     (x = new(); x.d1 = d1; x.e = e; x.dt = dt; x.fdt = fdt;
#       x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x)
#   iTfbd(d1 ::iTfbd,
#         d2 ::iTfbd,
#         e  ::Float64,
#         dt ::Float64,
#         fdt::Float64,
#         iμ ::Bool,
#         iψ ::Bool,
#         fx ::Bool,
#         lλ ::Array{Float64,1},
#         lμ ::Array{Float64,1}) =
#     new(d1, d2, e, dt, fdt, iμ, iψ, fx, lλ, lμ)
# end


# # pretty-printing
# function Base.show(io::IO, t::iTfbd)
#   nt = ntips(t)
#   nf = nfossils(t)

#   print(io, "insane gbm-bd fossil tree with ",
#     nt , " tip",  (isone(nt) ? "" : "s" ),
#     ", (", ntipsextinct(t)," extinct) and ",
#     nf," fossil", (isone(nf) ? "" : "s" ))
# end




# """
#     iTfbd(tree::iTfbd)

# Produce a new copy of `iTfbd`.
# """
# function iTfbd(tree::iTfbd)
#   if def1(tree)
#     if def2(tree)
#       iTfbd(iTfbd(tree.d1), iTfbd(tree.d2),
#         e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
#         isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
#     else
#       iTfbd(iTfbd(tree.d1),
#         e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
#         isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
#     end
#   else
#     iTfbd(e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
#       isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
#   end
# end




# """
#     iTfbd_wofe(tree::iTfbd)

# Creates a copy of a `iTfbd` tree without fossils extinct tips.
# """
# function iTfbd_wofe(tree::iTfbd)
#   if def1(tree) && def2(tree)
#     iTfbd(iTfbd_wofe(tree.d1), iTfbd_wofe(tree.d2),
#       e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
#       isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
#   else
#     iTfbd(e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
#       isfix(tree), copy(lλ(tree)), copy(lμ(tree)))
#   end
# end




# """
#     iTfbd(e0::Array{Int64,1},
#           e1::Array{Int64,1},
#           el::Array{Float64,1},
#           λs::Array{Array{Float64,1},1},
#           μs::Array{Array{Float64,1},1},
#           ea::Array{Int64,1},
#           ee::Array{Int64,1},
#           ef::Array{Int64,1},
#           ni::Int64,
#           ei::Int64,
#           δt::Float64)

# Transform edge structure to `iTfbd`.
# """
# function iTfbd(e0::Array{Int64,1},
#                e1::Array{Int64,1},
#                el::Array{Float64,1},
#                λs::Array{Array{Float64,1},1},
#                μs::Array{Array{Float64,1},1},
#                ea::Array{Int64,1},
#                ee::Array{Int64,1},
#                ef::Array{Int64,1},
#                ni::Int64,
#                ei::Int64,
#                δt::Float64)

#   # if tip
#   if in(ei, ea)
#     return iTfbd(el[ei], δt, δt, false, false, false, λs[ei], μs[ei])
  
#   # if extinct
#   elseif in(ei, ee)
#     return iTfbd(el[ei], δt, δt, true, false, false, λs[ei], μs[ei])

#   # if fossil
#   elseif in(ei, ef)
#     ei1 = findfirst(isequal(ni), e0)
#     n1  = e1[ei1]
#     return iTfbd(iTfbd(e0, e1, el, λs, μs, ea, ee, ef, n1, ei1, δt),
#                  el[ei], δt, δt, false, true, false, λs[ei], μs[ei])

#   # if internal
#   else
#     ei1, ei2 = findall(isequal(ni), e0)
#     n1, n2   = e1[ei1:ei2]
#     return iTfbd(iTfbd(e0, e1, el, λs, μs, ea, ee, ef, n1, ei1, δt),
#                  iTfbd(e0, e1, el, λs, μs, ea, ee, ef, n2, ei2, δt),
#                  el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt),
#                  false, false, false, λs[ei], μs[ei])
#   end
# end



