#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#






"""
    sTxs

A composite recursive type of supertype `Tx`
representing a binary phylogenetic tree with traits `xv` and their rates `lσ2`
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    edge
  dt:   choice of time lag
  fdt:  final `dt`
  xv:   array of a Brownian motion evolution of `X`.
  lσ2:  array of a Brownian motion evolution of `log(σ)`.

  sTxs(d1 ::sTxs,
       d2 ::sTxs,
       e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       xv ::Array{Float64,1},
       lσ2::Array{Float64,1})
"""
mutable struct sTxs <: sT
  d1 ::sTxs
  d2 ::sTxs
  e  ::Float64
  dt ::Float64
  fdt::Float64
  xv ::Array{Float64,1}
  lσ2::Array{Float64,1}

  sTxs() = new()
  sTxs(e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       xv ::Array{Float64,1},
       lσ2::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt; x.xv = xv; x.lσ2 = lσ2; x)
  sTxs(d1 ::sTxs,
       e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       xv ::Array{Float64,1},
       lσ2::Array{Float64,1}) =
    (x = new(); 
      x.d1 = d1; x.e = e; x.dt = dt; x.fdt = fdt; x.xv = xv; x.lσ2 = lσ2; x)
  sTxs(d1 ::sTxs,
       d2 ::sTxs,
       e  ::Float64,
       dt ::Float64,
       fdt::Float64,
       xv ::Array{Float64,1},
       lσ2::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, xv, lσ2)
end

# pretty-printing
Base.show(io::IO, t::sTxs) =
  print(io, "insane simple `xs` tree with ", ntips(t), " tips")




"""
    sTxs(tree::sTxs)

Produce a new copy of `sTxs`.
"""
function sTxs(tree::sTxs)
  if def1(tree)
    if def2(tree)
      sTxs(sTxs(tree.d1), sTxs(tree.d2),
        e(tree), dt(tree), fdt(tree), copy(xv(tree)), copy(lσ2(tree)))
    else
      sTxs(sTxs(tree.d1),
        e(tree), dt(tree), fdt(tree), copy(xv(tree)), copy(lσ2(tree)))
    end
  else
    sTxs(e(tree), dt(tree), fdt(tree), copy(xv(tree)), copy(lσ2(tree)))
  end
end





"""
    iTxb

A composite recursive type of supertype `iT`
representing a binary phylogenetic tree with no extinction
and `λ` evolving as a Geometric Brownian motion  for `insane` use,
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   edge
  fx:  if fix tree
  dt:  choice of time lag
  fdt: final `δt`
  lλ:  array of a Brownian motion evolution of `log(λ)`
  xv:   array of a Brownian motion evolution of `X`.
  lσ2:  array of a Brownian motion evolution of `log(σ)`.

    iTxb()

Constructs an empty `iTxb` object.

    iTxb(e::Float64)

Constructs an empty `iTxb` object with pendant edge `pe`.

    iTxb(d1::iTxb, d2::iTxb, e::Float64)

Constructs an `iTxb` object with two `iTxb` daughters and pendant edge `pe`.
"""
mutable struct iTxb <: iT
  d1 ::iTxb
  d2 ::iTxb
  e  ::Float64
  dt ::Float64
  fdt::Float64
  fx ::Bool
  lλ ::Array{Float64,1}
  xv ::Array{Float64,1}
  lσ2::Array{Float64,1}

  iTxb() = new()
  iTxb(e::Float64, dt::Float64, fdt::Float64, fx::Bool, lλ::Array{Float64,1},
    xv::Array{Float64,1}, lσ2::Array{Float64,1}) =
      (x = new(); x.e = e; x.dt = dt; x.fdt = fdt; x.fx = fx; x.lλ = lλ;
        x.xv = xv; x.lσ2 = lσ2; x)
  iTxb(d1::iTxb, d2::iTxb, e::Float64, dt::Float64, fdt::Float64, fx::Bool, 
    lλ::Array{Float64,1}, xv::Array{Float64,1}, lσ2::Array{Float64,1}) =
      new(d1, d2, e, dt, fdt, fx, lλ, xv, lσ2)
end


# pretty-printing
Base.show(io::IO, t::iTxb) =
  print(io, "insane trait pure-birth tree with ", ntips(t), " tips")




# """
#     iTxb(e0::Array{Int64,1},
#          e1::Array{Int64,1},
#          el::Array{Float64,1},
#          λs::Array{Array{Float64,1},1},
#          ea::Array{Int64,1},
#          ni::Int64,
#          ei::Int64,
#          δt::Float64)

# Transform edge structure to `iTxb`.
# """
# function iTxb(e0::Array{Int64,1},
#               e1::Array{Int64,1},
#               el::Array{Float64,1},
#               λs::Array{Array{Float64,1},1},
#               ea::Array{Int64,1},
#               ni::Int64,
#               ei::Int64,
#               δt::Float64)

#   # if tip
#   if in(ei, ea)
#     return iTxb(el[ei], δt, δt, true, λs[ei])
#   else
#     ei1, ei2 = findall(isequal(ni), e0)
#     n1, n2   = e1[ei1:ei2]
#     return iTxb(iTxb(e0, e1, el, λs, ea, n1, ei1, δt),
#                 iTxb(e0, e1, el, λs, ea, n2, ei2, δt),
#                 el[ei], δt, (el[ei] == 0.0 ? 0.0 : δt), true, λs[ei])
#   end
# end




"""
    iTxb(tree::iTxb)

Produce a new copy of `iTxb`.
"""
function iTxb(tree::iTxb)
  if def1(tree)
    iTxb(iTxb(tree.d1), iTxb(tree.d2),
      e(tree), dt(tree), fdt(tree), isfix(tree), 
      copy(lλ(tree)), copy(xv(tree)), copy(lσ2(tree)))
  else
    iTxb(e(tree), dt(tree), fdt(tree), isfix(tree),
         copy(lλ(tree)), copy(xv(tree)), copy(lσ2(tree)))
  end
end


