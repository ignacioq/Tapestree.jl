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






#=
Union Types
=#




"""
    Union type for label trees

Tlabel = Union{sT_label, sTf_label}
"""
Tlabel = Union{sT_label, sTf_label}




"""
    Union type for fossil data

iTf = Union{sTf_label, sTfbd, sTfbdX, iTfbd, iTfbdX}
"""
iTf = Union{sTf_label, sTfbd, iTfbd, sTfpe, sTxs}




"""
    Union type for unlabelled fossil data

uTf = Union{sTfbd, sTfbdX, iTfbd, iTfbdX}
"""
uTf = Union{sTfbd, iTfbd, sTxs}




"""
    Union type for gbm-bd data

iTbdUX = Union{iTbdX, iTfbdX}
"""
iTbdU = Union{iTbd, iTfbd}




"""
    Union type for simple trait data

Tx = Union{sTpbx, sTbdx, sTfbdx}
"""
Tx = Union{sTxs, sTpe, sTfpe}




"""
    Union type for trait and rate data

Txs = Union{sTxs}
"""
Txs = Union{sTxs}



"""
    Union type for trait and rate data

Txs = Union{sTxs}
"""
Tpe = Union{sTpe, sTfpe}


