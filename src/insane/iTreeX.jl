#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sTpbx

The simplest composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  fx: if fix
  xi: initial trait value
  xf: final trait value

    sTpbx()

Constructs an empty `sTpbx` object.

    sTpbx(e::Float64, fx::Bool, xi::Float64, xf::Float64)

Constructs an `sTpbx` object with two `sTpbx` daughters and edge `e`,
fix information `fx`, initial node trait `xi` and final `xf`.
"""
mutable struct sTpbx <: sT
  d1::sTpbx
  d2::sTpbx
  e ::Float64
  fx::Bool
  xi::Float64
  xf::Float64

  sTpbx() = new()
  sTpbx(e::Float64, fx::Bool, xi::Float64, xf::Float64) =
    (x = new(); x.e = e; x.fx = fx; x.xi = xi; x.xf = xf; x)
  sTpbx(d1::sTpbx, d2::sTpbx, e::Float64, fx::Bool, xi::Float64, xf::Float64) =
    new(d1, d2, e, fx, xi, xf)
end

# pretty-printing
Base.show(io::IO, t::sTpbx) =
  print(io, "insane trait pure-birth tree with ", ntips(t), " tips")




"""
    sTpbx(tree::sTpbx)

Creates a copy of `sTpbx`.
"""
function sTpbx(tree::sTpbx)
  if def1(tree)
    sTpbx(sTpbx(tree.d1), sTpbx(tree.d2),
      e(tree), isfix(tree), xi(tree), xf(tree))
  else
    sTpbx(e(tree), isfix(tree), xi(tree), xf(tree))
  end
end




"""
    sTbdx

The simplest composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: if extinct
  fx: if fix
  xi: initial trait value
  xf: final trait value


    sTbdx(e::Float64, iμ::Bool, fx::Bool, xi::Float64, xf::Float64)

Constructs an `sTbdx` object with two `sTbdx` daughters and edge `e`,
fix information `fx`, initial node trait `xi` and final `xf`.
"""
mutable struct sTbdx <: sT
  d1::sTbdx
  d2::sTbdx
  e ::Float64
  iμ::Bool
  fx::Bool
  xi::Float64
  xf::Float64

  sTbdx() = new()
  sTbdx(e::Float64, iμ::Bool, fx::Bool,  xi::Float64, xf::Float64) =
    (t = new(); t.e = e; t.iμ = iμ;t.fx = fx; t.xi = xi; t.xf = xf; t)
  sTbdx(d1::sTbdx, d2::sTbdx, e::Float64, iμ::Bool, fx::Bool,
    xi::Float64, xf::Float64) =
      new(d1, d2, e, iμ, fx, xi, xf)
end

# pretty-printing
Base.show(io::IO, t::sTbdx) =
  print(io, "insane trait birth-death tree with ", ntips(t), " tips (",
    ntipsextinct(t)," extinct)")




"""
    sTbdx(tree::sTbdx)

Creates a copy of `sTbdx`.
"""
function sTbdx(tree::sTbdx)
  if def1(tree)
    sTbdx(sTbdx(tree.d1), sTbdx(tree.d2),
      e(tree), isextinct(tree), isfix(tree), xi(tree), xf(tree))
  else
    sTbdx(e(tree), isextinct(tree), isfix(tree), xi(tree), xf(tree))
  end
end




"""
    sTfbdx

The simplest composite recursive type of supertype `sTf`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: is an extinction node
  iψ: is a fossil node
  fx: if it is fix
  xi: initial trait value
  xf: final trait value

    sTfbdx(e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

    sTfbdx(d1::sTfbdx, e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

    sTfbdx(d1::sTfbdx, d2::sTfbdx, e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

Constructs an `sTfbdx` object with one sampled ancestor, one `sTfbdx` daughter and
edge `e`.
"""
mutable struct sTfbdx <: sT
  d1::sTfbdx
  d2::sTfbdx
  e ::Float64
  iμ::Bool
  iψ::Bool
  fx::Bool
  xi::Float64
  xf::Float64

  sTfbdx() = new()
  sTfbdx(e::Float64, iμ::Bool, iψ::Bool, fx::Bool, xi::Float64, xf::Float64) =
    (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.xi = xi; x.xf = xf;
      x)
  sTfbdx(d1::sTfbdx, e::Float64, iμ::Bool, iψ::Bool, fx::Bool,
    xi::Float64, xf::Float64) =
    (x = new(); x.d1 = d1; x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx;
      x.xi = xi; x.xf = xf; x)
  sTfbdx(d1::sTfbdx, d2::sTfbdx, e::Float64, iμ::Bool, iψ::Bool, fx::Bool,
    xi::Float64, xf::Float64) = new(d1, d2, e, iμ, iψ, fx, xi, xf)
end

# pretty-printing
function Base.show(io::IO, t::sTfbdx)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane trait simple fossil tree with ",
    nt , " tip",  (isone(nt) ? "" : "s" ),
    ", (", ntipsextinct(t)," extinct) and ",
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    sTfbdx(tree::sTfbdx)

Produces a copy of `sTfbdx`.
"""
function sTfbdx(tree::sTfbdx)
  if def1(tree)
    d1 = sTfbdx(tree.d1)
    if def2(tree)
      sTfbdx(d1, sTfbdx(tree.d2), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree), xi(tree), xf(tree))
    else
      sTfbdx(d1, e(tree),
        isextinct(tree), isfossil(tree), isfix(tree), xi(tree), xf(tree))
    end
  else
    sTfbdx(e(tree), isextinct(tree), isfossil(tree), isfix(tree),
      xi(tree), xf(tree))
  end
end






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
       lσ2 ::Array{Float64,1})
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
    iTpbx

A composite recursive type of supertype `iT`
representing a binary phylogenetic tree with  `λ` evolving as a
Geometric Brownian motion and no extinction, with traits and their rates
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    edge
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`
  xv:   array of a Brownian motion evolution of `X`.
  lσ2:   array of a Brownian motion evolution of `log(σ)`.

  iTpbx(d1 ::iTpbx,
        d2 ::iTpbx,
        e  ::Float64,
        fx ::Bool,
        dt ::Float64,
        fdt::Float64,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1},
        lσ2 ::Array{Float64,1})
"""
mutable struct iTpbx <: iT
  d1 ::iTpbx
  d2 ::iTpbx
  e  ::Float64
  dt ::Float64
  fdt::Float64
  fx ::Bool
  lλ ::Array{Float64,1}
  xv ::Array{Float64,1}
  lσ2 ::Array{Float64,1}

  iTpbx() = new()
  iTpbx(e  ::Float64,
        fx ::Bool,
        dt ::Float64,
        fdt::Float64,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1},
        lσ2 ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.fx = fx; x.lλ = lλ; x.xv = xv; lσ2 = x.lσ2; x)
  iTpbx(d1 ::iTpbx,
        d2 ::iTpbx,
        e  ::Float64,
        fx ::Bool,
        dt ::Float64,
        fdt::Float64,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1},
        lσ2 ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, fx, lλ, xv, lσ2)
end

# pretty-printing
Base.show(io::IO, t::iTpbx) =
  print(io, "insane trait ipb tree with ", ntips(t), " tips")




"""
    iTpbx(tree::iTpbx)

Produce a new copy of `iTpbx`.
"""
function iTpbx(tree::iTpbx)
  if def1(tree)
    iTpbx(iTpbx(tree.d1), iTpbx(tree.d2),
      e(tree), isfix(tree), dt(tree), fdt(tree),
      copy(lλ(tree)), copy(xv(tree)), copy(lσ2(tree)))
  else
    iTpbx(e(tree), isfix(tree), dt(tree), fdt(tree),
      copy(lλ(tree)), copy(xv(tree)), copy(lσ2(tree)))
  end
end





"""
    iTcex

A composite recursive type of supertype `iT`
representing a binary phylogenetic tree with  `λ` evolving as a
Geometric Brownian motion and constant `μ` for `insane` use,
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    edge
  iμ:   if extinct node
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`
  xv:   array of a Brownian motion evolution of `X`.

  iTcex(d1 ::iTcex,
          d2 ::iTcex,
          e  ::Float64,
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1})
"""
mutable struct iTcex <: iT
  d1 ::iTcex
  d2 ::iTcex
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  xv ::Array{Float64,1}

  iTcex() = new()
  iTcex(e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.xv = xv; x)
  iTcex(d1 ::iTcex,
        d2 ::iTcex,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, iμ, fx, lλ, xv)
end

# pretty-printing
Base.show(io::IO, t::iTcex) =
  print(io, "insane trait gbm-ce tree with ", ntips(t), 
    " tips (", ntipsextinct(t)," extinct)")




"""
    iTcex(tree::iTcex)

Produce a new copy of `iTcex`.
"""
function iTcex(tree::iTcex)
  if def1(tree)
    iTcex(iTcex(tree.d1), iTcex(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(xv(tree)))
  else
    iTcex(e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(xv(tree)))
  end
end




"""
    iTbdx

A composite recursive type of supertype `iT`
representing a binary phylogenetic tree with  `λ` and `μ`
evolving as a Geometric Brownian motion with fossils for `insane` use,
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   pendant edge
  iμ:  if extinct node
  fx:  if fix (observed) node
  dt:  choice of time lag
  fdt: final `dt`
  lλ:  array of a Brownian motion evolution of `log(λ)`
  lμ:  array of a Brownian motion evolution of `log(μ)`
  xv:  array of a Brownian motion evolution of `X`.

  iTbdx(d1 ::iTbdx,
        d2 ::iTbdx,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1},
        xv ::Array{Float64,1})
"""
mutable struct iTbdx <: iT
  d1 ::iTbdx
  d2 ::iTbdx
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  lμ ::Array{Float64,1}
  xv ::Array{Float64,1}

  iTbdx() = new()
  iTbdx(e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x.xv = xv; x)
  iTbdx(d1 ::iTbdx,
        d2 ::iTbdx,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, iμ, fx, lλ, lμ, xv)
end


# pretty-printing
function Base.show(io::IO, t::iTbdx)
  nt = ntips(t)

  print(io, "insane trait gbm-bd tree with ", 
    nt , " tip",  (isone(nt) ? "" : "s" ), 
    " (", ntipsextinct(t)," extinct)")
end




"""
    iTbdx(tree::iTbdx)

Produce a new copy of `iTbdx`.
"""
function iTbdx(tree::iTbdx)
  if def1(tree)
    iTbdx(iTbdx(tree.d1), iTbdx(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
  else
    iTbdx(e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
  end
end



"""
    iTfbdx

A composite recursive type of supertype `iT`
representing a binary phylogenetic tree with  `λ` and `μ`
evolving as a Geometric Brownian motion with fossils for `insane` use,
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  e:   pendant edge
  iμ:  if extinct node
  iψ:  if is fossil
  fx:  if fix (observed) node
  dt:  choice of time lag
  fdt: final `dt`
  lλ:  array of a Brownian motion evolution of `log(λ)`
  lμ:  array of a Brownian motion evolution of `log(μ)`
  xv:  array of a Brownian motion evolution of `X`.

  iTfbdx(d1 ::iTfbdx,
         d2 ::iTfbdx,
         e  ::Float64,
         dt ::Float64,
         fdt::Float64,
         iμ ::Bool,
         iψ ::Bool,
         fx ::Bool,
         lλ ::Array{Float64,1},
         lμ ::Array{Float64,1},
         xv ::Array{Float64,1})
"""
mutable struct iTfbdx <: iT
  d1 ::iTfbdx
  d2 ::iTfbdx
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  iψ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  lμ ::Array{Float64,1}
  xv ::Array{Float64,1}


  iTfbdx() = new()
  iTfbdx(e  ::Float64,
         dt ::Float64,
         fdt::Float64,
         iμ ::Bool,
         iψ ::Bool,
         fx ::Bool,
         lλ ::Array{Float64,1},
         lμ ::Array{Float64,1},
         xv ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x.xv = xv; x)
  iTfbdx(d1 ::iTfbdx,
         e  ::Float64,
         dt ::Float64,
         fdt::Float64,
         iμ ::Bool,
         iψ ::Bool,
         fx ::Bool,
         lλ ::Array{Float64,1},
         lμ ::Array{Float64,1},
         xv ::Array{Float64,1}) =
    (x = new(); x.d1 = d1; x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x.xv = xv; x)
  iTfbdx(d1 ::iTfbdx,
         d2 ::iTfbdx,
         e  ::Float64,
         dt ::Float64,
         fdt::Float64,
         iμ ::Bool,
         iψ ::Bool,
         fx ::Bool,
         lλ ::Array{Float64,1},
         lμ ::Array{Float64,1},
         xv ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, iμ, iψ, fx, lλ, lμ, xv)
end


# pretty-printing
function Base.show(io::IO, t::iTfbdx)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane trait ibd fossil tree with ", 
    nt , " tip",  (isone(nt) ? "" : "s" ), 
    ", (", ntipsextinct(t)," extinct) and ", 
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    iTfbdx(tree::iTfbdx)

Produce a new copy of `iTfbdx`.
"""
function iTfbdx(tree::iTfbdx)
  if def1(tree)
    if def2(tree)
      iTfbdx(iTfbdx(tree.d1), iTfbdx(tree.d2),
        e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
        isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
    else
      iTfbdx(iTfbdx(tree.d1),
        e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
        isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
    end
  else
    iTfbdx(e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
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

iTf = Union{sTfbd, sTfbdX, iTfbd}
"""
iTf = Union{sTf_label, sTfbd, cTfbd, iTfbd, sTfpe, sTxs}




"""
    Union type for unlabelled fossil data

uTf = Union{sTfbd, sTfbdX, iTfbd}
"""
uTf = Union{sTfbd, cTfbd, iTfbd, sTxs}




"""
    Union type for gbm-bd data

iTbdU = Union{iTbd, iTfbd}
"""
cTbdU = Union{cTbd, cTfbd}




"""
    Union type for gbm-bd data

iTbdU = Union{iTbd, iTfbd}
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


