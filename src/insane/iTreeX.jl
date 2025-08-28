#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sTbX

The simplest composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  fx: if fix
  xi: initial trait value
  xf: final trait value

    sTbX()

Constructs an empty `sTbX` object.

    sTbX(e::Float64, fx::Bool, xi::Float64, xf::Float64)

Constructs an `sTbX` object with two `sTbX` daughters and edge `e`,
fix information `fx`, initial node trait `xi` and final `xf`.
"""
mutable struct sTbX <: sT
  d1::sTbX
  d2::sTbX
  e ::Float64
  fx::Bool
  xi::Float64
  xf::Float64

  sTbX() = new()
  sTbX(e::Float64, fx::Bool, xi::Float64, xf::Float64) =
    (x = new(); x.e = e; x.fx = fx; x.xi = xi; x.xf = xf; x)
  sTbX(d1::sTbX, d2::sTbX, e::Float64, fx::Bool, xi::Float64, xf::Float64) =
    new(d1, d2, e, fx, xi, xf)
end

# pretty-printing
Base.show(io::IO, t::sTbX) =
  print(io, "insane trait pure-birth tree with ", ntips(t), " tips")




"""
    sTbX(tree::sTbX)

Creates a copy of `sTbX`.
"""
function sTbX(tree::sTbX)
  if def1(tree)
    sTbX(sTbX(tree.d1), sTbX(tree.d2),
      e(tree), isfix(tree), xi(tree), xf(tree))
  else
    sTbX(e(tree), isfix(tree), xi(tree), xf(tree))
  end
end




"""
    sTbdX

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


    sTbdX(e::Float64, iμ::Bool, fx::Bool, xi::Float64, xf::Float64)

Constructs an `sTbdX` object with two `sTbdX` daughters and edge `e`,
fix information `fx`, initial node trait `xi` and final `xf`.
"""
mutable struct sTbdX <: sT
  d1::sTbdX
  d2::sTbdX
  e ::Float64
  iμ::Bool
  fx::Bool
  xi::Float64
  xf::Float64

  sTbdX() = new()
  sTbdX(e::Float64, iμ::Bool, fx::Bool,  xi::Float64, xf::Float64) =
    (t = new(); t.e = e; t.iμ = iμ;t.fx = fx; t.xi = xi; t.xf = xf; t)
  sTbdX(d1::sTbdX, d2::sTbdX, e::Float64, iμ::Bool, fx::Bool,
    xi::Float64, xf::Float64) =
      new(d1, d2, e, iμ, fx, xi, xf)
end

# pretty-printing
Base.show(io::IO, t::sTbdX) =
  print(io, "insane trait birth-death tree with ", ntips(t), " tips (",
    ntipsextinct(t)," extinct)")




"""
    sTbdX(tree::sTbdX)

Creates a copy of `sTbdX`.
"""
function sTbdX(tree::sTbdX)
  if def1(tree)
    sTbdX(sTbdX(tree.d1), sTbdX(tree.d2),
      e(tree), isextinct(tree), isfix(tree), xi(tree), xf(tree))
  else
    sTbdX(e(tree), isextinct(tree), isfix(tree), xi(tree), xf(tree))
  end
end




"""
    sTfbdX

The simplest composite recursive type of supertype `sT`
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

    sTfbdX(e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

    sTfbdX(d1::sTfbdX, e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

    sTfbdX(d1::sTfbdX, d2::sTfbdX, e::Float64, iμ::Bool, iψ::Bool, fx::Bool)

Constructs an `sTfbdX` object with one sampled ancestor, one `sTfbdX` daughter and
edge `e`.
"""
mutable struct sTfbdX <: sT
  d1::sTfbdX
  d2::sTfbdX
  e ::Float64
  iμ::Bool
  iψ::Bool
  fx::Bool
  xi::Float64
  xf::Float64

  sTfbdX() = new()
  sTfbdX(e::Float64, iμ::Bool, iψ::Bool, fx::Bool, xi::Float64, xf::Float64) =
    (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x.xi = xi; x.xf = xf;
      x)
  sTfbdX(d1::sTfbdX, e::Float64, iμ::Bool, iψ::Bool, fx::Bool,
    xi::Float64, xf::Float64) =
    (x = new(); x.d1 = d1; x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx;
      x.xi = xi; x.xf = xf; x)
  sTfbdX(d1::sTfbdX, d2::sTfbdX, e::Float64, iμ::Bool, iψ::Bool, fx::Bool,
    xi::Float64, xf::Float64) = new(d1, d2, e, iμ, iψ, fx, xi, xf)
end

# pretty-printing
function Base.show(io::IO, t::sTfbdX)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane trait simple fossil tree with ",
    nt , " tip",  (isone(nt) ? "" : "s" ),
    ", (", ntipsextinct(t)," extinct) and ",
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    sTfbdX(tree::sTfbdX)

Produces a copy of `sTfbdX`.
"""
function sTfbdX(tree::sTfbdX)
  if def1(tree)
    d1 = sTfbdX(tree.d1)
    if def2(tree)
      sTfbdX(d1, sTfbdX(tree.d2), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree), xi(tree), xf(tree))
    else
      sTfbdX(d1, e(tree),
        isextinct(tree), isfossil(tree), isfix(tree), xi(tree), xf(tree))
    end
  else
    sTfbdX(e(tree), isextinct(tree), isfossil(tree), isfix(tree),
      xi(tree), xf(tree))
  end
end





"""
    iTbX

A composite recursive type of supertype `iT`
representing a binary phylogenetic tree with  `λ` evolving as a
Geometric Brownian motion and no extinction,
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  e:    edge
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`
  xv:   array of a Brownian motion evolution of `X`.

  iTbX(d1 ::iTbX,
          d2 ::iTbX,
          e  ::Float64,
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1})
"""
mutable struct iTbX <: iT
  d1 ::iTbX
  d2 ::iTbX
  e  ::Float64
  dt ::Float64
  fdt::Float64
  fx ::Bool
  lλ ::Array{Float64,1}
  xv ::Array{Float64,1}

  iTbX() = new()
  iTbX(e  ::Float64,
        fx ::Bool,
        dt ::Float64,
        fdt::Float64,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.fx = fx; x.lλ = lλ; x.xv = xv; x)
  iTbX(d1 ::iTbX,
        d2 ::iTbX,
        e  ::Float64,
        fx ::Bool,
        dt ::Float64,
        fdt::Float64,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    new(d1, d2, e, dt, fdt, fx, lλ, xv)
end

# pretty-printing
Base.show(io::IO, t::iTbX) =
  print(io, "insane trait gbm-pb tree with ", ntips(t), " tips")




"""
    iTbX(tree::iTbX)

Produce a new copy of `iTbX`.
"""
function iTbX(tree::iTbX)
  if def1(tree)
    iTbX(iTbX(tree.d1), iTbX(tree.d2),
      e(tree), isfix(tree), dt(tree), fdt(tree),
      copy(lλ(tree)), copy(xv(tree)))
  else
    iTbX(e(tree), isfix(tree), dt(tree), fdt(tree),
      copy(lλ(tree)), copy(xv(tree)))
  end
end





"""
    iTceX

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

  iTceX(d1 ::iTceX,
          d2 ::iTceX,
          e  ::Float64,
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool,
          fx ::Bool,
          lλ ::Array{Float64,1})
"""
mutable struct iTceX <: iT
  d1 ::iTceX
  d2 ::iTceX
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  xv ::Array{Float64,1}

  iTceX() = new()
  iTceX(e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.xv = xv; x)
  iTceX(d1 ::iTceX,
        d2 ::iTceX,
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
Base.show(io::IO, t::iTceX) =
  print(io, "insane trait gbm-ce tree with ", ntips(t), 
    " tips (", ntipsextinct(t)," extinct)")




"""
    iTceX(tree::iTceX)

Produce a new copy of `iTceX`.
"""
function iTceX(tree::iTceX)
  if def1(tree)
    iTceX(iTceX(tree.d1), iTceX(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(xv(tree)))
  else
    iTceX(e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(xv(tree)))
  end
end




"""
    iTbdX

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

  iTbdX(d1 ::iTbdX,
        d2 ::iTbdX,
        e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1},
        xv ::Array{Float64,1})
"""
mutable struct iTbdX <: iT
  d1 ::iTbdX
  d2 ::iTbdX
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  lμ ::Array{Float64,1}
  xv ::Array{Float64,1}

  iTbdX() = new()
  iTbdX(e  ::Float64,
        dt ::Float64,
        fdt::Float64,
        iμ ::Bool,
        fx ::Bool,
        lλ ::Array{Float64,1},
        lμ ::Array{Float64,1},
        xv ::Array{Float64,1}) =
    (x = new(); x.e = e; x.dt = dt; x.fdt = fdt;
      x.iμ = iμ; x.fx = fx; x.lλ = lλ; x.lμ = lμ; x.xv = xv; x)
  iTbdX(d1 ::iTbdX,
        d2 ::iTbdX,
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
function Base.show(io::IO, t::iTbdX)
  nt = ntips(t)

  print(io, "insane trait gbm-bd tree with ", 
    nt , " tip",  (isone(nt) ? "" : "s" ), 
    " (", ntipsextinct(t)," extinct)")
end




"""
    iTbdX(tree::iTbdX)

Produce a new copy of `iTbdX`.
"""
function iTbdX(tree::iTbdX)
  if def1(tree)
    iTbdX(iTbdX(tree.d1), iTbdX(tree.d2),
      e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
  else
    iTbdX(e(tree), dt(tree), fdt(tree), isextinct(tree),
      isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
  end
end



"""
    iTfbdX

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

  iTfbdX(d1 ::iTfbdX,
         d2 ::iTfbdX,
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
mutable struct iTfbdX <: iT
  d1 ::iTfbdX
  d2 ::iTfbdX
  e  ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  iψ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  lμ ::Array{Float64,1}
  xv ::Array{Float64,1}


  iTfbdX() = new()
  iTfbdX(e  ::Float64,
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
  iTfbdX(d1 ::iTfbdX,
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
  iTfbdX(d1 ::iTfbdX,
         d2 ::iTfbdX,
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
function Base.show(io::IO, t::iTfbdX)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane trait gbm-bd fossil tree with ", 
    nt , " tip",  (isone(nt) ? "" : "s" ), 
    ", (", ntipsextinct(t)," extinct) and ", 
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    iTfbdX(tree::iTfbdX)

Produce a new copy of `iTfbdX`.
"""
function iTfbdX(tree::iTfbdX)
  if def1(tree)
    if def2(tree)
      iTfbdX(iTfbdX(tree.d1), iTfbdX(tree.d2),
        e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
        isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
    else
      iTfbdX(iTfbdX(tree.d1),
        e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
        isfix(tree), copy(lλ(tree)), copy(lμ(tree)), copy(xv(tree)))
    end
  else
    iTfbdX(e(tree), dt(tree), fdt(tree), isextinct(tree), isfossil(tree),
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
    Union type for simple trait data

sTX = Union{sTbX, sTbdX, sTfbdX}
"""
sTX = Union{sTbX, sTbdX, sTfbdX}



"""
    Union type for gbm trait data

iTX = Union{iTbX, iTceX, iTbdX, iTfbdX}
"""
iTX = Union{iTbX, iTceX, iTbdX, iTfbdX}




"""
    Union type for fossil data

iTf = Union{sTf_label, sTfbd, sTfbdX, iTfbd, iTfbdX}
"""
iTf = Union{sTf_label, sTfbd, sTfbdX, iTfbd, iTfbdX}



"""
    Union type for unlabelled fossil data

uTf = Union{sTfbd, sTfbdX, iTfbd, iTfbdX}
"""
uTf = Union{sTfbd, sTfbdX, iTfbd, iTfbdX}



"""
    Union type for gbm-bd data

iTbdU = Union{iTbd, iTfbd, iTpbd, iTbdX, iTfbdX}
"""
iTbdU = Union{iTbd, iTfbd, iTpbd, iTbdX, iTfbdX}




"""
    Union type for gbm-bd data

iTbdUX = Union{iTbdX, iTfbdX}
"""
iTbdUX = Union{iTbdX, iTfbdX}



#=
Type aliases (for compatibility with older versions)
=#


const iTgbmb = iTb

const iTgbmce = iTce

const iTgbmct = iTct

const iTgbmbd = iTbd
