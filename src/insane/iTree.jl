#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    iTree

An abstract type for all composite recursive types 
representing a binary phylogenetic tree for `insane` use
"""
abstract type iTree end




"""
    sT

An abstract type for all composite recursive types 
representing a simple binary phylogenetic tree for `insane` use
"""
abstract type sT <: iTree end




"""
    sTpb

The simplest composite recursive type of supertype `iTree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node

    sTpb()

Constructs an empty `sTpb` object.

    sTpb(pe::Float64)

Constructs an empty `sTpb` object with pendant edge `pe`.

    sTpb(d1::sTpb, d2::sTpb, pe::Float64)

Constructs an `sTpb` object with two `sTpb` daughters and pendant edge `pe`.
"""
mutable struct sTpb <: sT
  d1::Union{sTpb, Nothing}
  d2::Union{sTpb, Nothing}
  pe::Float64

  # inner constructor
  sTpb(d1::Union{sTpb, Nothing}, d2::Union{sTpb, Nothing}, pe::Float64) = 
    new(d1, d2, pe)
end

# outer constructors
sTpb() = sTpb(nothing, nothing, 0.0)

sTpb(pe::Float64) = sTpb(nothing, nothing, pe)


# pretty-printing
Base.show(io::IO, t::sTpb) = 
  print(io, "insane simple pure-birth tree with ", sntn(t), " tips")




"""
    sTbd

The simplest composite recursive type of supertype `sTbdree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node
  fx: if it is fix

    sTbd()

Constructs an empty `sTbd` object.

    sTbd(pe::Float64)

Constructs an empty `sTbd` object with pendant edge `pe`.

    sTbd(d1::sTbd, d2::sTbd, pe::Float64)

Constructs an `sTbd` object with two `sTbd` daughters and pendant edge `pe`.
"""
mutable struct sTbd <: sT
  d1::Union{sTbd, Nothing}
  d2::Union{sTbd, Nothing}
  pe::Float64
  iμ::Bool
  fx::Bool

  # inner constructor
  sTbd(d1::Union{sTbd, Nothing}, d2::Union{sTbd, Nothing}, 
    pe::Float64, iμ::Bool, fx::Bool) = new(d1, d2, pe, iμ, fx)
end

# outer constructors
sTbd() = sTbd(nothing, nothing, 0.0, false, false)

sTbd(pe::Float64) = sTbd(nothing, nothing, pe, false, false)

sTbd(pe::Float64, iμ::Bool) = sTbd(nothing, nothing, pe, iμ, false)

sTbd(d1::sTbd, d2::sTbd, pe::Float64) = sTbd(d1, d2, pe, false, false)

sTbd(d1::sTbd, d2::sTbd, pe::Float64, iμ::Bool) = 
  sTbd(d1, d2, pe, iμ, false)

# pretty-printing
Base.show(io::IO, t::sTbd) = 
  print(io, "insane simple birth-death tree with ", sntn(t), " tips (", snen(t)," extinct)")




"""
    iTgbm

An abstract type for all composite recursive types 
representing a binary phylogenetic tree with Geometric
Brownian motion rates for `insane` use
"""
abstract type iTgbm <: iTree end




"""
    iTgbmpb

A composite recursive type of supertype `iTgbm` 
representing a binary phylogenetic tree with no extinction
and `λ` evolving as a Geometric Brownian motion  for `insane` use, 
with the following fields:

  d1:  daughter tree 1
  d2:  daughter tree 2
  pe:  pendant edge
  dt:  choice of time lag
  fdt: final `δt`
  lλ:  array of a Brownian motion evolution of `log(λ)`

    iTgbmpb()

Constructs an empty `iTgbmpb` object.

    iTgbmpb(pe::Float64)

Constructs an empty `iTgbmpb` object with pendant edge `pe`.

    iTgbmpb(d1::iTgbmpb, d2::iTgbmpb, pe::Float64)

Constructs an `iTgbmpb` object with two `iTgbmpb` daughters and pendant edge `pe`.
"""
mutable struct iTgbmpb <: iTgbm
  d1 ::Union{iTgbmpb, Nothing}
  d2 ::Union{iTgbmpb, Nothing}
  pe ::Float64
  dt ::Float64
  fdt::Float64
  lλ ::Array{Float64,1}

  # inner constructor
  iTgbmpb(d1::Union{iTgbmpb, Nothing}, d2::Union{iTgbmpb, Nothing}, pe::Float64, 
    dt::Float64, fdt::Float64, lλ::Array{Float64,1}) = 
      new(d1, d2, pe, dt, fdt, lλ)
end

# outer constructors
iTgbmpb() = 
  iTgbmpb(nothing, nothing, 0.0, 0.0, 0.0, Float64[])

# pretty-printing
Base.show(io::IO, t::iTgbmpb) = 
  print(io, "insane pb-gbm tree with ", sntn(t), " tips")



"""
    iTgbmpb(tree::sTpb, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            σλ  ::Float64)

Promotes an `sTpb` to `iTgbmpb` according to some values for `λ` diffusion.
"""
function iTgbmpb(tree::sTpb, 
                 δt  ::Float64, 
                 srδt::Float64, 
                 lλa ::Float64, 
                 σλ  ::Float64)

  pet  = pe(tree)

  if iszero(pet)
    iTgbmpb(iTgbmpb(tree.d1, δt, srδt, lλa, σλ), 
            iTgbmpb(tree.d2, δt, srδt, lλa, σλ),
            pe(tree), δt, 0.0, Float64[lλa])

  else
    nt   = Int64(fld(pet,δt))
    fdti = mod(pet, δt)

    if iszero(fdti)
      fdti = δt
    end

    lλv = sim_bm(lλa, σλ, srδt, nt, fdti)
    l   = lastindex(lλv)

    iTgbmpb(iTgbmpb(tree.d1, δt, srδt, lλv[l], σλ), 
            iTgbmpb(tree.d2, δt, srδt, lλv[l], σλ),
            pe(tree), δt, fdti, lλv)

  end
end

"""
    iTgbmpb(::Nothing, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            σλ  ::Float64)

Promotes an `sTpb` to `iTgbmpb` according to some values for `λ` diffusion.
"""
iTgbmpb(::Nothing, δt::Float64, srδt::Float64, lλa::Float64, σλ::Float64) = 
  nothing




"""
    iTgbmce

A composite recursive type of supertype `iTgbm` 
representing a binary phylogenetic tree with  `λ` evolving as a 
Geometric Brownian motion and constant `μ` for `insane` use, 
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  pe:   pendant edge
  iμ:   if extinct node
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`

  iTgbmce(d1  ::Union{iTgbmce, Nothing}, 
          d2  ::Union{iTgbmce, Nothing}, 
          pe  ::Float64, 
          dt  ::Float64,
          fdt::Float64,
          iμ   ::Bool, 
          fx   ::Bool, 
          lλ   ::Array{Float64,1})

Constructs an `iTgbmce` object with two `iTgbmce` daughters and pendant edge `pe`.
"""
mutable struct iTgbmce <: iTgbm
  d1 ::Union{iTgbmce, Nothing}
  d2 ::Union{iTgbmce, Nothing}
  pe ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}

  # inner constructor
  iTgbmce(d1 ::Union{iTgbmce, Nothing}, 
          d2 ::Union{iTgbmce, Nothing}, 
          pe ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1}) = 
    new(d1, d2, pe, dt, fdt, iμ, fx, lλ)
end

# outer constructors
iTgbmce() = iTgbmce(nothing, nothing, 0.0, 0.0, 0.0, false, false, Float64[])

# outer constructors
iTgbmce(d1 ::Union{iTgbmce, Nothing}, 
        d2 ::Union{iTgbmce, Nothing}, 
        pe ::Float64, 
        dt ::Float64,
        fdt::Float64,
        lλ ::Array{Float64,1}) = 
  iTgbmce(d1, d2, pe, dt, fdt, false, false, lλ)

# pretty-printing
Base.show(io::IO, t::iTgbmce) = 
  print(io, "insane gbm-ce tree with ", sntn(t), " tips (", snen(t)," extinct)")




"""
    sTbd(tree::iTgbmce)

Demotes a tree of type `iTgbmce` to `sTbd`.
"""
sTbd(tree::iTgbmce) =
  sTbd(sTbd(tree.d1), sTbd(tree.d2), pe(tree), isextinct(tree), false)

"""
    sTbd(tree::Nothing) 

Demotes a tree of type `iTgbmce` to `sTbd`.
"""
sTbd(tree::Nothing) = nothing




"""
    iTgbmce(tree::sTbd, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            σλ  ::Float64,
            σμ  ::Float64)

Promotes an `sTbd` to `iTgbmce` according to some values for `λ` and `μ` 
diffusion.
"""
function iTgbmce(tree::sTbd, 
                 δt  ::Float64, 
                 srδt::Float64, 
                 lλa ::Float64, 
                 σλ  ::Float64,
                 σμ  ::Float64)

  pet  = pe(tree)

  # if crown root
  if iszero(pet)
    iTgbmce(iTgbmce(tree.d1, δt, srδt, lλa, σλ, σμ), 
            iTgbmce(tree.d2, δt, srδt, lλa σλ, σμ),
            pe(tree), δt, 0.0, isextinct(tree), isfix(tree), Float64[lλa])
  else
    nt   = Int64(fld(pet,δt))
    fdti = mod(pet, δt)

    lλv = sim_bm(lλa, σλ, srδt, nt, fdti)

    if iszero(fdti)
      fdti = δt
    end

    l = lastindex(lλv)

    iTgbmce(iTgbmce(tree.d1, δt, srδt, lλv[l], σλ, σμ), 
            iTgbmce(tree.d2, δt, srδt, lλv[l], σλ, σμ),
            pe(tree), δt, fdti, isextinct(tree), isfix(tree), lλv)
  end
end


"""
    iTgbmce(::Nothing, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            lμa ::Float64, 
            σλ  ::Float64,
            σμ  ::Float64)

Promotes an `sTbd` to `iTgbmce` according to some values for `λ` and `μ` 
diffusion.
"""
iTgbmce(::Nothing, 
        δt  ::Float64, 
        srδt::Float64, 
        lλa ::Float64, 
        σλ  ::Float64,
        σμ  ::Float64) = nothing




"""
    iTgbmbd

A composite recursive type of supertype `iTgbm` 
representing a binary phylogenetic tree with  `λ` and `μ` 
evolving as a Geometric Brownian motion  for `insane` use, 
with the following fields:

  d1:   daughter tree 1
  d2:   daughter tree 2
  pe:   pendant edge
  iμ:   if extinct node
  fx:   if fix (observed) node
  dt:   choice of time lag
  fdt:  final `dt`
  lλ:   array of a Brownian motion evolution of `log(λ)`
  lμ:   array of a Brownian motion evolution of `log(μ)`

  iTgbmbd(d1  ::Union{iTgbmbd, Nothing}, 
          d2  ::Union{iTgbmbd, Nothing}, 
          pe  ::Float64, 
          dt  ::Float64,
          fdt::Float64,
          iμ   ::Bool, 
          fx   ::Bool, 
          lλ   ::Array{Float64,1},
          lμ   ::Array{Float64,1})

Constructs an `iTgbmbd` object with two `iTgbmbd` daughters and pendant edge `pe`.
"""
mutable struct iTgbmbd <: iTgbm
  d1 ::Union{iTgbmbd, Nothing}
  d2 ::Union{iTgbmbd, Nothing}
  pe ::Float64
  dt ::Float64
  fdt::Float64
  iμ ::Bool
  fx ::Bool
  lλ ::Array{Float64,1}
  lμ ::Array{Float64,1}

  # inner constructor
  iTgbmbd(d1 ::Union{iTgbmbd, Nothing}, 
          d2 ::Union{iTgbmbd, Nothing}, 
          pe ::Float64, 
          dt ::Float64,
          fdt::Float64,
          iμ ::Bool, 
          fx ::Bool, 
          lλ ::Array{Float64,1},
          lμ ::Array{Float64,1}) = 
    new(d1, d2, pe, dt, fdt, iμ, fx, lλ, lμ)
end

# outer constructors
iTgbmbd() = iTgbmbd(nothing, nothing, 0.0, 0.0, 0.0, false, false, 
            Float64[], Float64[])

# outer constructors
iTgbmbd(d1 ::Union{iTgbmbd, Nothing}, 
        d2 ::Union{iTgbmbd, Nothing}, 
        pe ::Float64, 
        dt ::Float64,
        fdt::Float64,
        lλ ::Array{Float64,1}) = 
  iTgbmbd(d1, d2, pe, dt, fdt, false, false, lλ, Float64[])

# pretty-printing
Base.show(io::IO, t::iTgbmbd) = 
  print(io, "insane bd-gbm tree with ", sntn(t), " tips (", snen(t)," extinct)")




"""
    sTbd(tree::iTgbmbd)

Demotes a tree of type `iTgbmbd` to `sTbd`.
"""
sTbd(tree::iTgbmbd) =
  sTbd(sTbd(tree.d1), sTbd(tree.d2), pe(tree), isextinct(tree), false)

"""
    sTbd(tree::Nothing) 

Demotes a tree of type `iTgbmbd` to `sTbd`.
"""
sTbd(tree::Nothing) = nothing




"""
    iTgbmbd(tree::sTbd, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            lμa ::Float64, 
            σλ  ::Float64,
            σμ  ::Float64)

Promotes an `sTbd` to `iTgbmbd` according to some values for `λ` and `μ` 
diffusion.
"""
function iTgbmbd(tree::sTbd, 
                 δt  ::Float64, 
                 srδt::Float64, 
                 lλa ::Float64, 
                 lμa ::Float64, 
                 σλ  ::Float64,
                 σμ  ::Float64)

  pet  = pe(tree)

  # if crown root
  if iszero(pet)
    iTgbmbd(iTgbmbd(tree.d1, δt, srδt, lλa, lμa, σλ, σμ), 
            iTgbmbd(tree.d2, δt, srδt, lλa, lμa, σλ, σμ),
            pe(tree), δt, 0.0, isextinct(tree), isfix(tree), 
            Float64[lλa], Float64[lμa])
  else
    nt   = Int64(fld(pet,δt))
    fdti = mod(pet, δt)

    lλv = sim_bm(lλa, σλ, srδt, nt, fdti)
    lμv = sim_bm(lμa, σμ, srδt, nt, fdti)

    if iszero(fdti)
      fdti = δt
    end

    l = lastindex(lμv)

    iTgbmbd(iTgbmbd(tree.d1, δt, srδt, lλv[l], lμv[l], σλ, σμ), 
            iTgbmbd(tree.d2, δt, srδt, lλv[l], lμv[l], σλ, σμ),
            pe(tree), δt, fdti, isextinct(tree), isfix(tree), lλv, lμv)
  end
end


"""
    iTgbmbd(::Nothing, 
            δt  ::Float64, 
            srδt::Float64, 
            lλa ::Float64, 
            lμa ::Float64, 
            σλ  ::Float64,
            σμ  ::Float64)

Promotes an `sTbd` to `iTgbmbd` according to some values for `λ` and `μ` 
diffusion.
"""
iTgbmbd(::Nothing, 
        δt  ::Float64, 
        srδt::Float64, 
        lλa ::Float64, 
        lμa ::Float64, 
        σλ  ::Float64,
        σμ  ::Float64) = nothing

