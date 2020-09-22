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
    iT

An abstract type for all composite recursive types 
representing a simple binary phylogenetic tree for `insane` use
"""
abstract type iT <: iTree end




"""
    iTpb

The simplest composite recursive type of supertype `iTpbree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node

    iTpb()

Constructs an empty `iTpb` object.

    iTpb(pe::Float64)

Constructs an empty `iTpb` object with pendant edge `pe`.

    iTpb(d1::iTpb, d2::iTpb, pe::Float64)

Constructs an `iTpb` object with two `iTpb` daughters and pendant edge `pe`.
"""
mutable struct iTpb <: iT
  d1::Union{iTpb, Nothing}
  d2::Union{iTpb, Nothing}
  pe::Float64

  # inner constructor
  iTpb(d1::Union{iTpb, Nothing}, d2::Union{iTpb, Nothing}, pe::Float64) = 
    new(d1, d2, pe)
end

# outer constructors
iTpb() = iTpb(nothing, nothing, 0.0)

iTpb(pe::Float64) = iTpb(nothing, nothing, pe)


# pretty-printing
Base.show(io::IO, t::iTpb) = 
  print(io, "insane simple pure-birth tree with ", sntn(t), " tips")




"""
    iTbd

The simplest composite recursive type of supertype `iTbdree` 
representing a binary phylogenetic tree for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node

    iTbd()

Constructs an empty `iTbd` object.

    iTbd(pe::Float64)

Constructs an empty `iTbd` object with pendant edge `pe`.

    iTbd(d1::iTbd, d2::iTbd, pe::Float64)

Constructs an `iTbd` object with two `iTbd` daughters and pendant edge `pe`.
"""
mutable struct iTbd <: iT
  d1::Union{iTbd, Nothing}
  d2::Union{iTbd, Nothing}
  pe::Float64
  iμ::Bool
  fx::Bool

  # inner constructor
  iTbd(d1::Union{iTbd, Nothing}, d2::Union{iTbd, Nothing}, 
    pe::Float64, iμ::Bool, fx::Bool) = new(d1, d2, pe, iμ, fx)
end

# outer constructors
iTbd() = iTbd(nothing, nothing, 0.0, false, false)

iTbd(pe::Float64) = iTbd(nothing, nothing, pe, false, false)

iTbd(pe::Float64, iμ::Bool) = iTbd(nothing, nothing, pe, iμ, false)

iTbd(d1::iTbd, d2::iTbd, pe::Float64) = iTbd(d1, d2, pe, false, false)

iTbd(d1::iTbd, d2::iTbd, pe::Float64, iμ::Bool) = 
  iTbd(d1, d2, pe, iμ, false)

# pretty-printing
Base.show(io::IO, t::iTbd) = 
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

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  ts: array of time vectors concordant with `lλ`
  lλ: array of a Brownian motion evolution of `log(λ)`

    iTgbmpb()

Constructs an empty `iTgbmpb` object.

    iTgbmpb(pe::Float64)

Constructs an empty `iTgbmpb` object with pendant edge `pe`.

    iTgbmpb(d1::iTgbmpb, d2::iTgbmpb, pe::Float64)

Constructs an `iTgbmpb` object with two `iTgbmpb` daughters and pendant edge `pe`.
"""
mutable struct iTgbmpb <: iTgbm
  d1::Union{iTgbmpb, Nothing}
  d2::Union{iTgbmpb, Nothing}
  pe::Float64
  ts::Array{Float64,1}
  lλ::Array{Float64,1}

  # inner constructor
  iTgbmpb(d1::Union{iTgbmpb, Nothing}, d2::Union{iTgbmpb, Nothing}, 
    pe::Float64, ts::Array{Float64,1}, lλ::Array{Float64,1}) = 
      new(d1, d2, pe, ts, lλ)
end

# outer constructors
iTgbmpb() = 
  iTgbmpb(nothing, nothing, 0.0, Float64[], Float64[])

iTgbmpb(pe::Float64) = 
  iTgbmpb(nothing, nothing, pe, Float64[], Float64[])

iTgbmpb(d1::iTgbmpb, d2::iTgbmpb, pe::Float64) = 
  iTgbmpb(d1, d2, pe, Float64[], Float64[])

# pretty-printing
Base.show(io::IO, t::iTgbmpb) = 
  print(io, "insane pb-gbm tree with ", sntn(t), " tips")



"""
    iTgbmpb(tree::iTpb, δt::Float64, lλa::Float64, σλ::Float64)

Promotes an `iTpb` to `iTgbmpb` according to some values for `λ` diffusion.
"""
function iTgbmpb(tree::iTpb, 
                 δt::Float64, 
                 srδt::Float64, 
                 lλa::Float64, 
                 σλ::Float64)

  # make ts vector
  pet = pe(tree)
  tsv = [0.0:δt:pet...]
  if tsv[end] != pet
    push!(tsv, pet)
  end

  lλv = sim_bm(lλa, σλ, srδt, tsv)

  iTgbmpb(iTgbmpb(tree.d1, δt, srδt, lλv[end], σλ), 
          iTgbmpb(tree.d2, δt, srδt, lλv[end], σλ),
          pet, tsv, lλv)
end

"""
    iTgbmpb(::Nothing, δt::Float64, lλa::Float64, σλ::Float64)

Promotes an `iTpb` to `iTgbmpb` according to some values for `λ` diffusion.
"""
iTgbmpb(::Nothing, δt::Float64, srδt::Float64, lλa::Float64, σλ::Float64) = 
  nothing




"""
    iTgbmbd

A composite recursive type of supertype `iTgbm` 
representing a binary phylogenetic tree with  `λ` and `μ` 
evolving as a Geometric Brownian motion  for `insane` use, 
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  pe: pendant edge
  iμ: is an extinction node
  ts: array of time vectors concordant with `lλ` and `lμ`
  lλ: array of a Brownian motion evolution of `log(λ)`
  lμ: array of a Brownian motion evolution of `log(μ)`

    iTgbmbd()

Constructs an empty `iTgbmbd` object.

    iTgbmbd(pe::Float64)

Constructs an empty `iTgbmbd` object with pendant edge `pe`.

    iTgbmbd(d1::iTgbmbd, d2::iTgbmbd, pe::Float64)

Constructs an `iTgbmbd` object with two `iTgbmbd` daughters and pendant edge `pe`.
"""
mutable struct iTgbmbd <: iTgbm
  d1::Union{iTgbmbd, Nothing}
  d2::Union{iTgbmbd, Nothing}
  pe::Float64
  iμ::Bool
  fx::Bool
  ts::Array{Float64,1}
  lλ::Array{Float64,1}
  lμ::Array{Float64,1}

  # inner constructor
  iTgbmbd(d1::Union{iTgbmbd, Nothing}, d2::Union{iTgbmbd, Nothing}, 
    pe::Float64, iμ::Bool, fx::Bool, 
    ts::Array{Float64,1}, lλ::Array{Float64,1},lμ::Array{Float64,1}) = 
      new(d1, d2, pe, iμ, fx, ts, lλ, lμ)
end

# outer constructors
iTgbmbd() = iTgbmbd(nothing, nothing, 0.0, false, false, 
            Float64[], Float64[], Float64[])

iTgbmbd(pe::Float64) = iTgbmbd(nothing, nothing, pe, false, false,
                       Float64[], Float64[], Float64[])

iTgbmbd(pe::Float64, iμ::Bool) = iTgbmbd(nothing, nothing, pe, iμ, false,
                                 Float64[], Float64[], Float64[])

iTgbmbd(d1::iTgbmbd, d2::iTgbmbd, pe::Float64) = 
  iTgbmbd(d1, d2, pe, false, false,Float64[], Float64[], Float64[])

iTgbmbd(d1::iTgbmbd, d2::iTgbmbd, pe::Float64, iμ::Bool) = 
  iTgbmbd(d1, d2, pe, iμ, false, Float64[], Float64[], Float64[])

# pretty-printing
Base.show(io::IO, t::iTgbmbd) = 
  print(io, "insane bd-gbm tree with ", sntn(t), " tips (", snen(t)," extinct)")


