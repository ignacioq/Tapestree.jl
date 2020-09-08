#=

Geometric Brownian motion insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




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



# conversion
# iTgbmbd(tree::iT) = 
#   iTgbmbd(d1, d2, pe, iμ, false, Float64[], Float64[], Float64[])

# convert(::Type{iTgbmbd}, x::iT) = iTgbmbd(x)








