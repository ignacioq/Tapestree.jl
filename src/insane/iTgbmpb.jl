#=

Geometric Brownian motion insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




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



# conversion
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



