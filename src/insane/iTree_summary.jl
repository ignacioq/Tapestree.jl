#=

Summary functions for vector of `iTree`

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    extractp!(tree::T, p::Array{Float64,1}, lv::Function) where {T <: iTgbm}

Extract all the data augmented parameters returned by function `lv`.
"""
function extractp!(treev::Array{T,1}, vs::Array{Array{Float64,1},1}, lv::Function) where {T <: iTgbm}

  nt = lastindex(treev)

  for j in Base.OneTo(nt)
    append!(vs[j], lv(treev[j]))
  end

  if isnothing(treev[1].d1)
    treev1 = nothing
  else
    treev1 = iTgbmpb[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end 
  end

  if isnothing(treev[1].d2)
    treev2 = nothing
  else
    treev2 = iTgbmpb[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end 
  end

  extractp!(treev1, vs, lv)
  extractp!(treev2, vs, lv)
end

"""
    extractp!(::Nothing, p::Array{Float64,1}, lv::Function)

Extract all the data augmented parameters returned by function `lv`.
"""
extractp!(::Nothing, vs::Array{Array{Float64,1},1}, lv::Function) = nothing






"""
    extractp!(tree::T, p::Array{Float64,1}, lv::Function) where {T <: iTgbm}

Extract all the data augmented parameters returned by function `lv`.
"""
function extractp!(tree::T, p::Array{Float64,1}, lv::Function) where {T <: iTgbm}
  append!(p, lv(tree))

  extractp!(tree.d1, p, lv)
  extractp!(tree.d2, p, lv)
end

"""
    extractp!(::Nothing, p::Array{Float64,1}, lv::Function)

Extract all the data augmented parameters returned by function `lv`.
"""
extractp!(::Nothing, p::Array{Float64,1}, lv::Function) = nothing




"""
    iquantile(treev::Array{iTgbmpb,1}, p::Float64, lv::Function)


Make an `iTgbmpb` with the quantile specified by `p` in data specified in 
function `lv`.
"""
function iquantile(treev::Array{iTgbmpb,1}, p::Float64, lv::Function)

  tsv = ts(treev[1])
  nt = lastindex(treev)

  # make vector of lambdas
  vs = Array{Float64,1}[]
  for t in treev
    push!(vs, lv(t)) 
  end

  sv = Float64[]
  # make fill vector to estimate statistics
  v = Array{Float64,1}(undef, nt)
  for i in Base.OneTo(lastindex(tsv))
    for t in Base.OneTo(nt)
      v[t] = vs[t][i]
    end
    push!(sv, quantile(v, p))
  end

  if isnothing(treev[1].d1)
    treev1 = nothing
  else
    treev1 = iTgbmpb[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end 
  end

  if isnothing(treev[1].d2)
    treev2 = nothing
  else
    treev2 = iTgbmpb[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end 
  end

  iTgbmpb(iquantile(treev1, p, lv), iquantile(treev2, p, lv),
    tsv[end], tsv, sv)
end


"""
    iquantile(::Nothing, p::Float64, lv::Function)

Make an `iTgbmpb` with the quantile specified by `p` in data specified in 
function `lv`.
"""
iquantile(::Nothing, p::Float64, lv::Function) = 
  nothing