#=

Summary functions for vector of `iTree`

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    mcmc_array(treev::Array{T,1},
               δt   ::Float64,
               lv   ::Function) where {T <: iTgbm}

Return an Array with a row for each sampled tree for interpolated
parameters accessed by `lv` at times determined by `δt`.
"""
function mcmc_array(treev::Array{T,1},
                    δt   ::Float64,
                    lv   ::Function) where {T <: iTgbm}

  @inbounds begin

    nti = extractp(treev[1], δt, lv)
    vi = Float64[]
    linearize_gbm!(nti, lv, vi)

    ra = Array{Float64,2}(undef, lastindex(treev), lastindex(vi))
    ra[1,:] = vi

    for i in 2:lastindex(treev)
      nti = extractp(treev[i], δt, lv)
      vi  = Float64[]
      linearize_gbm!(nti, lv, vi)
      ra[i, :] = vi
    end
  end

  return ra
end




"""
    linearize_gbm!(tree::T, 
                   lv  ::Function,
                   v   ::Array{Float64,1}) where {T <: iTgbm}

Extract the parameters given by `lv` into a linear Array.
"""
function linearize_gbm!(tree::T, 
                        lv  ::Function,
                        v   ::Array{Float64,1}) where {T <: iTgbm}

  append!(v, lv(tree))
  if !istip(tree.d1)
    linearize_gbm!(tree.d1::iTgbmpb, lv, v)
  end
  if !istip(tree.d2)
    linearize_gbm!(tree.d2::iTgbmpb, lv, v)
  end
end




"""
    extractp(tree::T, δt::Float64, lv::Function) where {T <: iTgbm}

Log-linearly predict Geometric Brownian motion for `lv` at times given by `δt`.
"""
function extractp(tree::T, δt::Float64, lv::Function) where {T <: iTgbm}

  t = ts(tree)
  v = lv(tree)

  # make ts vector
  pet = pe(tree)
  tsv = [0.0:δt:pet...]
  if istip(tree) && tsv[end] != pet
    push!(tsv, pet)
  end

  pv = Float64[]

  # linearly predict λ
  for i in tsv
    ix, out = idxrange(t, i)
    if out
      push!(pv, linpred(i, t[ix], t[ix+1], v[ix], v[ix+1]))
    else
      push!(pv, v[ix])
    end
  end

  iTgbmpb(extractp(tree.d1, δt, lv), 
          extractp(tree.d2, δt, lv),
          pet, tsv, pv)
end

"""
    extractp(tree::Nothing, δt::Float64, lv::Function)

Log-linearly predict Geometric Brownian motion for `lv` at times given by `δt`.
"""
extractp(tree::Nothing, δt::Float64, lv::Function) = nothing




"""
    linpred(val::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64)

Estimate val according to linear interpolation for a range.
"""
linpred(val::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64) = 
  (y1 + (val - x1)*(y2 - y1)/(x2 - x1))::Float64




"""
    idxrange(x::Array{Float64,1}, val::Float64)

Get indexes in sorted vector `x` corresponding to the range in which 
`val` is in using a sort of uniroot algorithm.
"""
function idxrange(x::Array{Float64,1}, val::Float64)
  
  @inbounds begin

    a::Int64 = 1

    if x[a] > val
      return a, false
    end

    b::Int64 = lastindex(x)
  
    if x[b] < val
      return b, false
    end

    mid::Int64 = div(b,2)

    while b-a > 1
      val < x[mid] ? b = mid : a = mid
      mid = div(b + a, 2)
    end

    if x[a] == val 
      return a, false
    elseif x[b] == val
      return b, false
    else
      return a, true
    end

  end
end 




"""
    iquantile(treev::Array{T,1}, 
              p    ::Float64, 
              lv   ::Function) where {T <: iTgbm}

Make an `iTgbmpb` with the quantile specified by `p` in data specified in 
function `lv`.
"""
function iquantile(treev::Array{T,1}, 
                   p    ::Float64, 
                   lv   ::Function) where {T <: iTgbm}

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
iquantile(::Nothing, p::Float64, lv::Function) = nothing



