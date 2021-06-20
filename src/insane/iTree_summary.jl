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

    nti = extractp(treev[1], δt)
    vi = Float64[]
    linearize_gbm!(nti, lv, vi)

    ra = Array{Float64,2}(undef, lastindex(treev), lastindex(vi))
    ra[1,:] = vi

    for i in 2:lastindex(treev)
      nti = extractp(treev[i], δt)
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
    linearize_gbm!(tree.d1::T, lv, v)
  end
  if !istip(tree.d2)
    linearize_gbm!(tree.d2::T, lv, v)
  end
end




"""
    extractp(tree::T, nδt::Float64, lv::Function) where {T <: iTgbm}

Log-linearly predict Geometric Brownian motion for `lv` at times given by `nδt`.
"""
function extractp(tree::iTgbmpb, nδt::Float64)

  pet = pe(tree)
  δt  = dt(tree)
  v   = lλ(tree)
  n   = floor(Int64, pet/nδt)

  pv = Float64[]
  for i in Base.OneTo(n+1)
    iti = Float64(i-1)*nδt
    ix  = floor(Int64, iti/δt) + 1
    tts = δt*Float64(ix-1)
    ttf = δt*Float64(ix)
    push!(pv, linpred(iti, tts, ttf, v[ix], v[ix+1]))
  end

  if (Float64(n+1)*nδt) !== pet
    push!(pv, v[end])
    nfdt = mod(pet,nδt)
  else
    nfdt = nδt
  end

  iTgbmpb(extractp(tree.d1, nδt), 
          extractp(tree.d2, nδt), pet, nδt, nfdt, pv)
end



"""
    extractp(tree::Nothing, δt::Float64, lv::Function)

Log-linearly predict Geometric Brownian motion for `lv` at times given by `δt`.
"""
extractp(tree::Nothing, nδt::Float64) = nothing



"""
    extractp(tree::iTgbmbd, nδt::Float64)

Log-linearly predict Geometric Brownian motion for `lλ` and `lμ` 
at times given by `nδt`.
"""
function extractp(tree::iTgbmbd, nδt::Float64)

  pet = pe(tree)
  δt  = dt(tree)
  vλ  = lλ(tree)
  vμ  = lμ(tree)
  n   = floor(Int64, pet/nδt)

  pλv = Float64[]
  pμv = Float64[]
  for i in Base.OneTo(n+1)
    iti = Float64(i-1)*nδt
    ix  = floor(Int64, iti/δt) + 1
    tts = δt*Float64(ix-1)
    ttf = δt*Float64(ix)
    push!(pλv, linpred(iti, tts, ttf, vλ[ix], vλ[ix+1]))
    push!(pμv, linpred(iti, tts, ttf, vμ[ix], vμ[ix+1]))
  end

  if (Float64(n+1)*nδt) !== pet
    push!(pλv, vλ[end])
    push!(pμv, vμ[end])
    nfdt = mod(pet,nδt)
  else
    nfdt = nδt
  end

  iTgbmbd(extractp(tree.d1, nδt), 
          extractp(tree.d2, nδt), pet, nδt, nfdt, false, false, pλv, pμv)
end

"""
    extractp(tree::Nothing, δt::Float64)

Log-linearly predict Geometric Brownian motion for `lλ` and `lμ` 
at times given by `nδt`.
"""
extractp(tree::Nothing, δt::Float64) = nothing



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
function iquantile(treev::Array{iTgbmpb,1}, 
                   p    ::Float64)

  nt = lastindex(treev)
  t1 = treev[1]

  # make vector of lambdas
  vs = Array{Float64,1}[]
  for t in treev
    push!(vs, lλ(t)) 
  end

  sv = Float64[]
  # make fill vector to estimate statistics
  v = Array{Float64,1}(undef, nt)
  for i in Base.OneTo(lastindex(vs[1]))
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

  iTgbmpb(iquantile(treev1, p), iquantile(treev2, p),
    pe(t1), dt(t1), fdt(t1), sv)
end

"""
    iquantile(::Nothing, p::Float64, lv::Function)

Make an `iTgbmpb` with the quantile specified by `p` in data specified in 
function `lv`.
"""
iquantile(::Nothing, p::Float64) = nothing




"""
    iquantile(treev::Array{T,1}, 
              p    ::Float64, 
              lv   ::Function) where {T <: iTgbm}

Make an `iTgbmbd` with the quantile specified by `p` in data specified in 
function `lv`.
"""
function iquantile(treev::Array{iTgbmbd,1}, 
                   p    ::Float64)

  nt  = lastindex(treev)

  t1 = treev[1]

  # make vector of lambdas and mus
  vsλ = Array{Float64,1}[]
  vsμ = Array{Float64,1}[]
  for t in treev
    push!(vsλ, lλ(t)) 
    push!(vsμ, lμ(t)) 
  end

  svλ = Float64[]
  svμ = Float64[]
  # make fill vector to estimate statistics
  vλ = Array{Float64,1}(undef, nt)
  vμ = Array{Float64,1}(undef, nt)
  for i in Base.OneTo(lastindex(vsλ[1]))
    for t in Base.OneTo(nt)
      vλ[t] = vsλ[t][i]
      vμ[t] = vsμ[t][i]
    end
    push!(svλ, quantile(vλ, p))
    push!(svμ, quantile(vμ, p))
  end

  if isnothing(treev[1].d1)
    treev1 = nothing
  else
    treev1 = iTgbmbd[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end 
  end

  if isnothing(treev[1].d2)
    treev2 = nothing
  else
    treev2 = iTgbmbd[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end 
  end

  iTgbmbd(iquantile(treev1, p), 
          iquantile(treev2, p), 
          pe(t1), dt(t1), fdt(t1), false, false, svλ, svμ)
end

"""
    iquantile(::Nothing, p::Float64)

Make an `iTgbmbb` with the quantile specified by `p` in data specified in 
function `lv`.
"""
iquantile(::Nothing, p::Float64) = nothing


