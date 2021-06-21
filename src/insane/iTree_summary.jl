#=

Summary functions for vector of `iTree`

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    mcmc_array(treev::Array{iTgbmpb,1}, δt::Float64)

Return an Array with a row for each sampled tree for interpolated
`λ` parameters at times determined by `δt`.
"""
function mcmc_array(treev::Array{iTgbmpb,1}, δt::Float64)

  @inbounds begin

    vi = Float64[]
    extract_vector!(treev[1], vi, δt, 0.0)

    ra = Array{Float64,2}(undef, lastindex(treev), lastindex(vi))
    ra[1,:] = vi

    for i in 2:lastindex(treev)
      vi = Float64[]
      extract_vector!(treev[i], vi, δt, 0.0)
      ra[i, :] = vi
    end
  end

  return ra
end




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
    extract_tree(tree::iTgbmpb, nδt::Float64)

Log-linearly predict Geometric Brownian motion for `lv` at times given by `nδt`,
and return a tree.
"""
function extract_tree(tree::iTgbmpb, nδt::Float64)

  pet = pe(tree)
  δt  = dt(tree)
  v   = lλ(tree)
  n   = floor(Int64, pet/nδt)
  i1  = isapprox(pet, Float64(n)*nδt, atol = 1e-11) ? 0 : 1

  pv = Float64[]
  for i in Base.OneTo(n+i1)
    iti = Float64(i-1)*nδt
    ix  = floor(Int64, iti/δt) + 1
    tts = δt*Float64(ix-1)
    ttf = δt*Float64(ix)
    push!(pv, linpred(iti, tts, ttf, v[ix], v[ix+1]))
  end

  push!(pv, v[end])

  if isapprox(pet, Float64(n)*nδt, atol = 1e-11)
    nfdt = nδt
  else
    nfdt = mod(pet,nδt)
  end

  iTgbmpb(extract_tree(tree.d1, nδt), 
          extract_tree(tree.d2, nδt), pet, nδt, nfdt, pv)
end



"""
    extractp(tree::Nothing, δt::Float64, lv::Function)

Log-linearly predict Geometric Brownian motion for `lv` at times given by `δt`.
"""
extract_tree(tree::Nothing, nδt::Float64) = nothing




"""
    extract_vector!(tree::iTgbmpb,
                   v   ::Array{Float64,1}, 
                   nδt ::Float64, 
                   ct  ::Float64)

Log-linearly predict Geometric Brownian motion for `λ` at times given by `nδt`
and return a vector.
"""
function extract_vector!(tree::iTgbmpb, 
                         v   ::Array{Float64,1}, 
                         nδt ::Float64, 
                         ct  ::Float64)

  pet = pe(tree)
  δt  = dt(tree)
  vt  = lλ(tree)
  n   = floor(Int64, (pet - ct)/nδt)
  i1  = isapprox(pet - ct, Float64(n)*nδt, atol = 1e-11) ? 0 : 1

  pv = Float64[]
  iti = ct
  for i in Base.OneTo(n+i1)
    ix  = fld(iti,δt)
    tts = δt *  ix
    ttf = δt * (ix + 1.0)
    Ix  = Int64(ix)
    push!(pv, linpred(iti + 0.5*δt, tts, ttf, vt[Ix+1], vt[Ix+2]))
    iti = Float64(i)*nδt + ct
  end

  append!(v, pv)

  if !isnothing(tree.d1)
    extract_vector!(tree.d1::iTgbmpb, v, nδt, max(0.0, iti - pet))
  end
  if !isnothing(tree.d2)
    extract_vector!(tree.d2::iTgbmpb, v, nδt, max(0.0, iti - pet))
  end
end




"""
    extractp(tree::iTgbmbd, nδt::Float64)

Log-linearly predict Geometric Brownian motion for `lλ` and `lμ` 
at times given by `nδt`.
"""
function extract_tree(tree::iTgbmbd, nδt::Float64)

  pet = pe(tree)
  δt  = dt(tree)
  vλ  = lλ(tree)
  vμ  = lμ(tree)
  n   = floor(Int64, pet/nδt)

  i1  = isapprox(pet, Float64(n)*nδt, atol = 1e-11) ? 0 : 1
  pλv = Float64[]
  pμv = Float64[]
  for i in Base.OneTo(n+i1)
    iti = Float64(i-1)*nδt
    ix  = floor(Int64, iti/δt) + 1
    tts = δt*Float64(ix-1)
    ttf = δt*Float64(ix)
    push!(pλv, linpred(iti, tts, ttf, vλ[ix], vλ[ix+1]))
    push!(pμv, linpred(iti, tts, ttf, vμ[ix], vμ[ix+1]))
  end

  push!(pλv, vλ[end])
  push!(pμv, vμ[end])

  if isapprox(pet, Float64(n)*nδt, atol = 1e-11)
    nfdt = nδt
  else
    nfdt = mod(pet,nδt)
  end

  iTgbmbd(extract_tree(tree.d1, nδt), 
          extract_tree(tree.d2, nδt), pet, nδt, nfdt, false, false, pλv, pμv)
end

"""
    extract_tree(tree::Nothing, δt::Float64)

Log-linearly predict Geometric Brownian motion for `lλ` and `lμ` 
at times given by `nδt`.
"""
extract_tree(tree::Nothing, δt::Float64) = nothing



"""
    linpred(val::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64)

Estimate `y` at `x = val` according to linear interpolation for a range.
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


