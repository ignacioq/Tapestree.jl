#=

Summary functions for vector of `iTree`

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    time_quantile(r::Array{Array{Float64,1},1}, p::Array{Float64,1})

Return quantiles given by probabilities `p` vector.
"""
function time_quantile(r::Array{Array{Float64,1},1}, p::Array{Float64,1})

  fx = Array{Float64,2}(undef, lastindex(r), lastindex(p))
  for (i, ri) in enumerate(r)
    fx[i,:] = quantile(ri, p)
  end
  return fx
end




"""
    time_fun(f::Function, r::Array{Array{Float64,1},1})

Apply function `f` to each vector in `r` and return the value.
"""
function time_fun(f::Function, r::Array{Array{Float64,1},1})
  lr = lastindex(r)
  fx = Array{Float64,1}(undef, lr)
  for i in Base.OneTo(lr)
    fx[i] = f(r[i])
  end
  return fx
end




"""
    time_rate(tree::T, f::Function, δt::Float64) where {T <: iT}

Extract values from `f` function at times sampled every `δt` across the tree.
"""
function time_rate(tree::T, f::Function, δt::Float64) where {T <: iTree}

  th = treeheight(tree)

  # make time vector (present = 0.0)
  ts  = [0.0:δt:th...]
  reverse!(ts)

  # make vector of vector to push rates
  r = [Float64[] for i in Base.OneTo(lastindex(ts))]

  # do recursive
  _time_rate!(tree, ts, δt, r, 1, th, f)
  pop!(r)
  pop!(ts)

  return ts, r
end




"""
    _time_rate!(tree::T,
               ts  ::Array{Float64,1},
               tdt ::Float64,
               r   ::Array{Array{Float64,1},1},
               tii ::Int64,
               ct  ::Float64,
               f  ::Function) where {T <: iT}

Extract values from `f` function at times `ts` across the tree.
"""
@inline function _time_rate!(tree::T,
                             ts  ::Array{Float64,1},
                             tdt ::Float64,
                             r   ::Array{Array{Float64,1},1},
                             tii ::Int64,
                             ct  ::Float64,
                             f   ::Function) where {T <: iTree}
  if ct < accerr
    return nothing
  end

  et = e(tree)
  δt = dt(tree)
  vt = f(tree)

  nts, re = divrem(et - (ct - ts[tii]), tdt, RoundDown)
  nts = max(-1, Int64(nts) - Int64(re < accerr))
  tsr = tii:(tii + nts)

  # have to match ts times to f vector
  @simd for i in tsr
    tsi = ts[i]
    bt  = ct - tsi
    ix  = max(0.0, div(bt, δt, RoundDown))
    tts = δt *  ix
    ttf = δt * (ix + 1.0)
    Ix  = Int64(ix) + 1
    push!(r[i], linpred(bt, tts, ttf, vt[Ix], vt[Ix+1]))
  end

  if def1(tree)
    _time_rate!(tree.d1, ts, tdt, r, tii + nts + 1, ct - et, f)
    if def2(tree)
      _time_rate!(tree.d2, ts, tdt, r, tii + nts + 1, ct - et, f)
    end
  end

  return nothing
end



"""
    aggr_rates(tree::T,
               f   ::Function,
               δt  ::Float64;
               af = x -> quantile(x, 0.5)) where {T <: iT}

Return the aggregated rates from function `f` by `af`.
"""
function aggr_rates(tree::T,
                    f   ::Function,
                    δt  ::Float64;
                    af = x -> quantile(x, 0.5)) where {T <: iT}

  # prepare data
  ts, r = time_rate(tree, f, δt)
  lts = lastindex(ts)
  
  m = Array{Float64,1}(undef, lts)
  for (i, ri) in enumerate(r)
    m[i] = af(ri)
  end

  return ts, m
end




"""
    aggr_rates(trees::Vector{T},
               f    ::Function,
               δt   ::Float64;
               af = x -> quantile(x, 0.5)) where {T <: iT}

Return the aggregated rates from function `f` by `af`.
"""
function aggr_rates(trees::Vector{T},
                    f    ::Function,
                    δt   ::Float64;
                    af = x -> quantile(x, 0.5)) where {T <: iT}

  ntrees = lastindex(trees)
  riv = Vector{Float64}[]

  lts = 0
  ts  = Float64[]
  for t in trees

    # tree extracting function
    tsi, ri = time_rate(t, f, δt)

    # aggregating function
    ri = map(af, ri)

    if lastindex(tsi) > lts
      ts  = tsi
      reverse!(ts)
      lts = lastindex(tsi)
    end

    reverse!(ri)
    push!(riv, ri)
  end

  # estimate quantiles
  q = fill(NaN, lts, ntrees)
  # estimate quantiles
  for i in Base.OneTo(ntrees)
    ri = riv[i]
    lr = lastindex(ri)
    q[1:lr,i] = ri
  end

  m = Array{Float64}(undef, lts)
  for i in Base.OneTo(lts)
    qi = q[i,:]
    filter!(!isnan, qi)
    m[i] = quantile(qi, 0.5)
  end

  return ts, m
end




"""
    aggr_ltt(nts::Vector{Ltt}, 
                  tdt::Float64;
                  af = x -> quantile(x, 0.5))

Return the aggregated ltt by `af`.
"""
function aggr_ltt(nts::Vector{Ltt}, 
                  tdt::Float64;
                  af = x -> quantile(x, 0.5))

  n  = lastindex(nts)
  th = maximum(map(x -> maximum(x.t), nts))

  # make time vector (present = 0.0)
  ts  = [0.0:tdt:th...]
  lts = length(ts)

  q = zeros(Int64, lts, n)
  for j in Base.OneTo(n), i in Base.OneTo(lts)
    q[i,j] = nspt(nts[j], ts[i])
  end

  m = Array{Float64}(undef, lts)
  for i in Base.OneTo(lts)
    qi = q[i,:]
    filter!(!isnan, qi)
    m[i] = af(qi)
  end

  return ts, m
end




"""
    sample(treev::Array{T,1},
           f    ::Function,
           δt   ::Float64) where {T <: iT}

Return an Array with a row for each sampled tree for interpolated
parameters accessed by `f` at times determined by `δt`.
"""
function sample(tree::T,
                f   ::Function,
                δt  ::Float64) where {T <: iTree}
  vi = Float64[]
  extract_vector!(tree, vi, δt, 0.0, f)

  return vi
end




"""
    sample(treev::Array{T,},
           f    ::Function,
           δt   ::Float64) where {T <: cT}

Return an Array with a row for attributes accessed by `f` at each branch.
"""
function sample(tree::T, f::Function) where {T <: cT}

  vi = Float64[]
  extract_vector!(tree, vi, δt, 0.0, f)

  return vi
end




"""
    sample(treev::Array{T,1},
           f    ::Function,
           δt   ::Float64) where {T <: iTree}

Return an Array with a row for each sampled tree for interpolated
parameters accessed by `f` at times determined by `δt`.
"""
function sample(treev::Array{T,1},
                f    ::Function,
                δt   ::Float64) where {T <: iTree}

  @inbounds begin

    vi = Float64[]
    extract_vector!(treev[1], vi, δt, 0.0, f)

    ra = Array{Float64,2}(undef, lastindex(treev), lastindex(vi))
    ra[1,:] = vi

    for i in 2:lastindex(treev)
      vi = Float64[]
      extract_vector!(treev[i], vi, δt, 0.0, f)
      ra[i, :] = vi
    end
  end

  return ra
end




"""
    extract_vector!(tree::T,
                    v   ::Array{Float64,1},
                    nδt ::Float64,
                    ct  ::Float64,
                    f  ::Function) where {T <: iTree}

Log-linearly predict Geometric Brownian motion for `λ` at times given by `nδt`
and return a vector.
"""
function extract_vector!(tree::T,
                         v   ::Array{Float64,1},
                         nδt ::Float64,
                         ct  ::Float64,
                         f   ::Function) where {T <: iT}

  et = e(tree)
  δt = dt(tree)
  vt = f(tree)
  n  = floor(Int64, (et - ct)/nδt)
  i1 = isapprox(et - ct, Float64(n)*nδt, atol = 1e-11) ? 0 : 1

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

  if def1(tree)
    extract_vector!(tree.d1::T, v, nδt, max(0.0, iti - et), f)
    if def2(tree)
      extract_vector!(tree.d2::T, v, nδt, max(0.0, iti - et), f)
    end
  end
end





"""
    linearize_gbm(tree::T, f::Function) where {T <: iT}

Extract the parameters given by `f` into a linear Array.
"""
function linearize_gbm(tree::T, f::Function) where {T <: iT}
  v = Float64[]
  _linearize_gbm!(tree, f, v)
  return v
end




"""
    _linearize_gbm!(tree::T,
                    f  ::Function,
                    v   ::Array{Float64,1}) where {T <: iT}

Extract the parameters given by `f` into a linear Array, initialized with an
array `v`.
"""
function _linearize_gbm!(tree::T,
                         f  ::Function,
                         v   ::Array{Float64,1}) where {T <: iT}

  append!(v, f(tree))
  if def1(tree)
    _linearize_gbm!(tree.d1::T, f, v)
    if def2(tree)
      _linearize_gbm!(tree.d2::T, f, v)
    end
  end
end




"""
    extract_tree(tree::iTb, nδt::Float64)

Log-linearly predict Geometric Brownian motion for `f` at times given by `nδt`,
and return a tree.
"""
function extract_tree(tree::iTb, nδt::Float64)

  et = e(tree)
  δt = dt(tree)
  v  = lλ(tree)
  n  = floor(Int64, et/nδt)
  i1 = isapprox(et, Float64(n)*nδt, atol = 1e-11) ? 0 : 1

  pv = Float64[]
  for i in Base.OneTo(n+i1)
    iti = Float64(i-1)*nδt
    ix  = floor(Int64, iti/δt) + 1
    tts = δt*Float64(ix-1)
    ttf = δt*Float64(ix)
    push!(pv, linpred(iti, tts, ttf, v[ix], v[ix+1]))
  end

  push!(pv, v[end])

  if isapprox(et, Float64(n)*nδt, atol = 1e-11)
    nfdt = nδt
  else
    nfdt = mod(et,nδt)
  end

  if def1(tree)
    iTb(extract_tree(tree.d1, nδt),
         extract_tree(tree.d2, nδt), et, nδt, nfdt, pv)
  else
    iTb(et, nδt, nfdt, pv)
  end
end




"""
    extractp(tree::iTbd, nδt::Float64)

Log-linearly predict Geometric Brownian motion for `lλ` and `lμ`
at times given by `nδt`.
"""
function extract_tree(tree::iTbd, nδt::Float64)

  et = e(tree)
  δt = dt(tree)
  vλ = lλ(tree)
  vμ = lμ(tree)
  n  = floor(Int64, et/nδt)

  i1  = isapprox(et, Float64(n)*nδt, atol = 1e-11) ? 0 : 1
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

  if isapprox(et, Float64(n)*nδt, atol = 1e-11)
    nfdt = nδt
  else
    nfdt = mod(et,nδt)
  end

  if def1(tree)
    iTbd(extract_tree(tree.d1, nδt),
            extract_tree(tree.d2, nδt), et, nδt, nfdt, false, false, pλv, pμv)
  else
    iTbd(et, nδt, nfdt, false, false, pλv, pμv)
  end
end




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
    iquantile(treev::Array{iTce,1}, p::Float64)

Make an `iTce` with the quantile specified by `p` in data specified in
function `f`.
"""
function iquantile(treev::Array{iTce,1}, p::Float64)

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

  if isdefined(t1, :d1)
    treev1 = iTce[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = iTce[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    iTce(iquantile(treev1, p), iquantile(treev2, p),
      e(t1), dt(t1), fdt(t1), isextinct(t1), true, sv)
  else
    iTce(e(t1), dt(t1), fdt(t1), isextinct(t1), true, sv)
  end
end




"""
    iquantile(treev::Array{iTct,1}, p::Float64)

Make an `iTct` with the quantile specified by `p` in data specified in
function `f`.
"""
function iquantile(treev::Vector{iTct}, p::Float64)

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

  if isdefined(t1,:d1)
    treev1 = iTct[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = iTct[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    iTct(iquantile(treev1, p), iquantile(treev2, p),
      e(t1), dt(t1), fdt(t1), isextinct(t1), true, sv)
  else
    iTct(e(t1), dt(t1), fdt(t1), isextinct(t1), true, sv)
  end

end




"""
    iquantile(reev::Vector{iTb}, p::Float64)

Make an `iTb` with the quantile specified by `p`.
"""
function iquantile(treev::Vector{iTb}, p::Float64)

  nts = lastindex(treev)
  t1  = treev[1]

  # make vector of lambdas
  vs = Array{Float64,1}[]
  for t in treev
    push!(vs, lλ(t))
  end

  sv = Float64[]
  # make fill vector to estimate statistics
  v = Array{Float64,1}(undef, nts)
  for i in Base.OneTo(lastindex(vs[1]))
    for t in Base.OneTo(nts)
      v[t] = vs[t][i]
    end
    push!(sv, quantile(v, p))
  end

  if isdefined(t1, :d1)
    treev1 = iTb[]
    for t in Base.OneTo(nts)
        push!(treev1, treev[t].d1)
    end
    treev2 = iTb[]
    for t in Base.OneTo(nts)
        push!(treev2, treev[t].d2)
    end
    iTb(iquantile(treev1, p), iquantile(treev2, p),
         e(t1), dt(t1), fdt(t1), true, sv)
  else
    iTb(e(t1), dt(t1), fdt(t1), true, sv)
  end
end




"""
    iquantile(treev::Vector{iTbd}, p::Float64)

Make an `iTbd` with the quantile specified by `p` in data specified in
function `f`.
"""
function iquantile(treev::Vector{iTbd}, p::Float64)

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

  if def1(t1)
    treev1 = iTbd[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = iTbd[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    iTbd(iquantile(treev1, p),
         iquantile(treev2, p),
         e(t1), dt(t1), fdt(t1), isextinct(t1), true, svλ, svμ)
  else
    iTbd(e(t1), dt(t1), fdt(t1), isextinct(t1), true, svλ, svμ)
  end
end




"""
    iquantile(treev::Vector{iTfbd}, p::Float64)

Make an `iTfbd` with the quantile specified by `p` in data specified in
function `f`.
"""
function iquantile(treev::Vector{iTfbd}, p::Float64)

  nts = lastindex(treev)
  t1  = treev[1]

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
  vλ = Array{Float64,1}(undef, nts)
  vμ = Array{Float64,1}(undef, nts)
  for i in Base.OneTo(lastindex(vsλ[1]))
    for t in Base.OneTo(nts)
      if lastindex(vsλ[t]) != lastindex(vsλ[1])
        t2 = treev[t]
      end

      vλ[t] = vsλ[t][i]
      vμ[t] = vsμ[t][i]
    end
    push!(svλ, quantile(vλ, p))
    push!(svμ, quantile(vμ, p))
  end

  if def1(t1)
    treev1 = iTfbd[]
    for t in Base.OneTo(nts)
        push!(treev1, treev[t].d1)
    end

    if def2(t1)
      treev2 = iTfbd[]
      for t in Base.OneTo(nts)
          push!(treev2, treev[t].d2)
      end

      iTfbd(iquantile(treev1, p),
            iquantile(treev2, p),
            e(t1), dt(t1), fdt(t1), false, false, false, svλ, svμ)
    else
      iTfbd(iquantile(treev1, p),
            e(t1), dt(t1), fdt(t1), false, true, false, svλ, svμ)
    end
  else
    iTfbd(e(t1), dt(t1), fdt(t1), false, isfossil(t1), false, svλ, svμ)
  end
end




"""
    imean(treev::Vector{sTpe})

Make a mean `sTpe`.
"""
function imean(treev::Vector{T}) where {T <: Tpe}

  nts = lastindex(treev)
  t1  = treev[1]

  n   = _count_nodes!(t1, 0)
  xiv = zeros(Float64, n)
  xfv = zeros(Float64, n)

  for t in treev
    _sum_xv!(t, 0, xiv, xfv)
  end

  @turbo xiv ./= Float64(nts)
  @turbo xfv ./= Float64(nts)

  i, tree = make_tm(t1, 0, xiv, xfv)

  return tree
end



"""
    make_tm(tree::sTpe, 
            i   ::Int64, 
            xiv::Vector{Float64}, 
            xfv::Vector{Float64})

Make mean tree from node data.
"""
function make_tm(tree::sTpe, 
                 i   ::Int64, 
                 xiv::Vector{Float64}, 
                 xfv::Vector{Float64})

  i += 1
  xii, xfi = xiv[i], xfv[i]

  if def1(tree)
     i, d1 = make_tm(tree.d1, i, xiv, xfv)
     i, d2 = make_tm(tree.d2, i, xiv, xfv)

     return i, sTpe(d1, d2, e(tree), isextinct(tree), xii, xfi, sh(tree), true)
  else
    return i, sTpe(e(tree), isextinct(tree), xii, xfi, sh(tree), true)
  end
end




"""
    make_tm(tree::sfTpe, 
            i   ::Int64, 
            xiv::Vector{Float64}, 
            xfv::Vector{Float64})

Make mean tree from node data.
"""
function make_tm(tree::sTfpe, 
                 i   ::Int64, 
                 xiv::Vector{Float64}, 
                 xfv::Vector{Float64})

  i += 1
  xii, xfi = xiv[i], xfv[i]

  if def1(tree)
    i, d1 = make_tm(tree.d1, i, xiv, xfv)
    if def2(tree)
      i, d2 = make_tm(tree.d2, i, xiv, xfv)
      return i, sTfpe(d1, d2, e(tree), isextinct(tree), isfossil(tree), xii, xfi, sh(tree), true)
    else
      return i, sTfpe(d1, e(tree), isextinct(tree), isfossil(tree), xii, xfi, sh(tree), true)
    end
  else
    return i, sTfpe(e(tree), isextinct(tree), isfossil(tree), xii, xfi, sh(tree), true)
  end
end




"""
    _count_nodes!(tree::T, i::Int64) where {T <: Tpe}

Count number of nodes
"""
function _count_nodes!(tree::T, i::Int64) where {T <: Tpe}
  i += 1

  if def1(tree)
    i = _count_nodes!(tree.d1, i)
    if def2(tree)
      i = _count_nodes!(tree.d2, i)
    end
  end

  return i
end




"""
    _sum_xv!(tree::T, 
             i   ::Int64, 
             xiv ::Vector{Float64},
             xfv ::Vector{Float64}) where {T <: Tpe}

Make sum of xi and xf for each node.
"""
function _sum_xv!(tree::T, 
                  i   ::Int64, 
                  xiv ::Vector{Float64},
                  xfv ::Vector{Float64}) where {T <: Tpe}

  i += 1

  xiv[i] += xi(tree)
  xfv[i] += xf(tree)

  if def1(tree)
    i = _sum_xv!(tree.d1, i, xiv, xfv)
    if def2(tree)
      i = _sum_xv!(tree.d2, i, xiv, xfv)
    end
  end

  return i
end



"""
    imean(treev::Vector{cTb})

Make an `cTb` with the geometric mean.
"""
function imean(treev::Vector{cTb})

  nt  = lastindex(treev)

  t1 = treev[1]

  # make vector of lambdas
  ss = 0.0
  for t in treev
    ss += lλ(t)
  end
  mλ = ss/Float64(nt)

  if def1(t1)
    treev1 = cTb[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = cTb[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    cTb(imean(treev1),
        imean(treev2),
        e(t1), true, mλ)
  else
    cTb(e(t1), true, mλ)
  end
end




"""
    imean(treev::Vector{T}) where {T <: cT}

Make an `cTce` & `cTct` with the geometric mean.
"""
function imean(treev::Vector{T}) where {T <: cT}

  nt  = lastindex(treev)

  t1 = treev[1]

  # make vector of lambdas
  ss = 0.0
  for t in treev
    ss += lλ(t)
  end
  mλ = ss/Float64(nt)

  if def1(t1)
    treev1 = T[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = T[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    T(imean(treev1),
         imean(treev2),
         e(t1), isextinct(t1), true, mλ)
  else
    T(e(t1), isextinct(t1), true, mλ)
  end
end




"""
    imean(treev::Vector{cTbd})

Make an `cTbd` with the geometric mean.
"""
function imean(treev::Vector{cTbd})

  nt  = lastindex(treev)

  t1 = treev[1]

  # make vector of lambdas
  ssλ = ssμ = 0.0
  for t in treev
    ssλ += lλ(t)
    ssμ += lμ(t)
  end
  mλ = ssλ/Float64(nt)
  mμ = ssμ/Float64(nt)

  if def1(t1)
    treev1 = cTbd[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = cTbd[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    cTbd(imean(treev1), imean(treev2),
         e(t1), isextinct(t1), true, mλ, mμ)
  else
    cTbd(e(t1), isextinct(t1), true, mλ, mμ)
  end
end




"""
    imean(treev::Vector{cTfbd})

Make an `cTfbd` with the geometric mean.
"""
function imean(treev::Vector{cTfbd})

  nt  = lastindex(treev)

  t1 = treev[1]

  # make vector of lambdas
  ssλ = ssμ = 0.0
  for t in treev
    ssλ += lλ(t)
    ssμ += lμ(t)
  end
  mλ = ssλ/Float64(nt)
  mμ = ssμ/Float64(nt)

  if def1(t1)
    treev1 = cTfbd[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    if def2(t1)

      treev2 = cTfbd[]
      for t in Base.OneTo(nt)
          push!(treev2, treev[t].d2)
      end

      cTfbd(imean(treev1), imean(treev2),
           e(t1), isextinct(t1), isfossil(t1), true, mλ, mμ)
    else
      cTfbd(imean(treev1),
           e(t1), isextinct(t1), isfossil(t1), true, mλ, mμ)
    end
  else
    cTfbd(e(t1), isextinct(t1), isfossil(t1), true, mλ, mμ)
  end
end




"""
    imean(treev::Vector{iTb})

Make an `iT` with the geometric mean.
"""
function imean(treev::Vector{iTb})

  nt  = lastindex(treev)

  t1 = treev[1]

  # make vector of lambdas and mus
  vsλ = Array{Float64,1}[]
  for t in treev
    push!(vsλ, lλ(t))
  end

  svλ = Float64[]
  # make fill vector to estimate statistics
  vλ = Array{Float64,1}(undef, nt)
  for i in Base.OneTo(lastindex(vsλ[1]))
    for t in Base.OneTo(nt)
      vλ[t] = vsλ[t][i]
    end
    push!(svλ, mean(vλ))
  end

  if def1(t1)
    treev1 = iTb[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = iTb[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    iTb(imean(treev1),
         imean(treev2),
         e(t1), dt(t1), fdt(t1), true, svλ)
  else
    iTb(e(t1), dt(t1), fdt(t1), true, svλ)
  end
end





"""
    imean(treev::Vector{T})

Make an `iT` with the geometric mean.
"""
function imean(treev::Vector{T}) where {T <: iT}

  nt  = lastindex(treev)

  t1 = treev[1]

  # make vector of lambdas and mus
  vsλ = Array{Float64,1}[]
  for t in treev
    push!(vsλ, lλ(t))
  end

  svλ = Float64[]
  # make fill vector to estimate statistics
  vλ = Array{Float64,1}(undef, nt)
  for i in Base.OneTo(lastindex(vsλ[1]))
    for t in Base.OneTo(nt)
      vλ[t] = vsλ[t][i]
    end
    push!(svλ, mean(vλ))
  end

  if def1(t1)
    treev1 = T[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = T[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    T(imean(treev1),
      imean(treev2),
      e(t1), dt(t1), fdt(t1), isextinct(t1), true, svλ)
  else
    T(e(t1), dt(t1), fdt(t1), isextinct(t1), true, svλ)
  end
end



"""
    imean(treev::Vector{iTbd})

Make an `iTbd` with the geometric mean.
"""
function imean(treev::Vector{iTbd})

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
    push!(svλ, mean(vλ))
    push!(svμ, mean(vμ))
  end

  if def1(t1)
    treev1 = iTbd[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end
    treev2 = iTbd[]
    for t in Base.OneTo(nt)
        push!(treev2, treev[t].d2)
    end

    iTbd(imean(treev1),
         imean(treev2),
         e(t1), dt(t1), fdt(t1), isextinct(t1), true, svλ, svμ)
  else
    iTbd(e(t1), dt(t1), fdt(t1), isextinct(t1), true, svλ, svμ)
  end
end





"""
    imean(treev::Vector{iTfbd})

Make an `iTfbd` with the geometric mean.
"""
function imean(treev::Vector{iTfbd})

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
    push!(svλ, mean(vλ))
    push!(svμ, mean(vμ))
  end

  if def1(t1)
    treev1 = iTfbd[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end

    if def2(t1)
      treev2 = iTfbd[]
      for t in Base.OneTo(nt)
          push!(treev2, treev[t].d2)
      end

      iTfbd(imean(treev1),
            imean(treev2),
            e(t1), dt(t1), fdt(t1), isextinct(t1), isfossil(t1), true, svλ, svμ)
    else
      iTfbd(imean(treev1),
            e(t1), dt(t1), fdt(t1), isextinct(t1), isfossil(t1), true, svλ, svμ)
    end
  else
    iTfbd(e(t1), dt(t1), fdt(t1), isextinct(t1), isfossil(t1), true, svλ, svμ)
  end
end






"""
    imean(treev::Vector{sTxs})

Make an `sTxs` with the geometric mean.
"""
function imean(treev::Vector{sTxs})

  nt  = lastindex(treev)

  t1 = treev[1]

  # make vector of lambdas and mus
  vsx = Array{Float64,1}[]
  vsσ = Array{Float64,1}[]
  for t in treev
    push!(vsx, xv(t))
    push!(vsσ, lσ2(t))
  end

  svx = Float64[]
  svσ = Float64[]
  # make fill vector to estimate statistics
  vx = Array{Float64,1}(undef, nt)
  vσ = Array{Float64,1}(undef, nt)
  for i in Base.OneTo(lastindex(vsx[1]))
    for t in Base.OneTo(nt)
      vx[t] = vsx[t][i]
      vσ[t] = vsσ[t][i]
    end
    push!(svx, mean(vx))
    push!(svσ, mean(vσ))
  end

  if def1(t1)
    treev1 = sTxs[]
    for t in Base.OneTo(nt)
        push!(treev1, treev[t].d1)
    end

    if def2(t1)
      treev2 = sTxs[]
      for t in Base.OneTo(nt)
          push!(treev2, treev[t].d2)
      end

      sTxs(imean(treev1),
           imean(treev2),
           e(t1), dt(t1), fdt(t1), svx, svσ)
    else
      sTxs(imean(treev1),
           e(t1), dt(t1), fdt(t1), svx, svσ)
    end
  else
    sTxs(e(t1), dt(t1), fdt(t1), svx, svσ)
  end
end




