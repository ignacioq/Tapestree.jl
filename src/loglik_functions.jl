#=

likelihood functions tribe model

Ignacio Quintero Mächler

t(-_-t)

May 15 2017

=#





"""
    Eδx(μ::Float64, ωx::Float64, δt::Float64)

Return the expected value given the weighted average with sympatric species.
"""
Eδx(μ::Float64, ωx::Float64, δt::Float64) = @fastmath ωx * μ * δt





"""
    f_λ1(λ1::Float64, ω1::Float64, δx::Float64)

Estimate rates for area colonization based 
on the difference between lineage traits and area averages.
"""
function f_λ1(λ1::Float64, ω1::Float64, δx::Float64)
  if iszero(δx) 
    return λ1
  else
    return λ1 * (1.0 - ω1/(1.0 + δx))
  end
end





"""
    f_λ0(λ0::Float64, ω0::Float64, δx::Float64)

Estimate rates for area colonization/loss based 
on the difference between lineage traits and area averages.
"""
function f_λ0(λ0::Float64, ω0::Float64, δx::Float64) 
  if iszero(δx)
    return λ0
  else
    return λ0 * (1.0 + ω0/(1.0 + δx))
  end
end





"""
    makellf(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64, narea::Int64)

Make likelihood and likelihood ratio functions 
for all trait matrix and biogeography history.
"""
function makellf(δt   ::Array{Float64,1}, 
                 Y    ::Array{Int64,3}, 
                 ntip ::Int64, 
                 narea::Int64,
                 m    ::Int64,
                 nedge::Int64)

  # get initial range
  const wf23 = Int64[]
  for j = Base.OneTo(ntip)
    push!(wf23, findfirst(Y[:,j,1] .!= 23))
  end

  # number of normal evaluations
  n = 0
  for i = wf23
    n += length(i:(m-1))
  end

  # normal constant
  const normC = -0.5*log(2.0π)*n

  function llf(X     ::Array{Float64,2},
               Y     ::Array{Int64,3}, 
               LA    ::Array{Float64,2},
               LD    ::Array{Float64,3},
               ωx    ::Float64,
               ω1    ::Float64,
               ω0    ::Float64,
               λ1    ::Float64,
               λ0    ::Float64,
               stemev::Vector{Vector{Float64}},
               brs   ::Array{Int64,3},
               σ²    ::Float64)

    ll::Float64 = normC

    @inbounds begin

      # trait likelihood
      for j = Base.OneTo(ntip)
        @simd for i = wf23[j]:(m-1)
          ll += logdnorm_tc(X[(i+1),j], 
                            X[i,j] + Eδx(LA[i,j], ωx, δt[i]), 
                            δt[i]*σ²)::Float64
        end
      end

      # biogeograhic likelihood
      for k = Base.OneTo(narea)
        ll += brll(stemev[k], λ1, λ0, brs[nedge,1,k])::Float64
        for j = Base.OneTo(ntip)
          ll += bitvectorll(Y, λ1, λ0, ω1, ω0, LD, δt, 
                            j, k, wf23[j], m)::Float64
        end
      end

    end

    return ll::Float64
  end

  function llrf(Xc     ::Array{Float64,2},
                Xp     ::Array{Float64,2},
                Yc     ::Array{Int64,3}, 
                Yp     ::Array{Int64,3}, 
                LAc    ::Array{Float64,2},
                LAp    ::Array{Float64,2},
                LDc    ::Array{Float64,3},
                LDp    ::Array{Float64,3},
                ωx     ::Float64,
                ω1     ::Float64,
                ω0     ::Float64,
                λ1     ::Float64,
                λ0     ::Float64,
                stemevc::Vector{Vector{Float64}},
                stemevp::Vector{Vector{Float64}},
                brs    ::Array{Int64,3},
                brsp   ::Array{Int64,3},
                σ²     ::Float64)

    llr::Float64 = 0.0

    @inbounds begin

      # trait likelihood
      for j = Base.OneTo(ntip)
        @simd for i = wf23[j]:(m-1)
          llr += llrdnorm_xμ(Xp[(i+1),j], Xc[(i+1),j],
                             Xp[i,j] + Eδx(LAp[i,j], ωx, δt[i]), 
                             Xc[i,j] + Eδx(LAc[i,j], ωx, δt[i]),
                             δt[i]*σ²)::Float64
        end
      end

      # biogeograhic likelihood
      for k = Base.OneTo(narea)
        llr += brll(stemevp[k], λ1, λ0, brsp[nedge,1,k]) -
               brll(stemevc[k], λ1, λ0,  brs[nedge,1,k])
        for j = Base.OneTo(ntip)
          llr += bitvectorll(Yp, λ1, λ0, ω1, ω0, LDp, δt, j, k, wf23[j], m) -
                 bitvectorll(Yc, λ1, λ0, ω1, ω0, LDc, δt, j, k, wf23[j], m)
        end
      end

    end

    return llr::Float64
  end


  return llf, llrf
end





"""
    llr_bm(Xc ::Array{Float64,2},
           Xp ::Array{Float64,2},
           δt ::Array{Float64,1},
           σ² ::Float64, 
           idx::UnitRange)

Estimate the proposal probability of a given path according to Brownian Motion.
"""
function llr_bm(Xc ::Array{Float64,2},
                Xp ::Array{Float64,2},
                idx::UnitRange,
                δt ::Array{Float64,1},
                σ² ::Float64)

  @inbounds begin
    llr::Float64 = 0.0

    @simd for i in Base.OneTo(length(idx)-1)
      llr += llrdnorm_xμ(Xc[idx[i+1]], Xp[idx[i+1]], 
                         Xc[idx[i]],   Xp[idx[i]],
                         (δt[i+1] - δt[i])*σ²)
    end
  end

  return llr::Float64
end






"""
    bitvectorll(y ::Array{Int64,1},
                λ1::Float64,
                λ0::Float64,
                ω1::Float64,
                ω0::Float64,
                Δx::Array{Float64,1},
                δt::Array{Float64,1})

Return the likelihood for a
bit vector (composed of 0s and 1s).
"""
function bitvectorll(Y ::Array{Int64,3},
                     λ1::Float64,
                     λ0::Float64,
                     ω1::Float64,
                     ω0::Float64,
                     LD::Array{Float64,3},
                     δt::Array{Float64,1},
                     j ::Int64,
                     k ::Int64,
                     yi::Int64,
                     m ::Int64)

  ll::Float64 = 0.0
  i ::Int64   = yi

  @inbounds begin

    if iszero(Y[yi,j,k])

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          ll += nell(δt[i], f_λ1(λ1, ω1, LD[i,j,k]))::Float64
          i += 1
        end
        i >= m && break

        ll += evll(δt[i], f_λ1(λ1, ω1, LD[i,j,k]))::Float64
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          ll += nell(δt[i], f_λ0(λ0, ω0, LD[i,j,k]))::Float64
          i += 1
        end
        i >= m && break

        ll += evll(δt[i], f_λ0(λ0, ω0, LD[i,j,k]))::Float64
        i += 1
      end

    else

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          ll += nell(δt[i], f_λ0(λ0, ω0, LD[i,j,k]))::Float64
          i += 1
        end
        i >= m && break

        ll += evll(δt[i], f_λ0(λ0, ω0, LD[i,j,k]))::Float64
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          ll += nell(δt[i], f_λ1(λ1, ω1, LD[i,j,k]))::Float64
          i += 1
        end
        i >= m && break

        ll += evll(δt[i], f_λ1(λ1, ω1, LD[i,j,k]))::Float64
        i += 1
      end
    end

  end

  return ll::Float64
end





"""
    bitvectorllr_ω1(Y  ::Array{Int64,3},
                    λ0 ::Float64,
                    ω0c::Float64,
                    ω0p::Float64,
                    LD ::Array{Float64,3},
                    δt ::Array{Float64,1},
                    j  ::Int64,
                    k  ::Int64,
                    yi ::Int64,
                    m  ::Int64)

Return the likelihood ratio for a
bit vector when updating `ω1`.
"""
function bitvectorllr_ω1(Y  ::Array{Int64,3},
                         λ1 ::Float64,
                         ω1c::Float64,
                         ω1p::Float64,
                         LD ::Array{Float64,3},
                         δt ::Array{Float64,1},
                         j  ::Int64,
                         k  ::Int64,
                         yi ::Int64,
                         m  ::Int64)

  llr::Float64 = 0.0
  i  ::Int64   = yi

  @inbounds begin

    if iszero(Y[yi,j,k])

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          llr += nellr_ω1(δt[i], ω1c, ω1p, λ1, LD[i,j,k])
          i += 1
        end
        i >= m && break

        llr += evllr_ω1(δt[i], ω1c, ω1p, λ1, LD[i,j,k])::Float64
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          i += 1
        end
        i >= m && break
        i += 1
      end

    else

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          i += 1
        end
        i >= m && break
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          llr += nellr_ω1(δt[i], ω1c, ω1p, λ1, LD[i,j,k])::Float64
          i += 1
        end
        i >= m && break

        llr += evllr_ω1(δt[i], ω1c, ω1p, λ1, LD[i,j,k])::Float64
        i += 1
      end
    end

  end

  return llr::Float64
end





"""
    bitvectorllr_ω0(Y  ::Array{Int64,3},
                    λ0 ::Float64,
                    ω0c::Float64,
                    ω0p::Float64,
                    LD ::Array{Float64,3},
                    δt ::Array{Float64,1},
                    j  ::Int64,
                    k  ::Int64,
                    yi ::Int64,
                    m  ::Int64)

Return the likelihood ratio for a
bit vector when updating `ω0`.
"""
function bitvectorllr_ω0(Y  ::Array{Int64,3},
                         λ0 ::Float64,
                         ω0c::Float64,
                         ω0p::Float64,
                         LD ::Array{Float64,3},
                         δt ::Array{Float64,1},
                         j  ::Int64,
                         k  ::Int64,
                         yi ::Int64,
                         m  ::Int64)

  llr::Float64 = 0.0
  i ::Int64   = yi

  @inbounds begin

    if iszero(Y[yi,j,k])

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          i += 1
        end

        i >= m && break
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          llr += nellr_ω0(δt[i], ω0c, ω0p, λ0, LD[i,j,k])
          i += 1
        end
        i >= m && break

        llr += evllr_ω0(δt[i], ω0c, ω0p, λ0, LD[i,j,k])::Float64
        i += 1
      end

    else

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          llr += nellr_ω0(δt[i], ω0c, ω0p, λ0, LD[i,j,k])
          i += 1
        end
        i >= m && break

        llr += evllr_ω0(δt[i], ω0c, ω0p, λ0, LD[i,j,k])::Float64
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          i += 1
        end

        i >= m && break
        i += 1
      end
    end

  end

  return llr::Float64
end






"""
    bitvectorllr_λ1(Y  ::Array{Int64,3},
                    λ0 ::Float64,
                    ω0c::Float64,
                    ω0p::Float64,
                    LD ::Array{Float64,3},
                    δt ::Array{Float64,1},
                    j  ::Int64,
                    k  ::Int64,
                    yi ::Int64,
                    m  ::Int64)

Return the likelihood ratio for a
bit vector when updating `ω1`.
"""
function bitvectorllr_λ1(Y  ::Array{Int64,3},
                         λ1c::Float64,
                         λ1p::Float64,
                         ω1 ::Float64,
                         LD ::Array{Float64,3},
                         δt ::Array{Float64,1},
                         j  ::Int64,
                         k  ::Int64,
                         yi ::Int64,
                         m  ::Int64)

  llr = 0.0
  i   = yi

  @inbounds begin

    if iszero(Y[yi,j,k])

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          llr += nellr_λ1(δt[i], ω1, λ1c, λ1p, LD[i,j,k])
          i   += 1
        end
        i >= m && break

        llr += evllr_λ1(δt[i], ω1, λ1c, λ1p, LD[i,j,k])::Float64
        i   += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          i += 1
        end
        i >= m && break
        i += 1
      end

    else

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          i += 1
        end
        i >= m && break
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          llr += nellr_λ1(δt[i], ω1, λ1c, λ1p, LD[i,j,k])
          i += 1
        end
        i >= m && break

        llr += evllr_λ1(δt[i], ω1, λ1c, λ1p, LD[i,j,k])::Float64
        i   += 1
      end
    end

  end

  return llr::Float64
end






"""
    bitvectorllr_λ0(Y  ::Array{Int64,3},
                    λ0c::Float64,
                    λ0p::Float64,
                    ω0 ::Float64,
                    LD ::Array{Float64,3},
                    δt ::Array{Float64,1},
                    j  ::Int64,
                    k  ::Int64,
                    yi ::Int64,
                    m  ::Int64)

Return the likelihood ratio for a
bit vector when updating `ω1`.
"""
function bitvectorllr_λ0(Y  ::Array{Int64,3},
                         λ0c::Float64,
                         λ0p::Float64,
                         ω0 ::Float64,
                         LD ::Array{Float64,3},
                         δt ::Array{Float64,1},
                         j  ::Int64,
                         k  ::Int64,
                         yi ::Int64,
                         m  ::Int64)

  llr::Float64 = 0.0
  i ::Int64   = yi

  @inbounds begin

    if iszero(Y[yi,j,k])

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          i += 1
        end
        i >= m && break
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          llr += nellr_λ0(δt[i], ω0, λ0c, λ0p, LD[i,j,k])
          i += 1
        end
        i >= m && break
        llr += evllr_λ0(δt[i], ω0, λ0c, λ0p, LD[i,j,k])::Float64
        i += 1
      end

    else

      while i < m
        while i < m && Y[i,j,k] == Y[i+1,j,k]
          llr += nellr_λ0(δt[i], ω0, λ0c, λ0p, LD[i,j,k])
          i += 1
        end
        i >= m && break
        llr += evllr_λ0(δt[i], ω0, λ0c, λ0p, LD[i,j,k])::Float64
        i += 1

        while i < m && Y[i,j,k] == Y[i+1,j,k]
          i += 1
        end
        i >= m && break
        i += 1
      end
    end

  end

  return llr::Float64
end





"""
    bitbitll(y1::Int64, y2::Int64, λ1::Float64, λ0::Float64, ω1::Float64, ω0::Float64, Δx::Float64, δt::Float64)

Return the likelihood for a
bit vector of length 2 (composed of 0s and 1s).
"""
function bitbitll(y1::Int64,
                  y2::Int64,
                  λ1::Float64,
                  λ0::Float64,
                  ω1::Float64,
                  ω0::Float64,
                  Δx::Float64,
                  δt::Float64)

  # event or non-event
  if y1 == y2 
    return nell(δt, y1 == 0 ? f_λ1(λ1, ω1, Δx) : f_λ0(λ0, ω0, Δx))::Float64
  else
    return evll(δt, y1 == 0 ? f_λ1(λ1, ω1, Δx) : f_λ0(λ0, ω0, Δx))::Float64
  end
end





"""
    makellr_biogeo(Y    ::Array{Int64,3},
                   δt   ::Vector{Float64},
                   narea::Int64,
                   ntip ::Int64,
                   m    ::Int64)

Make likelihood function for when updating ω1 & ω0.
"""
function makellr_biogeo(Y    ::Array{Int64,3},
                        δt   ::Vector{Float64},
                        narea::Int64,
                        ntip ::Int64,
                        m    ::Int64,
                        nedge::Int64)

  # which is 23 (23 = NaN) in each column
  const wf23 = Int64[]
  for j = Base.OneTo(ntip)
    push!(wf23, findfirst(Y[:,j,1] .!= 23))
  end

  function fω1(Y  ::Array{Int64,3}, 
               λ1 ::Float64,
               ω1c::Float64,
               ω1p::Float64,
               LD ::Array{Float64,3})

    llr::Float64 = 0.0

    @inbounds begin

      for k = Base.OneTo(narea), j = Base.OneTo(ntip)
        llr +=  bitvectorllr_ω1(Y, λ1, ω1c, ω1p, LD, δt, 
                                j, k, wf23[j], m)
      end
    end

    return llr
  end


  function fω0(Y  ::Array{Int64,3}, 
               λ0 ::Float64,
               ω0c::Float64,
               ω0p::Float64,
               LD ::Array{Float64,3})

    llr::Float64 = 0.0

    @inbounds begin

      for k = Base.OneTo(narea), j = Base.OneTo(ntip)
        llr += bitvectorllr_ω0(Y, λ0, ω0c, ω0p, LD, δt, 
                               j, k, wf23[j], m)
      end
    end

    return llr
  end


  function fλ1(Y      ::Array{Int64,3}, 
               λ1c    ::Float64,
               λ1p    ::Float64,
               λ0     ::Float64,
               ω1     ::Float64,
               LD     ::Array{Float64,3},
               stemevc::Vector{Vector{Float64}},
               brs    ::Array{Int64,3})

    llr::Float64 = 0.0

    @inbounds begin

      for k = Base.OneTo(narea)
        llr += brll(stemevc[k], λ1p, λ0, brs[nedge,1,k])::Float64 -
               brll(stemevc[k], λ1c, λ0, brs[nedge,1,k])::Float64
        for j = Base.OneTo(ntip)
          llr += bitvectorllr_λ1(Y, λ1c, λ1p, ω1, LD, δt, 
                                 j, k, wf23[j], m)
        end
      end
    end

    return llr
  end


  function fλ0(Y      ::Array{Int64,3}, 
               λ1     ::Float64,
               λ0c    ::Float64,
               λ0p    ::Float64,
               ω0     ::Float64,
               LD     ::Array{Float64,3},
               stemevc::Vector{Vector{Float64}},
               brs    ::Array{Int64,3})

    llr::Float64 = 0.0

    @inbounds begin

      for k = Base.OneTo(narea)
        llr += brll(stemevc[k], λ1, λ0p, brs[nedge,1,k])::Float64 -
               brll(stemevc[k], λ1, λ0c, brs[nedge,1,k])::Float64
        for j = Base.OneTo(ntip)
          llr += bitvectorllr_λ0(Y, λ0c, λ0p, ω0, LD, δt, 
                                 j, k, wf23[j], m)
        end
      end
    end

    return llr
  end


  return fω1, fω0, fλ1, fλ0
end





"""
    makellf_bgiid(bridx_a::Array{Array{UnitRange{Int64},1},1},
                  δt     ::Array{Float64,1},
                  narea  ::Int64,
                  nedge  ::Int64,
                  m      ::Int64)

Make triad likelihood function for the mutual 
independence model (iid), the proposal density 
for data augmented biogeographic histories.
"""
function makellf_bgiid(bridx_a::Array{Array{UnitRange{Int64},1},1},
                       δt     ::Array{Float64,1},
                       narea  ::Int64,
                       nedge  ::Int64,
                       m      ::Int64)

  # prepare δts
  const δtA = Array{Float64,1}[]

  for j=bridx_a[1][1:(nedge-1)]
    inds = zeros(Int64,length(j) - 1)
    for i = eachindex(inds)
      inds[i] = rowind(j[i], m)
    end
    push!(δtA, δt[inds])
  end

  function fiid(Y      ::Array{Int64,3},
                stemev::Array{Array{Float64,1},1},
                brs    ::Array{Int64,3},
                triad  ::Array{Int64,1},
                λϕ1    ::Float64,
                λϕ0    ::Float64)

    ll::Float64 = 0.0

    @inbounds begin

      pr, d1, d2 = triad::Array{Int64,1}

      if pr < nedge 
        for k = Base.OneTo(narea)
          ll += bitvectorll_iid(Y, bridx_a[k][pr], λϕ1, λϕ0, δtA[pr]) +
                bitvectorll_iid(Y, bridx_a[k][d1], λϕ1, λϕ0, δtA[d1]) +
                bitvectorll_iid(Y, bridx_a[k][d2], λϕ1, λϕ0, δtA[d2])
        end
      else 
        for k = Base.OneTo(narea)
          ll += bitvectorll_iid(Y, bridx_a[k][d1], λϕ1, λϕ0, δtA[d1]) +
                bitvectorll_iid(Y, bridx_a[k][d2], λϕ1, λϕ0, δtA[d2]) +
                brll(stemev[k], λϕ1, λϕ0, brs[nedge,1,k])
        end
      end

    end

    return ll::Float64
  end


  function fiidbr(Y      ::Array{Int64,3}, 
             stemevc::Array{Array{Float64,1},1},
             brs    ::Array{Int64,3},
             br     ::Int64,
             λϕ1    ::Float64,
             λϕ0    ::Float64)

    ll::Float64 = 0.0

    @inbounds begin
      if br == nedge
        for k = Base.OneTo(narea)
          ll += brll(stemevc[k], λϕ1, λϕ0, brs[nedge,1,k])
        end
      else
        for k = Base.OneTo(narea)
          ll += bitvectorll_iid(Y, bridx_a[k][br], λϕ1, λϕ0,  δtA[br])
        end
      end
    end
    
    return ll::Float64
  end


  return fiid, fiidbr
end






"""
    stem_llr(λ1     ::Float64,
             λ0     ::Float64,
             stemc  ::Array{Int64,1},
             stemp  ::Array{Int64,1},
             stemevc::Array{Array{Float64,1},1},
             stemevp::Array{Array{Float64,1},1},
             narea  ::Int64)

Estimate likelihood ratio for stem branch.
"""
function stem_llr(λ1     ::Float64,
                  λ0     ::Float64,
                  brs    ::Array{Int64,3},
                  brsp   ::Array{Int64,3},
                  stemevc::Array{Array{Float64,1},1},
                  stemevp::Array{Array{Float64,1},1},
                  narea  ::Int64,
                  nedge  ::Int64)

  @inbounds begin

    ll::Float64 = 0.0
    for k in Base.OneTo(narea)
      ll += brll(stemevp[k], λ1, λ0, brsp[nedge,1,k]) - 
            brll(stemevc[k], λ1, λ0, brs[ nedge,1,k])
    end
  end

  return ll::Float64
end





"""
    stemiid_propr(λϕ1    ::Float64,
                  λϕ0    ::Float64,
                  stemc  ::Array{Int64,1},
                  stemp  ::Array{Int64,1},
                  stemevc::Array{Array{Float64,1},1},
                  stemevp::Array{Array{Float64,1},1},
                  narea  ::Int64) 

Proposal ratio for stem.
"""
function stemiid_propr(λϕ1    ::Float64,
                       λϕ0    ::Float64,
                       brs    ::Array{Int64,3,},
                       brsp   ::Array{Int64,3},
                       stemevc::Array{Array{Float64,1},1},
                       stemevp::Array{Array{Float64,1},1},
                       narea  ::Int64,
                       nedge  ::Int64)

  @inbounds begin

    ll::Float64 = 0.0
    for k in Base.OneTo(narea)
      ll += brll(stemevc[k], λϕ1, λϕ0, brs[ nedge,1,k]) - 
            brll(stemevp[k], λϕ1, λϕ0, brsp[nedge,1,k])
    end
  end

  return ll::Float64
end




"""
    bitvectorll_iid(y::Array{Int64,1}, λ1::Float64, λ0::Float64, δt::Array{Float64,1})

Return likelihood under the independence model 
for a bit vector.
"""
function bitvectorll_iid(Y  ::Array{Int64,3},
                         idx::UnitRange{Int64},
                         λ1 ::Float64,
                         λ0 ::Float64,
                         δt ::Array{Float64,1})

  @inbounds begin

    ll::Float64  = 0.0

    cur_s::Int64   = Y[idx[1]]
    cur_λ::Float64 = cur_s == 0 ? λ1 : λ0

    for i = Base.OneTo(length(δt))
      if Y[idx[i]] == Y[idx[i+1]]
        ll += nell(δt[i], cur_λ)::Float64
      else
        ll += evll(δt[i], cur_λ)::Float64
        cur_s = 1 - cur_s
        cur_λ = cur_s == 0 ? λ1 : λ0
      end
    end

  end

  return ll::Float64
end






"""
    makellr_ωxupd(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64)

Make likelihood ratio function when updating `ωx`.
"""
function makellr_ωxσupd(δt  ::Vector{Float64}, 
                        Y   ::Array{Int64,3}, 
                        ntip::Int64)

  # which is 23 (i.e., NaN) in each column
  const w23 = UnitRange{Int64}[]
  for i = Base.OneTo(ntip)
    non23 = find(Y[:,i,1] .!= 23)
    push!(w23,colon(non23[1],non23[end-1]))
  end

  function fωx(X  ::Array{Float64,2},
               LAp::Array{Float64,2},
               LAn::Array{Float64,2},
               ωxc::Float64,
               ωxp::Float64,
               σ² ::Float64)

    llr::Float64 = 0.0

    @inbounds begin
      # trait likelihood
      if ωxp >= 0.0
        if ωxc >= 0.0
          # if ωxp >= 0.0 & ωxc >= 0.0
          for j = Base.OneTo(ntip)
            @simd for i = w23[j]
              llr += llrdnorm_ωx(X[(i+1),j], X[i,j],
                                 Eδx(LAp[i,j], ωxp, δt[i]),
                                 Eδx(LAp[i,j], ωxc, δt[i]), 
                                 δt[i]*σ²)
            end
          end
        else
          # if ωxp >= 0.0 & ωxc < 0.0
          for j = Base.OneTo(ntip)
            @simd for i = w23[j]
              llr += llrdnorm_ωx(X[(i+1),j], X[i,j],
                                 Eδx(LAp[i,j], ωxp, δt[i]),
                                 Eδx(LAn[i,j], ωxc, δt[i]), 
                                 δt[i]*σ²)
            end
          end
        end
      else
        if ωxc >= 0.0
          # if ωxp < 0.0 & ωxc >= 0.0
          for j = Base.OneTo(ntip)
            @simd for i = w23[j]
              llr += llrdnorm_ωx(X[(i+1),j], X[i,j],
                                 Eδx(LAn[i,j], ωxp, δt[i]),
                                 Eδx(LAp[i,j], ωxc, δt[i]), 
                                 δt[i]*σ²)
            end
          end
        else
          # if ωxp < 0.0 & ωxc < 0.0
          for j = Base.OneTo(ntip)
            @simd for i = w23[j]
              llr += llrdnorm_ωx(X[(i+1),j], X[i,j], 
                                 Eδx(LAn[i,j], ωxp, δt[i]),
                                 Eδx(LAn[i,j], ωxc, δt[i]), 
                                 δt[i]*σ²)
            end
          end
        end
      end
    end

    return llr::Float64
  end


  function fσ(X  ::Array{Float64,2},
              LA ::Array{Float64,2},
              ωx ::Float64,
              σ²c::Float64,
              σ²p::Float64)

    llr::Float64 = 0.0

    @inbounds begin
      for j = Base.OneTo(ntip)
        @simd for i = w23[j]
          llr += llrdnorm_σ²(X[(i+1),j],  
                             X[i,j] + Eδx(LA[i,j], ωx, δt[i]), 
                             δt[i]*σ²p, δt[i]*σ²c)
        end
      end
    end

    return llr::Float64
  end

  return fωx, fσ
end




"""
    makellr_Xupd(δt::Vector{Float64}, narea::Int64)

Make likelihood function for an internal node update in `X`.
"""
function makellr_XRupd(δt   ::Vector{Float64}, 
                       narea::Int64,
                       wcol ::Array{Array{Int64,1},1})

  δt1 = δt[1]
  wci = wcol[1]
 
  function fx(xi  ::Int64,
              xpi ::Array{Float64,1},
              X   ::Array{Float64,2},
              lapi::Array{Float64,1},
              ldpi::Array{Float64,2},
              LA  ::Array{Float64,2},
              LD  ::Array{Float64,3},
              Y   ::Array{Int64,3},
              ωx  ::Float64,
              ω1  ::Float64,
              ω0  ::Float64,
              λ1  ::Float64,
              λ0  ::Float64,
              σ²  ::Float64)

    # normal likelihoods
    llr::Float64 = 0.0

    @inbounds begin

      # loop for parent nodes
      δxim1 = δt[xi-1]
      for j = wcol[xi-1]
        llr += llrdnorm_x(xpi[j], X[xi,j], 
                          X[xi-1,j] + Eδx(LA[xi-1,j], ωx, δxim1), 
                          δxim1*σ²)
      end

      # loop for daughter nodes
      δxi = δt[xi]
      for j = wcol[xi]
        llr += llrdnorm_μ(X[xi+1, j],
                          xpi[j]  + Eδx(lapi[j],  ωx, δxi),
                          X[xi,j] + Eδx(LA[xi,j], ωx, δxi),
                          δxi*σ²)

        for k = Base.OneTo(narea)
          llr += bitbitll(Y[xi,j,k], Y[xi+1,j,k], 
                          λ1, λ0, ω1, ω0, ldpi[j,k], δxi)::Float64 -
                 bitbitll(Y[xi,j,k], Y[xi+1,j,k], 
                          λ1, λ0, ω1, ω0, LD[xi,j,k], δxi)::Float64
        end
      end
    end

    return llr::Float64
  end

  function fr(xpi ::Array{Float64,1},
             X   ::Array{Float64,2},
             σ²  ::Float64)

    llr::Float64 = 0.0

    @inbounds begin

      # loop for daughter nodes
      for j = wci
        llr += llrdnorm_μ(X[2,j], xpi[j], X[1,j], δt1*σ²)
      end
    end

    return llr::Float64
  end

  return fx, fr
end





"""
    brll(brevs::Array{Float64,1}, λ1::Float64, λ0::Float64, si::Int64)

Return likelihood for a branch in continuous time.
"""
function brll(brevs::Array{Float64,1}, λ1::Float64, λ0::Float64, si::Int64)

  ll ::Float64 = 0.0

  if endof(brevs) > 1 
    for i = Base.OneTo(endof(brevs)-1)
      ll += evll(brevs[i], iszero(si) ? λ1 : λ0)::Float64
      si = 1 - si
    end
  end

  ll += nell(brevs[end], iszero(si) ? λ1 : λ0)::Float64
  
  return ll::Float64
end





"""
    evll(t::Float64, λ::Float64)

Return log-likelihood for events.
"""
evll(t::Float64, λ::Float64) = (Base.Math.JuliaLibm.log(λ) - (λ * t))::Float64





"""
    evllr_ω1(t  ::Float64,
             ω1c::Float64,
             ω1p::Float64,
             λ1 ::Float64,
             δx ::Float64)

Return log-likelihood ratio for events when updating `ω1`.
"""
function evllr_ω1(t  ::Float64,
                  ω1c::Float64,
                  ω1p::Float64,
                  λ1 ::Float64,
                  δx ::Float64)
  @fastmath begin
    return Base.Math.JuliaLibm.log((1.0 + δx - ω1p)/(1.0 + δx - ω1c)) +
           λ1*t*(ω1p - ω1c)/(1.0 + δx)
  end
end





"""
    evllr_ω0(t  ::Float64,
             ω0c::Float64,
             ω0p::Float64,
             λ0 ::Float64,
             δx ::Float64)

Return log-likelihood ratio for events when updating `ω1`.
"""
function evllr_ω0(t  ::Float64,
                  ω0c::Float64,
                  ω0p::Float64,
                  λ0 ::Float64,
                  δx ::Float64)
  @fastmath begin
    return Base.Math.JuliaLibm.log((1.0 + δx + ω0p)/(1.0 + δx + ω0c)) +
           λ0*t*(ω0c - ω0p)/(1 + δx)
  end
end





"""
    evllr_λ1(t  ::Float64,
             ω1::Float64,
             λ1c::Float64,
             λ1p::Float64,
             δx ::Float64)

Return log-likelihood ratio for non-events when updating `λ1`.
"""
function evllr_λ1(t  ::Float64,
                  ω1::Float64,
                  λ1c::Float64,
                  λ1p::Float64,
                  δx ::Float64)
  @fastmath begin
    return Base.Math.JuliaLibm.log(λ1p/λ1c) +
           t*(1.0 - ω1/(1 + δx))*(λ1c - λ1p)
  end
end




"""
    evllr_λ0(t::Float64, λ0c::Float64, λ0p::Float64)

Return log-likelihood ratio for non-events when updating `λ0`.
"""
function evllr_λ0(t  ::Float64,
                  ω0::Float64,
                  λ0c::Float64,
                  λ0p::Float64,
                  δx ::Float64)
  @fastmath begin
    return Base.Math.JuliaLibm.log(λ0p/λ0c) +
           t*(1.0 + ω0/(1.0 + δx))*(λ0c - λ0p)
  end
end





"""
    nell(t::Float64, λ::Float64)

Return log-likelihood for nonevents.
"""
nell(t::Float64, λ::Float64) = -λ*t






"""
    nellr_ω1(t  ::Float64,
             ω1c::Float64,
             ω1p::Float64,
             λ1 ::Float64,
             δx ::Float64)

Return log-likelihood ratio for non-events when updating `ω1`.
"""
function nellr_ω1(t  ::Float64,
                  ω1c::Float64,
                  ω1p::Float64,
                  λ1 ::Float64,
                  δx ::Float64)
  @fastmath begin
    return λ1*t*(ω1p - ω1c)/(1.0 + δx)
  end
end







"""
    nellr_ω0(t  ::Float64,
             ω0c::Float64,
             ω0p::Float64,
             λ0 ::Float64,
             δx ::Float64)

Return log-likelihood ratio for non-events when updating `ω0`.
"""
function nellr_ω0(t  ::Float64,
                  ω0c::Float64,
                  ω0p::Float64,
                  λ0 ::Float64,
                  δx ::Float64)
  @fastmath begin
    return λ0*t*(ω0c - ω0p)/(1 + δx)
  end
end






"""
    nellr_λ1(t  ::Float64,
             ω1::Float64,
             λ1c::Float64,
             λ1p::Float64,
             δx ::Float64)

Return log-likelihood ratio for non-events when updating `λ1`.
"""
function nellr_λ1(t  ::Float64,
                  ω1::Float64,
                  λ1c::Float64,
                  λ1p::Float64,
                  δx ::Float64)
  @fastmath begin
    return t*(1.0 - ω1/(1 + δx))*(λ1c - λ1p)
  end
end





"""
    nellr_λ0(t  ::Float64, 
             ω0 ::Float64, 
             λ0c::Float64, 
             λ0p::Float64, 
             δx ::Float64)

Return log-likelihood ratio for non-events when updating `λ0`.
"""
function nellr_λ0(t  ::Float64,
                  ω0 ::Float64,
                  λ0c::Float64,
                  λ0p::Float64,
                  δx ::Float64)
  @fastmath begin
    return t*(1.0 + ω0/(1.0 + δx))*(λ0c - λ0p)
  end
end






"""
    allλpr(λc::Array{Float64,2}, λprior::Float64)

Return log-prior for all areas 
"""
function allλpr(λ1    ::Float64,
                λ0    ::Float64,
                λprior::Float64)
  @fastmath 2.0*Base.Math.JuliaLibm.log(λprior) - λprior * (λ1 + λ0)
end





"""
    logdexp(x::Float64, λ::Float64)

Compute the logarithmic transformation of the 
**Exponential** density with mean `λ` for `x`.
"""
logdexp(x::Float64, λ::Float64) = 
  @fastmath Base.Math.JuliaLibm.log(λ) - λ * x




"""
    llrdexp_x(xp::Float64, xc::Float64, λ::Float64)

Compute the loglik ratio of the 
**Exponential** density for proposal 
`xp` given current `xc` both with mean `λ`.
"""
llrdexp_x(xp::Float64, xc::Float64, λ::Float64) = 
  @fastmath λ * (xc - xp)





"""
    logdbeta(x::Float64, α::Float64, β::Float64)

Compute the logarithmic transformation of the 
**Beta** density with shape `α` and shape `β` for `x`.
"""
logdbeta(x::Float64, α::Float64, β::Float64) = 
  @fastmath ((α-1.0) * Base.Math.JuliaLibm.log(x)                 +
             (β-1.0) * Base.Math.JuliaLibm.log(1.0 - x)           +
             Base.Math.JuliaLibm.log(gamma(α + β)/(gamma(α)*gamma(β))))





"""
    llrdbeta_x(xp::Float64, xc::Float64, α::Float64, β::Float64)

Compute the logarithmic ratio for the **Beta** density 
with shape `α` and shape `β` between `xp` and `xc`.
"""
llrdbeta_x(xp::Float64, xc::Float64, α::Float64, β::Float64) = 
  @fastmath ((α-1.0) * Base.Math.JuliaLibm.log(xp/xc) +
             (β-1.0) * Base.Math.JuliaLibm.log((1.0 - xp)/(1.0 - xc)))





"""
    logdnorm(x::Float64, μ::Float64, σ²::Float64)
  
Compute the logarithmic transformation of the 
**Normal** density with mean `μ` and variance `σ²` for `x`.
"""
logdnorm(x::Float64, μ::Float64, σ²::Float64) = 
  @fastmath -(0.5*Base.Math.JuliaLibm.log(2.0π) +
              0.5*Base.Math.JuliaLibm.log(σ²)   +
              (x - μ)^2/(2.0σ²))





"""
    logdnorm_tc(x::Float64, μ::Float64, σ²::Float64)

Compute the logarithmic transformation of the 
**Normal** density with mean `μ` and variance `σ²` for `x`, up to a constant
"""
logdnorm_tc(x::Float64, μ::Float64, σ²::Float64) =
  @fastmath -0.5*Base.Math.JuliaLibm.log(σ²) - 
            (x - μ)^2/(2.0σ²)::Float64





"""
    llrdnorm_ωx(x::Float64, xi::Float64, μp::Float64, μc::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `ωx` updates
"""
llrdnorm_ωx(x::Float64, xi::Float64, μp::Float64, μc::Float64, σ²::Float64) =
  @fastmath (-(x - xi - μp)^2 + (x - xi - μc)^2)/(2.0σ²)





"""
    llrdnorm_σ²(x::Float64, μ::Float64, σ²p::Float64, σ²c::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `σ²` updates
"""
llrdnorm_σ²(x::Float64, μ::Float64, σ²p::Float64, σ²c::Float64) = 
  @fastmath -0.5*(Base.Math.JuliaLibm.log(σ²p/σ²c) +
                  (x - μ)^2*(1.0/σ²p - 1.0/σ²c))





"""
    llrdnorm_μ(x::Float64, μp::Float64, μc::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `μ` updates
"""
llrdnorm_μ(x::Float64, μp::Float64, μc::Float64, σ²::Float64) =
  @fastmath ((x - μc)^2 - (x - μp)^2)/(2.0σ²)





"""
    llrdnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `x` updates
"""
llrdnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64) =
  @fastmath ((xc - μ)^2 - (xp - μ)^2)/(2.0σ²)




"""
    llrdnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `x` updates
"""
llrdnorm_xμ(xp::Float64, xc::Float64, μp::Float64, μc::Float64, σ²::Float64) =
  @fastmath ((xc - μc)^2 - (xp - μp)^2)/(2.0σ²)





"""
    logdhcau(x::Float64, scl::Float64)

Compute the logarithmic transformation of the 
**Half-Cauchy** density with scale `scl` for `x`.
"""
logdhcau(x::Float64, scl::Float64) = 
  @fastmath Base.Math.JuliaLibm.log(2.0 * scl/(π *(x * x + scl * scl)))





"""
    logdhcau1(x::Float64)
  
Compute the logarithmic transformation of the 
**Half-Cauchy** density with scale of 1 for `x`.
"""
logdhcau1(x::Float64) = 
  @fastmath Base.Math.JuliaLibm.log(2.0/(π * (x * x + 1.)))



