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
Eδx(μ::Float64, ωx::Float64, δt::Float64) = ωx * μ * δt





"""
    f_λ1(λ1::Float64, ω1::Float64, δx::Float64)

Estimate rates for area colonization based 
on the difference between lineage traits and area averages.
"""
function f_λ1(λ1::Float64, ω1::Float64, δx::Float64)
  if iszero(δx) 
    return λ1
  elseif ω1 < 0.0
    return λ1 * exp(-abs2(ω1)/δx)
  else
    return λ1 * exp(-abs2(ω1)*δx)
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
  elseif ω0 < 0.0
    return λ0 + Base.Math.JuliaLibm.log(1. + abs2(ω0)/δx)
  else
    return λ0 + Base.Math.JuliaLibm.log(1. + abs2(ω0)*δx)
  end
end





"""
    makellf(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64, narea::Int64)

Make likelihood function for all trait matrix 
and biogeography history.
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
  for i in wf23
    n += length(i:(m-1))
  end

  # normal constant
  const normC = -0.5*log(2.0π)*n

  function f(X      ::Array{Float64,2},
             Y      ::Array{Int64,3}, 
             LA     ::Array{Float64,2},
             LD     ::Array{Float64,3},
             ωx     ::Float64,
             ω1     ::Float64,
             ω0     ::Float64,
             λ1     ::Float64,
             λ0     ::Float64,
             stemevc::Vector{Vector{Float64}},
             brs    ::Array{Int64,3},
             σ²     ::Float64)

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
        ll += brll(stemevc[k], λ1, λ0, brs[nedge,1,k])::Float64
        for j = Base.OneTo(ntip)
          ll += bitvectorll(Y, λ1, λ0, ω1, ω0, LD, δt, 
                            j, k, wf23[j], m)::Float64
        end
      end

    end

    return ll::Float64
  end

  return f
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
      llr += logdnorm_tc(Xc[idx[i+1]], Xc[idx[i]], (δt[i+1] - δt[i])*σ²) -
             logdnorm_tc(Xp[idx[i+1]], Xp[idx[i]], (δt[i+1] - δt[i])*σ²)
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
    makellr_λ_upd(Y::Array{Int64,3}, δt::Vector{Float64}, narea::Int64)

Make likelihood ratio function for when updating λ.
"""
function makellr_λ_upd(Y    ::Array{Int64,3},
                       δt   ::Vector{Float64},
                       narea::Int64,
                       ntip ::Int64,
                       m    ::Int64, 
                       nedge::Int64)

  # get initial range
  const wf23 = Int64[]
  for j = Base.OneTo(ntip)
    push!(wf23, findfirst(Y[:,j,1] .!= 23))
  end

  function f(Y      ::Array{Int64,3}, 
             λ1c    ::Float64,
             λ0c    ::Float64,
             λ1p    ::Float64,
             λ0p    ::Float64,
             ω1     ::Float64,
             ω0     ::Float64,
             LD     ::Array{Float64,3},
             stemevc::Vector{Vector{Float64}},
             brs    ::Array{Int64,3})

    ll::Float64 = 0.0

    @inbounds begin

      for k = Base.OneTo(narea)
        ll += brll(stemevc[k], λ1p, λ0p, brs[nedge,1,k])::Float64 -
              brll(stemevc[k], λ1c, λ0c, brs[nedge,1,k])::Float64

        for j = Base.OneTo(ntip)
          ll += bitvectorll(Y, λ1p, λ0p, ω1, ω0, LD, δt, 
                            j, k, wf23[j], m)::Float64 -
                bitvectorll(Y, λ1c, λ0c, ω1, ω0, LD, δt, 
                            j, k, wf23[j], m)::Float64
        end
      end
    
    end

    return ll::Float64
  end

  return f
end





"""
    makellr_ω1μ_upd(Y::Array{Int64,3}, δt::Vector{Float64}, narea::Int64)

Make likelihood function for when updating ω1 & ω0.
"""
function makellr_ω10_upd(Y    ::Array{Int64,3},
                         δt   ::Vector{Float64},
                         narea::Int64,
                         ntip ::Int64,
                         m    ::Int64)

  # which is 23 (23 = NaN) in each column
  const wf23 = Int64[]
  for j = Base.OneTo(ntip)
    push!(wf23, findfirst(Y[:,j,1] .!= 23))
  end

  function f(Y  ::Array{Int64,3}, 
             λ1 ::Float64,
             λ0 ::Float64,
             ω1c::Float64,
             ω0c::Float64,
             ω1p::Float64,
             ω0p::Float64,
             LD ::Array{Float64,3})

    ll::Float64 = 0.0

    @inbounds begin

      for k = Base.OneTo(narea), j = Base.OneTo(ntip)
        ll +=  bitvectorll(Y, λ1, λ0, ω1p, ω0p, LD, δt, 
                           j, k, wf23[j], m)::Float64 -
               bitvectorll(Y, λ1, λ0, ω1c, ω0c, LD, δt, 
                           j, k, wf23[j], m)::Float64
      end
    end

    return ll
  end

  return f
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

  function f(Y      ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1},
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
                brll(stemevc[k], λϕ1, λϕ0, brs[nedge,1,k])
        end
      end

    end

    return ll::Float64
  end

  return f
end





"""
    makellf_biogeo_upd_iid_br(bridx_a::Array{Array{Array{Int64,1},1},1}, δt::Array{Float64,1}, narea::Int64, nedge::Int64, m::Int64)

Make single branch likelihood function for the mutual 
independence model (iid), the proposal density 
for data augmented biogeographic histories.
"""
function makellf_bgiid_br(bridx_a::Array{Array{UnitRange{Int64},1},1},
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

  function f(Y      ::Array{Int64,3}, 
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

  return f
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
function bitvectorll_iid(Y     ::Array{Int64,3},
                         idx   ::UnitRange{Int64},
                         λ1    ::Float64,
                         λ0    ::Float64,
                         δt    ::Array{Float64,1})

  @inbounds begin

    ll::Float64  = 0.0

    cur_s::Int64   = Y[idx[1]]
    cur_λ::Float64 = cur_s == 0 ? λ1 : λ0

    for i = Base.OneTo(endof(δt))
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
    makellr_σ²ωxupd(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64)

Make likelihood function for all trait matrix, `X`.
"""
function makellr_σ²ωxupd(δt  ::Vector{Float64}, 
                         Y   ::Array{Int64,3}, 
                         ntip::Int64)

  # which is 23 (i.e., NaN) in each column
  const w23 = UnitRange{Int64}[]
  for i = Base.OneTo(ntip)
    non23 = find(Y[:,i,1] .!= 23)
    push!(w23,colon(non23[1],non23[end-1]))
  end

  function f(X  ::Array{Float64,2},
             LA ::Array{Float64,2},
             ωxc::Float64,
             ωxp::Float64,
             σ²c::Float64,
             σ²p::Float64)

    llr::Float64 = 0.0

    @inbounds begin

      # trait likelihood
      for j = Base.OneTo(ntip)
        @simd for i = w23[j]
          llr += logdnorm_tc(X[(i+1),j], 
                             X[i,j] + Eδx(LA[i,j], ωxp, δt[i]), 
                             δt[i]*σ²p)::Float64 -
                 logdnorm_tc(X[(i+1),j], 
                             X[i,j] + Eδx(LA[i,j], ωxc, δt[i]), 
                             δt[i]*σ²c)::Float64
        end
      end

    end

    return llr::Float64
  end

  return f
end






"""
    makellr_Xupd(δt::Vector{Float64}, narea::Int64)

Make likelihood function for an internal node update in `X`.
"""
function makellr_Xupd(δt   ::Vector{Float64}, 
                      narea::Int64,
                      wcol ::Array{Array{Int64,1},1})

  function f(xi    ::Int64,
             xpi   ::Array{Float64,1},
             X     ::Array{Float64,2},
             lapi  ::Array{Float64,1},
             ldpi  ::Array{Float64,2},
             LA    ::Array{Float64,2},
             LD    ::Array{Float64,3},
             Y     ::Array{Int64,3},
             ωx    ::Float64,
             ω1    ::Float64,
             ω0    ::Float64,
             λ1    ::Float64,
             λ0    ::Float64,
             σ²    ::Float64)

    # normal likelihoods
    llr::Float64 = 0.0

    @inbounds begin

      # loop for parent nodes
      δxim1 = δt[xi-1]
      for j = wcol[xi-1]
        llr += logdnorm_tc(xpi[j], 
                           X[xi-1,j] + Eδx(LA[xi-1,j], ωx, δxim1), 
                           δxim1*σ²)::Float64 -
               logdnorm_tc(X[xi,j],
                           X[xi-1,j] + Eδx(LA[xi-1,j], ωx, δxim1), 
                           δxim1*σ²)::Float64
      end

      # loop for daughter nodes
      δxi = δt[xi]
      for j = wcol[xi]
        llr += logdnorm_tc(X[xi+1, j], 
                           xpi[j]  + Eδx(lapi[j], ωx, δxi), 
                           δxi*σ²)::Float64 -
               logdnorm_tc(X[xi+1, j], 
                           X[xi,j] + Eδx(LA[xi,j], ωx, δxi), 
                           δxi*σ²)::Float64

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

  return f
end





"""
    makellf_Rupd(δt::Vector{Float64}, narea::Int64)

Make likelihood function
for the root update in `X`.
"""
function makellr_Rupd(δt1  ::Float64, 
                      wci  ::Array{Int64,1})

  function f(xpi ::Array{Float64,1},
             X   ::Array{Float64,2},
             σ²  ::Float64)

    llr::Float64 = 0.0

    @inbounds begin

      # loop for daughter nodes
      for j = wci
        # trait likelihood
        llr += logdnorm_tc(X[2,j], xpi[j], δt1*σ²)::Float64 -
               logdnorm_tc(X[2,j], X[1,j], δt1*σ²)::Float64
      end
    end

    return llr::Float64
  end

  return f
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
    nell(t::Float64, λ::Float64)

Return log-likelihood for nonevents.
"""
nell(t::Float64, λ::Float64) = (-1 * λ * t)::Float64





"""
    allλpr(λc::Array{Float64,2}, λprior::Float64)

Return log-prior for all areas 
"""
function allλpr(λ1    ::Float64,
                λ0    ::Float64,
                λprior::Float64)
  return (logdexp(λ1, λprior) + logdexp(λ0, λprior))::Float64
end

