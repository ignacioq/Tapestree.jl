#=

Likelihood Functions for joint
Biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

May 15 2017

=#




"""
    E_sde(xi::Float64, μ::Float64, ωx::Float64, δt::Float64)

Return the expected value according to reworked competition model.
"""
function E_sde(xi::Float64, μ::Float64, ωx::Float64, δt::Float64)
  @fastmath begin
    if ωx < 0.0
      if μ - xi < 0.0
        return (-ωx * exp(μ - xi) * δt)::Float64
      elseif μ - xi > 0.0
        return ( ωx * exp(xi - μ) * δt)::Float64
      else
        return 0.0
      end
    else
       return (ωx * (μ - xi) * δt)::Float64
    end
  end
end





"""
    E_sde_old(xi::Float64, μ::Float64, ωx::Float64, δt::Float64)

Return the expected value according to reworked competition model.
"""
function E_sde_old(xi::Float64, μ::Float64, ωx::Float64, δt::Float64)
  return (ωx * (μ - xi) * δt)::Float64
end





"""
    makellf(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64, narea::Int64)

Make likelihood function for all trait matrix 
and biogeography history.
"""
function makellf(δt   ::Vector{Float64}, 
                 Y    ::Array{Int64, 3}, 
                 ntip ::Int64, 
                 narea::Int64)

  # which is 23 (23 = NaN) in each column
  const w23 = UnitRange{Int64}[]
  for i = Base.OneTo(ntip)
    non23 = find(Y[:,i,1] .!= 23)
    push!(w23,colon(non23[1],non23[end-1]))
  end

  # number of normal evaluations
  n = 0
  for i in w23
    n += length(i)
  end

  # normal constant
  const normC = -0.5*log(2.0π)*n

  function f(X      ::Array{Float64,2},
             Y      ::Array{Int64,3}, 
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             ωx     ::Float64,
             ω1     ::Float64,
             ω0     ::Float64,
             λ      ::Array{Float64,1},
             stemevc::Vector{Vector{Float64}},
             stemss ::Vector{Int64},
             σ²     ::Float64)

    ll::Float64 = normC

    @inbounds @fastmath begin

      # trait likelihood
      for j = Base.OneTo(ntip)
        @simd for i = w23[j]
          ll += logdnorm_tc(X[(i+1),j], 
                            X[i,j] + 
                            E_sde(X[i,j], linavg[i,j], ωx, δt[i]), 
                            δt[i]*σ²)::Float64
        end
      end

      # biogeograhic likelihood
      for k = Base.OneTo(narea)
        ll += brll(stemevc[k], λ[1], λ[2], stemss[k])::Float64
        for j = Base.OneTo(ntip)
          ll += bitvectorll(Y[w23[j][1]:(w23[j][end]+1),j,k], 
                            λ[1], λ[2], ω1, ω0, 
                            lindiff[w23[j],j,k], δt[w23[j]])::Float64
        end
      end

    end

    ll::Float64
  end

  return f
end







"""
    bitvectorll(y ::Array{Int64,1}, λ1::Float64, λ0::Float64, ω1::Float64, ω0::Float64, Δx::Array{Float64,1}, δt::Array{Float64,1})

Return the likelihood for a
bit vector (composed of 0s and 1s).
"""
function bitvectorll(y ::Array{Int64,1},
                     λ1::Float64,
                     λ0::Float64,
                     ω1::Float64,
                     ω0::Float64,
                     Δx::Array{Float64,1},
                     δt::Array{Float64,1})

  ll::Float64 = 0.0

  @inbounds begin

    cur_s::Int64 = y[1]

    if cur_s == 0 
      cur_λ, cur_ω = λ1, ω1
    else
      cur_λ, cur_ω = λ0, ω0
    end

    for i = Base.OneTo(endof(y)-1)

      if y[i] == y[i+1]
        ll += nell(δt[i], f_λ(cur_λ, cur_ω, Δx[i]))::Float64
      else
        ll += evll(δt[i], f_λ(cur_λ, cur_ω, Δx[i]))::Float64
        cur_s = 1 - cur_s

        if cur_s == 0 
          cur_λ, cur_ω = λ1, ω1
        else
          cur_λ, cur_ω = λ0, ω0
        end
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
  if y1::Int64 == y2::Int64 
    return nell(δt, y1::Int64 == 0 ? f_λ(λ1, ω1, Δx) : f_λ(λ0, ω0, Δx))::Float64
  else
    return evll(δt, y1::Int64 == 0 ? f_λ(λ1, ω1, Δx) : f_λ(λ0, ω0, Δx))::Float64
  end
end





"""
    f_λ(λ::Float64, ω::Float64, Δx::Float64)

Estimate rates for area colonization/loss based 
on the difference between lineage traits and area averages.
"""
f_λ(λ::Float64, ω::Float64, Δx::Float64) = @fastmath (λ * exp(ω*Δx))::Float64





"""
    makellf_λ_upd(Y::Array{Int64,3}, δt::Vector{Float64}, narea::Int64)

Make likelihood function for when updating λ.
"""
function makellf_λ_upd(Y    ::Array{Int64,3},
                       δt   ::Vector{Float64},
                       narea::Int64)

  const coloop = Base.OneTo(size(Y,2))

  # which is 23 (23 = NaN) in each column
  const w23 = UnitRange{Int64}[]
  for i=Base.OneTo(size(Y,2))
    non23 = find(Y[:,i,1] .!= 23)
    push!(w23,colon(non23[1],non23[end-1]))
  end


  function f(Y      ::Array{Int64,3}, 
             λ::Array{Float64,1},
             ω1     ::Float64,
             ω0     ::Float64,
             lindiff::Array{Float64,3},
             stemevc::Vector{Vector{Float64}},
             stemss ::Vector{Int64})

    ll::Float64 = 0.0

    @inbounds begin

      for k = Base.OneTo(narea)
        ll += brll(stemevc[k], λ[1], λ[2], stemss[k])::Float64
        for j = coloop
          ll += bitvectorll(Y[w23[j][1]:(w23[j][end]+1),j,k], 
                            λ[1], λ[2], ω1, ω0, 
                            lindiff[w23[j],j,k], δt[w23[j]])::Float64
        end
      end
    
    end

    return ll::Float64
  end

  return f
end





"""
    makellf_ω1μ_upd(Y::Array{Int64,3}, δt::Vector{Float64}, narea::Int64)

Make likelihood function for when updating ω1 & ω0.
"""
function makellf_ω10_upd(Y   ::Array{Int64,3},
                         δt   ::Vector{Float64},
                         narea::Int64)

  const coloop = Base.OneTo(size(Y,2))

  # which is 23 (23 = NaN) in each column
  const w23 = UnitRange{Int64}[]
  for i = Base.OneTo(size(Y,2))
    non23 = find(Y[:,i,1] .!= 23)
    push!(w23,colon(non23[1],non23[end-1]))
  end

  function f(Y      ::Array{Int64,3}, 
             λ      ::Array{Float64,1},
             ω1     ::Float64,
             ω0     ::Float64,
             lindiff::Array{Float64,3})

    ll::Float64 = 0.0

    @inbounds begin

      for k = Base.OneTo(narea), j = coloop
          ll += bitvectorll(Y[w23[j][1]:(w23[j][end]+1),j,k], 
                            λ[1], λ[2], ω1, ω0, 
                            lindiff[w23[j],j,k], δt[w23[j]])::Float64
      end

    end

    return ll
  end

  return f
end





"""
    makellf_biogeo_upd_iid(bridx_a::Array{Array{Array{Int64,1},1},1}, δt::Array{Float64,1}, narea::Int64, nedge::Int64, m::Int64)

Make triad likelihood function for the mutual 
independence model (iid), the proposal density 
for data augmented biogeographic histories.
"""
function makellf_biogeo_upd_iid(bridx_a::Array{Array{UnitRange{Int64},1},1},
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


  function f(Y     ::Array{Int64,3}, 
             λ     ::Array{Float64,1},
             ω1    ::Float64,
             ω0    ::Float64,
             avg_Δx::Array{Float64,2},
             triad ::Array{Int64,1})

    ll::Float64 = 0.0

    @inbounds begin

      pr, d1, d2 = triad::Array{Int64,1}

      if pr < nedge 
        for j=Base.OneTo(narea)
          ll += bitvectorll_iid(Y[bridx_a[j][pr]], 
                                λ[1], ω1, λ[2], ω0, avg_Δx[pr,j], δtA[pr]) +
                bitvectorll_iid(Y[bridx_a[j][d1]], 
                                λ[1], ω1, λ[2], ω0, avg_Δx[d1,j], δtA[d1]) +
                bitvectorll_iid(Y[bridx_a[j][d2]], 
                                λ[1], ω1, λ[2], ω0, avg_Δx[d2,j], δtA[d2])::Float64
        end
      else 
        for j=Base.OneTo(narea)
          ll += bitvectorll_iid(Y[bridx_a[j][d1]], 
                                λ[1], ω1, λ[2], ω0, avg_Δx[d1,j], δtA[d1]) +
                bitvectorll_iid(Y[bridx_a[j][d2]], 
                                λ[1], ω1, λ[2], ω0, avg_Δx[d2,j], δtA[d2])::Float64
        end
      end

    end

    return ll::Float64
  end

  return f
end




"""
    bitvectorll_iid(y::Array{Int64,1}, λ1::Float64, λ0::Float64, δt::Array{Float64,1})

Return likelihood under the independence model 
for a bit vector.
"""
function bitvectorll_iid(y     ::Array{Int64,1},
                         λ1    ::Float64,
                         ω1    ::Float64,
                         λ0    ::Float64,
                         ω0    ::Float64,
                         avg_Δx::Float64,
                         δt    ::Array{Float64,1})

  @inbounds begin

    ll::Float64  = 0.0
    λt1::Float64 = f_λ(λ1, ω1, avg_Δx)
    λt0::Float64 = f_λ(λ0, ω0, avg_Δx)

    cur_s::Int64   = y[1]
    cur_λ::Float64 = cur_s == 0 ? λt1 : λt0

    for i=Base.OneTo(endof(y)-1)
      if y[i] == y[i+1]
        ll += nell(δt[i], cur_λ)::Float64
      else
        ll += evll(δt[i], cur_λ)::Float64
        cur_s = 1 - cur_s
        cur_λ = cur_s == 0 ? λt1 : λt0
      end
    end

  end

  ll::Float64
end







"""
    makellf_σ²ωxupd(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64)

Make likelihood function for all trait matrix, `X`.
"""
function makellf_σ²ωxupd(δt  ::Vector{Float64}, 
                         Y   ::Array{Int64,3}, 
                         ntip::Int64)

  # which is 23 (i.e., NaN) in each column
  const w23 = UnitRange{Int64}[]
  for i = Base.OneTo(ntip)
    non23 = find(Y[:,i,1] .!= 23)
    push!(w23,colon(non23[1],non23[end-1]))
  end

  # number of normal evaluations
  n = 0
  for i = w23
    n += length(i)
  end

  # normal constant
  const normC = -0.5*log(2.0π)*n

  function f(X ::Array{Float64,2},
             la::Array{Float64,2},
             ωx::Float64,
             σ²::Float64)

    ll::Float64 = normC

    @inbounds @fastmath begin

      # trait likelihood
      for j = Base.OneTo(ntip)
        @simd for i = w23[j]
          ll += logdnorm_tc(X[(i+1),j], 
                          X[i,j] + 
                          E_sde(X[i,j], la[i,j], ωx, δt[i]), 
                          δt[i]*σ²)
        end
      end
    
    end

    return ll::Float64
  end

  return f
end




"""
    makellf_Xupd(δt::Vector{Float64}, narea::Int64)

Make likelihood function for an internal node update in `X`.
"""
function makellr_Xupd(δt   ::Vector{Float64}, 
                      narea::Int64)

  function f(i     ::Int64,
             wci   ::Array{Int64,1},
             wcim1 ::Array{Int64,1},
             xpi   ::Array{Float64,1},
             xci   ::Array{Float64,1},
             xcm1  ::Array{Float64,1},
             xcp1  ::Array{Float64,1},
             lapi  ::Array{Float64,1},
             ldpi  ::Array{Float64,2},
             laci  ::Array{Float64,1},
             lacim1::Array{Float64,1},
             ldci  ::Array{Float64,2},
             Y     ::Array{Int64,3},
             ωx    ::Float64,
             ω1    ::Float64,
             ω0    ::Float64,
             λ     ::Array{Float64,1},
             σ²    ::Float64)

    # normal likelihoods
    ll::Float64 = 0.0

    @inbounds begin

      # loop for parent nodes
      for j = wcim1
        ll += logdnorm_tc(xpi[j], 
                          xcm1[j] + 
                          E_sde(xcm1[j], lacim1[j], ωx, δt[i-1]), 
                          δt[i-1]*σ²)::Float64 -
              logdnorm_tc(xci[j],
                          xcm1[j] + 
                          E_sde(xcm1[j], lacim1[j], ωx, δt[i-1]), 
                          δt[i-1]*σ²)::Float64
      end

      # loop for daughter nodes
      for j = wci
        # trait likelihood
        ll += logdnorm_tc(xcp1[j], 
                          xpi[j] + 
                          E_sde(xpi[j], lapi[j], ωx, δt[i]), 
                          δt[i]*σ²)::Float64 -
              logdnorm_tc(xcp1[j], 
                          xci[j] + 
                          E_sde(xci[j], laci[j], ωx, δt[i]), 
                          δt[i]*σ²)::Float64
      end

        # biogeograhic likelihoods
      for k = Base.OneTo(narea), j = wci
        ll += bitbitll(Y[i,j,k], Y[i+1,j,k], 
                       λ[1], λ[2], ω1, ω0, ldpi[j,k], δt[i])::Float64 -
              bitbitll(Y[i,j,k], Y[i+1,j,k], 
                       λ[1], λ[2], ω1, ω0, ldci[j,k], δt[i])::Float64
      end
      
    end

    return ll::Float64
  end

  return f
end






"""
    makellf_Rupd(δt::Vector{Float64}, narea::Int64)

Make likelihood function
for the root update in `X`.
"""
function makellr_Rupd(δt1  ::Float64, 
                      narea::Int64)

  function f(wci   ::Array{Int64,1},
             xpi   ::Array{Float64,1},
             xci   ::Array{Float64,1},
             xcp1  ::Array{Float64,1},
             lapi  ::Array{Float64,1},
             ldpi  ::Array{Float64,2},
             laci  ::Array{Float64,1},
             ldci  ::Array{Float64,2},
             Y     ::Array{Int64,3},
             ωx    ::Float64,
             ω1    ::Float64,
             ω0    ::Float64,
             λ     ::Array{Float64,1},
             σ²    ::Float64)

    ll::Float64 = 0.0

    @inbounds @fastmath begin

      # loop for daughter nodes
      for j = eachindex(wci)
        # trait likelihood
        ll += logdnorm_tc(xcp1[j], 
                          xpi[j] + 
                          E_sde(xpi[j], lapi[j], ωx, δt1), 
                          δt1*σ²)::Float64 -
              logdnorm_tc(xcp1[j], 
                          xci[j] + 
                          E_sde(xci[j], laci[j], ωx, δt1), 
                          δt1*σ²)::Float64

        # biogeograhic likelihoods
        for k = Base.OneTo(narea)
          ll += bitbitll(Y[1,wci[j],k], Y[2,wci[j],k], 
                         λ[1], λ[2], ω1, ω0, ldpi[j,k], δt1)::Float64 -
                bitbitll(Y[1,wci[j],k], Y[2,wci[j],k], 
                         λ[1], λ[2], ω1, ω0, ldci[j,k], δt1)::Float64
        end

      end

    end

    ll::Float64
  end

  return f
end





"""
    brll(brevs::Array{Float64,1}, λ1::Float64, λ0::Float64, si::Int64)

Return likelihood for a branch in continuous time.
"""
function brll(brevs::Array{Float64,1}, λ1::Float64, λ0::Float64, si::Int64)

  cst::Int64   = si 
  lb ::Int64   = endof(brevs)
  ll ::Float64 = 0.0

  if lb > 1 
    for i = Base.OneTo(lb-1)
      ll += evll(brevs[i], cst == 0 ? λ1 : λ0)::Float64
      cst = 1 - cst
    end
  end

  ll += nell(brevs[lb], cst == 0 ? λ1 : λ0)::Float64
  
  return ll::Float64
end




"""
    evll(t::Float64, λ::Float64)

Return log-likelihood for events.
"""
evll(t::Float64, λ::Float64) = @fastmath (log(λ) - (λ * t))::Float64





"""
    nell(t::Float64, λ::Float64)

Return log-likelihood for nonevents.
"""
nell(t::Float64, λ::Float64) = (-1 * λ * t)::Float64





"""
    allλpr(λc::Array{Float64,2}, λprior::Float64)

Return log-prior for all areas 
"""
function allλpr(λc    ::Array{Float64,1},
                λprior::Float64)

  pr::Float64 = 0.0
  
  for j in λc
    pr += logdexp(j, λprior)::Float64
  end

  return pr::Float64
end

