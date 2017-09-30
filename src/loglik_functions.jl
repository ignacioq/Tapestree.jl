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
  if ωx < 0.0
    Δx::Float64 = μ - xi
    if Δx < 0.0
      return (-1 * ωx * exp(Δx) * δt)::Float64
    else
      return (ωx * exp(-Δx) * δt)::Float64
    end
  else
     return (ωx * (μ - xi) * δt)::Float64
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
    makellf(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64, wcol::Vector{Vector{Int64}}, narea::Int64)

Make likelihood function for all trait matrix 
and biogeography history.
"""
function makellf(δt   ::Vector{Float64}, 
                 Y    ::Array{Int64, 3}, 
                 ntip ::Int64, 
                 wcol ::Vector{Vector{Int64}},
                 narea::Int64)

  const coloop = Base.OneTo(size(Y,2))
  
  # which is 23 (23 = NaN) in each column
  const w23 = Array{Int64,1}[]
  for i = coloop
    push!(w23,find(Y[:,i,1] .!= 23))
  end

  # to avoid re-subsetting add 
  # a dummy at the end of δt
  const dδt = copy(δt)
  push!(dδt,0.)

  # number of normal evaluations
  n = 0
  for i in w23
    n += (length(i) - 1)
  end

  # normal constant
  const normC = -0.5*log(2.0π)*n

  function f(X      ::Array{Float64,2},
             Y      ::Array{Int64,3}, 
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             ωx     ::Float64,
             ωλ     ::Float64,
             ωμ     ::Float64,
             λ      ::Array{Float64,1},
             stemevc::Vector{Vector{Float64}},
             stemss,
             σ²     ::Float64)

    ll::Float64 = normC

    @inbounds begin

      # trait likelihood
      for j=Base.OneTo(ntip), i=w23[j][Base.OneTo(end-1)]

        ll += -0.5*log(δt[i]*σ²) -
              abs2(X[(i+1),j] -
                  (X[i,j] + E_sde(X[i,j], linavg[i,j], ωx, δt[i])))/
              (2.0*δt[i]*σ²)
      end

      # biogeograhic likelihood
      for j in Base.OneTo(narea)
        ll += brll(stemevc[j], λ[1], λ[2], stemss[j])
        for i = coloop
          wh = w23[i]
          ll += bitvectorll(Y[wh,i,j], λ[1], λ[2], ωλ, ωμ, 
                            lindiff[wh,i,j], dδt[wh]) 
        end
      end

    end

    ll
  end

  return f
end




"""
    bitvectorll(y ::Array{Int64,1}, λ1::Float64, λ0::Float64, ωλ::Float64, ωμ::Float64, Δx::Array{Float64,1}, δt::Array{Float64,1})

Return the likelihood for a
bit vector (composed of 0s and 1s).
"""
function bitvectorll(y ::Array{Int64,1},
                     λ1::Float64,
                     λ0::Float64,
                     ωλ::Float64,
                     ωμ::Float64,
                     Δx::Array{Float64,1},
                     δt::Array{Float64,1})

  ll::Float64 = 0.0

  @inbounds begin

    cur_s::Int64 = y[1]

    if cur_s == 0 
      cur_λ, cur_ω = λ1, ωλ
    else
      cur_λ, cur_ω = λ0, ωμ
    end

    for i=Base.OneTo(endof(y)-1)

      if y[i] == y[i+1]
        ll += nell(δt[i], ratest(cur_λ, cur_ω, Δx[i]))
      else
        ll += evll(δt[i], ratest(cur_λ, cur_ω, Δx[i]))
        cur_s = 1 - cur_s

        if cur_s == 0 
          cur_λ, cur_ω = λ1, ωλ
        else
          cur_λ, cur_ω = λ0, ωμ
        end
      end
    end

  end

  return ll
end



"""
    bitbitll(y1::Int64, y2::Int64, λ1::Float64, λ0::Float64, ωλ::Float64, ωμ::Float64, Δx::Float64, δt::Float64)

Return the likelihood for a
bit vector of length 2 (composed of 0s and 1s).
"""
function bitbitll(y1::Int64,
                  y2::Int64,
                  λ1::Float64,
                  λ0::Float64,
                  ωλ::Float64,
                  ωμ::Float64,
                  Δx::Float64,
                  δt::Float64)

  # event or non-event
  if y1 == y2 
    if y1 == 0 
      return nell(δt, ratest(λ1, ωλ, Δx))
    else
      return nell(δt, ratest(λ0, ωμ, Δx))
    end
  else
    if y1 == 0 
      return evll(δt, ratest(λ1, ωλ, Δx))
    else
      return evll(δt, ratest(λ0, ωμ, Δx))
    end
  end
end



"""
    ratest(λ::Float64, ω::Float64, absdiff::Float64)

Estimate λ colonization/extirpation rate
based on the absolute difference on X.
"""
ratest(λ::Float64, ω::Float64, absdiff::Float64) = 
  λ * exp(ω * absdiff)





"""
    makellf_λ_upd(Y::Array{Int64,3}, δt::Vector{Float64}, narea::Int64)

Make likelihood function for when updating λ.
"""
function makellf_λ_upd(Y    ::Array{Int64,3},
                       δt   ::Vector{Float64},
                       narea::Int64)

  const coloop = Base.OneTo(size(Y,2))

  # which is 23 (23 = NaN) in each column
  const w23 = Array{Int64,1}[]
  for i=coloop
    push!(w23,find(Y[:,i,1] .!= 23))
  end

  # to avoid re-subsetting add 
  # a dummy at the end of δt
  const dδt = copy(δt)
  push!(dδt,0.)

  function f(Y      ::Array{Int64,3}, 
             λ::Array{Float64,1},
             ωλ     ::Float64,
             ωμ     ::Float64,
             lindiff::Array{Float64,3},
             stemevc::Vector{Vector{Float64}},
             stemss ::Vector{Int64})

    ll::Float64 = 0.0

    @inbounds begin

      for j=Base.OneTo(narea)

        ll += brll(stemevc[j], λ[1], λ[2], stemss[j])

        for i=coloop
          wh = w23[i]
          ll += bitvectorll(Y[wh,i,j], λ[1], λ[2], ωλ, ωμ, 
                            lindiff[wh,i,j], dδt[wh]) 
        end
      end
    
    end

    return ll
  end

  return f
end




"""
    makellf_ωλμ_upd(Y::Array{Int64,3}, δt::Vector{Float64}, narea::Int64)

Make likelihood function for when updating ωλ & ωμ.
"""
function makellf_ωλμ_upd(Y   ::Array{Int64,3},
                         δt   ::Vector{Float64},
                         narea::Int64)

  const coloop = Base.OneTo(size(Y,2))

  # which is 23 (23 = NaN) in each column
  const w23 = Array{Int64,1}[]
  for i=coloop
    push!(w23,find(Y[:,i,1] .!= 23))
  end

  # to avoid re-subsetting add 
  # a dummy at the end of δt
  const dδt = copy(δt)
  push!(dδt,0.)

  function f(Y      ::Array{Int64,3}, 
             λ::Array{Float64,1},
             ωλ     ::Float64,
             ωμ     ::Float64,
             lindiff::Array{Float64,3})

    ll::Float64 = 0.0

    @inbounds begin

      for j=Base.OneTo(narea), i=coloop
        wh = w23[i]
        ll += bitvectorll(Y[wh,i,j], λ[1], λ[2], ωλ, ωμ, 
                          lindiff[wh,i,j], dδt[wh]) 
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
function makellf_biogeo_upd_iid(bridx_a::Array{Array{Array{Int64,1},1},1},
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


  function f(Y    ::Array{Int64,3}, 
             λ::Array{Float64,1},
             triad::Array{Int64,1})

    ll::Float64 = 0.0

    @inbounds begin

      pr, d1, d2 = triad

      if pr < nedge 
        for j=Base.OneTo(narea)
          ll += bitvectorll_iid(Y[bridx_a[j][pr]], λ[1], λ[2], δtA[pr]) +
                bitvectorll_iid(Y[bridx_a[j][d1]], λ[1], λ[2], δtA[d1]) +
                bitvectorll_iid(Y[bridx_a[j][d2]], λ[1], λ[2], δtA[d2])
        end
      else 
        for j=Base.OneTo(narea)
          ll += bitvectorll_iid(Y[bridx_a[j][d1]], λ[1], λ[2], δtA[d1]) +
                bitvectorll_iid(Y[bridx_a[j][d2]], λ[1], λ[2], δtA[d2])
        end
      end

    end

    return ll
  end

  return f
end



"""
    bitvectorll_iid(y::Array{Int64,1}, λ1::Float64, λ0::Float64, δt::Array{Float64,1})

Return likelihood under the independence model 
for a bit vector.
"""
function bitvectorll_iid(y ::Array{Int64,1},
                         λ1::Float64,
                         λ0::Float64,
                         δt::Array{Float64,1})

  ll::Float64 = 0.0

  @inbounds begin

    cur_s::Int64   = y[1]
    cur_λ::Float64 = cur_s == 0 ? λ1 : λ0

    for i=Base.OneTo(endof(y)-1)
      if y[i] == y[i+1]
        ll += nell(δt[i], cur_λ)
      else
        ll += evll(δt[i], cur_λ)
        cur_s = 1 - cur_s
        cur_λ = cur_s == 0 ? λ1 : λ0
      end
    end

  end

  ll
end





"""
    makellf_σ²ωxupd(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64)

Make likelihood function for all trait matrix, `X`.
"""
function makellf_σ²ωxupd(δt  ::Vector{Float64}, 
                         Y   ::Array{Int64, 3}, 
                         ntip::Int64)

  # which is 23 (i.e., NaN) in each column
  const w23 = Array{Int64}[]
  for i=Base.OneTo(size(Y,2))
    push!(w23,find(Y[:,i,1] .!= 23)[Base.OneTo(end-1)])
  end

  # number of normal evaluations
  n = 0
  for i=w23
    n += length(i)
  end

  # normal constant
  const normC = -0.5*log(2.0π)*n

  function f(X ::Array{Float64,2},
             la::Array{Float64,2},
             ωx::Float64,
             σ²::Float64)

    ll::Float64 = normC

    @inbounds begin

      # trait likelihood
      for j=Base.OneTo(ntip), i=w23[j]

        ll += -0.5*log(δt[i]*σ²) -
              abs2(X[(i+1),j] -
                  (X[i,j] + E_sde(X[i,j], la[i,j], ωx, δt[i])))/
              (2.0*δt[i]*σ²)

      end
    
    end

    ll
  end

  return f
end



"""
    makellf_Xupd(δt::Vector{Float64}, narea::Int64)

Make likelihood function for an internal node update in `X`.
"""
function makellf_Xupd(δt   ::Vector{Float64}, 
                      narea::Int64)

  function f(k    ::Int64,
             wck  ::Array{Int64,1},
             wckm1::Array{Int64,1},
             X    ::Array{Float64,2},
             Y    ::Array{Int64,3},
             lak  ,
             lakm1::Array{Float64,1},
             ldk,
             ωx   ::Float64,
             ωλ   ::Float64,
             ωμ   ::Float64,
             λ::Array{Float64,1},
             σ²   ::Float64)

    # normal likelihoods
    ll::Float64 = 0.0

    @inbounds begin

      # loop for parent nodes
      if k != 1               # if not the root
        for i in eachindex(wckm1)
          wci = wckm1[i]

          ll += -0.5*log(δt[k-1]*σ²) -
                abs2(X[k,wci] -
                    (X[k-1,wci] + E_sde(X[k-1,wci], lakm1[i], ωx, δt[k-1]))
                )/(2.0*δt[k-1]*σ²)
        end
      end

      # loop for daughter nodes
      for i=eachindex(wck)
        wci = wck[i]

        # trait likelihood
        ll += -0.5*log(δt[k]*σ²) -
              abs2(X[(k+1),wci] -
                  (X[k,wci] + E_sde(X[k,wci], lak[i], ωx, δt[k]))
              )/(2.0*δt[k]*σ²)

        # biogeograhic likelihoods
        for j=Base.OneTo(narea)
          ll += bitbitll(Y[k,wci,j], Y[k+1,wci,j], 
                          λ[1], λ[2], ωλ, ωμ, ldk[i,j], δt[k])
        end
      end

    end

    ll
  end

  return f
end




"""
    makellf_Rupd(δt::Vector{Float64}, narea::Int64)

Make likelihood function
for the root update in `X`.
"""
function makellf_Rupd(δt   ::Vector{Float64}, 
                      narea::Int64)

  function f(k    ::Int64,
             wck  ::Array{Int64,1},
             X    ::Array{Float64,2},
             Y    ::Array{Int64,3},
             lak,
             ldk,
             ωx   ::Float64,
             ωλ   ::Float64,
             ωμ   ::Float64,
             λ::Array{Float64,1},
             σ²   ::Float64)

    # normal likelihoods
    ll::Float64 = 0.0

    @inbounds begin

      # loop for daughter nodes
      for i=eachindex(wck)
        wci = wck[i]

        # trait likelihood
        ll += -0.5*log(δt[k]*σ²) -
              abs2(X[(k+1),wci] -
                  (X[k,wci] + E_sde(X[k,wci], lak[i], ωx, δt[k]))
              )/(2.0*δt[k]*σ²)

        # biogeograhic likelihoods
        for j = Base.OneTo(narea)
          ll += bitbitll(Y[k,wci,j], Y[k+1,wci,j], 
                         λ[1], λ[2], ωλ, ωμ,ldk[i,j], δt[k])
        end
      end

    end

    ll
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
    for i=Base.OneTo(lb-1)
      ll += evll(brevs[i], cst == 0 ? λ1 : λ0)
      cst = 1 - cst
    end
  end

  ll += nell(brevs[lb], cst == 0 ? λ1 : λ0)
  
  return ll
end




"""
    evll(t::Float64, λ::Float64)

Return log-likelihood for events.
"""
evll(t::Float64, λ::Float64) = @fastmath log(λ) - (λ * t)




"""
    nell(t::Float64, λ::Float64)

Return log-likelihood for nonevents.
"""
nell(t::Float64, λ::Float64) = -1 * λ * t




"""
    allλpr(λc::Array{Float64,2}, λprior::Float64)

Return log-prior for all areas 
"""
function allλpr(λc    ::Array{Float64,1},
                λprior::Float64)

  pr::Float64 = 0.0
  
  for j in λc
    pr += logdexp(j, λprior)
  end

  return pr
end

