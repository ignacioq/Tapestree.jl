#=
Utility functions for simulations in Compete

Ignacio Quintero

August 29 2017

t(-_-t)

=#



bs = 2.1     # branch start
bf = 4.2     # branch end
bt = bf - bs # branch time

bδt = bt/1000
λ1 = 0.3
λ0 = 0.1
ω1 = -0.1
ω0 = 0.1
δt = 0.01






Yt = rand(0:1, 4, 3)
Xt = rand(4)





# number of species and areas
n, k = size(Yt)


arav = zeros(k)   # area averages
liav = zeros(n)   # lineage averages
alλ1 = zeros(n,k)  # area specific lineage rates
alλ0 = zeros(n,k)  # area specific lineage rates




### start

arav, liav = ar_lin_avg(Xt, Yt, arav, liav, n, k)

alλ1, alλ0 = lin_rates(Xt, arav, alλ1, alλ0, λ1, λ0, ω1, ω0, n, k)



### step

# biogeographic step



Ytn = copy(Yt)
biogeosam_1step(alλ1, alλ0, δt, Yt, n, k)





  for i in eachindex(λ1v)

    vp = copy(vi)
    vr = biogeosam_1step(λ1v[i], λ0v[i], bδt[i], vp)

    while check_sam(vr)
      vr = biogeosam_1step(λ1v[i], λ0v[i], bδt[i], vp)
    end

    vi = vr[1]
  end


# trait step


# repeat



check_sam(Ytc::Array{Int64,2}, nch::Array{Int64,1}, n::Int64, k::Int64)



"""
    f_λ(λ::Float64, ω::Float64, Δx::Float64)

Estimate rates for area colonization/loss based on the difference between lineage traits and area averages.
"""
f_λ(λ::Float64, ω::Float64, Δx::Float64) = @fastmath λ * exp(ω*Δx)





"""
    lin_rates(Xt  ::Array{Float64,1}, arav::Array{Float64,1}, alλ1::Array{Float64,2}, alλ0::Array{Float64,2}, λ1::Float64, λ0::Float64, ω1::Float64, ω0::Float64, n::Int64, k::Int64)

Estimate lineage specific area colonization/loss rates based on differences between each each lineages trait and area averages.
"""
function lin_rates(Xt  ::Array{Float64,1}, 
                   arav::Array{Float64,1}, 
                   alλ1::Array{Float64,2},
                   alλ0::Array{Float64,2},
                   λ1  ::Float64,
                   λ0  ::Float64,
                   ω1  ::Float64,
                   ω0  ::Float64,
                   n   ::Int64,
                   k   ::Int64)
  @inbounds begin

    for j in Base.OneTo(k), i in Base.OneTo(n)
      Δx = abs(Xt[i] - arav[j])
      alλ1[i,j] = f_λ(λ1, ω1, Δx)
      alλ0[i,j] = f_λ(λ0, ω0, Δx)
    end

    return alλ1, alλ0
  end
end





"""
    ar_lin_avg(Xt::Array{Float64,1}, Yt::Array{Int64,2}, arav::Array{Float64,1}, 
              liav::Array{Float64,1}, n::Int64, k::Int64)

Estimate area and lineage specific averages given sympatry configuration.
"""
function ar_lin_avg(Xt  ::Array{Float64,1}, 
                    Yt  ::Array{Int64,2}, 
                    arav::Array{Float64,1},
                    liav::Array{Float64,1},
                    n   ::Int64,
                    k   ::Int64)

  @inbounds begin
    # estimate area averages
    for j in Base.OneTo(k)

      aa = 0.0
      na = 0
      for i in Base.OneTo(n)
        aa += Yt[i,j]*Xt[i]
        na += Yt[i,j]
      end

      arav[j] = aa/(na == 0 ? 1 : na) 
    end

    # estimate lineage averages
    for i in Base.OneTo(n)

      la = 0.0
      na = 0 
      for j in Base.OneTo(k)
        la += arav[j]*Yt[i,j]
        na += Yt[i,j]
      end
      
      liav[i] = la/na
    end

  end

  return arav, liav
end






"""
    traitsam_1step(xi::Float64, μ::Float64, δt::Float64, ωx::Float64, σ::Float64)

Sample one step for trait evolution history: X(t + δt).
"""
traitsam_1step(xi::Float64, μ::Float64, δt::Float64, ωx::Float64, σ::Float64) = 
  xi += ωx*(μ - xi)*δt + randn()*σ*sqrt(δt)





"""
    biogeosam_1step(λ1::Float64, λ0::Float64, δt::Float64, v1::Array{Int64,1})

Sample one step for biogeographic history Y(t + δt).
"""
function biogeosam_1step(alλ1 ::Array{Float64,2}, 
                         alλ0 ::Array{Float64,2}, 
                         δt   ::Float64, 
                         Ytn  ::Array{Int64,2},
                         n    ::Int64,
                         k    ::Int64)

  @inbounds begin

    nch = zeros(Int64, n)
    for j in Base.OneTo(k), i in Base.OneTo(n)
      if Ytn[i,j] == 0
        if rand() < alλ1[i,j]*δt
          setindex!(Ytn,1,i,j)
          nch[i] += 1
        end
      else 
        if rand() < alλ0[i,j]*δt
          setindex!(Ytn,0,i,j)
          nch[i] += 1 
        end
      end
    end

  end

  return Ytn, nch
end





"""
    check_sam(Ytc::Array{Int64,2}, nch::Array{Int64,1}, n::Int64, k::Int64)

Returns false if biogeographic step consists of *only one* change or if the species 
does *not* go globally extinct. 
"""
function check_sam(Ytc::Array{Int64,2}, nch::Array{Int64,1}, n::Int64, k::Int64)

  @fastmath @inbounds begin

    # check lineage did not go extinct
    for i in Base.OneTo(n)
      s = 0
      for j in Base.OneTo(k)
        s += Ytc[i,j]
      end
      if s == 0
        return true
      end
    end


    # check that, at most, only one event happened per lineage
    for i in nch
      if i > 1
        return true
      end
    end

    return false

  end
end









