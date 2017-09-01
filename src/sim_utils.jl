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

vi = [0,1,1,1,0]






Yt = rand(0:1, 4, 3)
Xt = rand(4)





# number of species and areas
n, k = size(Yt)


arav  = zeros(k)   # area averages
liav  = zeros(n)   # lineage averages
aldif = zeros(n,k) # area specific lineage differences




### start

arav, liav = ar_lin_avg(Xt, Yt, arav, liav, n, k)

ar_lin_dif(Xt, arav, aldif, n, k)

### step

# biogeographic step

# trait step


# repeat





"""
    f_λ(λ::Float64, ω::Float64, Δx::Float64)

Estimate rates for area colonization/loss based on the difference between lineage traits and area averages.
"""
f_λ(λ::Float64, ω::Float64, Δx::Float64) = @fastmath λ * exp(ω*Δx)





"""
    ar_lin_dif(Xt::Array{Float64,1}, arav::Array{Float64,1}, aldif::Array{Float64,2})

Estimate differences between each each lineages trait and area averages.
"""
function ar_lin_dif(Xt::Array{Float64,1}, 
                    arav::Array{Float64,1}, 
                    aldif::Array{Float64,2},
                    n   ::Int64,
                    k   ::Int64)
  @inbounds begin

    for j in Base.OneTo(k), i in Base.OneTo(n)
      aldif[i,j] = abs(Xt[i] - arav[j])
    end

    return aldif
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

      arav[j] = aa/na
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
function biogeosam_1step(λ1::Float64, λ0::Float64, δt::Float64, v1::Array{Int64,1})

  nch = 0
  for i in eachindex(v1)
    if v1[i] == 0
      if rand() < λ1*δt
        setindex!(v1,1,i)
        nch += 1
      end
    else 
      if rand() < λ0*δt
        setindex!(v1,0,i)
        nch += 1 
      end
    end
  end

  return v1, nch
end





"""
    check_sam(obj::Tuple{Array{Int64,1},Int64})

Returns false if biogeographic step consists of *only one* change or if the species 
does *not* go globally extinct. 
"""
function check_sam(obj::Tuple{Array{Int64,1},Int64})

  @fastmath  begin

    s = 0
    for a in obj[1]
      s += a
    end

    if obj[2] > 1 || s == 0
      return true
    else
      return false
    end
  
  end
end










# """
#     biogeosam_branch

# """
# function biogeosam_branch(λ1v::Array{Float64,1}, 
#                           λ0v::Array{Float64,1}, 
#                           bδt::Array{Float64,1}, 
#                           vi ::Array{Int64,1})

#   for i in eachindex(λ1v)

#     vp = copy(vi)
#     vr = biogeosam_1step(λ1v[i], λ0v[i], bδt[i], vp)

#     while check_sam(vr)
#       vr = biogeosam_1step(λ1v[i], λ0v[i], bδt[i], vp)
#     end

#     vi = vr[1]
#   end


# end



