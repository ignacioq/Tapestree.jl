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



transpose(Yt) * Xt ./ reducedim(+, Yt, 1)







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



