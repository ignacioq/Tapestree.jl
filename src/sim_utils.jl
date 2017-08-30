#=
Utility functions for simulations in Compete

Ignacio Quintero

August 29 2017

t(-_-t)

=#


"""
    biogeosam_1step(λ1::Float64, λ0::Float64, v1::Array{Int64,1})

Sample one step for biogeographic history.
"""
function biogeosam_1step(λ1::Float64, λ0::Float64, v1::Array{Int64,1})

  nch = 0
  for i in eachindex(v1)
    if v1[i] == 0
      if rand() < λ1
        setindex!(v1,1,i)
        nch += 1
      end
    else 
      if rand() < λ0
        setindex!(v1,0,i)
        nch += 1 
      end
    end
  end

  return v1, nch
end




"""
    biogeosam_1step(λ1::Float64, λ0::Float64, v1::Array{Int64,1})

Sample one step for biogeographic history.
"""
function check_sam(obj::Tuple{Array{Int64,1},Int64})
  
  s = 0
  for a in obj[1]
    s += a
  end

  if obj[2] > 1 || s == 0
    return false
  else
    return true
  end

end



v1
x = biogeosam_1step(0.5,0.5,v1); x
check_sam(x)









