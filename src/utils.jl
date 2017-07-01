using ProgressMeter

#=
Utility functions for turnover

Ignacio Quintero

February 6 2017

t(-_-t)

=#



function makescalef(obj_ar::Float64)
  nar = 1.0 - obj_ar

  function f(window::Float64, rate::Float64)
    if (rate > obj_ar)
      window *= (1. + (rate - obj_ar) / nar)
    else
      window /= (2. - rate / obj_ar)
    end
    window
  end

  return f
end



# uniform update
uniupt(p::Float64, tn::Float64) = abs(p + (rand()-0.5) * tn)



# additive update
addupt(p::Float64, tn::Float64) = p + randn() * tn



# additive update
absaddupt(p::Float64, tn::Float64) = abs(p + randn() * tn)



# multiplicative update
mulupt(p::Float64, tn::Float64) = p * exp((rand() - 0.5) * tn)




# logarithmic update
logupt(p::Float64, tn::Float64) = exp(log(p) + (rand() - 0.5) * tn)




# get row indexing from matrix indexing
rowind(x::Int64, nrow::Int64) = mod1(x,nrow)




# get column indexing from matrix indexing
colind(x::Int64, nrow::Int64) = cld(x, nrow)



# get vector indexing from column and row
vecind(row::Int64, col::Int64, nrow::Int64) = row + nrow*(col - 1)




# global scaling factor
function makeglobalscalef(obj_ar::Float64)

  function f(λ::Float64, globalrate::Float64, stepsize::Float64)
    newλ::Float64 = log(λ) + stepsize * (globalrate - obj_ar)
    exp(newλ)
  end

  return f
end




# paramete mean and Σ adaptive update
function adaptiveupd(Σ       ::Array{Float64,2},
                     psample ::Array{Float64,1},
                     pmean   ::Array{Float64,1},
                     stepsize::Float64)

  pdif  ::Array{Float64,1} = psample - pmean
  pmeanN::Array{Float64,1} = pmean + stepsize * pdif
  ΣN    ::Array{Float64,2} = Σ + stepsize *
                             (BLAS.ger!(1., pdif, pdif, zeros(2,2)) - Σ)

  pmeanN, ΣN
end




# make function for the stepsize for the adaptive update
function makestepsize(C::Float64, η::Float64)
  
  β::Float64 = rand(linspace((1./(1. + η)),1))

  function f(t::Int64)
    C/(t^β)
  end

  return f
end




# multivariate update
function makemvnproposal(Σ::Array{Float64,2})

  evc  = eigvecs(Σ)
  evl  = diagm(eigvals(Σ))
  spde = evc * sqrt(evl)

  function f(pvec::Vector{Float64})
    pvec + spde * randn(2)
  end

  return f
end





# return index for closest value 
# in non-sorted arrays 
function indmindif(x::Array{Float64,1}, val::Float64) 
  ibest  = start(eachindex(x)) 
  dxbest = abs(x[ibest]-val) 

  for j in eachindex(x) 
    dx = abs(x[j]-val) 
    if dx < dxbest 
        dxbest = dx 
        ibest  = j 
    end 
  end 
  ibest 
end 

