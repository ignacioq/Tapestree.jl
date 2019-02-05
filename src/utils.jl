#=

General utilities for ESSE

Ignacio Quintero Mächler

t(-_-t)

October 30 2017

=#



"""
    build_par_names(k::Int64, T::Bool)

Build dictionary for parameter names and indexes.
"""
function build_par_names(k::Int64, T::Bool = true)

  par_nams = String[]
  # build parameters name 
  for i in 0:(k-1)
    push!(par_nams, "lambda$i")
  end

  # add μ names
  for i in 0:(k-1)
    push!(par_nams, "mu$i")
  end

  # add q names
  for j in 0:(k-1), i in 0:(k-1)
    if i == j 
      continue
    end
    push!(par_nams, "q$j$i")
  end

  # add betas
  if T
    for i in 0:(k-1)
      push!(par_nams, "beta$i")
    end
  else
    # add betas
    for j in 0:(k-1), i in 0:(k-1)
      if i == j 
        continue
      end
      push!(par_nams, "beta$j$i")
    end
  end

  pardic = Dict(par_nams[i] => i for i in Base.OneTo(length(par_nams)))

  return pardic
end




"""
    build_par_names(k::Int64)

Build dictionary for parameter names and indexes.
"""
function build_par_names(k::Int64)

  par_nams = String[]
  # build parameters name 
  for i in 0:(k-1)
    push!(par_nams, "lambda$i")
  end

  # add μ names
  for i in 0:(k-1)
    push!(par_nams, "mu$i")
  end

  # add q names
  for j in 0:(k-1), i in 0:(k-1)
    if i == j 
      continue
    end
    push!(par_nams, "q$j$i")
  end

  pardic = Dict(par_nams[i] => i for i in Base.OneTo(length(par_nams)))

  return pardic
end






"""
    set_constraints(constraints::NTuple{endof(constraints),String},
                         pardic::Dict{String,Int64})

Make a Dictionary linking parameter that are to be the same.
"""
function set_constraints(constraints,
                         pardic     ::Dict{String,Int64})

  @inbounds begin
    conpar = Dict{Int64,Int64}()

    for c in constraints
      spl = split(c, '=')
      if length(spl) < 2
        continue
      end
      conpar[pardic[strip(spl[1])]] = pardic[strip(spl[2])]
    end
  end
  return conpar
end






"""
    states_to_values(tipst::Dict{Int64,Int64}, k::Int64)
    
Transform numbered tip_values to array with 1s and 0s
"""
function states_to_values(tipst::Dict{Int64,Int64}, k::Int64)
  tip_val = Dict{Int64,Array{Float64,1}}()
  for (key, val) in tipst
    push!(tip_val, key => zeros(k))
    setindex!(tip_val[key], 1.0, val)
  end

  return tip_val
end







