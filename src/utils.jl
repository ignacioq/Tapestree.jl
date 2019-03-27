#=

General utilities for ESSE

Ignacio Quintero Mächler

t(-_-t)

October 30 2017

=#



"""
  sets(x::Array{String,1})

Create all sets from each string (letters) in x.
"""
function sets(x::Array{String,1})
  ss = [""]
  for elem = x, j = eachindex(ss)
    push!(ss, ss[j]*elem)
  end

  return ss
end





"""
  sets(x::Array{SubString{String},1})

Create all subsets from each string (letters) in x.
"""
function sets(x::Array{SubString{String},1})
  ss = [""]
  for elem = x, j = eachindex(ss)
    push!(ss, ss[j]*elem)
  end

  return ss
end





"""
  sets(x::Vector{Set{Int64}})

Create all sets from a vector of sets of length `1`.
"""
function sets(x::Vector{Set{Int64}})

  ss = Vector{Set{Int64}}()
  push!(ss, Set{Int64}())
  for i = x, j = eachindex(ss)
    push!(ss, union(ss[j], i))
  end

  return ss
end




"""
   vicsubsets(x::String)

Create all vicariance subsets from a string of areas.
"""
function vicsubsets(x::String)
  lx  = lastindex(x)
  sx  = split(x, "")
  ss  = sets(sx)
  sort!(ss, by = length)
  pop!(ss)
  popfirst!(ss)
  ls = lastindex(ss)

  vss = Array{NTuple{2, String},1}()
  for i = Base.OneTo(ls)
    push!(vss, (ss[i], ss[ls-i+1]))
  end

  return vss
end





"""
    vicsubsets(x::Set{Int64})

Create all vicariance subsets from a Set of areas.
"""
function vicsubsets(x::Set{Int64})
  lx  = length(x)

  ss = Vector{Set{Int64}}()
  push!(ss, Set{Int64}())
  for i = x, j = eachindex(ss)
    push!(ss, union(ss[j], i))
  end
  sort!(ss, by = length)
  pop!(ss)
  popfirst!(ss)

  ls = lastindex(ss)
  vss = Array{NTuple{2, Set{Int64}},1}()
  for i = Base.OneTo(ls)
    push!(vss, (ss[i], ss[ls-i+1]))
  end

  return vss
end





"""
    build_par_names(k::Int64, h::Int64, ny::Int64, mod::NTuple{3,Bool})

Build dictionary for parameter names and indexes for EGeoHiSSE for
`k` areas, `h` hidden states and `ny` covariates.
"""
function build_par_names(k::Int64, h::Int64, ny::Int64, model::NTuple{3,Bool})

  # generate individual area names
  ia = String[]
  for i = 0:(k-1)
    push!(ia, string('A' + i))
  end

  ## build parameters name 
  par_nams = String[]

  # add λ names but only one between-regions speciation rate
  for j = 0:(h-1)
    for a = ia
      push!(par_nams, "lambda_$(a)_$j")
    end
    lastindex(ia) > 1 && push!(par_nams, "lambda_W_$j")
  end

  # add μ names for endemics
  for j = 0:(h-1), i = ia
    push!(par_nams, "mu_$(i)_$j")
  end

  # add area gains `g`
  # transitions can only through **one area** transition
  for i = 0:(h-1), a = ia, b = ia
    a == b && continue
    push!(par_nams, "gain_$(a)$(b)_$i")
  end

  # add area looses `l`
  # transitions can only through **one area** transition
  for i = 0:(h-1), a = ia
    push!(par_nams, "loss_$(a)_$i")
  end

  # add q between hidden states
  for j = 0:(h-1), i = 0:(h-1)
    j == i && continue 
    push!(par_nams, "q_$(j)$(i)")
  end

  ## add betas
  # if model is on Q
  if model[3]
    for i = 0:(h-1), a = ia, b = ia, l = Base.OneTo(div(ny,k*(k-1)))
      a == b && continue
      push!(par_nams, "beta$(l)_$(a)$(b)_$(i)")
    end
  # if model is speciation or extinction
  elseif model[1] || model[2]
    for i = 0:(h-1), a = ia, l = Base.OneTo(div(ny,k))
      push!(par_nams, "beta$(l)_$(a)_$i")
    end
  end

  pardic = Dict(par_nams[i] => i for i = Base.OneTo(lastindex(par_nams)))

  return pardic
end





"""
    build_par_names(k::Int64, T::Bool)

Build dictionary for parameter names and indexes for EHiSSE.
"""
function build_par_names(k::Int64, T::Bool)

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
    set_constraints(constraints::NTuple{length(constraints),String},
                    pardic     ::Dict{String,Int64})

Make a Dictionary linking parameter that are to be the same.
"""
function set_constraints(constraints::NTuple{length(constraints),String},
                         pardic     ::Dict{String,Int64})

  conpar = Dict{Int64,Int64}()
  zerov  = Int64[]

  for c in constraints
    spl = split(c, '=')
    if length(spl) < 2
      continue
    end
    if isequal(strip(spl[end]), "0")
      for i in Base.OneTo(length(spl)-1) 
        push!(zerov, pardic[strip(spl[i])])
      end
    else
      for i in Base.OneTo(length(spl)-1) 
        sp1 = strip(spl[i])
        sp2 = strip(spl[i+1])
        conpar[pardic[strip(spl[i])]] = pardic[strip(spl[i+1])]
      end
    end
  end

  return conpar, zerov
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








