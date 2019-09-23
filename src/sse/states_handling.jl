#=

handling of discrete states for ESSE and EGeoHiSSE

Ignacio Quintero Mächler

t(-_-t)

October 30 2017

=#


"""
    struct Sgh
      g::Set{Int64}
      h::Int64
    end

Composite type for Geographical State with `g` areas & `h` hidden states
"""
struct Sgh
  g::Set{Int64}
  h::Int64
end





"""
    isequal(x::Sgh, y::Sgh)

Compares equality between two `Sgh` types.
"""
isSghequal(x::Sgh, y::Sgh) = x.g == y.g && x.h == y.h 





"""
    create_states(k::Int64, h::Int64)

Create EGeoHiSSE states
"""
function create_states(k::Int64, h::Int64)

  # create individual areas subsets
  gs = Array{Set{Int64}, 1}(undef, k)
  for i = Base.OneTo(k)
    gs[i] = Set([i])
  end
  gS = sets(gs)
  popfirst!(gS)            # remove empty
  sort!(gS, by = length)   # arrange by range size

  # add hidden states and create Sgh objects
  S = Array{Sgh, 1}()
  for i = 0:(h-1), j = gS
    push!(S, Sgh(j, i))
  end

  return S
end




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
    build_par_names(k::Int64, h::Int64, ny::Int64, model::Array{Int64,1})

Build dictionary for parameter names and indexes for EGeoHiSSE for
`k` areas, `h` hidden states and `ny` covariates.
"""
function build_par_names(k    ::Int64, 
                         h    ::Int64, 
                         ny   ::Int64, 
                         model::NTuple{3,Bool})

  # generate individual area names
  ia = Array{Char,1}(undef,0)
  for i = 0:(k-1)
    push!(ia, ('A' + i)::Char)
  end

  ## build parameters name 
  par_nams = Array{String,1}(undef,0)

  # add λ names but only one between-regions speciation rate
  for j = 0:(h-1)
    for a::Char = ia
      push!(par_nams, "lambda_"*string(a)*"_"*string(j))
    end
    lastindex(ia) > 1 && push!(par_nams, "lambda_W_"*string(j))
  end

  # add μ names for endemics
  for j = 0:(h-1), i = ia
    push!(par_nams, "mu_"*string(i)*"_"*string(j))
  end

  # add area gains `g`
  # transitions can only through **one area** transition
  for i = 0:(h-1), a = ia, b = ia
    a == b && continue
    push!(par_nams, "gain_"*string(a)*string(b)*"_"*string(i))
  end

  # add area looses `l`
  # transitions can only through **one area** transition
  for i = 0:(h-1), a = ia
    push!(par_nams, "loss_"*string(a)*"_"*string(i))
  end

  # add q between hidden states
  for j = 0:(h-1), i = 0:(h-1)
    j == i && continue 
    push!(par_nams, "q_"*string(j)*string(i))
  end

  yppar = ny == 1 ? 1 : div(ny,
    model[1]*k + 
    model[2]*k + 
    model[3]*k*(k-1))

  # add betas
  if model[1]
    for i = 0:(h-1), a = ia, l = Base.OneTo(yppar)
      push!(par_nams, 
        "beta_lambda_"*string(l)*"_"*string(a)*"_"*string(i))
    end
  end
  if model[2]
    for i = 0:(h-1), a = ia, l = Base.OneTo(yppar)
      push!(par_nams, 
        "beta_mu_"*string(l)*"_"*string(a)*"_"*string(i))
    end
  end
  if model[3]
    for i = 0:(h-1), a = ia, b = ia, l = Base.OneTo(yppar)
      a == b && continue
      push!(par_nams, 
        "beta_q_"*string(l)*"_"*string(a)*string(b)*"_"*string(i))
    end
  end

  pardic::Dict{String, Int64} = Dict(par_nams[i]::String => i::Int64
    for i = Base.OneTo(lastindex(par_nams)))

  return pardic::Dict{String, Int64}
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
    set_constraints(constraints::NTuple{N,String},
                    pardic     ::Dict{String,Int64},
                    h          ::Int64)

Make a Dictionary linking parameter that are to be the same.
"""
function set_constraints(constraints::NTuple{N,String},
                         pardic     ::Dict{String,Int64},
                         k          ::Int64,
                         h          ::Int64,
                         ny         ::Int64,
                         model      ::NTuple{3,Bool}) where {N}

  # Dict of hidden state correspondence
  hsc = dict_hscor(k, h, ny, model)

  dcp  = Dict{Int64,Int64}()  # constrains dictionary for base hidden state 0
  dcfp = Dict{Int64,Int64}()  # constrains dictionary for other hidden states
  zp   = Int64[]              # zeros for base hidden state 0
  zfp  = Int64[]              # zeros for other hidden states
  chs  = Int64[]              # no hidden states left

  for c in constraints
    spl = map(x -> strip(x), split(c, '='))

    # if no equality
    if length(spl) < 2
      @warn "No equality found in $c"
      continue
    end

    # for parameters fixed to `0`
    if isequal(strip(spl[end]), "0")
      for i in Base.OneTo(length(spl)-1)
        spli = strip(spl[i])
        hsi  = split(spli, "_")[end]
        if hsi == 0 
          push!(zp, pardic[spli])
        else
          push!(zfp, pardic[spli])
          push!(zfp, pardic[spli])
        end
      end

    # for parameter equalities
    else
      for i in Base.OneTo(length(spl)-1) 
        # strip possible white spaces
        sp1 = spl[i]
        sp2 = spl[i+1]

        if haskey(hsc, pardic[sp1])
          # both have hidden states
          if haskey(hsc, pardic[sp2])
            dcp[hsc[pardic[sp2]]] = hsc[pardic[sp1]]
            dcfp[pardic[sp2]]     = pardic[sp1]
          # only 1 has hidden state
          else
            dcp[pardic[sp2]] = hsc[pardic[sp1]]
            push!(zfp, pardic[sp1])            
          end
        # only 2 has hidden state
        elseif haskey(hsc, pardic[sp2])
          dcp[pardic[sp1]] = hsc[pardic[sp2]]
          push!(zfp, pardic[sp2]) 
        # both are basal hidden state
        else
          dcp[pardic[sp2]] = pardic[sp1]
        end
      end
    end

    # no hidden states remaining for constraints?
    nhs = 0
    for i in 0:(h-1)
      re = Regex(".*_"*string(i)*"\$|.*_"*string(i)*" ")
      nhs += occursin(re, c) ? 1 : 0
    end
    if nhs == h
      for s in spl 
        push!(chs, pardic[s])
      end
    end
  end

  return dcp, dcfp, zp, zfp, chs
end





"""
    dict_hscor(k::Int64, h::Int64, ny::Int64)

Create dictionary of hidden states correspondence.
"""
function dict_hscor(k::Int64, h::Int64, ny::Int64, model::NTuple{3, Bool})

 # number of covariates
  yppar = ny == 1 ? 1 : div(ny,
    model[1]*k + 
    model[2]*k + 
    model[3]*k*(k-1))

  # starting indices for models 2 and 3
  m2s = model[1]*k*h
  m3s = m2s + model[2]*k*h

  hsc = Dict{Int64,Int64}()

  for j in 1:(h-1)
    for i in 1:k

      # speciation
      hsc[(k+1)*j + i] = (k+1)*(j-1) + i

      # extinction
      s = (k+1)*h 
      hsc[s + k*j + i] = s + k*(j-1) + i

      # gain
      s = (2k+1)*h
      for a in 1:(k-1)
        hsc[s + k*(k-1)*j + a + (k-1)*(i-1)] = s + k*(k-1)*(j-1) + a + (k-1)*(i-1)
      end

      # loss 
      s = (2k+1 + k*(k-1))*h
      hsc[s + k*j + i] = s + k*(j-1) + i

      # betas
      s = (3k+1+k*(k-1))*h + h*(h-1)

      if model[1]
        for l = Base.OneTo(yppar)
          hsc[s + k*j + i + (l-1)*(i-1)] = s + k*(j-1) + i + (l-1)*(i-1)
        end
      end
      if model[2]
        for l = Base.OneTo(yppar)
          hsc[s + k*j + m2s + i + (l-1)*(i-1)] = s + k*(j-1) + m2s + i + (l-1)*(i-1)
        end
      end

      if model[3]
        for a = 1:(k-1), l = Base.OneTo(yppar)
          hsc[s + k*(k-1)*j + m3s + a + (a-1)*yypar + (i-1)*(k-1)] = 
            s + k*(k-1)*(j-1) + m3s + a + (a-1)*yypar + (i-1)*(k-1)
        end
      end
    end

    # between-region speciation
    hsc[(k+1)*j + (k+1)] = (k+1)*(j-1) + (k+1)
  end

  return hsc
end





"""
    states_to_values(tipst::Dict{Int64,Int64}, S::Array{Sgh,1}, k::Int64, h::Int64)
    
Transform numbered tip_values to array with 1s and 0s
"""
function states_to_values(tipst::Dict{Int64,Int64}, S::Array{Sgh,1}, k::Int64)

  tip_val = Dict{Int64,Array{Float64,1}}()
  for (key, val) in tipst
    push!(tip_val, key => zeros(k))
    for s in S[val].g
      setindex!(tip_val[key], 1.0, s)
    end
  end

  return tip_val
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











