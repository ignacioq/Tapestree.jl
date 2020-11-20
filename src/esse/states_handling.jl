#=

Handling of discrete states and parameters for ESSE

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

Build dictionary for parameter names and indexes for `ESSE.g` for
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
    lastindex(ia) > 1 && push!(par_nams, "loss_"*string(a)*"_"*string(i))
  end

  # add q between hidden states
  for j = 0:(h-1), i = 0:(h-1)
    j == i && continue 
    push!(par_nams, "q_"*string(j)*string(i))
  end

  if any(model)

    yppar = ny == 1 ? 1 : ceil(Int64,ny/
      (model[1]*k + 
       model[2]*k + 
       model[3]*k*(k-1)))

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
  end

  pardic = Dict(v => i for (i,v) = enumerate(par_nams))

  return pardic::Dict{String, Int64}
end





"""
    build_par_names(k::Int64, T::Bool)

Build dictionary for parameter names and indexes for `ESSE.d`.
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

  dcp0  = Dict{Int64,Int64}()  # constrains dictionary for parameters
  zp    = Set(Int64[])         # zeros for parameters
  zfp   = Set(Int64[])         # zeros for λ hidden states

  for c in constraints

    spl = map(x -> strip(x), split(c, '='))

    # if no equality
    if length(spl) < 2
      @warn "No equality found in $c, no constraints applied for this expression"
      continue
    end

    # if hidden state rates `q_ij`
    if occursin("q", spl[1])
      for i in Base.OneTo(length(spl)-1) 
        dcp0[pardic[spl[i+1]]] = pardic[spl[i]]
      end

      continue
    end

    # for parameters fixed to `0`
    if isequal(spl[end], "0")
      for i in Base.OneTo(length(spl)-1)
        hsi  = split(spl[i], "_")[end]
        if hsi != "0" && occursin(r"^lambda_.*", spl[i])
          push!(zfp, pardic[spl[i]])
        else
          push!(zp, pardic[spl[i]])
        end
      end

    # for constraints
    else
      for i in Base.OneTo(length(spl)-1) 

        sp1 = spl[i]
        sp2 = spl[i+1]

        # if not λ
        if !occursin(r"^lambda_.*", sp1) && !occursin(r"^lambda_.*", sp2)
          dcp0[pardic[sp1]] = pardic[sp2]

        else
          spl1 = split(sp1, "_")
          spl2 = split(sp2, "_")

          # concatenate without hidden state
          spc1 = spc2 = ""
          for s in Base.OneTo(lastindex(spl1)-1)
            spc1 *= spl1[s]*"_"
          end
          for s in Base.OneTo(lastindex(spl2)-1)
            spc2 *= spl2[s]*"_"
          end

          # convert to hidden states
          hs1 = parse(Int64, spl1[end])
          hs2 = parse(Int64, spl2[end])

          # minmax for hidden states
          mnh, mxh = minmax(hs1, hs2)

          # which the maximum
          spmn, spmx = hs1 < hs2 ? (sp1, sp2) : (sp2, sp1)

          # set non-shared hidden states to 0
          wp = pardic[spmx]
          for j in (mnh+1):mxh
            push!(zfp, wp)
            wp = hsc[wp]
          end
        end
      end
    end
  end

  # check double feedback loops
  dcp = Dict{Int64,Int64}()
  for (key, value) in dcp0
    if in(value, keys(dcp)) && key == dcp[value]
      continue
    end
    dcp[key] = value
  end

  pnv = collect(values(pardic))
  pnk = collect(keys(pardic))
  pnk = pnk[sortperm(pnv)]

  if length(dcp) > 0
    ss = "Enforced parameter equalities: \n"
    for (k,v) in dcp
      ss *= "$(pnk[k]) = $(pnk[v]) \n"
    end
    @info ss
  end

  if length(zp) > 0 || length(zfp) > 0
    ss = "Parameters set to 0: \n"
    for k in zp
      ss *= "$(pnk[k]) \n"
    end
    for k in zfp
      ss *= "$(pnk[k]) \n"
    end

    @info ss
  end


  return dcp, zp, zfp
end





"""
    dict_hscor(k::Int64, h::Int64, ny::Int64)

Create dictionary of hidden states correspondence for λ.
"""
function dict_hscor(k::Int64, h::Int64, ny::Int64, model::NTuple{3, Bool})

  hsc = Dict{Int64,Int64}()

  if isone(k)

    for j in 1:(h-1)
      # speciation
      s = 0
      hsc[s + j + 1] = s + j
    end

  else

    for j in 1:(h-1)
      for i in 1:k
        # within-region speciation
        hsc[(k+1)*j + i] = (k+1)*(j-1) + i
      end

      # between-region speciation
      hsc[(k+1)*j + (k+1)] = (k+1)*(j-1) + (k+1)
    end

  end

  return hsc
end





"""
  set_multivariate(mvpars::NTuple{N,String},
                   pardic::Dict{String,Int64}) where {N}

Returns all vectors for multivariate updates given mvpars.
"""
function set_multivariate(mvpars::NTuple{N,String},
                          pardic::Dict{String,Int64}) where {N}

  mvps = Array{Int64,1}[]
  for c in mvpars
    spl = map(x -> strip(x), split(c, '='))
    si = spl[1]

    if si == ""
      break
    end

    for (k,v) in pardic
      # if matches a parameter
      if occursin(Regex("^$si"), k)
        mvp = [v]
        # get area of parameter
        ar = match(r"[A-Z]", k).match 
        occursin('W', ar) && continue
        # get hidden state
        hd = match(r"[0-9]$", k).match 
        # get other parameters in c that has this area and hidden state
        for sj in spl
          si == sj && continue
          for (kj, vj) in pardic
            if occursin(Regex("^$sj"), kj) && occursin(ar*"_"*hd, kj)
              push!(mvp, vj)
            end
          end
        end
        push!(mvps, sort!(mvp))
      end
    end
  end

  return mvps
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











