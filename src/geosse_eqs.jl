#=

GeoSSE equations

Ignacio Quintero Mächler

t(-_-t)

Created 15 02 2019
=#





"""
    make_geosse(k::Int64)

GeoSSE ODE equation for k areas.
"""
function make_geosse(k::Int64)

  # n states
  ns = 2^k - 1

  #create subsets
  sa = Array{String,1}(undef, k)
  for i in Base.OneTo(k)
    sa[i] = string('A' + (i-1))
  end
  S = sets(sa)
  popfirst!(S)            # remove empty
  sort!(S, by = length)   # arrange by range size

  # Metaprogram to generate GeoSSE equations
  eqs = quote end

  for ri = Base.OneTo(ns)

    # range
    r = S[ri]

    # length of r
    lr = lastindex(r)

    # which single areas occur in r
    ia = findall(x -> occursin(x, r), sa)
    oa = setdiff(1:k, ia)

    #= 
    likelihoods
    =#

    # no events
    nev = noevents_expr(ri, lr, ia, oa, k, ns, false)

    # local extinction
    # remove if !isone(lr)
    lex = localext_expr(r, ia, sa, S, k)

    # dispersal
    # remove if lr == k
    dis = dispersal_expr(r, lr, ia, oa, S, ns, k, false)

    # within-region speciation
    wrs = wrspec_expr(ri, ia, ns)

    # between-region speciation
    # remove if !isone(lr)
    brs = brspec_expr(r, S, ns, k)

    # if single area
    if ri <= k
      push!(eqs.args, 
        :(du[$ri] = $nev + $dis + $wrs))
    # if widespread
    elseif ri != ns
      push!(eqs.args, 
        :(du[$ri] = $nev + $lex + $dis + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$ri] = $nev + $lex + $wrs + $brs))
    end

    #= 
    extinctions
    =#

    # no events
    nev = noevents_expr(ri, lr, ia, oa, k, ns, true)

    # extinction
    ext = ext_expr(r, ia, sa, S, ns, k) 

    # dispersal and extinction
    dis = dispersal_expr(r, lr, ia, oa, S, ns, k, true)

    # within-region extinction
    wrs = wrsext_expr(ri, ia, ns)

    # between-region extinction
    brs = brsext_expr(r, S, ns, k)

    # if single area
    if ri <= k
      push!(eqs.args, 
        :(du[$(ri + ns)] = $nev + $ext + $dis + $wrs))
    # if widespread
    elseif ri != ns
      push!(eqs.args, 
        :(du[$(ri + ns)] = $nev + $ext + $dis + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$(ri + ns)] = $nev + $ext + $wrs + $brs))
    end

  end
  
  ## aesthetic touches
  # remove REPL comment
  popfirst!(eqs.args)

  # reorder with likelihoods first
  permute!(eqs.args, append!([1:2:2ns...],[2:2:2ns...]))

  return eqs
end







"""
    make_geosse(k::Int64)

GeoSSE ODE equation for k areas.
"""
function make_geosse(k::Int64)

  # n states
  ns = 2^k - 1

  #create subsets
  sa = Array{String,1}(undef, k)
  for i in Base.OneTo(k)
    sa[i] = string('A' + (i-1))
  end
  S = subsets(sa)
  popfirst!(S)            # remove empty
  sort!(S, by = length)   # arrange by range size

  # Metaprogram to generate GeoSSE equations
  eqs = quote end

  for ri = Base.OneTo(ns)

    # range
    r = S[ri]

    # length of r
    lr = lastindex(r)

    # which single areas occur in r
    ia = findall(x -> occursin(x, r), sa)
    oa = setdiff(1:k, ia)

    #= 
    likelihoods
    =#

    # no events
    nev = noevents_expr(ri, lr, ia, oa, k, ns, false)

    # local extinction
    # remove if !isone(lr)
    lex = localext_expr(r, ia, sa, S, k)

    # dispersal
    # remove if lr == k
    dis = dispersal_expr(r, lr, ia, oa, S, ns, k, false)

    # within-region speciation
    wrs = wrspec_expr(ri, ia, ns)

    # between-region speciation
    # remove if !isone(lr)
    brs = brspec_expr(r, S, ns, k)

    # if single area
    if ri <= k
      push!(eqs.args, 
        :(du[$ri] = $nev + $dis + $wrs))
    # if widespread
    elseif ri != ns
      push!(eqs.args, 
        :(du[$ri] = $nev + $lex + $dis + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$ri] = $nev + $lex + $wrs + $brs))
    end

    #= 
    extinctions
    =#

    # no events
    nev = noevents_expr(ri, lr, ia, oa, k, ns, true)

    # extinction
    ext = ext_expr(r, ia, sa, S, ns, k) 

    # dispersal and extinction
    dis = dispersal_expr(r, lr, ia, oa, S, ns, k, true)

    # within-region extinction
    wrs = wrsext_expr(ri, ia, ns)

    # between-region extinction
    brs = brsext_expr(r, S, ns, k)

    # if single area
    if ri <= k
      push!(eqs.args, 
        :(du[$(ri + ns)] = $nev + $ext + $dis + $wrs))
    # if widespread
    elseif ri != ns
      push!(eqs.args, 
        :(du[$(ri + ns)] = $nev + $ext + $dis + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$(ri + ns)] = $nev + $ext + $wrs + $brs))
    end

  end
  
  ## aesthetic touches
  # remove REPL comment
  popfirst!(eqs.args)

  eqs.args[:] = eqs.args[append!([1:2:end...],[2:2:end...])]

  return eqs
end





"""
    noevents_expr(ri::Int64,
                  lr::Int64,
                  ia::Array{Int64,1},
                  oa::Array{Int64,1},
                  k ::Int64)

Return expression for no events.
"""
function noevents_expr(ri ::Int64,
                       lr ::Int64,
                       ia ::Array{Int64,1},
                       oa ::Array{Int64,1},
                       k  ::Int64,
                       ns ::Int64,
                       ext::Bool)

  wu = ext ? ns : 0

  ts = isone(lr) ? 0 : (k*(k-1) + k)

  ex = :(+ ($(2^(lr-1) - 1.) * p[$(k+1)]))
  for (i, a) = enumerate(ia)
    push!(ex.args, :(p[$a] + p[$(a + k + ts + 1)]))
    for j = oa
      j -= a <= j ? 1 : 0
      push!(ex.args[i+2].args, :(p[$(2k + 1 + (k-1)*(a-1) + j)]))
    end
  end
  ex = :(-1.0 * $ex * u[$(ri + wu)])

  # remove 0 product if single area
  if isone(lr)
    ex.args[3] = ex.args[3].args[3]
  end

  return ex
end





"""
    localext_expr(r ::String,
                  ia::Array{Int64,1},
                  sa::Array{String,1},
                  S ::Array{String,1},
                  k ::Int64)

Return expression for local extinction.
"""
function localext_expr(r ::String,
                       ia::Array{Int64,1},
                       sa::Array{String,1},
                       S ::Array{String,1},
                       k ::Int64)

  ex = Expr(:call, :+)
  for a = ia
    push!(ex.args, :(p[$(k^2 + k + 1 + a)] * 
                     u[$(findfirst(x -> x == replace(r, sa[a] => ""), S))]))
  end

  return ex
end





"""
    dispersal_expr(r ::String,
                   lr::Int64,
                   ia::Array{Int64,1},
                   oa::Array{Int64,1},
                   S ::Array{String,1},
                   k ::Int64)

Return expression for dispersal.
"""
function dispersal_expr(r  ::String,
                        lr ::Int64,
                        ia ::Array{Int64,1},
                        oa ::Array{Int64,1},
                        S  ::Array{String,1},
                        ns ::Int64,
                        k  ::Int64,
                        ext::Bool)

  wu = ext ? ns : 0

  ida = findall(x -> all(occursin.(split(r,""),x)) && 
                     lastindex(x) == (lr + 1), 
                S)

  ex = Expr(:call, :+)
  for a = ia, (i, j) = enumerate(oa)
    j -= a <= j ? 1 : 0
    push!(ex.args, :(p[$(2k + 1 + (k-1)*(a-1) + j)] * u[$(ida[i] + wu)]))
  end

  return ex
end





"""
    wrspec_expr(ri::Int64,
                ia::Array{Int64,1},
                ns::Int64)

Return expression for within-region speciation.
"""
function wrspec_expr(ri ::Int64,
                     ia ::Array{Int64,1},
                     ns ::Int64)

  if isone(lastindex(ia)) 
    ex = :(2.0 * p[$ri] * u[$(ri + ns)] * u[$(ri)])
  else

    ex = Expr(:call, :+)
    for i = ia
      push!(ex.args, :(p[$i] * (u[$(i + ns)] * u[$(ri)] + u[$(ri + ns)] * u[$(i)])))
    end
  end

  return ex
end





"""
    brspec_expr(r::String,
                S ::Array{String,1},
                ns::Int64)

Return expression for between-region speciation.
"""
function brspec_expr(r  ::String,
                     S  ::Array{String,1},
                     ns ::Int64,
                     k  ::Int64)

  va  = vicsubsets(r)
  ex = Expr(:call, :+)
  for (la, ra) = va
    push!(ex.args,
      :(u[$(findfirst(isequal(ra), S) + ns)] *
        u[$(findfirst(isequal(la), S))]))
  end
  ex = :($(2^(lastindex(r)-1) - 1.0) * p[$(k+1)] * $ex)

  isone(ex.args[2]) && deleteat!(ex.args, 2)

  return ex
end





"""
    ext_expr(r ::String,
                  ia::Array{Int64,1},
                  sa::Array{String,1},
                  S ::Array{String,1},
                  k ::Int64)

Return expression for extinction.
"""
function ext_expr(r ::String,
                  ia::Array{Int64,1},
                  sa::Array{String,1},
                  S ::Array{String,1},
                  ns::Int64,
                  k ::Int64)

  if isone(lastindex(ia))
    ex = :(p[$(k+1 + ia[1])])
  else
    ex = Expr(:call, :+)
    for a = ia
      push!(ex.args, :(p[$(k^2 + k + 1 + a)] * 
                       u[$(findfirst(x -> x == replace(r, sa[a] => ""), S) + ns)]))
    end
  end

  return ex
end




"""
    wrsext_expr(ri::Int64,
                ia::Array{Int64,1},
                ns::Int64)

Return expression of extinction for within-region speciation.
"""
function wrsext_expr(ri ::Int64,
                     ia ::Array{Int64,1},
                     ns ::Int64)

  if isone(lastindex(ia)) 
    ex = :(p[$ri] * u[$(ri + ns)]^2)
  else
    ex = Expr(:call, :+)
    for i = ia
      push!(ex.args, :(p[$i] * u[$(i + ns)] * u[$(ri + ns)]))
    end
  end

  return ex
end






"""
    brsext_expr(r::String,
                S ::Array{String,1},
                ns::Int64)

Return expression for extinction for between-region speciation.
"""
function brsext_expr(r  ::String,
                     S  ::Array{String,1},
                     ns ::Int64,
                     k  ::Int64)

  va  = vicsubsets(r)[1:div(end,2)]
  ex = Expr(:call, :+)
  for (la, ra) = va
    push!(ex.args,
      :(u[$(findfirst(isequal(ra), S) + ns)] *
        u[$(findfirst(isequal(la), S) + ns)]))
  end
  ex = :($(2^(lastindex(r)-1) - 1.0) * p[$(k+1)] * $ex)

  isone(ex.args[2]) && deleteat!(ex.args, 2)

  return ex
end






#=
    1   2   3    4   5   6   7
p = sa, sb, sab, xa, xb, qa, qb

    1   2   3    4   5   6
u = da, db, dab, ea, eb, eab
=#
"""
    geosse_2k(du::Array{Float64,1}, 
              u::Array{Float64,1}, 
              p::Array{Float64,1}, 
              t::Float64)

GeoSSE ODE equation for 2 areas.
"""
function geosse_2k(du::Array{Float64,1}, 
                   u::Array{Float64,1}, 
                   p::Array{Float64,1}, 
                   t::Float64)

  @inbounds begin
    # probabilities
    du[1] = -(p[1] + p[4] + p[6])*u[1] + p[6]*u[3] + 2.0*p[1]*u[1]*u[4]
    du[2] = -(p[2] + p[5] + p[7])*u[2] + p[7]*u[3] + 2.0*p[2]*u[2]*u[5]
    du[3] = -(p[1] + p[2] + p[3] + p[4] + p[5])*u[3] + 
              p[4]*u[2] + p[5]*u[1]        +
              p[1]*(u[4]*u[3] + u[6]*u[1]) + 
              p[2]*(u[5]*u[3] + u[6]*u[2]) + 
              p[3]*(u[4]*u[2] + u[5]*u[1])

    # extinction
    du[4] = -(p[1] + p[4] + p[6])*u[4] + p[4] + p[6]*u[6] + p[1]*u[4]^2
    du[5] = -(p[2] + p[5] + p[7])*u[5] + p[5] + p[7]*u[6] + p[2]*u[5]^2
    du[6] = -(p[1] + p[2] + p[3] + p[4] + p[5])*u[6] + 
              p[4]*u[5] + p[5]*u[4] + p[1]*u[6]*u[4] + 
              p[2]*u[6]*u[5] + p[3]*u[4]*u[5]
  end
  return nothing
end




