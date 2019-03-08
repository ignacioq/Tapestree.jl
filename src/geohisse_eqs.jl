#=

GeoHiSSE equations

Ignacio Quintero Mächler

t(-_-t)


Created 07 03 2019
=#





"""
    struct ghs
      g::Set{Int64}
      h::Int64
    end

Composite type for Geographical `g` & Hidden `h` states
"""
struct ghs
  g::Set{Int64}
  h::Int64
end


"""
    isequal(x::ghs, y::ghs)

Compares equality between two `ghs` types.
"""
isghsequal(x::ghs, y::ghs) = x.g == y.g && x.h == y.h 




k = 2
h = 2


sort(collect(build_par_names(k,h,(true,false,false))), by = x -> x[2])



"""
    make_geohisse(k::Int64, h::Int64)

GeoSSE ODE equation for k areas and `h` hidden states.
"""
function make_geohisse(k::Int64, h::Int64)

  # n states
  ns = (2^k - 1)*h

  # create individual areas subsets
  gs = Array{Set{Int64}, 1}(undef, k)
  for i = Base.OneTo(k)
    gs[i] = Set([i])
  end
  gS = sets(gs)
  popfirst!(gS)            # remove empty
  sort!(gS, by = length)   # arrange by range size

  # add hidden states and create ghs objects
  S = Array{ghs, 1}()
  for i = 0:(h-1), j = gS
    push!(S, ghs(j, i))
  end

  # Metaprogram to generate equations
  eqs = quote end

  for si = Base.OneTo(ns)

    # state range
    s = S[si]

    # length of geographic range
    ls = length(s.g)

    # which single areas do not occur in r
    oa = setdiff(1:k, s.g)

    #= 
    likelihoods
    =#

    # no events
    nev = noevents_expr(si, s, ls, oa, k, h, ns, false)

    # local extinction
    # remove if !isone(lr)
    lex = localext_expr(s, S, k, h)

    # dispersal
    # remove if lr == k
    dispersal_expr(s, oa, S, ns, k, h, false)

    # within-region speciation
    wrs = wrspec_expr(si, s, ns, k)








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
function noevents_expr(si ::Int64,
                       s  ::ghs, 
                       ls ::Int64,
                       oa ::Array{Int64,1},
                       k  ::Int64,
                       h  ::Int64,
                       ns ::Int64,
                       ext::Bool)
  # if extinction
  wu = ext ? ns : 0

  ts = isone(ls) ? 0 : k*h*k

  # between-region speciation
  ex = :(+ ($(2^(ls-1) - 1.) * p[$(k+1+s.h*(k+1))]))

  for (i, v) = enumerate(s.g)
    # speciation and extinction
    push!(ex.args, :(p[$(v + s.h*(k+1))] + p[$(v + (k+1)*h + s.h*k + ts)]))
    # dispersal
    for j = oa
      j -= v <= j ? 1 : 0
      push!(ex.args[i+2].args, :(p[$(2k*h + h + (k-1)*(v-1) + j + s.h*k)]))
    end
  end

  # add hidden state shifts
  push!(ex.args[3].args, :(p[$(s.h + 1 + h*(1+k)^2)]))

  # multiply by u
  ex = :(-1.0 * $ex * u[$(si + wu)])

  # remove 0 product if single area
  if isone(ls)
    ex.args[3] = ex.args[3].args[3]
  end

  return ex
end





"""
    localext_expr(s ::ghs,
                  S ::Array{ghs,1},
                  k ::Int64,
                  h ::Int64)

Return expression for local extinction.
"""
function localext_expr(s ::ghs,
                       S ::Array{ghs,1},
                       k ::Int64,
                       h ::Int64)

  ex = Expr(:call, :+)
  for i = s.g
    push!(ex.args, :(p[$(i + (k+1)*h + s.h*k + k*h*k)] * 
                     u[$(findfirst(x -> isghsequal(x, ghs(setdiff(s.g, i),s.h)), S))]))
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
function dispersal_expr(s  ::ghs,
                        oa ::Array{Int64,1},
                        S  ::Array{ghs,1},
                        ns ::Int64,
                        k  ::Int64,
                        h  ::Int64,
                        ext::Bool)

  wu = ext ? ns : 0

  ida = findall(x -> length(union(s.g, x.g)) == 2 && 
                     length(intersect(s.g, x.g)) == 1, S)

  ex = Expr(:call, :+)
  for a = s.g, (i, j) = enumerate(oa)
    j -= a <= j ? 1 : 0
    push!(ex.args, :(p[$(2k*h + h + (k-1)*(a-1) + j + s.h*k)] * u[$(ida[i] + wu)]))
  end

  return ex
end





"""
    wrspec_expr(ri::Int64,
                ia::Array{Int64,1},
                ns::Int64)

Return expression for within-region speciation.
"""
function wrspec_expr(si ::Int64,
                     s  ::ghs,
                     ns ::Int64,
                     k  ::Int64)

  if isone(length(s.g)) 
    ex = :(2.0 * p[$si] * u[$(si + ns)] * u[$(si)])
  else

    ex = Expr(:call, :+)
    for i = s.g
      push!(ex.args, :(p[$(i + s.h*(k+1))] * 
        (u[$(i + s.h*(2^k-1) + ns)] * u[$(si)] + 
         u[$(si + ns)] * u[$(i + s.h*(2^k-1))])))
    end
  end

  return ex
end





"""
    brspec_expr(s ::ghs,
                S ::Array{ghs,1},
                ns::Int64,
                k ::Int64)

Return expression for between-region speciation.
"""
function brspec_expr(s ::ghs,
                     S ::Array{ghs,1},
                     ns::Int64,
                     k ::Int64)

  va  = vicsubsets(s.g)
  ex = Expr(:call, :+)
  for (la, ra) = va
    push!(ex.args,
      :(u[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h + ns)] *
        u[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)]))
  end
  ex = :($(2^(length(s.g)-1) - 1.0) * p[$(k+1 + (2^k-1)*s.h)] * $ex)

  isone(ex.args[2]) && deleteat!(ex.args, 2)

  return ex
end

 brspec_expr(s, S, ns, k)

#=
check for different hidden states
=#





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





"""
    egeohisse_2k_s(du::Array{Float64,1}, 
                   u::Array{Float64,1}, 
                   p::Array{Float64,1}, 
                   t::Float64)

Speciation EGeoHiSSE equation for 2 areas and 2 hidden states.
"""
function make_egeohisse_2k_s(af::Function)

  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)

    # UNIdimensional time function
    aft = af(t)

    sa1 = p[1]  * exp(aft*p[21])
    sb1 = p[2]  * exp(aft*p[22])
    sa2 = p[10] * exp(aft*p[23])
    sb2 = p[11] * exp(aft*p[24])

    @inbounds begin

      ### Hidden State 1
      ## probabilities 
      # area A
      du[1] = -(sa1 + p[4] + p[6] + p[19])*u[1]            + 
                p[6]*u[3] + p[19]*u[7] + 2.0*sa1*u[1]*u[4]
      # area B
      du[2] = -(sb1 + p[5] + p[7] + p[19])*u[2]            + 
                p[7]*u[3] + p[19]*u[8] + 2.0*sb1*u[2]*u[5]
      # area AB
      du[3] = -(sa1 + sb1 + p[3] + p[8] + p[9] + p[19])*u[3] + 
                p[8]*u[1] + p[9]*u[2]+ p[19]*u[9]            +
                sa1*(u[4]*u[3] + u[6]*u[1])                  + 
                sb1*(u[5]*u[3] + u[6]*u[2])                  + 
                p[3]*(u[4]*u[2] + u[5]*u[1])
      ## extinction
      # area A
      du[4] = -(sa1 + p[4] + p[6] + p[19])*u[4]             + 
                p[4] + p[19]*u[10] + p[6]*u[6] + sa1*u[4]^2
      # area B
      du[5] = -(sb1 + p[5] + p[7] + p[19])*u[5]             + 
                p[5] + p[19]*u[11] + p[7]*u[6] + sb1*u[5]^2
      # area AB
      du[6] = -(sa1 + sb1 + p[3] + p[8] + p[9] + p[19])*u[6]   +
                p[8]*u[4] + p[9]*u[5] + p[19]*u[12]              + 
                sa1*u[6]*u[4] + sb1*u[6]*u[5] + p[3]*u[4]*u[5]
      ### Hidden State 2
      ## probabilities 
      # area A
      du[7] = -(sa2 + p[13] + p[15] + p[20])*u[7]            + 
                p[15]*u[9] + p[20]*u[1] + 2.0*sa2*u[7]*u[10]
      # area B
      du[8] = -(sb2 + p[14] + p[16] + p[20])*u[8]            + 
                p[16]*u[9] + p[20]*u[2] + 2.0*sb2*u[8]*u[11]
      # area AB
      du[9] = -(sa2 + sb2 + p[12] + p[17] + p[18] + p[20])*u[9] + 
                p[17]*u[7] + p[18]*u[8]+ p[20]*u[3]                 +
                sa2*(u[10]*u[9] + u[12]*u[7])                     + 
                sb2*(u[11]*u[9] + u[12]*u[8])                     + 
                p[12]*(u[10]*u[8] + u[11]*u[7])
      ## extinction
      # area A
      du[10] = -(sa2 + p[13] + p[15] + p[20])*u[10]            + 
                 p[13] + p[20]*u[4] + p[15]*u[12] + sa2*u[10]^2
      # area B
      du[11] = -(sb2 + p[14] + p[16] + p[20])*u[11]            + 
                 p[14] + p[20]*u[5] + p[16]*u[12] + sb2*u[11]^2
      # area AB
      du[12] = -(sa2 + sb2 + p[12] + p[17] + p[18] + p[20])*u[6]       +
                 p[17]*u[10] + p[18]*u[11] + p[20]*u[6]                    + 
                 sa2*u[10]*u[12] + sb2*u[11]*u[12] + p[12]*u[10]*u[11]
    end

    return nothing
  end

  return f
end





"""
    geohisse_2k(du::Array{Float64,1}, 
                u::Array{Float64,1}, 
                p::Array{Float64,1}, 
                t::Float64)

GeoHiSSE + extinction most general ODE equation for 2 areas & 2 hidden states.
"""
function geohisse_2k(du::Array{Float64,1}, 
                     u::Array{Float64,1}, 
                     p::Array{Float64,1}, 
                     t::Float64)

  @inbounds begin

    ### Hidden State 1
    ## probabilities 
    # area A
    du[1] = -(p[1] + p[4] + p[6] + p[19])*u[1]            + 
              p[6]*u[3] + p[19]*u[7] + 2.0*p[1]*u[1]*u[4]
    # area B
    du[2] = -(p[2] + p[5] + p[7] + p[19])*u[2]            + 
              p[7]*u[3] + p[19]*u[8] + 2.0*p[2]*u[2]*u[5]
    # area AB
    du[3] = -(p[1] + p[2] + p[3] + p[8] + p[9] + p[19])*u[3] + 
              p[8]*u[1] + p[9]*u[2]+ p[19]*u[9]              +
              p[1]*(u[4]*u[3] + u[6]*u[1])                   + 
              p[2]*(u[5]*u[3] + u[6]*u[2])                   + 
              p[3]*(u[4]*u[2] + u[5]*u[1])
    ## extinction
    # area A
    du[4] = -(p[1] + p[4] + p[6] + p[19])*u[4]             + 
              p[4] + p[19]*u[10] + p[6]*u[6] + p[1]*u[4]^2
    # area B
    du[5] = -(p[2] + p[5] + p[7] + p[19])*u[5]             + 
              p[5] + p[19]*u[11] + p[7]*u[6] + p[2]*u[5]^2
    # area AB
    du[6] = -(p[1] + p[2] + p[3] + p[8] + p[9] + p[19])*u[6]   +
              p[8]*u[4] + p[9]*u[5] + p[19]*u[12]              + 
              p[1]*u[6]*u[4] + p[2]*u[6]*u[5] + p[3]*u[4]*u[5]
    ### Hidden State 2
    ## probabilities 
    # area A
    du[7] = -(p[10] + p[13] + p[15] + p[20])*u[7]            + 
              p[15]*u[9] + p[20]*u[1] + 2.0*p[10]*u[7]*u[10]
    # area B
    du[8] = -(p[11] + p[14] + p[16] + p[20])*u[8]            + 
              p[16]*u[9] + p[20]*u[2] + 2.0*p[11]*u[8]*u[11]
    # area AB
    du[9] = -(p[10] + p[11] + p[12] + p[17] + p[18] + p[20])*u[9] + 
              p[17]*u[7] + p[18]*u[8]+ p[20]*u[3]                 +
              p[10]*(u[10]*u[9] + u[12]*u[7])                     + 
              p[11]*(u[11]*u[9] + u[12]*u[8])                     + 
              p[12]*(u[10]*u[8] + u[11]*u[7])
    ## extinction
    # area A
    du[10] = -(p[10] + p[13] + p[15] + p[20])*u[10]            + 
               p[13] + p[20]*u[4] + p[15]*u[12] + p[10]*u[10]^2
    # area B
    du[11] = -(p[11] + p[14] + p[16] + p[20])*u[11]            + 
               p[14] + p[20]*u[5] + p[16]*u[12] + p[11]*u[11]^2
    # area AB
    du[12] = -(p[10] + p[11] + p[12] + p[17] + p[18] + p[20])*u[6]       +
               p[17]*u[10] + p[18]*u[11] + p[20]*u[6]                    + 
               p[10]*u[10]*u[12] + p[11]*u[11]*u[12] + p[12]*u[10]*u[11]
  end

  return nothing
end




