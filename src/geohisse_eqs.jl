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




#sort(collect(build_par_names(k,h,(true,false,false))), by = x -> x[2])



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
    dis = dispersal_expr(s, oa, S, ns, k, h, false)

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, false) : :0

    # within-region speciation
    wrs = wrspec_expr(si, s, ns, k)

    # between-region speciation
    # remove if !isone(lr)
    brs = brspec_expr(s, S, ns, k, false)

    # if single area
    if isone(length(s.g))
      push!(eqs.args, 
        :(du[$si] = $nev + $dis + $hid + $wrs))
    elseif length(s.g) != k
      push!(eqs.args, 
        :(du[$si] = $nev + $lex + $dis + $hid + $wrs + $brs))
    # if widespread
    else
      push!(eqs.args, 
        :(du[$si] = $nev + $lex + $hid + $wrs + $brs))
    end

    #= 
    extinctions
    =#

    # no events
    nev = noevents_expr(si, s, ls, oa, k, h, ns, true)

    # extinction
    ext = ext_expr(s, S, ns, k, h) 

    # dispersal and extinction
    dis = dispersal_expr(s, oa, S, ns, k, h, true)

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, true) : :0

    # within-region extinction
    wrs = wrsext_expr(si, s, ns, k)

    # between-region extinction
    brs = brspec_expr(s, S, ns, k, true)

    # if single area
    if isone(length(s.g))
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $dis + $hid + $wrs))
    elseif length(s.g) != k
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $dis + $hid + $wrs + $brs))
    # if widespread
    else
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $hid + $wrs + $brs))
    end

  end
  
  ## aesthetic touches
  # remove REPL comment
  popfirst!(eqs.args)

  eqs.args[:] = eqs.args[append!([1:2:end...],[2:2:end...])]

  return eqs
end





"""
    noevents_expr(si ::Int64,
                  s  ::ghs, 
                  ls ::Int64,
                  oa ::Array{Int64,1},
                  k  ::Int64,
                  h  ::Int64,
                  ns ::Int64,
                  ext::Bool)

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
  ex = :(+ ($(2.0^(ls-1) - 1.) * p[$(k+1+s.h*(k+1))]))

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
  for hi in setdiff(0:(h-1), s.h)
    hi -= s.h <= hi ? 1 : 0
    push!(ex.args[3].args, :(p[$(s.h + 1 + hi + h*(k+1)^2)]))
  end

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
    dispersal_expr(s  ::ghs,
                   oa ::Array{Int64,1},
                   S  ::Array{ghs,1},
                   ns ::Int64,
                   k  ::Int64,
                   h  ::Int64,
                   ext::Bool)

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

  ida = findall(x -> length(union(s.g, x.g))     == 2 && 
                     length(intersect(s.g, x.g)) == 1 &&
                     s.h == x.h, S)

  ex = Expr(:call, :+)
  for a = s.g, (i, j) = enumerate(oa)
    j -= a <= j ? 1 : 0
    push!(ex.args, :(p[$(2k*h + h + (k-1)*(a-1) + j + s.h*k)] * 
                     u[$(ida[i] + wu)]))
  end

  return ex
end





"""
    hidtran_expr(s  ::ghs,
                 S  ::Array{ghs,1},
                 ns ::Int64,
                 k  ::Int64,
                 h  ::Int64,
                 ext::Bool)

Return expression for hidden states transitions.
"""
function hidtran_expr(s  ::ghs,
                      S  ::Array{ghs,1},
                      ns ::Int64,
                      k  ::Int64,
                      h  ::Int64,
                      ext::Bool)

  wu = ext ? ns : 0

  hs = findall(x -> isequal(s.g, x.g) && 
                   !isequal(s.h, x.h), S)

  ex = Expr(:call, :+)
  for i in hs
    hi = S[i].h
    hi -= s.h <= hi ? 1 : 0
    push!(ex.args, :(p[$(s.h + 1 + hi + h*(k+1)^2)] * 
                     u[$(i + wu)]))
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
    brspec_expr(s  ::ghs,
                S  ::Array{ghs,1},
                ns ::Int64,
                k  ::Int64, 
                ext::Bool)

Return expression for between-region speciation.
"""
function brspec_expr(s  ::ghs,
                     S  ::Array{ghs,1},
                     ns ::Int64,
                     k  ::Int64, 
                     ext::Bool)

  if ext
    va = vicsubsets(s.g)[1:div(end,2)]
    wu = ns
  else
    va = vicsubsets(s.g)
    wu = 0
  end

  ex = Expr(:call, :+)
  for (la, ra) = va
    push!(ex.args,
      :(u[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h + ns)] *
        u[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h + wu)]))
  end
  ex = :($(2^(length(s.g)-1) - 1.0) * p[$(k+1 + (2^k-1)*s.h)] * $ex)

  isone(ex.args[2]) && deleteat!(ex.args, 2)

  return ex
end





"""
    ext_expr(s ::ghs,
             S ::Array{ghs,1},
             ns::Int64,
             k ::Int64,
             h ::Int64)

Return expression for extinction.
"""
function ext_expr(s ::ghs,
                  S ::Array{ghs,1},
                  ns::Int64,
                  k ::Int64,
                  h ::Int64)

  if isone(length(s.g))
    ex = Expr(:call, :+)
    for i in s.g push!(ex.args, :(p[$((k+1)*h + h*(s.h) + i)])) end
  else
    ex = Expr(:call, :+)
    for i = s.g
      push!(ex.args, :(p[$(i + (k+1)*h + s.h*k + k*h*k)] * 
                       u[$(findfirst(x -> isghsequal(x, ghs(setdiff(s.g, i),s.h)), S) + ns)]))
    end
  end

  return ex
end





"""
    wrsext_expr(si::Int64,
                s ::ghs,
                ns::Int64,
                k ::Int64)

Return expression of extinction for within-region speciation.
"""
function wrsext_expr(si::Int64,
                     s ::ghs,
                     ns::Int64,
                     k ::Int64)

  if isone(length(s.g)) 
    ex = :(p[$si] * u[$(si + ns)]^2)
  else
    ex = Expr(:call, :+)
    for i = s.g
      push!(ex.args, :(p[$(i + s.h*(k+1))] * 
                       u[$(i + s.h*(2^k-1) + ns)] * u[$(si + ns)]))
    end
  end

  return ex
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



