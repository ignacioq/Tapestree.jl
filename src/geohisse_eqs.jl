#=

GeoHiSSE equations

Ignacio Quintero Mächler

t(-_-t)


Created 07 03 2019
=#





"""
    make_geohisse(k::Int64, h::Int64, name::Symbol)

makes GeoSSE ODE equation function for `k` areas 
and `h` hidden states of name `:name`.
"""
function make_geohisse(k::Int64, h::Int64, name::Symbol)

  # n states
  ns = (2^k - 1)*h

  # create individual areas subsets
  S = create_states(k, h)

  # start Expression for ODE equations
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

    # remove if !isone(lr)
    if length(s.g) > 1
      # local extinction
      lex = localext_expr(s, S, k, h)
      # between-region speciation
      brs = brspec_expr(s, S, ns, k, false)
    end

    # dispersal
    # remove if lr == k
    if length(s.g) != k
      dis = dispersal_expr(s, oa, S, ns, k, h, false)
    end

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, false) : :0.0

    # within-region speciation
    wrs = wrspec_expr(si, s, ns, k)

    # push `D` equation to to eqs
    if isone(length(s.g))
      push!(eqs.args, 
        :(du[$si] = $nev + $dis + $hid + $wrs))
    elseif length(s.g) == k
      push!(eqs.args, 
        :(du[$si] = $nev + $lex + $hid + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$si] = $nev + $lex + $dis + $hid + $wrs + $brs))
    end

    #= 
    extinctions
    =#

    # no events
    nev = noevents_expr(si, s, ls, oa, k, h, ns, true)

    # extinction
    ext = ext_expr(s, S, ns, k, h) 

    # dispersal and extinction
    if length(s.g) != k
      dis = dispersal_expr(s, oa, S, ns, k, h, true)
    end

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, true) : :0.0

    # within-region extinction
    wrs = wrsext_expr(si, s, ns, k)

    # between-region extinction
    if length(s.g) > 1
      brs = brspec_expr(s, S, ns, k, true)
    end

    # push `E` equation to to eqs
    if isone(length(s.g))
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $dis + $hid + $wrs))
    elseif length(s.g) == k
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $hid + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $dis + $hid + $wrs + $brs))
    end

  end

  ## aesthetic touches
  # remove REPL comment
  popfirst!(eqs.args)

  eqs.args[:] = eqs.args[append!([1:2:end...],[2:2:end...])]

  # make function
  ex = 
    quote 
      function $name(du::Array{Float64,1}, u::Array{Float64,1}, p::Array{Float64,1}, t::Float64)
        @inbounds begin
          $eqs
        end
        return nothing
      end
    end

  return eval(ex)
end





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



"""
    create_states(k::Int64, h::Int64)

Create GeoHiSSE states
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

  # add hidden states and create ghs objects
  S = Array{ghs, 1}()
  for i = 0:(h-1), j = gS
    push!(S, ghs(j, i))
  end

  return S
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
      push!(ex.args[i+2].args, 
        :(p[$(2h*k + h + s.h*k*(k-1) + (k-1)*(v-1) + j)]))
    end
  end

  # add hidden state shifts
  for hi in setdiff(0:(h-1), s.h)
    hi -= s.h <= hi ? 0 : -1
    push!(ex.args[3].args, :(p[$(h*(3k + 1 + k*(k-1)) + s.h*(h-1) + hi)]))
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

  ida = findall(x -> length(union(s.g, x.g))     == length(s.g)+1 && 
                     length(intersect(s.g, x.g)) == length(s.g)   &&
                     s.h == x.h, S)

  ex = Expr(:call, :+)
  for a = s.g, (i, j) = enumerate(oa)
    j -= a <= j ? 1 : 0
    push!(ex.args, :(p[$(2h*k + h + s.h*k*(k-1) + (k-1)*(a-1) + j)] * 
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
    hi -= s.h <= hi ? 0 : -1
    push!(ex.args, :(p[$(h*(3k + 1 + k*(k-1)) + s.h*(h-1) + hi)] * 
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
  ex = :($(2^(length(s.g)-1) - 1.0) * p[$(k+1+s.h*(k+1))] * $ex)

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

GeoHiSSE + extinction general ODE equation for 2 areas & 2 hidden states.
"""
function geohisse_2k(du::Array{Float64,1}, 
                     u::Array{Float64,1}, 
                     p::Array{Float64,1}, 
                     t::Float64)

  @inbounds begin

    ### likelihoods 
    ## hidden state 1
    # area A
    du[1] = -(p[1] + p[7] + p[11] + p[19])*u[1]            + 
              p[11]*u[3] + p[19]*u[4] + 2.0*p[1]*u[1]*u[7]
    # area B
    du[2] = -(p[2] + p[8] + p[12] + p[19])*u[2]            + 
              p[12]*u[3] + p[19]*u[5] + 2.0*p[2]*u[2]*u[8]
    # area AB
    du[3] = -(p[1] + p[2] + p[3] + p[15] + p[16] + p[19])*u[3] + 
              p[16]*u[1] + p[15]*u[2]+ p[19]*u[6]              + 
              p[1]*(u[7]*u[3] + u[9]*u[1])                     + 
              p[2]*(u[8]*u[3] + u[9]*u[2])                     + 
              p[3]*(u[7]*u[2] + u[8]*u[1])
    ## hidden state 2
    # area A
    du[4] = -(p[4] + p[9] + p[13] + p[20])*u[4]            + 
              p[13]*u[6] + p[20]*u[1] + 2.0*p[4]*u[4]*u[10]
    # area B
    du[5] = -(p[5] + p[10] + p[14] + p[20])*u[5]            + 
              p[14]*u[6] + p[20]*u[2] + 2.0*p[5]*u[5]*u[11]
    # area AB
    du[6] = -(p[4] + p[5] + p[6] + p[17] + p[18] + p[20])*u[6] + 
              p[17]*u[5] + p[18]*u[4]+ p[20]*u[3]              + 
              p[4]*(u[10]*u[6] + u[12]*u[4])                   + 
              p[5]*(u[11]*u[6] + u[12]*u[5])                   + 
              p[6]*(u[10]*u[5] + u[11]*u[4])

    ### extinction
    ## hidden state 1
    # area A
    du[7] = -(p[1] + p[7] + p[11] + p[19])*u[7]             + 
              p[7] + p[19]*u[10] + p[11]*u[9] + p[1]*u[7]^2
    # area B
    du[8] = -(p[2] + p[8] + p[12] + p[19])*u[8]             + 
              p[8] + p[19]*u[11] + p[12]*u[9] + p[2]*u[8]^2
    # area AB
    du[9] = -(p[1] + p[2] + p[3] + p[15] + p[16] + p[19])*u[9]   +
              p[15]*u[8] + p[16]*u[7] + p[19]*u[12]              + 
              p[1]*u[9]*u[7] + p[2]*u[9]*u[8] + p[3]*u[7]*u[8]
    ## hidden state 2
    # area A
    du[10] = -(p[4] + p[9] + p[13] + p[20])*u[10]            + 
               p[9] + p[20]*u[7] + p[13]*u[12] + p[4]*u[10]^2
    # area B
    du[11] = -(p[5] + p[10] + p[14] + p[20])*u[11]            + 
               p[10] + p[20]*u[8] + p[14]*u[12] + p[5]*u[11]^2
    # area AB
    du[12] = -(p[4] + p[5] + p[6] + p[17] + p[18] + p[20])*u[12]       +
               p[17]*u[11] + p[18]*u[10] + p[20]*u[9]                    + 
               p[4]*u[10]*u[12] + p[5]*u[11]*u[12] + p[6]*u[10]*u[11]
  end

  return nothing
end




