#=

EGeoHiSSE equations

Ignacio Quintero Mächler

t(-_-t)

Created 18 03 2019
=#






"""
    make_egeohisse(k    ::Int64,
                   h    ::Int64,
                   ny   ::Int64,
                   af!  ::Function,
                   model::NTuple{3, Bool},
                   name ::Symbol)

Creates Covariate GeoHiSSE ODE equation function for `k` areas, 
`h` hidden states and `ny` covariates of name `:name`.
"""
function make_egeohisse(k    ::Int64,
                        h    ::Int64,
                        ny   ::Int64,
                        af!  ::Function,
                        model::NTuple{3, Bool})
  # n states
  ns = (2^k - 1)*h

  # create individual areas subsets
  S = create_states(k, h)

  # start Expression for ODE equations
  eqs = quote end
  popfirst!(eqs.args)

  # add environmental function
  push!(eqs.args, :(Base.invokelatest(af!,t, r)))

  # compute exponential for each β*covariable
  expex = exp_expr(k, h, ny, model)

  for ex in expex.args
    push!(eqs.args, ex)
  end

  """
  - Add closure for eaft vector!!
  """

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
    nev = noevents_expr(si, s, ls, oa, k, h, ns, model, false)

    # remove if !isone(lr)
    if ls > 1
      # local extinction
      lex = localext_expr(s, S, k, h)
      # between-region speciation
      brs = brspec_expr(s, S, ns, k, false)
    end

    # dispersal
    # remove if lr == k
    if ls != k
      dis = dispersal_expr(s, oa, S, ns, k, h, model[3], false)
    end

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, false) : :0.0

    # within-region speciation
    wrs = wrspec_expr(si, s, ns, k, model[1])

    # push `D` equation to to eqs
    if isone(ls)
      push!(eqs.args, 
        :(du[$si] = $nev + $dis + $hid + $wrs))
    elseif ls == k
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
    nev = noevents_expr(si, s, ls, oa, k, h, ns, model, true)

    # extinction
    ext = ext_expr(s, S, ns, k, h, model[2]) 

    # dispersal and extinction
    if ls != k
      dis = dispersal_expr(s, oa, S, ns, k, h, model[3], true)
    end

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, true) : :0.0

    # within-region extinction
    wrs = wrsext_expr(si, s, ns, k, model[1])

    # between-region extinction
    if ls > 1
      brs = brspec_expr(s, S, ns, k, true)
    end

    # push `E` equation to to eqs
    if isone(ls)
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $dis + $hid + $wrs))
    elseif ls == k
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $hid + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$(si + ns)] = $nev + $ext + $dis + $hid + $wrs + $brs))
    end

  end

  ## aesthetic touches
  # sort `du`s
  eqs.args[[(end-(ns*2)+1):end...]] = 
    eqs.args[append!([(end-(ns*2)+1):2:end...],
            [(end-(ns*2)+2):2:end...])]

  ex = quote 
    function ode_full(du   ::Array{Float64,1}, 
                      u    ::Array{Float64,1}, 
                      p    ::Array{Float64,1}, 
                      t    ::Float64,
                      r    ::Array{Float64,1},
                      eaft ::Array{Float64,1},
                      k    ::Int64,
                      h    ::Int64,
                      ny   ::Int64,
                      af!  ::Function,
                      model::NTuple{3, Bool})
      @inbounds begin
        $eqs
      end
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
    exp_expr(k    ::Int64,
             h    ::Int64,
             ny   ::Int64,
             model::NTuple{3,Bool})

Return exponential expression for one time evaluation of covariates.
"""
function exp_expr(k    ::Int64,
                  h    ::Int64,
                  ny   ::Int64,
                  model::NTuple{3,Bool})


  # if speciation
  if model[1] 
    pky = isone(ny) ? 1 : div(ny,k)
    ex = quote end
    pop!(ex.args)
    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(pky)
        rex = isone(ny) ? :(r[1]) : :(r[$(yi+pky*(i-1))])
        push!(coex.args,
          :(p[$(h*(3k+k*(k-1)+2)+yi+pky*(i-1)+pky*k*(j-1))] * 
            $rex))
      end
      push!(ex.args, :(eaft[$(i + k*(j-1))] = p[$(i + (k+1)*(j-1))]*exp($coex)))
    end

    return ex

  # if extinction
  elseif  model[2]
    pky = isone(ny) ? 1 : div(ny,k)
    ex = quote end
    pop!(ex.args)
    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(pky)
        rex = isone(ny) ? :(r[1]) : :(r[$(yi+pky*(i-1))])
        push!(coex.args,
          :(p[$(h*(3k+k*(k-1)+2)+yi+pky*(i-1)+pky*k*(j-1))] * 
            $rex))
      end
      push!(ex.args, :(eaft[$(i + k*(j-1))] = p[$(h*(k+1) + i + k*(j-1))]*
                       exp($coex)))
    end

    return ex

  # if dispersal
  elseif model[3]
    pky = isone(ny) ? 1 : div(ny,k*(k-1))
    ex = quote end
    pop!(ex.args)
    for j = Base.OneTo(h), i = Base.OneTo(k*(k-1))
      coex = Expr(:call, :+)
      for yi in Base.OneTo(pky)
        rex = isone(ny) ? :(r[1]) : :(r[$(yi+pky*(i-1))])
        push!(coex.args,
          :(p[$(h*(3k+k*(k-1)+2)+yi+pky*(i-1)+pky*k*(k-1)*(j-1))] * 
            $rex))
      end
      push!(ex.args, :(eaft[$(i + k*(k-1)*(j-1))] = p[$(h*(2k+1) + i + k*(k-1)*(j-1))]*
                       exp($coex)))
    end

    return ex
  end
end






"""
    noevents_expr(si   ::Int64,
                  s    ::ghs, 
                  ls   ::Int64,
                  oa   ::Array{Int64,1},
                  k    ::Int64,
                  h    ::Int64,
                  ns   ::Int64,
                  model::NTuple{3,Bool},
                  ext  ::Bool)

Return expression for no events.
"""
function noevents_expr(si   ::Int64,
                       s    ::ghs, 
                       ls   ::Int64,
                       oa   ::Array{Int64,1},
                       k    ::Int64,
                       h    ::Int64,
                       ns   ::Int64,
                       model::NTuple{3,Bool},
                       ext  ::Bool)

  # if extinction
  wu = ext ? ns : 0

  ts = isone(ls) ? 0 : k*h*k

  # between-region speciation
  ex = :(+ ($(2.0^(ls-1) - 1.) * p[$(k+1+s.h*(k+1))]))
  for (i, v) = enumerate(s.g)

    ##Covariates
    # if speciation model
    if model[1]
      push!(ex.args, :(eaft[$(v + s.h*k)] + p[$(v + (k+1)*h + s.h*k + ts)]))
    # if extinction model and only for endemic extinction
    elseif model[2] && isone(ls)
      push!(ex.args, :(p[$(v + s.h*(k+1))] + eaft[$(v + s.h*k)]))
    # if neither
    else
      push!(ex.args, :(p[$(v + s.h*(k+1))] + p[$(v + (k+1)*h + s.h*k + ts)]))
    end

    # dispersal
    if model[3]
      for j = oa
        j -= v <= j ? 1 : 0
        push!(ex.args[i+2].args, 
          :(eaft[$(s.h*k*(k-1) + (k-1)*(v-1) + j)]))
      end
    else
      for j = oa
        j -= v <= j ? 1 : 0
        push!(ex.args[i+2].args, 
          :(p[$(2h*k + h + s.h*k*(k-1) + (k-1)*(v-1) + j)]))
      end
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
                     u[$(findfirst(x -> isghsequal(x, 
                                        ghs(setdiff(s.g, i),s.h)), S))]))
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
                   mdQ::Bool,
                   ext::Bool)

Return expression for dispersal.
"""
function dispersal_expr(s  ::ghs,
                        oa ::Array{Int64,1},
                        S  ::Array{ghs,1},
                        ns ::Int64,
                        k  ::Int64,
                        h  ::Int64,
                        mdQ::Bool,
                        ext::Bool)

  wu = ext ? ns : 0

  ida = findall(x -> length(union(s.g, x.g))     == length(s.g)+1 && 
                     length(intersect(s.g, x.g)) == length(s.g)   &&
                     s.h == x.h, S)

  ex = Expr(:call, :+)

  # if dispersal covariate
  if mdQ
    for a = s.g, (i, j) = enumerate(oa)
      j -= a <= j ? 1 : 0
      push!(ex.args, :(eaft[$(s.h*k*(k-1) + (k-1)*(a-1) + j)] * 
                       u[$(ida[i] + wu)]))
    end
  else
    for a = s.g, (i, j) = enumerate(oa)
      j -= a <= j ? 1 : 0
      push!(ex.args, :(p[$(2h*k + h + s.h*k*(k-1) + (k-1)*(a-1) + j)] * 
                       u[$(ida[i] + wu)]))
    end
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
                     k  ::Int64,
                     mdS::Bool)
  # if model speciation
  if mdS
    if isone(length(s.g))
      ex = Expr(:call, :*, 2.0)
      for i = s.g
        push!(ex.args, :(eaft[$(i + s.h*k)] * 
               u[$(i + s.h*(2^k-1) + ns)] * u[$(i + s.h*(2^k-1))]))
      end
    else
      ex = Expr(:call, :+)
      for i = s.g
        push!(ex.args, :(eaft[$(i + s.h*k)] * 
          (u[$(i + s.h*(2^k-1) + ns)] * u[$(si)] + 
           u[$(si + ns)] * u[$(i + s.h*(2^k-1))])))
      end
    end

  # if model *NOT* speciation
  else
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
function ext_expr(s  ::ghs,
                  S  ::Array{ghs,1},
                  ns ::Int64,
                  k  ::Int64,
                  h  ::Int64,
                  mdE::Bool)

  if isone(length(s.g))
    if mdE
      ex = Expr(:call, :+)
      for i in s.g push!(ex.args, :(eaft[$(k*s.h + i)])) end
    else
      ex = Expr(:call, :+)
      for i in s.g push!(ex.args, :(p[$((k+1)*h + k*s.h + i)])) end
    end
  else
    ex = Expr(:call, :+)
    for i = s.g
      push!(ex.args, :(p[$(i + (k+1)*h + s.h*k + k*h*k)] * 
                       u[$(findfirst(x -> isghsequal(x, 
                                          ghs(setdiff(s.g, i),s.h)), S) + ns)]))
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
function wrsext_expr(si ::Int64,
                     s  ::ghs,
                     ns ::Int64,
                     k  ::Int64,
                     mdS::Bool)

  # if speciation model
  if mdS
    if isone(length(s.g))
      ex = Expr(:call, :*, 2.0)
      for i = s.g
        push!(ex.args, :(eaft[$(i + s.h*k)] * u[$(i + s.h*(2^k-1) + ns)]^2))
      end
    else
      ex = Expr(:call, :+)
      for i = s.g
        push!(ex.args, :(eaft[$(i + s.h*k)] * 
                         u[$(i + s.h*(2^k-1) + ns)] * u[$(si + ns)]))
      end
    end
  # if *NOT* speciation model
  else
    if isone(length(s.g)) 
      ex = :(p[$si] * u[$(si + ns)]^2)
    else
      ex = Expr(:call, :+)
      for i = s.g
        push!(ex.args, :(p[$(i + s.h*(k+1))] * 
                         u[$(i + s.h*(2^k-1) + ns)] * u[$(si + ns)]))
      end
    end
  end

  return ex
end





"""
    egeohisse_2k_gen_s(k  ::Int64,
                       h  ::Int64,
                       ny ::Int64,
                       af!::Function)

EGeoHiSSE + extinction ODE equation for `2` areas, `2` hidden states 
and `2` ny for `speciation` model.
"""
function egeohisse_2k_gen_s(k  ::Int64,
                            h  ::Int64,
                            ny ::Int64,
                            af!::Function)

  r    = Array{Float64,1}(undef, 4)
  eaft = Array{Float64,1}(undef, 4)

  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)
    @inbounds begin
      af!(t, r)
      eaft[1] = p[1] * exp(p[21] * r[1] + p[22] * r[2])
      eaft[2] = p[2] * exp(p[23] * r[3] + p[24] * r[4])
      eaft[3] = p[4] * exp(p[25] * r[1] + p[26] * r[2])
      eaft[4] = p[5] * exp(p[27] * r[3] + p[28] * r[4])
      du[1] = -1.0 * (eaft[1] + p[7] + p[11] + p[19]) * u[1] + +(p[11] * u[3]) + +(p[19] * u[4]) + 2.0 * (eaft[1] * u[7] * u[1])
      du[2] = -1.0 * (eaft[2] + p[8] + p[12] + p[19]) * u[2] + +(p[12] * u[3]) + +(p[19] * u[5]) + 2.0 * (eaft[2] * u[8] * u[2])
      du[3] = -1.0 * (1.0 * p[3] + (eaft[2] + p[16] + p[19]) + (eaft[1] + p[15])) * u[3] + (p[16] * u[1] + p[15] * u[2]) + +(p[19] * u[6]) + (eaft[2] * (u[8] * u[3] + u[9] * u[2]) + eaft[1] * (u[7] * u[3] + u[9] * u[1])) + p[3] * (u[7] * u[2] + u[8] * u[1])
      du[4] = -1.0 * (eaft[3] + p[9] + p[13] + p[20]) * u[4] + +(p[13] * u[6]) + +(p[20] * u[1]) + 2.0 * (eaft[3] * u[10] * u[4])
      du[5] = -1.0 * (eaft[4] + p[10] + p[14] + p[20]) * u[5] + +(p[14] * u[6]) + +(p[20] * u[2]) + 2.0 * (eaft[4] * u[11] * u[5])
      du[6] = -1.0 * (1.0 * p[6] + (eaft[4] + p[18] + p[20]) + (eaft[3] + p[17])) * u[6] + (p[18] * u[4] + p[17] * u[5]) + +(p[20] * u[3]) + (eaft[4] * (u[11] * u[6] + u[12] * u[5]) + eaft[3] * (u[10] * u[6] + u[12] * u[4])) + p[6] * (u[10] * u[5] + u[11] * u[4])
      du[7] = -1.0 * (eaft[1] + p[7] + p[11] + p[19]) * u[7] + +(p[7]) + +(p[11] * u[9]) + +(p[19] * u[10]) + 2.0 * (eaft[1] * u[7] ^ 2)
      du[8] = -1.0 * (eaft[2] + p[8] + p[12] + p[19]) * u[8] + +(p[8]) + +(p[12] * u[9]) + +(p[19] * u[11]) + 2.0 * (eaft[2] * u[8] ^ 2)
      du[9] = -1.0 * (1.0 * p[3] + (eaft[2] + p[16] + p[19]) + (eaft[1] + p[15])) * u[9] + (p[16] * u[7] + p[15] * u[8]) + +(p[19] * u[12]) + (eaft[2] * u[8] * u[9] + eaft[1] * u[7] * u[9]) + p[3] * +(u[7] * u[8])
      du[10] = -1.0 * (eaft[3] + p[9] + p[13] + p[20]) * u[10] + +(p[9]) + +(p[13] * u[12]) + +(p[20] * u[7]) + 2.0 * (eaft[3] * u[10] ^ 2)
      du[11] = -1.0 * (eaft[4] + p[10] + p[14] + p[20]) * u[11] + +(p[10]) + +(p[14] * u[12]) + +(p[20] * u[8]) + 2.0 * (eaft[4] * u[11] ^ 2)
      du[12] = -1.0 * (1.0 * p[6] + (eaft[4] + p[18] + p[20]) + (eaft[3] + p[17])) * u[12] + (p[18] * u[10] + p[17] * u[11]) + +(p[20] * u[9]) + (eaft[4] * u[11] * u[12] + eaft[3] * u[10] * u[12]) + p[6] * +(u[10] * u[11])
    end
    return nothing
  end
end





"""
    egeohisse_2k_gen_e(k  ::Int64,
                       h  ::Int64,
                       ny ::Int64,
                       af!::Function)

EGeoHiSSE + extinction ODE equation for `2` areas, `2` hidden states and 
`2` ny for `extinction` model.
"""
function egeohisse_2k_gen_e(k  ::Int64,
                            h  ::Int64,
                            ny ::Int64,
                            af!::Function)

  r    = Array{Float64,1}(undef, 4)
  eaft = Array{Float64,1}(undef, 4)

  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)
    @inbounds begin
      af!(t, r)
      eaft[1] = p[7] * exp(p[21] * r[1] + p[22] * r[2])
      eaft[2] = p[8] * exp(p[23] * r[3] + p[24] * r[4])
      eaft[3] = p[9] * exp(p[25] * r[1] + p[26] * r[2])
      eaft[4] = p[10] * exp(p[27] * r[3] + p[28] * r[4])
      du[1] = -1.0 * (p[1] + eaft[1] + p[11] + p[19]) * u[1] + +(p[11] * u[3]) + +(p[19] * u[4]) + 2.0 * p[1] * u[7] * u[1]
      du[2] = -1.0 * (p[2] + eaft[2] + p[12] + p[19]) * u[2] + +(p[12] * u[3]) + +(p[19] * u[5]) + 2.0 * p[2] * u[8] * u[2]
      du[3] = -1.0 * (1.0 * p[3] + (p[2] + p[16] + p[19]) + (p[1] + p[15])) * u[3] + (p[16] * u[1] + p[15] * u[2]) + +(p[19] * u[6]) + (p[2] * (u[8] * u[3] + u[9] * u[2]) + p[1] * (u[7] * u[3] + u[9] * u[1])) + p[3] * (u[7] * u[2] + u[8] * u[1])
      du[4] = -1.0 * (p[4] + eaft[3] + p[13] + p[20]) * u[4] + +(p[13] * u[6]) + +(p[20] * u[1]) + 2.0 * p[4] * u[10] * u[4]
      du[5] = -1.0 * (p[5] + eaft[4] + p[14] + p[20]) * u[5] + +(p[14] * u[6]) + +(p[20] * u[2]) + 2.0 * p[5] * u[11] * u[5]
      du[6] = -1.0 * (1.0 * p[6] + (p[5] + p[18] + p[20]) + (p[4] + p[17])) * u[6] + (p[18] * u[4] + p[17] * u[5]) + +(p[20] * u[3]) + (p[5] * (u[11] * u[6] + u[12] * u[5]) + p[4] * (u[10] * u[6] + u[12] * u[4])) + p[6] * (u[10] * u[5] + u[11] * u[4])
      du[7] = -1.0 * (p[1] + eaft[1] + p[11] + p[19]) * u[7] + +(eaft[1]) + +(p[11] * u[9]) + +(p[19] * u[10]) + p[1] * u[7] ^ 2
      du[8] = -1.0 * (p[2] + eaft[2] + p[12] + p[19]) * u[8] + +(eaft[2]) + +(p[12] * u[9]) + +(p[19] * u[11]) + p[2] * u[8] ^ 2
      du[9] = -1.0 * (1.0 * p[3] + (p[2] + p[16] + p[19]) + (p[1] + p[15])) * u[9] + (p[16] * u[7] + p[15] * u[8]) + +(p[19] * u[12]) + (p[2] * u[8] * u[9] + p[1] * u[7] * u[9]) + p[3] * +(u[7] * u[8])
      du[10] = -1.0 * (p[4] + eaft[3] + p[13] + p[20]) * u[10] + +(eaft[3]) + +(p[13] * u[12]) + +(p[20] * u[7]) + p[4] * u[10] ^ 2
      du[11] = -1.0 * (p[5] + eaft[4] + p[14] + p[20]) * u[11] + +(eaft[4]) + +(p[14] * u[12]) + +(p[20] * u[8]) + p[5] * u[11] ^ 2
      du[12] = -1.0 * (1.0 * p[6] + (p[5] + p[18] + p[20]) + (p[4] + p[17])) * u[12] + (p[18] * u[10] + p[17] * u[11]) + +(p[20] * u[9]) + (p[5] * u[11] * u[12] + p[4] * u[10] * u[12]) + p[6] * +(u[10] * u[11])
    end
    return nothing
  end
end





"""
    egeohisse_2k_gen_q(k  ::Int64,
                       h  ::Int64,
                       ny ::Int64,
                       af!::Function)

EGeoHiSSE + extinction ODE equation for `2` areas, `2` hidden states and 
`2` ny for `transition` model.
"""
function egeohisse_2k_gen_q(k  ::Int64,
                            h  ::Int64,
                            ny ::Int64,
                            af!::Function)

  r    = Array{Float64,1}(undef, 4)
  eaft = Array{Float64,1}(undef, 4)

  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)
    @inbounds begin
      af!(t, r)
      eaft[1] = p[11] * exp(p[21] * r[1] + p[22] * r[2])
      eaft[2] = p[12] * exp(p[23] * r[3] + p[24] * r[4])
      eaft[3] = p[13] * exp(p[25] * r[1] + p[26] * r[2])
      eaft[4] = p[14] * exp(p[27] * r[3] + p[28] * r[4])
      du[1] = -1.0 * (p[1] + p[7] + eaft[1] + p[19]) * u[1] + +(eaft[1] * u[3]) + +(p[19] * u[4]) + 2.0 * p[1] * u[7] * u[1]
      du[2] = -1.0 * (p[2] + p[8] + eaft[2] + p[19]) * u[2] + +(eaft[2] * u[3]) + +(p[19] * u[5]) + 2.0 * p[2] * u[8] * u[2]
      du[3] = -1.0 * (1.0 * p[3] + (p[2] + p[16] + p[19]) + (p[1] + p[15])) * u[3] + (p[16] * u[1] + p[15] * u[2]) + +(p[19] * u[6]) + (p[2] * (u[8] * u[3] + u[9] * u[2]) + p[1] * (u[7] * u[3] + u[9] * u[1])) + p[3] * (u[7] * u[2] + u[8] * u[1])
      du[4] = -1.0 * (p[4] + p[9] + eaft[3] + p[20]) * u[4] + +(eaft[3] * u[6]) + +(p[20] * u[1]) + 2.0 * p[4] * u[10] * u[4]
      du[5] = -1.0 * (p[5] + p[10] + eaft[4] + p[20]) * u[5] + +(eaft[4] * u[6]) + +(p[20] * u[2]) + 2.0 * p[5] * u[11] * u[5]
      du[6] = -1.0 * (1.0 * p[6] + (p[5] + p[18] + p[20]) + (p[4] + p[17])) * u[6] + (p[18] * u[4] + p[17] * u[5]) + +(p[20] * u[3]) + (p[5] * (u[11] * u[6] + u[12] * u[5]) + p[4] * (u[10] * u[6] + u[12] * u[4])) + p[6] * (u[10] * u[5] + u[11] * u[4])
      du[7] = -1.0 * (p[1] + p[7] + eaft[1] + p[19]) * u[7] + +(p[7]) + +(eaft[1] * u[9]) + +(p[19] * u[10]) + p[1] * u[7] ^ 2
      du[8] = -1.0 * (p[2] + p[8] + eaft[2] + p[19]) * u[8] + +(p[8]) + +(eaft[2] * u[9]) + +(p[19] * u[11]) + p[2] * u[8] ^ 2
      du[9] = -1.0 * (1.0 * p[3] + (p[2] + p[16] + p[19]) + (p[1] + p[15])) * u[9] + (p[16] * u[7] + p[15] * u[8]) + +(p[19] * u[12]) + (p[2] * u[8] * u[9] + p[1] * u[7] * u[9]) + p[3] * +(u[7] * u[8])
      du[10] = -1.0 * (p[4] + p[9] + eaft[3] + p[20]) * u[10] + +(p[9]) + +(eaft[3] * u[12]) + +(p[20] * u[7]) + p[4] * u[10] ^ 2
      du[11] = -1.0 * (p[5] + p[10] + eaft[4] + p[20]) * u[11] + +(p[10]) + +(eaft[4] * u[12]) + +(p[20] * u[8]) + p[5] * u[11] ^ 2
      du[12] = -1.0 * (1.0 * p[6] + (p[5] + p[18] + p[20]) + (p[4] + p[17])) * u[12] + (p[18] * u[10] + p[17] * u[11]) + +(p[20] * u[9]) + (p[5] * u[11] * u[12] + p[4] * u[10] * u[12]) + p[6] * +(u[10] * u[11])
    end
    return nothing
  end
end


