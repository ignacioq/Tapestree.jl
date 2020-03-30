#=

EGeoHiSSE equations

Ignacio Quintero Mächler

t(-_-t)

Created 18 03 2019
=#





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

  # last non β parameters
  bbase = h*(k+1) + 2*k*h + k*(k-1)*h + h*(h-1)

  # ncov
  ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

  # y per parameter
  yppar = isone(ny) ? 1 : div(ny,ncov)

  # starting indices for models 2 and 3
  m2s = model[1]*k*h
  m3s = m2s + model[2]*k*h

  y2s = model[1]*yppar*k
  y3s = y2s + model[2]*yppar*k

  # start expression 
  ex = quote end
  pop!(ex.args)

  isone(ny) && push!(ex.args, :(r1 = r[1]))

  # speciation 
  if model[1]
    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(yppar)
        rex = isone(ny) ? :(r1) : :(r[$(yi+yppar*(i-1))])
        push!(coex.args, :(p[$(bbase + yi + yppar*((i-1) + k*(j-1)))] * $rex))
      end
      push!(ex.args, :(eaft[$(i + k*(j-1))] = p[$(i + (k+1)*(j-1))]*exp($coex)))
    end
  end 
  # extinction 
  if model[2]
    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(yppar)
        rex = isone(ny) ? :(r1) : :(r[$(y2s + yi+yppar*(i-1))])
        push!(coex.args, 
          :(p[$(bbase + yi + m2s*yppar + yppar*((i-1) + k*(k-1)*(j-1)))] * $rex))
      end
      push!(ex.args, :(eaft[$(i + k*(j-1) + m2s)] = 
        p[$((k+1)*h + i + k*(j-1))]*exp($coex)))
    end
  end
  # transition 
  if model[3]
    for j = Base.OneTo(h), i = Base.OneTo(k*(k-1))
      coex = Expr(:call, :+)
      for yi in Base.OneTo(yppar)
        rex = isone(ny) ? :(r1) : :(r[$(y3s + yi+yppar*(i-1))])
        push!(coex.args, 
          :(p[$(bbase + yi + m3s*yppar + yppar*((i-1) + k*(k-1)*(j-1)))] * $rex))
      end
      push!(ex.args, :(eaft[$(i + k*(j-1) + m3s)] = 
        p[$(h*(2k+1) + i + k*(k-1)*(j-1))]*exp($coex)))
    end
  end

  return ex
end





"""
    noevents_expr(si   ::Int64,
                  s    ::Sgh, 
                  ls   ::Int64,
                  oa   ::Array{Int64,1},
                  k    ::Int64,
                  h    ::Int64,
                  ns   ::Int64,
                  model::NTuple{3,Bool},
                  ext  ::Bool)

Return expression for no events in Vector form.
"""
function noevents_expr(si   ::Int64,
                       s    ::Sgh, 
                       ls   ::Int64,
                       oa   ::Array{Int64,1},
                       k    ::Int64,
                       h    ::Int64,
                       ns   ::Int64,
                       model::NTuple{3,Bool},
                       ext  ::Bool)

  # starting indices for models 2 and 3
  m2s = model[1]*k*h
  m3s = m2s + model[2]*k*h

  # if extinction ODE
  wu = ext ? ns : 0

  ts = isone(ls) ? 0 : k*h*k

  # between-region speciation
  ex = Expr(:call, :+)
  ls > 1 && push!(ex.args, :($(2.0^(ls-1) - 1.0) * p[$(k+1+s.h*(k+1))]))

  for v in s.g

    #speciation
    if model[1]
      push!(ex.args, :(eaft[$(v + s.h*k)]))
    else
      push!(ex.args, :(p[$(v + s.h*(k+1))]))
    end

    # extinction 
    if model[2] 
      push!(ex.args, :(eaft[$(m2s + v + s.h*k)]))
    else
      push!(ex.args, :(p[$(v + (k+1)*h + s.h*k + ts)]))
    end

    # dispersal
    if model[3]
      for j = oa
        j -= v <= j ? 1 : 0
        push!(ex.args, 
          :(eaft[$(m3s + s.h*k*(k-1) + (k-1)*(v-1) + j)]))
      end
    else
      for j = oa
        j -= v <= j ? 1 : 0
        push!(ex.args, 
          :(p[$(2h*k + h + s.h*k*(k-1) + (k-1)*(v-1) + j)]))
      end
    end
  end

  # add hidden state shifts
  for hi in setdiff(0:(h-1), s.h)
    hi -= s.h <= hi ? 0 : -1
    push!(ex.args, :(p[$(h*(3k + 1 + k*(k-1)) + s.h*(h-1) + hi)]))
  end

  # multiply by u
  ex = :(-1.0 * $ex * u[$(si + wu)])

  return ex
end





"""
    noevents_expr(si     ::Int64,
                  s      ::Sgh, 
                  ls     ::Int64,
                  oa     ::Array{Int64,1},
                  k      ::Int64,
                  h      ::Int64,
                  ns     ::Int64,
                  i      ::Int64,
                  model  ::NTuple{3,Bool},
                  ext    ::Bool)

Return expression for no events in Matrix form.
"""
function noevents_expr(si     ::Int64,
                       s      ::Sgh, 
                       ls     ::Int64,
                       oa     ::Array{Int64,1},
                       k      ::Int64,
                       h      ::Int64,
                       ns     ::Int64,
                       i      ::Int64,
                       model  ::NTuple{3,Bool})

  # starting indices for models 2 and 3
  m2s = model[1]*k*h
  m3s = m2s + model[2]*k*h


  ts = isone(ls) ? 0 : k*h*k

  # between-region speciation
  ex = Expr(:call, :+)
  ls > 1 && push!(ex.args, :($(2.0^(ls-1) - 1.0) * p[$(k+1+s.h*(k+1))]))

  for v in s.g

    #speciation
    if model[1]
      push!(ex.args, :(eaft[$(v + s.h*k)]))
    else
      push!(ex.args, :(p[$(v + s.h*(k+1))]))
    end

    # extinction 
    if model[2] 
      push!(ex.args, :(eaft[$(m2s + v + s.h*k)]))
    else
      push!(ex.args, :(p[$(v + (k+1)*h + s.h*k + ts)]))
    end

    # dispersal
    if model[3]
      for j = oa
        j -= v <= j ? 1 : 0
        push!(ex.args, 
          :(eaft[$(m3s + s.h*k*(k-1) + (k-1)*(v-1) + j)]))
      end
    else
      for j = oa
        j -= v <= j ? 1 : 0
        push!(ex.args, 
          :(p[$(2h*k + h + s.h*k*(k-1) + (k-1)*(v-1) + j)]))
      end
    end
  end

  # add hidden state shifts
  for hi in setdiff(0:(h-1), s.h)
    hi -= s.h <= hi ? 0 : -1
    push!(ex.args, :(p[$(h*(3k + 1 + k*(k-1)) + s.h*(h-1) + hi)]))
  end

  # multiply by u
  ex = :(-1.0 * $ex * u[$si, $i])

  return ex
end





"""
    localext_expr(s ::Sgh,
                  S ::Array{Sgh,1},
                  k ::Int64,
                  h ::Int64,
                  model::NTuple{3,Bool})

Return expression for local extinction in Vector form.
"""
function localext_expr(s ::Sgh,
                       S ::Array{Sgh,1},
                       k ::Int64,
                       h ::Int64,
                       model::NTuple{3,Bool})

  ex = Expr(:call, :+)

  if model[2]
    # starting indices for models 2 and 3
    m2s = model[1]*k*h

    for i = s.g
      push!(ex.args, :(eaft[$(m2s + i + s.h*k)] * 
                       u[$(findfirst(x -> isSghequal(x, 
                                          Sgh(setdiff(s.g, i),s.h)), S))]))
    end
  else
    for i = s.g
      push!(ex.args, :(p[$(i + (k+1)*h + s.h*k + k*h*k)] * 
                       u[$(findfirst(x -> isSghequal(x, 
                                          Sgh(setdiff(s.g, i),s.h)), S))]))
    end
  end

  return ex
end





"""
    localext_expr(s     ::Sgh,
                  S     ::Array{Sgh,1},
                  k     ::Int64,
                  h     ::Int64,
                  i     ::Int64
                  model ::NTuple{3,Bool})

Return expression for local extinction in Matrix form.
"""
function localext_expr(s     ::Sgh,
                       S     ::Array{Sgh,1},
                       k     ::Int64,
                       h     ::Int64,
                       i     ::Int64,
                       model ::NTuple{3,Bool})

  ex = Expr(:call, :+)

  if model[2]
    # starting indices for models 2 and 3
    m2s = model[1]*k*h

    for j = s.g
      push!(ex.args, :(eaft[$(m2s + j + s.h*k)] * 
                       u[$(findfirst(x -> isSghequal(x, 
                                          Sgh(setdiff(s.g, j),s.h)), S)),$i]))
    end
  else
    for j = s.g
      push!(ex.args, :(p[$(j + (k+1)*h + s.h*k + k*h*k)] * 
                       u[$(findfirst(x -> isSghequal(x, 
                                          Sgh(setdiff(s.g, j),s.h)), S)),$i]))
    end
  end

  return ex
end





"""
    dispersal_expr(s    ::Sgh,
                   oa   ::Array{Int64,1},
                   S    ::Array{Sgh,1},
                   ns   ::Int64,
                   k    ::Int64,
                   h    ::Int64,
                   model::NTuple{3,Bool},
                   ext  ::Bool)

Return expression for dispersal in Vector form.
"""
function dispersal_expr(s    ::Sgh,
                        oa   ::Array{Int64,1},
                        S    ::Array{Sgh,1},
                        ns   ::Int64,
                        k    ::Int64,
                        h    ::Int64,
                        model::NTuple{3,Bool},
                        ext  ::Bool)

  m3s = model[1]*k*h + model[2]*k*h

  wu = ext ? ns : 0

  ida = findall(x -> length(union(s.g, x.g))     == length(s.g)+1 && 
                     length(intersect(s.g, x.g)) == length(s.g)   &&
                     s.h == x.h, S)

  ex = Expr(:call, :+)

  # if dispersal covariate
  if model[3]
    for a = s.g, (i, j) = enumerate(oa)
      j -= a <= j ? 1 : 0
      push!(ex.args, :(eaft[$(m3s + s.h*k*(k-1) + (k-1)*(a-1) + j)] * 
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
    dispersal_expr(s    ::Sgh,
                   oa   ::Array{Int64,1},
                   S    ::Array{Sgh,1},
                   ns   ::Int64,
                   k    ::Int64,
                   h    ::Int64,
                   i    ::Int64,
                   model::NTuple{3,Bool})

Return expression for dispersal in Matrix form.
"""
function dispersal_expr(s    ::Sgh,
                        oa   ::Array{Int64,1},
                        S    ::Array{Sgh,1},
                        ns   ::Int64,
                        k    ::Int64,
                        h    ::Int64,
                        i    ::Int64,
                        model::NTuple{3,Bool})

  m3s = model[1]*k*h + model[2]*k*h

  ida = findall(x -> length(union(s.g, x.g))     == length(s.g)+1 && 
                     length(intersect(s.g, x.g)) == length(s.g)   &&
                     s.h == x.h, S)

  ex = Expr(:call, :+)

  # if dispersal covariate
  if model[3]
    for a = s.g, (l, j) = enumerate(oa)
      j -= a <= j ? 1 : 0
      push!(ex.args, :(eaft[$(m3s + s.h*k*(k-1) + (k-1)*(a-1) + j)] * 
                       u[$(ida[l]), $i]))
    end
  else
    for a = s.g, (l, j) = enumerate(oa)
      j -= a <= j ? 1 : 0
      push!(ex.args, :(p[$(2h*k + h + s.h*k*(k-1) + (k-1)*(a-1) + j)] * 
                       u[$(ida[l]), $i]))
    end
  end

  return ex
end





"""
    hidtran_expr(s  ::Sgh,
                 S  ::Array{Sgh,1},
                 ns ::Int64,
                 k  ::Int64,
                 h  ::Int64,
                 i  ::Int64)

Return expression for hidden states transitions in Matrix form.
"""
function hidtran_expr(s  ::Sgh,
                      S  ::Array{Sgh,1},
                      ns ::Int64,
                      k  ::Int64,
                      h  ::Int64,
                      i  ::Int64)

  hs = findall(x -> isequal(s.g, x.g) && 
                   !isequal(s.h, x.h), S)

  ex = Expr(:call, :+)
  for j in hs
    hi = S[j].h
    hi -= s.h <= hi ? 0 : -1
    push!(ex.args, :(p[$(h*(3k + 1 + k*(k-1)) + s.h*(h-1) + hi)] * 
                     u[$j, $i]))
  end
  return ex
end





"""
    hidtran_expr(s  ::Sgh,
                 S  ::Array{Sgh,1},
                 ns ::Int64,
                 k  ::Int64,
                 h  ::Int64,
                 ext::Bool)

Return expression for hidden states transitions in Vector form.
"""
function hidtran_expr(s  ::Sgh,
                      S  ::Array{Sgh,1},
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
    wrspec_expr(si   ::Int64,
                s    ::Sgh,
                ns   ::Int64,
                k    ::Int64,
                i    ::Int64,
                model::NTuple{3,Bool})

Return expression for within-region speciation in Matrix form.
"""
function wrspec_expr(si   ::Int64,
                     s    ::Sgh,
                     ns   ::Int64,
                     k    ::Int64,
                     i    ::Int64,
                     model::NTuple{3,Bool})

  # if model speciation
  if model[1]
    if isone(length(s.g))
      ex = Expr(:call, :*, 2.0)
      for j = s.g
        push!(ex.args, :(eaft[$(j + s.h*k)] * 
               rE[$(j + s.h*(2^k-1))] * u[$(j + s.h*(2^k-1)), $i]))
      end
    else
      ex = Expr(:call, :+)
      for j = s.g
        push!(ex.args, :(eaft[$(j + s.h*k)] * 
          (rE[$(j + s.h*(2^k-1))] * u[$(si), $i] + 
           rE[$(si)] * u[$(j + s.h*(2^k-1)), $i])))
      end
    end

  # if model *NOT* speciation
  else
    if isone(length(s.g)) 
      ex = :(2.0 * p[$si] * rE[$(si)] * u[$(si)])
    else
      ex = Expr(:call, :+)
      for j = s.g
        push!(ex.args, :(p[$(j + s.h*(k+1))] * 
          (rE[$(j + s.h*(2^k-1))] * u[$(si), $i] + 
           rE[$(si)] * u[$(j + s.h*(2^k-1)), $i])))
      end
    end
  end

  return ex
end





"""
    wrspec_expr(si   ::Int64,
                s    ::Sgh,
                ns   ::Int64,
                k    ::Int64,
                model::NTuple{3,Bool})

Return expression for within-region speciation in Vector form.
"""
function wrspec_expr(si   ::Int64,
                     s    ::Sgh,
                     ns   ::Int64,
                     k    ::Int64,
                     model::NTuple{3,Bool})

  # if model speciation
  if model[1]
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
    brspec_expr(s   ::Sgh,
                S   ::Array{Sgh,1},
                ns  ::Int64,
                k   ::Int64,
                i   ::Int64, 
                ext ::Bool, 
                flow::Bool)

Return expression for between-region speciation in Vector form.
"""
function brspec_expr(s   ::Sgh,
                     S   ::Array{Sgh,1},
                     ns  ::Int64,
                     k   ::Int64,
                     ext ::Bool, 
                     flow::Bool)

  if ext
    va = vicsubsets(s.g)[1:div(end,2)]
    wu1 = flow ? 0 : ns
    wu2 = flow ? 0 : ns
  else
    va = vicsubsets(s.g)
    wu1 = flow ? 0 : ns
    wu2 = 0
  end

  ex = Expr(:call, :+)
  for (la, ra) = va
  push!(ex.args,
    :(u[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h + wu1)] *
      u[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h + wu2)]))
  end

  ex = :(p[$(k+1+s.h*(k+1))] * $ex)

  return ex
end





"""
    brspec_expr(s   ::Sgh,
                S   ::Array{Sgh,1},
                ns  ::Int64,
                k   ::Int64,
                i   ::Int64)


Return expression for between-region speciation in Matrix form.
"""
function brspec_expr(s   ::Sgh,
                     S   ::Array{Sgh,1},
                     ns  ::Int64,
                     k   ::Int64,
                     i   ::Int64)

  va = vicsubsets(s.g)
  
  ex = Expr(:call, :+)

  for (la, ra) = va
  push!(ex.args,
    :(rE[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)] *
       u[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h), $i]))
  end

  ex = :(p[$(k+1+s.h*(k+1))] * $ex)

  return ex
end





"""
    ext_expr(s    ::Sgh,
             S    ::Array{Sgh,1},
             ns   ::Int64,
             k    ::Int64,
             h    ::Int64,
             model::NTuple{3,Bool})

Return expression for extinction.
"""
function ext_expr(s    ::Sgh,
                  S    ::Array{Sgh,1},
                  ns   ::Int64,
                  k    ::Int64,
                  h    ::Int64,
                  model::NTuple{3,Bool},
                  flow ::Bool)

  wu = flow ? 0 : ns

  m2s = model[1]*k*h

  if isone(length(s.g))
    if model[2]
      ex = Expr(:call, :+)
      for i in s.g push!(ex.args, :(eaft[$(m2s + k*s.h + i)])) end
    else
      ex = Expr(:call, :+)
      for i in s.g push!(ex.args, :(p[$((k+1)*h + k*s.h + i)])) end
    end
  else
    if model[2]
      ex = Expr(:call, :+)
      for i = s.g
        push!(ex.args, :(eaft[$(m2s + k*s.h + i)] * 
                         u[$(findfirst(x -> isSghequal(x, 
                                            Sgh(setdiff(s.g, i),s.h)), S) + wu)]))
      end
    else
      ex = Expr(:call, :+)
      for i = s.g
        push!(ex.args, :(p[$(i + (k+1)*h + s.h*k + k*h*k)] * 
                         u[$(findfirst(x -> isSghequal(x, 
                                            Sgh(setdiff(s.g, i),s.h)), S) + wu)]))
      end
    end
  end

  return ex
end





"""
    wrsext_expr(si   ::Int64,
                s    ::Sgh,
                ns   ::Int64,
                k    ::Int64,
                model::NTuple{3,Bool},
                flow ::Bool)

Return expression of extinction for within-region speciation.
"""
function wrsext_expr(si   ::Int64,
                     s    ::Sgh,
                     ns   ::Int64,
                     k    ::Int64,
                     model::NTuple{3,Bool},
                     flow ::Bool)

  wu = flow ? 0 : ns

  # if speciation model
  if model[1]
    if isone(length(s.g))
      for i = s.g
        ex = :(eaft[$(i + s.h*k)] * u[$(i + s.h*(2^k-1) + wu)]^2)
      end
    else
      ex = Expr(:call, :+)
      for i = s.g
        push!(ex.args, :(eaft[$(i + s.h*k)] * 
                         u[$(i + s.h*(2^k-1) + wu)] * u[$(si + wu)]))
      end
    end
  # if *NOT* speciation model
  else
    if isone(length(s.g)) 
      ex = :(p[$si] * u[$(si + wu)]^2)
    else
      ex = Expr(:call, :+)
      for i = s.g
        push!(ex.args, :(p[$(i + s.h*(k+1))] * 
                         u[$(i + s.h*(2^k-1) + wu)] * u[$(si + wu)]))
      end
    end
  end

  return ex
end





"""
    geohisse_full(du  ::Array{Float64,1}, 
                  u   ::Array{Float64,1}, 
                  p   ::Array{Float64,1}, 
                  t   ::Float64,
                  r   ::Array{Float64,1},
                  eaft::Array{Float64,1},
                  af! ::Function,
                  ::Val{k},
                  ::Val{h},
                  ::Val{ny},
                  ::Val{model}) where {k, h, ny, model}

Creates Covariate GeoHiSSE ODE equation function for `k` areas, 
`h` hidden states and `ny` covariates.
"""
@generated function geohisse_full(du  ::Array{Float64,1}, 
                                  u   ::Array{Float64,1}, 
                                  p   ::Array{Float64,1}, 
                                  t   ::Float64,
                                  r   ::Array{Float64,1},
                                  eaft::Array{Float64,1},
                                  af! ::Function,
                                  ::Val{k},
                                  ::Val{h},
                                  ::Val{ny},
                                  ::Val{model}) where {k, h, ny, model}

  # n states
  ns = (2^k - 1)*h

  # create individual areas subsets
  S = create_states(k, h)

  # start Expression for ODE equations
  eqs = quote end
  popfirst!(eqs.args)

  # add environmental function
  push!(eqs.args, :(af!(t, r)))

  # compute exponential for each β*covariable
  expex = exp_expr(k, h, ny, model)

  for ex in expex.args
    push!(eqs.args, ex)
  end

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
      lex = localext_expr(s, S, k, h, model)

      # between-region speciation
      brs = brspec_expr(s, S, ns, k, false, false)
    end

    # dispersal
    # remove if lr == k
    if ls != k
      dis = dispersal_expr(s, oa, S, ns, k, h, model, false)
    end

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, false) : :0.0

    # within-region speciation
    wrs = wrspec_expr(si, s, ns, k, model)

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
    ext = ext_expr(s, S, ns, k, h, model, false) 

    # dispersal and extinction
    if ls != k
      dis = dispersal_expr(s, oa, S, ns, k, h, model, true)
    end

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, true) : :0.0

    # within-region extinction
    wrs = wrsext_expr(si, s, ns, k, model, false)

    # between-region extinction
    if ls > 1
      brs = brspec_expr(s, S, ns, k, true, false)
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

  #Core.println(eqs)

  return quote 
    @inbounds begin
      $eqs
    end
    return nothing
  end

end






"""
    make_egeohisse(::Val{k},
                   ::Val{h},
                   ::Val{ny},
                   ::Val{model},
                   af!::Function) where {k, h, ny, model}

Make closure for EGeoHiSSE.
"""
function make_egeohisse(::Val{k},
                        ::Val{h},
                        ::Val{ny},
                        ::Val{model},
                        af!::Function) where {k, h, ny, model}

  r    = Array{Float64,1}(undef, ny)
  eaft = Array{Float64,1}(undef,
    model[1]*k*h + model[2]*k*h + model[3]*k*(k-1)*h)

  # make ode function with closure
  ode_fun = (du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64) -> 
    begin
      geohisse_full(du::Array{Float64,1},
                    u ::Array{Float64,1},
                    p ::Array{Float64,1},
                    t ::Float64,
                    r,eaft,af!,Val(k),Val(h),Val(ny),Val(model))
      return nothing
    end

  return ode_fun
end






"""
    geohisse_E(du  ::Array{Float64,1}, 
               u   ::Array{Float64,1}, 
               p   ::Array{Float64,1}, 
               t   ::Float64,
               r   ::Array{Float64,1},
               eaft::Array{Float64,1},
               af! ::Function,
               ::Val{k},
               ::Val{h},
               ::Val{ny},
               ::Val{model}) where {k, h, ny, model}

Creates Covariate GeoHiSSE Extinction ODE equation for `k` areas, 
`h` hidden states and `ny` covariates for the specific `model`.
"""
@generated function geohisse_E(du  ::Array{Float64,1}, 
                               u   ::Array{Float64,1}, 
                               p   ::Array{Float64,1}, 
                               t   ::Float64,
                               r   ::Array{Float64,1},
                               eaft::Array{Float64,1},
                               af! ::Function,
                               ::Val{k},
                               ::Val{h},
                               ::Val{ny},
                               ::Val{model}) where {k, h, ny, model}

  # n states
  ns = (2^k - 1)*h

  # create individual areas subsets
  S = create_states(k, h)

  # start Expression for ODE equations
  eqs = quote end
  popfirst!(eqs.args)

  # add environmental function
  push!(eqs.args, :(af!(t, r)))

  # compute exponential for each β*covariable
  expex = exp_expr(k, h, ny, model)

  for ex in expex.args
    push!(eqs.args, ex)
  end

  for si = Base.OneTo(ns)

    # state range
    s = S[si]

    # length of geographic range
    ls = length(s.g)

    # which single areas do not occur in r
    oa = setdiff(1:k, s.g)

    # no events
    nev = noevents_expr(si, s, ls, oa, k, h, ns, model, false)

    # extinction
    ext = ext_expr(s, S, ns, k, h, model, true) 

    # dispersal and extinction
    if ls != k
      dis = dispersal_expr(s, oa, S, ns, k, h, model, false)
    end

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns, k, h, false) : :0.0

    # within-region extinction
    wrs = wrsext_expr(si, s, ns, k, model, true)

    # between-region extinction
    if ls > 1
      brs = brspec_expr(s, S, ns, k, true, true)
    end

    # push `E` equation to to eqs
    if isone(ls)
      push!(eqs.args, 
        :(du[$(si)] = $nev + $ext + $dis + $hid + $wrs))
    elseif ls == k
      push!(eqs.args, 
        :(du[$(si)] = $nev + $ext + $hid + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$(si)] = $nev + $ext + $dis + $hid + $wrs + $brs))
    end

  end

  # print for checking
  #Core.println(eqs)

  return quote 
    @inbounds begin
      $eqs
    end
    return nothing
  end

end





"""
    make_egeohisse_E(::Val{k},
                     ::Val{h},
                     ::Val{ny},
                     ::Val{model},
                     af!::Function) where {k, h, ny, model}

Make closure for Extinction function of EGeoHiSSE.
"""
function make_egeohisse_E(::Val{k},
                        ::Val{h},
                        ::Val{ny},
                        ::Val{model},
                        af!::Function) where {k, h, ny, model}

  r    = Array{Float64,1}(undef, ny)
  eaft = Array{Float64,1}(undef,
    model[1]*k*h + model[2]*k*h + model[3]*k*(k-1)*h)

  # make ode function with closure
  ode_fun = (du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64) -> 
    begin
      geohisse_E(du::Array{Float64,1}, 
                 u ::Array{Float64,1}, 
                 p ::Array{Float64,1}, 
                 t ::Float64,
                 r, eaft, af!, Val(k), Val(h), Val(ny), Val(model))
      return nothing
    end

  return ode_fun
end





"""
    geohisse_M(du  ::Array{Float64,1}, 
               u   ::Array{Float64,1}, 
               p   ::Array{Float64,1}, 
               t   ::Float64,
               r   ::Array{Float64,1},
               eaft::Array{Float64,1},
               af! ::Function,
               ::Val{k},
               ::Val{h},
               ::Val{ny},
               ::Val{model}) where {k, h, ny, model}

Creates Covariate GeoHiSSE Likelihood ODE equation function for `k` areas, 
`h` hidden states and `ny` covariates in Matrix form.
"""
@generated function geohisse_M(du  ::Array{Float64,2}, 
                               u   ::Array{Float64,2}, 
                               pp  ::Array{Array{Float64,1},1}, 
                               t   ::Float64,
                               r   ::Array{Float64,1},
                               eaft::Array{Float64,1},
                               rE  ::Array{Float64,1},
                               af! ::Function,
                               afE!::Function,
                               idxl::Int64,
                               ::Val{k},
                               ::Val{h},
                               ::Val{ny},
                               ::Val{model}) where {k, h, ny, model}

  # n states
  ns = (2^k - 1)*h

  # create individual areas subsets
  S = create_states(k, h)

  # start Expression for ODE equations
  eqs = quote end
  popfirst!(eqs.args)

  # add parameter subset
  push!(eqs.args, :(p = pp[idxl]))

  # add extinction function
  push!(eqs.args, :(afE!(t, rE, pp)))

  # add environmental function
  push!(eqs.args, :(af!(t, r)))

  # compute exponential for each β*covariable
  expex = exp_expr(k, h, ny, model)

  for ex in expex.args
    push!(eqs.args, ex)
  end

  for i in Base.OneTo(ns), si = Base.OneTo(ns)

    # state range
    s = S[si]

    # length of geographic range
    ls = length(s.g)

    # which single areas do not occur in r
    oa = setdiff(1:k, s.g)

    # no events
    nev = noevents_expr(si, s, ls, oa, k, h, ns, i, model)

    # remove if !isone(lr)
    if ls > 1
      # local extinction
      lex = localext_expr(s, S, k, h, i, model)

      # between-region speciation
      brs = brspec_expr(s, S, ns, k, i)
    end

    # dispersal
    # remove if lr == k
    if ls != k
      dis = dispersal_expr(s, oa, S, ns, k, h, i, model)
    end

    # hidden states transitions
    hid = h > 1 ? hidtran_expr(s, S, ns,k, h, i) : :0.0

    # within-region speciation
    wrs = wrspec_expr(si, s, ns, k, i, model)

    # push `D` equation to to eqs
    if isone(ls)
      push!(eqs.args, 
        :(du[$si, $i] = $nev + $dis + $hid + $wrs))
    elseif ls == k
      push!(eqs.args, 
        :(du[$si, $i] = $nev + $lex + $hid + $wrs + $brs))
    else
      push!(eqs.args, 
        :(du[$si, $i] = $nev + $lex + $dis + $hid + $wrs + $brs))
    end

  end

  ## aesthetic touches
  # sort `du`s

  #Core.println(eqs)

  return quote 
    @inbounds begin
      $eqs
    end
    return nothing
  end

end





"""
    make_egeohisse_M(::Val{k},
                     ::Val{h},
                     ::Val{ny},
                     ::Val{model},
                     af!::Function) where {k, h, ny, model}

Make closure for EGeoHiSSE given `E(t)` and in matrix form.
"""
function make_egeohisse_M(::Val{k},
                          ::Val{h},
                          ::Val{ny},
                          ::Val{model},
                          af! ::Function,
                          afE!::Function,
                          idxl::Int64) where {k, h, ny, model}

  ns = (2^k - 1)*h

  r    = Array{Float64,1}(undef, ny)
  eaft = Array{Float64,1}(undef,
    model[1]*k*h + model[2]*k*h + model[3]*k*(k-1)*h)
  rE  = Array{Float64,1}(undef,ns)

  # make ode function with closure
  ode_fun = (du::Array{Float64,2}, 
             u ::Array{Float64,2}, 
             pp::Array{Array{Float64,1},1}, 
             t ::Float64) -> 
    begin
      geohisse_M(du::Array{Float64,2}, 
                 u ::Array{Float64,2}, 
                 pp::Array{Array{Float64,1},1}, 
                 t ::Float64,
                 r, eaft, rE, af!, afE!, idxl,
                 Val(k), Val(h), Val(ny), Val(model))
      return nothing
    end

  return ode_fun
end



