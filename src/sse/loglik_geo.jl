#=

Log-likelihood for ODEs, integrating along the tree

Ignacio Quintero Mächler

t(-_-t)

Created 20 09 2017
Updated 26 03 2019

=#





"""
    make_lhf(llf, prf)

Make log posterior function with the likelihood, **llf**, 
and prior, **lpf**, functions.
"""
function make_lhf(llf::Function, 
                  lpf::Function, 
                  conp::Dict{Int64,Int64})

  function f(p::Array{Float64,1}, wp::Int64)
    while haskey(conp, wp)
      tp = conp[wp]
      p[tp] = p[wp]
      wp = tp
    end
    return llf(p) + lpf(p)
  end

  return f
end





"""
    make_llf(tip_val::Dict{Int64,Array{Float64,1}},
             ed     ::Array{Int64,2},
             el     ::Array{Float64,1},
             ode_fun,
             af!    ::Function,
             p      ::Array{Float64,1},
             h      ::Int64,
             ny     ::Int64,
             model  ::Int64)

Make likelihood function for a tree given an ODE function.
"""
function make_llf(tip_val::Dict{Int64,Array{Float64,1}},
                  ed     ::Array{Int64,2},
                  el     ::Array{Float64,1},
                  ode_solve,
                  af!    ::Function,
                  ::Val{k},
                  ::Val{h},
                  ::Val{ny},
                  ::Val{model}) where {k, h, ny, model}

  ns   = h*(2^k-1)
  ntip = length(tip_val)

  # add root of length 0
  ed = cat(ed, [2*ntip ntip + 1], dims = 1)
  push!(el, 0.0)

  # get absolute times of branches as related to z(t)
  elrt = abs_time_branches(el, ed, ntip)

  elrt1 = elrt[:,1]
  elrt2 = elrt[:,2]

  ne    = size(ed,1)
  child = ed[:,2]
  wtp   = findall(child .<= ntip)

  # make trios
  trios = maketriads(ed, rev = true)

  # preallocate tip likelihoods
  led = Array{Float64,1}[]
  for i in Base.OneTo(ne)
    push!(led, zeros(Float64, 2*ns))
  end

  # create states 
  S = create_states(k, h)

  # assign states to terminal branches
  for wi in wtp 
    wig = Set(findall(map(x -> isone(x), tip_val[child[wi]])))
    led[wi][findall(map(x -> isequal(x.g, wig), S))] .= 1.0
  end

  # partial log-likelihoods
  llik = Array{Float64,1}(undef,ns)

  # log lambdas 
  lλs = Array{Float64,1}(undef,h*(k+1))

  # log lambdas covariates
  lλts = Array{Float64,1}(undef,h*k)

  # preallocate output of y covariates
  r = Array{Float64,1}(undef,ny)

  # preallocate weights
  w = Array{Float64,1}(undef,ns)

  # preallocate extinction probabilities at root
  extp = Array{Float64,1}(undef, ns)

  # make speciation events and closure
  λevent! = (t   ::Float64, 
             llik::Array{Float64,1}, 
             ud1 ::Array{Float64,1}, 
             ud2 ::Array{Float64,1},
             lλs ::Array{Float64,1},
             lλts::Array{Float64,1},
             p   ::Array{Float64,1}) ->
    begin
      λevent_full(t, llik, ud1, ud2, lλs, lλts, p, r, af!,
        Val(k), Val(h), Val(ny), Val(model))
      return nothing
    end

  # make root full likelihood estimation and closure
  rootll = (t   ::Float64,
            llik::Array{Float64,1},
            extp::Array{Float64,1},
            w   ::Array{Float64,1},
            p   ::Array{Float64,1},
            lλs ::Array{Float64,1},
            lλts::Array{Float64,1}) -> 
    begin
      rootll_full(t, llik, extp, w, p, lλs, lλts, r, af!,
        Val(k), Val(h), Val(ny), Val(model))
    end

  function f(p::Array{Float64,1})

    @inbounds begin

      llxtra = 0.0

      for i in Base.OneTo(h*(k+1))
        lλs[i] = log(p[i])
      end

      # loop for integrating over internal branches
      for triad in trios

        pr, d1, d2 = triad::Array{Int64,1}

        ud1 = ode_solve(led[d1], p, elrt2[d1], elrt1[d1])::Array{Float64,1}

        check_negs(ud1, ns) && return -Inf

        ud2 = ode_solve(led[d2], p, elrt2[d2], elrt1[d2])::Array{Float64,1}

        check_negs(ud2, ns) && return -Inf

        # update likelihoods with speciation event
        λevent!(elrt2[pr], llik, ud1, ud2, lλs, lλts, p)

        # loglik to sum for integration
        tosum   = minimum(llik)/2.0
        llxtra -= tosum

        # assign the remaining likelihoods &
        # assign extinction probabilities and 
        # check for extinction of `1.0`
        for i in Base.OneTo(ns)
          isone(ud1[i+ns]) && return -Inf
          led[pr][i+ns] = ud1[i+ns]
          led[pr][i] = exp(llik[i] - tosum)
        end
      end

      # assign root likelihood in non log terms &
      # assign root extinction probabilities
      for i in Base.OneTo(ns)
        llik[i] = led[ne][i]
        extp[i] = led[ne][i+ns]
      end

      # estimate likelihood weights
      normbysum!(llik, w, ns)

      # combine root likelihoods
      ll = rootll(elrt1[ne], llik, extp, w, p, lλs, lλts)

      return (log(ll) - llxtra)::Float64
    end
  end

  return f
end






"""
    make_rootll(t   ::Float64,
                llik::Array{Float64,1},
                extp::Array{Float64,1},
                w   ::Array{Float64,1},
                p   ::Array{Float64,1},
                lλs ::Array{Float64,1},
                lλts::Array{Float64,1},
                r   ::Array{Float64,1},
                ::Val{k}, 
                ::Val{h}, 
                ::Val{ny}, 
                ::Val{mdS}) where {k,h,ny,mdS}
Generated function for full tree likelihood at the root.
"""
@generated function rootll_full(t   ::Float64,
                                llik::Array{Float64,1},
                                extp::Array{Float64,1},
                                w   ::Array{Float64,1},
                                p   ::Array{Float64,1},
                                lλs ::Array{Float64,1},
                                lλts::Array{Float64,1},
                                r   ::Array{Float64,1},
                                af! ::Function,
                                ::Val{k},
                                ::Val{h},
                                ::Val{ny},
                                ::Val{model})::Float64 where {k, h, ny, model}

  S = create_states(k, h)

  eqs = quote end
  popfirst!(eqs.args)

  if model == 1
    # add environmental function
    push!(eqs.args, :(af!(t, r)))

    # estimate covariate lambdas
    pky = isone(ny) ? 1 : div(ny,k)
    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(pky)
        rex = isone(ny) ? :(r[1]) : :(r[$(yi+pky*(i-1))])
        push!(coex.args,
          :(p[$(h*(k+1) + 2*k*h + k*(k-1)*h + h*(h-1) + yi + pky* ((i-1) + k*(j-1)))] * 
            $rex))
      end
      push!(eqs.args, 
        :(lλts[$(i+k*(j-1))] = lλs[$(i+(k+1)*(j-1))] + $coex))
    end

    eq = Expr(:call, :+)
    for j in Base.OneTo(h)

      # for single areas
      for i in Base.OneTo(k)
        push!(eq.args,
          :(llik[$(i + (2^k-1)*(j-1))]*w[$(i + (2^k-1)*(j-1))] / 
            (exp(lλts[$(i + k*(j-1))])*(1.0-extp[$(i + (j-1)*(2^k-1))])^2)))
      end

      # for widespread
      for i in k+1:2^k-1
        wl = :(llik[$(i + (2^k-1)*(j-1))]*w[$(i + (2^k-1)*(j-1))] / 
          (l*(1.0-extp[$(i + (j-1)*(2^k-1))])^2))

        lams = Expr(:call, :+, 
          :($(2.0^(length(S[i + (2^k-1)*(j-1)].g)) - 2)*
            p[$((k+1)*(1+(j-1)))]))
        for a in S[i + (2^k-1)*(j-1)].g
          push!(lams.args, :(2.0 * exp(lλts[$(a + k*(j-1))])))
        end

        # change l in wl
        wl.args[3].args[2] = lams

        push!(eq.args, wl)
      end
    end

  else

    eq = Expr(:call, :+)
    for j in Base.OneTo(h)
      # for single areas
      for i in Base.OneTo(k)
        push!(eq.args,
          :(llik[$(i + (2^k-1)*(j-1))]*w[$(i + (2^k-1)*(j-1))] / 
            (p[$(i + (k+1)*(j-1))]*(1.0-extp[$(i + (j-1)*(2^k-1))])^2)))
      end
      # for widespread
      for i in k+1:2^k-1
        wl = :(llik[$(i + (2^k-1)*(j-1))]*w[$(i + (2^k-1)*(j-1))] / 
          (l*(1.0-extp[$(i + (j-1)*(2^k-1))])^2))

        lams = Expr(:call, :+, 
          :($(2.0^(length(S[i + (2^k-1)*(j-1)].g)) - 2)*
            p[$((k+1)*(1+(j-1)))]))
        for a in S[i + (2^k-1)*(j-1)].g
          push!(lams.args, :(2.0 * p[$(a + (k+1)*(j-1))]))
        end

        # change l in wl
        wl.args[3].args[2] = lams

        push!(eq.args, wl)
      end
    end

  end

  push!(eqs.args, :($eq))

  return quote
    @inbounds begin
      $eqs
    end
  end
end










"""
    λevent_full(t   ::Float64, 
                llik::Array{Float64,1}, 
                ud1 ::Array{Float64,1}, 
                ud2 ::Array{Float64,1},
                lλs ::Array{Float64,1},
                lλts::Array{Float64,1},
                p   ::Array{Float64,1},
                r   ::Array{Float64,1},
                ::Val{k}, 
                ::Val{h}, 
                ::Val{ny}, 
                ::Val{S},
                ::Val{mdS})
Generated function for speciation event likelihoods
"""
@generated function λevent_full(t   ::Float64, 
                                llik::Array{Float64,1}, 
                                ud1 ::Array{Float64,1}, 
                                ud2 ::Array{Float64,1},
                                lλs ::Array{Float64,1},
                                lλts::Array{Float64,1},
                                p   ::Array{Float64,1},
                                r   ::Array{Float64,1},
                                af! ::Function,
                                ::Val{k}, 
                                ::Val{h}, 
                                ::Val{ny}, 
                                ::Val{model}) where {k, h, ny, model}

  eqs = quote end
  popfirst!(eqs.args)

  S = create_states(k, h)

  # if speciation model
  if model == 1 
    # add environmental function
    push!(eqs.args, :(af!(t, r)))

    # estimate covariate lambdas
    pky = isone(ny) ? 1 : div(ny,k)
    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(pky)
        rex = isone(ny) ? :(r[1]) : :(r[$(yi+pky*(i-1))])
        push!(coex.args,
          :(p[$(h*(k+1) + 2*k*h + k*(k-1)*h + h*(h-1) + yi + pky* ((i-1) + k*(j-1)))] * 
            $rex))
      end
      push!(eqs.args, 
        :(lλts[$(i+k*(j-1))] = lλs[$(i+(k+1)*(j-1))] + $coex))
    end

    # likelihood for individual areas states
    for j = Base.OneTo(h), i = Base.OneTo(k)
      push!(eqs.args, 
          :(llik[$(i+(2^k-1)*(j-1))] = 
            log(ud1[$(i+(2^k-1)*(j-1))] * ud2[$(i+(2^k-1)*(j-1))]) + 
            lλts[$(i+k*(j-1))]))
    end

    # likelihood for widespread states
    for j = Base.OneTo(h), i = k+1:2^k-1
      s = S[i + (2^k-1)*(j-1)]

      # within region speciation
      ex = Expr(:call, :+)
      for a in s.g
        push!(ex.args, :((ud1[$(a + s.h*(2^k-1))]    *  ud2[$(i + (2^k-1)*(j-1))] + 
                          ud1[$(i + (2^k-1)*(j-1))] * ud2[$(a + s.h*(2^k-1))]) *
                         exp(lλts[$(a + s.h*(k))])))
      end

      # between region speciation
      for (la, ra) = vicsubsets(s.g)[1:div(end,2)]
        push!(ex.args, 
         :((ud1[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] + 
            ud1[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)]) *
            p[$((k+1)*(1+s.h))]))
      end

      push!(eqs.args, :(llik[$(i + (2^k-1)*(j-1))] = log($ex)))
    end

  # *not* speciation model
  else

    # likelihood for individual areas states
    for j = Base.OneTo(h), i = Base.OneTo(k)
      push!(eqs.args, 
          :(llik[$(i+(2^k-1)*(j-1))] = 
            log(ud1[$(i+(2^k-1)*(j-1))] * ud2[$(i+(2^k-1)*(j-1))]) + 
            lλs[$(i+(k+1)*(j-1))]))
    end

    # likelihood for widespread states
    for j = Base.OneTo(h), i = k+1:2^k-1
      s = S[i + (2^k-1)*(j-1)]

      # within region speciation
      ex = Expr(:call, :+)
      for a in s.g
        push!(ex.args, :((ud1[$(a + s.h*(k+1))]*ud2[$(i + (2^k-1)*(j-1))] + 
                          ud1[$(i + (2^k-1)*(j-1))]*ud2[$(a + s.h*(k+1))]) *
                         p[$(a + s.h*(k+1))]))
      end

      # between region speciation
      for (la, ra) = vicsubsets(s.g)[1:div(end,2)]
        push!(ex.args, 
         :((ud1[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] + 
            ud1[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)]) *
            p[$((k+1)*(1+s.h))]))
      end

      push!(eqs.args, :(llik[$(i + (2^k-1)*(j-1))] = log($ex)))
    end
  end

  return quote 
    @inbounds begin
      $eqs
    end
    return nothing
  end
end





"""
    normbysum!(v::Array{Float64,1}, r::Array{Float64,1}, ns::Int64)

Return weights by the sum of all elements of the array to sum 1. 
"""
function normbysum!(v::Array{Float64,1}, r::Array{Float64,1}, ns::Int64)
  s = sum(v)
  @inbounds begin
    @simd for i in Base.OneTo(ns)
      r[i] = v[i]/s
    end
  end
end






"""
    lpf_full(p      ::Array{Float64,1},
             λpriors::Float64,
             μpriors::Float64,
             lpriors::Float64,
             gpriors::Float64,
             qpriors::Float64,
             βpriors::Tuple{Float64,Float64},
             k      ::Int64,
             h      ::Int64,
             mdQ    ::Int)
`@generated` log-prior function.
"""
@generated function lpf_full(p::Array{Float64,1},
                             ::Val{λpriors},
                             ::Val{μpriors},
                             ::Val{lpriors},
                             ::Val{gpriors},
                             ::Val{qpriors},
                             ::Val{βpriors},
                             ::Val{k},
                             ::Val{h},
                             ::Val{ny},
                             ::Val{model}) where {λpriors,
                                                  μpriors,
                                                  lpriors,
                                                  gpriors,
                                                  qpriors,
                                                  βpriors,
                                                  k, h, ny, model}

  eq = Expr(:call, :+)
  # speciation priors
  for i in Base.OneTo(h*(k+1))
    push!(eq.args, :(logdexp(p[$i], $λpriors)))
  end
  # global extinction priors
  for i in (h*(k+1)+1):(h*(k+1)+k*h)
    push!(eq.args, :(logdexp(p[$i], $μpriors)))
  end
  # area colonization priors
  for i in (h*(k+1)+k*h+1):(h*(k+1)+k*h+k*(k-1)*h)
    push!(eq.args, :(logdexp(p[$i], $gpriors)))
  end
  # area loss priors
  for i in (h*(k+1)+k*h+k*(k-1)*h+1):(h*(k+1)+2k*h+k*(k-1)*h)
    push!(eq.args, :(logdexp(p[$i], $lpriors)))
  end
  # hidden states transition
  for i in (h*(k+1)+2k*h+k*(k-1)*h+1):(h*(k+1)+2k*h+k*(k-1)*h+h*(h-1))
    push!(eq.args, :(logdexp(p[$i], $qpriors)))
  end
  # betas
  if isone(ny)
    nb = model == 3 ? 
      (h*(k+1)+2k*h+k*(k-1)*h+h*(h-1) + k*(k-1)) :
      (h*(k+1)+2k*h+k*(k-1)*h+h*(h-1) + k*h)
    for i in (h*(k+1)+2k*h+k*(k-1)*h+h*(h-1)+1):nb
      push!(eq.args, :(logdnorm(p[$i], $(βpriors[1]), $(βpriors[2]))))
    end
  else
    nb = model == 3 ? 
      (h*(k+1)+2k*h+k*(k-1)*h+h*(h-1) + k*(k-1)*div(ny,k*(k-1))*h) :
      (h*(k+1)+2k*h+k*(k-1)*h+h*(h-1) + k*h*div(ny,k))
    for i in (h*(k+1)+2k*h+k*(k-1)*h+h*(h-1)+1):nb
      push!(eq.args, :(logdnorm(p[$i], $(βpriors[1]), $(βpriors[2]))))
    end
  end

  return quote 
    @inbounds begin
      lq = $eq
    end
    return lq
  end
end






""" 
    make_lpf(::Val{λpriors}, 
                  ::Val{μpriors}, 
                  ::Val{lpriors}, 
                  ::Val{gpriors}, 
                  ::Val{qpriors}, 
                  ::Val{βpriors}, 
                  ::Val{k}, 
                  ::Val{h}, 
                  ::Val{ny}, 
                  ::Val{model}) where {λpriors,
                                       μpriors,
                                       lpriors,
                                       gpriors,
                                       qpriors,
                                       βpriors,
                                       k, h, ny, model}

Make log prior function
"""
function make_lpf(::Val{λpriors}, 
                  ::Val{μpriors}, 
                  ::Val{lpriors}, 
                  ::Val{gpriors}, 
                  ::Val{qpriors}, 
                  ::Val{βpriors}, 
                  ::Val{k}, 
                  ::Val{h}, 
                  ::Val{ny}, 
                  ::Val{model}) where {λpriors,
                                       μpriors,
                                       lpriors,
                                       gpriors,
                                       qpriors,
                                       βpriors,
                                       k, h, ny, model}

  # create prior function
  lpf = (p::Array{Float64,1}) ->
    lpf_full(p, Val(λpriors), Val(μpriors), Val(lpriors), 
      Val(gpriors), Val(qpriors), Val(βpriors), 
      Val(k), Val(h), Val(ny), Val(model))

  return lpf
end





"""
    check_negs(x::Array{Float64,1}, k::Int64)

Return `true` if any in x is less or equal than 0.  
"""
function check_negs(x::Array{Float64,1}, k::Int64)
  @inbounds begin
    @simd for i in Base.OneTo(k)
      if x[i] < 0.0 
        return true
      end
    end
    return false
  end
end















