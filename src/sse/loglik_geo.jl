#=

Log-likelihood for ODEs, integrating along the tree

Ignacio Quintero Mächler

t(-_-t)

Created 20 09 2017
Updated 26 03 2019

=#



"""
    make_lhf(llf ::Function, 
             lpf ::Function, 
             assign_hidfacs!::Function,
             dcp ::Dict{Int64,Int64},
             dcfp::Dict{Int64,Int64})

Make log posterior function with the likelihood, **llf**, 
and prior, **lpf**, functions.
"""
function make_lhf(llf            ::Function, 
                  lpf            ::Function, 
                  assign_hidfacs!::Function,
                  dcp            ::Dict{Int64,Int64},
                  dcfp           ::Dict{Int64,Int64})

  ks  = keys(dcp)
  ksf = keys(dcfp)

  function f(p ::Array{Float64,1}, 
             fp::Array{Float64,1})

    # constraints
    for wp in ks
      while haskey(dcp, wp)
        tp = dcp[wp]
        p[tp] = p[wp]
        wp = tp
      end
    end

    for wp in ksf
      while haskey(dcfp, wp)
        tp = dcfp[wp]
        fp[tp] = fp[wp]
        wp = tp
      end
    end

   # factors
    assign_hidfacs!(p, fp)

    return llf(p) + lpf(p, fp)
  end

  return f
end





"""
    make_loglik(X        ::Array{Array{Float64,1},1},
                abts1    ::Array{Float64,1},
                abts2    ::Array{Float64,1},
                trios    ::Array{Array{Int64,1},1},
                ode_solve::Function,
                λevent!  ::Function, 
                rootll   ::Function,
                k        ::Int64,
                h        ::Int64,
                ns       ::Int64)

Make likelihood function for a tree given an ODE function.
"""
function make_loglik(X        ::Array{Array{Float64,1},1},
                     abts1    ::Array{Float64,1},
                     abts2    ::Array{Float64,1},
                     trios    ::Array{Array{Int64,1},1},
                     ode_solve::Function,
                     λevent!  ::Function, 
                     rootll   ::Function,
                     k        ::Int64,
                     h        ::Int64,
                     ns       ::Int64,
                     ned      ::Int64)

  # preallocate vectors
  llik = Array{Float64,1}(undef,ns)
  w    = Array{Float64,1}(undef,ns)
  extp = Array{Float64,1}(undef, ns)

  function f(p::Array{Float64,1})

    @inbounds begin

      llxtra = 0.0

      # loop for integrating over internal branches
      for triad in trios

        pr, d1, d2 = triad::Array{Int64,1}

        ud1 = ode_solve(X[d1], p, abts2[d1], abts1[d1])::Array{Float64,1}

        check_negs(ud1, ns) && return -Inf

        ud2 = ode_solve(X[d2], p, abts2[d2], abts1[d2])::Array{Float64,1}

        check_negs(ud2, ns) && return -Inf

        # update likelihoods with speciation event
        λevent!(abts2[pr], llik, ud1, ud2, p)

        # loglik to sum for integration
        tosum = normbysum!(llik, ns)

        llxtra += log(tosum)

        # assign the remaining likelihoods &
        # assign extinction probabilities and 
        # check for extinction of `1.0`
        for i in Base.OneTo(ns)
          ud1[i+ns] > 1.0 && return -Inf
          X[pr][i+ns] = ud1[i+ns]
          X[pr][i]    = llik[i]
        end
      end

      # assign root likelihood in non log terms &
      # assign root extinction probabilities
      for i in Base.OneTo(ns)
        llik[i] = X[ned][i]
        extp[i] = X[ned][i+ns]
      end

      # estimate likelihood weights
      normbysum!(llik, w, ns)

      # combine root likelihoods
      ll = rootll(abts1[ned], llik, extp, w, p)

      return (log(ll) + llxtra)::Float64
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
                λts::Array{Float64,1},
                r   ::Array{Float64,1},
                ::Val{k}, 
                ::Val{h}, 
                ::Val{ny}, 
                ::Val{mdS}) where {k,h,ny,mdS}

Generated function for full tree likelihood at the root for normal likelihoods.
"""
@generated function rootll_full(t   ::Float64,
                                llik::Array{Float64,1},
                                extp::Array{Float64,1},
                                w   ::Array{Float64,1},
                                p   ::Array{Float64,1},
                                λts ::Array{Float64,1},
                                r   ::Array{Float64,1},
                                af! ::Function,
                                ::Val{k},
                                ::Val{h},
                                ::Val{ny},
                                ::Val{model})::Float64 where {k, h, ny, model}

  S = create_states(k, h)

  eqs = quote end
  popfirst!(eqs.args)

  if model[1]
   # last non β parameters
    bbase = h*(k+1) + 2*k*h + k*(k-1)*h + h*(h-1)

    # ncov
    ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

    # y per parameter
    yppar = isone(ny) ? 1 : div(ny,ncov)

    # add environmental function
    push!(eqs.args, :(af!(t, r)))

    isone(ny) && push!(eqs.args, :(r1 = r[1]))

    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(yppar)
        rex = isone(ny) ? :(r1) : :(r[$(yi+yppar*(i-1))])
        push!(coex.args, :(p[$(bbase + yi + yppar*((i-1) + k*(j-1)))] * $rex))
      end
      push!(eqs.args, :(λts[$(i + k*(j-1))] = p[$(i + (k+1)*(j-1))]*exp($coex)))
    end

    eq = Expr(:call, :+)
    for j in Base.OneTo(h)
      # for single areas
      for i in Base.OneTo(k)
        push!(eq.args,
          :(llik[$(i + (2^k-1)*(j-1))]*w[$(i + (2^k-1)*(j-1))] / 
            (λts[$(i + k*(j-1))]*(1.0-extp[$(i + (j-1)*(2^k-1))])^2)))
      end
      # for widespread
      for i in k+1:2^k-1
        wl = :(llik[$(i + (2^k-1)*(j-1))]*w[$(i + (2^k-1)*(j-1))] /(le))

        lams = Expr(:call, :+)
        # single area speciation
        for a in S[i + (2^k-1)*(j-1)].g
          push!(lams.args, :(
            λts[$(a + k*(j-1))] * 
           (1.0-extp[$(a + (j-1)*(2^k-1))]) * (1.0-extp[$(i + (j-1)*(2^k-1))])))
        end

        # for allopatric speciation
        for (la, ra) = vicsubsets(S[i + (2^k-1)*(j-1)].g)[1:(div(end,2))]
          push!(lams.args, 
            :(p[$((k+1)*(1+(j-1)))] *
             (1.0-extp[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*(j-1))]) * 
             (1.0-extp[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*(j-1))])))
        end

        # change l in wl
        wl.args[3] = lams

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
        wl = :(llik[$(i + (2^k-1)*(j-1))]*w[$(i + (2^k-1)*(j-1))] /(le))

        lams = Expr(:call, :+)
        # single area speciation
        for a in S[i + (2^k-1)*(j-1)].g
          push!(lams.args, :(
            p[$(a + (k+1)*(j-1))] * 
           (1.0-extp[$(a + (j-1)*(2^k-1))]) * (1.0-extp[$(i + (j-1)*(2^k-1))])))
        end

        # for allopatric speciation
        for (la, ra) = vicsubsets(S[i + (2^k-1)*(j-1)].g)[1:(div(end,2))]
          push!(lams.args, 
            :(p[$((k+1)*(1+(j-1)))] *
             (1.0-extp[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*(j-1))]) * 
             (1.0-extp[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*(j-1))])))
        end

        # change l in wl
        wl.args[3] = lams

        push!(eq.args, wl)
      end
    end
  end

  push!(eqs.args, :($eq))

  #Core.println(eqs)

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
                λts::Array{Float64,1},
                p   ::Array{Float64,1},
                r   ::Array{Float64,1},
                ::Val{k}, 
                ::Val{h}, 
                ::Val{ny}, 
                ::Val{S},
                ::Val{mdS})

Generated function for speciation event likelihoods for normal likelihood.
"""
@generated function λevent_full(t   ::Float64, 
                                llik::Array{Float64,1}, 
                                ud1 ::Array{Float64,1}, 
                                ud2 ::Array{Float64,1},
                                p   ::Array{Float64,1},
                                λts ::Array{Float64,1},
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
  if model[1] 
 
    # last non β parameters
    bbase = h*(k+1) + 2*k*h + k*(k-1)*h + h*(h-1)

    # ncov
    ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

    # y per parameter
    yppar = isone(ny) ? 1 : div(ny,ncov)

    # add environmental function
    push!(eqs.args, :(af!(t, r)))

    isone(ny) && push!(eqs.args, :(r1 = r[1]))

    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(yppar)
        rex = isone(ny) ? :(r1) : :(r[$(yi+yppar*(i-1))])
        push!(coex.args, :(p[$(bbase + yi + yppar*((i-1) + k*(j-1)))] * $rex))
      end
      push!(eqs.args, :(λts[$(i + k*(j-1))] = p[$(i + (k+1)*(j-1))]*exp($coex)))
    end

    # likelihood for individual areas states
    for j = Base.OneTo(h), i = Base.OneTo(k)
      push!(eqs.args, 
          :(llik[$(i+(2^k-1)*(j-1))] = 
            ud1[$(i+(2^k-1)*(j-1))] * ud2[$(i+(2^k-1)*(j-1))] * λts[$(i+k*(j-1))]))
    end

    # likelihood for widespread states
    for j = Base.OneTo(h), i = k+1:2^k-1
      s = S[i + (2^k-1)*(j-1)]

      # within region speciation
      ex = Expr(:call, :+)
      for a in s.g
        push!(ex.args, :((ud1[$(a + s.h*(2^k-1))]   * ud2[$(i + (2^k-1)*(j-1))] + 
                          ud1[$(i + (2^k-1)*(j-1))] * ud2[$(a + s.h*(2^k-1))]) *
                          0.5 * λts[$(a + s.h*(k))]))
      end

      # between region speciation
      for (la, ra) = vicsubsets(s.g)[1:div(end,2)]
        push!(ex.args, 
         :((ud1[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] + 
            ud1[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)]) *
            0.5 * p[$((k+1)*(1+s.h))]))
      end

      push!(eqs.args, :(llik[$(i + (2^k-1)*(j-1))] = $ex))
    end

  # *not* speciation model
  else

    # likelihood for individual areas states
    for j = Base.OneTo(h), i = Base.OneTo(k)
      push!(eqs.args, 
          :(llik[$(i+(2^k-1)*(j-1))] = 
            ud1[$(i+(2^k-1)*(j-1))] * ud2[$(i+(2^k-1)*(j-1))] * p[$(i+(k+1)*(j-1))]))
    end

    # likelihood for widespread states
    for j = Base.OneTo(h), i = k+1:2^k-1
      s = S[i + (2^k-1)*(j-1)]

      # within region speciation
      ex = Expr(:call, :+)
      for a in s.g
        push!(ex.args, :((ud1[$(a + s.h*(k+1))]*ud2[$(i + (2^k-1)*(j-1))] + 
                          ud1[$(i + (2^k-1)*(j-1))]*ud2[$(a + s.h*(k+1))]) *
                          0.5 * p[$(a + s.h*(k+1))]))
      end

      # between region speciation
      for (la, ra) = vicsubsets(s.g)[1:div(end,2)]
        push!(ex.args, 
         :((ud1[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] + 
            ud1[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)]) *
            0.5 * p[$((k+1)*(1+s.h))]))
      end

      push!(eqs.args, :(llik[$(i + (2^k-1)*(j-1))] = $ex))
    end
  end

  #Core.println(eqs)


  return quote 
    @inbounds begin
      $eqs
    end
    return nothing
  end
end






"""
    λevent_full(t   ::Float64, 
                llik::Array{Array{Float64,1},1}, 
                ud1 ::Array{Float64,1}, 
                ud2 ::Array{Float64,1},
                p   ::Array{Float64,1},
                pr  ::Int64,
                λts ::Array{Float64,1},
                r   ::Array{Float64,1},
                af! ::Function,
                ::Val{k}, 
                ::Val{h}, 
                ::Val{ny}, 
                ::Val{model}) where {k, h, ny, model}

Generated function for speciation event likelihoods for flow algorithm.
"""
@generated function λevent_full(t   ::Float64, 
                                llik::Array{Array{Float64,1},1}, 
                                ud1 ::Array{Float64,1}, 
                                ud2 ::Array{Float64,1},
                                p   ::Array{Float64,1},
                                pr  ::Int64,
                                λts ::Array{Float64,1},
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
  if model[1] 
 
    # last non β parameters
    bbase = h*(k+1) + 2*k*h + k*(k-1)*h + h*(h-1)

    # ncov
    ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

    # y per parameter
    yppar = isone(ny) ? 1 : div(ny,ncov)

    # add environmental function
    push!(eqs.args, :(af!(t, r)))

    isone(ny) && push!(eqs.args, :(r1 = r[1]))

    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(yppar)
        rex = isone(ny) ? :(r1) : :(r[$(yi+yppar*(i-1))])
        push!(coex.args, :(p[$(bbase + yi + yppar*((i-1) + k*(j-1)))] * $rex))
      end
      push!(eqs.args, :(λts[$(i + k*(j-1))] = p[$(i + (k+1)*(j-1))]*exp($coex)))
    end

    # likelihood for individual areas states
    for j = Base.OneTo(h), i = Base.OneTo(k)
      push!(eqs.args, 
          :(llik[pr][$(i+(2^k-1)*(j-1))] = 
            ud1[$(i+(2^k-1)*(j-1))] * ud2[$(i+(2^k-1)*(j-1))] * λts[$(i+k*(j-1))]))
    end

    # likelihood for widespread states
    for j = Base.OneTo(h), i = k+1:2^k-1
      s = S[i + (2^k-1)*(j-1)]

      # within region speciation
      ex = Expr(:call, :+)
      for a in s.g
        push!(ex.args, :((ud1[$(a + s.h*(2^k-1))]   * ud2[$(i + (2^k-1)*(j-1))] + 
                          ud1[$(i + (2^k-1)*(j-1))] * ud2[$(a + s.h*(2^k-1))]) *
                          0.5 * λts[$(a + s.h*(k))]))
      end

      # between region speciation
      for (la, ra) = vicsubsets(s.g)[1:div(end,2)]
        push!(ex.args, 
         :((ud1[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] + 
            ud1[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)]) *
            0.5 * p[$((k+1)*(1+s.h))]))
      end

      push!(eqs.args, :(llik[pr][$(i + (2^k-1)*(j-1))] = $ex))
    end

  # *not* speciation model
  else

    # likelihood for individual areas states
    for j = Base.OneTo(h), i = Base.OneTo(k)
      push!(eqs.args, 
          :(llik[pr][$(i+(2^k-1)*(j-1))] = 
            ud1[$(i+(2^k-1)*(j-1))] * ud2[$(i+(2^k-1)*(j-1))] * p[$(i+(k+1)*(j-1))]))
    end

    # likelihood for widespread states
    for j = Base.OneTo(h), i = k+1:2^k-1
      s = S[i + (2^k-1)*(j-1)]

      # within region speciation
      ex = Expr(:call, :+)
      for a in s.g
        push!(ex.args, :((ud1[$(a + s.h*(k+1))]*ud2[$(i + (2^k-1)*(j-1))] + 
                          ud1[$(i + (2^k-1)*(j-1))]*ud2[$(a + s.h*(k+1))]) *
                          0.5 * p[$(a + s.h*(k+1))]))
      end

      # between region speciation
      for (la, ra) = vicsubsets(s.g)[1:div(end,2)]
        push!(ex.args, 
         :((ud1[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] + 
            ud1[$(findfirst(x -> isequal(la,x.g), S) + (2^k-1)*s.h)] *
            ud2[$(findfirst(x -> isequal(ra,x.g), S) + (2^k-1)*s.h)]) *
            0.5 * p[$((k+1)*(1+s.h))]))
      end

      push!(eqs.args, :(llik[pr][$(i + (2^k-1)*(j-1))] = $ex))
    end
  end

  #Core.println(eqs)


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
    normbysum!(v::Array{Float64,1}, ns::Int64)

Return weights by the sum of all elements of the array to sum 1 and overwrites. 
"""
function normbysum!(v::Array{Float64,1}, ns::Int64)
  s = sum(v)
  @inbounds begin
    @simd for i in Base.OneTo(ns)
      v[i] /= s
    end
  end
  return s
end






"""
    make_lpf(pupd   ::Array{Int64,1},
             phid   ::Array{Int64,1},
             λpriors::Float64,
             μpriors::Float64,
             gpriors::Float64,
             lpriors::Float64,
             qpriors::Float64,
             βpriors::Tuple{Float64,Float64},
             hpriors::Float64,
             k      ::Int64,
             h      ::Int64,
             model  ::Tuple{Bool,Bool,Bool})

Make log-prior function.
"""
function make_lpf(pupd   ::Array{Int64,1},
                  phid   ::Array{Int64,1},
                  λpriors::Float64,
                  μpriors::Float64,
                  gpriors::Float64,
                  lpriors::Float64,
                  qpriors::Float64,
                  βpriors::Tuple{Float64,Float64},
                  hpriors::Float64,
                  k      ::Int64,
                  h      ::Int64,
                  ny     ::Int64,
                  model  ::Tuple{Bool,Bool,Bool})

  λupds = intersect(1:(h*(k+1)), pupd)
  μupds = intersect((h*(k+1)+1):(h*(k+1)+k*h), pupd)
  gupds = intersect((h*(k+1)+k*h+1):(h*(k+1)+k*h+k*(k-1)*h), pupd)
  lupds = intersect((h*(k+1)+k*h+k*(k-1)*h+1):(h*(k+1)+2k*h+k*(k-1)*h), pupd)
  qupds = intersect((h*(k+1)+2k*h+k*(k-1)*h+1):(h*(k+1)+2k*h+k*(k-1)*h+h*(h-1)), 
            pupd)
  bbase = (h*(k+1)+2k*h+k*(k-1)*h+h*(h-1))
  yppar = ny == 1 ? 1 : div(ny,
    model[1]*k + 
    model[2]*k + 
    model[3]*k*(k-1))

  ncov  = model[1]*k*yppar + model[2]*k*yppar + model[3]*k*(k-1)*yppar
  βupds = intersect((bbase+1):(bbase+ncov), pupd)

  βp_m, βp_v = βpriors

  function f(p ::Array{Float64,1}, 
             fp::Array{Float64,1})
    lq = 0.0

    # speciation priors
    for i in λupds
      lq += logdexp(p[i], λpriors)
    end

    # global extinction priors
    for i in μupds
      lq += logdexp(p[i], μpriors)
    end

    # area colonization priors
    for i in gupds
      lq += logdexp(p[i], gpriors)
    end

    # area loss priors
    for i in lupds
      lq += logdexp(p[i], lpriors)
    end

    # hidden states transition
    for i in qupds
      lq += logdexp(p[i], qpriors)
    end

    # betas
    for i in βupds
      lq += logdnorm(p[i], βp_m, βp_v)
    end
  
    # hidden states factors
    for i in phid
      lq += logdexp(fp[i], hpriors)
    end

    return lq
  end

  return f
end





"""
    make_assign_hidfacs(::Val{k},
                        ::Val{h},
                        ::Val{ny},
                        ::Val{model}) where {k, h, ny, model})

Generated function to assign factors to parameters given 
factor+parameter `fp` vector.
"""
function make_assign_hidfacs(::Val{k},
                             ::Val{h},
                             ::Val{ny},
                             ::Val{model}) where {k, h, ny, model}

  assign_hidfacs! = (p::Array{Float64,1}, fp::Array{Float64,1}) -> 
    assign_hidfacs_full(p, fp, Val(k), Val(h), Val(ny), Val(model))

  return assign_hidfacs!
end





"""
    assign_hidfacs(p ::Array{Float64,1},
                   fp::Array{Float64,1},
                   ::Val{k},
                   ::Val{h},
                   ::Val{ny},
                   ::Val{model}) where {k, h, ny, model}

Generated function to assign factors to parameters given 
factor+parameter `fp` vector.
"""
@generated function assign_hidfacs_full(p ::Array{Float64,1},
                                        fp::Array{Float64,1},
                                        ::Val{k},
                                        ::Val{h},
                                        ::Val{ny},
                                        ::Val{model}) where {k, h, ny, model}

  # number of covariates
  yppar = ny == 1 ? 1 : div(ny,
    model[1]*k + 
    model[2]*k + 
    model[3]*k*(k-1))

  # starting indices for models 2 and 3
  m2s = model[1]*k*h*yppar
  m3s = m2s + model[2]*k*h*yppar

  ex = quote end
  pop!(ex.args)

  for j in 1:(h-1)
    for i in 1:k

      # speciation
      push!(ex.args, :(p[$((k+1)*j + i)] = 
        p[$((k+1)*(j-1) + i)] + fp[$((k+1)*j + i)]))

      # extinction
      s = (k+1)*h 
      push!(ex.args, 
        :(p[$(s + k*j + i)] = 
            p[$(s + k*(j-1) + i)] + fp[$(s + k*j + i)]))

      # gain
      s = (2k+1)*h
      for a in 1:(k-1)
        push!(ex.args, 
          :(p[$(s + k*(k-1)*j + a + (k-1)*(i-1))] = 
              p[$(s + k*(k-1)*(j-1) + a + (k-1)*(i-1))] + 
              fp[$(s + k*(k-1)*j + a + (k-1)*(i-1))]))
      end

      # loss 
      s = (2k+1 + k*(k-1))*h
      push!(ex.args, 
        :(p[$(s + k*j + i)] = 
            p[$(s + k*(j-1) + i)] + fp[$(s + k*j + i)]))

      # betas
      s = (3k+1+k*(k-1))*h + h*(h-1)
      if model[1]
        for l = Base.OneTo(yppar)
          push!(ex.args, 
            :(p[$(s + l + yppar*(i-1) + k*j*yppar)] = 
              p[$(s + l + yppar*(i-1) + k*(j-1)*yppar)] + 
              fp[$(s + l + yppar*(i-1) + k*j*yppar)]))
        end
      end
      if model[2]
        for l = Base.OneTo(yppar)
          push!(ex.args,
            :(p[$(s + m2s + l + yppar*(i-1) + k*j*yppar)] = 
              p[$(s + m2s + l + yppar*(i-1) + k*(j-1)*yppar)] + 
              fp[$(s + m2s + l + yppar*(i-1) + k*j*yppar)]))
        end
      end
      if model[3]
        for a = 1:(k-1), l = Base.OneTo(yppar)
          push!(ex.args, 
            :(p[$(s + m3s + l + (a-1)*yypar + yppar*(i-1) + k*j*yppar)] = 
              p[$(s + m3s + l + (a-1)*yypar + yppar*(i-1) + k*(j-1)*yppar)] + 
              fp[$(s + m3s + l + (a-1)*yypar + yppar*(i-1) + k*j*yppar)]))
        end
      end
    end
    # between-region speciation
    push!(ex.args, :(p[$((k+1)*j + (k+1))] = 
      (p[$((k+1)*(j-1) + (k+1))] + fp[$((k+1)*j + (k+1))])))
  end

  #Core.println(ex)

  return quote
    @inbounds begin
      $ex
    end
  end
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





"""
    make_λevent(h    ::Int64, 
                k    ::Int64, 
                ny   ::Int64, 
                flow ::Bool,
                model::NTuple{3,Bool},
                af!  ::Function)

Make function for λevent.
"""
function make_λevent(h    ::Int64, 
                     k    ::Int64, 
                     ny   ::Int64, 
                     flow ::Bool,
                     model::NTuple{3,Bool},
                     af!  ::Function)

  λts = Array{Float64,1}(undef,h*k)
  r   = Array{Float64,1}(undef,ny)

  # make speciation events and closure
  if flow
    λevent! = (t   ::Float64, 
               llik::Array{Array{Float64,1},1},
               ud1 ::Array{Float64,1},
               ud2 ::Array{Float64,1},
               p   ::Array{Float64,1},
               pr  ::Int64) ->
    begin
      λevent_full(t, llik, ud1, ud2, p, pr, λts, r, af!,
        Val(k), Val(h), Val(ny), Val(model))
      return nothing
    end

  else
    λevent! = (t   ::Float64, 
               llik::Array{Float64,1},
               ud1 ::Array{Float64,1},
               ud2 ::Array{Float64,1},
               p   ::Array{Float64,1}) ->
    begin
      λevent_full(t, llik, ud1, ud2, p, λts, r, af!,
        Val(k), Val(h), Val(ny), Val(model))
      return nothing
    end
  end

  return λevent!
end





"""
    make_rootll(h    ::Int64, 
                k    ::Int64, 
                ny   ::Int64, 
                model::NTuple{3,Bool},
                af!  ::Function)

Make root conditioning likelihood function.
"""
function make_rootll(h    ::Int64, 
                     k    ::Int64, 
                     ny   ::Int64, 
                     model::NTuple{3,Bool},
                     af!  ::Function)

  λts = Array{Float64,1}(undef,h*k)
  r   = Array{Float64,1}(undef,ny)

  rootll = (t   ::Float64,
            llik::Array{Float64,1},
            extp::Array{Float64,1},
            w   ::Array{Float64,1},
            p   ::Array{Float64,1}) -> 
    begin
      rootll_full(t, llik, extp, w, p, λts, r, af!,
        Val(k), Val(h), Val(ny), Val(model))::Float64
    end

  return rootll
end










