#=

ESSE marginal ancestral estate probabilities

Ignacio Quintero Mächler

t(-_-t)

September 26 2017

=#




"""
    sample_node_ps(R       ::Array{Float64,2},
                   spf     ::Function,
                   nsamples::Int64, 
                   ns      ::Int64, 
                   ned     ::Int64)

Run node marginal probabilities estimation.
"""
function sample_node_ps(R       ::Array{Float64,2},
                        A       ::Array{Array{Float64,1},1},
                        spf     ::Function,
                        nsamples::Int64, 
                        ns      ::Int64, 
                        ned     ::Int64)

  nlit = size(R, 1)
  npar = size(R, 2) - 2

  nsamples = nsamples > nlit ? nlit : nsamples
  lever = ceil(Int64, nlit/nsamples)

  #prellocate results
  S = Array{Float64,2}(undef, nsamples, ns*ned+2)

  pv = zeros(npar)

  ss = 1:lever:nlit

  for s in Base.OneTo(lastindex(ss))
    spf(R[ss[s],3:(npar+2)])
    S[s,1] = R[ss[s],1]
    S[s,2] = R[ss[s],2]

    # copy to fix
    @simd for i in Base.OneTo(ned)
      i_i = (i + (i-1)*(ns-1)) + 2
      i_f = (i + (i-1)*(ns-1) + (ns-1)) + 2
      S[s,i_i:i_f] = copy(A[i])
    end
  end

  return S
end




"""
    make_state_posteriors(llf  ::Function,
                          lpf  ::Function, 
                          llfnj::Function,
                          X    ::Array{Array{Float64,1},1},
                          U    ::Array{Array{Float64,1},1},
                          M    ::Array{Array{Float64,1},1},
                          ns   ::Int64,
                          ned  ::Int64,
                          k    ::Int64,
                          h    ::Int64)

Make function to estimate node states marginal posterior probabilities.
"""
function make_state_posteriors(llf  ::Function,
                               lpf  ::Function, 
                               llfnj::Function,
                               X    ::Array{Array{Float64,1},1},
                               U    ::Array{Array{Float64,1},1},
                               A    ::Array{Array{Float64,1},1},
                               ns   ::Int64,
                               ned  ::Int64,
                               k    ::Int64,
                               h    ::Int64)

  Xtl = deepcopy(X)::Array{Array{Float64,1},1}

  # estimate fp for prior
  hdic = Dict{Int64, Int64}()
  if isone(k)
    for i in Base.OneTo(h-1)
      for j in Base.OneTo(k)
        push!(hdic, j + k*(i-1) => j + k*i)
      end
    end
  else
    for i in Base.OneTo(h-1)
      for j in Base.OneTo(k+1)
        push!(hdic, j + (k+1)*(i-1) => j + (k+1)*i)
      end
    end
  end

  fp = isone(k) ? zeros(k*h) : zeros((k+1)*h)

  function f(p::Array{Float64,1})

    @inbounds begin
      # estimate total likelihood
      tll  = llf(p)

      # estimate fp from p
      for (k,v) in hdic
        fp[v] = p[v] - p[k]
      end

      # estimate prior
      tlp  = lpf(p, fp)
      # estimate posterior
      tlpp =  tlp + tll

      # copy to fix
      @simd for i in Base.OneTo(ned)
        unsafe_copyto!(Xtl[i], 1, X[i], 1, ns)
      end

      for n in Base.OneTo(ned)
        An = A[n]
        for i in Base.OneTo(ns)
          An[i] = exp(llfnj(p, n, i, Xtl, U))
          if isnan(An[i])
            An[i] = 0.0
          end
        end
        sA = sum(An)
        for i in Base.OneTo(ns)
          An[i] /= sA
        end
      end
    end
  end

  return f
end




"""
    nodetoroot_triads(ed   ::Array{Int64,2},
                      trios::Array{Array{Int64,1},1})

Create dictionary of trios for all nodes such that it only goes from
node to tip.
"""
function nodetoroot_triads(ed   ::Array{Int64,2},
                           trios::Array{Array{Int64,1},1})

  ed2 = ed[:,2]

  tdic = Array{Tuple{Array{Int64,1},Int64},1}[]

  for i in Base.OneTo(size(ed,1))

    # make node-to-root path
    ndic = Int64[] 

    n = i
    while true
      push!(ndic, n)
      pr, da = ed[n,:]
      n = findfirst(isequal(pr), ed2)
      if isnothing(n)
        break
      end
    end

    # determine trios where the nodes in ndic are
    tdicn = Tuple{Array{Int64,1},Int64}[]
    for n in ndic, t in trios
      @views pos = findfirst(isequal(n), t[2:3])
      if !isnothing(pos)
        push!(tdicn, (t, pos))
      end
    end

    push!(tdic, tdicn)
  end

  return tdic
end




"""
    make_loglik_nj(X        ::Array{Array{Float64,1},1},
                   abts1    ::Array{Float64,1},
                   abts2    ::Array{Float64,1},
                   trios    ::Array{Array{Int64,1},1},
                   int      ::DiffEqBase.DEIntegrator,
                   λevent!  ::Function, 
                   rootll   ::Function,
                   ns       ::Int64,
                   ned      ::Int64)

Make likelihood function for tree and tip states following an ODE function 
where node `n` has state `j`.
"""
function make_loglik_nj(X         ::Array{Array{Float64,1},1},
                        tdic      ::Array{Array{Tuple{Array{Int64,1},Int64},1},1},
                        abts1     ::Array{Float64,1},
                        abts2     ::Array{Float64,1},
                        trios     ::Array{Array{Int64,1},1},
                        int       ::DiffEqBase.DEIntegrator,
                        λevent!   ::Function, 
                        rootll_nj!::Function,
                        ns        ::Int64,
                        ned       ::Int64)

  # preallocate vectors
  llik    = Array{Float64,1}(undef, ns)
  llik_tl = Array{Float64,1}(undef, ns)
  w       = Array{Float64,1}(undef, ns)
  w_tl    = Array{Float64,1}(undef, ns)
  extp    = Array{Float64,1}(undef, ns)

  function f(p   ::Array{Float64,1}, 
             n   ::Int64, 
             j   ::Int64, 
             Xtl::Array{Array{Float64,1},1},
             U   ::Array{Array{Float64,1},1})

    @inbounds begin

      int.p = p

      llxtra = 0.0

      tdicn = tdic[n]

      # copy from fix
      @simd for i in Base.OneTo(ned)
        unsafe_copyto!(X[i], 1, Xtl[i], 1, ns)
      end

      # make node likelihood `0.0` to all other states
      for i in Base.OneTo(ns)
        i == j && continue
        X[n][i] = 0.0
      end

      # loop for integrating over internal branches from node-to-root
      for (triad, pos) in tdicn

        pr, d1, d2 = triad

        if isone(pos)
          ud1 = solvef(int, X[d1], abts2[d1], abts1[d1])::Array{Float64,1}
          check_negs(ud1, ns) && return -Inf
          ud2 = U[d2]
        else
          ud1 = U[d1]
          ud2 = solvef(int, X[d2], abts2[d2], abts1[d2])::Array{Float64,1}
          check_negs(ud2, ns) && return -Inf
        end

        # update marginal likelihoods with speciation event
        λevent!(abts2[pr], llik, ud1, ud2, p)

        # update total likelihoods with speciation event
        λevent!(abts2[pr], llik_tl, U[d1], U[d2], p)

        # loglik for marginal likelihood of n
        tosum = normbysum!(llik, ns)

        # loglik for total likelihood
        tosum_tl = normbysum!(llik_tl, ns)

        # loglik to sum for fixed node
        llxtra += log(tosum) - log(tosum_tl)

        # assign the remaining likelihoods &
        # assign extinction probabilities and 
        # check for extinction of `1.0`
        Xpr = X[pr]
        for i in Base.OneTo(ns)
          ud1[i+ns] >= 1.0 && return -Inf
          Xpr[i+ns] = ud1[i+ns]
          Xpr[i]    = llik[i]
        end
      end

      # assign root likelihood in non log terms &
      # assign root extinction probabilities
      Xned    = X[ned]
      Xned_tl = Xtl[ned]
      for i in Base.OneTo(ns)
        Xned[i+ns] >= 1.0 && return -Inf
        llik[i]    = Xned[i]
        extp[i]    = Xned[i+ns]
        llik_tl[i] = Xned_tl[i]
      end

      if n == ned
        # estimate likelihood weights
        normbysum!(llik_tl, w_tl, ns)

        # combine root likelihoods for marginal likelihood
        rootll_nj!(abts1[ned], llik_tl, extp, w_tl, p)

        # for root, weights and extinction have already been estimated
        return (log(llik_tl[j]) - log(sum(llik_tl)) + llxtra)
      else
        # estimate likelihood weights
        normbysum!(llik,    w,    ns)
        normbysum!(llik_tl, w_tl, ns)

        # combine root likelihoods for marginal likelihood
        rootll_nj!(abts1[ned], llik,    extp, w,    p)
        rootll_nj!(abts1[ned], llik_tl, extp, w_tl, p)

        return (log(sum(llik)) - log(sum(llik_tl)) + llxtra)::Float64
      end
    end
  end

  return f
end





"""
    rootll_full_nj(t   ::Float64,
                   llik::Array{Float64,1},
                   extp::Array{Float64,1},
                   w   ::Array{Float64,1},
                   p   ::Array{Float64,1},
                   λts ::Array{Float64,1},
                   r   ::Array{Float64,1},
                   af! ::Function,
                   ::Type{Val{k}},
                   ::Type{Val{h}},
                   ::Type{Val{ny}},
                   ::Type{Val{model}})::Float64 where {k, h, ny, model}

Generated function for full tree likelihood at the root for 
**pruning** likelihoods that returns a vector for each state.
"""
@generated function rootll_full_nj(t   ::Float64,
                                   llik::Array{Float64,1},
                                   extp::Array{Float64,1},
                                   w   ::Array{Float64,1},
                                   p   ::Array{Float64,1},
                                   λts ::Array{Float64,1},
                                   r   ::Array{Float64,1},
                                   af! ::Function,
                                   ::Type{Val{k}},
                                   ::Type{Val{h}},
                                   ::Type{Val{ny}},
                                   ::Type{Val{model}})::Float64 where {k, h, ny, model}

  S = create_states(k, h)

  eqs = quote end
  popfirst!(eqs.args)

  if model[1]

    # last non β parameters
    bbase = isone(k) ? 
      (2h + h*(h-1)) : (h*(k+1) + 2*k*h + k*(k-1)*h + h*(h-1))

    # ncov
    ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

    # y per parameter
    yppar = isone(ny) ? 1 : ceil(Int64,ny/ncov)

    # add environmental function
    push!(eqs.args, :(af!(t, r)))

    isone(ny) && push!(eqs.args, :(r1 = r[1]))

    for j = Base.OneTo(h), i = Base.OneTo(k)
      coex = Expr(:call, :+)
      for yi in Base.OneTo(yppar)
        rex = isone(ny) ? :(r1) : :(r[$(yi+yppar*(i-1))])
        push!(coex.args, :(p[$(bbase + yi + yppar*((i-1) + k*(j-1)))] * $rex))
      end
      if isone(k)
        push!(eqs.args, :(λts[$(i + (j-1))] = p[$(i + (j-1))]*exp($coex)))
      else
        push!(eqs.args, :(λts[$(i + k*(j-1))] = p[$(i + (k+1)*(j-1))]*exp($coex)))
      end
    end

    eq = quote end
    popfirst!(eq.args)
    for j in Base.OneTo(h)
      # for single areas
      for i in Base.OneTo(k)
        push!(eq.args,
          :(llik[$(i + (2^k-1)*(j-1))] *= 
            w[$(i + (2^k-1)*(j-1))] / 
            (λts[$(i + k*(j-1))]*(1.0-extp[$(i + (j-1)*(2^k-1))])^2)))
      end
      # for widespread
      for i in k+1:2^k-1
        wl = :(llik[$(i + (2^k-1)*(j-1))] *= 
               w[$(i + (2^k-1)*(j-1))] /(le))

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

        # change le in wl
        wl.args[2].args[3] = lams

        push!(eq.args, wl)
      end
    end

  # not speciation model
  else

    eq = quote end
    popfirst!(eq.args)
    for j in Base.OneTo(h)
      # for single areas
      for i in Base.OneTo(k)
        if isone(k)
          push!(eq.args,
            :(llik[$(i + (2^k-1)*(j-1))] *= 
              w[$(i + (2^k-1)*(j-1))] / 
              (p[$(i + (j-1))]*(1.0-extp[$(i + (j-1)*(2^k-1))])^2)))
        else
          push!(eq.args,
            :(llik[$(i + (2^k-1)*(j-1))] *= 
              w[$(i + (2^k-1)*(j-1))] / 
              (p[$(i + (k+1)*(j-1))]*(1.0-extp[$(i + (j-1)*(2^k-1))])^2)))
        end
      end
      # for widespread
      for i in k+1:2^k-1
        wl = :(llik[$(i + (2^k-1)*(j-1))] *= 
               w[$(i + (2^k-1)*(j-1))] /(le))

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
        wl.args[2].args[3] = lams

        push!(eq.args, wl)
      end
    end
  end

  push!(eqs.args, :($eq))

  @debug eqs

  return quote
    @inbounds begin
      $eqs
    end
  end
end






"""
    make_rootll_nj(::Type{Val{h}}, 
                   ::Type{Val{k}}, 
                   ::Type{Val{ny}}, 
                   ::Type{Val{model}},
                   af!  ::Function) where {h, k, ny, model}

Make root conditioning likelihood function.
"""
function make_rootll_nj(::Type{Val{h}}, 
                        ::Type{Val{k}}, 
                        ::Type{Val{ny}}, 
                        ::Type{Val{model}},
                        af!  ::Function) where {h, k, ny, model}

  λts = Array{Float64,1}(undef,h*k)
  r   = Array{Float64,1}(undef,ny)

  rootll = let r = r, λts = λts, af! = af!
      (t   ::Float64,
              llik::Array{Float64,1},
              extp::Array{Float64,1},
              w   ::Array{Float64,1},
              p   ::Array{Float64,1}) -> 
      begin
        rootll_full_nj(t, llik, extp, w, p, λts, r, af!,
          Val{k}, Val{h}, Val{ny}, Val{model})::Float64
      end
  end
  
  return rootll
end



