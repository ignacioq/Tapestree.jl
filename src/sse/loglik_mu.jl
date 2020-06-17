#=

Log-likelihood for MuSSE type ODEs, integrating along the tree

Ignacio Quintero Mächler

t(-_-t)

Created 20 09 2017

=#




"""
    make_llf(tip_val::Dict{Int64,Array{Float64,1}},
             ed     ::Array{Int64,2},
             el     ::Array{Float64,1},
             ode_fun,
             af     ::Function,
             p      ::Array{Float64,1};
             sbrlen ::Float64 = 5.0)

Make likelihood function for a tree given an ODE function.
"""
function make_llf(tip_val::Dict{Int64,Array{Float64,1}},
                  ed     ::Array{Int64,2},
                  el     ::Array{Float64,1},
                  ode_fun,
                  af     ::Function,
                  p      ::Array{Float64,1},
                  md     ::Bool, 
                  ws     ::Bool;
                  sbrlen ::Float64 = 5.0)

  k    = length(tip_val[1])::Int64 
  k2   = 2k::Int64
  ke   = (k+1:2k)::UnitRange{Int64}
  ntip = length(tip_val)::Int64
  bk   = k + (k*k)
  r    = Array{Float64}(undef, k)

  # add long root branch
  ed = cat(ed, [2*ntip ntip + 1], dims = 1)
  push!(el, sbrlen)

  # get absolute times of branches as related to z(t)
  elrt = abs_time_branches(el, ed, ntip)

  ne    = size(ed,1)::Int64
  child = ed[:,2]::Array{Int64,1}
  wtp   = findall(child .<= ntip)::Array{Int64,1}

  # make trios
  trios = maketriads(ed, rev = true)::Array{Array{Int64,1},1}

  # preallocate tip likelihoods
  led = Array{Float64,1}[]
  for i in Base.OneTo(ne)
    push!(led, zeros(Float64, k2))
  end

  # assign states to terminal branches
  for wi in wtp 
     led[wi][1:k] = tip_val[child[wi]]
  end

  # partial log-likelihoods
  llik = Array{Float64}(undef,k)

  # root likelihood
  rll = Array{Float64}(undef,k2)

  # log lambdas 
  lλs = Array{Float64}(undef,k)

  # make ode solver
  ode_solve = make_solver(ode_fun, p)

  λevent! = make_λevent(af, k, ws, md)

  function f(p::Array{Float64})

    @inbounds begin

      ll     = 0.0::Float64
      llxtra = 0.0::Float64

      for i in Base.OneTo(k)
        lλs[i] = log(p[i])::Float64
      end

      # loop for integrating over internal branches
      for triad in trios

        pr, d1, d2 = triad::Array{Int64,1}

        ud1 = ode_solve(led[d1], p, elrt[d1,2], elrt[d1,1])::Array{Float64,1}

        if check_negs(ud1, k)
          return -Inf
        end

        ud2 = ode_solve(led[d2], p, elrt[d2,2], elrt[d2,1])::Array{Float64,1}

        if check_negs(ud2, k)
          return -Inf
        end

        # update likelihoods with speciation event
        λevent!(elrt[pr,2], llik, ud1, ud2, lλs, p)

        # loglik to sum for integration
        tosum   = minimum(llik)
        llxtra -= tosum

        # assign the remaining likelihoods
        for i in Base.OneTo(k)
          led[pr][i] = exp(llik[i] - tosum)
        end

        # assign extinction probabilities
        @views led[pr][ke] = ud1[ke]
      end

      # stem branch
      rll = ode_solve(led[ne], p, elrt[ne,2], elrt[ne,1])::Array{Float64,1}

      if check_negs(rll, k)
          return -Inf
      end

      for i in Base.OneTo(k)
        ll += rll[i]::Float64
      end

      return (log(ll) - llxtra)::Float64
    end
  end

end





"""
    make_λevent(af, k::Int64, ws::Bool, md::Bool)

Make function for speciation event likelihoods
"""
function make_λevent(af, k::Int64, ws::Bool, md::Bool)

  r  = Array{Float64}(undef, k)
  bk = k + (k*k)

  # speciation ESSE model with MULTIvariate z(t)
  function f1(t   ::Float64, 
              llik::Array{Float64,1}, 
              ud1 ::Array{Float64,1}, 
              ud2 ::Array{Float64,1}, 
              lλs ::Array{Float64,1}, 
              p   ::Array{Float64,1}) 

    af(t, r)

    # estimate ll to sum,
    for i in Base.OneTo(k)
      llik[i] = log(ud1[i] * ud2[i]) + lλs[i] + p[bk+i] * r[i]
    end

    return nothing
  end

  # speciation ESSE model with UNIvariate z(t)
  function f2(t   ::Float64, 
              llik::Array{Float64,1}, 
              ud1 ::Array{Float64,1}, 
              ud2 ::Array{Float64,1}, 
              lλs ::Array{Float64,1}, 
              p   ::Array{Float64,1}) 

    r = af(t)

    # estimate ll to sum,
    for i in Base.OneTo(k)
      llik[i] = log(ud1[i] * ud2[i]) + lλs[i] + p[bk+i] * r
    end

    return nothing
  end

  # other ESSE model
  function f3(t   ::Float64, 
              llik::Array{Float64,1}, 
              ud1 ::Array{Float64,1}, 
              ud2 ::Array{Float64,1}, 
              lλs ::Array{Float64,1}, 
              p   ::Array{Float64,1}) 

    # estimate ll to sum,
    for i in Base.OneTo(k)
      llik[i] = log(ud1[i] * ud2[i]) + lλs[i]
    end

    return nothing
  end

  if ws
    if md 
      return f1
    else
      return f2
    end
  else
    return f3
  end 
end





"""
    make_lpf(λpriors::Float64, μpriors::Float64, qpriors::Float64, βpriors, k::Int64)

Make log-prior function.
"""
function make_lpf(λpriors::Float64,
                  μpriors::Float64,
                  qpriors::Float64,
                  βpriors::Tuple{Float64,Float64},
                  k      ::Int64,
                  T      ::Bool)

  l2::UnitRange{Int64} = (k+1):2k      
  l3::UnitRange{Int64} = (2k+1):(k+k*k)
  lp::Float64          = 0.0 

  if T
    l4 = ((k+k*k+1):(2k+k*k))::UnitRange{Int64} 
  else
    l4 = ((k+k*k+1):(2k*k))::UnitRange{Int64} 
  end

  function f(p::Array{Float64,1})
    lp = 0.0::Float64

    for i in Base.OneTo(k)
      lp += logdexp(p[i], λpriors)::Float64
    end

    for i in l2
      lp += logdexp(p[i], μpriors)::Float64
    end

    for i in l3
      lp += logdexp(p[i], qpriors)::Float64
    end

    for i in l4
      lp += logdnorm(p[i], βpriors[1], βpriors[2])::Float64
    end

    return lp::Float64
  end

  return f
end






"""
    make_lpf(λpriors::Float64, μpriors::Float64, qpriors::Float64, k::Int64)

Make log-prior function.
"""
function make_lpf(λpriors::Float64,
                  μpriors::Float64,
                  qpriors::Float64,
                  k      ::Int64)

  l2::UnitRange{Int64} = (k+1):2k      
  l3::UnitRange{Int64} = (2k+1):(k+k*k)
  lp::Float64          = 0.0 

  function f(p::Array{Float64,1})
    lp = 0.0::Float64

    for i in Base.OneTo(k)
      lp += logdexp(p[i], λpriors)::Float64
    end

    for i in l2
      lp += logdexp(p[i], μpriors)::Float64
    end

    for i in l3
      lp += logdexp(p[i], qpriors)::Float64
    end

    return lp::Float64
  end

  return f
end











