#=

constant birth-death simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#





"""
    sim_cbd(t::Float64, λ::Float64, μ::Float64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ` conditioned on extinction before time `t`.
"""
function sim_cbd_ext(t::Float64, λ::Float64, μ::Float64)

  tw, λorμ = cbd_wait_ext(λ, μ, t)

  if tw > t
    return sTbd(t, false, false)
  end

  if λorμ
    return sTbd(sim_cbd_ext(t - tw, λ, μ), 
                sim_cbd_ext(t - tw, λ, μ), tw, false, false)
  else
    return sTbd(tw, true, false)
  end
end




"""
    cbd_wait_ext(λ::Float64, μ::Float64, t::Float64)

Sample a per-lineage waiting time for constant birth-death species
with speciation rate `λ` and extinction rate `μ` conditioned on 
extinction before time `t`.
"""
function cbd_wait_ext(λ::Float64, μ::Float64, t::Float64)
  λt = rexp(λ)
  μt = rexpt(μ, t)

  if λt < μt
    return λt, true
  else
    return μt, false
  end
end




"""
    sim_cbd(t::Float64, λ::Float64, μ::Float64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ`.
"""
function sim_cbd(t::Float64,
                 λ::Float64,
                 μ::Float64)

  tw = cbd_wait(λ, μ)

  if tw > t
    return sTbd(t, false, false)
  end

  if λorμ(λ, μ)
    return sTbd(sim_cbd(t - tw, λ, μ), sim_cbd(t - tw, λ, μ), tw, false, false)
  else
    return sTbd(tw, true, false)
  end
end



"""
    sim_cbd(t ::Float64,
            λ ::Float64,
            μ ::Float64,
            na::Int64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ`.
"""
function sim_cbd(t ::Float64,
                 λ ::Float64,
                 μ ::Float64,
                 na::Int64)

  tw = cbd_wait(λ, μ)

  if tw > t
    na += 1
    return sTbd(t), na
  end

  if λorμ(λ, μ)
    d1, na = sim_cbd(t - tw, λ, μ, na)
    d2, na = sim_cbd(t - tw, λ, μ, na)

    return sTbd(d1, d2, tw), na
  else
    return sTbd(tw, true), na
  end
end




"""
    _sim_cbd_t(t   ::Float64,
               λ   ::Float64,
               μ   ::Float64,
               lr  ::Float64,
               lU  ::Float64,
               Iρi ::Float64,
               na  ::Int64,
               nn ::Int64,
               nlim::Int64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ` for terminal branches.
"""
function _sim_cbd_t(t   ::Float64,
                    λ   ::Float64,
                    μ   ::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    Iρi ::Float64,
                    na  ::Int64,
                    nn ::Int64,
                    nlim::Int64)

  if isfinite(lr) && nn < nlim

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(Iρi * Float64(na)/Float64(na-1))
      end
      
      if nlr < lr && lU >= nlr
        return sTbd(), na, nn, NaN
      else
        return sTbd(t, false, false), na, nn, nlr
      end
    else
      if λorμ(λ, μ)
        nn += 1
        d1, na, nn, lr = _sim_cbd_t(t - tw, λ, μ, lr, lU, Iρi, na, nn, nlim)
        d2, na, nn, lr = _sim_cbd_t(t - tw, λ, μ, lr, lU, Iρi, na, nn, nlim)

        return sTbd(d1, d2, tw, false, false), na, nn, lr
      else
        return sTbd(tw, true, false), na, nn, lr
      end
    end
  end

  return sTbd(), na, nn, NaN
end




"""
    _sim_cbd_i(t   ::Float64,
               λ   ::Float64,
               μ   ::Float64,
               na  ::Int64,
               nn  ::Int64,
               nlim::Int64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ` for internal branches.
"""
function _sim_cbd_i(t   ::Float64,
                    λ   ::Float64,
                    μ   ::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

  if nn < nlim

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      return sTbd(t, false, false), na, nn
    end

    if λorμ(λ, μ)
      nn += 1
      d1, na, nn = _sim_cbd_i(t - tw, λ, μ, na, nn, nlim)
      d2, na, nn = _sim_cbd_i(t - tw, λ, μ, na, nn, nlim)

      return sTbd(d1, d2, tw, false, false), na, nn
    else
      return sTbd(tw, true, false), na, nn
    end
  end

  return sTbd(), na, nn
end




"""
    _sim_cbd_it(t   ::Float64,
                λ   ::Float64,
                μ   ::Float64,
                lr  ::Float64,
                lU  ::Float64,
                Iρi ::Float64,
                nn ::Int64,
                nlim::Int64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ` for continuing internal branches.
"""
function _sim_cbd_it(t   ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     lr  ::Float64,
                     lU  ::Float64,
                     Iρi ::Float64,
                     na  ::Int64,
                     nn ::Int64,
                     nlim::Int64)

  if lU < lr && nn < nlim

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      lr += log(Iρi)
      return sTbd(t, false, false), na, nn, lr
    end

    if λorμ(λ, μ)
      nn += 1
      d1, na, nn, lr = _sim_cbd_it(t - tw, λ, μ, lr, lU, Iρi, na, nn, nlim)
      d2, na, nn, lr = _sim_cbd_it(t - tw, λ, μ, lr, lU, Iρi, na, nn, nlim)

      return sTbd(d1, d2, tw, false, false), na, nn, lr
    else
      return sTbd(tw, true, false), na, nn, lr
    end

  end

  return sTbd(), na, nn, NaN
end




"""
    sim_cbd_surv(t::Float64, λ::Float64, μ::Float64, surv::Bool, nn::Int64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ` until it goes extinct or survives.
"""
function sim_cbd_surv(t   ::Float64,
                      λ   ::Float64,
                      μ   ::Float64,
                      surv::Bool,
                      nn ::Int64)

  if !surv && nn < 500

    tw = cbd_wait(λ, μ)

    if tw > t
      return true, nn
    end

    if λorμ(λ, μ)
      nn += 1
      surv, nn = sim_cbd_surv(t - tw, λ, μ, surv, nn)
      surv, nn = sim_cbd_surv(t - tw, λ, μ, surv, nn)

      return surv, nn
    else
      return surv, nn
    end
  end

  return true, nn
end




"""
   sim_cbd_b(n::Int64, λ::Float64, μ::Float64)

Simulate constant birth-death in backward time.
"""
function sim_cbd_b(n::Int64,
                   λ::Float64,
                   μ::Float64)

  nF = Float64(n)
  nI = n

  # disjoint trees vector
  tv = sTbd[]
  for i in Base.OneTo(nI)
    push!(tv, sTbd(0.0))
  end

  # start simulation
  while true
    w = cbd_wait(nF, λ, μ)

    for t in tv
      adde!(t, w)
    end

    # if speciation
    if λorμ(λ, μ)
      if isone(nI)
        return tv[nI]
      else
        j, k = samp2(Base.OneTo(nI))
        tv[j] = sTbd(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    # if extinction
    else
      nI += 1
      nF += 1.0
      push!(tv, sTbd(0.0, true))
    end
  end
end




"""
    sim_cbd_b(λ::Float64,
              μ::Float64,
              mxth::Float64,
              maxn::Int64)

Simulate constant birth-death in backward time conditioned on 1 survival
and not having a greater tree height than `mxth`.
"""
function sim_cbd_b(λ::Float64,
                   μ::Float64,
                   mxth::Float64,
                   maxn::Int64)

  nF = 1.0
  nI = 1

  # disjoint trees vector
  tv = [sTbd(0.0, false)]

  th = 0.0

  # start simulation
  while true
    w   = cbd_wait(nF, λ, μ)

    # track backward time
    th += w

    if nI > maxn
      return tv[nI], (mxth + 0.1)
    end

    if th > mxth
     return tv[nI], th
    end

    for t in tv
      adde!(t, w)
    end

    # if speciation
    if λorμ(λ, μ)
      if isone(nI)
        return tv[nI], th
      else
        j, k = samp2(Base.OneTo(nI))
        tv[j] = sTbd(tv[j], tv[k], 0.0)
        deleteat!(tv,k)
        nI -= 1
        nF -= 1.0
      end
    # if extinction
    else
      nI += 1
      nF += 1.0
      push!(tv, sTbd(0.0, true))
    end
  end
end




"""
    samp2(o::Base.OneTo{Int64})

Sample `2` without replacement from `o`.
"""
function samp2(o::Base.OneTo{Int64})
  j = rand(o)
  k = rand(o)
  while k == j
    k = rand(o)
  end
  return j, k
end




"""
    cbd_wait(n::Float64, λ::Float64, μ::Float64)

Sample a waiting time for constant birth-death when `n` species
are alive with speciation rate `λ` and extinction rate `μ`.
"""
cbd_wait(n::Float64, λ::Float64, μ::Float64) = rexp(n*(λ + μ))




"""
    cbd_wait(λ::Float64, μ::Float64)

Sample a per-lineage waiting time for constant birth-death species
with speciation rate `λ` and extinction rate `μ`.
"""
cbd_wait(λ::Float64, μ::Float64) = rexp(λ + μ)





"""
    rexp(r::Float64)

Generate an exponential sample with rate `r`.
"""
rexp(r::Float64) = @fastmath randexp()/r



"""
    rexpt(r::Float64)

Generate an exponential sample with rate `r` truncated to be smaller than `t`.
"""
rexpt(r::Float64, t::Float64) = 
  @fastmath μt = -log(1.0 - rand()*(1.0 - exp(- r*t)))/r



"""
    λorμ(λ::Float64, μ::Float64)

Return `true` if speciation event
"""
λorμ(λ::Float64, μ::Float64) = rand() < (λ/(λ + μ))



