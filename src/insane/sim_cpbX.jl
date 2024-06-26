#=

constant pure-birth simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#



"""
  sim_cpb(t::Float64, λ::Float64, x0::Float64, σx::Float64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cpb(t::Float64, λ::Float64, x0::Float64, σx::Float64)

  tw = cpb_wait(λ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    return sTpbX(t, false, x0, x1)
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  sTpbX(sim_cpb(t - tw, λ, x1, σx),
        sim_cpb(t - tw, λ, x1, σx),
        tw, false, x0, x1)
end




"""
    sim_cpb(t::Float64, λ::Float64, x0::Float64, σx::Float64, na::Int64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cpb(t::Float64, λ::Float64, x0::Float64, σx::Float64, na::Int64)

  tw = cpb_wait(λ)

  if tw > t
    na += 1
    return sTpbX(t, false, x0, rnorm(x0, sqrt(t) * σx)), na
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  d1, na = sim_cpb(t - tw, λ, x1, σx, na)
  d2, na = sim_cpb(t - tw, λ, x1, σx, na)

  sTpbX(d1, d2, tw, false, x0, x1), na
end




"""
    _sim_cpb_t(t    ::Float64,
               λ    ::Float64,
               x0   ::Float64,
               σx   ::Float64,
               lr   ::Float64,
               lU   ::Float64,
               Iρi  ::Float64,
               na   ::Int64,
               nn   ::Int64,
               nlim ::Int64,
               xfist::Vector{Float64},
               xfst ::Vector{Float64},
               est  ::Vector{Float64})

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function _sim_cpb_t(t   ::Float64,
                    λ   ::Float64,
                    x0  ::Float64,
                    σx  ::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    Iρi ::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64,
                    xist::Vector{Float64},
                    xfst::Vector{Float64},
                    est ::Vector{Float64})

  if isfinite(lr) && nn < nlim
    tw = cpb_wait(λ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(Iρi * Float64(na)/Float64(na-1))
      end
      if nlr < lr && lU >= nlr
        return sTpbX(), na, nn, NaN
      else
        x1 = rnorm(x0, sqrt(t) * σx)
        push!(xist, x0)
        push!(xfst, x1)
        push!(est, t)
        return sTpbX(t, false, x0, x1), na, nn, nlr
      end
    end

    nn += 1
    x1 = rnorm(x0, sqrt(tw) * σx)
    d1, na, nn, lr = 
      _sim_cpb_t(t - tw, λ, x1, σx, lr, lU, Iρi, na, nn, nlim, xist, xfst, est)
    d2, na, nn, lr = 
      _sim_cpb_t(t - tw, λ, x1, σx, lr, lU, Iρi, na, nn, nlim, xist, xfst, est)

    return sTpbX(d1, d2, tw, false, x0, x1), na, nn, lr
  end

  return sTpbX(), na, nn, NaN
end




"""
    _sim_cpb_i(t   ::Float64,
               λ   ::Float64,
               x0  ::Float64,
               σx  ::Float64,
               nn  ::Int64,
               nlim::Int64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function _sim_cpb_i(t   ::Float64,
                    λ   ::Float64,
                    x0  ::Float64,
                    σx  ::Float64,
                    nn  ::Int64,
                    nlim::Int64)

  if nn < nlim

    tw = cpb_wait(λ)

    if tw > t
        x1 = rnorm(x0, sqrt(t) * σx)

      return sTpbX(t, false, x0, x1), nn
    end

    nn += 1
    x1 = rnorm(x0, sqrt(tw) * σx)

    d1, nn = _sim_cpb_i(t - tw, λ, x1, σx, nn, nlim)
    d2, nn = _sim_cpb_i(t - tw, λ, x1, σx, nn, nlim)

    return sTpbX(d1, d2, tw, false, x0, x1), nn
  end

  sTpbX(), nn
end




"""
    _sim_cpb_it(t   ::Float64,
                λ   ::Float64,
                x0  ::Float64,
                σx  ::Float64,
                lr  ::Float64,
                lU  ::Float64,
                Iρi ::Float64,
                nn  ::Int64,
                nlim::Int64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function _sim_cpb_it(t   ::Float64,
                     λ   ::Float64,
                     x0  ::Float64,
                     σx  ::Float64,
                     lr  ::Float64,
                     lU  ::Float64,
                     Iρi ::Float64,
                     nn  ::Int64,
                     nlim::Int64)

  if lU < lr && nn < nlim

    tw = cpb_wait(λ)

    if tw > t
      lr += log(Iρi)
      return sTpbX(t, false, x0, rnorm(x0, sqrt(t) * σx)), nn, lr
    end

    nn += 1
    x1  = rnorm(x0, sqrt(tw) * σx)
    d1, nn, lr = _sim_cpb_it(t - tw, λ, x1, σx, lr, lU, Iρi, nn, nlim)
    d2, nn, lr = _sim_cpb_it(t - tw, λ, x1, σx, lr, lU, Iρi, nn, nlim)

    return sTpbX(d1, d2, tw, false, x0, x1), nn, lr
  end

  return sTpbX(), nn, NaN
end



