#=

constant pure-birth simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#



"""
  sim_cb(t::Float64, λ::Float64, x0::Float64, σx::Float64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cb(t::Float64, λ::Float64, x0::Float64, σx::Float64)

  tw = cb_wait(λ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    return sTbX(t, false, x0, x1)
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  sTbX(sim_cb(t - tw, λ, x1, σx),
        sim_cb(t - tw, λ, x1, σx),
        tw, false, x0, x1)
end




"""
    sim_cb(t::Float64, λ::Float64, x0::Float64, σx::Float64, na::Int64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cb(t::Float64, λ::Float64, x0::Float64, σx::Float64, na::Int64)

  tw = cb_wait(λ)

  if tw > t
    na += 1
    return sTbX(t, false, x0, rnorm(x0, sqrt(t) * σx)), na
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  d1, na = sim_cb(t - tw, λ, x1, σx, na)
  d2, na = sim_cb(t - tw, λ, x1, σx, na)

  sTbX(d1, d2, tw, false, x0, x1), na
end




"""
    _sim_cb_t(t    ::Float64,
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
function _sim_cb_t(t   ::Float64,
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
    tw = cb_wait(λ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(Iρi * Float64(na)/Float64(na-1))
      end
      if nlr < lr && lU >= nlr
        return sTbX(), na, nn, NaN
      else
        x1 = rnorm(x0, sqrt(t) * σx)
        push!(xist, x0)
        push!(xfst, x1)
        push!(est, t)
        return sTbX(t, false, x0, x1), na, nn, nlr
      end
    end

    nn += 1
    x1 = rnorm(x0, sqrt(tw) * σx)
    d1, na, nn, lr = 
      _sim_cb_t(t - tw, λ, x1, σx, lr, lU, Iρi, na, nn, nlim, xist, xfst, est)
    d2, na, nn, lr = 
      _sim_cb_t(t - tw, λ, x1, σx, lr, lU, Iρi, na, nn, nlim, xist, xfst, est)

    return sTbX(d1, d2, tw, false, x0, x1), na, nn, lr
  end

  return sTbX(), na, nn, NaN
end




"""
    _sim_cb_i(t   ::Float64,
               λ   ::Float64,
               x0  ::Float64,
               σx  ::Float64,
               nn  ::Int64,
               nlim::Int64)

Simulate a constant pure-birth `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function _sim_cb_i(t   ::Float64,
                    λ   ::Float64,
                    x0  ::Float64,
                    σx  ::Float64,
                    nn  ::Int64,
                    nlim::Int64)

  if nn < nlim

    tw = cb_wait(λ)

    if tw > t
        x1 = rnorm(x0, sqrt(t) * σx)

      return sTbX(t, false, x0, x1), nn
    end

    nn += 1
    x1 = rnorm(x0, sqrt(tw) * σx)

    d1, nn = _sim_cb_i(t - tw, λ, x1, σx, nn, nlim)
    d2, nn = _sim_cb_i(t - tw, λ, x1, σx, nn, nlim)

    return sTbX(d1, d2, tw, false, x0, x1), nn
  end

  sTbX(), nn
end




"""
    _sim_cb_it(t   ::Float64,
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
function _sim_cb_it(t   ::Float64,
                     λ   ::Float64,
                     x0  ::Float64,
                     σx  ::Float64,
                     lr  ::Float64,
                     lU  ::Float64,
                     Iρi ::Float64,
                     nn  ::Int64,
                     nlim::Int64)

  if lU < lr && nn < nlim

    tw = cb_wait(λ)

    if tw > t
      lr += log(Iρi)
      return sTbX(t, false, x0, rnorm(x0, sqrt(t) * σx)), nn, lr
    end

    nn += 1
    x1  = rnorm(x0, sqrt(tw) * σx)
    d1, nn, lr = _sim_cb_it(t - tw, λ, x1, σx, lr, lU, Iρi, nn, nlim)
    d2, nn, lr = _sim_cb_it(t - tw, λ, x1, σx, lr, lU, Iρi, nn, nlim)

    return sTbX(d1, d2, tw, false, x0, x1), nn, lr
  end

  return sTbX(), nn, NaN
end



