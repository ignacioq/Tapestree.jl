#=

constant birth-death simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    sim_cbd(t::Float64, λ::Float64, μ::Float64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ`.
"""
function sim_cbd(t ::Float64,
                 λ ::Float64,
                 μ ::Float64,
                 x0::Float64,
                 σx::Float64)

  tw = cbd_wait(λ, μ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    return sTbdX(t, false, false, x0, x1)
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  if λorμ(λ, μ)
    return sTbdX(sim_cbd(t - tw, λ, μ, x1, σx),
                 sim_cbd(t - tw, λ, μ, x1, σx),
                 tw, false, false, x0, x1)
  else
    return sTbdX(tw, true, false, x0, x1)
  end
end





"""
    sim_cbd(t ::Float64,
            λ ::Float64,
            μ ::Float64,
            x0::Float64,
            σx::Float64
            na::Int64)

Simulate a constant birth-death `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cbd(t ::Float64,
                 λ ::Float64,
                 μ ::Float64,
                 x0::Float64,
                 σx::Float64,
                 na::Int64)

  tw = cbd_wait(λ, μ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    na += 1
    return sTbdX(t, false, false, x0, x1), na
  end

  x1 = rnorm(x0, sqrt(tw) * σx)


  if λorμ(λ, μ)
    d1, na = sim_cbd(t - tw, λ, μ, x1, σx, na)
    d2, na = sim_cbd(t - tw, λ, μ, x1, σx, na)

    return sTbdX(d1, d2, tw, false, false, x0, x1), na
  else
    return sTbdX(tw, true, false, x0, x1), na
  end
end





"""
    _sim_cbd_t(t   ::Float64,
               λ   ::Float64,
               μ   ::Float64,
               x0  ::Float64,
               σx  ::Float64,
               lr  ::Float64,
               lU  ::Float64,
               Iρi ::Float64,
               na  ::Int64,
               nn  ::Int64,
               nlim::Int64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ` for terminal branches.
"""
function _sim_cbd_t(t   ::Float64,
                    λ   ::Float64,
                    μ   ::Float64,
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

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(Iρi * Float64(na)/Float64(na-1))
      end
      if nlr < lr && lU >= nlr
        return sTbdX(), na, nn, NaN
      else
        x1 = rnorm(x0, sqrt(t) * σx)
        push!(xist, x0)
        push!(xfst, x1)
        push!(est, t)
        return sTbdX(t, false, false, x0, x1), na, nn, nlr
      end
    else
      if λorμ(λ, μ)
        nn += 1
        x1 = rnorm(x0, sqrt(tw) * σx)

        d1, na, nn, lr = 
          _sim_cbd_t(t - tw, λ, μ, x1, σx, lr, lU, Iρi, na, nn, nlim, 
            xist, xfst, est)
        d2, na, nn, lr = 
          _sim_cbd_t(t - tw, λ, μ, x1, σx, lr, lU, Iρi, na, nn, nlim, 
            xist, xfst, est)

        return sTbdX(d1, d2, tw, false, false, x0, x1), na, nn, lr
      else
        return sTbdX(tw, true, false, x0, rnorm(x0, sqrt(tw) * σx)), na, nn, lr
      end
    end
  end

  return sTbdX(), na, nn, NaN
end




"""
    _sim_cbd_i(t   ::Float64,
               λ   ::Float64,
               μ   ::Float64,
               x0  ::Float64,
               σx  ::Float64,
               na  ::Int64,
               nn  ::Int64,
               nlim::Int64)

Simulate a constant birth-death `iTree` of height `t` with speciation rate `λ`
and extinction rate `μ` for internal branches.
"""
function _sim_cbd_i(t   ::Float64,
                    λ   ::Float64,
                    μ   ::Float64,
                    x0  ::Float64,
                    σx  ::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

  if nn < nlim

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      return sTbdX(t, false, false, x0, rnorm(x0, sqrt(t) * σx)), na, nn
    end

    if λorμ(λ, μ)
      nn += 1
      x1 = rnorm(x0, sqrt(tw) * σx)
      d1, na, nn = _sim_cbd_i(t - tw, λ, μ, x1, σx, na, nn, nlim)
      d2, na, nn = _sim_cbd_i(t - tw, λ, μ, x1, σx, na, nn, nlim)

      return sTbdX(d1, d2, tw, false, false, x0, x1), na, nn
    else
      return sTbdX(tw, true, false, x0, rnorm(x0, sqrt(tw) * σx)), na, nn
    end
  end

  return sTbdX(), na, nn
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
                     x0  ::Float64,
                     σx  ::Float64,
                     lr  ::Float64,
                     lU  ::Float64,
                     Iρi ::Float64,
                     na  ::Int64,
                     nn  ::Int64,
                     nlim::Int64)

  if lU < lr && nn < nlim

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      lr += log(Iρi)
      return sTbdX(t, false, false), na, nn, lr
    end

    if λorμ(λ, μ)
      nn += 1
      d1, na, nn, lr = _sim_cbd_it(t - tw, λ, μ, lr, lU, Iρi, na, nn, nlim)
      d2, na, nn, lr = _sim_cbd_it(t - tw, λ, μ, lr, lU, Iρi, na, nn, nlim)

      return sTbdX(d1, d2, tw, false, false), na, nn, lr
    else
      return sTbdX(tw, true, false), na, nn, lr
    end
  end

  return sTbdX(), na, nn, NaN
end


