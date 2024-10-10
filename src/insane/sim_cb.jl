#=

constant pure-birth simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    cb_wait(λ::Float64)

Sample a per-lineage waiting time for pure-birth species
with speciation rate `λ`.
"""
cb_wait(λ::Float64) = rexp(λ)




"""
    sim_cb(t::Float64, λ::Float64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`.
"""
function sim_cb(t::Float64, λ::Float64)

  tw = cb_wait(λ)

  if tw > t
    return sTb(t)
  end

  sTb(sim_cb(t - tw, λ), sim_cb(t - tw, λ), tw, false)
end




"""
    sim_cb(t::Float64, λ::Float64, na::Int64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`.
"""
function sim_cb(t::Float64, λ::Float64, na::Int64)

  tw = cb_wait(λ)

  if tw > t
    na += 1
    return sTb(t), na
  end
  d1, na = sim_cb(t - tw, λ, na)
  d2, na = sim_cb(t - tw, λ, na)

  sTb(d1, d2, tw, false), na
end




"""
    _sim_cb_t(t   ::Float64,
               λ   ::Float64,
               lr  ::Float64,
               lU  ::Float64,
               Iρi ::Float64,
               na  ::Int64,
               nn  ::Int64,
               nlim::Int64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`
for terminal branches.
"""
function _sim_cb_t(t   ::Float64,
                    λ   ::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    Iρi ::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

  if isfinite(lr) && nn < nlim

    tw = cb_wait(λ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(Iρi * Float64(na)/Float64(na-1))
      end
      if nlr < lr && lU >= nlr
        return sTb(), na, nn, NaN
      else
        return sTb(t, false), na, nn, nlr
      end
    else
      nn += 1
      d1, na, nn, lr = _sim_cb_t(t - tw, λ, lr, lU, Iρi, na, nn, nlim)
      d2, na, nn, lr = _sim_cb_t(t - tw, λ, lr, lU, Iρi, na, nn, nlim)

      return sTb(d1, d2, tw, false), na, nn, lr
    end
  end

  return sTb(), na, nn, NaN
end




"""
    _sim_cb_i(t   ::Float64,
               λ   ::Float64,
               nn ::Int64,
               nlim::Int64,
               rej ::Bool)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`
for interal branches.
"""
function _sim_cb_i(t   ::Float64,
                    λ   ::Float64,
                    nn ::Int64,
                    nlim::Int64)

  if nn < nlim

    tw = cb_wait(λ)

    if tw > t
      return sTb(t, false), nn
    else
      nn += 1
      d1, nn = _sim_cb_i(t - tw, λ, nn, nlim)
      d2, nn = _sim_cb_i(t - tw, λ, nn, nlim)

      return sTb(d1, d2, tw, false), nn
    end
  end

  return sTb(), nn
end




"""
    _sim_cb_it(t   ::Float64,
                λ   ::Float64,
                lr  ::Float64,
                lU  ::Float64,
                Iρi ::Float64,
                nn ::Int64,
                nlim::Int64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`
for continuing internal branches.
"""
function _sim_cb_it(t   ::Float64,
                     λ   ::Float64,
                     lr  ::Float64,
                     lU  ::Float64,
                     Iρi ::Float64,
                     nn ::Int64,
                     nlim::Int64)

  if lU < lr && nn < nlim

    tw = cb_wait(λ)

    if tw > t
      lr += log(Iρi)
      return sTb(t, false), nn, lr
    else
      nn += 1
      d1, nn, lr = _sim_cb_it(t - tw, λ, lr, lU, Iρi, nn, nlim)
      d2, nn, lr = _sim_cb_it(t - tw, λ, lr, lU, Iρi, nn, nlim)

      return sTb(d1, d2, tw, false), nn, lr
    end
  end

  return sTb(), nn, NaN
end




