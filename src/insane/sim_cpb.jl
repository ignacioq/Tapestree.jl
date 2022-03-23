#=

constant pure-birth simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    cpb_wait(λ::Float64)

Sample a per-lineage waiting time for pure-birth species
with speciation rate `λ`.
"""
cpb_wait(λ::Float64) = rexp(λ)




"""
    sim_cpb(t::Float64, λ::Float64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`.
"""
function sim_cpb(t::Float64, λ::Float64)

  tw = cpb_wait(λ)

  if tw > t
    return sTpb(t)
  end

  sTpb(sim_cpb(t - tw, λ), sim_cpb(t - tw, λ), tw, false)
end




"""
    sim_cpb(t::Float64, λ::Float64, na::Int64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`.
"""
function sim_cpb(t::Float64, λ::Float64, na::Int64)

  tw = cpb_wait(λ)

  if tw > t
    na += 1
    return sTpb(t), na
  end
  d1, na = sim_cpb(t - tw, λ, na)
  d2, na = sim_cpb(t - tw, λ, na)

  sTpb(d1, d2, tw, false), na
end




"""
    _sim_cpb_t(t   ::Float64,
               λ   ::Float64,
               lr  ::Float64,
               lU  ::Float64,
               Iρi ::Float64,
               na  ::Int64,
               nsp ::Int64,
               nlim::Int64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`
for terminal branches.
"""
function _sim_cpb_t(t   ::Float64,
                    λ   ::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    Iρi ::Float64,
                    na  ::Int64,
                    nsp ::Int64,
                    nlim::Int64)

  if isfinite(lr) && nsp < nlim

    tw = cpb_wait(λ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(Iρi * Float64(na)/Float64(na-1))
      end
      if nlr >= lr
        return sTpb(t, false), na, nsp, nlr
      elseif lU < nlr
        return sTpb(t, false), na, nsp, nlr
      else
        return sTpb(0.0, false), na, nsp, NaN
      end
    else
      nsp += 1
      d1, na, nsp, lr = _sim_cpb_t(t - tw, λ, lr, lU, Iρi, na, nsp, nlim)
      d2, na, nsp, lr = _sim_cpb_t(t - tw, λ, lr, lU, Iρi, na, nsp, nlim)

      return sTpb(d1, d2, tw, false), na, nsp, lr
    end
  end

  return sTpb(0.0, false), na, nsp, NaN
end




"""
    _sim_cpb_i(t   ::Float64,
              λ   ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nsp ::Int64,
              nlim::Int64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`
for interal branches.
"""
function _sim_cpb_i(t   ::Float64,
                    λ   ::Float64,
                    nsp ::Int64,
                    nlim::Int64)

  if nsp < nlim

    tw = cpb_wait(λ)

    if tw > t
      return sTpb(t, false), nsp
    else
      nsp += 1
      d1, nsp = _sim_cpb_i(t - tw, λ, nsp, nlim)
      d2, nsp = _sim_cpb_i(t - tw, λ, nsp, nlim)

      return sTpb(d1, d2, tw, false), nsp
    end
  end

  return sTpb(0.0, false), nsp
end




"""
    _sim_cpb_it(t   ::Float64,
                λ   ::Float64,
                lr  ::Float64,
                lU  ::Float64,
                Iρi ::Float64,
                nsp ::Int64,
                nlim::Int64)

Simulate a constant pure-birth `iTree` of height `t` with speciation rate `λ`
for continuing internal branches.
"""
function _sim_cpb_it(t   ::Float64,
                     λ   ::Float64,
                     lr  ::Float64,
                     lU  ::Float64,
                     Iρi ::Float64,
                     nsp ::Int64,
                     nlim::Int64)

  if lU < lr && nsp < nlim

    tw = cpb_wait(λ)

    if tw > t
      lr += log(Iρi)
      return sTpb(t, false), nsp, lr
    else
      nsp += 1
      d1, nsp, lr = _sim_cpb_it(t - tw, λ, lr, lU, Iρi, nsp, nlim)
      d2, nsp, lr = _sim_cpb_it(t - tw, λ, lr, lU, Iρi, nsp, nlim)

      return sTpb(d1, d2, tw, false), nsp, lr
    end
  end

  return sTpb(0.0, false), nsp, NaN
end




