#=

constant fossilized birth-death simulation

Jérémy Andréoletti
Ignacio Quintero Mächler

v(°-°v)
t(-_-t)

Created 07 10 2021
=#




"""
    sim_cfbd(t::Float64, λ::Float64, μ::Float64, ψ::Float64)

Simulate a constant fossilized birth-death `iTree` of height `t` with speciation
rate `λ`, extinction rate `μ` and fossilization rate `ψ`.
"""
function sim_cfbd(t::Float64, λ::Float64, μ::Float64, ψ::Float64)

  tw = cfbd_wait(λ, μ, ψ)

  # if reached the present
  if tw > t
    return sTfbd(t, false, false, false)
  end

  # speciation
  if λevent(λ, μ, ψ)
    return sTfbd(sim_cfbd(t - tw, λ, μ, ψ),
                 sim_cfbd(t - tw, λ, μ, ψ),
                 tw, false, false, false)
  # extinction
  elseif μevent(μ, ψ)
    return sTfbd(tw, true, false, false)
  # fossil sampling
  else
    return sTfbd(sim_cfbd(t - tw, λ, μ, ψ), 
             tw, false, true, false)
  end
end



"""
    sim_cfbd(t  ::Float64,
             λ  ::Float64,
             μ  ::Float64,
             ψ  ::Vector{Float64},
             ψts::Vector{Float64})

Simulate a constant fossilized birth-death `iTree` of height `t` with speciation
rate `λ`, extinction rate `μ` and a vector of fossilization rates `ψ`.
"""
function sim_cfbd(t  ::Float64,
                  λ  ::Float64,
                  μ  ::Float64,
                  ψ  ::Vector{Float64},
                  ψts::Vector{Float64})
  return _sim_cfbd(t, λ, μ, ψ, ψts, 1, lastindex(ψ))
end



"""
    _sim_cfbd(t  ::Float64,
              λ  ::Float64,
              μ  ::Float64,
              ψ  ::Vector{Float64},
              ψts::Vector{Float64},
              ix ::Int64,
              nep::Int64)

Simulate a constant fossilized birth-death `iTree` of height `t` with speciation
rate `λ`, extinction rate `μ` and a vector of fossilization rates `ψ`.
"""
function _sim_cfbd(t  ::Float64,
                   λ  ::Float64,
                   μ  ::Float64,
                   ψ  ::Vector{Float64},
                   ψts::Vector{Float64},
                   ix ::Int64,
                   nep::Int64)

  @inbounds ψi = ψ[ix]

  tw = cfbd_wait(λ, μ, ψi)

  # ψ epoch change
  if ix < nep
    @inbounds ψti = ψts[ix]
    if t - tw < ψti
      e0 = t - ψti
      t0 = _sim_cfbd(ψti, λ, μ, ψ, ψts, ix + 1, nep)
      sete!(t0, e(t0) + e0)
      return t0
    end
  end

  # if reached the present
  if tw > t
    return sTfbd(t, false, false, false)
  end

  # speciation
  if λevent(λ, μ, ψi)
    return sTfbd(_sim_cfbd(t - tw, λ, μ, ψ, ψts, ix, nep),
                 _sim_cfbd(t - tw, λ, μ, ψ, ψts, ix, nep),
                 tw, false, false, false)
  # extinction
  elseif μevent(μ, ψi)
    return sTfbd(tw, true, false, false)
  # fossil sampling
  else
    return sTfbd(_sim_cfbd(t - tw, λ, μ, ψ, ψts, ix, nep), 
             tw, false, true, false)
  end
end




"""
    _sim_cfbd_t(t   ::Float64,
                λ   ::Float64,
                μ   ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ix  ::Int64,
                nep ::Int64,
                lr  ::Float64,
                lU  ::Float64,
                Iρi ::Float64,
                na  ::Int64,
                nn  ::Int64,
                nlim::Int64)

Simulate a constant fossilized birth-death `iTree` of height `t` with speciation
rate `λ`, extinction rate `μ` and fossilization rate `ψ` for terminal branches,
conditioned on no fossilizations.
"""
function _sim_cfbd_t(t   ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     ψ   ::Vector{Float64},
                     ψts ::Vector{Float64},
                     ix  ::Int64,
                     nep ::Int64,
                     lr  ::Float64,
                     lU  ::Float64,
                     Iρi ::Float64,
                     na  ::Int64,
                     nn  ::Int64,
                     nlim::Int64)

  if isfinite(lr) && nn < nlim

    @inbounds ψi  = ψ[ix]

    tw = cfbd_wait(λ, μ, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        t0, na, nn, lr = 
          _sim_cfbd_t(ψti, λ, μ, ψ, ψts, ix + 1, nep, lr, lU, Iρi, na, nn, nlim)
        sete!(t0, e(t0) + e0)
        return t0, na, nn, lr
      end
    end

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(Iρi * Float64(na)/Float64(na-1))
      end
      if nlr < lr && lU >= nlr
        return sTfbd(), na, nn, NaN
      else
        return sTfbd(t, false, false, false), na, nn, nlr
      end
    end

    # speciation
    if λevent(λ, μ, ψi)
      nn += 1
      d1, na, nn, lr =
        _sim_cfbd_t(t - tw, λ, μ, ψ, ψts, ix, nep, lr, lU, Iρi, na, nn, nlim)
      d2, na, nn, lr =
        _sim_cfbd_t(t - tw, λ, μ, ψ, ψts, ix, nep, lr, lU, Iρi, na, nn, nlim)

      return sTfbd(d1, d2, tw, false, false, false), na, nn, lr
    # extinction
    elseif μevent(μ, ψi)

      return sTfbd(tw, true, false, false), na, nn, lr
    # fossil sampling
    else
      return sTfbd(t, false, false, false), na, nn, NaN
    end
  end

  return sTfbd(), na, nn, NaN
end




"""
    _sim_cfbd_i(t   ::Float64,
                te  ::Float64,
                λ   ::Float64,
                μ   ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ix  ::Int64,
                nep ::Int64,
                na  ::Int64,
                nf  ::Int64,
                nn  ::Int64,
                nlim::Int64)

Simulate a constant fossilized birth-death `iTree` of height `t` with
speciation rate `λ`, extinction rate `μ` and fossilization rate `ψ`
for internal branches, conditioned on no fossilizations.
"""
function _sim_cfbd_i(t   ::Float64,
                     te  ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     ψ   ::Vector{Float64},
                     ψts ::Vector{Float64},
                     ix  ::Int64,
                     nep ::Int64,
                     na  ::Int64,
                     nf  ::Int64,
                     nn  ::Int64,
                     nlim::Int64)

  if iszero(nf) && nn < nlim

    @inbounds ψi  = ψ[ix]

    tw = cfbd_wait(λ, μ, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti > te
        e0 = t - ψti
        t0, na, nf, nn  = 
          _sim_cfbd_i(ψti, te, λ, μ, ψ, ψts, ix + 1, nep, na, nf, nn, nlim)
        sete!(t0, e(t0) + e0)
        return t0, na, nf, nn
      end
    end

    if tw > (t - te)
      na += 1
      return sTfbd(t - te, false, false, false), na, nf, nn
    end

    # speciation
    if λevent(λ, μ, ψi)
      nn += 1
      d1, na, nf, nn = 
        _sim_cfbd_i(t - tw, te, λ, μ, ψ, ψts, ix, nep, na, nf, nn, nlim)
      d2, na, nf, nn = 
        _sim_cfbd_i(t - tw, te, λ, μ, ψ, ψts, ix, nep, na, nf, nn, nlim)

      return sTfbd(d1, d2, tw, false, false, false), na, nf, nn
    # extinction
    elseif μevent(μ, ψi)

      return sTfbd(tw, true, false, false), na, nf, nn
    # fossil sampling
    else
      return sTfbd(), na, 1, nn
    end
  end

  return sTfbd(), na, nf, nn
end




"""
    _sim_cfbd_it(t   ::Float64,
                 λ   ::Float64,
                 μ   ::Float64,
                 ψ   ::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 Iρi ::Float64,
                 na  ::Int64,
                 nf  ::Int64,
                 nn  ::Int64,
                 nlim::Int64)

Simulate a constant fossilized birth-death `iTree` of height `t` with
speciation rate `λ`, extinction rate `μ` and fossilization rate `ψ`
for continuing internal branches, conditioned on no fossilizations.
"""
function _sim_cfbd_it(t   ::Float64,
                      λ   ::Float64,
                      μ   ::Float64,
                      ψ   ::Vector{Float64},
                      ψts ::Vector{Float64},
                      ix  ::Int64,
                      nep ::Int64,
                      lr  ::Float64,
                      lU  ::Float64,
                      Iρi ::Float64,
                      na  ::Int64,
                      nn  ::Int64,
                      nlim::Int64)

  if isfinite(lr) && nn < nlim

    @inbounds ψi  = ψ[ix]

    tw = cfbd_wait(λ, μ, ψi)

    # ψ epoch change
    if ix < nep
      ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        t0, na, nn, lr = 
          _sim_cfbd_it(ψti, λ, μ, ψ, ψts, ix + 1, nep, lr, lU, Iρi, na, nn, nlim)
        sete!(t0, e(t0) + e0)
        return t0, na, nn, lr
      end
    end

    if tw > t
      na += 1
      lr += log(Iρi)
      return sTfbd(t, false, false, false), na, nn, lr
    end

    # speciation
    if λevent(λ, μ, ψi)
      nn += 1
      d1, na, nn, lr =
        _sim_cfbd_it(t - tw, λ, μ, ψ, ψts, ix, nep, lr, lU, Iρi, na, nn, nlim)
      d2, na, nn, lr =
        _sim_cfbd_it(t - tw, λ, μ, ψ, ψts, ix, nep, lr, lU, Iρi, na, nn, nlim)

      return sTfbd(d1, d2, tw, false, false, false), na, nn, lr
    # extinction
    elseif μevent(μ, ψi)

      return sTfbd(tw, true, false, false), na, nn, lr
    # fossil sampling
    else
      return sTfbd(), na, nn, NaN
    end
  end

  return sTfbd(), na, nn, NaN
end




"""
    cfbd_wait(n::Float64, λ::Float64, μ::Float64, ψ::Float64)

Sample a waiting time for constant fossilized birth-death when `n` species
are alive with speciation rate `λ` and extinction rate `μ`.
"""
cfbd_wait(n::Float64, λ::Float64, μ::Float64, ψ::Float64) = rexp(n*(λ + μ + ψ))




"""
    cfbd_wait(λ::Float64, μ::Float64, ψ::Float64)

Sample a per-lineage waiting time for constant fossilized birth-death
with speciation rate `λ` and extinction rate `μ`.
"""
cfbd_wait(λ::Float64, μ::Float64, ψ::Float64) = rexp(λ + μ + ψ)




"""
    λevent(λ::Float64, μ::Float64, ψ::Float64)

Return `true` if speciation event
"""
λevent(λ::Float64, μ::Float64, ψ::Float64) = (λ/(λ + μ + ψ)) > rand()




"""
    λevent(μ::Float64, ψ::Float64)

Return `true` if extinction event, conditioned on "not a speciation event"
"""
μevent(μ::Float64, ψ::Float64) = (μ/(μ + ψ)) > rand()



