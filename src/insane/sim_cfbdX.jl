#=

constant birth-death simulation

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    sim_cfbd(t ::Float64,
             λ ::Float64,
             μ ::Float64,
             ψ ::Float64,
             x0::Float64,
             σx::Float64)

Simulate a constant fossilized birth-death `iTree` of height `t` 
with speciation rate `λ`, extinction rate `μ` and fossil sampling rate `ψ`.
"""
function sim_cfbd(t ::Float64,
                  λ ::Float64,
                  μ ::Float64,
                  ψ ::Float64,
                  x0::Float64,
                  σx::Float64)

  tw = cfbd_wait(λ, μ, ψ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    return sTfbdx(t, false, false, false, x0, x1)
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  if λevent(λ, μ, ψ)
    # speciation
    return sTfbdx(sim_cfbd(t - tw, λ, μ, ψ, x1, σx), 
                  sim_cfbd(t - tw, λ, μ, ψ, x1, σx), 
                  tw, false, false, false, x0, x1)
  elseif μevent(μ, ψ)
    # extinction
    return sTfbdx(tw, true, false, false, x0, x1)
  else
    # fossil sampling
    return sTfbdx(sim_cfbd(t - tw, λ, μ, ψ, x1, σx), 
                  tw, false, true, false, x0, x1)
  end
end





"""
    sim_cfbd(t ::Float64,
             λ ::Float64,
             μ ::Float64,
             ψ ::Float64,
             x0::Float64,
             σx::Float64,
             na::Int64)

Simulate a constant  fossilized birth-death `iTree` with traits of height `t`
with speciation rate `λ` starting at trait `x0` with rate `σx`.
"""
function sim_cfbd(t ::Float64,
                  λ ::Float64,
                  μ ::Float64,
                  ψ ::Float64,
                  x0::Float64,
                  σx::Float64,
                  na::Int64,
                  nf::Int64)

  tw = cfbd_wait(λ, μ, ψ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σx)
    na += 1
    return sTfbdx(t, false, false, false, x0, x1), na, nf
  end

  x1 = rnorm(x0, sqrt(tw) * σx)

  # speciation
  if λevent(λ, μ, ψ)
    d1, na, nf = sim_cfbd(t - tw, λ, μ, ψ, x1, σx, na, nf)
    d2, na, nf = sim_cfbd(t - tw, λ, μ, ψ, x1, σx, na, nf)

    return sTfbdx(d1, d2, tw, false, false, false, x0, x1), na, nf
  # extinction
  elseif μevent(μ, ψ)
    return sTfbdx(tw, true, false, false, x0, x1), na, nf
  # fossil sampling
  else
    nf += 1
    d1, na, nf = sim_cfbd(t - tw, λ, μ, ψ, x1, σx, na, nf)
    return sTfbdx(d1, tw, false, true, false, x0, x1), na, nf
  end
end




"""
     _sim_cfbd_t(t   ::Float64,
                 λ   ::Float64,
                 μ   ::Float64,
                 ψ   ::Float64,
                 x0  ::Float64,
                 σx  ::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 iρi ::Float64,
                 na  ::Int64,
                 nn  ::Int64,
                 nlim::Int64,
                 xist::Vector{Float64},
                 xfst::Vector{Float64},
                 est ::Vector{Float64})

Simulate a constant fossilized birth-death `iTree` of height `t` with speciation
rate `λ`, extinction rate `μ` and fossilization rate `ψ` for terminal branches, 
conditioned on no fossilizations.
"""
function _sim_cfbd_t(t   ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     ψ   ::Float64,
                     x0  ::Float64,
                     σx  ::Float64,
                     lr  ::Float64,
                     lU  ::Float64,
                     iρi ::Float64,
                     na  ::Int64,
                     nn  ::Int64,
                     nlim::Int64,
                     xist::Vector{Float64},
                     xfst::Vector{Float64},
                     est ::Vector{Float64})

  if isfinite(lr) && nn < nlim

    tw = cfbd_wait(λ, μ, ψ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(iρi * Float64(na)/Float64(na-1))
      end
      if nlr < lr && lU >= nlr
        return sTfbdx(), na, nn, NaN
      else
        x1 = rnorm(x0, sqrt(t) * σx)
        push!(xist, x0)
        push!(xfst, x1)
        push!(est, t)
        return sTfbdx(t, false, false, false, x0, x1), na, nn, nlr
      end
    end

    # speciation
    if λevent(λ, μ, ψ)
      nn += 1
      x1 = rnorm(x0, sqrt(tw) * σx)

      d1, na, nn, lr = 
        _sim_cfbd_t(t - tw, λ, μ, ψ, x1, σx, lr, lU, iρi, na, nn, nlim, 
          xist, xfst, est)
      d2, na, nn, lr = 
        _sim_cfbd_t(t - tw, λ, μ, ψ, x1, σx, lr, lU, iρi, na, nn, nlim, 
          xist, xfst, est)

      return sTfbdx(d1, d2, tw, false, false, false, x0, x1), na, nn, lr
    # extinction
    elseif μevent(μ, ψ)

      return sTfbdx(tw, true, false, false, x0, rnorm(x0, sqrt(tw) * σx)), 
             na, nn, lr
    # fossil sampling
    else
      return sTfbdx(), na, nn, NaN
    end
  end

  return sTfbdx(), na, nn, NaN
end




"""
    _sim_cfbd_i(t   ::Float64,
                λ   ::Float64,
                μ   ::Float64,
                ψ   ::Float64,
                x0  ::Float64,
                σx  ::Float64,
                na  ::Int64,
                nf  ::Int64,
                nn  ::Int64,
                nlim::Int64,
                xist::Vector{Float64},
                xfst::Vector{Float64},
                est ::Vector{Float64})

Simulate a constant fossilized birth-death `iTree` of height `t` with 
speciation rate `λ`, extinction rate `μ` and fossilization rate `ψ` 
for internal branches, conditioned on no fossilizations.
"""
function _sim_cfbd_i(t   ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     ψ   ::Float64,
                     x0  ::Float64,
                     σx  ::Float64,
                     na  ::Int64,
                     nf  ::Int64,
                     nn  ::Int64,
                     nlim::Int64,
                     xist::Vector{Float64},
                     xfst::Vector{Float64},
                     est ::Vector{Float64})


  if iszero(nf) && nn < nlim

    tw = cfbd_wait(λ, μ, ψ)

    if tw > t
      na += 1
      x1 = rnorm(x0, sqrt(t) * σx)
      push!(xist, x0)
      push!(xfst, x1)
      push!(est, t)
      return sTfbdx(t, false, false, false, x0, x1), na, nf, nn
    end

    # speciation
    if λevent(λ, μ, ψ)
      nn += 1
      x1 = rnorm(x0, sqrt(tw) * σx)
      d1, na, nf, nn = 
        _sim_cfbd_i(t - tw, λ, μ, ψ, x1, σx, na, nf, nn, nlim, xist, xfst, est)
      d2, na, nf, nn = 
        _sim_cfbd_i(t - tw, λ, μ, ψ, x1, σx, na, nf, nn, nlim, xist, xfst, est)

      return sTfbdx(d1, d2, tw, false, false, false, x0, x1), na, nf, nn
    # extinction
    elseif μevent(μ, ψ)

      return sTfbdx(tw, true, false, false, x0, rnorm(x0, sqrt(tw) * σx)), 
             na, nf, nn
    # fossil sampling
    else
      nf += 1
      return sTfbdx(), na, nf, nn
    end
  end

  return sTfbdx(), na, nf, nn
end




"""
    _sim_cfbd_i(t   ::Float64,
                λ   ::Float64,
                μ   ::Float64,
                ψ   ::Float64,
                x0  ::Float64,
                σx  ::Float64,
                na  ::Int64,
                nf  ::Int64,
                nn  ::Int64,
                nlim::Int64,
                xist::Vector{Float64},
                xfst::Vector{Float64},
                est ::Vector{Float64})

Simulate a constant fossilized birth-death `iTree` of height `t` with 
speciation rate `λ`, extinction rate `μ` and fossilization rate `ψ` 
for internal branches, conditioned on no fossilizations.
"""
function _sim_cfbd_i(t   ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     ψ   ::Float64,
                     x0  ::Float64,
                     σx  ::Float64,
                     na  ::Int64,
                     nf  ::Int64,
                     nn  ::Int64,
                     nlim::Int64)

  if iszero(nf) && nn < nlim

    tw = cfbd_wait(λ, μ, ψ)

    if tw > t
      na += 1
      return sTfbdx(t, false, false, false, x0, rnorm(x0, sqrt(t) * σx)), 
             na, nf, nn
    end

    # speciation
    if λevent(λ, μ, ψ)
      nn += 1
      x1 = rnorm(x0, sqrt(tw) * σx)
      d1, na, nf, nn = _sim_cfbd_i(t - tw, λ, μ, ψ, x1, σx, na, nf, nn, nlim)
      d2, na, nf, nn = _sim_cfbd_i(t - tw, λ, μ, ψ, x1, σx, na, nf, nn, nlim)

      return sTfbdx(d1, d2, tw, false, false, false, x0, x1), na, nf, nn
    # extinction
    elseif μevent(μ, ψ)

      return sTfbdx(tw, true, false, false, x0, rnorm(x0, sqrt(tw) * σx)), 
             na, nf, nn
    # fossil sampling
    else
      nf += 1
      return sTfbdx(), na, nf, nn
    end
  end

  return sTfbdx(), na, nf, nn
end





"""
    _sim_cfbd_it(t   ::Float64,
                 λ   ::Float64,
                 μ   ::Float64,
                 ψ   ::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 iρi ::Float64,
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
                      ψ   ::Float64,
                      x0  ::Float64,
                      σx  ::Float64,
                      lr  ::Float64,
                      lU  ::Float64,
                      iρi ::Float64,
                      na  ::Int64,
                      nn  ::Int64,
                      nlim::Int64)

  if lU < lr && nn < nlim

    tw = cfbd_wait(λ, μ, ψ)

    if tw > t
      na += 1
      lr += log(iρi)
      return sTfbdx(t, false, false, false, x0, rnorm(x0, sqrt(t) * σx)), 
             na, nn, lr
    end

    # speciation
    if λevent(λ, μ, ψ)
      nn += 1
      x1 = rnorm(x0, sqrt(tw) * σx)
      d1, na, nn, lr = 
        _sim_cfbd_it(t - tw, λ, μ, ψ, x1, σx, lr, lU, iρi, na, nn, nlim)
      d2, na, nn, lr = 
        _sim_cfbd_it(t - tw, λ, μ, ψ, x1, σx, lr, lU, iρi, na, nn, nlim)

      return sTfbdx(d1, d2, tw, false, false, false, x0, x1), na, nn, lr
    # extinction
    elseif μevent(μ, ψ)

      return sTfbdx(tw, true, false, false, x0, rnorm(x0, sqrt(tw) * σx)), 
             na, nn, lr
    # fossil sampling
    else
      return sTfbdx(), na, nn, NaN
    end
  end

  return sTfbdx(), na, nn, NaN
end

