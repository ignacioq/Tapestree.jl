#=

constant fossil punkeek simulation

Ignacio Quintero Mächler

t(-_-t)

Created 20 11 2024
=#




"""
    sim_cfpe(t ::Float64,
             λ ::Float64,
             μ ::Float64,
             ψ ::Float64,
             x0::Float64,
             σa::Float64,
             σk::Float64)

Simulate a constant fossil punkeek model of height `t` with speciation rate `λ`,
extinction rate `μ`, fossilization rate `ψ`, starting trait value `x0`, 
and anagenetic and cladogenetic variance `σa` and `σk`.
"""
function sim_cfpe(t ::Float64,
                  λ ::Float64,
                  μ ::Float64,
                  ψ ::Float64,
                  x0::Float64,
                  σa::Float64,
                  σk::Float64)

  tw = cfbd_wait(λ, μ, ψ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σa)
    return sTfpe(t, false, false, x0, x1, false, false)
  end

  x1 = rnorm(x0, sqrt(tw) * σa)

  # speciation
  if λevent(λ, μ, ψ)
    xk  = rnorm(x1, σk)
    shi = rand(Bool)
    xl, xr = if shi xk, x1 else x1, xk end

    return sTfpe(sim_cfpe(t - tw, λ, μ, ψ, xl, σa, σk), 
                 sim_cfpe(t - tw, λ, μ, ψ, xr, σa, σk), 
                 tw, false, false, x0, x1, shi, false)

  # extinction
  elseif μevent(μ, ψ)
    return sTfpe(tw, true, false, x0, x1, false, false)

  # fossil
  else
    return sTfpe(sim_cfpe(t - tw, λ, μ, ψ, x1, σa, σk), 
                 tw, false, true, x0, x1, false, false)
  end
end



"""
here
"""


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
    sim_cfpe(t ::Float64,
            λ ::Float64,
            μ ::Float64,
            x0::Float64,
            σa::Float64,
            σk::Float64)

Simulate a constant punkeek model of height `t` with speciation rate `λ`,
extinction rate `μ`, starting trait value `x0`, and anagenetic and 
cladogenetic variance `σa` and `σk`.
"""
function sim_cfpe(t ::Float64,
                  λ ::Float64,
                  μ ::Float64,
                  x0::Float64,
                  σa::Float64,
                  σk::Float64)

  tw = cbd_wait(λ, μ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σa)
    return sTfpe(t, false, x0, x1, false, false)
  end

  x1 = rnorm(x0, sqrt(tw) * σa)

  if λorμ(λ, μ)
    xk = rnorm(x1, σk)
    shi = rand(Bool)
    xl, xr = if shi xk, x1 else x1, xk end

    return sTfpe(sim_cfpe(t - tw, λ, μ, xl, σa, σk), 
                 sim_cfpe(t - tw, λ, μ, xr, σa, σk), 
                 tw, false, x0, x1, shi, false)
  else
    return sTfpe(tw, true, x0, x1, false, false)
  end
end




"""
    _sim_cfpe_t(t   ::Float64,
               λ   ::Float64,
               μ   ::Float64,
               x0  ::Float64,
               σa  ::Float64,
               σk  ::Float64,
               lr  ::Float64,
               lU  ::Float64,
               iρi ::Float64,
               na  ::Int64,
               nn  ::Int64,
               nlim::Int64,
               xist::Vector{Float64},
               xfst::Vector{Float64},
               est ::Vector{Float64})

Simulate a constant punkeek model of height `t` with speciation rate `λ`,
extinction rate `μ`, starting trait value `x0`, and anagenetic and 
cladogenetic variance `σa` and `σk` for terminal branches.
"""
function _sim_cfpe_t(t   ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     ψ   ::Vector{Float64},
                     x0  ::Float64,
                     σa  ::Float64,
                     σk  ::Float64,
                     ψts ::Vector{Float64},
                     ix  ::Int64,
                     nep ::Int64,
                     lr  ::Float64,
                     iρi ::Float64,
                     na  ::Int64,
                     nn  ::Int64,
                     nlim::Int64,
                     xist::Vector{Float64},
                     xfst::Vector{Float64},
                     est ::Vector{Float64})

  if isfinite(lr) && nn < nlim

    @inbounds ψi  = ψ[ix]

    tw = cfbd_wait(λ, μ, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        x1 = rnorm(x0, sqrt(e0) * σa)
        t0, na, nn, lr = 
          _sim_cfpe_t(ψti, λ, μ, ψ, x1, σa, σk, ψts, ix + 1, nep, 
            lr, iρi, na, nn, nlim, xist, xfst, est)
        sete!(t0, e(t0) + e0)
        setxi!(t0, x0)
        return t0, na, nn, lr
      end
    end

    # if tip
    if tw > t
      na += 1
      if na > 1
        lr += log(iρi * Float64(na)/Float64(na-1))
      end
      if isfinite(lr)
        x1 = rnorm(x0, sqrt(t) * σa)
        push!(xist, x0)
        push!(xfst, x1)
        push!(est, t)
        return sTfpe(t, false, false, x0, x1, false, false), na, nn, lr
      else
        return sTfpe(), na, nn, NaN
      end
    end

    x1 = rnorm(x0, sqrt(tw) * σa)

    # if speciation
    if λevent(λ, μ, ψi)
      nn += 1
      xk  = rnorm(x1, σk)
      shi = rand(Bool)
      xl, xr = if shi xk, x1 else x1, xk end

      d1, na, nn, lr = 
        _sim_cfpe_t(t - tw, λ, μ, ψ, xl, σa, σk, ψts, ix, nep, lr, iρi, 
          na, nn, nlim, xist, xfst, est)
      d2, na, nn, lr = 
        _sim_cfpe_t(t - tw, λ, μ, ψ, xr, σa, σk, ψts, ix, nep, lr, iρi, 
          na, nn, nlim, xist, xfst, est)

      return sTfpe(d1, d2, tw, false, false, x0, x1, shi, false), na, nn, lr
    # if extinction
    elseif μevent(μ, ψi)

      return sTfpe(tw, true, false, x0, x1, false, false), na, nn, lr
    # if fossil sample
    else
      return sTfpe(), na, nn, NaN
    end
  end

  return sTfpe(), na, nn, NaN
end




"""
    _sim_cfpe_i(t   ::Float64,
                te  ::Float64,
                λ   ::Float64,
                μ   ::Float64,
                ψ   ::Vector{Float64},
                x0  ::Float64,
                σa  ::Float64,
                σk  ::Float64,
                ψts ::Vector{Float64},
                ix  ::Int64,
                nep ::Int64,
                na  ::Int64,
                nf  ::Int64,
                nn  ::Int64,
                nlim::Int64,
                xfst::Vector{Float64})

Simulate a constant punkeek model of height `t` with speciation rate `λ`,
extinction rate `μ`, starting trait value `x0`, and anagenetic and 
cladogenetic variance `σa` and `σk` for terminal branches.
"""
function _sim_cfpe_i(t   ::Float64,
                     te  ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     ψ   ::Vector{Float64},
                     x0  ::Float64,
                     σa  ::Float64,
                     σk  ::Float64,
                     ψts ::Vector{Float64},
                     ix  ::Int64,
                     nep ::Int64,
                     na  ::Int64,
                     nf  ::Int64,
                     nn  ::Int64,
                     nlim::Int64,
                     xfst::Vector{Float64})

  if iszero(nf) && nn < nlim

    @inbounds ψi  = ψ[ix]

    tw = cfbd_wait(λ, μ, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti > te
        e0 = t - ψti
        x1 = rnorm(x0, sqrt(e0) * σa)
        t0, na, nf, nn = 
          _sim_cfpe_i(ψti, te, λ, μ, ψ, x1, σa, σk, ψts, ix + 1, nep, 
            na, nf, nn, nlim, xfst)
        sete!(t0, e(t0) + e0)
        setxi!(t0, x0)
        return t0, na, nf, nn
      end
    end

    # if tip
    if tw > (t - te)
      na += 1
      x1 = rnorm(x0, sqrt(t) * σa)
      push!(xfst, x1)
      return sTfpe(t - te, false, false, x0, x1, false, false), na, nf, nn
    end

    x1 = rnorm(x0, sqrt(tw) * σa)

    # if speciation
    if λevent(λ, μ, ψi)
      nn += 1
      xk  = rnorm(x1, σk)
      shi = rand(Bool)
      xl, xr = if shi xk, x1 else x1, xk end

      d1, na, nf, nn = 
        _sim_cfpe_i(t - tw, te, λ, μ, ψ, xl, σa, σk, ψts, ix, nep, 
          na, nf, nn, nlim, xfst)
      d2, na, nf, nn = 
        _sim_cfpe_i(t - tw, te, λ, μ, ψ, xr, σa, σk, ψts, ix, nep, 
          na, nf, nn, nlim, xfst)

      return sTfpe(d1, d2, tw, false, false, x0, x1, shi, false), na, nf, nn
    # if extinction
    elseif μevent(μ, ψi)

      return sTfpe(tw, true, false, x0, x1, false, false), na, nf, nn
    # if fossil sample
    else
      return sTfpe(), na, 1, nn
    end
  end

  return sTfpe(), na, nf, nn
end




"""
    _sim_cfpe_i(t   ::Float64,
                te  ::Float64,
                λ   ::Float64,
                μ   ::Float64,
                ψ   ::Vector{Float64},
                x0  ::Float64,
                σa  ::Float64,
                σk  ::Float64,
                ψts ::Vector{Float64},
                ix  ::Int64,
                nep ::Int64,
                na  ::Int64,
                nf  ::Int64,
                nn  ::Int64,
                nlim::Int64,
                xist::Vector{Float64},
                xfst::Vector{Float64},
                est ::Vector{Float64})

Simulate a constant punkeek model of height `t` with speciation rate `λ`,
extinction rate `μ`, starting trait value `x0`, and anagenetic and 
cladogenetic variance `σa` and `σk` for terminal branches.
"""
function _sim_cfpe_i(t   ::Float64,
                     te  ::Float64,
                     λ   ::Float64,
                     μ   ::Float64,
                     ψ   ::Vector{Float64},
                     x0  ::Float64,
                     σa  ::Float64,
                     σk  ::Float64,
                     ψts ::Vector{Float64},
                     ix  ::Int64,
                     nep ::Int64,
                     na  ::Int64,
                     nf  ::Int64,
                     nn  ::Int64,
                     nlim::Int64,
                     xist::Vector{Float64},
                     xfst::Vector{Float64},
                     est ::Vector{Float64})

  if iszero(nf) && nn < nlim

    @inbounds ψi  = ψ[ix]

    tw = cfbd_wait(λ, μ, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti > te
        e0 = t - ψti
        x1 = rnorm(x0, sqrt(e0) * σa)
        t0, na, nf, nn = 
          _sim_cfpe_i(ψti, te, λ, μ, ψ, x1, σa, σk, ψts, ix + 1, nep, 
            na, nf, nn, nlim, xist, xfst, est)
        sete!(t0, e(t0) + e0)
        setxi!(t0, x0)
        return t0, na, nf, nn
      end
    end

    # if tip
    if tw > (t - te)
      na += 1
      x1 = rnorm(x0, sqrt(t) * σa)
      push!(xist, x0)
      push!(xfst, x1)
      push!(est, t)
      return sTfpe(t - te, false, false, x0, x1, false, false), na, nf, nn
    end

    x1 = rnorm(x0, sqrt(tw) * σa)

    # if speciation
    if λevent(λ, μ, ψi)
      nn += 1
      xk  = rnorm(x1, σk)
      shi = rand(Bool)
      xl, xr = if shi xk, x1 else x1, xk end

      d1, na, nf, nn = 
        _sim_cfpe_i(t - tw, te, λ, μ, ψ, xl, σa, σk, ψts, ix, nep, 
          na, nf, nn, nlim, xist, xfst, est)
      d2, na, nf, nn = 
        _sim_cfpe_i(t - tw, te, λ, μ, ψ, xr, σa, σk, ψts, ix, nep, 
          na, nf, nn, nlim, xist, xfst, est)

      return sTfpe(d1, d2, tw, false, false, x0, x1, shi, false), na, nf, nn
    # if extinction
    elseif μevent(μ, ψi)

      return sTfpe(tw, true, false, x0, x1, false, false), na, nf, nn
    # if fossil sample
    else
      return sTfpe(), na, 1, nn
    end
  end

  return sTfpe(), na, nf, nn
end




"""
    _sim_cfpe_it(t   ::Float64,
                λ   ::Float64,
                μ   ::Float64,
                x0  ::Float64,
                σa  ::Float64,
                σk  ::Float64,
                lr  ::Float64,
                lU  ::Float64,
                iρi ::Float64,
                na  ::Int64,
                nn  ::Int64,
                nlim::Int64)

Simulate a constant punkeek model of height `t` with speciation rate `λ`,
extinction rate `μ`, starting trait value `x0`, and anagenetic and 
cladogenetic variance `σa` and `σk` for internal-terminal branches.
"""
function _sim_cfpe_it(t   ::Float64,
                      λ   ::Float64,
                      μ   ::Float64,
                      ψ   ::Vector{Float64},
                      x0  ::Float64,
                      σa  ::Float64,
                      σk  ::Float64,
                      ψts ::Vector{Float64},
                      ix  ::Int64,
                      nep ::Int64,
                      lr  ::Float64,
                      lU  ::Float64,
                      iρi ::Float64,
                      na  ::Int64,
                      nn  ::Int64,
                      nlim::Int64)

  if lU < lr && nn < nlim

    @inbounds ψi  = ψ[ix]

    tw = cfbd_wait(λ, μ, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        x1 = rnorm(x0, sqrt(e0) * σa)
        t0, na, nn, lr = 
          _sim_cfpe_it(ψti, λ, μ, ψ, x1, σa, σk, ψts, ix + 1, nep, 
            lr, lU, iρi, na, nn, nlim)
        sete!(t0, e(t0) + e0)
        setxi!(t0, x0)
        return t0, na, nn, lr
      end
    end

    # if tip
    if tw > t
      na += 1
      lr += log(iρi)
      if isfinite(lr)
        x1 = rnorm(x0, sqrt(t) * σa)
        return sTfpe(t, false, false, x0, x1, false, false), na, nn, lr
      else
        return sTfpe(), na, nn, NaN
      end
    end

    x1 = rnorm(x0, sqrt(tw) * σa)

    # if speciation
    if λevent(λ, μ, ψi)
      nn += 1
      xk  = rnorm(x1, σk)
      shi = rand(Bool)
      xl, xr = if shi xk, x1 else x1, xk end

      d1, na, nn, lr = 
        _sim_cfpe_it(t - tw, λ, μ, ψ, xl, σa, σk, ψts, ix, nep, lr, lU, iρi, 
          na, nn, nlim)
      d2, na, nn, lr = 
        _sim_cfpe_it(t - tw, λ, μ, ψ, xr, σa, σk, ψts, ix, nep, lr, lU, iρi, 
          na, nn, nlim)

      return sTfpe(d1, d2, tw, false, false, x0, x1, shi, false), na, nn, lr
    # if extinction
    elseif μevent(μ, ψi)

      return sTfpe(tw, true, false, x0, x1, false, false), na, nn, lr
    # if fossil sample
    else
      return sTfpe(), na, nn, NaN
    end
  end

  return sTfpe(), na, nn, NaN
end



