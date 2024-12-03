#=

constant punkeek simulation

Ignacio Quintero Mächler

t(-_-t)

Created 20 11 2024
=#


"""
    sim_cpe(t ::Float64,
            λ ::Float64,
            μ ::Float64,
            x0::Float64,
            σa::Float64,
            σk::Float64)

Simulate a constant punkeek model of height `t` with speciation rate `λ`,
extinction rate `μ`, starting trait value `x0`, and anagenetic and 
cladogenetic variance `σa` and `σk`.
"""
function sim_cpe(t ::Float64,
                 λ ::Float64,
                 μ ::Float64,
                 x0::Float64,
                 σa::Float64,
                 σk::Float64)

  tw = cbd_wait(λ, μ)

  if tw > t
    x1 = rnorm(x0, sqrt(t) * σa)
    return sTpe(t, false, x0, x1, false, false)
  end

  x1 = rnorm(x0, sqrt(tw) * σa)

  if λorμ(λ, μ)
    xk = rnorm(x1, σk)
    xl, xr = if rand() < 0.5 xk, x1 else x1, xk end

    return sTpe(sim_cpe(t - tw, λ, μ, xl, σa, σk), 
                sim_cpe(t - tw, λ, μ, xr, σa, σk), 
                tw, false, x0, x1, true, false)
  else
    return sTpe(tw, true, x0, x1, false, false)
  end
end




"""
    _sim_cpe_t(t   ::Float64,
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
function _sim_cpe_t(t   ::Float64,
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

  if isfinite(lr) && nn < nlim

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(iρi * Float64(na)/Float64(na-1))
      end
      if nlr < lr && lU >= nlr
        return sTpe(), na, nn, NaN
      else
        x1 = rnorm(x0, sqrt(t) * σa)
        push!(xist, x0)
        push!(xfst, x1)
        push!(est, t)
        return sTpe(t, false, x0, x1, false, false), na, nn, nlr
      end
    end

    x1 = rnorm(x0, sqrt(tw) * σa)

    if λorμ(λ, μ)
      nn += 1
      xk = rnorm(x1, σk)
      xl, xr = if rand() < 0.5 xk, x1 else x1, xk end

      d1, na, nn, lr = 
        _sim_cpe_t(t - tw, λ, μ, xl, σa, σk, lr, lU, iρi, na, nn, nlim, 
          xist, xfst, est)
      d2, na, nn, lr = 
        _sim_cpe_t(t - tw, λ, μ, xr, σa, σk, lr, lU, iρi, na, nn, nlim, 
          xist, xfst, est)

      return sTpe(d1, d2, tw, false, x0, x1, true, false), na, nn, lr
    else
      return sTpe(tw, true, x0, x1, false, false), na, nn, lr
    end
  end

  return sTpe(), na, nn, NaN
end





"""
    _sim_cpe_i(t   ::Float64,
               λ   ::Float64,
               μ   ::Float64,
               x0  ::Float64,
               σa  ::Float64,
               σk  ::Float64,
               na  ::Int64,
               nn  ::Int64,
               nlim::Int64)

Simulate a constant punkeek model of height `t` with speciation rate `λ`,
extinction rate `μ`, starting trait value `x0`, and anagenetic and 
cladogenetic variance `σa` and `σk` for internal branches.
"""
function _sim_cpe_i(t   ::Float64,
                    λ   ::Float64,
                    μ   ::Float64,
                    x0  ::Float64,
                    σa  ::Float64,
                    σk  ::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

  if nn < nlim

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      x1 = rnorm(x0, sqrt(t) * σa)
      return sTpe(t, false, x0, x1, false, false), na, nn
    end

    x1 = rnorm(x0, sqrt(tw) * σa)

    if λorμ(λ, μ)
      nn += 1
      xk = rnorm(x1, σk)
      xl, xr = if rand() < 0.5 xk, x1 else x1, xk end

      d1, na, nn = _sim_cpe_i(t - tw, λ, μ, xl, σa, σk, na, nn, nlim)
      d2, na, nn = _sim_cpe_i(t - tw, λ, μ, xr, σa, σk, na, nn, nlim)

      return sTpe(d1, d2, tw, false, x0, x1, true, false), na, nn
    else
      return sTpe(tw, true, x0, x1, false, false), na, nn
    end
  end

  return sTpe(), na, nn
end




"""
    _sim_cpe_it(t   ::Float64,
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
function _sim_cpe_it(t   ::Float64,
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

  if lU < lr && nn < nlim

    tw = cbd_wait(λ, μ)

    if tw > t
      na += 1
      lr += log(iρi)
      x1  = rnorm(x0, sqrt(t) * σa)
      return sTpe(t, false, x0, x1, false, false), na, nn, lr
    end

    x1 = rnorm(x0, sqrt(tw) * σa)

    if λorμ(λ, μ)
      nn += 1
      xk = rnorm(x1, σk)
      xl, xr = if rand() < 0.5 xk, x1 else x1, xk end

      d1, na, nn, lr = 
        _sim_cpe_it(t - tw, λ, μ, xl, σa, σk, lr, lU, iρi, na, nn, nlim)
      d2, na, nn, lr = 
        _sim_cpe_it(t - tw, λ, μ, xr, σa, σk, lr, lU, iρi, na, nn, nlim)

      return sTpe(d1, d2, tw, false, x0, x1, true, false), na, nn, lr
    else
      return sTpe(tw, true, x0, x1, false, false), na, nn, lr
    end
  end

  return sTpe(), na, nn, NaN
end



