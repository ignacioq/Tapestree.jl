#=

Diffused Brownian motion simulation on fixed trees

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#






"""
"""
function sim_dbm(tree::iTree, 
                 x0  ::Float64,
                 σ0  ::Float64,
                 γ   ::Float64,
                 δt  ::Float64)
  _sim_dbm(tree, x0, log(σ0), γ, δt, sqrt(δt))
end


"""
"""
function _sim_dbm(tree::iTree, 
                  x0  ::Float64,
                  lσ0 ::Float64,
                  γ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  et = e(tree)

  # simulate dbm
  if iszero(et)
    xv   = Float64[x0,  x0]
    lσ   = Float64[lσ0, lσ0]
    fdti = 0.0
    l    = 2
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    xv, lσ  = dbm(x0, lσ0, γ, fdti, srδt, nt)
    l   = nt + 2
  end

  if def1(tree)
    if def2(tree)
      x0  = xv[end]
      lσ0 = lσ[end]
      iTxs(_sim_dbm(tree.d1, x0, lσ0, γ, δt, srδt), 
           _sim_dbm(tree.d2, x0, lσ0, γ, δt, srδt), 
           et, δt, fdti, xv, lσ)
    else
      iTxs(_sim_dbm(tree.d1, xv[end], lσ[end], γ, δt, srδt), 
           et, δt, fdti, xv, lσ)
    end
  else
    iTxs(et, δt, fdti, xv, lσ)
  end
end



































"""
    sim_gbmpb(t   ::Float64;
              λ0  ::Float64 = 1.0,
              α   ::Float64 = 0.0,
              σλ  ::Float64 = 0.1,
              δt  ::Float64 = 1e-3,
              nlim::Int64   = 10_000,
              init::Symbol  = :crown)

Simulate `iTpb` according to a pure-birth geometric Brownian motion
conditional in stopping at time `t`.
"""
function sim_gbmpb(t   ::Float64;
                   λ0  ::Float64 = 1.0,
                   α   ::Float64 = 0.0,
                   σλ  ::Float64 = 0.1,
                   δt  ::Float64 = 1e-3,
                   nlim::Int64   = 10_000,
                   init::Symbol  = :crown)

  if init === :crown
    lλ0 = log(λ0)
    d1, nn = _sim_gbmpb(t, lλ0, α, σλ, δt, sqrt(δt), 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    d2, nn = _sim_gbmpb(t, lλ0, α, σλ, δt, sqrt(δt), 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    tree = iTpb(d1, d2, 0.0, δt, 0.0, false, [lλ0, lλ0])
   elseif init === :stem
    tree, nn = _sim_gbmpb(t, log(λ0), α, σλ, δt, sqrt(δt), nn + 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_gbmpb(t   ::Float64,
               λt  ::Float64,
               α   ::Float64,
               σλ  ::Float64,
               δt  ::Float64,
               srδt::Float64,
               nn ::Int64,
               nlim::Int64)

Simulate `iTpb` according to a pure-birth geometric Brownian motion.
"""
function _sim_gbmpb(t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    nn ::Int64,
                    nlim::Int64)

  if nn < nlim

    λv = Float64[λt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, t)
          nn += 1
          return iTpb(iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                      iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                      bt, δt, t, false, λv), nn
        end

        return iTpb(bt, δt, t, false, λv), nn
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nn += 1
        td1, nn = _sim_gbmpb(t, λt1, α, σλ, δt, srδt, nn, nlim)
        td2, nn = _sim_gbmpb(t, λt1, α, σλ, δt, srδt, nn, nlim)

        return iTpb(td1, td2, bt, δt, δt, false, λv), nn
      end

      λt = λt1
    end
  end

  return iTpb(), nn
end





"""
    _sim_gbmpb_t(t   ::Float64,
                 λt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 Iρi ::Float64,
                 na  ::Int64,
                 nn ::Int64,
                 nlim::Int64)

Simulate `iTpb` according to a pure-birth geometric Brownian motion for
terminal branches.
"""
function _sim_gbmpb_t(t   ::Float64,
                      λt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      lr  ::Float64,
                      lU  ::Float64,
                      Iρi ::Float64,
                      na  ::Int64,
                      nn ::Int64,
                      nlim::Int64)

  if isfinite(lr) && nn < nlim

    λv = Float64[λt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, t)
          nn += 1
          na += 2
          if na === 2
            nlr = lr + log(Iρi*2.0)
          else
            nlr = lr + log(Iρi * Iρi * Float64(na)/Float64(na-2))
          end
          if nlr < lr && lU >= nlr
            return iTpb(), na, nn, NaN
          else
            return iTpb(iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                        iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                        bt, δt, t, false, λv), na, nn, nlr
          end
        else
          na += 1
          nlr = lr
          if na > 1
            nlr += log(Iρi * Float64(na)/Float64(na-1))
          end
          if nlr >= lr
            return iTpb(bt, δt, t, false, λv), na, nn, nlr
          elseif lU < nlr
            return iTpb(bt, δt, t, false, λv), na, nn, nlr
          else
            return iTpb(), na, nn, NaN
          end
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nn += 1
        td1, na, nn, lr =
          _sim_gbmpb_t(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, na, nn, nlim)
        td2, na, nn, lr =
          _sim_gbmpb_t(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, na, nn, nlim)

        return iTpb(td1, td2, bt, δt, δt, false, λv), na, nn, lr
      end

      λt = λt1
    end
  end

  return iTpb(), na, nn, NaN
end




"""
    _sim_gbmpb_it(nsδt::Float64,
                  t   ::Float64,
                  λt  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  lr  ::Float64,
                  lU  ::Float64,
                  Iρi ::Float64,
                  na  ::Int64,
                  nn ::Int64,
                  nlim::Int64)
Simulate `iTpb` according to a pure-birth geometric Brownian motion,
starting with a non-standard δt with a limit in the number of species.
"""
function _sim_gbmpb_it(nsδt::Float64,
                       t   ::Float64,
                       λt  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       Iρi ::Float64,
                       nn  ::Int64,
                       nlim::Int64)

  λv = Float64[λt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt + accerr
    t   = isapprox(t, 0.0) ? 0.0 : isapprox(t, nsδt) ? nsδt : t
    bt += t
    λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
    λm  = exp(0.5*(λt + λt1))
    push!(λv, λt1)

    if divev(λm, t)
      nn += 1
      lr += 2.0*log(Iρi)
      return iTpb(iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                  iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                  bt, δt, t, false, λv), nn, lr
    else
      lr += log(Iρi)
      return iTpb(bt, δt, t, false, λv), nn, lr
    end
  end

  t  -= nsδt
  bt += nsδt

  λt1 = rnorm(λt + α*nsδt, sqrt(nsδt)*σλ)
  λm  = exp(0.5*(λt + λt1))
  push!(λv, λt1)

  if divev(λm, nsδt)
    nn += 1
    td1, nn, lr =
      _sim_gbmpb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)
    td2, nn, lr =
      _sim_gbmpb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)

    return iTpb(td1, td2, bt, δt, nsδt, false, λv), nn, lr
  end

  λt = λt1

  if lU < lr && nn < nlim

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, t)
          nn += 1
          lr  += 2.0*log(Iρi)
          return iTpb(iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                      iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                      bt, δt, t, false, λv), nn, lr
        else
          lr += log(Iρi)
          return iTpb(bt, δt, t, false, λv), nn, lr
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nn += 1
        td1, nn, lr =
          _sim_gbmpb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)
        td2, nn, lr =
          _sim_gbmpb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)

        return iTpb(td1, td2, bt, δt, δt, false, λv), nn, lr
      end

      λt = λt1
    end
  end

  return iTpb(), nn, NaN
end




"""
    _sim_gbmpb_it(t   ::Float64,
                 λt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 Iρi ::Float64,
                 nn ::Int64,
                 nlim::Int64)

Simulate `iTpb` according to a pure-birth geometric Brownian motion for
terminal branches.
"""
function _sim_gbmpb_it(t   ::Float64,
                       λt  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       Iρi ::Float64,
                       nn ::Int64,
                       nlim::Int64)

  if lU < lr && nn < nlim

    λv = Float64[λt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, t)
          nn += 1
          lr  += 2.0*log(Iρi)
          return iTpb(iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                      iTpb(0.0, δt, 0.0, false, [λt1, λt1]),
                      bt, δt, t, false, λv), nn, lr
        end

        lr += log(Iρi)
        return iTpb(bt, δt, t, false, λv), nn, lr
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nn += 1
        td1, nn, lr =
          _sim_gbmpb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)
        td2, nn, lr =
          _sim_gbmpb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)

        return iTpb(td1, td2, bt, δt, δt, false, λv), nn, lr
      end

      λt = λt1
    end
  end

  return iTpb(), nn, NaN
end




"""
    divev(λ::Float64, δt::Float64)

Return true if diversification event.
"""
divev(λ::Float64, δt::Float64) = @fastmath rand() < λ*δt


