#=

Anagenetic GBM birth-death Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sim_gbmct(t   ::Float64;
              λ0  ::Float64 = 1.0,
              α   ::Float64 = 0.0,
              σλ  ::Float64 = 0.1,
              ϵ   ::Float64 = 0.2,
              δt  ::Float64 = 1e-3,
              nlim::Int64   = 10_000,
              init::Symbol  = :crown)

Simulate `iTgbmpb` according to a pure-birth geometric Brownian motion.
"""
function sim_gbmct(t   ::Float64;
                   λ0  ::Float64 = 1.0,
                   α   ::Float64 = 0.0,
                   σλ  ::Float64 = 0.1,
                   ϵ   ::Float64 = 0.2,
                   δt  ::Float64 = 1e-3,
                   nlim::Int64   = 10_000,
                   init::Symbol  = :crown)

  if init === :crown
    lλ0 = log(λ0)
    d1, nsp = _sim_gbmct(t, lλ0, α, σλ, ϵ, δt, sqrt(δt), 1, nlim)
    if nsp >= nlim 
      @warn "maximum number of lineages surpassed"
    end

    d2, nsp = _sim_gbmct(t, lλ0, α, σλ, ϵ, δt, sqrt(δt), 1, nlim)
    if nsp >= nlim 
      @warn "maximum number of lineages surpassed"
    end

    tree = iTgbmct(d1, d2, 0.0, δt, 0.0, false, false, Float64[lλ0, lλ0])
  elseif init === :stem
    tree, nsp = _sim_gbmct(t, log(λ0), α, σλ, ϵ, δt, sqrt(δt), 1, nlim)

    if nsp >= nlim 
      @warn "maximum number of lineages surpassed"
    end
  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_gbmct(t   ::Float64,
               λt  ::Float64,
               α   ::Float64,
               σλ  ::Float64,
               ϵ   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Simulate `iTgbmct` according to a geometric Brownian motion.
"""
function _sim_gbmct(t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  λv = Float64[λt]
  bt = 0.0

  while true

    if t <= δt
      bt  += t

      t = max(0.0,t)
      srt = sqrt(t)
      λt1 = rnorm(λt + α*t, srt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, t)
        # if speciation
        if λorμ(λm, ϵ*λm)
          return iTgbmct(iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                         iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                  bt, δt, t, false, false, λv)
        # if extinction
        else
          return iTgbmct(bt, δt, t, true, false, λv)
        end
      end

      return iTgbmct(bt, δt, t, false, false, λv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt + α*δt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, δt)
      # if speciation
      if λorμ(λm, ϵ*λm)
        return iTgbmct(_sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt), 
                       _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt), 
                bt, δt, δt, false, false, λv)
      # if extinction
      else
        return iTgbmct(bt, δt, δt, true, false, λv)
      end
    end

    λt = λt1
  end
end




"""
    _sim_gbmct(t   ::Float64,
               λt  ::Float64,
               α   ::Float64,
               σλ  ::Float64,
               ϵ   ::Float64,
               δt  ::Float64,
               srδt::Float64,
               nsp ::Int64,
               nlim::Int64)

Simulate `iTgbmct` according to a geometric Brownian motion with a limit
on the number lineages allowed to reach.
"""
function _sim_gbmct(t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    nsp ::Int64,
                    nlim::Int64)

  if nsp < nlim

    λv = Float64[λt]
    bt = 0.0

    while true

      if t <= δt
        bt  += t

        t = max(0.0,t)
        srt = sqrt(t)
        λt1 = rnorm(λt + α*t, srt*σλ)

        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divevϵ(λm, ϵ, t)
          # if speciation
          if λorμ(λm, ϵ*λm)
            nsp += 1
            return iTgbmct(
                    iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                    iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                    bt, δt, t, false, false, λv), nsp
          # if extinction
          else
            return iTgbmct(bt, δt, t, true, false, λv), nsp
          end
        end

        return iTgbmct(bt, δt, t, false, false, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(λm, ϵ*λm)
          nsp += 1
          td1, nsp = _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt, nsp, nlim)
          td2, nsp = _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt, nsp, nlim)

          return iTgbmct(td1, td2, bt, δt, δt, false, false, λv), nsp
        # if extinction
        else
          return iTgbmct(bt, δt, δt, true, false, λv), nsp
        end
      end

      λt = λt1
    end

  else
    return iTgbmct(), nsp
  end
end




"""
    _sim_gbmct(nsδt::Float64,
               t   ::Float64,
               λt  ::Float64,
               α   ::Float64,
               σλ  ::Float64,
               ϵ   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Simulate `iTgbmct` according to a geometric Brownian motion starting 
with a non-standard δt.
"""
function _sim_gbmct(nsδt::Float64,
                    t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  λv = Float64[λt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt  += t

    t   = max(0.0,t)
    srt = sqrt(t)
    λt1 = rnorm(λt + α*t, srt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, t)
      # if speciation
      if λorμ(λm, ϵ*λm)
        return iTgbmct(iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                       iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                bt, δt, t, false, false, λv)
      # if extinction
      else
        return iTgbmct(bt, δt, t, true, false, λv)
      end
    end

    return iTgbmct(bt, δt, t, false, false, λv)
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt + α*nsδt, srnsδt*σλ)

  push!(λv, λt1)

  λm = exp(0.5*(λt + λt1))

  if divevϵ(λm, ϵ, nsδt)
    # if speciation
    if λorμ(λm, ϵ*λm)
      return iTgbmct(_sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt), 
                     _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt), 
              bt, δt, nsδt, false, false, λv)
    # if extinction
    else
      return iTgbmct(bt, δt, nsδt, true, false, λv)
    end
  end

  λt = λt1

  ## second: standard δt
  while true

    if t <= δt
      bt  += t

      t   = max(0.0,t)
      srt = sqrt(t)
      λt1 = rnorm(λt + α*t, srt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, t)
        # if speciation
        if λorμ(λm, ϵ*λm)
          return iTgbmct(iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                         iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                  bt, δt, t, false, false, λv)
        # if extinction
        else
          return iTgbmct(bt, δt, t, true, false, λv)
        end
      end

      return iTgbmct(bt, δt, t, false, false, λv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt + α*δt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, δt)
      # if speciation
      if λorμ(λm, ϵ*λm)
        return iTgbmct(_sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt), 
                       _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt), 
                bt, δt, δt, false, false, λv)
      # if extinction
      else
        return iTgbmct(bt, δt, δt, true, false, λv)
      end
    end

    λt = λt1
  end
end





"""
    _sim_gbmct(nsδt::Float64,
               t   ::Float64,
               λt  ::Float64,
               α   ::Float64,
               σλ  ::Float64,
               ϵ   ::Float64,
               δt  ::Float64,
               srδt::Float64, 
               nsp ::Int64,
               nlim::Int64)

Simulate `iTgbmct` according to a geometric Brownian motion starting 
with a non-standard δt with a limit in the number of species.
"""
function _sim_gbmct(nsδt::Float64,
                    t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64, 
                    nsp ::Int64,
                    nlim::Int64)

  λv = Float64[λt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt  += t

    t   = max(0.0, t)
    srt = sqrt(t)
    λt1 = rnorm(λt + α*t, srt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, t)
      # if speciation
      if λorμ(λm, ϵ*λm)
        nsp += 1

        return iTgbmct(
                 iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                 iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                 bt, δt, t, false, false, λv), nsp
      # if extinction
      else
        return iTgbmct(bt, δt, t, true, false, λv), nsp
      end
    end

    return iTgbmct(bt, δt, t, false, false, λv), nsp
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt + α*nsδt, srnsδt*σλ)

  push!(λv, λt1)

  λm = exp(0.5*(λt + λt1))

  if divevϵ(λm, ϵ, nsδt)
    # if speciation
    if λorμ(λm, ϵ*λm)
      nsp += 1
      td1, nsp = _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt, nsp, nlim)
      td2, nsp = _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt, nsp, nlim)

      return iTgbmct(td1, td2, bt, δt, nsδt, false, false, λv), nsp
    else
    # if extinction
      return iTgbmct(bt, δt, nsδt, true, false, λv), nsp
    end
  end

  λt = λt1

  if nsp < nlim

    ## second: standard δt
    while true

      if t <= δt
        bt  += t

        t   = max(0.0,t)
        srt = sqrt(t)
        λt1 = rnorm(λt + α*t, srt*σλ)

        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divevϵ(λm, ϵ, t)
          # if speciation
          if λorμ(λm, ϵ*λm)
            nsp += 1
            return iTgbmct(
                      iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                      iTgbmct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                           bt, δt, t, false, false, λv), nsp
          # if extinction
          else
            return iTgbmct(bt, δt, t, true, false, λv), nsp
          end
        end

        return iTgbmct(bt, δt, t, false, false, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(λm, ϵ*λm)
          nsp += 1
          td1, nsp = _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt, nsp, nlim)
          td2, nsp = _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt, nsp, nlim)

          return iTgbmct(td1, td2, bt, δt, δt, false, false, λv), nsp
        # if extinction
        else
          return iTgbmct(bt, δt, δt, true, false, λv), nsp
        end
      end

      λt = λt1
    end
  else
    return iTgbmct(), nsp
  end
end





"""
    divevϵ(λ::Float64, ϵ::Float64, δt::Float64)

Return true if diversification event for `ϵ` parametization.
"""
divevϵ(λ::Float64, ϵ::Float64, δt::Float64) = @fastmath rand() < (1.0 + ϵ)*λ*δt 




