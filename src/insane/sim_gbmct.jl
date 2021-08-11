#=

Anagenetic GBM birth-death Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    sim_gbmct(t   ::Float64,
              λt  ::Float64,
              ϵ   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64)

Simulate `iTgbmct` according to a geometric Brownian motion.
"""
function sim_gbmct(t   ::Float64,
                   λt  ::Float64,
                   ϵ   ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

  λv = Float64[λt]
  bt = 0.0

  while true

    if t <= δt
      bt  += t

      t = max(0.0,t)
      srt = sqrt(t)
      λt1 = rnorm(λt, srt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, t)
        # if speciation
        if λorμ(λm, ϵ*λm)
          return iTgbmct(sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt), 
                         sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt), 
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

    λt1 = rnorm(λt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, δt)
      # if speciation
      if λorμ(λm, ϵ*λm)
        return iTgbmct(sim_gbmct(t, λt1, ϵ, σλ, δt, srδt), 
                       sim_gbmct(t, λt1, ϵ, σλ, δt, srδt), 
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
    sim_gbmct(t   ::Float64,
            λt  ::Float64,
            ϵ   ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64,
            nsp ::Int64,
            nlim::Int64)

Simulate `iTgbmct` according to a geometric Brownian motion with a limit
on the number lineages allowed to reach.
"""
function sim_gbmct(t   ::Float64,
                   λt  ::Float64,
                   ϵ   ::Float64,
                   σλ  ::Float64,
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
        λt1 = rnorm(λt, srt*σλ)

        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divevϵ(λm, ϵ, t)
          # if speciation
          if λorμ(λm, ϵ*λm)
            nsp += 1
            td1, nsp = sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt, nsp, nlim)
            td2, nsp = sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt, nsp, nlim)

            return iTgbmct(td1, td2, bt, δt, t, false, false, λv), nsp
          # if extinction
          else
            return iTgbmct(bt, δt, t, true, false, λv), nsp
          end
        end

        return iTgbmct(bt, δt, t, false, false, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(λm, ϵ*λm)
          nsp += 1
          td1, nsp = sim_gbmct(t, λt1, ϵ, σλ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbmct(t, λt1, ϵ, σλ, δt, srδt, nsp, nlim)

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
    sim_gbmct(nsδt::Float64,
            t   ::Float64,
            λt  ::Float64,
            ϵ   ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Simulate `iTgbmct` according to a geometric Brownian motion starting 
with a non-standard δt.
"""
function sim_gbmct(nsδt::Float64,
                   t   ::Float64,
                   λt  ::Float64,
                   ϵ   ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

  λv = Float64[λt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt  += t

    t   = max(0.0,t)
    srt = sqrt(t)
    λt1 = rnorm(λt, srt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, t)
      # if speciation
      if λorμ(λm, ϵ*λm)
        return iTgbmct(sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt), 
                       sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt), 
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

  λt1 = rnorm(λt, srnsδt*σλ)

  push!(λv, λt1)

  λm = exp(0.5*(λt + λt1))

  if divevϵ(λm, ϵ, nsδt)
    # if speciation
    if λorμ(λm, ϵ*λm)
      return iTgbmct(sim_gbmct(t, λt1, ϵ, σλ, δt, srδt), 
                     sim_gbmct(t, λt1, ϵ, σλ, δt, srδt), 
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
      λt1 = rnorm(λt, srt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, t)
        # if speciation
        if λorμ(λm, ϵ*λm)
          return iTgbmct(sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt), 
                         sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt), 
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

    λt1 = rnorm(λt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, δt)
      # if speciation
      if λorμ(λm, ϵ*λm)
        return iTgbmct(sim_gbmct(t, λt1, ϵ, σλ, δt, srδt), 
                       sim_gbmct(t, λt1, ϵ, σλ, δt, srδt), 
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
    sim_gbmct(nsδt::Float64,
            t   ::Float64,
            λt  ::Float64,
            ϵ   ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64, 
            nsp ::Int64,
            nlim::Int64)

Simulate `iTgbmct` according to a geometric Brownian motion starting 
with a non-standard δt with a limit in the number of species.
"""
function sim_gbmct(nsδt::Float64,
                   t   ::Float64,
                   λt  ::Float64,
                   ϵ   ::Float64,
                   σλ  ::Float64,
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
    λt1 = rnorm(λt, srt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, t)
      # if speciation
      if λorμ(λm, ϵ*λm)
        nsp += 1
        td1, nsp = sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt, nsp, nlim)
        td2, nsp = sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt, nsp, nlim)

        return iTgbmct(td1, td2, bt, δt, t, false, false, λv), nsp
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

  λt1 = rnorm(λt, srnsδt*σλ)

  push!(λv, λt1)

  λm = exp(0.5*(λt + λt1))

  if divevϵ(λm, ϵ, nsδt)
    # if speciation
    if λorμ(λm, ϵ*λm)
      nsp += 1
      td1, nsp = sim_gbmct(t, λt1, ϵ, σλ, δt, srδt, nsp, nlim)
      td2, nsp = sim_gbmct(t, λt1, ϵ, σλ, δt, srδt, nsp, nlim)

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
        λt1 = rnorm(λt, srt*σλ)

        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divevϵ(λm, ϵ, t)
          # if speciation
          if λorμ(λm, ϵ*λm)
            nsp += 1
            td1, nsp = sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt, nsp, nlim)
            td2, nsp = sim_gbmct(0.0, λt1, ϵ, σλ, δt, srδt, nsp, nlim)

            return iTgbmct(td1, td2, bt, δt, t, false, false, λv), nsp
          # if extinction
          else
            return iTgbmct(bt, δt, t, true, false, λv), nsp
          end
        end

        return iTgbmct(bt, δt, t, false, false, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(λm, ϵ*λm)
          nsp += 1
          td1, nsp = sim_gbmct(t, λt1, ϵ, σλ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbmct(t, λt1, ϵ, σλ, δt, srδt, nsp, nlim)

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




"""
    rnorm(μ::Float64, σ::Float64)

Generate a normal variable with mean `μ` and standard deviation `σ`.
"""
rnorm(μ::Float64, σ::Float64) = @fastmath randn()*σ + μ



