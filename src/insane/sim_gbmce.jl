#=

Anagenetic GBM birth-death Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    sim_gbmce(t   ::Float64,
              λt  ::Float64,
              μ   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64)

Simulate `iTgbmce` according to a geometric Brownian motion.
"""
function sim_gbmce(t   ::Float64,
                   λt  ::Float64,
                   μ   ::Float64,
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

      if divev(λm, μ, t)
        # if speciation
        if λorμ(λm, μ)
          return iTgbmce(sim_gbmce(0.0, λt1, μ, σλ, δt, srδt), 
                         sim_gbmce(0.0, λt1, μ, σλ, δt, srδt), 
                  bt, δt, t, false, false, λv)
        # if extinction
        else
          return iTgbmce(bt, δt, t, true, false, λv)
        end
      end

      return iTgbmce(bt, δt, t, false, false, λv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divev(λm, μ, δt)
      # if speciation
      if λorμ(λm, μ)
        return iTgbmce(sim_gbmce(t, λt1, μ, σλ, δt, srδt), 
                       sim_gbmce(t, λt1, μ, σλ, δt, srδt), 
                bt, δt, δt, false, false, λv)
      # if extinction
      else
        return iTgbmce(bt, δt, δt, true, false, λv)
      end
    end

    λt = λt1
  end
end




"""
    sim_gbmce(t   ::Float64,
            λt  ::Float64,
            μ   ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64,
            nsp ::Int64,
            nlim::Int64)

Simulate `iTgbmce` according to a geometric Brownian motion with a limit
on the number lineages allowed to reach.
"""
function sim_gbmce(t   ::Float64,
                 λt  ::Float64,
                 μ   ::Float64,
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

        if divev(λm, μ, t)
          # if speciation
          if λorμ(λm, μ)
            nsp += 1
            td1, nsp = sim_gbmce(0.0, λt1, μ, σλ, δt, srδt, nsp, nlim)
            td2, nsp = sim_gbmce(0.0, λt1, μ, σλ, δt, srδt, nsp, nlim)

            return iTgbmce(td1, td2, bt, δt, t, false, false, λv), nsp
          # if extinction
          else
            return iTgbmce(bt, δt, t, true, false, λv), nsp
          end
        end

        return iTgbmce(bt, δt, t, false, false, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, μ, δt)
        # if speciation
        if λorμ(λm, μ)
          nsp += 1
          td1, nsp = sim_gbmce(t, λt1, μ, σλ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbmce(t, λt1, μ, σλ, δt, srδt, nsp, nlim)

          return iTgbmce(td1, td2, bt, δt, δt, false, false, λv), nsp
        # if extinction
        else
          return iTgbmce(bt, δt, δt, true, false, λv), nsp
        end
      end

      λt = λt1
    end

  else
    return iTgbmce(), nsp
  end
end





"""
    sim_gbmce(nsδt::Float64,
            t   ::Float64,
            λt  ::Float64,
            μ   ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Simulate `iTgbmce` according to a geometric Brownian motion starting 
with a non-standard δt.
"""
function sim_gbmce(nsδt::Float64,
                 t   ::Float64,
                 λt  ::Float64,
                 μ   ::Float64,
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

    if divev(λm, μ, t)
      # if speciation
      if λorμ(λm, μ)
        return iTgbmce(sim_gbmce(0.0, λt1, μ, σλ, δt, srδt), 
                       sim_gbmce(0.0, λt1, μ, σλ, δt, srδt), 
                bt, δt, t, false, false, λv)
      # if extinction
      else
        return iTgbmce(bt, δt, t, true, false, λv)
      end
    end

    return iTgbmce(bt, δt, t, false, false, λv)
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt, srnsδt*σλ)

  push!(λv, λt1)

  λm = exp(0.5*(λt + λt1))

  if divev(λm, μ, nsδt)
    # if speciation
    if λorμ(λm, μ)
      return iTgbmce(sim_gbmce(t, λt1, μ, σλ, δt, srδt), 
                     sim_gbmce(t, λt1, μ, σλ, δt, srδt), 
              bt, δt, nsδt, false, false, λv)
    # if extinction
    else
      return iTgbmce(bt, δt, nsδt, true, false, λv)
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

      if divev(λm, μ, t)
        # if speciation
        if λorμ(λm, μ)
          return iTgbmce(sim_gbmce(0.0, λt1, μ, σλ, δt, srδt), 
                         sim_gbmce(0.0, λt1, μ, σλ, δt, srδt), 
                  bt, δt, t, false, false, λv)
        # if extinction
        else
          return iTgbmce(bt, δt, t, true, false, λv)
        end
      end

      return iTgbmce(bt, δt, t, false, false, λv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divev(λm, μ, δt)
      # if speciation
      if λorμ(λm, μ)
        return iTgbmce(sim_gbmce(t, λt1, μ, σλ, δt, srδt), 
                       sim_gbmce(t, λt1, μ, σλ, δt, srδt), 
                bt, δt, δt, false, false, λv)
      # if extinction
      else
        return iTgbmce(bt, δt, δt, true, false, λv)
      end
    end

    λt = λt1
  end
end





"""
    sim_gbmce(nsδt::Float64,
            t   ::Float64,
            λt  ::Float64,
            μ   ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64, 
            nsp ::Int64,
            nlim::Int64)

Simulate `iTgbmce` according to a geometric Brownian motion starting 
with a non-standard δt with a limit in the number of species.
"""
function sim_gbmce(nsδt::Float64,
                   t   ::Float64,
                   λt  ::Float64,
                   μ   ::Float64,
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

    if divev(λm, μ, t)
      # if speciation
      if λorμ(λm, μ)
        nsp += 1
        td1, nsp = sim_gbmce(0.0, λt1, μ, σλ, δt, srδt, nsp, nlim)
        td2, nsp = sim_gbmce(0.0, λt1, μ, σλ, δt, srδt, nsp, nlim)

        return iTgbmce(td1, td2, bt, δt, t, false, false, λv), nsp
      # if extinction
      else
        return iTgbmce(bt, δt, t, true, false, λv), nsp
      end
    end

    return iTgbmce(bt, δt, t, false, false, λv), nsp
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt, srnsδt*σλ)

  push!(λv, λt1)

  λm = exp(0.5*(λt + λt1))

  if divev(λm, μ, nsδt)
    # if speciation
    if λorμ(λm, μ)
      nsp += 1
      td1, nsp = sim_gbmce(t, λt1, μ, σλ, δt, srδt, nsp, nlim)
      td2, nsp = sim_gbmce(t, λt1, μ, σλ, δt, srδt, nsp, nlim)

      return iTgbmce(td1, td2, bt, δt, nsδt, false, false, λv), nsp
    else
    # if extinction
      return iTgbmce(bt, δt, nsδt, true, false, λv), nsp
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

        if divev(λm, μ, t)
          # if speciation
          if λorμ(λm, μ)
            nsp += 1
            td1, nsp = sim_gbmce(0.0, λt1, μ, σλ, δt, srδt, nsp, nlim)
            td2, nsp = sim_gbmce(0.0, λt1, μ, σλ, δt, srδt, nsp, nlim)

            return iTgbmce(td1, td2, bt, δt, t, false, false, λv), nsp
          # if extinction
          else
            return iTgbmce(bt, δt, t, true, false, λv), nsp
          end
        end

        return iTgbmce(bt, δt, t, false, false, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, μ, δt)
        # if speciation
        if λorμ(λm, μ)
          nsp += 1
          td1, nsp = sim_gbmce(t, λt1, μ, σλ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbmce(t, λt1, μ, σλ, δt, srδt, nsp, nlim)

          return iTgbmce(td1, td2, bt, δt, δt, false, false, λv), nsp
        # if extinction
        else
          return iTgbmce(bt, δt, δt, true, false, λv), nsp
        end
      end

      λt = λt1
    end
  else
    return iTgbmce(), nsp
  end
end





"""
    divev(λ::Float64, μ::Float64, δt::Float64)
Return true if diversification event.
"""
divev(λ::Float64, μ::Float64, δt::Float64) = @fastmath rand() < (λ + μ)*δt 




"""
    rnorm(μ::Float64, σ::Float64)
Generate a normal variable with mean `μ` and variance `σ`.
"""
rnorm(μ::Float64, σ::Float64) = @fastmath randn()*σ + μ



