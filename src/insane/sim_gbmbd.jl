#=

Anagenetic GBM birth-death Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    sim_gbmbd(t   ::Float64,
              λt  ::Float64,
              μt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64)

Simulate `iTgbmbd` according to a `gbmbd`.
"""
function sim_gbmbd(t   ::Float64,
                   λt  ::Float64,
                   μt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

  λv = Float64[λt]
  μv = Float64[μt]
  bt = 0.0

  while true

    if t <= δt
      bt  += t

      t = max(0.0,t)
      srt = sqrt(t)
      λt1 = rnorm(λt + α*t, srt*σλ)
      μt1 = rnorm(μt, srt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, t)
        # if speciation
        if λorμ(λm, μm)
          return iTgbmbd(iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                             Float64[μt1, μt1]), 
                         iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                             Float64[μt1, μt1]), 
                  bt, δt, t, false, false, λv, μv)
        # if extinction
        else
          return iTgbmbd(bt, δt, t, true, false, λv, μv)
        end
      end

      return iTgbmbd(bt, δt, t, false, false, λv, μv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt + α*δt, srδt*σλ)
    μt1 = rnorm(μt, srδt*σμ)

    push!(λv, λt1)
    push!(μv, μt1)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, δt)
      # if speciation
      if λorμ(λm, μm)
        return iTgbmbd(sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt), 
                       sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt), 
                bt, δt, δt, false, false, λv, μv)
      # if extinction
      else
        return iTgbmbd(bt, δt, δt, true, false, λv, μv)
      end
    end

    λt = λt1
    μt = μt1
  end
end




"""
    sim_gbmbd(t   ::Float64,
              λt  ::Float64,
              μt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              nsp ::Int64,
              nlim::Int64)

Simulate `iTgbmbd` according to a geometric Brownian motion with a limit
on the number lineages allowed to reach.
"""
function sim_gbmbd(t   ::Float64,
                   λt  ::Float64,
                   μt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   nsp ::Int64,
                   nlim::Int64)

  if nsp < nlim

    λv = Float64[λt]
    μv = Float64[μt]
    bt = 0.0

    while true

      if t <= δt
        bt  += t

        t = max(0.0,t)
        srt = sqrt(t)
        λt1 = rnorm(λt + α*t, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, t)
          # if speciation
          if λorμ(λm, μm)
            nsp += 1
            return iTgbmbd(iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                               Float64[μt1, μt1]), 
                           iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                               Float64[μt1, μt1]), 
                           bt, δt, t, false, false, λv, μv), nsp
          # if extinction
          else
            return iTgbmbd(bt, δt, t, true, false, λv, μv), nsp
          end
        end

        return iTgbmbd(bt, δt, t, false, false, λv, μv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        # if speciation
        if λorμ(λm, μm)
          nsp += 1
          td1, nsp = sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt, nsp, nlim)

          return iTgbmbd(td1, td2, bt, δt, δt, false, false, λv, μv), nsp
        # if extinction
        else
          return iTgbmbd(bt, δt, δt, true, false, λv, μv), nsp
        end
      end

      λt = λt1
      μt = μt1
    end

  else
    return iTgbmbd(), nsp
  end
end





"""
    sim_gbmbd(nsδt::Float64,
              t   ::Float64,
              λt  ::Float64,
              μt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64)

Simulate `iTgbmbd` according to a geometric Brownian motion starting 
with a non-standard δt.
"""
function sim_gbmbd(nsδt::Float64,
                   t   ::Float64,
                   λt  ::Float64,
                   μt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

  λv = Float64[λt]
  μv = Float64[μt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt  += t

    t   = max(0.0,t)
    srt = sqrt(t)
    λt1 = rnorm(λt + α*t, srt*σλ)
    μt1 = rnorm(μt, srt*σμ)

    push!(λv, λt1)
    push!(μv, μt1)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, t)
      # if speciation
      if λorμ(λm, μm)
        return iTgbmbd(iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                           Float64[μt1, μt1]), 
                       iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                           Float64[μt1, μt1]), 
                       bt, δt, t, false, false, λv, μv)
      # if extinction
      else
        return iTgbmbd(bt, δt, t, true, false, λv, μv)
      end
    end

    return iTgbmbd(bt, δt, t, false, false, λv, μv)
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt + α*nsδt, srnsδt*σλ)
  μt1 = rnorm(μt, srnsδt*σμ)

  push!(λv, λt1)
  push!(μv, μt1)

  λm = exp(0.5*(λt + λt1))
  μm = exp(0.5*(μt + μt1))

  if divev(λm, μm, nsδt)
    # if speciation
    if λorμ(λm, μm)
      return iTgbmbd(sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt), 
                     sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt), 
              bt, δt, nsδt, false, false, λv, μv)
    # if extinction
    else
      return iTgbmbd(bt, δt, nsδt, true, false, λv, μv)
    end
  end

  λt = λt1
  μt = μt1

  ## second: standard δt
  while true

    if t <= δt
      bt  += t

      t   = max(0.0,t)
      srt = sqrt(t)
      λt1 = rnorm(λt + α*t, srt*σλ)
      μt1 = rnorm(μt, srt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, t)
        # if speciation
        if λorμ(λm, μm)
          return iTgbmbd(iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                             Float64[μt1, μt1]), 
                         iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                             Float64[μt1, μt1]), 
                  bt, δt, t, false, false, λv, μv)
        # if extinction
        else
          return iTgbmbd(bt, δt, t, true, false, λv, μv)
        end
      end

      return iTgbmbd(bt, δt, t, false, false, λv, μv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt + α*δt, srδt*σλ)
    μt1 = rnorm(μt, srδt*σμ)

    push!(λv, λt1)
    push!(μv, μt1)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, δt)
      # if speciation
      if λorμ(λm, μm)
        return iTgbmbd(sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt), 
                       sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt), 
                bt, δt, δt, false, false, λv, μv)
      # if extinction
      else
        return iTgbmbd(bt, δt, δt, true, false, λv, μv)
      end
    end

    λt = λt1
    μt = μt1
  end
end





"""
    sim_gbmbd(nsδt::Float64,
              t   ::Float64,
              λt  ::Float64,
              μt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64, 
              nsp ::Int64,
              nlim::Int64)

Simulate `iTgbmbd` according to a geometric Brownian motion starting 
with a non-standard δt with a limit in the number of species.
"""
function sim_gbmbd(nsδt::Float64,
                   t   ::Float64,
                   λt  ::Float64,
                   μt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64, 
                   nsp ::Int64,
                   nlim::Int64)

  λv = Float64[λt]
  μv = Float64[μt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt  += t

    t   = max(0.0,t)
    srt = sqrt(t)
    λt1 = rnorm(λt + α*t, srt*σλ)
    μt1 = rnorm(μt, srt*σμ)

    push!(λv, λt1)
    push!(μv, μt1)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, t)
      # if speciation
      if λorμ(λm, μm)
        nsp += 1

        return iTgbmbd(iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                           Float64[μt1, μt1]), 
                       iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                           Float64[μt1, μt1]), 
                       bt, δt, t, false, false, λv, μv), nsp
      # if extinction
      else
        return iTgbmbd(bt, δt, t, true, false, λv, μv), nsp
      end
    end

    return iTgbmbd(bt, δt, t, false, false, λv, μv), nsp
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt + α*nsδt, srnsδt*σλ)
  μt1 = rnorm(μt, srnsδt*σμ)

  push!(λv, λt1)
  push!(μv, μt1)

  λm = exp(0.5*(λt + λt1))
  μm = exp(0.5*(μt + μt1))

  if divev(λm, μm, nsδt)
    # if speciation
    if λorμ(λm, μm)
      nsp += 1
      td1, nsp = sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt, nsp, nlim)
      td2, nsp = sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt, nsp, nlim)

      return iTgbmbd(td1, td2, bt, δt, nsδt, false, false, λv, μv), nsp
    else
    # if extinction
      return iTgbmbd(bt, δt, nsδt, true, false, λv, μv), nsp
    end
  end

  λt = λt1
  μt = μt1

  if nsp < nlim

    ## second: standard δt
    while true

      if t <= δt
        bt  += t

        t   = max(0.0,t)
        srt = sqrt(t)
        λt1 = rnorm(λt + α*t, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, t)
          # if speciation
          if λorμ(λm, μm)
            nsp += 1

            return iTgbmbd(iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                               Float64[μt1, μt1]), 
                           iTgbmbd(0.0, δt, 0.0, false, false, Float64[λt1, λt1], 
                                                               Float64[μt1, μt1]),
                           bt, δt, t, false, false, λv, μv), nsp
          # if extinction
          else
            return iTgbmbd(bt, δt, t, true, false, λv, μv), nsp
          end
        end

        return iTgbmbd(bt, δt, t, false, false, λv, μv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        # if speciation
        if λorμ(λm, μm)
          nsp += 1
          td1, nsp = sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbmbd(t, λt1, μt1, σλ, σμ, δt, srδt, nsp, nlim)

          return iTgbmbd(td1, td2, bt, δt, δt, false, false, λv, μv), nsp
        # if extinction
        else
          return iTgbmbd(bt, δt, δt, true, false, λv, μv), nsp
        end
      end

      λt = λt1
      μt = μt1
    end
  else
    return iTgbmbd(), nsp
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



