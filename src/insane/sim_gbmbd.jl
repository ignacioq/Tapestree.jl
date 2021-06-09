#=

Anagenetic GBM birth-death Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    sim_gbm(t   ::Float64,
            λt  ::Float64,
            μt  ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Simulate `iTgbmbd` according to a geometric Brownian motion.
"""
function sim_gbm(t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

  λv = Float64[λt]
  μv = Float64[μt]
  bt = 0.0

  while true

    if t <= δt
      bt += t

      if t > 0.0
        srt = sqrt(t)
        λt1 = rnorm(λt, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, t)
          push!(λv, λm)
          push!(μv, μm)

          # if speciation
          if λorμ(λm, μm)
            return iTgbmbd(sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt), 
                           sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt), 
                    bt, δt, t, false, false, λv, μv)
          # if extinction
          else
            return iTgbmbd(nothing, nothing, bt, δt, 0.5 * t, true, false, λv, μv)
          end
        end

        push!(λv, λt1)
        push!(μv, μt1)
      end

      return iTgbmbd(nothing, nothing, bt, δt, t, false, false, λv, μv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt, srδt*σλ)
    μt1 = rnorm(μt, srδt*σμ)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, δt)
      push!(λv, λm)
      push!(μv, μm)

      # if speciation
      if λorμ(λm, μm)
        return iTgbmbd(sim_gbm(t + 0.5*δt, λm, μm, σλ, σμ, δt, srδt), 
                       sim_gbm(t + 0.5*δt, λm, μm, σλ, σμ, δt, srδt), 
                bt, δt, 0.5*δt, false, false, λv, μv)
      # if extinction
      else
        return iTgbmbd(nothing, nothing, bt, δt, 0.5*δt, true, false, λv, μv)
      end
    end

    push!(λv, λt1)
    push!(μv, μt1)

    λt = λt1
    μt = μt1
  end
end




"""
    sim_gbm(t   ::Float64,
            λt  ::Float64,
            μt  ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            δt  ::Float64,
            srδt::Float64,
            nsp ::Int64,
            nlim::Int64)

Simulate `iTgbmbd` according to a geometric Brownian motion with a limit
on the number lineages allowed to reach.
"""
function sim_gbm(t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
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

        if t > 0.0
          srt = sqrt(t)
          λt1 = rnorm(λt, srt*σλ)
          μt1 = rnorm(μt, srt*σμ)

          λm = exp(0.5*(λt + λt1))
          μm = exp(0.5*(μt + μt1))

          if divev(λm, μm, t)
            push!(λv, λm)
            push!(μv, μm)

            # if speciation
            if λorμ(λm, μm)
              nsp += 1
              td1, nsp = sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)
              td2, nsp = sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)

              return iTgbmbd(td1, td2, bt, δt, 0.5 * t, false, false, λv, μv), nsp
            else
            # if extinction
              return iTgbmbd(nothing, nothing, bt, δt, 0.5 * t, true, false, λv, μv), nsp
            end
          end

          push!(λv, λt1)
          push!(μv, μt1)
        end

        return iTgbmbd(nothing, nothing, bt, δt, t, false, false, λv, μv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        push!(λv, λm)
        push!(μv, μm)

        # if speciation
        if λorμ(λm, μm)
          nsp += 1
          td1, nsp = sim_gbm(t + 0.5*δt, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbm(t + 0.5*δt, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)

          return iTgbmbd(td1, td2, bt, δt, 0.5*δt, false, false, λv, μv), nsp
        # if extinction
        else
          return iTgbmbd(nothing, nothing, bt, δt, 0.5*δt, true, false, λv, μv), nsp
        end
      end
      push!(λv, λt1)
      push!(μv, μt1)

      λt = λt1
      μt = μt1
    end

  else
    return iTgbmbd(), nsp
  end
end





"""
    sim_gbm(nsδt::Float64,
            t   ::Float64,
            λt  ::Float64,
            μt  ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Simulate `iTgbmbd` according to a geometric Brownian motion starting 
with a non-standard δt.
"""
function sim_gbm(nsδt::Float64,
                 t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

  λv = Float64[λt]
  μv = Float64[μt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt += t

    if t > 0.0
      srt = sqrt(t)
      λt1 = rnorm(λt, srt*σλ)
      μt1 = rnorm(μt, srt*σμ)
      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, t)
        push!(λv, λm)
        push!(μv, μm)

        # if speciation
        if λorμ(λm, μm)
          return iTgbmbd(sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt), 
                         sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt), 
                  bt, δt, 0.5 * t, false, false, λv, μv)
        # if extinction
        else
          return iTgbmbd(nothing, nothing, bt, δt, 0.5 * t, true, false, λv, μv)
        end
      end
      push!(λv, λt1)
      push!(μv, μt1)
    end

    return iTgbmbd(nothing, nothing, bt, δt, t, false, false, λv, μv)
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt, srnsδt*σλ)
  μt1 = rnorm(μt, srnsδt*σμ)


  λm = exp(0.5*(λt + λt1))
  μm = exp(0.5*(μt + μt1))

  if divev(λm, μm, nsδt)
    push!(λv, λm)
    push!(μv, μm)

    # if speciation
    if λorμ(λm, μm)
      return iTgbmbd(sim_gbm(0.5 * nsδt, λm, μm, σλ, σμ, δt, srδt), 
                     sim_gbm(0.5 * nsδt, λm, μm, σλ, σμ, δt, srδt), 
              bt, 0.5 * nsδt, false, false, λv, μv)
    # if extinction
    else
      return iTgbmbd(nothing, nothing, bt, δt, 0.5 * nsδt, true, false, λv, μv)
    end
  end

  push!(λv, λt1)
  push!(μv, μt1)

  λt = λt1
  μt = μt1

  ## second: standard δt
  while true

    if t <= δt
      bt += t

      if t > 0.0
        srt = sqrt(t)
        λt1 = rnorm(λt, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, δt)
          push!(λv, λm)
          push!(μv, μm)

          # if speciation
          if λorμ(λm, μm)
            return iTgbmbd(sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt), 
                           sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt), 
                    bt, δt, 0.5 * t, false, false, λv, μv)
          # if extinction
          else
            return iTgbmbd(nothing, nothing, bt, δt, 0.5 * t, true, false, λv, μv)
          end
        end

        push!(λv, λt1)
        push!(μv, μt1)
      end

      return iTgbmbd(nothing, nothing, bt, δt, t, false, false, λv, μv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt, srδt*σλ)
    μt1 = rnorm(μt, srδt*σμ)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, δt)
      push!(λv, λm)
      push!(μv, μm)

      # if speciation
      if λorμ(λm, μm)
        return iTgbmbd(sim_gbm(t + 0.5*δt, λm, μm, σλ, σμ, δt, srδt), 
                       sim_gbm(t + 0.5*δt, λm, μm, σλ, σμ, δt, srδt), 
                bt, δt, 0.5*δt, false, false, λv, μv)
      # if extinction
      else
        return iTgbmbd(nothing, nothing, bt, δt, 0.5*δt, true, false, λv, μv)
      end
    end

    push!(λv, λt1)
    push!(μv, μt1)

    λt = λt1
    μt = μt1
  end
end





"""
    sim_gbm(nsδt::Float64,
            t   ::Float64,
            λt  ::Float64,
            μt  ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            δt  ::Float64,
            srδt::Float64, 
            nsp ::Int64,
            nlim::Int64)

Simulate `iTgbmbd` according to a geometric Brownian motion starting 
with a non-standard δt with a limit in the number of species.
"""
function sim_gbm(nsδt::Float64,
                 t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
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

    if t > 0.0
      srt = sqrt(t)
      λt1 = rnorm(λt, srt*σλ)
      μt1 = rnorm(μt, srt*σμ)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, t)
        push!(λv, λm)
        push!(μv, μm)

        # if speciation
        if λorμ(λm, μm)
          nsp += 1
          td1, nsp = sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)

          return iTgbmbd(td1, td2, bt, δt, 0.5 * t, false, false, λv, μv), nsp
        else
        # if extinction
          return iTgbmbd(nothing, nothing, bt, δt, 0.5 * t, true, false, λv, μv), nsp
        end
      end

      push!(λv, λt1)
      push!(μv, μt1)
    end

    return iTgbmbd(nothing, nothing, bt, δt, t, false, false, λv, μv), nsp
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt, srnsδt*σλ)
  μt1 = rnorm(μt, srnsδt*σμ)

  λm = exp(0.5*(λt + λt1))
  μm = exp(0.5*(μt + μt1))

  if divev(λm, μm, nsδt)
    push!(λv, λm)
    push!(μv, μm)

    # if speciation
    if λorμ(λm, μm)
      nsp += 1
      td1, nsp = sim_gbm(0.5 * nsδt, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)
      td2, nsp = sim_gbm(0.5 * nsδt, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)

      return iTgbmbd(td1, td2, bt, δt, 0.5 * nsδt, false, false, λv, μv), nsp
    else
    # if extinction
      return iTgbmbd(nothing, nothing, bt, δt, 0.5 * nsδt, true, false, λv, μv), nsp
    end
  end

  push!(λv, λt1)
  push!(μv, μt1)

  λt = λt1
  μt = μt1

  if nsp < nlim

    ## second: standard δt
    while true

      if t <= δt
        bt  += t
        if t > 0.0
          srt = sqrt(t)
          λt1 = rnorm(λt, srt*σλ)
          μt1 = rnorm(μt, srt*σμ)

          λm = exp(0.5*(λt + λt1))
          μm = exp(0.5*(μt + μt1))

          if divev(λm, μm, t)
            push!(λv, λm)
            push!(μv, μm)

            # if speciation
            if λorμ(λm, μm)
              nsp += 1
              td1, nsp = sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)
              td2, nsp = sim_gbm(0.5 * t, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)

              return iTgbmbd(td1, td2, bt, δt, 0.5 * t, false, false, λv, μv), nsp
            # if extinction
            else
              return iTgbmbd(nothing, nothing, bt, δt, 0.5 * t, true, false, λv, μv), nsp
            end
          end

          push!(λv, λt1)
          push!(μv, μt1)
        end

        return iTgbmbd(nothing, nothing, bt, δt, t, false, false, λv, μv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        push!(λv, λm)
        push!(μv, μm)

        # if speciation
        if λorμ(λm, μm)
          nsp += 1
          td1, nsp = sim_gbm(t + 0.5*δt, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)
          td2, nsp = sim_gbm(t + 0.5*δt, λm, μm, σλ, σμ, δt, srδt, nsp, nlim)

          return iTgbmbd(td1, td2, bt, δt, 0.5*δt, false, false, λv, μv), nsp
        # if extinction
        else
          return iTgbmbd(nothing, nothing, bt, δt, 0.5*δt, true, false, λv, μv), nsp
        end
      end

      push!(λv, λt1)
      push!(μv, μt1)

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




