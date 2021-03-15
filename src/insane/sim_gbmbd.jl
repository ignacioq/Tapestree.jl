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
  ts = Float64[bt]

  while true

    tc = t - δt

    if tc < 0.0
      if t > 1e-11
        bt  += t
        srt = sqrt(t)
        λt1 = rnorm(λt, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)
        push!(ts, bt)
      end

      return iTgbmbd(nothing, nothing, bt, false, false, ts, λv, μv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt, srδt*σλ)
    μt1 = rnorm(μt, srδt*σμ)

    push!(λv, λt1)
    push!(μv, μt1)
    push!(ts, bt)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, δt)
      # if speciation
      if λorμ(λm, μm)
        return iTgbmbd(sim_gbm(t, λt1, μt1, σλ, σμ, δt, srδt), 
                       sim_gbm(t, λt1, μt1, σλ, σμ, δt, srδt), 
                bt, false, false, ts, λv, μv)
      # if extinction
      else
        return iTgbmbd(nothing, nothing, bt, true, false, ts, λv, μv)
      end
    end

    λt = λt1
    μt = μt1
  end
end





"""
    sim_ov_gbm(t   ::Float64,
               ii  ::Int64,
               tl  ::Int64,
               bbiλ::Array{Float64,1},
               bbiμ::Array{Float64,1},
               tsi ::Array{Float64,1},
               σλ  ::Float64,
               σμ  ::Float64,
               δt  ::Float64,
               srδt::Float64)

Simulate `iTgbmbd` on top of a Brownian bridge for a branch.
"""
function sim_ov_gbm(t   ::Float64,
                    ii  ::Int64,
                    tl  ::Int64,
                    bbiλ::Array{Float64,1},
                    bbiμ::Array{Float64,1},
                    tsi ::Array{Float64,1},
                    σλ  ::Float64,
                    σμ  ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  λt = bbiλ[ii]
  μt = bbiμ[ii]
  λv = Float64[λt]
  μv = Float64[μt]
  bt = 0.0
  tv = Float64[bt]

  for i in ii:(tl-2)

    t  -= δt
    bt += δt

    λt1 = bbiλ[i+1]
    μt1 = bbiμ[i+1]

    push!(λv, λt1)
    push!(μv, μt1)
    push!(tv, bt)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, δt)
      # if speciation
      if λorμ(λm, μm)
        if rand() < 0.5
          return iTgbmbd(
                   sim_ov_gbm(t, i+1, tl, bbiλ, bbiμ, tsi, σλ, σμ, δt, srδt), 
                   sim_gbm(t, λt1, μt1, σλ, σμ, δt, srδt), 
                 bt, false, true, tv, λv, μv)
        else
          return iTgbmbd(
                   sim_gbm(t, λt1, μt1, σλ, σμ, δt, srδt), 
                   sim_ov_gbm(t, i+1, tl, bbiλ, bbiμ, tsi, σλ, σμ, δt, srδt), 
                 bt, false, true, tv, λv, μv)
        end
      # if extinction
      else
        return iTgbmbd(nothing, nothing, bt, true, true, tv, λv, μv)
      end
    end

    λt = λt1
    μt = μt1
  end

  bt += (tsi[tl] - tsi[tl-1])
  λt1 = bbiλ[tl]
  μt1 = bbiμ[tl]

  push!(λv, λt1)
  push!(μv, μt1)
  push!(tv, bt)

  return iTgbmbd(nothing, nothing, bt, false, true, tv, λv, μv)
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
  ts = Float64[bt]

  ## first: non-standard δt
  tc = t - nsδt

  if tc < 0.0
    if t > 1e-11
      bt  += t
      srt = sqrt(t)
      λt1 = rnorm(λt, srt*σλ)
      μt1 = rnorm(μt, srt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)
      push!(ts, bt)
    end

    return iTgbmbd(nothing, nothing, bt, false, false, ts, λv, μv)
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt, srnsδt*σλ)
  μt1 = rnorm(μt, srnsδt*σμ)

  push!(λv, λt1)
  push!(μv, μt1)
  push!(ts, bt)

  λm = exp(0.5*(λt + λt1))
  μm = exp(0.5*(μt + μt1))

  if divev(λm, μm, nsδt)
    # if speciation
    if λorμ(λm, μm)
      return iTgbmbd(sim_gbm(t, λt1, μt1, σλ, σμ, δt, srδt), 
                     sim_gbm(t, λt1, μt1, σλ, σμ, δt, srδt), 
              bt, false, false, ts, λv, μv)
    # if extinction
    else
      return iTgbmbd(nothing, nothing, bt, true, false, ts, λv, μv)
    end
  end

  λt = λt1
  μt = μt1

  ## second: standard δt
  while true

    tc = t - δt

    if tc < 0.0
      if t > 1e-11
        bt  += t
        srt = sqrt(t)
        λt1 = rnorm(λt, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)
        push!(ts, bt)
      end

      return iTgbmbd(nothing, nothing, bt, false, false, ts, λv, μv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt, srδt*σλ)
    μt1 = rnorm(μt, srδt*σμ)

    push!(λv, λt1)
    push!(μv, μt1)
    push!(ts, bt)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, δt)
      # if speciation
      if λorμ(λm, μm)
        return iTgbmbd(sim_gbm(t, λt1, μt1, σλ, σμ, δt, srδt), 
                       sim_gbm(t, λt1, μt1, σλ, σμ, δt, srδt), 
                bt, false, false, ts, λv, μv)
      # if extinction
      else
        return iTgbmbd(nothing, nothing, bt, true, false, ts, λv, μv)
      end
    end

    λt = λt1
    μt = μt1
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




