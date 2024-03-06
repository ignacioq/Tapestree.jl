#=

Diffused Brownian motion utilities

Ignacio Quintero Mächler

t(-_-t)

Created 26 01 2024
=#




"""
    function dbm(xa  ::Float64,
                 σa  ::Float64,
                 α   ::Float64,
                 γ   ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64,
                 n   ::Int64)

Returns a diffused Brownian motion vectors starting with rates `σa` and 
trait `xa`.
"""
@inline function dbm(xa  ::Float64,
                     σa  ::Float64,
                     α   ::Float64,
                     γ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64,
                     n   ::Int64)
  @inbounds begin
    l = n + 2
    x = randn(l)
    σ = randn(l)

    # rates
    σ[1] = σa
    @turbo for i in Base.OneTo(n)
      σ[i+1] *= srδt*γ
      σ[i+1] += α*δt
    end
    σ[l] *= sqrt(fdt)*γ
    σ[l] += α*fdt
    cumsum!(σ, σ)

    # values
    x[1] = xa
    @turbo for i in Base.OneTo(n)
      x[i+1] *= srδt*exp(0.5*(σ[i] + σ[i+1]))
    end
    x[l] *= sqrt(fdt)*exp(0.5*(σ[l-1] + σ[l]))
    cumsum!(x, x)
  end

  return x, σ
end




"""
    dbm!(x   ::Vector{Float64},
         xa  ::Float64,
         σ   ::Vector{Float64},
         σa  ::Float64,
         α   ::Float64,
         γ   ::Float64,
         δt  ::Float64
         fdt ::Float64,
         srδt::Float64)

Returns a diffused Brownian motion in place starting with rates `σa` and 
trait `xa`.
"""
@inline function dbm!(x   ::Vector{Float64},
                      xa  ::Float64,
                      σ   ::Vector{Float64},
                      σa  ::Float64,
                      α   ::Float64,
                      γ   ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(x)
    randn!(σ)

    # rates
    σ[1] = σa
    @turbo for i in Base.OneTo(l-2)
      σ[i+1] *= srδt*γ
      σ[i+1] += α*δt
    end
    σ[l] *= sqrt(fdt)*γ
    σ[l] += α*fdt
    cumsum!(σ, σ)

    # values
    x[1] = xa
    @turbo for i in Base.OneTo(l-2)
      x[i+1] *= srδt*exp(0.5*(σ[i] + σ[i+1]))
    end
    x[l] *= sqrt(fdt)*exp(0.5*(σ[l-1] + σ[l]))
    cumsum!(x, x)
  end

  return nothing
end




"""
    dbm!(x   ::Vector{Float64},
         xa  ::Float64,
         σ   ::Vector{Float64},
         fdt ::Float64,
         srδt::Float64)

Returns a diffused Brownian motion in place conditional on rate path `σ` and 
initial trait `xa`.
"""
@inline function dbm!(x   ::Vector{Float64},
                      xa  ::Float64,
                      σ   ::Vector{Float64},
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(x)
    # values
    x[1] = xa
    @turbo for i in Base.OneTo(l-2)
      x[i+1] *= srδt*exp(0.5*(σ[i] + σ[i+1]))
    end
    x[l] *= sqrt(fdt)*exp(0.5*(σ[l-1] + σ[l]))
    cumsum!(x, x)
  end

  return nothing
end




"""
    dbm(xa  ::Float64,
        σ   ::Vector{Float64},
        fdt ::Float64,
        srδt::Float64)

Returns a diffused Brownian motion conditional on rate path `σ` and 
initial trait `xa`.
"""
@inline function dbm(xa  ::Float64,
                     σ   ::Vector{Float64},
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin
    l = lastindex(σ)
    x = randn(l)
    # values
    x[1] = xa
    @turbo for i in Base.OneTo(l-2)
      x[i+1] *= srδt*exp(0.5*(σ[i] + σ[i+1]))
    end
    x[l] *= sqrt(fdt)*exp(0.5*(σ[l-1] + σ[l]))
    cumsum!(x, x)
  end

  return x
end




"""
    dbb(xi  ::Float64,
        xf  ::Float64,
        σi  ::Float64,
        σf  ::Float64,
        γ   ::Float64,
        δt  ::Float64,
        fdt ::Float64,
        srδt::Float64,
        n   ::Int64)

Diffused Brownian bridge simulation.
"""
@inline function dbb(xi  ::Float64,
                     xf  ::Float64,
                     σi  ::Float64,
                     σf  ::Float64,
                     γ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64,
                     n   ::Int64)

  @inbounds begin
    l = n + 2
    σ = randn(l)
    x = randn(l)
    σ[1] = σi
    x[1] = xi
    if l > 2
      # rates
      for i in Base.OneTo(n)
        σ[i+1] *= srδt*γ
        σ[i+1] += σ[i]
      end
      σ[l] *= sqrt(fdt)*γ
      σ[l] += σ[l-1]

      # make rates bridge
      ite = 1.0/(Float64(n) * δt + fdt)
      σdf = (σ[l] - σf)
      @turbo for i = Base.OneTo(l-1)
        σ[i] -= (Float64(i-1) * δt * ite * σdf)
      end

      # values
      @turbo for i = Base.OneTo(n)
        x[i+1] *= srδt*exp(0.5*(σ[i] + σ[i+1]))
      end
      x[l] *= sqrt(fdt)*exp(0.5*(σ[l-1] + σ[l]))
      cumsum!(x, x)

      # make bridge
      xdf = (x[l] - xf)
      @turbo for i = Base.OneTo(l-1)
        x[i] -= (Float64(i-1) * δt * ite * xdf)
      end
    end

    σ[l] = σf
    x[l] = xf
  end

  return x, σ
end




"""
    dbb!(x   ::Vector{Float64},
         xi  ::Float64,
         xf  ::Float64,
         σ   ::Vector{Float64},
         σi  ::Float64,
         σf  ::Float64,
         γ   ::Float64,
         δt  ::Float64,
         fdt ::Float64,
         srδt::Float64)

Diffused Brownian bridge simulation.
"""
@inline function dbb!(x   ::Vector{Float64},
                      xi  ::Float64,
                      xf  ::Float64,
                      σ   ::Vector{Float64},
                      σi  ::Float64,
                      σf  ::Float64,
                      γ   ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(σ)
    randn!(x)
    σ[1] = σi
    x[1] = xi

    if l > 2
      # rates
      for i in Base.OneTo(l-2)
        σ[i+1] *= srδt*γ
        σ[i+1] += σ[i]
      end
      σ[l] *= sqrt(fdt)*γ
      σ[l] += σ[l-1]

      # make rates bridge
      ite = 1.0/(Float64(l-2) * δt + fdt)
      σdf = (σ[l] - σf)
      @turbo for i = Base.OneTo(l-1)
        σ[i] -= (Float64(i-1) * δt * ite * σdf)
      end

      # values
      @turbo for i = Base.OneTo(l-2)
        x[i+1] *= srδt*exp(0.5*(σ[i] + σ[i+1]))
      end
      x[l] *= sqrt(fdt)*exp(0.5*(σ[l-1] + σ[l]))
      cumsum!(x, x)

      # make values bridge
      xdf = (x[l] - xf)
      @turbo for i = Base.OneTo(l-1)
        x[i] -= (Float64(i-1) * δt * ite * xdf)
      end
    end

    σ[l] = σf
    x[l] = xf
  end
end



"""
    dbb!(x   ::Vector{Float64},
         xi  ::Float64,
         xf  ::Float64,
         σ   ::Vector{Float64},
         δt  ::Float64,
         fdt ::Float64,
         srδt::Float64)

Diffused Brownian bridge simulation conditional on a rate path `σ`.
"""
@inline function dbb!(x   ::Vector{Float64},
                      xi  ::Float64,
                      xf  ::Float64,
                      σ   ::Vector{Float64},
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(x)

    # values
    x[1] = xi
    @turbo for i = Base.OneTo(l-2)
      x[i+1] *= srδt*exp(0.5*(σ[i] + σ[i+1]))
    end
    x[l] *= sqrt(fdt)*exp(0.5*(σ[l-1] + σ[l]))
    cumsum!(x, x)

    # make values bridge
    ite = 1.0/(Float64(l-2) * δt + fdt)
    xdf = (x[l] - xf)
    @turbo for i = Base.OneTo(l-1)
      x[i] -= (Float64(i-1) * δt * ite * xdf)
    end
    x[l] = xf
  end
end





"""
    duoprop(x1::Float64,
            x2::Float64,
            σ1::Float64,
            σ2::Float64)

Proposal for a duo of Gaussians.
"""
function duoprop(x1::Float64,
                 x2::Float64,
                 σ1::Float64,
                 σ2::Float64)

  iσ = 1.0/(σ1^2 + σ2^2)
  return rnorm((x1*σ2^2 + x2*σ1^2) * iσ, 
                σ1*σ2*sqrt(iσ))
end




"""
    trioprop(xpr::Float64,
             xd1::Float64,
             xd2::Float64,
             σpr::Float64,
             σd1::Float64,
             σd2::Float64)

Proposal for a trio of Gaussians.
"""
function trioprop(xpr::Float64,
                  xd1::Float64,
                  xd2::Float64,
                  σpr::Float64,
                  σd1::Float64,
                  σd2::Float64)

  s = 1.0/(1.0/(σpr^2) + 1.0/(σd1^2) + 1.0/(σd2^2))
  return rnorm((xpr/(σpr^2) + xd1/(σd1^2) + xd2/(σd2^2))*s,
                sqrt(s))
end




"""
    intσ(σ  ::Vector{Float64},
         δt ::Float64,
         fdt::Float64)

Returns the integral for a vector of log-rates `σ`.
"""
function intσ(σ  ::Vector{Float64},
              δt ::Float64,
              fdt::Float64)

  l = lastindex(σ)

  ss = 0.0
  @turbo for i = Base.OneTo(l-2)
    ss += exp(σ[i] + σ[i+1])
  end
  ss *= δt
  ss += exp(σ[l-1] + σ[l]) * fdt

  return sqrt(ss)
end




