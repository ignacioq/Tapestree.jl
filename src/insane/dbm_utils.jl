#=

Diffused Brownian motion utilities

Ignacio Quintero Mächler

t(-_-t)

Created 26 01 2024
=#




"""
    function dbm(xa  ::Float64,
                 αx  ::Float64,
                 lσ2a::Float64,
                 ασ  ::Float64,
                 γ   ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64,
                 n   ::Int64)

Returns a diffused Brownian motion vectors starting with rates 
`lσ2a` (`ln(σ²(t))`) and trait `xa`. 
"""
@inline function dbm(xa  ::Float64,
                     αx  ::Float64,
                     lσ2a::Float64,
                     ασ  ::Float64,
                     γ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64,
                     n   ::Int64)
  @inbounds begin
    l   = n + 2
    x   = randn(l)
    lσ2 = randn(l)

    # rates
    lσ2[1] = lσ2a
    if n > 0
      @turbo for i in Base.OneTo(n)
        lσ2[i+1] *= srδt*γ
        lσ2[i+1] += ασ*δt
      end
    end
    lσ2[l] *= sqrt(fdt)*γ
    lσ2[l] += ασ*fdt
    cumsum!(lσ2, lσ2)

    # values
    x[1] = xa
    if n > 0
      @turbo for i in Base.OneTo(n)
        x[i+1] *= srδt*exp(0.25*(lσ2[i] + lσ2[i+1]))
        x[i+1] += αx*δt
      end
    end
    x[l] *= sqrt(fdt)*exp(0.25*(lσ2[l-1] + lσ2[l]))
    x[l] += αx*fdt
    cumsum!(x, x)
  end

  return x, lσ2
end




"""
    dbm!(x   ::Vector{Float64},
         xa  ::Float64,
         αx  ::Float64,
         lσ2 ::Vector{Float64},
         lσ2a::Float64,
         ασ  ::Float64,
         γ   ::Float64,
         δt  ::Float64,
         fdt ::Float64,
         srδt::Float64)

Returns a diffused Brownian motion in place starting with rates `lσ2a` and 
trait `xa`.
"""
@inline function dbm!(x   ::Vector{Float64},
                      xa  ::Float64,
                      αx  ::Float64,
                      lσ2 ::Vector{Float64},
                      lσ2a::Float64,
                      ασ  ::Float64,
                      γ   ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(x)
    randn!(lσ2)

    # rates
    lσ2[1] = lσ2a
    if l > 2
      @turbo for i in Base.OneTo(l-2)
        lσ2[i+1] *= srδt*γ
        lσ2[i+1] += ασ*δt
      end
    end
    lσ2[l] *= sqrt(fdt)*γ
    lσ2[l] += ασ*fdt
    cumsum!(lσ2, lσ2)

    # values
    x[1] = xa
    if l > 2
      @turbo for i in Base.OneTo(l-2)
        x[i+1] *= srδt*exp(0.25*(lσ2[i] + lσ2[i+1]))
        x[i+1] += αx*δt
      end
    end
    x[l] *= sqrt(fdt)*exp(0.25*(lσ2[l-1] + lσ2[l]))
    x[l] += αx*fdt
    cumsum!(x, x)
  end

  return nothing
end




"""
    dbm!(x   ::Vector{Float64},
         xa  ::Float64,
         αx  ::Float64,
         lσ2 ::Vector{Float64},
         δt  ::Float64,
         fdt ::Float64,
         srδt::Float64)

Returns a diffused Brownian motion in place conditional on rate path `lσ2` and 
initial trait `xa`.
"""
@inline function dbm!(x   ::Vector{Float64},
                      xa  ::Float64,
                      αx  ::Float64,
                      lσ2 ::Vector{Float64},
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(x)
    # values
    x[1] = xa
    if l > 2
      @turbo for i in Base.OneTo(l-2)
        x[i+1] *= srδt*exp(0.25*(lσ2[i] + lσ2[i+1]))
        x[i+1] += αx*δt
      end
    end
    x[l] *= sqrt(fdt)*exp(0.25*(lσ2[l-1] + lσ2[l]))
    x[l] += αx*fdt
   cumsum!(x, x)
  end

  return nothing
end




"""
    dbm(xa  ::Float64,
        αx  ::Float64,
        lσ2 ::Vector{Float64},
        δt  ::Float64
        fdt ::Float64,
        srδt::Float64)

Returns a diffused Brownian motion conditional on rate path `lσ2` and 
initial trait `xa`.
"""
@inline function dbm(xa  ::Float64,
                     αx  ::Float64,
                     lσ2 ::Vector{Float64},
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64)

  @inbounds begin
    l = lastindex(lσ2)
    x = randn(l)
    # values
    x[1] = xa
    if l > 2
      @turbo for i in Base.OneTo(l-2)
        x[i+1] *= srδt*exp(0.25*(lσ2[i] + lσ2[i+1]))
        x[i+1] += αx*δt
      end
    end
    x[l] *= sqrt(fdt)*exp(0.25*(lσ2[l-1] + lσ2[l]))
    x[l] += αx*fdt
    cumsum!(x, x)
  end

  return x
end




"""
    dbb(xi  ::Float64,
        xf  ::Float64,
        lσ2i::Float64,
        lσ2f::Float64,
        γ   ::Float64,
        δt  ::Float64,
        fdt ::Float64,
        srδt::Float64,
        n   ::Int64)

Diffused Brownian bridge simulation.
"""
@inline function dbb(xi  ::Float64,
                     xf  ::Float64,
                     lσ2i::Float64,
                     lσ2f::Float64,
                     γ   ::Float64,
                     δt  ::Float64,
                     fdt ::Float64,
                     srδt::Float64,
                     n   ::Int64)

  @inbounds begin
    l      = n + 2
    lσ2    = randn(l)
    x      = randn(l)
    lσ2[1] = lσ2i
    x[1]   = xi
    if l > 2
      # rates
      for i in Base.OneTo(n)
        lσ2[i+1] *= srδt*γ
        lσ2[i+1] += lσ2[i]
      end
      lσ2[l] *= sqrt(fdt)*γ
      lσ2[l] += lσ2[l-1]

      # make rates bridge
      ite = 1.0/(Float64(n) * δt + fdt)
      lσ2df = (lσ2[l] - lσ2f)
      @turbo for i = Base.OneTo(l-1)
        lσ2[i] -= (Float64(i-1) * δt * ite * lσ2df)
      end

      # values
      @turbo for i = Base.OneTo(n)
        x[i+1] *= srδt*exp(0.25*(lσ2[i] + lσ2[i+1]))
      end
      x[l] *= sqrt(fdt)*exp(0.25*(lσ2[l-1] + lσ2[l]))
      cumsum!(x, x)

      # make bridge
      xdf = (x[l] - xf)
      @turbo for i = Base.OneTo(l-1)
        x[i] -= (Float64(i-1) * δt * ite * xdf)
      end
    end

    lσ2[l] = lσ2f
    x[l]   = xf
  end

  return x, lσ2
end




"""
    dbb!(x   ::Vector{Float64},
         xi  ::Float64,
         xf  ::Float64,
         lσ2 ::Vector{Float64},
         lσ2i::Float64,
         lσ2f::Float64,
         γ   ::Float64,
         δt  ::Float64,
         fdt ::Float64,
         srδt::Float64)

Diffused Brownian bridge simulation.
"""
@inline function dbb!(x   ::Vector{Float64},
                      xi  ::Float64,
                      xf  ::Float64,
                      lσ2 ::Vector{Float64},
                      lσ2i::Float64,
                      lσ2f::Float64,
                      γ   ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(lσ2)
    randn!(x)
    lσ2[1] = lσ2i
    x[1]   = xi

    if l > 2
      # rates
      for i in Base.OneTo(l-2)
        lσ2[i+1] *= srδt*γ
        lσ2[i+1] += lσ2[i]
      end
      lσ2[l] *= sqrt(fdt)*γ
      lσ2[l] += lσ2[l-1]

      # make rates bridge
      ite = 1.0/(Float64(l-2) * δt + fdt)
      lσ2df = (lσ2[l] - lσ2f)
      @turbo for i = Base.OneTo(l-1)
        lσ2[i] -= (Float64(i-1) * δt * ite * lσ2df)
      end

      # values
      @turbo for i = Base.OneTo(l-2)
        x[i+1] *= srδt*exp(0.25*(lσ2[i] + lσ2[i+1]))
      end
      x[l] *= sqrt(fdt)*exp(0.25*(lσ2[l-1] + lσ2[l]))
      cumsum!(x, x)

      # make values bridge
      xdf = (x[l] - xf)
      @turbo for i = Base.OneTo(l-1)
        x[i] -= (Float64(i-1) * δt * ite * xdf)
      end
    end

    lσ2[l] = lσ2f
    x[l]   = xf
  end
end



"""
    dbb!(x   ::Vector{Float64},
         xi  ::Float64,
         xf  ::Float64,
         lσ2 ::Vector{Float64},
         δt  ::Float64,
         fdt ::Float64,
         srδt::Float64)

Diffused Brownian bridge simulation conditional on a rate path `lσ2`.
"""
@inline function dbb!(x   ::Vector{Float64},
                      xi  ::Float64,
                      xf  ::Float64,
                      lσ2 ::Vector{Float64},
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)
    randn!(x)

    # values
    x[1] = xi
    if l > 2
      @turbo for i = Base.OneTo(l-2)
        x[i+1] *= srδt*exp(0.25*(lσ2[i] + lσ2[i+1]))
      end
    end
    x[l] *= sqrt(fdt)*exp(0.25*(lσ2[l-1] + lσ2[l]))
    cumsum!(x, x)

    # make values bridge
    ite = 1.0/(Float64(l-2) * δt + fdt)
    xdf = (x[l] - xf)
    if l > 2
      @turbo for i = Base.OneTo(l-1)
        x[i] -= (Float64(i-1) * δt * ite * xdf)
      end
    end
    x[l] = xf
  end
end




"""
    duoprop(x1 ::Float64,
            x2 ::Float64,
            σ21::Float64,
            σ22::Float64)

Proposal for a duo of Gaussians.
"""
function duoprop(x1 ::Float64,
                 x2 ::Float64,
                 σ21::Float64,
                 σ22::Float64)

  iσ2 = 1.0/(σ21 + σ22)
  return rnorm((x1*σ22 + x2*σ21) * iσ2, sqrt(σ21*σ22*iσ2))
end




"""
    trioprop(xp ::Float64,
             x1 ::Float64,
             x2 ::Float64,
             σ2p::Float64,
             σ21::Float64,
             σ22::Float64)

Proposal for a trio of Gaussians.
"""
function trioprop(xp ::Float64,
                  x1 ::Float64,
                  x2 ::Float64,
                  σ2p::Float64,
                  σ21::Float64,
                  σ22::Float64)

  s = 1.0/(1.0/σ2p + 1.0/σ21 + 1.0/σ22)
  return rnorm((xp/σ2p + x1/σ21 + x2/σ22)*s, sqrt(s))
end




"""
    intσ2(lσ2::Vector{Float64},
          δt ::Float64,
          fdt::Float64)

Returns the integral for a vector of log-rates `lσ2`.
"""
function intσ2(lσ2::Vector{Float64},
               δt ::Float64,
               fdt::Float64)

  l = lastindex(lσ2)

  ss = 0.0
  if l > 2
    @turbo for i = Base.OneTo(l-2)
      ss += exp(0.5*(lσ2[i] + lσ2[i+1]))
    end
    ss *= δt
    ss += exp(0.5*(lσ2[l-1] + lσ2[l])) * fdt
  end

  return ss
end


