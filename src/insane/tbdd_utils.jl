#=

Correlated diffusion utilities

Ignacio Quintero Mächler

t(-_-t)

Created 12 02 2026
=#



"""
    cbb!(lλv ::Vector{Float64},
         lλi ::Float64,
         lλf ::Float64,
         x   ::Vector{Float64},
         βλ  ::Float64,
         σλ  ::Float64,
         δt  ::Float64,
         fdt ::Float64,
         srδt::Float64)

Correlated Brownian bridge simulation conditional on `x`.
"""
@inline function cbb!(lλv ::Vector{Float64},
                      lλi ::Float64,
                      lλf ::Float64,
                      x   ::Vector{Float64},
                      βλ  ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)
  @inbounds begin
    l = lastindex(x)
    randn!(lλv)
    lλv[1] = lλi

    if l > 2
      # speciation rates conditional con x
      srδtσλ = srδt*σλ
      for i in Base.OneTo(l-2)
        lλv[i+1] *= srδtσλ
        lλv[i+1] += lλv[i] + βλ*(x[i+1] - x[i])
      end
      lλv[l] *= sqrt(fdt)*σλ
      lλv[l] += lλv[l-1] + βλ*(x[l] - x[l-1])

      # make rates bridge
      ite = (lλv[l] - lλf) * δt/(Float64(l-2) * δt + fdt)
      @turbo for i = Base.OneTo(l-1)
        lλv[i] -= Float64(i-1) * ite
      end
    end

    lλv[l] = lλf
  end
end




"""
    cbb!(x   ::Vector{Float64},
        xi  ::Float64,
        xf  ::Float64,
        lσ2 ::Vector{Float64},
        lλv ::Vector{Float64},
        lλi ::Float64,
        lλf ::Float64,
        βλ  ::Float64,
        σλ  ::Float64,
        δt  ::Float64,
        fdt ::Float64,
        srδt::Float64)

Correlated Brownian bridge simulation for `x` and `lλ` conditional 
on trait evolutionary rates `lσ2`.
"""
@inline function cbb!(x   ::Vector{Float64},
                      xi  ::Float64,
                      xf  ::Float64,
                      lσ2 ::Vector{Float64},
                      lλv ::Vector{Float64},
                      lλi ::Float64,
                      lλf ::Float64,
                      βλ  ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64)

  @inbounds begin
    l = lastindex(x)

    # generate random normals
    randn!(x)
    randn!(lλv)

    # initial values
    x[1]   = xi
    lλv[1] = lλi
    if l > 2
      srδtσλ = srδt*σλ
      for i = Base.OneTo(l-2)
        x[i+1]   *= srδt*exp(0.25*(lσ2[i] + lσ2[i+1]))
        xi        = x[i]
        x[i+1]   += xi

        lλv[i+1] *= srδtσλ
        lλv[i+1] += lλv[i] + βλ*(x[i+1] - xi)
      end
      srfdt   = sqrt(fdt)
      x[l]   *= srfdt*exp(0.25*(lσ2[l-1] + lσ2[l]))
      xlm1    = x[l-1]
      x[l]   += xlm1

      lλv[l] *= srfdt*σλ
      lλv[l] += lλv[l-1] + βλ*(x[l] - xlm1)

      # make values bridge
      ite  = δt/(Float64(l-2) * δt + fdt)
      itex = (x[l]   - xf)  * ite
      iteλ = (lλv[l] - lλf) * ite
      if l > 2
        @turbo for i = Base.OneTo(l-1)
          x[i]   -= Float64(i-1) * itex
          lλv[i] -= Float64(i-1) * iteλ
        end
      end
    end
    x[l]   = xf
    lλv[l] = lλf
  end
end


