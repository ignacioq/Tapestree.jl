#=

trait pure-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_xb(Ξ  ::Vector{iTxb},
            idf::Vector{iBffs},
            ασ ::Float64,
            σσ ::Float64,
            αλ ::Float64,
            βλ ::Float64,
            σλ ::Float64,
            δt ::Float64)

Returns the log-likelihood for a `iTxb` according to trait pure-birth diffusion.
"""
function llik_xb(Ξ  ::Vector{iTxb},
                 idf::Vector{iBffs},
                 ασ ::Float64,
                 σσ ::Float64,
                 αλ ::Float64,
                 βλ ::Float64,
                 σλ ::Float64,
                 δt ::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      ll += llik_xb(Ξ[i], ασ, σσ, αλ, βλ, σλ, δt)
      if d2(idf[i]) > 0
        ll += λt(Ξ[i])
      end
    end
  end

  return ll
end




"""
    llik_xb(tree::iTxb,
            ασ  ::Float64,
            σσ  ::Float64,
            αλ  ::Float64,
            βλ  ::Float64,
            σλ  ::Float64,
            δt  ::Float64)

Returns the log-likelihood for a `iTxb` according to trait pure-birth diffusion.
"""
function llik_xb(tree::iTxb,
                 ασ  ::Float64,
                 σσ  ::Float64,
                 αλ  ::Float64,
                 βλ  ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64)

  if istip(tree)
    ll_xb_b(xv(tree), lσ2(tree), lλ(tree), 
            ασ, σσ, αλ, βλ, σλ, δt, fdt(tree), false)
  else
    ll_xb_b(xv(tree), lσ2(tree), lλ(tree), 
            ασ, σσ, αλ, βλ, σλ, δt, fdt(tree), true) +
    llik_xb(tree.d1::iTxb, ασ, σσ, αλ, βλ, σλ, δt)   +
    llik_xb(tree.d2::iTxb, ασ, σσ, αλ, βλ, σλ, δt)
  end
end




"""
    ll_xb_b(vx  ::Array{Float64,1},
            vlσ2::Array{Float64,1},
            vlλ ::Array{Float64,1},
            ασ  ::Float64,
            σσ  ::Float64,
            αλ  ::Float64,
            βλ  ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            fdt ::Float64,
            λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_xb_b(vx  ::Array{Float64,1},
                 vlσ2::Array{Float64,1},
                 vlλ ::Array{Float64,1},
                 ασ  ::Float64,
                 σσ  ::Float64,
                 αλ  ::Float64,
                 βλ  ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 λev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(vlλ)-2
    n  = Float64(nI)

    ll = llx = llσ2 = llλ = llpb = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        lσ2i   = vlσ2[i]
        lσ2i1  = vlσ2[i+1]
        dxi    = vx[i+1] - vx[i]
        lλi    = vlλ[i]
        lλi1   = vlλ[i+1]
        llx   += -0.5*dxi^2/(exp(0.5*(lσ2i1 + lσ2i))*δt) - 
                  0.25*(lσ2i1 + lσ2i)
        llσ2  += (lσ2i1 - lσ2i - ασ*δt)^2
        llλ   += (lλi1 - lλi - αλ*δt - βλ*dxi)^2
        llpb  += exp(0.5*(lλi + lλi1))
      end

      # add to global likelihood
      ll += llx - n*(1.5*log(δt) + log(σσ*σλ))                      +
            llσ2*(-0.5/(σσ^2*δt)) + llλ*(-0.5/(σλ^2*δt))            -
            n*2.756815599614018008622906563687138259410858154296875 - # 1.5 * log(2.0π)
            llpb*δt
    end

    lλi1 = vlλ[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lσ2i   = vlσ2[nI+1]
      lσ2i1  = vlσ2[nI+2]
      dxi    = vx[nI+2] - vx[nI+1]
      lλi    = vlλ[nI+1]
      # add to global likelihood
      ll += -0.5*dxi^2/(exp(0.5*(lσ2i1 + lσ2i))*fdt)              - 
            0.25*(lσ2i1 + lσ2i) - 1.5*log(fdt) - log(σσ*σλ)       +
            (lσ2i1 - lσ2i - ασ*fdt)^2 * (-0.5/(σσ^2*fdt))         + 
            (lλi1 - lλi - αλ*fdt - βλ*dxi)^2 * (-0.5/(σλ^2*fdt))  -
            2.756815599614018008622906563687138259410858154296875 - # 1.5 * log(2.0π)
            exp(0.5*(lλi + lλi1))*fdt
    end
    if λev
      ll += lλi1
    end
  end

  return ll
end




"""
    llr_xb_σ(vx  ::Array{Float64,1},
             ασ  ::Float64,
             σσ  ::Float64,
             lσ2p::Array{Float64,1},
             lσ2c::Array{Float64,1},
             δt  ::Float64,
             fdt ::Float64)

Returns the acceptance ratio and changes in gibbs quanta for a `σ²(t)` 
path proposal (the likelihood for the GBM for `σ²` cancels out).
"""
function llr_xb_σ(vx  ::Array{Float64,1},
                  ασ  ::Float64,
                  σσ  ::Float64,
                  lσ2p::Array{Float64,1},
                  lσ2c::Array{Float64,1},
                  δt  ::Float64,
                  fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(vx)-2

    llσxr = llσσr = ssσr = 0.0
    if nI > 0
      dσxr = 0.0
      @turbo for i in Base.OneTo(nI)
        lσ2ci  = lσ2c[i]
        lσ2ci1 = lσ2c[i+1]
        dlσ2c  = lσ2ci1 - lσ2ci
        lσ2pi  = lσ2p[i]
        lσ2pi1 = lσ2p[i+1]
        dlσ2p  = lσ2pi1 - lσ2pi
        llσσr += (dlσ2p - ασ*δt)^2 - (dlσ2c - ασ*δt)^2
        llσxr += -0.5*(vx[i+1] - vx[i])^2/δt *
                 (1.0/exp(0.5*(lσ2pi + lσ2pi1)) - 
                  1.0/exp(0.5*(lσ2ci + lσ2ci1))) + 
                 0.25*(lσ2ci + lσ2ci1 - lσ2pi - lσ2pi1)
        dσxr  += dlσ2c   - dlσ2p
        ssσr  += dlσ2p^2 - dlσ2c^2
      end
      llσσr *= -0.5/(σσ^2*δt)
      ssσr  /= 2.0 * δt
      ssσr  += (ασ * dσxr)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      lσ2ci  = lσ2c[nI+1]
      lσ2ci1 = lσ2c[nI+2]
      dlσ2c  = lσ2ci1 - lσ2ci
      lσ2pi  = lσ2p[nI+1]
      lσ2pi1 = lσ2p[nI+2]
      dlσ2p  = lσ2pi1 - lσ2pi
      llσσr += -0.5*((dlσ2p - ασ*fdt)^2 - (dlσ2c - ασ*fdt)^2)/(σσ^2*fdt)
      llσxr += -(0.5*(vx[nI+2] - vx[nI+1])^2/fdt) *
                (1.0/exp(0.5*(lσ2pi + lσ2pi1)) - 
                 1.0/exp(0.5*(lσ2ci + lσ2ci1))) + 
                 0.25*(lσ2ci + lσ2ci1 - lσ2pi - lσ2pi1)
      ssσr  += (dlσ2p^2 - dlσ2c^2)/(2.0*fdt) + ασ*(dlσ2c - dlσ2p)
    end
  end

  return llσxr, llσσr, ssσr
end




"""
    llr_xb_b_sep(vxp ::Array{Float64,1},
                 vxc ::Array{Float64,1},
                 vlσ2::Array{Float64,1},
                 lλp ::Array{Float64,1},
                 lλc ::Array{Float64,1},
                 ασ  ::Float64,
                 σσ  ::Float64,
                 αλ  ::Float64,
                 βλ  ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64,
                 λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth
separately for the Brownian motion and the pure-birth
"""
function llr_xb_b_sep(vxp ::Array{Float64,1},
                      vxc ::Array{Float64,1},
                      vlσ2::Array{Float64,1},
                      lλp ::Array{Float64,1},
                      lλc ::Array{Float64,1},
                      ασ  ::Float64,
                      σσ  ::Float64,
                      αλ  ::Float64,
                      βλ  ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64,
                      λev ::Bool)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλp)-2

    llbmr = llbr = dxsr = dxlr = ssλr = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        dxpi   = vxp[i+1] - vxp[i]
        dxci   = vxc[i+1] - vxc[i]
        lλpi   = lλp[i]
        lλpi1  = lλp[i+1]
        dlλpi  = lλpi1 - lλpi
        lλci   = lλc[i]
        lλci1  = lλc[i+1]
        dlλci  = lλci1 - lλci
        llbmr += 0.5 * (dxci^2 - dxpi^2)/(exp(0.5*(vlσ2[i] + vlσ2[i+1]))*δt) 
        llbr  += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
        dxsr  += dxpi^2 - dxci^2
        dxlr  += dxpi * dlλpi - dxci * dlλci
        ssλr  += (dlλpi - αλ*δt - βλ*dxpi)^2 - (dlλci - αλ*δt - βλ*dxci)^2
      end

      llbmr += ssλr*(-0.5/(σλ^2*δt)) 
      llbr  *= -δt
      dxsr  /= δt
      dxlr  /= δt
      ssλr  /= 2.0*δt
    end

    lλpi1 = lλp[nI+2]
    lλci1 = lλc[nI+2]

   # add final non-standard `δt`
    if fdt > 0.0
      dxpi   = vxp[nI+2] - vxp[nI+1]
      dxci   = vxc[nI+2] - vxc[nI+1]
      lλpi   = lλp[nI+1]
      dlλpi  = lλpi1 - lλpi
      lλci   = lλc[nI+1]
      dlλci  = lλci1 - lλci
      ssλ0r  = (dlλpi - αλ*fdt - βλ*dxpi)^2 - (dlλci - αλ*fdt - βλ*dxci)^2
      llbmr += 0.5*(dxci^2 - dxpi^2)/(exp(0.5*(vlσ2[nI+1] + vlσ2[nI+2]))*fdt) +
               ssλ0r*(-0.5/(σλ^2*fdt))
      llbr  += -fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
      dxsr  += (dxpi^2 - dxci^2)/fdt
      dxlr  += (dxpi * dlλpi - dxci * dlλci)/fdt
      ssλr  += ssλ0r/(2.0*fdt)
    end

    irλr = -llbr 

    #if speciation
    if λev
      llbr  += lλpi1 - lλci1
    end
  end

  return llbmr, llbr, dxsr, dxlr, ssλr, irλr
end





"""
    llr_xb_b_sep(vxp  ::Array{Float64,1},
                 vxc  ::Array{Float64,1},
                 vlσ2p::Array{Float64,1},
                 vlσ2c::Array{Float64,1},
                 lλp  ::Array{Float64,1},
                 lλc  ::Array{Float64,1},
                 ασ   ::Float64,
                 σσ   ::Float64,
                 αλ   ::Float64,
                 βλ   ::Float64,
                 σλ   ::Float64,
                 δt   ::Float64,
                 fdt  ::Float64,
                 srδt ::Float64,
                 λev  ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth
separately for the Brownian motion and the pure-birth
"""
function llr_xb_b_sep(vxp  ::Array{Float64,1},
                      vxc  ::Array{Float64,1},
                      lσ2p::Array{Float64,1},
                      lσ2c::Array{Float64,1},
                      lλp  ::Array{Float64,1},
                      lλc  ::Array{Float64,1},
                      ασ   ::Float64,
                      σσ   ::Float64,
                      αλ   ::Float64,
                      βλ   ::Float64,
                      σλ   ::Float64,
                      δt   ::Float64,
                      fdt  ::Float64,
                      srδt ::Float64,
                      λev  ::Bool)

  """
  here: test this
  """

  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(lλp)-2

    llbmr = llbr = dxsr = dxlr = ssσr = ssλr = 0.0
    if nI > 0
      dσxr = 0.0
      @turbo for i in Base.OneTo(nI)
        dxpi   = vxp[i+1] - vxp[i]
        dxci   = vxc[i+1] - vxc[i]
        lσ2ci  = lσ2c[i]
        lσ2ci1 = lσ2c[i+1]
        dlσ2c  = lσ2ci1 - lσ2ci
        lσ2pi  = lσ2p[i]
        lσ2pi1 = lσ2p[i+1]
        dlσ2p  = lσ2pi1 - lσ2pi
        lλpi   = lλp[i]
        lλpi1  = lλp[i+1]
        dlλpi  = lλpi1 - lλpi
        lλci   = lλc[i]
        lλci1  = lλc[i+1]
        dlλci  = lλci1 - lλci
        llbmr += 0.5*dxci^2/(exp(0.5*(lσ2ci + lσ2ci1))*δt) -
                 0.5*dxpi^2/(exp(0.5*(lσ2pi + lσ2pi1))*δt) +
                 0.25*(lσ2ci + lσ2ci1 - lσ2pi - lσ2pi1)
        llσσr += (dlσ2p - ασ*δt)^2 - (dlσ2c - ασ*δt)^2
        llbr  += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
        dxsr  += dxpi^2 - dxci^2
        dxlr  += dxpi * dlλpi - dxci * dlλci
        dσxr  += dlσ2c   - dlσ2p
        ssσr  += dlσ2p^2 - dlσ2c^2
        ssλr  += (dlλpi - αλ*δt - βλ*dxpi)^2 - (dlλci - αλ*δt - βλ*dxci)^2
      end
      llbmr += ssλr*(-0.5/(σλ^2*δt)) + llσσr*(-0.5/(σσ^2*δt))
      llbr  *= -δt
      dxsr  /= δt
      dxlr  /= δt
      ssσr  /= 2.0 * δt
      ssσr  += (ασ * dσxr)
      ssλr  /= 2.0*δt
    end

    lλpi1 = lλp[nI+2]
    lλci1 = lλc[nI+2]

   # add final non-standard `δt`
    if fdt > 0.0
      dxpi   = vxp[nI+2] - vxp[nI+1]
      dxci   = vxc[nI+2] - vxc[nI+1]
      lσ2ci  = lσ2c[nI+1]
      lσ2ci1 = lσ2c[nI+2]
      dlσ2c  = lσ2ci1 - lσ2ci
      lσ2pi  = lσ2p[nI+1]
      lσ2pi1 = lσ2p[nI+2]
      dlσ2p  = lσ2pi1 - lσ2pi
      lλpi   = lλp[nI+1]
      dlλpi  = lλpi1 - lλpi
      lλci   = lλc[nI+1]
      dlλci  = lλci1 - lλci
      ssλ0r  = (dlλpi - αλ*fdt - βλ*dxpi)^2 - (dlλci - αλ*fdt - βλ*dxci)^2
      llbmr += 0.5*dxci^2/(exp(0.5*(lσ2ci + lσ2ci1))*δt)                -
               0.5*dxpi^2/(exp(0.5*(lσ2pi + lσ2pi1))*δt)                +
               0.25*(lσ2ci + lσ2ci1 - lσ2pi - lσ2pi1)                   +
               0.5*((dlσ2p - ασ*fdt)^2 - (dlσ2c - ασ*fdt)^2)/(σσ^2*fdt) +
               ssλ0r*(-0.5/(σλ^2*fdt))
      llbr  += -fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
      dxsr  += (dxpi^2 - dxci^2)/fdt
      dxlr  += (dxpi * dlλpi - dxci * dlλci)/fdt
      ssσr  += (dlσ2p^2 - dlσ2c^2)/(2.0*fdt) + ασ*(dlσ2c - dlσ2p)
      ssλr  += ssλ0r/(2.0*fdt)
    end

    irλr = -llbr 

    #if speciation
    if λev
      llbr  += lλpi1 - lλci1
    end
  end

  return llbmr, llbr, dxsr, dxlr, ssσr, ssλr, irλr
end





"""
    ll_gibbs_xb!(tree::iTxb,
                 ασ  ::Float64,
                 σσ  ::Float64,
                 αλ  ::Float64,
                 βλ  ::Float64,
                 σλ  ::Float64,
                 ll  ::Float64,
                 dxs ::Float64,
                 dxl ::Float64,
                 ddx ::Float64,
                 ddσ ::Float64,
                 ssσ ::Float64,
                 ddλ ::Float64,
                 ssλ ::Float64,
                 nλ  ::Float64,
                 irλ ::Float64)

Returns likelihood and quantities for Gibbs sampling 
for trait driven speciation `iTxb`.
"""
function ll_gibbs_xb!(tree::iTxb,
                      ασ  ::Float64,
                      σσ  ::Float64,
                      αλ  ::Float64,
                      βλ  ::Float64,
                      σλ  ::Float64,
                      ll  ::Float64,
                      dxs ::Float64,
                      dxl ::Float64,
                      ddx ::Float64,
                      ddσ ::Float64,
                      ssσ ::Float64,
                      ddλ ::Float64,
                      ssλ ::Float64,
                      nλ  ::Float64,
                      irλ ::Float64,
                      ns  ::Float64,
                      L   ::Float64)

  id1 = def1(tree)

  ll0, dxs0, dxl0, ddx0, ddσ0, ssσ0, ddλ0, ssλ0, nλ0, irλ0 = 
    ll_gibbs_xb_b(xv(tree), lσ2(tree), lλ(tree), 
      ασ, σσ, αλ, βλ, σλ, dt(tree), fdt(tree), id1)

  ll  += ll0
  dxs += dxs0
  dxl += dxl0
  ddx += ddx0
  ddσ += ddσ0
  ssσ += ssσ0
  ddλ += ddλ0
  ssλ += ssλ0
  nλ  += nλ0
  irλ += irλ0
  L   += e(tree)

  if id1
    ll, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L = 
      ll_gibbs_xb!(tree.d1, ασ, σσ, αλ, βλ, σλ,
                   ll, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L)
    if def2(tree)
        ns += 1.0
        ll, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L = 
         ll_gibbs_xb!(tree.d2, ασ, σσ, αλ, βλ, σλ,
                      ll, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L)
    end
  end

  return ll, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L
end




"""
    ll_gibbs_xb_b(vx  ::Array{Float64,1},
                  vlσ2::Array{Float64,1},
                  vlλ ::Array{Float64,1},
                  ασ  ::Float64,
                  σσ  ::Float64,
                  αλ  ::Float64,
                  βλ  ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth
separately for the Brownian motion and the pure-birth
"""
function ll_gibbs_xb_b(vx  ::Array{Float64,1},
                       vlσ2::Array{Float64,1},
                       vlλ ::Array{Float64,1},
                       ασ  ::Float64,
                       σσ  ::Float64,
                       αλ  ::Float64,
                       βλ  ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       λev ::Bool)
  @inbounds begin

    # estimate standard `δt` likelihood
    nI = lastindex(vx) - 2 
    n  = Float64(nI)

    ll = llx = llσ2 = llλ = llpb = dxs = dxl = ssσ = ssλ = nλ = irλ = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        lσ2i  = vlσ2[i]
        lσ2i1 = vlσ2[i+1]
        dxi   = vx[i+1] - vx[i]
        dxi2  = dxi^2
        lλi   = vlλ[i]
        lλi1  = vlλ[i+1]
        dλi   = lλi1 - lλi
        dxs  += dxi2
        dxl  += dxi * dλi
        llx  += -0.5*dxi2/(exp(0.5*(lσ2i1 + lσ2i))*δt) - 
                  0.25*(lσ2i1 + lσ2i)
        llσ2 += (lσ2i1 - lσ2i - ασ*δt)^2
        llλ  += (dλi - αλ*δt - βλ*dxi)^2
        llpb += exp(0.5*(lλi + lλi1))
      end
      llpb *= δt
      dxs  /= δt
      dxl  /= δt
      ssσ  = llσ2/(2.0*δt)
      ssλ  = llλ /(2.0*δt)
      irλ  = llpb
      nλ  += Float64(nI)
      # add to global likelihood
      ll += llx - n*(1.5*log(δt) + log(σσ*σλ))                      +
            llσ2*(-0.5/(σσ^2*δt)) + llλ*(-0.5/(σλ^2*δt))            -
            n*2.756815599614018008622906563687138259410858154296875 - # 1.5 * log(2.0π)
            llpb
    end

    lλi1 = vlλ[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lσ2i  = vlσ2[nI+1]
      lσ2i1 = vlσ2[nI+2]
      dlσ2  = (lσ2i1 - lσ2i - ασ*fdt)^2
      dxi   = vx[nI+2] - vx[nI+1]
      dxi2  = dxi^2
      lλi   = vlλ[nI+1]
      dλi   = lλi1 - lλi
      dxs  += dxi2/fdt
      dxl  += dxi*dλi/fdt
      ssσ  += dlσ2/(2.0*fdt)
      dlλ2  = (dλi - αλ*fdt - βλ*dxi)^2
      ssλ  += dlλ2/(2.0*fdt)
      irλi  = exp(0.5*(lλi + lλi1))*fdt
      irλ  += irλi
      nλ   += 1.0

      # add to global likelihood
      ll += -0.5*dxi2/(exp(0.5*(lσ2i1 + lσ2i))*fdt)         - 
            0.25*(lσ2i1 + lσ2i) - 1.5*log(fdt) - log(σσ*σλ) +
            dlσ2 * (-0.5/(σσ^2*fdt))                        + 
            dlλ2 * (-0.5/(σλ^2*fdt))                        -
            2.756815599614018008622906563687138259410858154296875 - # 1.5 * log(2.0π)
            irλi
    end
    if λev
      ll += lλi1
    end
  end

  return ll, dxs, dxl, (vx[nI+2] -  vx[1]), (vlσ2[nI+2] - vlσ2[1]), 
         ssσ, (vlλ[nI+2] - vlλ[1]), ssλ, nλ, irλ
end




"""
    _gibbs_quanta!(tree::iTxb,
                   ασ  ::Float64,
                   αλ  ::Float64,
                   βλ  ::Float64,
                   dxs ::Float64,
                   dxl ::Float64,
                   ddx ::Float64,
                   ddσ ::Float64,
                   ssσ ::Float64,
                   ddλ ::Float64,
                   ssλ ::Float64,
                   nλ  ::Float64,
                   irλ ::Float64)

Returns the quantities for Gibbs sampling for trait driven speciation `iTxb`.
"""
function _gibbs_quanta!(tree::iTxb,
                        ασ  ::Float64,
                        αλ  ::Float64,
                        βλ  ::Float64,
                        dxs ::Float64,
                        dxl ::Float64,
                        ddx ::Float64,
                        ddσ ::Float64,
                        ssσ ::Float64,
                        ddλ ::Float64,
                        ssλ ::Float64,
                        nλ  ::Float64,
                        irλ ::Float64)

  dxs0, dxl0, ddx0, ddσ0, ssσ0, ddλ0, ssλ0, nλ0, irλ0 = 
    _gibbs_quanta(xv(tree), lσ2(tree), lλ(tree), 
                  ασ, αλ, βλ, dt(tree), fdt(tree))

  dxs += dxs0
  dxl += dxl0
  ddx += ddx0
  ddσ += ddσ0
  ssσ += ssσ0
  ddλ += ddλ0
  ssλ += ssλ0
  nλ  += nλ0
  irλ += irλ0

  if def1(tree)
      dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ = 
        _gibbs_quanta!(tree.d1, ασ, αλ, βλ, 
                       dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ)
    if def2(tree)
        dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ = 
         _gibbs_quanta!(tree.d2, ασ, αλ, βλ, 
                        dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ)
    end
  end

  return dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ
end




"""
     _gibbs_quanta(vx  ::Vector{Float64},
                   vlσ2::Vector{Float64},
                   vlλ ::Vector{Float64},
                   ασ  ::Float64,
                   αλ  ::Float64,
                   βλ  ::Float64,
                   δt  ::Float64,
                   fdt ::Float64)

Returns the quantities for Gibbs sampling for trait driven speciation `iTxb`.
"""
function _gibbs_quanta(vx  ::Vector{Float64},
                       vlσ2::Vector{Float64},
                       vlλ ::Vector{Float64},
                       ασ  ::Float64,
                       αλ  ::Float64,
                       βλ  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64)
  @inbounds begin

    nI = lastindex(vx) - 2
    dxs = dxl = ssσ = ssλ = nλ = irλ = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        dxi  = vx[i+1] - vx[i]
        lλi  = vlλ[i]
        lλi1 = vlλ[i+1]
        dλi  = lλi1 - lλi

        dxs += dxi^2
        dxl += dxi * dλi
        ssσ += (vlσ2[i+1] - vlσ2[i] - ασ*δt)^2
        ssλ += (dλi - αλ*δt - βλ*dxi)^2
        irλ += exp(0.5*(lλi + lλi1))
      end

      # standardize
      dxs /= δt
      dxl /= δt
      ssσ /= 2.0*δt
      ssλ /= 2.0*δt
      irλ *= δt
      nλ  += Float64(nI)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      dxi  = vx[nI+2] - vx[nI+1]
      lλi  = vlλ[nI+1]
      lλi1 = vlλ[nI+2]
      dλi  = lλi1 - lλi

      dxs += dxi^2/fdt
      dxl += dxi*dλi/fdt
      ssσ += (vlσ2[nI+2] - vlσ2[nI+1] - ασ*fdt)^2/(2.0*fdt)
      ssλ += (dλi - αλ*fdt - βλ*dxi)^2/(2.0*fdt)
      irλ += exp(0.5*(lλi + lλi1))*fdt
      nλ  += 1.0
    end
  end

  return dxs, dxl, 
        (vx[nI+2] -  vx[1]), (vlσ2[nI+2] - vlσ2[1]), ssσ, 
        (vlλ[nI+2] - vlλ[1]), ssλ, nλ, irλ
end



