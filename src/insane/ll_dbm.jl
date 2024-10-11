#=

Diffused Brownian motion likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#



"""
    llik_dbm(Ξ   ::Vector{sTxs},
             αx  ::Float64,
             ασ  ::Float64,
             γ   ::Float64,
             δt  ::Float64)

Returns the log-likelihood for `sTxs` according to diffused Brownian motion.
"""
function llik_dbm(Ξ   ::Vector{sTxs},
                  αx  ::Float64,
                  ασ  ::Float64,
                  γ   ::Float64,
                  δt  ::Float64)
  ll = 0.0
  for ξ in Ξ
    ll += llik_dbm(ξ, αx, ασ, γ, δt)
  end

  return ll
end




"""
    llik_dbm_v!(ll  ::Vector{Float64},
                 Ξ   ::Vector{sTxs},
                 αx  ::Float64,
                 ασ  ::Float64,
                 γ   ::Float64,
                 δt  ::Float64)

Returns the log-likelihood for `sTxs` according to diffused Brownian motion.
"""
function llik_dbm_v!(ll  ::Vector{Float64},
                     Ξ   ::Vector{sTxs},
                     αx  ::Float64,
                     ασ  ::Float64,
                     γ   ::Float64,
                     δt  ::Float64)

  @inbounds begin
    for i in Base.OneTo(lastindex(Ξ))
      ll[i] = llik_dbm(Ξ[i], αx, ασ, γ, δt)
    end
  end

  return nothing
end




"""
    llik_dbm(tree::sTxs,
             αx  ::Float64,
             ασ  ::Float64,
             γ   ::Float64,
             δt  ::Float64)

Returns the log-likelihood for a `sTxs` according to diffused Brownian motion.
"""
function llik_dbm(tree::sTxs,
                  αx  ::Float64,
                  ασ  ::Float64,
                  γ   ::Float64,
                  δt  ::Float64)

  if def1(tree)
    if def2(tree)
      llik_dbm(tree.d1, αx, ασ, γ, δt)                        +
      llik_dbm(tree.d2, αx, ασ, γ, δt)                        +
      ll_dbm_b(xv(tree), αx, lσ2(tree), ασ, γ, δt, fdt(tree))
    else
      llik_dbm(tree.d1, αx, ασ, γ, δt)                        +
      ll_dbm_b(xv(tree), αx, lσ2(tree), ασ, γ, δt, fdt(tree))
    end
  else
    ll_dbm_b(xv(tree), αx, lσ2(tree), ασ, γ, δt, fdt(tree))
  end
end



"""
    ll_dbm_b(x   ::Array{Float64,1},
             lσ2 ::Array{Float64,1},
             αx  ::Float64, 
             ασ  ::Float64,
             γ   ::Float64,
             δt  ::Float64,
             fdt ::Float64)

Returns the log-likelihood for a branch according to diffused Brownian motion.
"""
function ll_dbm_b(x   ::Array{Float64,1},
                  αx  ::Float64, 
                  lσ2 ::Array{Float64,1},
                  ασ  ::Float64,
                  γ   ::Float64,
                  δt  ::Float64,
                  fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2
    n  = Float64(nI)

    ll = llx = llσ2 = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        lσ2i   = lσ2[i]
        lσ2i1  = lσ2[i+1]
        llσ2  += (lσ2i1 - lσ2i - ασ*δt)^2
        llx   += -0.5*(x[i+1] - x[i] - αx*δt)^2/(exp(0.5*(lσ2i1 + lσ2i))*δt) - 
                  0.25*(lσ2i1 + lσ2i)
      end
      # estimate global likelihood
      ll += llx - 0.5*n*log(δt) + 
            llσ2*(-0.5/(γ^2*δt)) - 0.5*n*log(γ^2*δt) - 
            n*1.83787706640934533908193770912475883960723876953125 # log(2.0π)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      lσ2i  = lσ2[nI+1]
      lσ2i1 = lσ2[nI+2]
      ll   += -0.5*(x[nI+2] - x[nI+1] - αx*fdt)^2/(exp(0.5*(lσ2i1 + lσ2i))*fdt) - 
               0.25*(lσ2i1 + lσ2i) - 0.5*log(fdt) + 
               (lσ2i1 - lσ2i - ασ*fdt)^2*(-0.5/(γ^2*fdt)) - 0.5*log(γ^2*fdt) - 
               1.83787706640934533908193770912475883960723876953125
    end
  end

  return ll
end



"""
    ll_dbm_ss_dd_b(x   ::Array{Float64,1},
                   αx  ::Float64, 
                   lσ2 ::Array{Float64,1},
                   ασ  ::Float64,
                   γ   ::Float64,
                   δt  ::Float64,
                   fdt ::Float64)

Returns the log-likelihood for a branch according to diffused Brownian motion.
"""
function ll_dbm_ss_dd_b(x   ::Array{Float64,1},
                        αx  ::Float64, 
                        lσ2 ::Array{Float64,1},
                        ασ  ::Float64,
                        γ   ::Float64,
                        δt  ::Float64,
                        fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2
    n  = Float64(nI)

    ll = llx = llσ2 = Ls = Xs = ss = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        lσ2i  = lσ2[i]
        lσ2i1 = lσ2[i+1]
        llσ2 += (lσ2i1 - lσ2i - ασ*δt)^2
        mσ    = exp(0.5*(lσ2i1 + lσ2i))
        dx    = x[i+1] - x[i]
        Ls   += 1.0/mσ
        Xs   += dx/mσ
        llx  += -0.5*(dx - αx*δt)^2/(mσ*δt) - 0.25*(lσ2i1 + lσ2i)
      end
      Ls *= δt
      # estimate standard squares
      ss += llσ2/(2.0*δt)
      # estimate global likelihood
      ll += llx - 0.5*n*log(δt) + 
            llσ2*(-0.5/(γ^2*δt)) - 0.5*n*log(γ^2*δt) - 
            n*1.83787706640934533908193770912475883960723876953125 # log(2.0π)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      lσ2i  = lσ2[nI+1]
      lσ2i1 = lσ2[nI+2]
      dlσ2  = (lσ2i1 - lσ2i - ασ*fdt)^2
      mσ    = exp(0.5*(lσ2i1 + lσ2i))
      dx    = x[nI+2] - x[nI+1]
      Ls   += fdt/mσ
      Xs   += dx/mσ
      ll   += -0.5*(dx - αx*fdt)^2/(mσ*fdt) - 0.25*(lσ2i1 + lσ2i) - 
               0.5*log(fdt) + dlσ2*(-0.5/(γ^2*fdt)) - 0.5*log(γ^2*fdt) - 
               1.83787706640934533908193770912475883960723876953125
      ss += dlσ2/(2.0*fdt)
    end
  end

  return ll, Ls, Xs, (lσ2[nI+2] - lσ2[1]), ss
end




"""
    llr_dbm_σ(x   ::Array{Float64,1},
              αx  ::Float64, 
              lσ2p::Array{Float64,1},
              lσ2c::Array{Float64,1},
              δt  ::Float64,
              fdt ::Float64)

Returns the acceptance ratio for a `σ²(t)` path proposal (the likelihood for 
the GBM for `σ²` cancels out).
"""
function llr_dbm_σ(x   ::Array{Float64,1},
                   αx  ::Float64, 
                   lσ2p::Array{Float64,1},
                   lσ2c::Array{Float64,1},
                   δt  ::Float64,
                   fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(x)-2

    acr = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        lσ2ci  = lσ2c[i]
        lσ2ci1 = lσ2c[i+1]
        lσ2pi  = lσ2p[i]
        lσ2pi1 = lσ2p[i+1]
        acr += -(0.5*(x[i+1] - x[i] - αx*δt)^2/δt) *
                (1.0/exp(0.5*(lσ2pi + lσ2pi1)) - 1.0/exp(0.5*(lσ2ci + lσ2ci1))) + 
                 0.25*(lσ2ci + lσ2ci1 - lσ2pi - lσ2pi1)
      end
    end
    # add final non-standard `δt`
    if fdt > 0.0
      lσ2ci  = lσ2c[nI+1]
      lσ2ci1 = lσ2c[nI+2]
      lσ2pi  = lσ2p[nI+1]
      lσ2pi1 = lσ2p[nI+2]
      acr += -(0.5*(x[nI+2] - x[nI+1] - αx*fdt)^2/fdt) *
              (1.0/exp(0.5*(lσ2pi + lσ2pi1)) - 1.0/exp(0.5*(lσ2ci + lσ2ci1))) + 
               0.25*(lσ2ci + lσ2ci1 - lσ2pi - lσ2pi1)
    end
  end

  return acr
end




"""
  llr_scale(Ξ   ::Vector{sTxs},
            αx  ::Float64,
            s   ::Vector{Float64},
            δt  ::Float64)

Returns the log-likelihood ratio for a `σ(t)` scaled proposal.
"""
function llr_scale(Ξ   ::Vector{sTxs},
                   αx  ::Float64,
                   s   ::Float64,
                   δt  ::Float64)

  @inbounds begin
    xn  = lastindex(Ξ)
    llr = zeros(xn)
    Lsr = zeros(xn)
    Xsr = zeros(xn)
    for i in Base.OneTo(xn)
      ξ = Ξ[i]
      llr[i], Lsr[i], Xsr[i] = llr_scale(xv(ξ), αx, lσ2(ξ), s, δt, fdt(ξ))
    end
  end

  return llr, Lsr, Xsr
end




"""
    llr_scale(x   ::Array{Float64,1},
              αx  ::Float64, 
              lσ2c::Array{Float64,1},
              s   ::Float64,
              δt  ::Float64,
              fdt ::Float64)

Returns the log-likelihood ratio for a `σ(t)` scaled proposal.
"""
function llr_scale(x   ::Array{Float64,1},
                   αx  ::Float64, 
                   lσ2c::Array{Float64,1},
                   s   ::Float64,
                   δt  ::Float64,
                   fdt ::Float64)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI  = lastindex(x)-2
    es = exp(s)

    llr = Lsr = Xsr = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        σca  = exp(0.5*(lσ2c[i] + lσ2c[i+1]))
        σpa  = σca * es
        iσd  = (1.0/σpa - 1.0/σca)
        dxi  = x[i+1] - x[i] 
        Lsr += iσd
        Xsr += dxi * iσd
        llr -= 0.5*(dxi - αx*δt)^2/δt * iσd
      end
      Lsr *= δt
      llr -= 0.5*Float64(nI)*s
    end
    # add final non-standard `δt`
    if fdt > 0.0
      σca  = exp(0.5*(lσ2c[nI+1] + lσ2c[nI+2]))
      σpa  = σca * es
      iσd  = (1.0/σpa - 1.0/σca)
      dxi  = x[nI+2] - x[nI+1] 
      Lsr += fdt * iσd
      Xsr += dxi * iσd
      llr += -(0.5*(dxi - αx*fdt)^2/fdt)*iσd - 0.5*s
    end
  end

  return llr, Lsr, Xsr
end




"""
    _ss_dd(tree::T,
           fx  ::Function,
           fσ  ::Function,
           ασ  ::Float64,
           Ls  ::Float64,
           Xs  ::Float64,
           ddσ ::Float64,
           ss  ::Float64,
           n   ::Float64) where {T <: iTree}

Returns the standardized sum of squares for rate `fσ`, the path number `n`, 
and the delta drift `ddσ`, as well as Gibbs variables for drift in `fx`.
"""
function _ss_dd(tree::T,
                fx  ::Function,
                fσ  ::Function,
                ασ  ::Float64,
                Ls  ::Float64,
                Xs  ::Float64,
                ddσ ::Float64,
                ss  ::Float64,
                n   ::Float64) where {T <: iTree}

  Ls0, Xs0, ddσ0, ss0, n0 = 
    _ss_dd_b(fx(tree), fσ(tree), ασ, dt(tree), fdt(tree))

  Ls  += Ls0
  Xs  += Xs0
  ddσ += ddσ0
  ss  += ss0
  n   += n0

  if def1(tree)
    Ls, Xs, ddσ, ss, n = _ss_dd(tree.d1, fx, fσ, ασ, Ls, Xs, ddσ, ss, n)
    if def2(tree)
      Ls, Xs, ddσ, ss, n = _ss_dd(tree.d2, fx, fσ, ασ, Ls, Xs, ddσ, ss, n)
    end
  end

  return Ls, Xs, ddσ, ss, n
end




"""
    _ss_dd_b(vx ::Array{Float64,1},
             vσ ::Array{Float64,1},
             αx ::Float64,
             ασ ::Float64,
             δt ::Float64,
             fdt::Float64)

Returns the standardized sum of squares for vectors `vx` & `vσ`, 
the path number `n`, and the delta drift `ddx` & `ddσ`.
"""
function _ss_dd_b(vx ::Array{Float64,1},
                  vσ ::Array{Float64,1},
                  ασ ::Float64,
                  δt ::Float64,
                  fdt::Float64)


    # estimate standard `δt` likelihood
    nI = lastindex(vσ)-2

    Ls = Xs = ss = n = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        vi  = vσ[i]
        vi1 = vσ[i+1]
        ivm = 1.0/exp(0.5*(vi + vi1))
        Ls += ivm
        Xs += (vx[i+1] - vx[i]) * ivm
        ss += (vi1 - vi - ασ*δt)^2
      end

      # standardize by time
      Ls *= δt
      ss *= 1.0/(2.0*δt)
      n  += Float64(nI)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      vi  = vσ[nI+1]
      vi1 = vσ[nI+2]
      ivm = 1.0/exp(0.5*(vi + vi1))
      Ls += fdt * ivm
      Xs += (vx[nI+2] - vx[nI+1]) * ivm
      ss += (vi1 - vi - ασ*fdt)^2/(2.0*fdt)
      n  += 1.0
    end

  return Ls, Xs, (vσ[nI+2] - vσ[1]), ss, n
end


