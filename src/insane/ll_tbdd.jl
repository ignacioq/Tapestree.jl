#=

trait pure-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_tb(Ξ   ::Vector{iTxb},
            idf ::Vector{iBffs},
            αλ  ::Float64,
            βλ  ::Float64,
            σλ  ::Float64,
            αx  ::Float64,
            ασ  ::Float64,
            σσ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Returns the log-likelihood for a `iTxb` according to trait pure-birth diffusion.
"""
function llik_tb(Ξ   ::Vector{iTxb},
                 idf ::Vector{iBffs},
                 αλ  ::Float64,
                 βλ  ::Float64,
                 σλ  ::Float64,
                 αx  ::Float64,
                 ασ  ::Float64,
                 σσ  ::Float64,
                 δt  ::Float64)

  @inbounds begin
    ll = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      ll += llik_tb(Ξ[i], αλ, βλ, σλ, αx, ασ, σσ, δt)
      if d2(idf[i]) > 0
        ll += λt(Ξ[i])
      end
    end
  end

  return ll
end




"""
    llik_tb(tree::iTxb,
            αλ  ::Float64,
            βλ  ::Float64,
            σλ  ::Float64,
            αx  ::Float64,
            ασ  ::Float64,
            σσ  ::Float64,
            δt  ::Float64,
            srδt::Float64)

Returns the log-likelihood for a `iTxb` according to trait pure-birth diffusion.
"""
function llik_tb(tree::iTxb,
                 αλ  ::Float64,
                 βλ  ::Float64,
                 σλ  ::Float64,
                 αx  ::Float64,
                 ασ  ::Float64,
                 σσ  ::Float64,
                 δt  ::Float64)

  if istip(tree)
    ll_tb_b(xv(tree), lσ2(tree), lλ(tree), 
            αλ, βλ, σλ, αx, ασ, σσ, δt, fdt(tree), false)
  else
    ll_tb_b(xv(tree), lσ2(tree), lλ(tree), 
            αλ, βλ, σλ, αx, ασ, σσ, δt, fdt(tree), true) +
    llik_tb(tree.d1::iTxb, αλ, βλ, σλ, αx, ασ, σσ, δt)   +
    llik_tb(tree.d2::iTxb, αλ, βλ, σλ, αx, ασ, σσ, δt)
  end
end




"""
    ll_tb_b(x   ::Array{Float64,1},
            lσ2v::Array{Float64,1},
            lλv ::Array{Float64,1},
            αλ  ::Float64,
            βλ  ::Float64,
            σλ  ::Float64,
            αx  ::Float64,
            ασ  ::Float64,
            σσ  ::Float64,
            δt  ::Float64,
            fdt ::Float64,
            λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_tb_b(x   ::Array{Float64,1},
                 lσ2v::Array{Float64,1},
                 lλv ::Array{Float64,1},
                 αλ  ::Float64,
                 βλ  ::Float64,
                 σλ  ::Float64,
                 αx  ::Float64,
                 ασ  ::Float64,
                 σσ  ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 λev ::Bool)

  @inbounds begin
    # estimate standard `δt` likelihood
    nI = lastindex(lλv)-2
    n  = Float64(nI)

    ll = llx = llσ2 = llλ = llpb = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        lσ2i   = lσ2v[i]
        lσ2i1  = lσ2v[i+1]
        dxi    = x[i+1] - x[i]
        lλi    = lλv[i]
        lλi1   = lλv[i+1]
        llx   += -0.5*(dxi - αx*δt)^2/(exp(0.5*(lσ2i1 + lσ2i))*δt) - 
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

    lλi1 = lλv[nI+2]

    # add final non-standard `δt`
    if fdt > 0.0
      lσ2i   = lσ2v[nI+1]
      lσ2i1  = lσ2v[nI+2]
      dxi    = x[nI+2] - x[nI+1]
      lλi    = lλv[nI+1]
      # add to global likelihood
      ll += -0.5*(dxi - αx*fdt)^2/(exp(0.5*(lσ2i1 + lσ2i))*fdt)   - 
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
    llik_tb_ssλ(tree::iTxb,
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

Returns the log-likelihood for a `iTxb` according to trait pure-birth diffusion.
"""
function llik_tb_ssλ(tree::iTxb,
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ns  ::Float64)

  if istip(tree)
    ll, dλ, ssλ, nλ, irλ = 
      ll_tb_b_ssλ(lλ(tree), α, σλ, δt, fdt(tree), srδt, false)
  else
    ns += 1.0

    ll, dλ, ssλ, nλ, irλ = 
      ll_tb_b_ssλ(lλ(tree), α, σλ, δt, fdt(tree), srδt, true)

    ll1, dλ1, ssλ1, nλ1, irλ1, ns = 
      llik_tb_ssλ(tree.d1, α, σλ, δt, srδt, ns)
    ll2, dλ2, ssλ2, nλ2, irλ2, ns = 
      llik_tb_ssλ(tree.d2, α, σλ, δt, srδt, ns)

    ll  += ll1  + ll2
    dλ  += dλ1  + dλ2
    ssλ += ssλ1 + ssλ2
    nλ  += nλ1  + nλ2
    irλ += irλ1 + irλ2
  end

  return ll, dλ, ssλ, nλ, irλ, ns
end




"""
    ll_tb_b_ssλ(lλv ::Array{Float64,1},
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 fdt ::Float64,
                 srδt::Float64,
                 λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth.
"""
function ll_tb_b_ssλ(lλv ::Array{Float64,1},
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      fdt ::Float64,
                      srδt::Float64,
                      λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλv)-2
  nλ = Float64(nI)

  ll = llbm = llb = ssλ = irλ = 0.0
  if nI > 0
    @turbo for i in Base.OneTo(nI)
      lλvi  = lλv[i]
      lλvi1 = lλv[i+1]
      llbm += (lλvi1 - lλvi - α*δt)^2
      llb += exp(0.5*(lλvi + lλvi1))
    end

    # standardized sum of squares
    ssλ += llbm/(2.0*δt)
    irλ += llb*δt

    # add to global likelihood
    ll += llbm*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) - 
          llb*δt
  end

  lλvi1 = lλv[nI+2]
  dλ    = lλvi1 - lλv[1]

  # add final non-standard `δt`
  if fdt > 0.0
    lλvi = lλv[nI+1]
    iri  = fdt*exp(0.5*(lλvi + lλvi1))
    ll  += ldnorm_bm(lλvi1, lλvi + α*fdt, sqrt(fdt)*σλ) -
           iri
    ssλ += (lλvi1 - lλvi - α*fdt)^2/(2.0*fdt)
    nλ  += 1.0
    irλ += iri
  end
  if λev
    ll += lλvi1
  end

  return ll, dλ, ssλ, nλ, irλ
end



"""
    llr_tb_b_sep(lλp ::Array{Float64,1},
                  lλc ::Array{Float64,1},
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  fdt ::Float64,
                  srδt::Float64,
                  λev ::Bool)

Returns the log-likelihood for a branch according to GBM pure-birth
separately for the Brownian motion and the pure-birth
"""
function llr_tb_b_sep(lλp ::Array{Float64,1},
                       lλc ::Array{Float64,1},
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       fdt ::Float64,
                       srδt::Float64,
                       λev ::Bool)

  # estimate standard `δt` likelihood
  nI = lastindex(lλp)-2

  llrbm = llrb = ssrλ = 0.0
  if nI > 0
    @turbo for i in Base.OneTo(nI)
      lλpi   = lλp[i]
      lλci   = lλc[i]
      lλpi1  = lλp[i+1]
      lλci1  = lλc[i+1]
      llrbm += (lλpi1 - lλpi - α*δt)^2 - (lλci1 - lλci - α*δt)^2
      llrb += exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1))
    end

    # standardized sum of squares
    ssrλ  += llrbm/(2.0*δt)
    # add to global likelihood
    llrbm *= (-0.5/((σλ*srδt)^2))
    llrb *= (-δt)
  end

  lλpi1 = lλp[nI+2]
  lλci1 = lλc[nI+2]

 # add final non-standard `δt`
  if fdt > 0.0
    lλpi   = lλp[nI+1]
    lλci   = lλc[nI+1]
    ssrλ  += ((lλpi1 - lλpi - α*fdt)^2 - (lλci1 - lλci - α*fdt)^2)/(2.0*fdt)
    llrbm += lrdnorm_bm_x(lλpi1, lλpi + α*fdt,
                          lλci1, lλci + α*fdt, sqrt(fdt)*σλ)
    llrb  -= fdt*(exp(0.5*(lλpi + lλpi1)) - exp(0.5*(lλci + lλci1)))
  end

  irrλ = -llrb 
  #if speciation
  if λev
    llrb  += lλpi1 - lλci1
  end

  return llrbm, llrb , ssrλ, irrλ
end




"""
    _ss_ir_dd(tree::T,
              f   ::Function,
              α   ::Float64,
              dd  ::Float64,
              ss  ::Float64,
              n   ::Float64,
              ir  ::Float64) where {T <: iTree}

Returns the standardized sum of squares for rate `v`, the path number `n`,
the integrated rate `ir` and the delta drift `dd`.
"""
function _ss_ir_dd(tree::T,
                   f   ::Function,
                   α   ::Float64,
                   dd  ::Float64,
                   ss  ::Float64,
                   n   ::Float64,
                   ir  ::Float64) where {T <: iTree}

  dd0, ss0, n0, ir0 = _ss_ir_dd_b(f(tree), α, dt(tree), fdt(tree))

  dd += dd0
  ss += ss0
  n  += n0
  ir += ir0

  if def1(tree)
    dd, ss, n, ir = _ss_ir_dd(tree.d1, f, α, dd, ss, n, ir)
    if def2(tree)
      dd, ss, n, ir = _ss_ir_dd(tree.d2, f, α, dd, ss, n, ir)
    end
  end

  return dd, ss, n, ir
end



"""
    _ss_ir_dd_b(v  ::Array{Float64,1},
                α  ::Float64,
                δt ::Float64,
                fdt::Float64)

Returns the standardized sum of squares for rate `v`, the path number `n`,
the integrated rate `ir` and the delta drift `dd`.
"""
function _ss_ir_dd_b(v  ::Array{Float64,1},
                     α  ::Float64,
                     δt ::Float64,
                     fdt::Float64)


    # estimate standard `δt` likelihood
    nI = lastindex(v)-2

    ss = ir = n = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        vi  = v[i]
        vi1 = v[i+1]
        ss += (vi1 - vi - α*δt)^2
        ir += exp(0.5*(vi + vi1))
      end
    
      # standardize
      ss *= 1.0/(2.0*δt)
      ir *= δt
      n  += Float64(nI)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      vi  = v[nI+1]
      vi1 = v[nI+2]
      ss += (vi1 - vi - α*fdt)^2/(2.0*fdt)
      n  += 1.0
      ir += fdt*exp(0.5*(vi + vi1))
    end

  return (v[nI+2] - v[1]), ss, n, ir
end



"""
    _ss(tree::T, f::Function, α::Float64) where {T <: iTree}

Returns the standardized sum of squares for rate `v`.
"""
function _ss(tree::T, f::Function, α::Float64) where {T <: iTree}

  ss = _ss_b(f(tree), α, dt(tree), fdt(tree))

  if def1(tree)
    ss += _ss(tree.d1, f, α)
    if def2(tree)
      ss += _ss(tree.d2, f, α)
    end
  end

  return ss
end




"""
    _ss_b(v::Array{Float64,1},
          α  ::Float64,
          δt ::Float64,
          fdt::Float64)

Returns the standardized sum of squares for rate `v`.
"""
function _ss_b(v::Array{Float64,1},
               α  ::Float64,
               δt ::Float64,
               fdt::Float64)

    # estimate standard `δt` likelihood
    nI = lastindex(v)-2

    ss = 0.0
    if nI > 0
      @turbo for i in Base.OneTo(nI)
        ss += (v[i+1] - v[i] - α*δt)^2
      end

      # standardize
      ss *= 1.0/(2.0*δt)
    end

    # add final non-standard `δt`
    if fdt > 0.0
      ss += (v[nI+2] - v[nI+1] - α*fdt)^2/(2.0*fdt)
    end

  return ss
end




"""
    int_rate(tree::iTree, f::Function)

Integrate rate given by `f`.
"""
function int_rate(tree::iTree, f::Function)

  if def1(tree)
    if def2(tree)
      int_rate(f(tree), dt(tree), fdt(tree)) +
      int_rate(tree.d1, f)         +
      int_rate(tree.d2, f)
    else
      int_rate(f(tree), dt(tree), fdt(tree)) +
      int_rate(tree.d1, f)
    end
  else
    int_rate(f(tree), dt(tree), fdt(tree))
  end
end




"""
    int_rate(v::Vector{Float64}, δt::Float64, fdt::Float64)

Integrate `v` rate.
"""
function int_rate(v::Vector{Float64}, δt::Float64, fdt::Float64)

  # estimate standard `δt` likelihood
  nI = lastindex(v)-2

  ir = 0.0
  if nI > 0
    @turbo for i in Base.OneTo(nI)
      ir += exp(0.5*(v[i] + v[i+1]))
    end
    ir *= δt
  end

  # add final non-standard `δt`
  if fdt > 0.0
    ir += fdt*exp(0.5*(v[nI+1] + v[nI+2]))
  end

  return ir
end



