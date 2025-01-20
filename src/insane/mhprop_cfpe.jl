#=

punkeek proposals

Ignacio Quintero Mächler

t(-_-t)

Created 25 11 2024
=#





"""
    _stem_update!(ξi       ::sTpe,
                  αx       ::Float64,
                  ασ       ::Float64,
                  γ        ::Float64,
                  δt       ::Float64,
                  srδt     ::Float64)

Perform punkeek stem update.
"""
function _fstem_update!(ξi ::sTfpe,
                        ξ1 ::sTfpe,
                        σa ::Float64,
                        ll ::Float64,
                        sσa::Float64)


  ei  = e(ξ1)
  xc  = xi(ξi)
  xfi = xf(ξ1)
  eσ  = σa*sqrt(ei)

  # sample mrca
  xp   = rnorm(xfi, eσ)
  ll  += lrdnorm_bm_x(xp, xc, xfi, eσ)
  sσa += (xfi - xp)^2/ei - (xfi - xc)^2/ei
  setxi!(ξi, xp)
  setxf!(ξi, xp)
  setxi!(ξ1, xp)

  return ll, sσa
end




"""
    _update_duo_x!(ξi  ::T,
                   ξ1  ::T,
                   σa  ::Float64,
                   ll  ::Float64,
                   sσa ::Float64) where {T <: Tpe}

Perform punkeek for **fixed** node.
"""
function _update_duo_x!(ξi  ::T,
                        ξ1  ::T,
                        xavg::Float64,
                        xstd::Float64,
                        σa  ::Float64,
                        ll  ::Float64,
                        sσa ::Float64) where {T <: Tpe}

  σa2 = σa^2
  xa, xic, x1 = xi(ξi), xf(ξi), xf(ξ1)
  ei, e1 = e(ξi), e(ξ1)

  # sample
  xip = trioprop(xavg, xa, x1, xstd^2, ei*σa2, e1*σa2)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa, ei*σa2) + llrdnorm_μ(x1, xip, xic, e1*σa2)
  sσa += ((xip - xa)^2 - (xic - xa)^2)/ei + ((x1 - xip)^2 - (x1 - xic)^2)/e1
  setxf!(ξi, xip)
  setxi!(ξ1, xip)

  return ll, sσa
end











