#=

punkeek proposals

Ignacio Quintero Mächler

t(-_-t)

Created 25 11 2024
=#



"""
     _stem_update!(ξi ::sTpe,
                   σa ::Float64,
                   ll ::Float64,
                   sσa::Float64)

Perform punkeek stem update.
"""
function _stem_update!(ξi ::sTpe,
                       σa ::Float64,
                       ll ::Float64,
                       sσa::Float64)
  ei  = e(ξi)
  xc  = xi(ξi)
  xfi = xf(ξi)
  eσ  = σa*sqrt(ei)

  # sample mrca
  xp   = rnorm(xfi, eσ)
  ll  += lrdnorm_bm_x(xp, xc, xfi, eσ)
  sσa += ((xfi - xp)^2 - (xfi - xc)^2)/(2.0*ei)
  setxi!(ξi, xp)

  return ll, sσa
end




"""
    _crown_update!(ξi ::sTpe,
                   ξ1 ::sTpe,
                   ξ2 ::sTpe,
                   σa ::Float64,
                   σk ::Float64,
                   ll ::Float64,
                   sσa::Float64,
                   sσk::Float64) 

Perform punkeek crown update.
"""
function _crown_update!(ξi ::sTpe,
                        ξ1 ::sTpe,
                        ξ2 ::sTpe,
                        σa ::Float64,
                        σk ::Float64,
                        ll ::Float64,
                        sσa::Float64,
                        sσk::Float64)

  σa2, σk2 = σa^2, σk^2
  xic, x1, x2 = xf(ξi), xf(ξ1), xf(ξ2)
  e1, e2 = e(ξ1), e(ξ2)

  ## which cladogenetic (unidentifiable at the root)
  # if d1 cladogenetic
  pk1 = llik_cpe_dyad(xic, x2, x1, e2, e1, σa2, σk2)
  # d2 cladogenetic
  pk2 = llik_cpe_dyad(xic, x1, x2, e1, e2, σa2, σk2)
  dpk = pk1 - pk2

  p1 = if dpk > 37.0
    1.0
  else
    o12 = exp(dpk)  # odds
    o12/(1.0 + o12) # probability
  end

  # p1 if x1 cladogenetic
  shp, xadp, xkdp, eadp, ekdp = 
    rand() < p1 ? (true, x2, x1, e2, e1) : (false, x1, x2, e1, e2)

  # sample new values
  xip = duoprop(xadp, xkdp, eadp*σa2, ekdp*σa2 + σk2)
  xkp = duoprop(xip, xkdp, σk2, ekdp*σa2)

  # update likelihood and gibbs quanta
  xkc, xadc, xkdc, eadc, ekdc = 
    sh(ξi) ? (xi(ξ1), x2, x1, e2, e1) : (xi(ξ2), x1, x2, e1, e2)

  ll += llik_cpe_trio(xip, xkp, xadp, xkdp, eadp, ekdp, σa2, σk2) -
        llik_cpe_trio(xic, xkc, xadc, xkdc, eadc, ekdc, σa2, σk2)

  sσa += 0.5*((xadp - xip)^2/eadp + (xkdp - xkp)^2/ekdp -
              (xadc - xic)^2/eadc - (xkdc - xkc)^2/ekdc)
  sσk += 0.5*((xkp - xip)^2 - (xkc - xic)^2)

  setsh!(ξi, shp)
  setxi!(ξi, xip)
  setxf!(ξi, xip)
  ξk, ξa = shp ? (ξ1, ξ2) : (ξ2, ξ1)
  setxi!(ξk, xkp) 
  setxi!(ξa, xip)

  return ll, sσa, sσk
end




"""
    _update_node!(tree::sTpe,
                  xavg::Float64,
                  xstd::Float64,
                  σa  ::Float64,
                  σk  ::Float64,
                  ll  ::Float64,
                  sσa ::Float64,
                  sσk ::Float64)

Perform punkeek internal node updates.
"""
function _update_node!(tree::sTpe,
                       xavg::Float64,
                       xstd::Float64,
                       σa  ::Float64,
                       σk  ::Float64,
                       ll  ::Float64,
                       sσa ::Float64,
                       sσk ::Float64,
                       ter ::Bool)
  if def1(tree)
    ll, sσa, sσk = _update_quartet!(tree, σa, σk, ll, sσa, sσk)

    ll, sσa, sσk = 
      _update_node!(tree.d1, xavg, xstd, σa, σk, ll, sσa, sσk, ter)
    ll, sσa, sσk = 
      _update_node!(tree.d2, xavg, xstd, σa, σk, ll, sσa, sσk, ter)
  else
    if !isfix(tree)
      ll, sσa = _update_tip!(tree, NaN, NaN, σa, ll, sσa)
    else
      if ter && isnan(xavg)
        ll, sσa = _update_tip!(tree, xavg, xstd, σa, ll, sσa)
      end
    end
  end

  return ll, sσa, sσk
end




"""
    _update_tip!(tree::sTpe,
                  σa  ::Float64, 
                  σk  ::Float64, 
                  ll  ::Float64, 
                  sσa ::Float64, 
                  sσk ::Float64) 

Perform punkeek tip updates.
"""
function _update_tip!(tree::sTpe,
                      xavg::Float64,
                      xstd::Float64,
                      σa  ::Float64, 
                      ll  ::Float64, 
                      sσa ::Float64)

  xa, xfc = xi(tree), xf(tree)
  ei = e(tree)

  # trait proposal
  xfp = NaN
  if xstd > 0.0
    xfp = duoprop(xavg, xa, xstd^2, ei*σa^2)
  else
    xfp = rnorm(xa, sqrt(ei)*σa)
  end

  ## update trackers
  ll  += llrdnorm_x(xfp, xfc, xa, ei*σa^2)
  sσa += ((xfp - xa)^2 - (xfc - xa)^2)/(2.0*ei)
  setxf!(tree, xfp)

  return ll, sσa
end




"""
    _update_duo!(ξi  ::sTpe,
                   ξ1  ::sTpe,
                   σa  ::Float64,
                   ll  ::Float64,
                   sσa ::Float64)

Perform punkeek for **unfixed** node.
"""
function _update_duo!(ξi  ::sTpe,
                      ξ1  ::sTpe,
                      σa  ::Float64,
                      ll  ::Float64,
                      sσa ::Float64)

  σa2 = σa^2
  xa, xic, x1 = xi(ξi), xf(ξi), xf(ξ1)
  ei, e1 = e(ξi), e(ξ1)

  # sample
  xip = duoprop(xa, x1, ei*σa2, e1*σa2)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa, ei*σa2) + llrdnorm_μ(x1, xip, xic, e1*σa2)
  sσa += 0.5*(((xip - xa)^2 - (xic - xa)^2)/ei + 
              ((x1 - xip)^2 - (x1 - xic)^2)/e1)
  setxf!(ξi, xip)
  setxi!(ξ1, xip)

  return ll, sσa
end




"""
    _update_quartet!(ξi ::sTpe,
                       σa ::Float64,
                       σk ::Float64,
                       ll ::Float64,
                       sσa::Float64,
                       sσk::Float64)

Make a punkeek quartet proposal.
"""
function _update_quartet!(ξi ::sTpe,
                          σa ::Float64,
                          σk ::Float64,
                          ll ::Float64,
                          sσa::Float64,
                          sσk::Float64)

  ll, sσa, sσk = _update_quartet!(ξi, ξi.d1, ξi.d2, σa, σk, ll, sσa, sσk)

  return ll, sσa, sσk
end




"""
    _update_quartet!(ξi ::sTpe,
                     ξ1 ::sTpe,
                     ξ2 ::sTpe,
                     σa ::Float64,
                     σk ::Float64,
                     ll ::Float64,
                     sσa::Float64,
                     sσk::Float64)

Perform a punkeek quartet update.
"""
function _update_quartet!(ξi ::sTpe,
                          ξ1 ::sTpe,
                          ξ2 ::sTpe,
                          σa ::Float64,
                          σk ::Float64,
                          ll ::Float64,
                          sσa::Float64,
                          sσk::Float64)

  σa2, σk2 = σa^2, σk^2
  xa, xic, x1, x2 = xi(ξi), xf(ξi), xf(ξ1), xf(ξ2)
  ei, e1, e2 = e(ξi), e(ξ1), e(ξ2)

  ## probabilities
  # which is cladogenetic
  pk1 = llik_cpe_triad(xa, x2, x1, ei, e2, e1, σa2, σk2) # d1 cladogenetic
  pk2 = llik_cpe_triad(xa, x1, x2, ei, e1, e2, σa2, σk2) # d2 cladogenetic
  dpk = pk1 - pk2

  p1 = if dpk > 37.0
    1.0
  else
    o12 = exp(dpk)  # odds
    o12/(1.0 + o12) # probability
  end

  # p1 if x1 cladogenetic
  shp, xadp, xkdp, eadp, ekdp = 
    rand() < p1 ? (true, x2, x1, e2, e1) : (false, x1, x2, e1, e2)

  # sample new values
  xip = trioprop(xa, xadp, xkdp, ei*σa2, eadp*σa2, ekdp*σa2 + σk2)
  xkp = duoprop(xip, xkdp, σk2, ekdp*σa2)

  # update likelihood and gibbs quanta
  xkc, xadc, xkdc, eadc, ekdc = 
    sh(ξi) ? (xi(ξ1), x2, x1, e2, e1) : (xi(ξ2), x1, x2, e1, e2)

  # look if need to control for the proposal probability
  ll += llik_cpe_quartet(xa, xip, xkp, xadp, xkdp, ei, eadp, ekdp, σa2, σk2) -
        llik_cpe_quartet(xa, xic, xkc, xadc, xkdc, ei, eadc, ekdc, σa2, σk2)

  sσa += 0.5*((xip - xa)^2/ei + (xadp - xip)^2/eadp + (xkdp - xkp)^2/ekdp -
              (xic - xa)^2/ei - (xadc - xic)^2/eadc - (xkdc - xkc)^2/ekdc)
  sσk += 0.5*((xkp - xip)^2 - (xkc - xic)^2)

  setsh!(ξi, shp)
  setxf!(ξi, xip)

  ξk, ξa = shp ? (ξ1, ξ2) : (ξ2, ξ1)
  setxi!(ξk, xkp) 
  setxi!(ξa, xip)

  return ll, sσa, sσk
end





