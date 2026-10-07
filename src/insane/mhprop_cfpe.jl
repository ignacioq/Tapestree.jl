#=

fossil punkeek proposals

Ignacio Quintero Mächler

t(-_-t)

Created 25 11 2024
=#




"""
     _stem_update!(ξi ::sTfpe,
                   α  ::Float64,
                   σa ::Float64,
                   ll ::Float64,
                   dα ::Float64,
                   sσa::Float64)

Perform punkeek stem update.
"""
function _stem_update!(ξi ::sTfpe,
                       α  ::Float64,
                       σa ::Float64,
                       ll ::Float64,
                       dα ::Float64,
                       sσa::Float64)
  ei  = e(ξi)
  xc  = xi(ξi)
  xfi = xf(ξi)
  eσ  = σa*sqrt(ei)

  # sample mrca
  xp   = rnorm(xfi - α*ei, eσ)
  ll  += lrdnorm_bm_x(xp, xc - α*ei, xfi, eσ)
  dα  += xc - xp
  sσa += ((xfi - xp - α*ei)^2 - (xfi - xc - α*ei)^2)/(2.0*ei)
  setxi!(ξi, xp)

  return ll, dα, sσa
end




"""
    _fstem_update!(ξi  ::sTfpe,
                   ξ1  ::sTfpe,
                   α   ::Float64,
                   σa  ::Float64,
                   ll  ::Float64,
                   dα  ::Float64,
                   sσa ::Float64,
                   ifxi::Bool)

Perform punkeek stem fossil update.
"""
function _fstem_update!(ξi  ::sTfpe,
                        ξ1  ::sTfpe,
                        α   ::Float64,
                        σa  ::Float64,
                        ll  ::Float64,
                        dα  ::Float64,
                        sσa ::Float64)
  e1  = e(ξ1)
  xc  = xi(ξi)
  xf1 = xf(ξ1)
  eσ  = σa*sqrt(e1)

  # sample mrca
  xp   = rnorm(xf1 - α*e1, eσ)
  ll  += lrdnorm_bm_x(xp, xc, xf1 - α*e1, eσ)
  dα  += xc - xp

  sσa += 0.5*((xf1 - xp - α*e1)^2 - (xf1 - xc - α*e1)^2)/e1
  setxi!(ξi, xp)
  setxf!(ξi, xp)
  setxi!(ξ1, xp)

  return ll, dα, sσa
end




"""
    _crown_update!(ξi ::sTfpe,
                   ξ1 ::sTfpe,
                   ξ2 ::sTfpe,
                   α  ::Float64,
                   σa ::Float64,
                   σk ::Float64,
                   ll ::Float64,
                   dα ::Float64,
                   sσa::Float64,
                   sσk::Float64)

Perform punkeek crown update.
"""
function _crown_update!(ξi ::sTfpe,
                        ξ1 ::sTfpe,
                        ξ2 ::sTfpe,
                        α  ::Float64,
                        σa ::Float64,
                        σk ::Float64,
                        ll ::Float64,
                        dα ::Float64,
                        sσa::Float64,
                        sσk::Float64)

  σa2, σk2 = σa^2, σk^2
  xic, x1, x2 = xf(ξi), xf(ξ1), xf(ξ2)
  e1, e2 = e(ξ1), e(ξ2)

  ## which cladogenetic (unidentifiable at the root)
  # d1 cladogenetic
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
  xip = duoprop(xadp - α*eadp, xkdp - α*ekdp, eadp*σa2, ekdp*σa2 + σk2)
  xkp = duoprop(xip, xkdp - α*ekdp, σk2, ekdp*σa2)

  # update likelihood and gibbs quanta
  xkc, xadc, xkdc, eadc, ekdc = 
    sh(ξi) ? (xi(ξ1), x2, x1, e2, e1) : (xi(ξ2), x1, x2, e1, e2)

  ll += llik_cpe_trio(xip, xkp, xadp - α*eadp, xkdp - α*ekdp, 
          eadp, ekdp, σa2, σk2) -
        llik_cpe_trio(xic, xkc, xadc - α*eadc, xkdc - α*ekdc, 
          eadc, ekdc, σa2, σk2)

  dα  += (xadp - xip) + (xkdp - xkp) - 
         (xadc - xic) - (xkdc - xkc)
  sσa += 0.5*((xadp - xip - α*eadp)^2/eadp + 
              (xkdp - xkp - α*ekdp)^2/ekdp -
              (xadc - xic - α*eadc)^2/eadc - 
              (xkdc - xkc - α*ekdc)^2/ekdc)
  sσk += 0.5*((xkp - xip)^2 - (xkc - xic)^2)

  setsh!(ξi, shp)
  setxi!(ξi, xip)
  setxf!(ξi, xip)
  ξk, ξa = shp ? (ξ1, ξ2) : (ξ2, ξ1)
  setxi!(ξk, xkp) 
  setxi!(ξa, xip)

  return ll, dα, sσa, sσk
end




"""
    _update_node!(tree::sTfpe,
                  xavg::Float64,
                  xstd::Float64,
                  α   ::Float64,
                  σa  ::Float64,
                  σk  ::Float64,
                  ll  ::Float64,
                  dα  ::Float64,
                  sσa ::Float64,
                  sσk ::Float64,
                  ter ::Bool)

Perform punkeek node updates.
"""
function _update_node!(tree::sTfpe,
                       xavg::Float64,
                       xstd::Float64,
                       α   ::Float64,
                       σa  ::Float64,
                       σk  ::Float64,
                       ll  ::Float64,
                       dα  ::Float64,
                       sσa ::Float64,
                       sσk ::Float64,
                       ter ::Bool)

  if def1(tree)
    if def2(tree)
      ll, dα, sσa, sσk = _update_quartet!(tree, α, σa, σk, ll, dα, sσa, sσk)

      ll, dα, sσa, sσk = 
        _update_node!(tree.d1, xavg, xstd, α, σa, σk, ll, dα, sσa, sσk, ter)
      ll, dα, sσa, sσk = 
        _update_node!(tree.d2, xavg, xstd, α, σa, σk, ll, dα, sσa, sσk, ter)
    else
      if xstd > 0.0
        ll, sσa = _update_duo!(tree, xavg, xstd, α, σa, ll, sσa)
      end
      ll, dα, sσa, sσk = 
        _update_node!(tree.d1, xavg, xstd, α, σa, σk, ll, dα, sσa, sσk, ter)
    end
  else
    if !isfix(tree)
      ll, dα, sσa = _update_tip!(tree, NaN, NaN, α, σa, ll, dα, sσa)
    else
      if ter && isnan(xavg)
        ll, dα, sσa = _update_tip!(tree, xavg, xstd, α, σa, ll, dα, sσa)
      end
    end
  end
  return ll, dα, sσa, sσk
end




"""
    _update_tip!(tree::sTfpe,
                 xavg::Float64,
                 xstd::Float64,
                 α   ::Float64,
                 σa  ::Float64, 
                 ll  ::Float64, 
                 dα  ::Float64,
                 sσa ::Float64)

Perform punkeek tip updates.
"""
function _update_tip!(tree::sTfpe,
                      xavg::Float64,
                      xstd::Float64,
                      α   ::Float64,
                      σa  ::Float64, 
                      ll  ::Float64, 
                      dα  ::Float64,
                      sσa ::Float64)

  xa, xfc = xi(tree), xf(tree)
  ei = e(tree)

  # trait proposal
  xfp = NaN
  if xstd > 0.0
    xfp = duoprop(xavg, xa + α*ei, xstd^2, ei*σa^2)
  else
    xfp = rnorm(xa + α*ei, sqrt(ei)*σa)
  end

  ## update trackers
  ll  += llrdnorm_x(xfp, xfc, xa + α*ei, ei*σa^2)
  dα  += xfp - xfc
  sσa += 0.5*((xfp - xa - α*ei)^2 - (xfc - xa - α*ei)^2)/ei
  setxf!(tree, xfp)

  return ll, dα, sσa
end




"""
    _update_duo!(ξi ::sTfpe,
                 xavg::Float64,
                 xstd::Float64,
                 α  ::Float64,
                 σa ::Float64,
                 ll ::Float64,
                 sσa::Float64)

Make a punkeek dup proposal.
"""
function _update_duo!(ξi ::sTfpe,
                      xavg::Float64,
                      xstd::Float64,
                      α  ::Float64,
                      σa ::Float64,
                      ll ::Float64,
                      sσa::Float64)

  ll, sσa = _update_duo!(ξi, ξi.d1, xavg, xstd, α, σa, ll, sσa)

  return ll, sσa
end





"""
    _update_duo!(ξi  ::sTfpe,
                 ξ1  ::sTfpe,
                 xavg::Float64,
                 xstd::Float64,
                 α   ::Float64,
                 σa  ::Float64,
                 ll  ::Float64,
                 sσa ::Float64)

Perform punkeek update for duo node.
"""
function _update_duo!(ξi  ::sTfpe,
                      ξ1  ::sTfpe,
                      xavg::Float64,
                      xstd::Float64,
                      α   ::Float64,
                      σa  ::Float64,
                      ll  ::Float64,
                      sσa ::Float64)

  σa2 = σa^2
  xa, xic, x1 = xi(ξi), xf(ξi), xf(ξ1)
  ei, e1 = e(ξi), e(ξ1)

  # sample
  xip = NaN
  if xstd > 0.0
    xip = trioprop(xavg, xa + α*ei, x1 - α*e1, xstd^2, ei*σa2, e1*σa2)
  else
    xip = duoprop(       xa + α*ei, x1 - α*e1,         ei*σa2, e1*σa2)
  end

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa + α*ei, ei*σa2) + 
         llrdnorm_μ(x1 - α*e1, xip, xic, e1*σa2)
  sσa += 0.5*(((xip - xa - α*ei)^2 - (xic - xa - α*ei)^2)/ei + 
              ((x1 - xip - α*e1)^2 - (x1 - xic - α*e1)^2)/e1)
  setxf!(ξi, xip)
  setxi!(ξ1, xip)

  return ll, sσa
end




"""
    _update_quartet!(ξi ::sTfpe,
                     α  ::Float64,
                     σa ::Float64,
                     σk ::Float64,
                     ll ::Float64,
                     dα ::Float64,
                     sσa::Float64,
                     sσk::Float64)

Make a punkeek quartet proposal.
"""
function _update_quartet!(ξi ::sTfpe,
                          α  ::Float64,
                          σa ::Float64,
                          σk ::Float64,
                          ll ::Float64,
                          dα ::Float64,
                          sσa::Float64,
                          sσk::Float64)

  ll, dα, sσa, sσk = 
    _update_quartet!(ξi, ξi.d1, ξi.d2, α, σa, σk, ll, dα, sσa, sσk)

  return ll, dα, sσa, sσk
end




"""
    _update_quartet!(ξi ::sTfpe,
                     ξ1 ::sTfpe,
                     ξ2 ::sTfpe,
                     α  ::Float64,
                     σa ::Float64,
                     σk ::Float64,
                     ll ::Float64,
                     dα ::Float64,
                     sσa::Float64,
                     sσk::Float64)

Perform a punkeek quartet update.
"""
function _update_quartet!(ξi ::sTfpe,
                          ξ1 ::sTfpe,
                          ξ2 ::sTfpe,
                          α  ::Float64,
                          σa ::Float64,
                          σk ::Float64,
                          ll ::Float64,
                          dα ::Float64,
                          sσa::Float64,
                          sσk::Float64)

  σa2, σk2 = σa^2, σk^2
  xa, xic, x1, x2 = xi(ξi), xf(ξi), xf(ξ1), xf(ξ2)
  ei, e1, e2 = e(ξi), e(ξ1), e(ξ2)

  ### probabilities
  ## which is cladogenetic
  # d1 cladogenetic
  pk1 = llik_cpe_triad(xa + α*ei, x2 - α*e2, x1 - α*e1, ei, e2, e1, σa2, σk2)
  # d2 cladogenetic
  pk2 = llik_cpe_triad(xa + α*ei, x1 - α*e1, x2 - α*e2, ei, e1, e2, σa2, σk2)
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
  xip = trioprop(xa + α*ei, xadp - α*eadp,  xkdp - α*ekdp, 
                    ei*σa2,      eadp*σa2, ekdp*σa2 + σk2)
  xkp = duoprop(xip, xkdp - α*ekdp, σk2, ekdp*σa2)

  # update likelihood and gibbs quanta
  xkc, xadc, xkdc, eadc, ekdc = 
    sh(ξi) ? (xi(ξ1), x2, x1, e2, e1) : (xi(ξ2), x1, x2, e1, e2)

  ll += llik_cpe_quartet(xa + α*ei, xip, xkp, xadp - α*eadp, xkdp - α*ekdp, 
          ei, eadp, ekdp, σa2, σk2) -
        llik_cpe_quartet(xa + α*ei, xic, xkc, xadc - α*eadc, xkdc - α*ekdc, 
          ei, eadc, ekdc, σa2, σk2)

  dα  += xip - xic                    + 
         (xadp - xip) + (xkdp - xkp)  - 
         (xadc - xic) - (xkdc - xkc)
  sσa += 0.5*(((xip - xa - α*ei)^2 - (xic - xa - α*ei)^2)/ei +
              (xadp - xip - α*eadp)^2/eadp                   -
              (xadc - xic - α*eadc)^2/eadc                   +
              (xkdp - xkp - α*ekdp)^2/ekdp                   -
              (xkdc - xkc - α*ekdc)^2/ekdc)
  sσk += 0.5*((xkp - xip)^2 - (xkc - xic)^2)

  setsh!(ξi, shp)
  setxf!(ξi, xip)

  ξk, ξa = shp ? (ξ1, ξ2) : (ξ2, ξ1)
  setxi!(ξk, xkp) 
  setxi!(ξa, xip)

  return ll, dα, sσa, sσk
end




