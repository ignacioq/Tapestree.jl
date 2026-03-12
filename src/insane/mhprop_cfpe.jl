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
  sσa += (xfi - xp - α*ei)^2/ei - (xfi - xc - α*ei)^2/ei
  setxi!(ξi, xp)

  return ll, dα, sσa
end




"""
    _fstem_update!(ξi ::sTfpe,
                   ξ1 ::sTfpe,
                   α  ::Float64,
                   σa ::Float64,
                   ll ::Float64,
                   dα ::Float64,
                   sσa::Float64)

Perform punkeek stem fossil update.
"""
function _fstem_update!(ξi ::sTfpe,
                        ξ1 ::sTfpe,
                        α  ::Float64,
                        σa ::Float64,
                        ll ::Float64,
                        dα ::Float64,
                        sσa::Float64)
  ei  = e(ξ1)
  xc  = xi(ξi)
  xfi = xf(ξ1)
  eσ  = σa*sqrt(ei)

  # sample mrca
  xp   = rnorm(xfi - α*ei, eσ)
  ll  += lrdnorm_bm_x(xp, xc - α*ei, xfi, eσ)
  dα  += xc - xp
  sσa += (xfi - xp - α*ei)^2/ei - (xfi - xc - α*ei)^2/ei
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
  # if x1 cladogenetic
  shp, xadp, xkdp, eadp, ekdp = 
    rand(Bool) ? (true, x2, x1, e2, e1) : (false, x1, x2, e1, e2)

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
         (xadc - xic) + (xkdc - xkc)
  sσa += (xadp - xip - α*eadp)^2/eadp + (xkdp - xkp - α*ekdp)^2/ekdp -
         (xadc - xic - α*eadc)^2/eadc - (xkdc - xkc - α*ekdc)^2/ekdc
  sσk += (xkp - xip)^2 - (xkc - xic)^2

  setsh!(ξi, shp)
  setxi!(ξi, xip)
  setxf!(ξi, xip)
  ξk, ξa = shp ? (ξ1, ξ2) : (ξ2, ξ1)
  setxi!(ξk, xkp) 
  setxi!(ξa, xip)

  return ll, dα, sσa, sσk
end




"""
    _update_node_x!(tree::sTfpe,
                    α   ::Float64,
                    σa  ::Float64,
                    σk  ::Float64,
                    ll  ::Float64,
                    dα  ::Float64,
                    sσa ::Float64,
                    sσk ::Float64)

Perform punkeek internal node updates.
"""
function _update_node_x!(tree::sTfpe,
                         α   ::Float64,
                         σa  ::Float64,
                         σk  ::Float64,
                         ll  ::Float64,
                         dα  ::Float64,
                         sσa ::Float64,
                         sσk ::Float64)

  if def1(tree)
    if def2(tree)
      ll, dα, sσa, sσk = _update_quartet_x!(tree, α, σa, σk, ll, dα, sσa, sσk)

      ll, dα, sσa, sσk = _update_node_x!(tree.d1, α, σa, σk, ll, dα, sσa, sσk)
      ll, dα, sσa, sσk = _update_node_x!(tree.d2, α, σa, σk, ll, dα, sσa, sσk)
    else
      ll, dα, sσa, sσk = _update_node_x!(tree.d1, α, σa, σk, ll, dα, sσa, sσk)
    end
  elseif !isfix(tree)
    ll, dα, sσa = _update_tip_x!(tree, α, σa, ll, dα, sσa)
  end

  return ll, dα, sσa, sσk
end





"""
    _update_leaf_x!(tree::sTfpe,
                    xavg::Float64,
                    xstd::Float64,
                    α   ::Float64,
                    σa  ::Float64,
                    σk  ::Float64,
                    ll  ::Float64,
                    dα  ::Float64,
                    sσa ::Float64,
                    sσk ::Float64)

Perform punkeek **fixed** leaf (terminal reconstructed edge) updates.
"""
function _update_leaf_x!(tree::sTfpe,
                         xavg::Float64,
                         xstd::Float64,
                         α   ::Float64,
                         σa  ::Float64,
                         σk  ::Float64,
                         ll  ::Float64,
                         dα  ::Float64,
                         sσa ::Float64,
                         sσk ::Float64)

  if def1(tree)
    if def2(tree)
      ll, dα, sσa, sσk = _update_quartet_x!(tree, α, σa, σk, ll, dα, sσa, sσk)

      ll, dα, sσa, sσk = 
        _update_leaf_x!(tree.d1, xavg, xstd, α, σa, σk, ll, dα, sσa, sσk)
      ll, dα, sσa, sσk = 
        _update_leaf_x!(tree.d2, xavg, xstd, α, σa, σk, ll, dα, sσa, sσk)
    else
      ll, dα, sσa, sσk = 
        _update_leaf_x!(tree.d1, xavg, xstd, α, σa, σk, ll, dα, sσa, sσk)
    end
  elseif isfix(tree)
    if !iszero(xstd)
      ll, dα, sσa = _update_tip_x!(tree, xavg, xstd, α, σa, ll, dα, sσa)
    end
  else
    ll, dα, sσa = _update_tip_x!(tree, α, σa, ll, dα, sσa)
  end

  return ll, dα, sσa, sσk
end




"""
    _update_leaf_x!(tree::sTfpe,
                    α   ::Float64,
                    σa  ::Float64,
                    σk  ::Float64,
                    ll  ::Float64,
                    dα  ::Float64,
                    sσa ::Float64,
                    sσk ::Float64)

Perform punkeek **unfixed** leaf (terminal reconstructed edge) updates.
"""
function _update_leaf_x!(tree::sTfpe,
                         α   ::Float64,
                         σa  ::Float64,
                         σk  ::Float64,
                         ll  ::Float64,
                         dα  ::Float64,
                         sσa ::Float64,
                         sσk ::Float64)

  if def1(tree)
    if def2(tree)
      ll, dα, sσa, sσk = _update_quartet_x!(tree, α, σa, σk, ll, dα, sσa, sσk)

      ll, dα, sσa, sσk = _update_leaf_x!(tree.d1, α, σa, σk, ll, dα, sσa, sσk)
      ll, dα, sσa, sσk = _update_leaf_x!(tree.d2, α, σa, σk, ll, dα, sσa, sσk)
    else
      ll, dα, sσa, sσk = _update_leaf_x!(tree.d1, α, σa, σk, ll, dα, sσa, sσk)
    end
  else
    ll, dα, sσa = _update_tip_x!(tree, α, σa, ll, dα, sσa)
  end

  return ll, dα, sσa, sσk
end



"""
    _update_tip_x!(tree::sTfpe,
                   α   ::Float64,
                   σa  ::Float64, 
                   ll  ::Float64, 
                   dα  ::Float64,
                   sσa ::Float64)

Perform punkeek **unfixed** tip updates.
"""
function _update_tip_x!(tree::sTfpe,
                        α   ::Float64,
                        σa  ::Float64, 
                        ll  ::Float64, 
                        dα  ::Float64,
                        sσa ::Float64)

  xa, xic = xi(tree), xf(tree)
  ei = e(tree)

  # proposal
  xip = rnorm(xa + α*ei, sqrt(ei)*σa)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa + α*ei, ei*σa^2)
  dα  += xip - xic
  sσa += ((xip - xa - α*ei)^2 - (xic - xa - α*ei)^2)/ei
  setxf!(tree, xip)

  return ll, dα, sσa
end




"""
    _update_tip_x!(tree::sTfpe,
                   xavg::Float64,
                   xstd::Float64,
                   α   ::Float64,
                   σa  ::Float64, 
                   ll  ::Float64, 
                   dα  ::Float64,
                   sσa ::Float64)

Perform punkeek **fixed** tip updates.
"""
function _update_tip_x!(tree::sTfpe,
                        xavg::Float64,
                        xstd::Float64,
                        α   ::Float64,
                        σa  ::Float64, 
                        ll  ::Float64, 
                        dα  ::Float64,
                        sσa ::Float64)

  xa, xic = xi(tree), xf(tree)
  ei = e(tree)

  xip = duoprop(xavg, xic + α*ei, xstd^2, ei*σa^2)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa + α*ei, ei*σa^2)
  dα  += xip - xic
  sσa += ((xip - xa - α*ei)^2 - (xic - xa - α*ei)^2)/ei
  setxf!(tree, xip)

  return ll, dα, sσa
end




"""
    _update_duo_x!(ξi  ::sTfpe,
                   ξ1  ::sTfpe,
                   α   ::Float64,
                   σa  ::Float64,
                   ll  ::Float64,
                   sσa ::Float64)

Perform punkeek for **unfixed** node.
"""
function _update_duo_x!(ξi  ::sTfpe,
                        ξ1  ::sTfpe,
                        α   ::Float64,
                        σa  ::Float64,
                        ll  ::Float64,
                        sσa ::Float64)

  σa2 = σa^2
  xa, xic, x1 = xi(ξi), xf(ξi), xf(ξ1)
  ei, e1 = e(ξi), e(ξ1)

  # sample
  xip = duoprop(xa + α*ei, x1 - α*e1, ei*σa2, e1*σa2)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa + α*ei, ei*σa2) + 
         llrdnorm_μ(x1 - α*e1, xip, xic, e1*σa2)
  sσa += ((xip - xa - α*ei)^2 - (xic - xa - α*ei)^2)/ei + 
         ((x1 - xip - α*e1)^2 - (x1 - xic - α*e1)^2)/e1
  setxf!(ξi, xip)
  setxi!(ξ1, xip)

  return ll, sσa
end




"""
    _update_duo_x!(ξi  ::sTfpe,
                   ξ1  ::sTfpe,
                   xavg::Float64,
                   xstd::Float64,
                   α   ::Float64,
                   σa  ::Float64,
                   ll  ::Float64,
                   sσa ::Float64)

Perform punkeek for **fixed** node.
"""
function _update_duo_x!(ξi  ::sTfpe,
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
  xip = trioprop(xavg, xa + α*ei, x1 - α*e1, xstd^2, ei*σa2, e1*σa2)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa + α*ei, ei*σa2) + 
         llrdnorm_μ(x1 - α*e1, xip, xic, e1*σa2)
  sσa += ((xip - xa - α*ei)^2 - (xic - xa - α*ei)^2)/ei + 
         ((x1 - xip - α*e1)^2 - (x1 - xic - α*e1)^2)/e1
  setxf!(ξi, xip)
  setxi!(ξ1, xip)

  return ll, sσa
end




"""
    _update_quartet_x!(ξi ::sTfpe,
                       α  ::Float64,
                       σa ::Float64,
                       σk ::Float64,
                       ll ::Float64,
                       dα ::Float64,
                       sσa::Float64,
                       sσk::Float64)

Make a punkeek quartet proposal.
"""
function _update_quartet_x!(ξi ::sTfpe,
                            α  ::Float64,
                            σa ::Float64,
                            σk ::Float64,
                            ll ::Float64,
                            dα ::Float64,
                            sσa::Float64,
                            sσk::Float64)

  ll, dα, sσa, sσk = 
    _update_node_x!(ξi, ξi.d1, ξi.d2, α, σa, σk, ll, dα, sσa, sσk)

  return ll, dα, sσa, sσk
end




"""
    _update_node_x!(ξi ::sTfpe,
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
function _update_node_x!(ξi ::sTfpe,
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
  o12 = exp(pk1 - pk2)  # odds
  p1  = o12/(1.0 + o12) # probability

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
  sσa += (xip - xa - α*ei)^2/ei       -
         (xic - xa - α*ei)^2/ei       +
         (xadp - xip - α*eadp)^2/eadp - 
         (xadc - xic - α*eadc)^2/eadc +
         (xkdp - xkp - α*ekdp)^2/ekdp -
         (xkdc - xkc - α*ekdc)^2/ekdc
  sσk += (xkp - xip)^2 - (xkc - xic)^2

  setsh!(ξi, shp)
  setxf!(ξi, xip)

  ξk, ξa = shp ? (ξ1, ξ2) : (ξ2, ξ1)
  setxi!(ξk, xkp) 
  setxi!(ξa, xip)

  return ll, dα, sσa, sσk
end




