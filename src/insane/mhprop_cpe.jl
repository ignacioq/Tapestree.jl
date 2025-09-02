#=

punkeek proposals

Ignacio Quintero Mächler

t(-_-t)

Created 25 11 2024
=#



"""
     _stem_update!(ξi ::T,
                   σa ::Float64,
                   ll ::Float64,
                   sσa::Float64) where {T <: Tpe}

Perform punkeek stem update.
"""
function _stem_update!(ξi ::T,
                       σa ::Float64,
                       ll ::Float64,
                       sσa::Float64) where {T <: Tpe}
  ei  = e(ξi)
  xc  = xi(ξi)
  xfi = xf(ξi)
  eσ  = σa*sqrt(ei)

  # sample mrca
  xp   = rnorm(xfi, eσ)
  ll  += lrdnorm_bm_x(xp, xc, xfi, eσ)
  sσa += (xfi - xp)^2/ei - (xfi - xc)^2/ei
  setxi!(ξi, xp)

  return ll, sσa
end




"""
    _crown_update!(ξi ::T,
                   ξ1 ::T,
                   ξ2 ::T,
                   σa ::Float64,
                   σk ::Float64,
                   ll ::Float64,
                   sσa::Float64,
                   sσk::Float64)  where {T <: Tpe}

Perform punkeek crown update.
"""
function _crown_update!(ξi ::T,
                        ξ1 ::T,
                        ξ2 ::T,
                        σa ::Float64,
                        σk ::Float64,
                        ll ::Float64,
                        sσa::Float64,
                        sσk::Float64) where {T <: Tpe}

  σa2, σk2 = σa^2, σk^2
  xic, x1, x2 = xf(ξi), xf(ξ1), xf(ξ2)
  e1, e2 = e(ξ1), e(ξ2)

  ## proposals
  # if d1 is cladogenetic
  xip1 = duoprop(x1, x2, σk2 + e1*σa2, e2*σa2)
  xkp1 = duoprop(xip1, x1, σk2, e1*σa2)

  # if d2 is cladogenetic
  xip2 = duoprop(x1, x2, e1*σa2, σk2 + e2*σa2)
  xkp2 = duoprop(xip2, x2, σk2, e2*σa2)

  ## probabilities
  # decide which with relative prob of K1 and K2
  pk1 = llik_trio(xip1, xkp1, x2, x1, e2, e1, σa2, σk2)
  pk2 = llik_trio(xip2, xkp2, x1, x2, e1, e2, σa2, σk2)
  o12 = exp(pk1 - pk2)  # odds
  p1  = o12/(1.0 + o12) # probability

  ## update trackers
  # current likelihood and ss
  if sh(ξi)
    xk   = xi(ξ1)
    ll  -= llik_trio(xic, xk, x2, x1, e2, e1, σa2, σk2)
    sσa -= (x2 - xic)^2/e2 + (x1 - xk)^2/e1
    sσk -= (xk - xic)^2
  else
    xk   = xi(ξ2)
    ll  -= llik_trio(xic, xk, x1, x2, e1, e2, σa2, σk2)
    sσa -= (x1 - xic)^2/e1 + (x2 - xk)^2/e2
    sσk -= (xk - xic)^2
  end

  # if x1 cladogenetic
  if rand() < p1
    ll  += pk1
    sσa += (x2 - xip1)^2/e2 + (x1 - xkp1)^2/e1
    sσk += (xkp1 - xip1)^2
    setsh!(ξi, true)
    setxi!(ξi, xip1)
    setxf!(ξi, xip1)
    setxi!(ξ1, xkp1)
    setxi!(ξ2, xip1)
  # if x2 cladogenetic
  else
    ll  += pk2
    sσa += (x1 - xip2)^2/e1 + (x2 - xkp2)^2/e2
    sσk += (xkp2 - xip2)^2
    setsh!(ξi, false)
    setxi!(ξi, xip2)
    setxf!(ξi, xip2)
    setxi!(ξ1, xip2)
    setxi!(ξ2, xkp2)
  end

  return ll, sσa, sσk
end



"""
    _update_node_x!(tree::T,
                    σa  ::Float64,
                    σk  ::Float64,
                    ll  ::Float64,
                    sσa ::Float64,
                    sσk ::Float64) where {T <: Tpe}

Perform punkeek internal node updates.
"""
function _update_node_x!(tree::T,
                         σa  ::Float64,
                         σk  ::Float64,
                         ll  ::Float64,
                         sσa ::Float64,
                         sσk ::Float64) where {T <: Tpe}

  if def1(tree)
    if def2(tree)
      ll, sσa, sσk = _update_quartet_x!(tree, σa, σk, ll, sσa, sσk)

      ll, sσa, sσk = _update_node_x!(tree.d1, σa, σk, ll, sσa, sσk)
      ll, sσa, sσk = _update_node_x!(tree.d2, σa, σk, ll, sσa, sσk)
    else
      ll, sσa, sσk = _update_node_x!(tree.d1, σa, σk, ll, sσa, sσk)
    end
  elseif !isfix(tree)
    ll, sσa = _update_tip_x!(tree, σa, ll, sσa)
  end

  return ll, sσa, sσk
end



"""
    _update_leaf_x!(tree::T,
                    xavg::Float64,
                    xstd::Float64,
                    σa  ::Float64,
                    σk  ::Float64,
                    ll  ::Float64,
                    sσa ::Float64,
                    sσk ::Float64) where {T <: Tpe}

Perform punkeek **fixed** leaf (terminal reconstructed edge) updates.
"""
function _update_leaf_x!(tree::T,
                         xavg::Float64,
                         xstd::Float64,
                         σa  ::Float64,
                         σk  ::Float64,
                         ll  ::Float64,
                         sσa ::Float64,
                         sσk ::Float64) where {T <: Tpe}

  if def1(tree)
    if def2(tree)
      ll, sσa, sσk = _update_quartet_x!(tree, σa, σk, ll, sσa, sσk)

      ll, sσa, sσk = _update_leaf_x!(tree.d1, xavg, xstd, σa, σk, ll, sσa, sσk)
      ll, sσa, sσk = _update_leaf_x!(tree.d2, xavg, xstd, σa, σk, ll, sσa, sσk)
    else
      ll, sσa, sσk = _update_leaf_x!(tree.d1, xavg, xstd, σa, σk, ll, sσa, sσk)
    end
  elseif isfix(tree)
    if !iszero(xstd)
      ll, sσa = _update_tip_x!(tree, xavg, xstd, σa, ll, sσa)
    end
  else
    ll, sσa = _update_tip_x!(tree, σa, ll, sσa)
  end

  return ll, sσa, sσk
end




"""
    _update_leaf_x!(tree::T,
                    σa  ::Float64,
                    σk  ::Float64,
                    ll  ::Float64,
                    sσa ::Float64,
                    sσk ::Float64) where {T <: Tpe}

Perform punkeek **unfixed** leaf (terminal reconstructed edge) updates.
"""
function _update_leaf_x!(tree::T,
                         σa  ::Float64,
                         σk  ::Float64,
                         ll  ::Float64,
                         sσa ::Float64,
                         sσk ::Float64) where {T <: Tpe}

  if def1(tree)
    if def2(tree)
      ll, sσa, sσk = _update_quartet_x!(tree, σa, σk, ll, sσa, sσk)

      ll, sσa, sσk = _update_leaf_x!(tree.d1, σa, σk, ll, sσa, sσk)
      ll, sσa, sσk = _update_leaf_x!(tree.d2, σa, σk, ll, sσa, sσk)
    else
      ll, sσa, sσk = _update_leaf_x!(tree.d1, σa, σk, ll, sσa, sσk)
    end
  else
    ll, sσa = _update_tip_x!(tree, σa, ll, sσa)
  end

  return ll, sσa, sσk
end



"""
    _update_tip_x!(tree::T,
                   σa  ::Float64, 
                   ll  ::Float64, 
                   sσa ::Float64)  where {T <: Tpe}

Perform punkeek **unfixed** tip updates.
"""
function _update_tip_x!(tree::T,
                        σa  ::Float64, 
                        ll  ::Float64, 
                        sσa ::Float64) where {T <: Tpe}

  xa, xic = xi(tree), xf(tree)
  ei = e(tree)

  # proposal
  xip = rnorm(xa, sqrt(ei)*σa)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa, ei*σa^2)
  sσa += ((xip - xa)^2 - (xic - xa)^2)/ei
  setxf!(tree, xip)

  return ll, sσa
end




"""
    _update_tip_x!(tree::T,
                    σa  ::Float64, 
                    σk  ::Float64, 
                    ll  ::Float64, 
                    sσa ::Float64, 
                    sσk ::Float64)  where {T <: Tpe}

Perform punkeek **fixed** tip updates.
"""
function _update_tip_x!(tree::T,
                        xavg::Float64,
                        xstd::Float64,
                        σa  ::Float64, 
                        ll  ::Float64, 
                        sσa ::Float64) where {T <: Tpe}

  xa, xic = xi(tree), xf(tree)
  ei = e(tree)

  xip = duoprop(xavg, xic, xstd^2, ei*σa^2)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa, ei*σa^2)
  sσa += ((xip - xa)^2 - (xic - xa)^2)/ei
  setxf!(tree, xip)

  return ll, sσa
end




"""
    _update_duo_x!(ξi  ::T,
                   ξ1  ::T,
                   σa  ::Float64,
                   ll  ::Float64,
                   sσa ::Float64) where {T <: Tpe}

Perform punkeek for **unfixed** node.
"""
function _update_duo_x!(ξi  ::T,
                        ξ1  ::T,
                        σa  ::Float64,
                        ll  ::Float64,
                        sσa ::Float64) where {T <: Tpe}

  σa2 = σa^2
  xa, xic, x1 = xi(ξi), xf(ξi), xf(ξ1)
  ei, e1 = e(ξi), e(ξ1)

  # sample
  xip = duoprop(xa, x1, ei*σa2, e1*σa2)

  ## update trackers
  ll  += llrdnorm_x(xip, xic, xa, ei*σa2) + llrdnorm_μ(x1, xip, xic, e1*σa2)
  sσa += ((xip - xa)^2 - (xic - xa)^2)/ei + ((x1 - xip)^2 - (x1 - xic)^2)/e1
  setxf!(ξi, xip)
  setxi!(ξ1, xip)

  return ll, sσa
end




"""
    _update_quartet_x!(ξi ::T,
                       σa ::Float64,
                       σk ::Float64,
                       ll ::Float64,
                       sσa::Float64,
                       sσk::Float64)

Make a punkeek quartet proposal.
"""
function _update_quartet_x!(ξi ::T,
                            σa ::Float64,
                            σk ::Float64,
                            ll ::Float64,
                            sσa::Float64,
                            sσk::Float64) where {T <: Tpe}

  ξ1 = ξi.d1
  ξ2 = ξi.d2
  σa2, σk2 = σa^2, σk^2
  xa, xic, x1, x2 = xi(ξi), xf(ξi), xf(ξ1), xf(ξ2)
  ei, e1, e2 = e(ξi), e(ξ1), e(ξ2)

  ## proposals
  # if d1 is cladogenetic
  xip1 = trioprop(xa, x1, x2, ei*σa2, σk2 + e1*σa2, e2*σa2)
  xkp1 = duoprop(xip1, x1, σk2, e1*σa2)

  # if d2 is cladogenetic
  xip2 = trioprop(xa, x1, x2, ei*σa2, e1*σa2, σk2 + e2*σa2)
  xkp2 = duoprop(xip2, x2, σk2, e2*σa2)

  ## probabilities
  # decide which with relative prob of K1 and K2
  pk1 = llik_quartet(xa, xip1, xkp1, x2, x1, ei, e2, e1, σa2, σk2)
  pk2 = llik_quartet(xa, xip2, xkp2, x1, x2, ei, e1, e2, σa2, σk2)
  o12 = exp(pk1 - pk2)   # odds
  p1  = o12/(1.0 + o12) # probability

  ## update trackers
  # current likelihood and 
  if sh(ξi)
    xk   = xi(ξ1)
    ll  -= llik_quartet(xa, xic, xk, x2, x1, ei, e2, e1, σa2, σk2)
    sσa -= (xic - xa)^2/ei + (x2 - xic)^2/e2 + (x1 - xk)^2/e1
    sσk -= (xk - xic)^2
  else
    xk   = xi(ξ2)
    ll  -= llik_quartet(xa, xic, xk, x1, x2, ei, e1, e2, σa2, σk2)
    sσa -= (xic - xa)^2/ei + (x1 - xic)^2/e1 + (x2 - xk)^2/e2
    sσk -= (xk - xic)^2
  end

  # if x1 cladogenetic
  if rand() < p1
    ll  += pk1
    sσa += (xip1 - xa)^2/ei + (x2 - xip1)^2/e2 + (x1 - xkp1)^2/e1
    sσk += (xkp1 - xip1)^2
    setsh!(ξi, true)
    setxf!(ξi, xip1)
    setxi!(ξ1, xkp1)
    setxi!(ξ2, xip1)
  # if x2 cladogenetic
  else
    ll  += pk2
    sσa += (xip2 - xa)^2/ei + (x1 - xip2)^2/e1 + (x2 - xkp2)^2/e2
    sσk += (xkp2 - xip2)^2
    setsh!(ξi, false)
    setxf!(ξi, xip2)
    setxi!(ξ1, xip2)
    setxi!(ξ2, xkp2)
  end

  return ll, sσa, sσk
end




"""
    _update_node_x!(ξi ::T,
                    ξ1 ::T,
                    ξ2 ::T,
                    σa ::Float64,
                    σk ::Float64,
                    ll ::Float64,
                    sσa::Float64,
                    sσk::Float64) where {T <: Tpe}

Perform a punkeek quartet update.
"""
function _update_node_x!(ξi ::T,
                         ξ1 ::T,
                         ξ2 ::T,
                         σa ::Float64,
                         σk ::Float64,
                         ll ::Float64,
                         sσa::Float64,
                         sσk::Float64) where {T <: Tpe}

  σa2, σk2 = σa^2, σk^2
  xa, xic, x1, x2 = xi(ξi), xf(ξi), xf(ξ1), xf(ξ2)
  ei, e1, e2 = e(ξi), e(ξ1), e(ξ2)

  ## proposals
  # if d1 is cladogenetic
  xip1 = trioprop(xa, x1, x2, ei*σa2, σk2 + e1*σa2, e2*σa2)
  xkp1 = duoprop(xip1, x1, σk2, e1*σa2)

  # if d2 is cladogenetic
  xip2 = trioprop(xa, x1, x2, ei*σa2, e1*σa2, σk2 + e2*σa2)
  xkp2 = duoprop(xip2, x2, σk2, e2*σa2)

  ## probabilities
  # decide which with relative prob of K1 and K2
  pk1 = llik_quartet(xa, xip1, xkp1, x2, x1, ei, e2, e1, σa2, σk2)
  pk2 = llik_quartet(xa, xip2, xkp2, x1, x2, ei, e1, e2, σa2, σk2)
  o12 = exp(pk1 - pk2)   # odds
  p1  = o12/(1.0 + o12) # probability

  ## update trackers
  # current likelihood and 
  if sh(ξi)
    xk   = xi(ξ1)
    ll  -= llik_quartet(xa, xic, xk, x2, x1, ei, e2, e1, σa2, σk2)
    sσa -= (xic - xa)^2/ei + (x2 - xic)^2/e2 + (x1 - xk)^2/e1
    sσk -= (xk - xic)^2
  else
    xk   = xi(ξ2)
    ll  -= llik_quartet(xa, xic, xk, x1, x2, ei, e1, e2, σa2, σk2)
    sσa -= (xic - xa)^2/ei + (x1 - xic)^2/e1 + (x2 - xk)^2/e2
    sσk -= (xk - xic)^2
  end

  # if x1 cladogenetic
  if rand() < p1
    ll  += pk1
    sσa += (xip1 - xa)^2/ei + (x2 - xip1)^2/e2 + (x1 - xkp1)^2/e1
    sσk += (xkp1 - xip1)^2
    setsh!(ξi, true)
    setxf!(ξi, xip1)
    setxi!(ξ1, xkp1)
    setxi!(ξ2, xip1)
  # if x2 cladogenetic
  else
    ll  += pk2
    sσa += (xip2 - xa)^2/ei + (x1 - xip2)^2/e1 + (x2 - xkp2)^2/e2
    sσk += (xkp2 - xip2)^2
    setsh!(ξi, false)
    setxf!(ξi, xip2)
    setxi!(ξ1, xip2)
    setxi!(ξ2, xkp2)
  end

  return ll, sσa, sσk
end

















