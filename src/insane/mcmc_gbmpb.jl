#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#





"""
    lλupdate!(tra ::iTgbmpb,
              trap::iTgbmpb,
              llc ::Float64, 
              prc ::Float64,
              σ²λ ::Float64, 
              srδt::Float64, 
              λa_prior::Tuple{Float64,Float64},
              dri  ::BitArray{1},
              ldr  ::Int64,
              ix   ::Int64)

Return subtree given BitVector `dri`.
"""
function lλupdate!(tra ::iTgbmpb,
                   trap::iTgbmpb,
                   llc ::Float64, 
                   prc ::Float64,
                   σ²λ ::Float64, 
                   srδt::Float64, 
                   λa_prior::Tuple{Float64,Float64},
                   dri  ::BitArray{1},
                   ldr  ::Int64,
                   ix   ::Int64)

  if ix == ldr 
    # if root
    if ldr == 0
      llc, prc = triad_lλupdate_root!(tra::iTgbmpb, trap::iTgbmpb, 
                   llc, prc, σ²λ, srδt, λa_prior)
    else
      if ter[1]
        if ter[2]
          # if both are terminal
          llc = triad_lλupdate_noded12!(tra::iTgbmpb, trap::iTgbmpb, 
                        llc, σ²λ, srδt)
        else
          # if d1 is terminal
          llc = triad_lλupdate_noded1!(tra::iTgbmpb, trap::iTgbmpb, 
                       llc, σ²λ, srδt)
        end
      elseif ter[2]
        # if d2 is terminal
        llc = triad_lλupdate_noded2!(tra::iTgbmpb, trap::iTgbmpb, 
                     llc, σ²λ, srδt)
      else
        # if no terminal branches involved
        llc = triad_lλupdate_node!(tra::iTgbmpb, trap::iTgbmpb, 
                     llc, σ²λ, srδt)
      end
    end

  elseif ix < ldr
    ix += 1
    if dri[ix]
      llc, prc = 
        lλupdate!(tra.d1::iTgbmpb, trap.d1::iTgbmpb, 
          llc, prc, σ²λ, srδt, λa_prior, dri, ldr, ix)
    else
      llc, prc = 
        lλupdate!(tra.d2::iTgbmpb, trap.d2::iTgbmpb, 
          llc, prc, σ²λ, srδt, λa_prior, dri, ldr, ix)
    end
  end

  return llc, prc
end






"""
    triad_lλupdate_noded12!(tra     ::iTgbmpb, 
                            trap    ::iTgbmpb,
                            llc     ::Float64,
                            prc     ::Float64,
                            σ²λ     ::Float64, 
                            srδt    ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
both daughters are terminal.
"""
function triad_lλupdate_noded12!(tra     ::iTgbmpb, 
                                 trap    ::iTgbmpb,
                                 llc     ::Float64,
                                 σ²λ     ::Float64, 
                                 srδt    ::Float64)

  # get reference vectors
  λprv_p = lλ(trap)
  λd1v_p = lλ(trap.d1)
  λd2v_p = lλ(trap.d2)
  λprv_c = lλ(tra)
  λd1v_c = lλ(tra.d1)
  λd2v_c = lλ(tra.d2)
  λpr = λprv_c[1]

  # time vectors
  tprv = ts(tra)
  td1v = ts(tra.d1)
  td2v = ts(tra.d2)

  # fill with Brownian motion
  bm!(λprv_p, λpr, tprv, σ²λ, srδt)
  lλp = λprv_p[end]
  bm!(λd1v_p, lλp, td1v, σ²λ, srδt)
  bm!(λd2v_p, lλp, td2v, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tprv, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(td1v, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(td2v, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tprv, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(td1v, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(td2v, λd2v_c, σ²λ, δt)

  # speciation likelihood ratio
  λlr = 2.0*lλp - 2.0*λd1v_c[1]

  # likelihood ratio
  llr = llbmprp + llpbprp + llbmd1p + llpbd1p + llbmd2p + llpbd2p -
        llbmprc - llpbprc - llbmd1c - llpbd1c - llbmd2c - llpbd2c +
        λlr
  # acceptance ratio
  acr = llpbprp + llpbd1p + llpbd2p -
        llpbprc - llpbd1c - llpbd2c + λlr

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_noded1!(tra     ::iTgbmpb, 
                           trap    ::iTgbmpb,
                           llc     ::Float64,
                           σ²λ     ::Float64, 
                           srδt    ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 1 is terminal.
"""
function triad_lλupdate_noded1!(tra     ::iTgbmpb, 
                                trap    ::iTgbmpb,
                                llc     ::Float64,
                                σ²λ     ::Float64, 
                                srδt    ::Float64)

  # get reference vectors
  λprv_p = lλ(trap)
  λd1v_p = lλ(trap.d1)
  λd2v_p = lλ(trap.d2)
  λprv_c = lλ(tra)
  λd1v_c = lλ(tra.d1)
  λd2v_c = lλ(tra.d2)
  λpr = λprv_c[1]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(tra)
  td1v = ts(tra.d1)
  td2v = ts(tra.d2)

  # node proposal
  lλp = duoprop(λpr, λd2, pe(tra), pe(tra.d2), σ²λ)

  bb!(λprv_p, λpr, lλp, tprv, σ²λ, srδt)
  bm!(λd1v_p, lλp, td1v, σ²λ, srδt)
  bb!(λd2v_p, lλp, λd2, td2v, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tprv, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(td1v, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(td2v, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tprv, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(td1v, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(td2v, λd2v_c, σ²λ, δt)

  # speciation likelihood ratio
  λlr = 2.0*lλp - 2.0*λd1v_c[1]

  # likelihood ratio
  llr = llbmprp + llpbprp + llbmd1p + llpbd1p + llbmd2p + llpbd2p -
        llbmprc - llpbprc - llbmd1c - llpbd1c - llbmd2c - llpbd2c +
        λlr

  # acceptance ratio
  acr = llpbprp + llpbd1p + llpbd2p -
        llpbprc - llpbd1c - llpbd2c + λlr

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_noded2!(tra     ::iTgbmpb, 
                           trap    ::iTgbmpb,
                           llc     ::Float64,
                           prc     ::Float64,
                           σ²λ     ::Float64, 
                           srδt    ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 2 is terminal.
"""
function triad_lλupdate_noded2!(tra     ::iTgbmpb, 
                                trap    ::iTgbmpb,
                                llc     ::Float64,
                                σ²λ     ::Float64, 
                                srδt    ::Float64)

  # get reference vectors
  λprv_p = lλ(trap)
  λd1v_p = lλ(trap.d1)
  λd2v_p = lλ(trap.d2)
  λprv_c = lλ(tra)
  λd1v_c = lλ(tra.d1)
  λd2v_c = lλ(tra.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]

  # time vectors
  tprv = ts(tra)
  td1v = ts(tra.d1)
  td2v = ts(tra.d2)

  # node proposal
  lλp = duoprop(λpr, λd1, pe(tra), pe(tra.d1), σ²λ)

  bb!(λprv_p, λpr, lλp, tprv, σ²λ, srδt)
  bb!(λd1v_p, lλp, λd1, td1v, σ²λ, srδt)
  bm!(λd2v_p, lλp, td2v, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tprv, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(td1v, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(td2v, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tprv, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(td1v, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(td2v, λd2v_c, σ²λ, δt)

  # speciation likelihood ratio
  λlr = 2.0*lλp - 2.0*λd1v_c[1]

  # likelihood ratio
  llr = llbmprp + llpbprp + llbmd1p + llpbd1p + llbmd2p + llpbd2p -
        llbmprc - llpbprc - llbmd1c - llpbd1c - llbmd2c - llpbd2c +
        λlr

  # acceptance ratio
  acr = llpbprp + llpbd1p + llpbd2p -
        llpbprc - llpbd1c - llpbd2c + λlr

  if -randexp() < acr 
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_node!(tra     ::iTgbmpb, 
                          trap    ::iTgbmpb,
                          llc     ::Float64,
                          σ²λ     ::Float64, 
                          srδt    ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
no daughters are terminal.
"""
function triad_lλupdate_node!(tra     ::iTgbmpb, 
                              trap    ::iTgbmpb,
                              llc     ::Float64,
                              σ²λ     ::Float64, 
                              srδt    ::Float64)

  # get reference vectors
  λprv_p = lλ(trap)
  λd1v_p = lλ(trap.d1)
  λd2v_p = lλ(trap.d2)
  λprv_c = lλ(tra)
  λd1v_c = lλ(tra.d1)
  λd2v_c = lλ(tra.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(tra)
  td1v = ts(tra.d1)
  td2v = ts(tra.d2)

  # node proposal
  lλp = trioprop(λpr, λd1, λd2, pe(trap), pe(trap.d1), pe(trap.d1), σ²λ)

  bb!(λprv_p, λpr, lλp, tprv, σ²λ, srδt)
  bb!(λd1v_p, lλp, λd1, td1v, σ²λ, srδt)
  bb!(λd2v_p, lλp, λd2, td2v, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tprv, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(td1v, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(td2v, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tprv, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(td1v, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(td2v, λd2v_c, σ²λ, δt)

  # speciation likelihood ratio
  λlr = 2.0*lλp - 2.0*λd1v_c[1]

  # likelihood ratio
  llr = llbmprp + llpbprp + llbmd1p + llpbd1p + llbmd2p + llpbd2p -
        llbmprc - llpbprc - llbmd1c - llpbd1c - llbmd2c - llpbd2c +
        λlr

  # acceptance ratio
  acr = llpbprp + llpbd1p + llpbd2p -
        llpbprc - llpbd1c - llpbd2c + λlr

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_root!(tra     ::iTgbmpb, 
                         trap    ::iTgbmpb,
                         llc     ::Float64,
                         prc     ::Float64,
                         σ²λ     ::Float64, 
                         srδt    ::Float64,
                         λa_prior::Tuple{Float64, Float64})

Make a trio of Brownian motion MCMC updates when the root is involved.
"""
function triad_lλupdate_root!(tra     ::iTgbmpb, 
                              trap    ::iTgbmpb,
                              llc     ::Float64,
                              prc     ::Float64,
                              σ²λ     ::Float64, 
                              srδt    ::Float64,
                              λa_prior::Tuple{Float64, Float64})

  # get reference vectors
  λprv_p = lλ(trap)
  λd1v_p = lλ(trap.d1)
  λd2v_p = lλ(trap.d2)
  λprv_c = lλ(tra)
  λd1v_c = lλ(tra.d1)
  λd2v_c = lλ(tra.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(tra)
  td1v = ts(tra.d1)
  td2v = ts(tra.d2)

  # proposal given daughters
  lλp = duoprop(λd1, λd2, pe(tra.d1), pe(tra.d2), σ²λ)

  # propose for root
  lλrp = rnorm(lλp, pe(tra)*σ²λ)

  # make Brownian bridge proposals
  bb!(λprv_p, lλrp, lλp, tprv, σ²λ, srδt)
  bb!(λd1v_p, lλp,  λd1, td1v, σ²λ, srδt)
  bb!(λd2v_p, lλp,  λd2, td2v, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tprv, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(td1v, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(td2v, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tprv, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(td1v, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(td2v, λd2v_c, σ²λ, δt)

  # speciation likelihood ratio
  λlr = 2.0*lλp - 2.0*λd1v_c[1]

  # likelihood ratio
  llr = llbmprp + llpbprp + llbmd1p + llpbd1p + llbmd2p + llpbd2p -
        llbmprc - llpbprc - llbmd1c - llpbd1c - llbmd2c - llpbd2c +
        λlr

  # prior ratio
  prr = llrdnorm_x(lλrp, λpr, λa_prior[1], λa_prior[2])

  # acceptance ratio
  acr = llpbprp + llpbd1p + llpbd2p -
        llpbprc - llpbd1c - llpbd2c + λlr

  if -randexp() < (acr + prr)
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc, prc
end




"""
    update_σ²λ!(σ²λc ::Float64,
                llc  ::Float64,
                prc  ::Float64,
                σ²λtn::Float64,
                δt   ::Float64,
                σ²λprior::Float64)

MCMC update for σ²λ.
"""
function update_σ²λ!(σ²λc ::Float64,
                     tree ::iTgbmpb,
                     llc  ::Float64,
                     prc  ::Float64,
                     σ²λtn::Float64,
                     δt   ::Float64,
                     σ²λprior::Float64)

  # parameter proposals
  σ²λp = mulupt(σ²λc, σ²λtn)::Float64

  # one could make a ratio likelihood function
  llp = llik_gbm(tree, σ²λp, δt)
  prr = llrdexp_x(σ²λp, σ²λc, σ²λprior)

  if -randexp() < (llp - llc + prr + log(σ²λp/σ²λc))
    llc  = llp::Float64
    prc += prr::Float64
    σ²λc = σ²λp::Float64
  end

  return llc, prc, σ²λc
end



