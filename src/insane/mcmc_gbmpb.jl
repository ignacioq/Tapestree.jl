#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#





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
                                 prc     ::Float64,
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
  tvpr = ts(tra)
  tvd1 = ts(tra.d1)
  tvd2 = ts(tra.d2)

  # fill with Brownian motion
  bm!(λvpr_p, λpr, tvpr, σ²λ, srδt)
  bm!(λvd1_p, λvpr_p[end], tvd1, σ²λ, srδt)
  bm!(λvd2_p, λvpr_p[end], tvd2, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tvpr, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(tvd1, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(tvd2, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tvpr, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(tvd1, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(tvd2, λd2v_c, σ²λ, δt)

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

  if -randexp() < (ar + prr)
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc, prc
end




"""
    triad_lλupdate_noded1!(tra     ::iTgbmpb, 
                           trap    ::iTgbmpb,
                           llc     ::Float64,
                           prc     ::Float64,
                           σ²λ     ::Float64, 
                           srδt    ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 1 is terminal.
"""
function triad_lλupdate_noded1!(tra     ::iTgbmpb, 
                                trap    ::iTgbmpb,
                                llc     ::Float64,
                                prc     ::Float64,
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
  tvpr = ts(tra)
  tvd1 = ts(tra.d1)
  tvd2 = ts(tra.d2)

  # node proposal
  lλp = duoprop(λpr, λd2, pe(tra), pe(tra.d2), σ²λ)

  bb!(λvpr_p, λpr, lλp, tvpr, σ²λ, srδt)
  bm!(λvd1_p, lλp, tvd1, σ²λ, srδt)
  bb!(λvd2_p, lλp, λd2, tvd2, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tvpr, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(tvd1, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(tvd2, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tvpr, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(tvd1, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(tvd2, λd2v_c, σ²λ, δt)

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

  if -randexp() < (ar + prr)
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc, prc
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
                                prc     ::Float64,
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
  tvpr = ts(tra)
  tvd1 = ts(tra.d1)
  tvd2 = ts(tra.d2)

  # node proposal
  lλp = duoprop(λpr, λd1, pe(tra), pe(tra.d1), σ²λ)

  bb!(λvpr_p, λpr, lλp, tvpr, σ²λ, srδt)
  bb!(λvd1_p, lλp, λd1, tvd1, σ²λ, srδt)
  bm!(λvd2_p, lλp, tvd2, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tvpr, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(tvd1, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(tvd2, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tvpr, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(tvd1, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(tvd2, λd2v_c, σ²λ, δt)

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

  if -randexp() < (ar + prr)
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc, prc
end




"""
    triad_lλupdate_node!(tra     ::iTgbmpb, 
                          trap    ::iTgbmpb,
                          llc     ::Float64,
                          prc     ::Float64,
                          σ²λ     ::Float64, 
                          srδt    ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
no daughters are terminal.
"""
function triad_lλupdate_node!(tra     ::iTgbmpb, 
                               trap    ::iTgbmpb,
                               llc     ::Float64,
                               prc     ::Float64,
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
  tvpr = ts(tra)
  tvd1 = ts(tra.d1)
  tvd2 = ts(tra.d2)

  # node proposal
  lλp = trioprop(λpr, λd1, λd2, pe(trap), pe(trap.d1), pe(trap.d1), σ²λ)

  bb!(λprv_p, λpr, lλp, tvpr, σ²λ, srδt)
  bb!(λd1v_p, lλp, λd1, tvd1, σ²λ, srδt)
  bb!(λd2v_p, lλp, λd2, tvd2, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tvpr, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(tvd1, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(tvd2, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tvpr, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(tvd1, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(tvd2, λd2v_c, σ²λ, δt)

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

  if -randexp() < (ar + prr)
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc, prc
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
  tvpr = ts(tra)
  tvd1 = ts(tra.d1)
  tvd2 = ts(tra.d2)

  # proposal given daughters
  lλp = duoprop(λd1, λd2, pe(tra.d1), pe(tra.d2), σ²λ)

  # propose for root
  lλrp = rnorm(lλp, pe(tra)*σ²λ)

  # make Brownian bridge proposals
  bb!(λprv_p, lλrp, lλp, tvpr, σ²λ, srδt)
  bb!(λd1v_p, lλp,  λd1, tvd1, σ²λ, srδt)
  bb!(λd2v_p, lλp,  λd2, tvd2, σ²λ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods separately
  llbmprp,llpbprp = ll_gbm_b_sep(tvpr, λprv_p, σ²λ, δt)
  llbmd1p,llpbd1p = ll_gbm_b_sep(tvd1, λd1v_p, σ²λ, δt)
  llbmd2p,llpbd2p = ll_gbm_b_sep(tvd2, λd2v_p, σ²λ, δt)
  llbmprc,llpbprc = ll_gbm_b_sep(tvpr, λprv_c, σ²λ, δt)
  llbmd1c,llpbd1c = ll_gbm_b_sep(tvd1, λd1v_c, σ²λ, δt)
  llbmd2c,llpbd2c = ll_gbm_b_sep(tvd2, λd2v_c, σ²λ, δt)

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

  if -randexp() < (ar + prr)
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc, prc
end



