#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#





"""
    insane_gbmpb(tree    ::iTpb, 
                 out_file::String;
                 λa_prior::Tuple{Float64,Float64} = (0.0,10.0),
                 σλprior ::Float64  = 0.1,
                 δt      ::Float64  = 1e-2,
                 niter   ::Int64    = 1_000,
                 nthin   ::Int64    = 10,
                 nburn   ::Int64    = 200,
                 tune_int::Int64    = 100,
                 σλi     ::Float64  = 0.5,
                 σλtni   ::Float64  = 1.0,
                 obj_ar  ::Float64  = 0.234,
                 pupdp   ::Tuple{Float64,Float64} = (0.9, 0.1),
                 prints  ::Int64    = 5)

Run insane for constant pure-birth.
"""
function insane_gbmpb(tree    ::iTpb, 
                      out_file::String;
                      σλprior ::Float64  = 0.1,
                      δt      ::Float64  = 1e-2,
                      niter   ::Int64    = 1_000,
                      nthin   ::Int64    = 10,
                      nburn   ::Int64    = 200,
                      tune_int::Int64    = 100,
                      σλi     ::Float64  = 0.5,
                      σλtni   ::Float64  = 1.0,
                      obj_ar  ::Float64  = 0.234,
                      prints  ::Int64    = 5,
                      pupdp   ::Tuple{Float64,Float64} = (0.9, 0.1),
                      λa_prior::Tuple{Float64,Float64} = (0.0,10.0))

  δt  *= treeheight(tree)
  srδt = sqrt(δt)

  # lλ root node
  lλa = log(λmle_cpb(tree))

  # make Ψ current and proposal parameters
  Ψc = iTgbmpb(tree, δt, srδt, lλa, σλi)
  Ψp = deepcopy(Ψc)

  # make fix Ψ directory
  idv = iDir[]
  bit = BitArray{1}()
  makeiDir!(Ψc, idv, bit)

  # make parent node directory to `iDir`
  inodes, terminus = make_inodes(idv)

  # parameter update vector
  nin = lastindex(inodes)

  pup = make_pup(pupdp, nin)

  # make scaling function
  scalef = makescalef(obj_ar)

  # burn-in phase
  llc, prc, σλc, σλtn =
    mcmc_burn_gbmpb(Ψp, Ψc, λa_prior, σλprior, nburn, tune_int, σλi, σλtni, 
      δt, srδt, idv, inodes, terminus, pup, prints, scalef)

  # mcmc
  R, Ψv = mcmc_gbmpb(Ψp, Ψc, llc, prc, σλc, λa_prior, σλprior, 
        niter, nthin, σλtn, δt, srδt, idv, inodes, terminus, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "sigma_lambda" => 2))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmpb(Ψp      ::iTgbmpb,
                    Ψc      ::iTgbmpb,
                    λa_prior::Tuple{Float64,Float64},
                    σλprior ::Float64,
                    nburn   ::Int64,
                    tune_int::Int64,
                    σλi     ::Float64,
                    σλtni   ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idv     ::Array{iDir,1},
                    inodes  ::Array{Int64,1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    prints  ::Int64,
                    scalef  ::Function)

MCMC burn-in chain for GBM pure-birth.
"""
function mcmc_burn_gbmpb(Ψp      ::iTgbmpb,
                         Ψc      ::iTgbmpb,
                         λa_prior::Tuple{Float64,Float64},
                         σλprior ::Float64,
                         nburn   ::Int64,
                         tune_int::Int64,
                         σλi     ::Float64,
                         σλtni   ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         idv     ::Array{iDir,1},
                         inodes  ::Array{Int64,1},
                         terminus::Array{BitArray{1}},
                         pup     ::Array{Int64,1},
                         prints  ::Int64,
                         scalef  ::Function)

  # initialize acceptance log
  ltn  = 0
  lup  = 0.0
  lac  = 0.0
  σλtn = σλtni

  # starting parameters
  σλc = σλi
  llc = llik_gbm(Ψc, σλc, δt, srδt)
  prc = logdexp(σλc, σλprior) + 
        logdnorm_tc(lλ(Ψc)[1], λa_prior[1], λa_prior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup
      ## parameter updates
      if iszero(pupi)
        llc, prc, σλc, ltn, lup, lac = 
          update_σλ!(σλc, Ψc, llc, prc, 
            σλtn, ltn, lup, lac, δt, srδt, σλprior)
      else
        ter  = terminus[pupi]
        inod = inodes[pupi]

        dri = dr(idv[inod])
        ldr = lastindex(dri)

        llc, prc = 
          lλupdate!(Ψp, Ψc, llc, prc, σλc, δt, srδt, 
            λa_prior, dri, ldr, ter, 0)
      end

      if ltn == tune_int
        σλtn = scalef(σλtn,lac/lup)
        ltn = 0
      end

    end

    next!(pbar)
  end

  return llc, prc, σλc, σλtn
end




"""
    mcmc_gbmpb(Ψp      ::iTgbmpb,
               Ψc      ::iTgbmpb,
               llc     ::Float64,
               prc     ::Float64,
               σλc     ::Float64,
               λa_prior::Tuple{Float64,Float64},
               σλprior ::Float64,
               niter   ::Int64,
               nthin   ::Int64,
               σλtn    ::Float64,
               δt      ::Float64,
               srδt    ::Float64,
               idv     ::Array{iDir,1},
               inodes  ::Array{Int64,1},
               terminus::Array{BitArray{1}},
               pup     ::Array{Int64,1},
               prints  ::Int64)

MCMC chain for GBM pure-birth.
"""
function mcmc_gbmpb(Ψp      ::iTgbmpb,
                    Ψc      ::iTgbmpb,
                    llc     ::Float64,
                    prc     ::Float64,
                    σλc     ::Float64,
                    λa_prior::Tuple{Float64,Float64},
                    σλprior ::Float64,
                    niter   ::Int64,
                    nthin   ::Int64,
                    σλtn    ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idv     ::Array{iDir,1},
                    inodes  ::Array{Int64,1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  R = Array{Float64,2}(undef, nlogs, 5)

  # make Ψ vector
  Ψv = iTgbmpb[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      if iszero(pupi)
        # `λ` diffusion rate updates
        llc, prc, σλc = 
          update_σλ!(σλc, Ψc, llc, prc, σλtn, δt, srδt, σλprior)
      else
        # gbm updates
        ter  = terminus[pupi]
        inod = inodes[pupi]

        dri = dr(idv[inod])
        ldr = lastindex(dri)

        llc, prc = 
          lλupdate!(Ψp, Ψc, llc, prc, σλc, δt, srδt, 
            λa_prior, dri, ldr, ter, 0)
      end
    end

    # log parameters
    lthin += 1
    if lthin == nthin
      lit += 1
      @inbounds begin
        R[lit,1] = Float64(lit)
        R[lit,2] = llc
        R[lit,3] = prc
        R[lit,4] = exp(lλ(Ψc)[1])
        R[lit,5] = σλc
        push!(Ψv, deepcopy(Ψc))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ψv
end




"""
    lλupdate!(Ψp      ::iTgbmpb,
              Ψc      ::iTgbmpb,
              llc     ::Float64, 
              prc     ::Float64,
              σλ      ::Float64, 
              δt      ::Float64, 
              srδt    ::Float64, 
              λa_prior::Tuple{Float64,Float64},
              dri     ::BitArray{1},
              ldr     ::Int64,
              ter     ::BitArray{1},
              ix      ::Int64)

Make a λgbm update for a triad.
"""
function lλupdate!(Ψp      ::iTgbmpb,
                   Ψc      ::iTgbmpb,
                   llc     ::Float64, 
                   prc     ::Float64,
                   σλ      ::Float64, 
                   δt      ::Float64, 
                   srδt    ::Float64, 
                   λa_prior::Tuple{Float64,Float64},
                   dri     ::BitArray{1},
                   ldr     ::Int64,
                   ter     ::BitArray{1},
                   ix      ::Int64)

  if ix == ldr 
    # if root
    if ldr == 0
      llc, prc = triad_lλupdate_root!(Ψp::iTgbmpb, Ψc::iTgbmpb, 
                   llc, prc, σλ, δt, srδt, λa_prior)
    else
      if ter[1]
        if ter[2]
          # if both are terminal
          llc = triad_lλupdate_noded12!(Ψp::iTgbmpb, Ψc::iTgbmpb, 
                        llc, σλ, δt, srδt)
        else
          # if d1 is terminal
          llc = triad_lλupdate_noded1!(Ψp::iTgbmpb, Ψc::iTgbmpb, 
                       llc, σλ, δt, srδt)
        end
      elseif ter[2]
        # if d2 is terminal
        llc = triad_lλupdate_noded2!(Ψp::iTgbmpb, Ψc::iTgbmpb, 
                     llc, σλ, δt, srδt)
      else
        # if no terminal branches involved
        llc = triad_lλupdate_node!(Ψp::iTgbmpb, Ψc::iTgbmpb, 
                     llc, σλ, δt, srδt)
      end
    end

  elseif ix < ldr
    ix += 1
    if dri[ix]
      llc, prc = 
        lλupdate!(Ψp.d1::iTgbmpb, Ψc.d1::iTgbmpb, 
          llc, prc, σλ, δt, srδt, λa_prior, dri, ldr, ter, ix)
    else
      llc, prc = 
        lλupdate!(Ψp.d2::iTgbmpb, Ψc.d2::iTgbmpb, 
          llc, prc, σλ, δt, srδt, λa_prior, dri, ldr, ter, ix)
    end
  end

  return llc, prc
end




"""
    triad_lλupdate_noded12!(Ψc  ::iTgbmpb, 
                            Ψp  ::iTgbmpb,
                            llc ::Float64,
                            σλ  ::Float64, 
                            δt  ::Float64)
                            srδt::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
both daughters are terminal.
"""
function triad_lλupdate_noded12!(Ψp  ::iTgbmpb, 
                                 Ψc  ::iTgbmpb,
                                 llc ::Float64,
                                 σλ  ::Float64, 
                                 δt  ::Float64,
                                 srδt::Float64)

  # get reference vectors
  λprv_p = lλ(Ψp)
  λd1v_p = lλ(Ψp.d1)
  λd2v_p = lλ(Ψp.d2)
  λprv_c = lλ(Ψc)
  λd1v_c = lλ(Ψc.d1)
  λd2v_c = lλ(Ψc.d2)
  λpr = λprv_c[1]

  # time vectors
  tprv = ts(Ψc)
  td1v = ts(Ψc.d1)
  td2v = ts(Ψc.d2)

  # fill with Brownian motion
  bm!(λprv_p, λpr, tprv, σλ, srδt)
  lλp = λprv_p[end]
  bm!(λd1v_p, lλp, td1v, σλ, srδt)
  bm!(λd2v_p, lλp, td2v, σλ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_noded1!(Ψc::iTgbmpb, 
                           Ψp::iTgbmpb,
                           llc  ::Float64,
                           σλ   ::Float64, 
                           δt   ::Float64)
                           srδt ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 1 is terminal.
"""
function triad_lλupdate_noded1!(Ψp  ::iTgbmpb, 
                                Ψc  ::iTgbmpb,
                                llc ::Float64,
                                σλ  ::Float64, 
                                δt  ::Float64,
                                srδt::Float64)

  # get reference vectors
  λprv_p = lλ(Ψp)
  λd1v_p = lλ(Ψp.d1)
  λd2v_p = lλ(Ψp.d2)
  λprv_c = lλ(Ψc)
  λd1v_c = lλ(Ψc.d1)
  λd2v_c = lλ(Ψc.d2)
  λpr = λprv_c[1]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(Ψc)
  td1v = ts(Ψc.d1)
  td2v = ts(Ψc.d2)

  # node proposal
  lλp = duoprop(λpr, λd2, pe(Ψc), pe(Ψc.d2), σλ)

  bb!(λprv_p, λpr, lλp, tprv, σλ, srδt)
  bm!(λd1v_p, lλp, td1v, σλ, srδt)
  bb!(λd2v_p, lλp, λd2, td2v, σλ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_noded2!(Ψc::iTgbmpb, 
                           Ψp::iTgbmpb,
                           llc  ::Float64,
                           σλ   ::Float64, 
                           δt   ::Float64,
                           srδt ::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
daughter 2 is terminal.
"""
function triad_lλupdate_noded2!(Ψp  ::iTgbmpb, 
                                Ψc  ::iTgbmpb,
                                llc ::Float64,
                                σλ  ::Float64, 
                                δt  ::Float64,
                                srδt::Float64)

  # get reference vectors
  λprv_p = lλ(Ψp)
  λd1v_p = lλ(Ψp.d1)
  λd2v_p = lλ(Ψp.d2)
  λprv_c = lλ(Ψc)
  λd1v_c = lλ(Ψc.d1)
  λd2v_c = lλ(Ψc.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]

  # time vectors
  tprv = ts(Ψc)
  td1v = ts(Ψc.d1)
  td2v = ts(Ψc.d2)

  # node proposal
  lλp = duoprop(λpr, λd1, pe(Ψc), pe(Ψc.d1), σλ)

  bb!(λprv_p, λpr, lλp, tprv, σλ, srδt)
  bb!(λd1v_p, lλp, λd1, td1v, σλ, srδt)
  bm!(λd2v_p, lλp, td2v, σλ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  if -randexp() < acr 
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_node!(Ψp  ::iTgbmpb, 
                         Ψc  ::iTgbmpb,
                         llc ::Float64,
                         σλ  ::Float64,
                         δt  ::Float64, 
                         srδt::Float64)

Make a trio of Brownian motion MCMC updates when node is internal and 
no daughters are terminal.
"""
function triad_lλupdate_node!(Ψp  ::iTgbmpb, 
                              Ψc  ::iTgbmpb,
                              llc ::Float64,
                              σλ  ::Float64,
                              δt  ::Float64, 
                              srδt::Float64)

  # get reference vectors
  λprv_p = lλ(Ψp)
  λd1v_p = lλ(Ψp.d1)
  λd2v_p = lλ(Ψp.d2)
  λprv_c = lλ(Ψc)
  λd1v_c = lλ(Ψc.d1)
  λd2v_c = lλ(Ψc.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(Ψc)
  td1v = ts(Ψc.d1)
  td2v = ts(Ψc.d2)

  # node proposal
  lλp = trioprop(λpr, λd1, λd2, pe(Ψc), pe(Ψc.d1), pe(Ψc.d2), σλ)

  bb!(λprv_p, λpr, lλp, tprv, σλ, srδt)
  bb!(λd1v_p, lλp, λd1, td1v, σλ, srδt)
  bb!(λd2v_p, lλp, λd2, td2v, σλ, srδt)

  ## make acceptance ratio 
  # estimate likelihoods
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  if -randexp() < acr
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    triad_lλupdate_root!(Ψc      ::iTgbmpb, 
                         Ψp      ::iTgbmpb,
                         llc     ::Float64,
                         prc     ::Float64,
                         σλ      ::Float64, 
                         δt      ::Float64,
                         σλ      ::Float64, 
                         δt      ::Float64,
                         srδt    ::Float64,
                         λa_prior::Tuple{Float64, Float64})

Make a trio of Brownian motion MCMC updates when the root is involved.
"""
function triad_lλupdate_root!(Ψp      ::iTgbmpb, 
                              Ψc      ::iTgbmpb,
                              llc     ::Float64,
                              prc     ::Float64,
                              σλ      ::Float64, 
                              δt      ::Float64,
                              srδt    ::Float64,
                              λa_prior::Tuple{Float64, Float64})

  # get reference vectors
  λprv_p = lλ(Ψp)
  λd1v_p = lλ(Ψp.d1)
  λd2v_p = lλ(Ψp.d2)
  λprv_c = lλ(Ψc)
  λd1v_c = lλ(Ψc.d1)
  λd2v_c = lλ(Ψc.d2)
  λpr = λprv_c[1]
  λd1 = λd1v_c[end]
  λd2 = λd2v_c[end]

  # time vectors
  tprv = ts(Ψc)
  td1v = ts(Ψc.d1)
  td2v = ts(Ψc.d2)

  # proposal given daughters
  lλp = duoprop(λd1, λd2, pe(Ψc.d1), pe(Ψc.d2), σλ)

  # propose for root
  lλrp = rnorm(lλp, sqrt(pe(Ψc))*σλ)

  # make Brownian bridge proposals
  bb!(λprv_p, lλrp, lλp, tprv, σλ, srδt)
  bb!(λd1v_p, lλp,  λd1, td1v, σλ, srδt)
  bb!(λd2v_p, lλp,  λd2, td2v, σλ, srδt)

  ## make acceptance ratio 
  llr, acr = llr_propr(tprv, td1v, td2v, λprv_p, λd1v_p, λd2v_p, 
    λprv_c, λd1v_c, λd2v_c, σλ, δt, srδt, lλp, λd1v_c[1])

  # prior ratio
  prr = llrdnorm_x(lλrp, λpr, λa_prior[1], λa_prior[2])

  # acceptance ratio
  acr += prr

  if -randexp() < acr 
    llc += llr
    prc += prr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc, prc
end




"""
    llr_propr(tprv  ::Array{Float64,1},
              td1v  ::Array{Float64,1},
              td2v  ::Array{Float64,1},
              λprv_p::Array{Float64,1},
              λd1v_p::Array{Float64,1},
              λd2v_p::Array{Float64,1},
              λprv_c::Array{Float64,1},
              λd1v_c::Array{Float64,1},
              λd2v_c::Array{Float64,1},
              σλ    ::Float64,
              σλ   ::Float64,
              δt    ::Float64,
              srδt  ::Float64,
              lλp   ::Float64,
              lλc   ::Float64)

Return the likelihood and proposal ratio for pure-birth gbm.
"""
function llr_propr(tprv  ::Array{Float64,1},
                   td1v  ::Array{Float64,1},
                   td2v  ::Array{Float64,1},
                   λprv_p::Array{Float64,1},
                   λd1v_p::Array{Float64,1},
                   λd2v_p::Array{Float64,1},
                   λprv_c::Array{Float64,1},
                   λd1v_c::Array{Float64,1},
                   λd2v_c::Array{Float64,1},
                   σλ    ::Float64,
                   δt    ::Float64,
                   srδt  ::Float64,
                   lλp   ::Float64,
                   lλc   ::Float64)

  # log likelihood ratio functions
  llrbm_pr, llrpb_pr = llr_gbm_b_sep(tprv, λprv_p, λprv_c, σλ, δt, srδt)
  llrbm_d1, llrpb_d1 = llr_gbm_b_sep(td1v, λd1v_p, λd1v_c, σλ, δt, srδt)
  llrbm_d2, llrpb_d2 = llr_gbm_b_sep(td2v, λd2v_p, λd2v_c, σλ, δt, srδt)

  llr = llrbm_pr + llrpb_pr +
        llrbm_d1 + llrpb_d1 +
        llrbm_d2 + llrpb_d2 +
        (lλp - lλc)

  acr = llr - llrbm_pr - llrbm_d1 - llrbm_d2 

  return llr, acr
end




"""
    update_σλ!(σλc    ::Float64,
               Ψ      ::iTgbmpb,
               llc    ::Float64,
               prc    ::Float64,
               σλtn   ::Float64,
               δt     ::Float64,
               srδt   ::Float64,
               σλprior::Float64)

MCMC update for σλ.
"""
function update_σλ!(σλc    ::Float64,
                    Ψ      ::iTgbmpb,
                    llc    ::Float64,
                    prc    ::Float64,
                    σλtn   ::Float64,
                    δt     ::Float64,
                    srδt   ::Float64,
                    σλprior::Float64)

  # parameter proposals
  σλp = mulupt(σλc, σλtn)::Float64

  # log likelihood and prior ratio
  llr = llr_gbm_bm(Ψ, σλp, σλc, srδt)
  prr = llrdexp_x(σλp, σλc, σλprior)

  if -randexp() < (llr + prr + log(σλp/σλc))
    σλc  = σλp
    llc += llr
    prc += prr
  end

  return llc, prc, σλc
end




"""
    update_σλ!(σλc    ::Float64,
               Ψ      ::iTgbmpb,
               llc    ::Float64,
               prc    ::Float64,
               σλtn   ::Float64,
               ltn    ::Int64,
               lup    ::Float64,
               lac    ::Float64,
               δt     ::Float64,
               srδt   ::Float64,
               σλprior::Float64)

MCMC update for σλ with acceptance log.
"""
function update_σλ!(σλc    ::Float64,
                    Ψ      ::iTgbmpb,
                    llc    ::Float64,
                    prc    ::Float64,
                    σλtn   ::Float64,
                    ltn    ::Int64,
                    lup    ::Float64,
                    lac    ::Float64,
                    δt     ::Float64,
                    srδt   ::Float64,
                    σλprior::Float64)

  # parameter proposals
  σλp = mulupt(σλc, σλtn)::Float64

  # log likelihood and prior ratio
  llr = llr_gbm_bm(Ψ, σλp, σλc, srδt)
  prr = llrdexp_x(σλp, σλc, σλprior)

  ltn += 1
  lup += 1.0

  if -randexp() < (llr + prr + log(σλp/σλc))
    σλc  = σλp
    llc += llr
    prc += prr
    lac += 1.0
  end

  return llc, prc, σλc, ltn, lup, lac
end




"""
    ncrep!(Ψp::iTgbmpb, 
           Ψc::iTgbmpb,
           σ    ::Float64)

Non-centered reparametization of data augmentation for `σ`.
"""
function ncrep!(Ψp::iTgbmpb, 
                Ψc::iTgbmpb,
                σ    ::Float64)

  ncrep!(lλ(Ψp), lλ(Ψc), ts(Ψc), σ)

  if !istip(Ψc.d1)
    ncrep!(Ψp.d1::iTgbmpb, Ψc.d1::iTgbmpb, σ)
  end
  if !istip(Ψc.d2)
    ncrep!(Ψp.d2::iTgbmpb, Ψc.d2::iTgbmpb, σ)
  end
end




"""
    make_pup(pupdp::NTuple{N,Float64}, 
             nin  ::Int64) where {N}

Make the weighted parameter update vector according to probabilities `pupdp`.
"""
function make_pup(pupdp::NTuple{N,Float64}, 
                  nin  ::Int64) where {N}

  # standardize pr vector
  pups = Array{Float64,1}(undef,N)
  spupdp = sum(pupdp)
  for i in Base.OneTo(N)
    pups[i] = pupdp[i]/spupdp
  end

  pup = Int64[]
  # da parameters
  if pups[1] > 0.0
    append!(pup,[1:nin...])
  end

  if pups[2] > 0.0
    if pups[1] > 0.0
      append!(pup, 
        fill(0, ceil(Int64, pups[2]*Float64(nin)/pups[1])))
    else 
      push!(pup,0)
    end
  end

  return return pup
end


