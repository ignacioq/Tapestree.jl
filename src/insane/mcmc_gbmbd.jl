#=

Anagenetic GBM birth-death MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#

function insane_gbmbd(tree    ::sTbd, 
                      out_file::String;
                      λprior  ::Float64           = 0.1,
                      μprior  ::Float64           = 0.1,
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      tune_int::Int64             = 100,
                      ϵi      ::Float64           = 0.2,
                      λi      ::Float64           = NaN,
                      μi      ::Float64           = NaN,
                      λtni    ::Float64           = 1.0,
                      μtni    ::Float64           = 1.0,
                      obj_ar  ::Float64           = 0.4,
                      pupdp   ::NTuple{4,Float64} = (0.4,0.4,0.1,0.1),
                      prints  ::Int64              = 5)

  # fix tree
  fixtree!(tree)

  # `n` tips, `th` treeheight define δt
  n    = sntn(tree)
  th   = treeheight(tree)
  δt  *= th
  srδt = sqrt(δt)

   # starting parameters (using method of moments)
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end

  # make Ψ current and proposal parameters
  Ψc = iTgbmbd(tree, δt, srδt, log(λc), log(μc), σλi, σμi)
  Ψp = deepcopy(Ψc)

  # make fix Ψ directory
  idv = iDir[]
  bit = BitArray{1}()
  makeiDir!(Ψc, idv, bit)

  # make parent node directory to `iDir`
  inodes, terminus = make_inodes(idv)

  # create wbr vector
  wbr = falses(lastindex(idv))

  # create da branches vector
  dabr = Int64[]

  # make survival conditioning function (stem or crown)
  # svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
  #                          stem_prob_surv_cbd

  scalef = makescalef(obj_ar)


  nin = lastindex(inodes)

  pupdp = (0.9, 0.1)
  
  pup = make_pup(pupdp, nin)




  # initialize acceptance log
  ltn  = 0
  lup  = 0.0
  lac  = 0.0
  σλtn = σλtni
  σμtn = σμtni

  # starting parameters
  σλc = σλi
  σμc = σμi

  llc = llik_gbm(Ψc, σλc, σμc, δt, srδt)
  prc = logdexp(σλc, σλprior)                            +
        logdexp(σμc, σμprior)                            +
        logdnorm_tc(lλ(Ψc)[1], λa_prior[1], λa_prior[2]) +
        logdnorm_tc(lμ(Ψc)[1], μa_prior[1], μa_prior[2])

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










"""
    update_σλ!(σλc    ::Float64,
               Ψ      ::iTgbmbd,
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
                    Ψ      ::iTgbmbd,
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


"""
here: start with log likelihood ratio for speciation rate variance
"""




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





