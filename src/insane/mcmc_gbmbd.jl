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
                      pupdp   ::NTuple{3,Float64} = (0.8,0.2,0.1),
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

  # number of internal nodes
  nin = lastindex(inodes)

  # parameter updates (gbm, graft, σλ & σμ)
  pup = make_pup(pupdp, nin)

  # initialize acceptance log
  lλtn = 0
  lλup = lλac =0.0
  σλtn = σλtni
  lμtn = 0
  lμup = lμac =0.0
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
      # graft
      if pupi == 0

      # prune
      elseif pupi == -1

      # update σλ
      elseif pupi == -2

        llc, prc, σλc, lλtn, lλup, lλac = 
          update_σ!(σλc, Ψc, llc, prc, 
            σλtn, lλtn, lλup, lλac, δt, srδt, σλprior, lλ)

      # update σμ
      elseif pupi == -3

        llc, prc, σμc, lμtn, lμup, lμac = 
          update_σ!(σμc, Ψc, llc, prc, 
            σμtn, lμtn, lμup, lμac, δt, srδt, σμprior, lμ)

      # update GBM
      else
        """
         here
        """

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
    update_σ!(σc    ::Float64,
              Ψ     ::iTgbmbd,
              llc   ::Float64,
              prc   ::Float64,
              σtn   ::Float64,
              ltn   ::Int64,
              lup   ::Float64,
              lac   ::Float64,
              δt    ::Float64,
              srδt  ::Float64,
              σprior::Float64)

MCMC update for `σ` with acceptance log.
"""
function update_σ!(σc    ::Float64,
                   Ψ     ::iTgbmbd,
                   llc   ::Float64,
                   prc   ::Float64,
                   σtn   ::Float64,
                   ltn   ::Int64,
                   lup   ::Float64,
                   lac   ::Float64,
                   δt    ::Float64,
                   srδt  ::Float64,
                   σprior::Float64,
                   lf    ::Function)

  # parameter proposals
  σp = mulupt(σc, σtn)::Float64

  # log likelihood and prior ratio
  llr = llr_gbm_bm(Ψ, σp, σc, srδt, lf)
  prr = llrdexp_x(σp, σc, σprior)

  ltn += 1
  lup += 1.0

  if -randexp() < (llr + prr + log(σp/σc))
    σc   = σp
    llc += llr
    prc += prr
    lac += 1.0
  end

  return llc, prc, σc, ltn, lup, lac
end







"""
    make_pup(pupdp::NTuple{3,Float64}, 
             nin  ::Int64)

Make the weighted parameter update vector according to probabilities `pupdp`.
"""
function make_pup(pupdp::NTuple{3,Float64}, 
                  nin  ::Int64)

  # standardize pr vector
  pups = Array{Float64,1}(undef,3)
  spupdp = sum(pupdp)
  for i in Base.OneTo(3)
    pups[i] = pupdp[i]/spupdp
  end

  pup = Int64[]
  # `gbm` parameters
  if pups[1] > 0.0
    append!(pup,[1:nin...])
  end

  # probability `nin` given 
  prnin = pups[1] > 0.0 ? Float64(nin)/pups[1] : 0.0

  # `graft/prune` parameters (0,-1)
  if pups[2] > 0.0
    if prnin > 0.0
      append!(pup, 
        fill(-1, ceil(Int64, pups[2]*prnin)))
      append!(pup, 
        fill(0, ceil(Int64, pups[2]*prnin)))
    else
      push!(pup, -1, 0)
    end
  end

  # `σλ & σμ` parameters (-2,-3)
  if pups[3] > 0.0
    if prnin > 0.0
      append!(pup, 
        fill(-2, ceil(Int64, pups[3]*prnin)))
      append!(pup, 
        fill(-3, ceil(Int64, pups[3]*prnin)))
    else
      push!(pup, -2, -3)
    end
  end

  return return pup
end


