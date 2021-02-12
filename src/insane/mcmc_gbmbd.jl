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
  idf = iBf[]
  bit = BitArray{1}()
  makeiBf!(Ψc, idf, bit)

  # make parent node directory to `iBf`
  inodes, terminus = make_inodes(idf)

  # create da branches vector
  ida = iBa[]

  # create wbf and wba vector
  wbf = falses(lastindex(idf))
  wba = falses(lastindex(ida))


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
  lλup = lλac = 0.0
  σλtn = σλtni
  lμtn = 0
  lμup = lμac = 0.0
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
      # update gbm
      if pupi > 0

        """
        do graft/prune first
        """

        ter  = terminus[pupi]
        inod = inodes[pupi]

        dri = dr(idv[inod])
        ldr = lastindex(dri)

        llc, prc = 
          lvupdate!(Ψp, Ψc, llc, prc, σλc, δt, srδt, 
            λa_prior, dri, ldr, ter, 0)

      # graft
      elseif pupi == 0

        """
        first graph and then update gbm (to account for grafted gbm)
        """





      # prune
      elseif pupi == -1


      # update σλ
      elseif pupi == -2
        llc, prc, σλc, lλtn, lλup, lλac = 
          update_σ!(σλc, Ψc, llc, prc, 
            σλtn, lλtn, lλup, lλac, δt, srδt, σλprior, lλ)

      # update σμ
      else
        llc, prc, σμc, lμtn, lμup, lμac = 
          update_σ!(σμc, Ψc, llc, prc, 
            σμtn, lμtn, lμup, lμac, δt, srδt, σμprior, lμ)
      end


      if ltn == tune_int
        σλtn = scalef(σλtn,lac/lup)
        ltn = 0
      end

    end

    next!(pbar)
  end







"""
    graftp(tree::sTbd,
           llc ::Float64,
           λc  ::Float64,
           μc  ::Float64,
           th  ::Float64,
           idv ::Array{iBf,1}, 
           wbr ::BitArray{1},
           dabr::Array{Int64,1},
           pupdp::NTuple{4,Float64})

Graft proposal function for constant birth-death.
"""
function graftp(Ψ    ::iTgbmbd,
                llc  ::Float64,
                λc   ::Float64,
                μc   ::Float64,
                th   ::Float64,
                idf  ::Array{iBf,1}, 
                ida  ::Array{iBa,1}, 
                wbf  ::BitArray{1},
                wba  ::BitArray{1},
                dabr ::Array{Int64,1},
                pupdp::NTuple{4,Float64},
                δt   ::Float64)
                

#=
check if we can choose the tree length first to simulate from the
lambda, extinction
=#

  λn, μn = λμprop() 

  #simulate extinct lineage
  t0, t0h = sim_cbd_b(λn, μn, th, 100)

  # if useful simulation
  if t0h < th

    # get height and number of intersecting branches
    h, nf, na, rn = randbranch(th, t0h, wbf, wba, idf, ida)

    # if branch is fixed
    if rn <= nf
      bf, i = getbranch(rn, wbf, idf)

      dri = dr(bf)
      ldr = lastindex(dri)
      dai = da(bf)

      # check `δt` to graft to (and index `hi`) and it's current 
      # `λ`, `μ` at `nh`
      λh, μh, nh, hi = λμath(Ψ, h, th, dri, ldr, 0)

    # if branch is from da
    else
      ba, i = getbranch(rn - nf, wba, ida)



    end




    """
    here
    """

    # make likelihood for GBM but knowing all rates are constant
    # also for proposal ratio




    # proposal ratio
    lpr = log(2.0 * μc * (th - t0h) * Float64(nbh) * pupdp[4]) - 
          log((Float64(lastindex(dabr)) + 1.0) * pupdp[3])

    # likelihood ratio
    llr = llik_cbd(t0, λc, μc) + log(2.0*λc) 

    if -randexp() < lpr #+ llr
      llc += llr



      # gbm augment tree with constant rates
      # but do the augmentation only if grafted
      t0 = iTgbmbd(t0, δt, srδt, log(0.1), log(0.1), 0.0, 0.0)



      # graft branch
      dri  = dr(br)
      tree = graftree!(tree, t0, dri, h, lastindex(dri), th, 0)
      # add n graft to branch
      addda!(br)
      # log branch as being data augmented
      push!(dabr, bri)
    end
  end

  return tree, llc
end










"""
    lvupdate!(Ψp      ::iTgbmbd,
              Ψc      ::iTgbmbd,
              llc     ::Float64, 
              prc     ::Float64,
              σ      ::Float64, 
              δt      ::Float64, 
              srδt    ::Float64, 
              a_prior::Tuple{Float64,Float64},
              dri     ::BitArray{1},
              ldr     ::Int64,
              ter     ::BitArray{1},
              ix      ::Int64)

Make a gbm update for a triad.
"""
function lvupdate!(Ψp      ::iTgbmbd,
                   Ψc      ::iTgbmbd,
                   llc     ::Float64, 
                   prc     ::Float64,
                   σ       ::Float64, 
                   δt      ::Float64, 
                   srδt    ::Float64, 
                   a_prior::Tuple{Float64,Float64},
                   dri     ::BitArray{1},
                   ldr     ::Int64,
                   ter     ::BitArray{1},
                   ix      ::Int64)

  """
  here
  """

  if ix == ldr 
    # if root
    if ldr == 0
      llc, prc = triad_lupdate_root!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                   llc, prc, σ, δt, srδt, a_prior)
    else
      if ter[1]
        if ter[2]
          # if both are terminal
          llc = triad_lupdate_noded12!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                        llc, σ, δt, srδt)
        else
          # if d1 is terminal
          llc = triad_lupdate_noded1!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                       llc, σ, δt, srδt)
        end
      elseif ter[2]
        # if d2 is terminal
        llc = triad_lupdate_noded2!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                     llc, σ, δt, srδt)
      else
        # if no terminal branches involved
        llc = triad_lupdate_node!(Ψp::iTgbmbd, Ψc::iTgbmbd, 
                     llc, σ, δt, srδt)
      end
    end

  elseif ix < ldr
    ix += 1
    if dri[ix]
      llc, prc = 
        lvupdate!(Ψp.d1::iTgbmbd, Ψc.d1::iTgbmbd, 
          llc, prc, σ, δt, srδt, a_prior, dri, ldr, ter, ix)
    else
      llc, prc = 
        lvupdate!(Ψp.d2::iTgbmbd, Ψc.d2::iTgbmbd, 
          llc, prc, σ, δt, srδt, a_prior, dri, ldr, ter, ix)
    end
  end

  return llc, prc
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




