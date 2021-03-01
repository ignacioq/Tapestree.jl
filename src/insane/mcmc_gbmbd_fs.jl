#=

Anagenetic GBM birth-death MCMC using forward simulation

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


  # make survival conditioning function (stem or crown)
  # svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
  #                          stem_prob_surv_cbd

  scalef = makescalef(obj_ar)

  # parameter updates (1: σλ & σμ, 2: gbm, 3: forward simulation,)
  pup = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, Int64(100.0 * pupdp[i])))
  end

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
      # update σλ or σμ
      if pupi === 1

        llc, prc, σλc, lλtn, lλup, lλac = 
          update_σ!(σλc, Ψc, llc, prc, 
            σλtn, lλtn, lλup, lλac, δt, srδt, σλprior, lλ)

        llc, prc, σμc, lμtn, lμup, lμac = 
          update_σ!(σμc, Ψc, llc, prc, 
            σμtn, lμtn, lμup, lμac, δt, srδt, σμprior, lμ)


      # gbm update
      elseif pupi === 2

        # ter  = terminus[pupi]
        # inod = inodes[pupi]

        # dri = dr(idv[inod])
        # ldr = lastindex(dri)

        # llc, prc = 
        #   lvupdate!(Ψp, Ψc, llc, prc, σλc, δt, srδt, 
        #     λa_prior, dri, ldr, ter, 0)


      # forward simulation update
      else

        """
        here
        """

        tree, llc = fsp(tree, rand(idf), llc, λc, μc, ntry)

      end


      if ltn == tune_int
        σλtn = scalef(σλtn,lac/lup)
        ltn = 0
      end

    end

    next!(pbar)
  end

  return llc
end











tree, llc = fsp(tree, rand(idf), llc, λc, μc, ntry)




bi = rand(idf)

# get branch start λ and μ
dri = dr(bi)
ldr = length(dri)
λt, μt = λμi(Ψc, dri, ldr, 0)




"""
    fsbi(bi::iBf, λc::Float64, μc::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi(bi  ::iBf, 
              λc  ::Float64, 
              μc  ::Float64, 
              ntry::Int64)

  # times
  tfb = tf(bi)

  # simulate tree
  t0  = sim_gbm(ti(bi) - tfb, λt, μt, σλ, σμ, δt, srδt)

  




  ne = snen(t0)
  nt = sntn(t0)

  ret = true

  # goes extinct
  if ne === nt
    ret = false
  else
    # ntry per unobserved branch to go extinct
    for i in Base.OneTo(nt - ne - 1)
      for j in Base.OneTo(ntry)
        st0 = sim_cbd(tfb, λc, μc)
        th0 = treeheight(st0)
        # if goes extinct before the present
        if (th0 + 1e-10) < tfb
          #graft to tip
          add1(t0, st0, 1, 0)
          break
        end
        if j === ntry
          ret = false
        end
      end
    end
  end

  return t0, ret
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







