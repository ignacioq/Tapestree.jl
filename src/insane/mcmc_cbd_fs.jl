#=

constant birth-death MCMC forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#



"""
    insane_cbd(tree    ::sTbd, 
               out_file::String;
               λprior  ::Float64           = 0.1,
               μprior  ::Float64           = 0.1,
               niter   ::Int64             = 1_000,
               nthin   ::Int64             = 10,
               nburn   ::Int64             = 200,
               tune_int::Int64             = 100,
               ϵi      ::Float64           = 0.4,
               λi      ::Float64           = NaN,
               μi      ::Float64           = NaN,
               λtni    ::Float64           = 1.0,
               μtni    ::Float64           = 1.0,
               obj_ar  ::Float64           = 0.4,
               pupdp   ::NTuple{4,Float64} = (0.4,0.4,0.1,0.1),
               prints  ::Int64              = 5)

Run insane for constant pure-birth.
"""
function insane_cbd(tree    ::sTbd, 
                    out_file::String;
                    λprior  ::Float64           = 0.1,
                    μprior  ::Float64           = 0.1,
                    niter   ::Int64             = 1_000,
                    nthin   ::Int64             = 10,
                    nburn   ::Int64             = 200,
                    tune_int::Int64             = 100,
                    ϵi      ::Float64           = 0.4,
                    λi      ::Float64           = NaN,
                    μi      ::Float64           = NaN,
                    λtni    ::Float64           = 1.0,
                    μtni    ::Float64           = 1.0,
                    obj_ar  ::Float64           = 0.4,
                    pupdp   ::NTuple{4,Float64} = (0.4,0.4,0.1,0.1),
                    prints  ::Int64              = 5)

  # tree characters
  th = treeheight(tree)
  n  = sntn(tree)

  fixtree!(tree)

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # make parameter updates scaling function for tuning
  pup = Int64[]
  for i in Base.OneTo(4) 
    append!(pup, fill(i, Int64(100.0 * pupdp[i])))
  end

  # make fix tree directory
  idf = iBf[]
  bit = BitArray{1}()
  makeiBf!(tree, idf, bit)

  # create wbf vector
  wbf = falses(lastindex(idf))

  # make survival conditioning function (stem or crown)
  svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
                           stem_prob_surv_cbd

  # adaptive phase
  llc, prc, tree, λc, μc, λtn, μtn, idf, dabr = 
      mcmc_burn_cbd(tree, n, th, tune_int, λprior, μprior, 
        nburn, ϵi, λi, μi, λtni, μtni, scalef, idf, wbf, dabr, pup, pupdp, 
        prints, svf)

  # mcmc
  R, tree = mcmc_cbd(tree, llc, prc, λc, μc, λprior, μprior,
        niter, nthin, λtn, μtn, th, idf, wbf, dabr, pup, pupdp, prints, svf)

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2), 
                ("n_extinct"   => 3),
                ("tree_length" => 4))

  write_ssr(R, pardic, out_file)

  return R, tree
end




"""
    mcmc_burn_cbd(tree    ::sTbd,
                  tl      ::Float64,
                  nt      ::Int64,
                  th      ::Float64,
                  tune_int::Int64,
                  λprior  ::Float64,
                  μprior  ::Float64,
                  nburn   ::Int64,
                  λtni    ::Float64, 
                  μtni    ::Float64, 
                  scalef  ::Function,
                  λprior  ::Float64,
                  μprior  ::Float64,
                  idf     ::Array{iBf,1},
                  wbf     ::BitArray{1},
                  dabr    ::Array{Int64,1})

adaptive MCMC phase for da chain for constant birth-death.
"""
function mcmc_burn_cbd(tree    ::sTbd,
                       n       ::Int64,
                       th      ::Float64,
                       tune_int::Int64,
                       λprior  ::Float64,
                       μprior  ::Float64,
                       nburn   ::Int64,
                       ϵi      ::Float64,
                       λi      ::Float64,
                       μi      ::Float64,
                       λtni    ::Float64, 
                       μtni    ::Float64, 
                       scalef  ::Function,
                       idf     ::Array{iBf,1},
                       wbf     ::BitArray{1},
                       dabr    ::Array{Int64,1},
                       pup     ::Array{Int64,1}, 
                       pupdp   ::NTuple{4,Float64},
                       prints  ::Int64,
                       svf     ::Function)

  # initialize acceptance log
  ltn = 0
  lup = Float64[0.0,0.0]
  lac = Float64[0.0,0.0]
  λtn = λtni
  μtn = μtni

  # starting parameters
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end

  llc = llik_cbd(tree, λc, μc) - svf(λc, μc, th)
  prc = logdexp(λc, λprior) + logdexp(μc, μprior)

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p == 1
        llc, prc, λc = λp(tree, llc, prc, λc, lac, λtn, μc, λprior, th, svf)
        lup[1] += 1.0
      end

      # μ proposal
      if p == 2
        llc, prc, μc = μp(tree, llc, prc, μc, lac, μtn, λc, μprior, th, svf)
        lup[2] += 1.0
      end
      
      # forward simulation proposal proposal

      """
      here
      """

      # log tuning parameters
      ltn += 1
      if ltn == tune_int
        λtn = scalef(λtn,lac[1]/lup[1])
        μtn = scalef(μtn,lac[2]/lup[2])
        ltn = 0
      end
    end

    next!(pbar)
  end

  return llc, prc, tree, λc, μc, λtn, μtn, idf, dabr
end





"""
    mcmc_cbd(tree  ::sTbd,
             llc   ::Float64,
             prc   ::Float64,
             λc    ::Float64,
             μc    ::Float64,
             λprior::Float64,
             μprior::Float64,
             niter ::Int64,
             nthin ::Int64,
             λtn   ::Float64,
             μtn   ::Float64, 
             th    ::Float64,
             idf   ::Array{iBf,1},
             wbf   ::BitArray{1},
             dabr  ::Array{Int64,1})

MCMC da chain for constant birth-death.
"""
function mcmc_cbd(tree  ::sTbd,
                  llc   ::Float64,
                  prc   ::Float64,
                  λc    ::Float64,
                  μc    ::Float64,
                  λprior::Float64,
                  μprior::Float64,
                  niter ::Int64,
                  nthin ::Int64,
                  λtn   ::Float64,
                  μtn   ::Float64, 
                  th    ::Float64,
                  idf   ::Array{iBf,1},
                  wbf   ::BitArray{1},
                  dabr  ::Array{Int64,1},
                  pup   ::Array{Int64,1}, 
                  pupdp ::NTuple{4,Float64},
                  prints::Int64,
                  svf   ::Function)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  R = Array{Float64,2}(undef, nlogs, 7)

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p == 1
        llc, prc, λc = λp(tree, llc, prc, λc, λtn, μc, λprior, th, svf)
      end

      # μ proposal
      if p == 2
        llc, prc, μc = μp(tree, llc, prc, μc, μtn, λc, μprior, th, svf)
      end
      
      # graft proposal
      if p == 3
        tree, llc = graftp(tree, llc, λc, μc, th, idf, wbf, dabr, pupdp)
      end
      
      # prune proposal
      if p == 4
        tree, llc = prunep(tree, llc, λc, μc, th, idf, wbf, dabr, pupdp)
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
        R[lit,4] = λc
        R[lit,5] = μc
        R[lit,6] = Float64(snen(tree))
        R[lit,7] = treelength(tree)
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, tree
end




bi = idf[1]
t0, ret = fsbi(bi, λc, μc, ntry)



# exchange t0 in branch bi



"""
    swapbranch!(tree ::sTbd,
                nbtr::sTbd,
                dri  ::BitArray{1}, 
                ldr  ::Int64,
                ix   ::Int64)

Prune tree at branch given by `dri` and grafted `wpr`.
"""
function swapbranch!(tree ::sTbd,
                     nbtr::sTbd,
                     dri  ::BitArray{1}, 
                     ldr  ::Int64,
                     ix   ::Int64)

  if ix == ldr
    addtree(nbtr, tree, 1, 0) 

    return nbtr
  elseif ix < ldr
    ix += 1
    if dri[ix]
      tree.d1 = 
        swapbranch!(tree.d1::sTbd, nbtr::sTbd, dri, ldr, ix)
    else
      tree.d2 = 
        swapbranch!(tree.d2::sTbd, nbtr::sTbd, dri, ldr, ix)
    end
  end

  return tree
end






"""
    fsp(tree::sTbd,
        llc ::Float64,
        λc  ::Float64,
        μc  ::Float64,
        th  ::Float64,
        idf ::Array{iBf,1}, 
        wbf ::BitArray{1},
        dabr::Array{Int64,1},
        pupdp::NTuple{4,Float64})

Forward simulation proposal function for constant birth-death.
"""
function fsp(tree::sTbd,
             bi  ::iBf,
             λc::Float64, 
             μc::Float64,
             ntry::Int64)

  # simulate an internal branch
  t0, ret = fsbi(bi, λc, μc, ntry)

  # if retain simulation
  if ret 

    dri = dr(bi)
    ldr = lastindex(dri)

    tree = swapbranch!(tree, nbtr, dri, ldr, 0)

    # proposal ratio
    lpr = - log(λc) 

    """
    here
    """

    # likelihood ratio
    llr = llik_cbd(t0, λc, μc) + log(2.0*λc) 


    llik_cbd(treex, λc, μc) - llik_cbd(tree, λc, μc)

    llik_cbd(t0, λc, μc) - log(λc)




    if -randexp() < lpr #+ llr
      llc += llr
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
    fsbi(bi::iBf, λc::Float64, μc::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi(bi::iBf, λc::Float64, μc::Float64, ntry::Int64)

  # times
  tfb = tf(bi)

  t0 = sim_cbd(ti(bi) - tfb, λc, μc)

  ne = snen(t0)
  nt = sntn(t0)

  ret = true

  # goes extinct
  if ne == nt
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
    addtree(tree::sTbd, dtree::sTbd, it::Int64, ix::Int64) 

Add `dtree` to not extinct tip in `tree` as speciation event.
"""
function addtree(tree::sTbd, dtree::sTbd, it::Int64, ix::Int64) 

  if istip(tree) && !isextinct(tree)
    ix += 1
    if ix === it
      tree.d1 = dtree.d1
      tree.d2 = dtree.d2
    end
    return ix 
  end

  if !isnothing(tree.d1) && ix <= it
    ix = addtree(tree.d1::sTbd, dtree::sTbd, it, ix)
  end
  if !isnothing(tree.d2) && ix <= it
    ix = addtree(tree.d2::sTbd, dtree::sTbd, it, ix)
  end

  return ix
end







"""
    add1(tree::sTbd, stree::sTbd, it::Int64, ix::Int64)

Add `stree` to tip in `tree` given by `it` in `tree.d1` order.
"""
function add1(tree::sTbd, stree::sTbd, it::Int64, ix::Int64) 

  if istip(tree) && !isextinct(tree)
    ix += 1

    if ix === it
      npe = pe(tree) + pe(stree)
      setpe!(tree, npe)
      setproperty!(tree, :iμ, stree.iμ)
      tree.d1 = stree.d1
      tree.d2 = stree.d2
    end
    return ix 
  end

  if !isnothing(tree.d1) && ix <= it
    ix = add1(tree.d1::sTbd, stree, it, ix)
  end
  if !isnothing(tree.d2) && ix <= it
    ix = add1(tree.d2::sTbd, stree, it, ix)
  end

  return ix
end




"""
    add2(tree::sTbd, stree::sTbd, it::Int64, ix::Int64) 

Add `stree` to tip in `tree` given by `it` in `tree.d2` order.
"""
function add2(tree::sTbd, stree::sTbd, it::Int64, ix::Int64) 

  if istip(tree) && !isextinct(tree)
    ix += 1

    if ix === it
      npe = pe(tree) + pe(stree)
      setpe!(tree, npe)
      setproperty!(tree, :iμ, stree.iμ)
      tree.d1 = stree.d1
      tree.d2 = stree.d2
    end
    return ix 
  end

  if !isnothing(tree.d2) && ix <= it
    ix = add2(tree.d2::sTbd, stree, it, ix)
  end
  if !isnothing(tree.d1) && ix <= it
    ix = add2(tree.d1::sTbd, stree, it, ix)
  end

  return ix
end













"""
    λp(tree  ::sTbd,
       llc   ::Float64,
       prc   ::Float64,
       λc    ::Float64,
       lac   ::Array{Float64,1},
       λtn   ::Float64,
       μc    ::Float64,
       λprior::Float64)

`λ` proposal function for constant birth-death in adaptive phase.
"""
function λp(tree  ::sTbd,
            llc   ::Float64,
            prc   ::Float64,
            λc    ::Float64,
            lac   ::Array{Float64,1},
            λtn   ::Float64,
            μc    ::Float64,
            λprior::Float64,
            th    ::Float64,
            svf   ::Function)

    λp = mulupt(λc, λtn)::Float64

    llp = llik_cbd(tree, λp, μc) - svf(λp, μc, th)
    prr = llrdexp_x(λp, λc, λprior)

    if -randexp() < (llp - llc + prr + log(λp/λc))
      llc     = llp::Float64
      prc    += prr::Float64
      λc      = λp::Float64
      lac[1] += 1.0
    end

    return llc, prc, λc
end




"""
    λp(tree  ::sTbd,
       llc   ::Float64,
       prc   ::Float64,
       λc    ::Float64,
       λtn   ::Float64,
       μc    ::Float64,
       λprior::Float64)

`λ` proposal function for constant birth-death.
"""
function λp(tree  ::sTbd,
            llc   ::Float64,
            prc   ::Float64,
            λc    ::Float64,
            λtn   ::Float64,
            μc    ::Float64,
            λprior::Float64,
            th    ::Float64,
            svf   ::Function)

    λp = mulupt(λc, λtn)::Float64

    llp = llik_cbd(tree, λp, μc) - svf(λp, μc, th)
    prr = llrdexp_x(λp, λc, λprior)

    if -randexp() < (llp - llc + prr + log(λp/λc))
      llc     = llp::Float64
      prc    += prr::Float64
      λc      = λp::Float64
    end

    return llc, prc, λc 
end




"""
    μp(tree  ::sTbd,
       llc   ::Float64,
       prc   ::Float64,
       μc    ::Float64,
       lac   ::Array{Float64,1},
       μtn   ::Float64,
       λc    ::Float64,
       μprior::Float64)

`μ` proposal function for constant birth-death in adaptive phase.
"""
function μp(tree  ::sTbd,
            llc   ::Float64,
            prc   ::Float64,
            μc    ::Float64,
            lac   ::Array{Float64,1},
            μtn   ::Float64,
            λc    ::Float64,
            μprior::Float64,
            th    ::Float64,
            svf   ::Function)

    μp = mulupt(μc, rand() < 0.3 ? μtn : 4.0*μtn)::Float64

    # one could make a ratio likelihood function
    sc  = svf(λc, μp, th)
    llp = sc == -Inf ? -Inf : llik_cbd(tree, λc, μp) - sc

    prr = llrdexp_x(μp, μc, μprior)

    if -randexp() < (llp - llc + prr + log(μp/μc))
      llc  = llp::Float64
      prc += prr::Float64
      μc   = μp::Float64
      lac[2] += 1.0
    end

    return llc, prc, μc 
end




"""
    μp(tree  ::sTbd,
       llc   ::Float64,
       prc   ::Float64,
       μc    ::Float64,
       μtn   ::Float64,
       λc    ::Float64,
       μprior::Float64)

`μ` proposal function for constant birth-death.
"""
function μp(tree  ::sTbd,
            llc   ::Float64,
            prc   ::Float64,
            μc    ::Float64,
            μtn   ::Float64,
            λc    ::Float64,
            μprior::Float64,
            th    ::Float64,
            svf   ::Function)

    μp = mulupt(μc, rand() < 0.3 ? μtn : 4.0*μtn)::Float64

    # one could make a ratio likelihood function
    sc  = svf(λc, μp, th)
    llp = sc == -Inf ? -Inf : llik_cbd(tree, λc, μp) - sc

    prr = llrdexp_x(μp, μc, μprior)

    if -randexp() < (llp - llc + prr + log(μp/μc))
      llc  = llp::Float64
      prc += prr::Float64
      μc   = μp::Float64
    end

    return llc, prc, μc 
end




"""
    λμprop()

Generate proportional proposals for `λ` and `μ`
using random samples from **LogNormal** distributions. 
"""
function λμprop() 

  lg = lnr()

  return lg*exp(randn()*0.3 - 0.044),
         lg*exp(randn()*0.3 - 0.044) 
end




"""
    lnr()

**LogNormal** random samples with median of `1`. 
"""
lnr() = @fastmath exp(randn())


