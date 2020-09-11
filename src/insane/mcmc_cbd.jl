#=

birth-death MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cbd(tree    ::iTbd, 
               out_file::String;
               λprior  ::Float64 = 0.1,
               μprior  ::Float64 = 0.1,
               niter   ::Int64   = 1_000,
               nthin   ::Int64   = 10,
               nburn   ::Int64   = 200,
               tune_int::Int64   = 100,
               λtni    ::Float64 = 1.0,
               μtni    ::Float64 = 1.0,
               obj_ar  ::Float64 = 0.4)

Run insane for constant pure-birth.
"""
function insane_cbd(tree    ::iTbd, 
                    out_file::String;
                    λprior  ::Float64           = 0.1,
                    μprior  ::Float64           = 0.1,
                    niter   ::Int64             = 1_000,
                    nthin   ::Int64             = 10,
                    nburn   ::Int64             = 200,
                    tune_int::Int64             = 100,
                    ϵi      ::Float64           = 0.4,
                    λi      ::Float64           = 100.,
                    μi      ::Float64           = 100.,
                    λtni    ::Float64           = 1.0,
                    μtni    ::Float64           = 1.0,
                    obj_ar  ::Float64           = 0.4,
                    pupdp   ::NTuple{4,Float64} = (0.4,0.4,0.1,0.1),
                    prints  ::Int64              = 5)

  # tree characters
  th = treeheight(tree)
  nt = sntn(tree)

  fixtree!(tree)

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # make parameter updates scaling function for tuning
  pup = Int64[]
  for i in Base.OneTo(4) 
    append!(pup, fill(i, Int64(100.0 * pupdp[i])))
  end

  # make fix tree directory
  idv = iDir[]
  bit = BitArray{1}()
  makeiDir!(tree, idv, bit)

  # create wbr vector
  wbr = falses(lastindex(idv))

  # create da branches vector
  dabr = Int64[]

  # make survival conditioning function (stem or crown)
  svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
                           stem_prob_surv_cbd

  # adaptive phase
  llc, prc, tree, λc, μc, λtn, μtn, idv, dabr = 
      mcmc_burn_cbd(tree, nt, th, tune_int, λprior, μprior, 
        nburn, ϵi, λi, μi, λtni, μtni, scalef, idv, wbr, dabr, pup, pupdp, 
        prints, svf)

  # mcmc
  R, tree = mcmc_cbd(tree, llc, prc, λc, μc, λprior, μprior,
        niter, nthin, λtn, μtn, th, idv, wbr, dabr, pup, pupdp, prints, svf)

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2), 
                ("n_extinct"   => 3),
                ("tree_length" => 4))

  write_ssr(R, pardic, out_file)

  return R, tree
end




"""
    mcmc_burn_cbd(tree    ::iTbd,
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
                  idv     ::Array{iDir,1},
                  wbr     ::BitArray{1},
                  dabr    ::Array{Int64,1})

MCMC da chain for constant birth-death.
"""
function mcmc_burn_cbd(tree    ::iTbd,
                       nt      ::Int64,
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
                       idv     ::Array{iDir,1},
                       wbr     ::BitArray{1},
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

  # starting parameters (using method of moments)
  if λi == 100.0 && μi == 100.0
    δ   = 1.0/th * log(Float64(nt)*(1.0 - ϵi) + ϵi)
    λc  = δ/(1.0 - ϵi)
    μc  = λc - δ
  else
    λc = λi
    μc = μi
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
      
      # graft proposal
      if p == 3
        tree, llc = graftp(tree, llc, λc, μc, th, idv, wbr, dabr, pupdp)
      end
      
      # prune proposal
      if p == 4
        tree, llc = prunep(tree, llc, λc, μc, th, idv, wbr, dabr, pupdp)
      end

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

  return llc, prc, tree, λc, μc, λtn, μtn, idv, dabr
end





"""
    mcmc_cbd(tree  ::iTbd,
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
             idv   ::Array{iDir,1},
             wbr   ::BitArray{1},
             dabr  ::Array{Int64,1})

MCMC da chain for constant birth-death.
"""
function mcmc_cbd(tree  ::iTbd,
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
                  idv   ::Array{iDir,1},
                  wbr   ::BitArray{1},
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
        tree, llc = graftp(tree, llc, λc, μc, th, idv, wbr, dabr, pupdp)
      end
      
      # prune proposal
      if p == 4
        tree, llc = prunep(tree, llc, λc, μc, th, idv, wbr, dabr, pupdp)
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




"""
    graftp(tree::iTbd,
           llc ::Float64,
           λc  ::Float64,
           μc  ::Float64,
           th  ::Float64,
           idv ::Array{iDir,1}, 
           wbr ::BitArray{1},
           dabr::Array{Int64,1},
           pupdp::NTuple{4,Float64})

Graft proposal function for constant birth-death.
"""
function graftp(tree::iTbd,
                llc ::Float64,
                λc  ::Float64,
                μc  ::Float64,
                th  ::Float64,
                idv ::Array{iDir,1}, 
                wbr ::BitArray{1},
                dabr::Array{Int64,1},
                pupdp::NTuple{4,Float64})

  #λn, μn = λμprop() 

  #simulate extinct lineage
  t0, t0h = sim_cbd_b(λc, μc, th, 100)

      # if useful simulation
  if t0h < th

    # randomly select branch to graft
    h, br, bri, nbh  = randbranch(th, t0h, idv, wbr)

    # proposal ratio
    lpr = log(2.0 * μc * (th - t0h) * Float64(nbh) * pupdp[4]) - 
          log((Float64(lastindex(dabr)) + 1.0) * pupdp[3])

    # likelihood ratio
    llr = llik_cbd(t0, λc, μc) + log(2.0*λc) 

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
    prunep(tree::iTbd,
           llc ::Float64,
           λc  ::Float64,
           μc  ::Float64,
           th  ::Float64,
           idv ::Array{iDir,1}, 
           wbr ::BitArray{1},
           dabr::Array{Int64,1},
           pupdp::NTuple{4,Float64})

Prune proposal function for constant birth-death.
"""
function prunep(tree::iTbd,
                llc ::Float64,
                λc  ::Float64,
                μc  ::Float64,
                th  ::Float64,
                idv ::Array{iDir,1}, 
                wbr ::BitArray{1},
                dabr::Array{Int64,1},
                pupdp::NTuple{4,Float64})

  ng = lastindex(dabr)

  if ng > 0
    dabri = rand(Base.OneTo(ng))
    br    = idv[dabr[dabri]]
    dri   = dr(br)
    ldr   = lastindex(dri)
    wpr   = rand(Base.OneTo(da(br)))

    # return tree height
    h, th0 = streeheight(tree, th, 0.0, dri, ldr, wpr, 0, 1)

    # get how many branches are cut at `h`
    nbh = branchescut!(wbr, h, idv)

    # proposal ratio
    lpr = log(Float64(ng) * pupdp[3]) -
          log(2.0 * μc * (th - th0) * Float64(nbh) * pupdp[4])

    # likelihood ratio
    llr = - stree_ll_cbd(tree, 0.0, λc, μc, dri, ldr, wpr, 0, 1) - 
            log(2.0 * λc)

    if -randexp() < lpr #+ llr
      llc += llr
      tree = prunetree!(tree, dri, ldr, wpr, 0, 1)
      # remove graft from branch
      # add graft to branch
      rmda!(br)
      # log branch as being data augmented
      deleteat!(dabr, dabri)
    end
  end

  return tree, llc
end





"""
    λp(tree  ::iTbd,
       llc   ::Float64,
       prc   ::Float64,
       λc    ::Float64,
       lac   ::Array{Float64,1},
       λtn   ::Float64,
       μc    ::Float64,
       λprior::Float64)

`λ` proposal function for constant birth-death in adaptive phase.
"""
function λp(tree  ::iTbd,
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
    λp(tree  ::iTbd,
       llc   ::Float64,
       prc   ::Float64,
       λc    ::Float64,
       λtn   ::Float64,
       μc    ::Float64,
       λprior::Float64)

`λ` proposal function for constant birth-death.
"""
function λp(tree  ::iTbd,
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
    μp(tree  ::iTbd,
       llc   ::Float64,
       prc   ::Float64,
       μc    ::Float64,
       lac   ::Array{Float64,1},
       μtn   ::Float64,
       λc    ::Float64,
       μprior::Float64)

`μ` proposal function for constant birth-death in adaptive phase.
"""
function μp(tree  ::iTbd,
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
    μp(tree  ::iTbd,
       llc   ::Float64,
       prc   ::Float64,
       μc    ::Float64,
       μtn   ::Float64,
       λc    ::Float64,
       μprior::Float64)

`μ` proposal function for constant birth-death.
"""
function μp(tree  ::iTbd,
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









