#=

constant birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cbd_fs(tree    ::sTbd, 
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
                  ntry    ::Int64             = 1,
                  λtni    ::Float64           = 1.0,
                  μtni    ::Float64           = 1.0,
                  obj_ar  ::Float64           = 0.4,
                  pupdp   ::NTuple{3,Float64} = (0.4,0.4,0.1),
                  prints  ::Int64              = 5)

Run insane for constant pure-birth.
"""
function insane_cbd_fs(tree    ::sTbd, 
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
                       ntry    ::Int64             = 2,
                       λtni    ::Float64           = 1.0,
                       μtni    ::Float64           = 1.0,
                       obj_ar  ::Float64           = 0.4,
                       pupdp   ::Array{Float64,1}  = [0.4,0.4,0.1],
                       prints  ::Int64             = 5)

  # tree characters
  th = treeheight(tree)
  n  = sntn(tree)

  fixtree!(tree)

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # make parameter updates scaling function for tuning
  pup = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, Int64(100.0 * pupdp[i])))
  end

  # make fix tree directory
  idf = iBf[]
  bit = BitArray{1}()
  makeiBf!(tree, idf, bit)

  # make survival conditioning function (stem or crown)
  svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
                           stem_prob_surv_cbd

  # adaptive phase
  llc, prc, tree, λc, μc, λtn, μtn = 
      mcmc_burn_cbd(tree, n, th, tune_int, λprior, μprior, ntry,
        nburn, ϵi, λi, μi, λtni, μtni, scalef, idf, pup, prints, svf)

  # mcmc
  R, tree = mcmc_cbd(tree, llc, prc, λc, μc, λprior, μprior, ntry,
        niter, nthin, λtn, μtn, th, idf, pup, prints, svf)

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2), 
                ("n_extinct"   => 3),
                ("tree_length" => 4))

  write_ssr(R, pardic, out_file)

  return R, tree
end




"""
    mcmc_burn_cbd(tree    ::sTbd,
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
                  pup     ::Array{Int64,1}, 
                  prints  ::Int64,
                  svf     ::Function)

Adaptive MCMC phase for da chain for constant birth-death using forward
simulation.
"""
function mcmc_burn_cbd(tree    ::sTbd,
                       n       ::Int64,
                       th      ::Float64,
                       tune_int::Int64,
                       λprior  ::Float64,
                       μprior  ::Float64,
                       ntry    ::Int64,
                       nburn   ::Int64,
                       ϵi      ::Float64,
                       λi      ::Float64,
                       μi      ::Float64,
                       λtni    ::Float64, 
                       μtni    ::Float64, 
                       scalef  ::Function,
                       idf     ::Array{iBf,1},
                       pup     ::Array{Int64,1}, 
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
      if p === 1
        llc, prc, λc = λp(tree, llc, prc, λc, lac, λtn, μc, λprior, th, svf)
        lup[1] += 1.0
      end

      # μ proposal
      if p === 2
        llc, prc, μc = μp(tree, llc, prc, μc, lac, μtn, λc, μprior, th, svf)
        lup[2] += 1.0
      end
      
      # forward simulation proposal proposal
      if p === 3
        tree, llc = fsp(tree, rand(idf), llc, λc, μc, ntry)
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

  return llc, prc, tree, λc, μc, λtn, μtn
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
             pup   ::Array{Int64,1}, 
             prints::Int64,
             svf   ::Function)

MCMC da chain for constant birth-death using forward simulation.
"""
function mcmc_cbd(tree  ::sTbd,
                  llc   ::Float64,
                  prc   ::Float64,
                  λc    ::Float64,
                  μc    ::Float64,
                  λprior::Float64,
                  μprior::Float64,
                  ntry  ::Int64,
                  niter ::Int64,
                  nthin ::Int64,
                  λtn   ::Float64,
                  μtn   ::Float64, 
                  th    ::Float64,
                  idf   ::Array{iBf,1},
                  pup   ::Array{Int64,1}, 
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

      # forward simulation proposal proposal
      if p === 3
        tree, llc = fsp(tree, rand(idf), llc, λc, μc, ntry)
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
    llik_cbd_f(tree::sTbd, λ::Float64, μ::Float64)

Estimate constant birth-death likelihood for the tree in a branch.
"""
function llik_cbd_f(tree::sTbd, λ::Float64, μ::Float64)

  if istip(tree)
    ll = - pe(tree)*(λ + μ)
  else
    ll = log(2.0*λ) - pe(tree)*(λ + μ)
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      return ll
    elseif ifx1
      ll += llik_cbd_f(tree.d1::sTbd, λ, μ) +
            llik_cbd(  tree.d2::sTbd, λ, μ)
    else
      ll += llik_cbd(  tree.d1::sTbd, λ, μ) + 
            llik_cbd_f(tree.d2::sTbd, λ, μ)
    end
  end

  return ll
end





"""
    br_ll_cbd(tree::sTbd,
              λc  ::Float64, 
              μc  ::Float64,
              dri ::BitArray{1}, 
              ldr ::Int64,
              ix  ::Int64)

Returns constant birth-death likelihood for whole branch `br`.
"""
function br_ll_cbd(tree::sTbd,
                   λc  ::Float64, 
                   μc  ::Float64,
                   dri ::BitArray{1}, 
                   ldr ::Int64,
                   ix  ::Int64)

  if ix == ldr
    return llik_cbd_f(tree, λc, μc)
  elseif ix < ldr
    ifx1 = isfix(tree.d1::sTbd)
    if ifx1 && isfix(tree.d2::sTbd)
      ix += 1
      if dri[ix]
        ll = br_ll_cbd(tree.d1::sTbd, λc, μc, dri, ldr, ix)
      else
        ll = br_ll_cbd(tree.d2::sTbd, λc, μc, dri, ldr, ix)
      end
    elseif ifx1
      ll = br_ll_cbd(tree.d1::sTbd, λc, μc, dri, ldr, ix)
    else
      ll = br_ll_cbd(tree.d2::sTbd, λc, μc, dri, ldr, ix)
    end
  end

  return ll
end





"""
    fsp(tree::sTbd,
        bi  ::iBf,
        llc ::Float64,
        λc  ::Float64, 
        μc  ::Float64,
        ntry::Int64)

Forward simulation proposal function for constant birth-death.
"""
function fsp(tree::sTbd,
             bi  ::iBf,
             llc ::Float64,
             λc  ::Float64, 
             μc  ::Float64,
             ntry::Int64)

  # forward simulate an internal branch
  t0, ret = fsbi(bi, λc, μc, ntry)

  fixalive!(t0)

  # if retain simulation
  if ret
    # get branch information
    dri = dr(bi)
    ldr = lastindex(dri)
    itb = it(bi)

    # if speciation (if branch is internal)
    iλ = itb ? 0.0 : log(2.0*λc)

    # likelihood ratio
    llr = llik_cbd(t0, λc, μc) + iλ - 
          br_ll_cbd(tree, λc, μc, dri, ldr, 0)

    if -randexp() <= 0.0
      llc += llr
      # swap branch
      tree = swapbranch!(tree, t0, dri, ldr, itb, 0)
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
    swapbranch!(tree ::sTbd,
                nbtr::sTbd,
                dri  ::BitArray{1}, 
                ldr  ::Int64,
                ix   ::Int64)

Swap branch given by `dri` by `nbtr` and return the tree.
"""
function swapbranch!(tree::sTbd,
                     nbtr::sTbd,
                     dri ::BitArray{1}, 
                     ldr ::Int64,
                     it  ::Bool,
                     ix  ::Int64)

  if ix == ldr
    if !it
      addtree(nbtr, tree) 
    end
    return nbtr
  elseif ix < ldr
    ifx1 = isfix(tree.d1::sTbd)
    if ifx1 && isfix(tree.d2::sTbd)
      ix += 1
      if dri[ix]
        tree.d1 = 
          swapbranch!(tree.d1::sTbd, nbtr::sTbd, dri, ldr, it, ix)
      else
        tree.d2 = 
          swapbranch!(tree.d2::sTbd, nbtr::sTbd, dri, ldr, it, ix)
      end
    elseif ifx1
      tree.d1 = 
          swapbranch!(tree.d1::sTbd, nbtr::sTbd, dri, ldr, it, ix)
    else
      tree.d2 = 
          swapbranch!(tree.d2::sTbd, nbtr::sTbd, dri, ldr, it, ix)
    end
  end

  return tree
end




"""
    addtree(tree::sTbd, dtree::sTbd) 

Add `dtree` to not extinct tip in `tree` as speciation event, making
sure that the daughters of `dtree` are fixed.
"""
function addtree(tree::sTbd, dtree::sTbd) 

  if istip(tree) && !isextinct(tree)

    dtree = fixds(dtree)

    tree.d1 = dtree.d1
    tree.d2 = dtree.d2

    return tree
  end

  if !isnothing(tree.d1)
    tree.d1 = addtree(tree.d1::sTbd, dtree::sTbd)
  end
  if !isnothing(tree.d2)
    tree.d2 = addtree(tree.d2::sTbd, dtree::sTbd)
  end

  return tree
end




"""
    fixds(tree::sTbd)

Returns the first tree with both daughters fixed.
"""
function fixds(tree::sTbd)

  ifx1 = isfix(tree.d1::sTbd)
  if ifx1 && isfix(tree.d2::sTbd)
    return tree
  elseif ifx1
    tree = fixds(tree.d1::sTbd)
  else
    tree = fixds(tree.d2::sTbd)
  end

  return tree
end




"""
    fixalive!(tree::sTbd)

Fixes the the path from root to the only species alive.
"""
function fixalive!(tree::sTbd)

  if istip(tree::sTbd) && !isextinct(tree::sTbd)
    fix!(tree::sTbd)
    return true
  end

  if !isnothing(tree.d2)
    f = fixalive!(tree.d2::sTbd)
    if f 
      fix!(tree)
      return true
    end
  end

  if !isnothing(tree.d1)
    f = fixalive!(tree.d1::sTbd)
    if f 
      fix!(tree)
      return true
    end
  end

  return false
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




