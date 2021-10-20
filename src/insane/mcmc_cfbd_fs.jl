#=

constant fossilized birth-death MCMC using forward simulation

Jérémy Andréoletti
Adapted from birth-death MCMC by Ignacio Quintero Mächler

v(°-°v)

Created 07 10 2021
=#




"""
    insane_cfbd_fs(tree    ::sTfbd, 
                   out_file::String,
                   λprior  ::Float64,
                   μprior  ::Float64,
                   ψprior  ::Float64,
                   niter   ::Int64,
                   nthin   ::Int64,
                   nburn   ::Int64,
                   tune_int::Int64,
                   ϵi      ::Float64,
                   λi      ::Float64,
                   μi      ::Float64,
                   ψi      ::Float64,
                   λtni    ::Float64,
                   μtni    ::Float64,
                   ψtni    ::Float64,
                   obj_ar  ::Float64,
                   pupdp   ::NTuple{4,Float64},
                   prints  ::Int64)

Run insane for constant pure-birth.
"""
function insane_cfbd_fs(tree    ::sTfbd, 
                        out_file::String,
                        λprior  ::Float64,
                        μprior  ::Float64,
                        ψprior  ::Float64,
                        niter   ::Int64,
                        nthin   ::Int64,
                        nburn   ::Int64,
                        tune_int::Int64,
                        ϵi      ::Float64,
                        λi      ::Float64,
                        μi      ::Float64,
                        ψi      ::Float64,
                        λtni    ::Float64,
                        μtni    ::Float64,
                        ψtni    ::Float64,
                        obj_ar  ::Float64,
                        pupdp   ::NTuple{4,Float64},
                        prints  ::Int64)

  # tree characters
  th = treeheight(tree)
  n  = ntips(tree)

  fixtree!(tree)

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # make fix tree directory
  idf = iBfffs[]
  bit = BitArray{1}()
  makeiBf!(tree, idf, bit)

  # make survival conditioning function (stem or crown) and identify branches
  # from `idf`
  # if iszero(e(tree))
  #    svf = crown_prob_surv_cfbd
  #    cb = findall(x -> isone(lastindex(dr(x))), idf)
  # else
  #    svf = stem_prob_surv_cfbd
  #    cb = Int64[findfirst(x -> iszero(lastindex(dr(x))), idf)]
  # end

  if iszero(e(tree))
     svf = cond_surv_crown
     cb = findall(x -> isone(lastindex(dr(x))), idf)
  else
     svf = cond_surv_stem
     cb = Int64[findfirst(x -> iszero(lastindex(dr(x))), idf)]
  end

  @info "Running constant fossilized birth-death with forward simulation"
  @info "Burnin adaptive phase"

  # adaptive phase
  llc, prc, tree, λc, μc, ψc, λtn, μtn, ψtn = 
      mcmc_burn_cfbd(tree, n, th, tune_int, λprior, μprior, ψprior,nburn, ϵi, 
              λi, μi, ψi, λtni, μtni, ψtni, scalef, idf, pup, prints, svf, cb)

  @info "Final MCMC"
  # mcmc
  R, tree = mcmc_cfbd(tree, llc, prc, λc, μc, ψc, λprior, μprior, ψprior,
              niter, nthin, λtn, μtn, ψtn, th, idf, pup, prints, svf, cb)

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2),
                ("psi"         => 3), 
                ("n_extinct"   => 4),
                ("tree_length" => 5))

  println("size(R)=$(size(R))")

  write_ssr(R, pardic, out_file)

  return R, tree
end





"""
    mcmc_burn_cfbd(tree    ::sTfbd,
                   n       ::Int64,
                   th      ::Float64,
                   tune_int::Int64,
                   λprior  ::Float64,
                   μprior  ::Float64,
                   ψprior  ::Float64,
                   nburn   ::Int64,
                   ϵi      ::Float64,
                   λi      ::Float64,
                   μi      ::Float64,
                   ψi      ::Float64,
                   λtni    ::Float64, 
                   μtni    ::Float64, 
                   ψtni    ::Float64, 
                   scalef  ::Function,
                   idf     ::Array{iBfffs,1},
                   pup     ::Array{Int64,1}, 
                   prints  ::Int64,
                   svf     ::Function)

Adaptive MCMC phase for da chain for constant fossilized birth-death using forward
simulation.
"""
function mcmc_burn_cfbd(tree    ::sTfbd,
                        n       ::Int64,
                        th      ::Float64,
                        tune_int::Int64,
                        λprior  ::Float64,
                        μprior  ::Float64,
                        ψprior  ::Float64,
                        nburn   ::Int64,
                        ϵi      ::Float64,
                        λi      ::Float64,
                        μi      ::Float64,
                        ψi      ::Float64, 
                        λtni    ::Float64, 
                        μtni    ::Float64, 
                        ψtni    ::Float64,
                        scalef  ::Function,
                        idf     ::Array{iBfffs,1},
                        pup     ::Array{Int64,1}, 
                        prints  ::Int64,
                        svf     ::Function,
                        cb      ::Array{Int64,1})

  # initialize acceptance log
  ltn = 0
  lup = Float64[0.0,0.0]
  lac = Float64[0.0,0.0]
  λtn = λtni
  μtn = μtni
  ψtn = ψtni

  # starting parameters
  if isnan(λi) && isnan(μi) && isnan(ψi)
    λc, μc = moments(Float64(n), th, ϵi)
    ψc = nfossils(tree)/treelength(tree)
  else
    λc, μc, ψc = λi, μi, ψi
  end

  # length idf
  lidf = lastindex(idf)

  # add simulated subtrees to all fossil tips
  for bi in filter(x -> it(x) && ifos(x), idf)
    dri = dr(bi)
    ldr = lastindex(dri)
    ret = false
    while !ret
      t0, ret = fsψtip(bi, λc, μc, ψc)
    end
    @show t0, e(t0)
    swapfossil!(tree, t0, dri, ldr, 0)
  end

  # likelihood
  llc = llik_cfbd(tree, λc, μc, ψc) + svf(tree, λc, μc)
  # llc = llik_cfbd(tree, λc, μc, ψc) - svf(λc, μc, th)
  prc = logdexp(λc, λprior) + logdexp(μc, μprior) + logdexp(ψc, ψprior)

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1
        llc, prc, λc = λp(tree, llc, prc, λc, lac, λtn, μc, ψc, λprior, th, svf)
        lup[1] += 1.0
      end

      # μ proposal
      if p === 2
        llc, prc, μc = μp(tree, llc, prc, μc, lac, μtn, λc, ψc, μprior, th, svf)
        lup[2] += 1.0
      end

      # ψ proposal
      if p === 3
        llc, prc, ψc = μp(tree, llc, prc, ψc, lac, μtn, λc, μc, ψprior, th, svf)
        lup[2] += 1.0
      end
      
      # forward simulation proposal proposal
      if p === 4
        bix = fIrand(lidf) + 1
        tree, llc = fsp(tree, idf[bix], llc, λc, μc, ψc, in(bix, cb))
      end

      # log tuning parameters
      ltn += 1
      if ltn == tune_int
        λtn = scalef(λtn,lac[1]/lup[1])
        μtn = scalef(μtn,lac[2]/lup[2])
        ψtn = scalef(ψtn,lac[2]/lup[2])
        ltn = 0
      end
    end

    next!(pbar)
  end

  return llc, prc, tree, λc, μc, ψc, λtn, μtn, ψtn
end





"""
    mcmc_cfbd(tree  ::sTfbd,
              llc   ::Float64,
              prc   ::Float64,
              λc    ::Float64,
              μc    ::Float64,
              ψc    ::Float64,
              λprior::Float64,
              μprior::Float64,
              ψprior::Float64,
              niter ::Int64,
              nthin ::Int64,
              λtn   ::Float64,
              μtn   ::Float64, 
              ψtn   ::Float64,
              th    ::Float64,
              idf   ::Array{iBfffs,1},
              pup   ::Array{Int64,1}, 
              prints::Int64,
              svf   ::Function)

MCMC da chain for constant fossilized birth-death using forward simulation.
"""
function mcmc_cfbd(tree  ::sTfbd,
                   llc   ::Float64,
                   prc   ::Float64,
                   λc    ::Float64,
                   μc    ::Float64,
                   ψc    ::Float64,
                   λprior::Float64,
                   μprior::Float64,
                   ψprior::Float64,
                   niter ::Int64,
                   nthin ::Int64,
                   λtn   ::Float64,
                   μtn   ::Float64, 
                   ψtn   ::Float64,
                   th    ::Float64,
                   idf   ::Array{iBfffs,1},
                   pup   ::Array{Int64,1}, 
                   prints::Int64,
                   svf   ::Function,
                   cb    ::Array{Int64,1})

  lidf = lastindex(idf)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 8)

  # make Ψ vector
  treev = sTfbd[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1
        llc, prc, λc = λp(tree, llc, prc, λc, λtn, μc, ψc, λprior, th, svf)

        # llci = llik_cfbd(tree, λc, μc, ψc) + svf(tree, λc, μc)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, i, 1
        #    return 
        # end
      end

      # μ proposal
      if p === 2
        llc, prc, μc = μp(tree, llc, prc, μc, μtn, λc, ψc, μprior, th, svf)
      
        # llci = llik_cfbd(tree, λc, μc, ψc) + svf(tree, λc, μc)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, i, 2
        #    return 
        # end
      end

      # ψ proposal
      if p === 3
        llc, prc, ψc = ψp(tree, llc, prc, ψc, ψtn, λc, μc, ψprior, th, svf)
      
        # llci = llik_cfbd(tree, λc, μc, ψc) + svf(tree, λc, μc)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, i, 2
        #    return 
        # end
      end

      # forward simulation proposal
      if p === 4
        bix = fIrand(lidf) + 1
        tree, llc = fsp(tree, idf[bix], llc, λc, μc, ψc, in(bix, cb))
      
        # llci = llik_cfbd(tree, λc, μc, ψc) + svf(tree, λc, μc)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, i, 3
        #    return 
        # end
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
        R[lit,6] = ψc
        R[lit,7] = Float64(ntipsextinct(tree))
        R[lit,8] = treelength(tree)
        push!(treev, deepcopy(tree))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev
end




"""
    llik_cfbd_f(tree::sTfbd, λ::Float64, μ::Float64)

Estimate constant fossilized birth-death likelihood for the tree in a branch.
"""
function llik_cfbd_f(tree::sTfbd, λ::Float64, μ::Float64, ψ::Float64)

  ll = - e(tree)*(λ + μ + ψ)

  if issampledancestor(tree)
    # end of the fixed branch
    return ll + log(ψ)
  end

  # birth
  if isdefined(tree, :d1) && isdefined(tree, :d2)
    ll += log(λ)
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      # end of the fixed branch
      return ll
    elseif ifx1
      ll += llik_cfbd_f(tree.d1::sTfbd, λ, μ, ψ) +
            llik_cfbd(  tree.d2::sTfbd, λ, μ, ψ)
    else
      ll += llik_cfbd(  tree.d1::sTfbd, λ, μ, ψ) + 
            llik_cfbd_f(tree.d2::sTfbd, λ, μ, ψ)
    end
  end

  return ll
end




"""
    br_ll_cfbd(tree::sTfbd,
               λ   ::Float64, 
               μ   ::Float64,
               ψ   ::Float64,
               dri ::BitArray{1}, 
               ldr ::Int64,
               ix  ::Int64)

Returns constant fossilized birth-death likelihood for whole branch at `dri`.
"""
function br_ll_cfbd(tree::sTfbd,
                    λ   ::Float64, 
                    μ   ::Float64,
                    ψ   ::Float64,
                    dri ::BitArray{1}, 
                    ldr ::Int64,
                    ix  ::Int64)
  if ix === ldr
    return llik_cfbd_f(tree, λ, μ, ψ)
  end
  
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  # Sampled ancestors
  if defd1 && !defd2 
    ix += 1
    return br_ll_cfbd(tree.d1::sTfbd, λ, μ, ψ, dri, ldr, ix)
  end
  if defd2 && !defd1
    ix += 1
    return br_ll_cfbd(tree.d2::sTfbd, λ, μ, ψ, dri, ldr, ix)
  end

  if ix < ldr
    ifx1 = isfix(tree.d1::sTfbd)
    if ifx1 && isfix(tree.d2::sTfbd)
      ix += 1
      if dri[ix]
        ll = br_ll_cfbd(tree.d1::sTfbd, λ, μ, ψ, dri, ldr, ix)
      else
        ll = br_ll_cfbd(tree.d2::sTfbd, λ, μ, ψ, dri, ldr, ix)
      end
    elseif ifx1
      ll = br_ll_cfbd(tree.d1::sTfbd, λ, μ, ψ, dri, ldr, ix)
    else
      ll = br_ll_cfbd(tree.d2::sTfbd, λ, μ, ψ, dri, ldr, ix)
    end
  end

  return ll
end




"""
    tipψ_ll_cfbd(tree::sTfbd,
                 λ   ::Float64, 
                 μ   ::Float64,
                 ψ   ::Float64,
                 dri ::BitArray{1}, 
                 ldr ::Int64,
                 ix  ::Int64)

Returns constant fossilized birth-death likelihood for augmented subtree 
following the observed fossil tip at `dri`.
"""
function tipψ_ll_cfbd(tree::sTfbd,
                      λ   ::Float64, 
                      μ   ::Float64,
                      ψ   ::Float64,
                      dri ::BitArray{1}, 
                      ldr ::Int64,
                      ix  ::Int64)
  # branch reached
  if ix === ldr
    # end of branch reached
    if isfossil(tree)
      return llik_cfbd(tree.d1, λ, μ, ψ)
    else
      if isfix(tree.d1::sTfbd)
        return br_ll_cfbd(tree.d1::sTfbd, λ, μ, ψ, dri, ldr, ix)
      else
        return br_ll_cfbd(tree.d2::sTfbd, λ, μ, ψ, dri, ldr, ix)
      end
    end
    
  end
  
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  # Sampled ancestors
  if defd1 && !defd2 
    ix += 1
    return br_ll_cfbd(tree.d1::sTfbd, λ, μ, ψ, dri, ldr, ix)
  end
  if defd2 && !defd1
    ix += 1
    return br_ll_cfbd(tree.d2::sTfbd, λ, μ, ψ, dri, ldr, ix)
  end

  if ix < ldr
    ifx1 = isfix(tree.d1::sTfbd)
    if ifx1 && isfix(tree.d2::sTfbd)
      ix += 1
      if dri[ix]
        ll = br_ll_cfbd(tree.d1::sTfbd, λ, μ, ψ, dri, ldr, ix)
      else
        ll = br_ll_cfbd(tree.d2::sTfbd, λ, μ, ψ, dri, ldr, ix)
      end
    elseif ifx1
      ll = br_ll_cfbd(tree.d1::sTfbd, λ, μ, ψ, dri, ldr, ix)
    else
      ll = br_ll_cfbd(tree.d2::sTfbd, λ, μ, ψ, dri, ldr, ix)
    end
  end

  return ll
end




"""
    fsp(tree::sTfbd,
        bi  ::iBfffs,
        llc ::Float64,
        λ   ::Float64, 
        μ   ::Float64,
        ψ   ::Float64,
        scb ::Bool)

Forward simulation proposal function for constant fossilized birth-death.
"""
function fsp(tree::sTfbd,
             bi  ::iBfffs,
             llc ::Float64,
             λ   ::Float64, 
             μ   ::Float64,
             ψ   ::Float64,
             scb ::Bool)

  # forward simulate an internal branch
  t0, ret = fsbi(bi, λ, μ, ψ)

  # if retain simulation
  if ret

    dri = dr(bi)
    ldr = lastindex(dri)
    itb = it(bi)
    iψb = ifos(bi)

    # if speciation (if branch is internal)
    iλ = itb||iψb ? 0.0 : log(λ)
    # if fossil sampling (if branch is internal)
    iψ = iψb ? 0.0 : log(ψ)

    ## likelihood ratio
    # if conditioning branch
    if scb 
      # if one of crown branches
      if isone(lastindex(dri))
        if dri[1]
          llr = llik_cfbd(   t0, λ, μ, ψ) + iλ + iψ    - 
                br_ll_cfbd(tree, λ, μ, ψ, dri, ldr, 0) +
                cond_surv_stem_p(t0,    λ, μ)      - 
                cond_surv_stem(tree.d1, λ, μ)
        else
          llr = llik_cfbd(   t0, λ, μ, ψ) + iλ + iψ    -
                br_ll_cfbd(tree, λ, μ, ψ, dri, ldr, 0) +
                cond_surv_stem_p(t0,    λ, μ)      -
                cond_surv_stem(tree.d2, λ, μ)
        end
      else
        llr = llik_cfbd(t0,    λ, μ, ψ) + iλ + iψ    -
              br_ll_cfbd(tree, λ, μ, ψ, dri, ldr, 0) +
              cond_surv_stem_p(t0, λ, μ)         -
              cond_surv_stem(tree, λ, μ)
      end
    else
      llr = llik_cfbd(t0, λ, μ, ψ) + iλ + iψ - 
            br_ll_cfbd(tree, λ, μ, ψ, dri, ldr, 0)
    end

    llc += llr

    # swap branch
    tree = swapbranch!(tree, t0, dri, ldr, itb, iψb, 0)
  end

  # if branch ending with fossil tip: proposal for following the simulation
  if it(bi) && ifos(bi)
    # forward simulate the prolongation of a fossil tip
    t0, ret = fsψtip(bi, λ, μ, ψ)
    if ret
      dri = dr(bi)
      ldr = lastindex(dri)
      llr = llik_cfbd(t0, λ, μ, ψ) - tipψ_ll_cfbd(tree, λ, μ, ψ, dri, ldr, 0)
      llc += llr

      # swap fossil tip
      swapfossil!(tree, t0, dri, ldr, 0)
    end
  end

  return tree, llc
end




"""
    fsbi(bi::iBfffs, λ::Float64, μ::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi(bi::iBfffs, λ::Float64, μ::Float64, ψ::Float64)

  # retain?
  ret = true

  # times
  tfb = tf(bi)

  # forward simulation during branch length
  t0 = sim_cfbd(ti(bi) - tfb, λ, μ, ψ)
  na = ntipsalive(t0)

  if iszero(na) || nfossils(t0)>0
    ret = false
  elseif isone(na)
    fixalive!(t0)
  elseif na > 1
    if it(bi)
      ret = false
    else
      # fix random tip
      _fixrtip!(t0, na)

      for j in Base.OneTo(na - 1)
        for i in Base.OneTo(2)
          st0 = sim_cfbd(tfb, λ, μ, ψ)
          # if there is any sampled fossil
          if nfossils(st0)>0
            ret = false
            break
          end
          # if goes extinct before the present
          if iszero(ntipsalive(st0))
            #graft to tip
            addtotip(t0, st0, false)
            break
          end
          if i === 2
            ret = false
          end
        end
        !ret && break
      end
    end
  end

  return t0, ret
end




"""
    fsψtip(bi::iBfffs, λ::Float64, μ::Float64, ψ::Float64)

Forward simulation for a fossil tip at branch `bi`
"""
function fsψtip(bi::iBfffs, λ::Float64, μ::Float64, ψ::Float64)

  # retain?
  ret = true

  # forward simulation during branch length
  t0 = sim_cfbd(tf(bi), λ, μ, ψ)

  # do not retain if any extant or fossil sampling
  if ntipsalive(t0)>0 || nfossils(t0)>0
    ret = false
  end

  return t0, ret
end




"""
    addtotip(tree::sTfbd, stree::sTfbd, ix::Bool) 

Add `stree` to tip in `tree` given by `it` in `tree.d1` order.
"""
function addtotip(tree::sTfbd, stree::sTfbd, ix::Bool) 

  if istip(tree)
    if isalive(tree) && !isfix(tree)

      sete!(tree, e(tree) + e(stree))
      setproperty!(tree, :iμ, isextinct(stree))

      if isdefined(stree, :d1)
        tree.d1 = stree.d1
        tree.d2 = stree.d2
      end

      ix = true
    end

    return ix 
  end

  if !ix
    ix = addtotip(tree.d1::sTfbd, stree, ix)
  end
  if !ix
    ix = addtotip(tree.d2::sTfbd, stree, ix)
  end

  return ix
end




"""
    λp(tree  ::sTfbd,
       llc   ::Float64,
       prc   ::Float64,
       λc    ::Float64,
       lac   ::Array{Float64,1},
       λtn   ::Float64,
       μc    ::Float64,
       ψc    ::Float64,
       λprior::Float64)

`λ` proposal function for constant fossilized birth-death in adaptive phase.
"""
function λp(tree  ::sTfbd,
            llc   ::Float64,
            prc   ::Float64,
            λc    ::Float64,
            lac   ::Array{Float64,1},
            λtn   ::Float64,
            μc    ::Float64,
            ψc    ::Float64,
            λprior::Float64,
            th    ::Float64,
            svf   ::Function)

    λp = mulupt(λc, λtn)::Float64

    llp = llik_cfbd(tree, λp, μc, ψc) + svf(tree, λp, μc)
    # llp = llik_cfbd(tree, λp, μc) - svf(λp, μc, th)

    prr = llrdexp_x(λp, λc, λprior)

    #println("λc=$(round(λc; digits=2)), λtn=$(round(λtn; digits=2)), λp=$(round(λp; digits=2)), ll=$(round(llp - llc + prr + log(λp/λc)))")

    if -randexp() < (llp - llc + prr + log(λp/λc))
      llc     = llp::Float64
      prc    += prr::Float64
      λc      = λp::Float64
      lac[1] += 1.0
    end

    return llc, prc, λc
end




"""
    λp(tree  ::sTfbd,
       llc   ::Float64,
       prc   ::Float64,
       λc    ::Float64,
       λtn   ::Float64,
       μc    ::Float64,
       ψc    ::Float64,
       λprior::Float64)

`λ` proposal function for constant fossilized birth-death.
"""
function λp(tree  ::sTfbd,
            llc   ::Float64,
            prc   ::Float64,
            λc    ::Float64,
            λtn   ::Float64,
            μc    ::Float64,
            ψc    ::Float64,
            λprior::Float64,
            th    ::Float64,
            svf   ::Function)
  
    λp = mulupt(λc, rand() < 0.3 ? λtn : 4.0*λtn)::Float64

    llp = llik_cfbd(tree, λp, μc, ψc) + svf(tree, λp, μc)
    # llp = llik_cfbd(tree, λp, μc) - svf(λp, μc, th)

    prr = llrdexp_x(λp, λc, λprior)

    #println("λc=$(round(λc; digits=2)), λtn=$(round(λtn; digits=2)), λp=$(round(λp; digits=2)), ll=$(round(llp - llc + prr + log(λp/λc)))")

    if -randexp() < (llp - llc + prr + log(λp/λc))
      llc     = llp::Float64
      prc    += prr::Float64
      λc      = λp::Float64
    end

    return llc, prc, λc 
end




"""
    μp(tree  ::sTfbd,
       llc   ::Float64,
       prc   ::Float64,
       μc    ::Float64,
       lac   ::Array{Float64,1},
       μtn   ::Float64,
       λc    ::Float64,
       ψc    ::Float64,
       μprior::Float64)

`μ` proposal function for constant fossilized birth-death in adaptive phase.
"""
function μp(tree  ::sTfbd,
            llc   ::Float64,
            prc   ::Float64,
            μc    ::Float64,
            lac   ::Array{Float64,1},
            μtn   ::Float64,
            λc    ::Float64,
            ψc    ::Float64,
            μprior::Float64,
            th    ::Float64,
            svf   ::Function)

    μp = mulupt(μc, μtn)::Float64

    # one could make a ratio likelihood function
    sc  = svf(tree, λc, μp)
    llp = isinf(sc) ? -Inf : llik_cfbd(tree, λc, μp, ψc) + sc
    # sc  = svf(λc, μp, th)
    # llp = isinf(sc) ? -Inf : llik_cfbd(tree, λc, μp, ψc) - sc

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
    μp(tree  ::sTfbd,
       llc   ::Float64,
       prc   ::Float64,
       μc    ::Float64,
       μtn   ::Float64,
       λc    ::Float64,
       ψc    ::Float64,
       μprior::Float64)

`μ` proposal function for constant fossilized birth-death.
"""
function μp(tree  ::sTfbd,
            llc   ::Float64,
            prc   ::Float64,
            μc    ::Float64,
            μtn   ::Float64,
            λc    ::Float64,
            ψc    ::Float64,
            μprior::Float64,
            th    ::Float64,
            svf   ::Function)

    μp = mulupt(μc, rand() < 0.3 ? μtn : 4.0*μtn)::Float64

    # one could make a ratio likelihood function
    sc  = svf(tree, λc, μp)
    llp = isinf(sc) ? -Inf : llik_cfbd(tree, λc, μp, ψc) + sc
    # sc  = svf(λc, μp, th)
    # llp = isinf(sc) ? -Inf : llik_cfbd(tree, λc, μp, ψc) - sc

    prr = llrdexp_x(μp, μc, μprior)

    if -randexp() < (llp - llc + prr + log(μp/μc))
      llc  = llp::Float64
      prc += prr::Float64
      μc   = μp::Float64
    end

    return llc, prc, μc 
end




"""
    ψp(tree  ::sTfbd,
       llc   ::Float64,
       prc   ::Float64,
       ψc    ::Float64,
       lac   ::Array{Float64,1},
       ψtn   ::Float64,
       λc    ::Float64,
       μc    ::Float64,
       ψprior::Float64,
       th    ::Float64,
       svf   ::Function)

`ψ` proposal function for constant fossilized birth-death in adaptive phase.
"""
function ψp(tree  ::sTfbd,
            llc   ::Float64,
            prc   ::Float64,
            ψc    ::Float64,
            lac   ::Array{Float64,1},
            ψtn   ::Float64,
            λc    ::Float64,
            μc    ::Float64,
            ψprior::Float64,
            th    ::Float64,
            svf   ::Function)

    ψp = mulupt(ψc, ψtn)::Float64

    # one could make a ratio likelihood function
    sc  = svf(tree, λc, μc)
    llp = isinf(sc) ? -Inf : llik_cfbd(tree, λc, μc, ψp) + sc
    # sc  = svf(λc, μc, th)
    # llp = isinf(sc) ? -Inf : llik_cfbd(tree, λc, μc, ψp) - sc

    prr = llrdexp_x(ψp, ψc, ψprior)

    if -randexp() < (llp - llc + prr + log(ψp/ψc))
      llc  = llp::Float64
      prc += prr::Float64
      ψc   = ψp::Float64
      lac[2] += 1.0
    end

    return llc, prc, ψc 
end




"""
    ψp(tree  ::sTfbd,
       llc   ::Float64,
       prc   ::Float64,
       ψc    ::Float64,
       ψtn   ::Float64,
       λc    ::Float64,
       μc    ::Float64,
       ψprior::Float64,
       th    ::Float64,
       svf   ::Function)

`ψ` proposal function for constant fossilized birth-death.
"""
function ψp(tree  ::sTfbd,
            llc   ::Float64,
            prc   ::Float64,
            ψc    ::Float64,
            ψtn   ::Float64,
            λc    ::Float64,
            μc    ::Float64,
            ψprior::Float64,
            th    ::Float64,
            svf   ::Function)

    ψp = mulupt(ψc, rand() < 0.3 ? ψtn : 4.0*ψtn)::Float64

    # one could make a ratio likelihood function
    sc  = svf(tree, λc, μc)
    llp = isinf(sc) ? -Inf : llik_cfbd(tree, λc, μc, ψp) + sc
    # sc  = svf(λc, μc, th)
    # llp = isinf(sc) ? -Inf : llik_cfbd(tree, λc, μc, ψp) - sc

    prr = llrdexp_x(ψp, ψc, ψprior)

    if -randexp() < (llp - llc + prr + log(ψp/ψc))
      llc  = llp::Float64
      prc += prr::Float64
      ψc   = ψp::Float64
    end

    return llc, prc, ψc 
end




"""
    λμψprop()

Generate proportional proposals for `λ`, `μ` and `ψ`
using random samples from **LogNormal** distributions. 
"""
function λμψprop() 

  lg = lnr()

  return lg*exp(randn()*0.3 - 0.044),
         lg*exp(randn()*0.3 - 0.044),
         lg*exp(randn()*0.3 - 0.044)
end




"""
    lnr()

**LogNormal** random samples with median of `1`. 
"""
lnr() = @fastmath exp(randn())





