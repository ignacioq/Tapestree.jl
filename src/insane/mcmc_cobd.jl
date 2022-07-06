#=

constant occurrence birth-death MCMC using forward simulation

Jérémy Andréoletti

v(^-^v)

Created 11 02 2022
=#




"""
    insane_cobd(tree     ::sTf_label, 
                out_file ::String;
                λ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                μ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                ψ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                niter    ::Int64                 = 1_000,
                nthin    ::Int64                 = 10,
                nburn    ::Int64                 = 200,
                marginal ::Bool                  = false,
                nitpp    ::Int64                 = 100, 
                nthpp    ::Int64                 = 10,
                K        ::Int64                 = 10,
                ϵi       ::Float64               = 0.4,
                λi       ::Float64               = NaN,
                μi       ::Float64               = NaN,
                ψi       ::Float64               = NaN,
                pupdp    ::NTuple{4,Float64}     = (0.2,0.2,0.2,0.2),
                prints   ::Int64                 = 5,
                tρ       ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant fossilized birth-death.
"""
function insane_cobd(tree     ::sTf_label, 
                     out_file ::String;
                     λ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                     μ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                     ψ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                     niter    ::Int64                 = 1_000,
                     nthin    ::Int64                 = 10,
                     nburn    ::Int64                 = 200,
                     marginal ::Bool                  = false,
                     nitpp    ::Int64                 = 100, 
                     nthpp    ::Int64                 = 10,
                     K        ::Int64                 = 11,
                     ϵi       ::Float64               = 0.4,
                     λi       ::Float64               = NaN,
                     μi       ::Float64               = NaN,
                     ψi       ::Float64               = NaN,
                     pupdp    ::NTuple{4,Float64}     = (0.2,0.2,0.2,0.2),
                     prints   ::Int64                 = 5,
                     tρ       ::Dict{String, Float64} = Dict("" => 1.0))

  n  = ntips(tree)
  th   = treeheight(tree)
  stem = !iszero(e(tree))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ)

  # starting parameters
  if isnan(λi) && isnan(μi) && isnan(ψi)
    # if only one tip
    if isone(n)
      λc = prod(λ_prior)
      μc = prod(μ_prior)
    else
      λc, μc = moments(Float64(n), ti(idf[1]), ϵi)
    end
    # if no sampled fossil
    if iszero(nfossils(tree))
      ψc = prod(ψ_prior)
    else
      ψc = nfossils(tree)/treelength(tree)
    end
  else
    λc, μc, ψc = λi, μi, ψi
  end
  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 1_000, stem)

  # make a decoupled tree and fix it
  Ξ = sTfbd[]
  sTfbd!(Ξ, tree)

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant fossilized birth-death with forward simulation"

  # adaptive phase
  llc, prc, λc, μc = mcmc_burn_cobd(Ξ, idf, λ_prior, μ_prior, ψ_prior, nburn, 
                                    λc, μc, ψc, mc, th, stem, pup, prints)

  # mcmc
  r, treev, λc, μc = mcmc_cobd(Ξ, idf, llc, prc, λc, μc, ψc, λ_prior, μ_prior, 
                               ψ_prior, mc, th, stem, niter, nthin, pup, prints)

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2),
                ("psi"         => 3))

  write_ssr(r, pardic, out_file)

  if marginal

    #= # reference distribution
    βs = [range(0.0, 1.0, K)...]
    reverse!(βs)

    # make reference posterior for `λ`
    @views p = r[:,4]
    m     = mean(p)
    v     = var(p)
    λ_rdist = (m^2/v, m/v)

    # make reference posterior for `μ`
    @views p = r[:,5]
    m  = mean(p)
    sd = std(p)

    if sum(x -> x < 0.2, p) > sum(x -> 0.2 < x < 0.4, p)
      μ0 = 0.0
    else
      μ0 = m
    end

    σ0 = max(0.5, sd)

    x1 = run_newton(μ0, σ0, m, sd)

    μ_rdist = (x1[1], x1[2])

    # marginal likelihood
    pp = ref_posterior(Ξ, idf, λc, μc, v, mc, th, stem, λ_prior, μ_prior, 
                       λ_rdist, μ_rdist, nitpp, nthpp, βs, pup)

    # process with reference distribution the posterior
    p1 = Vector{Float64}(undef, size(r,1))
    for i in Base.OneTo(size(r,1))
      p1[i] = r[i,2] + r[i,3] - 
              logdgamma(r[i,4], λ_rdist[1], λ_rdist[2]) -
              logdtnorm(r[i,5], μ_rdist[1], μ_rdist[2])
    end
    pp[1] = p1

    reverse!(pp)
    reverse!(βs)

    ml = gss(pp, βs)=#
  else
    ml = NaN
  end

  return r, treev, ml
end




"""
    mcmc_burn_cobd(Ξ        ::Vector{sTfbd},
                   idf      ::Array{iBffs,1},
                   λ_prior  ::NTuple{2,Float64},
                   μ_prior  ::NTuple{2,Float64},
                   ψ_prior  ::NTuple{2,Float64},
                   nburn    ::Int64,
                   λc       ::Float64,
                   μc       ::Float64,
                   ψc       ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Bool,
                   pup      ::Array{Int64,1}, 
                   prints   ::Int64)

Adaptive MCMC phase for da chain for constant fossilized birth-death using 
forward simulation.
"""
function mcmc_burn_cobd(Ξ        ::Vector{sTfbd},
                        idf      ::Array{iBffs,1},
                        λ_prior  ::NTuple{2,Float64},
                        μ_prior  ::NTuple{2,Float64},
                        ψ_prior  ::NTuple{2,Float64},
                        nburn    ::Int64,
                        λc       ::Float64,
                        μc       ::Float64,
                        ψc       ::Float64,
                        mc       ::Float64,
                        th       ::Float64,
                        stem     ::Bool,
                        pup      ::Array{Int64,1}, 
                        prints   ::Int64)

  el = lastindex(idf)                  # number of branches
  L  = treelength(Ξ)                   # tree length
  nfos = Float64(nfossils(Ξ))          # number of fossilization events
  ns = Float64(nnodesbifurcation(Ξ))   # number of speciation events
  ne = 0.0                             # number of extinction events

  # likelihood
  llc = llik_cobd(Ξ, λc, μc, ψc) + log(mc) + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2]) + 
        logdgamma(μc, μ_prior[1], μ_prior[2]) + 
        logdgamma(ψc, ψ_prior[1], ψ_prior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1
        llc, prc, λc, mc = 
                update_λ!(llc, prc, λc, ns, L, μc, mc, th, stem, λ_prior)
  
      # μ proposal
      elseif p ===2  
        llc, prc, μc, mc = 
                update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

      # ψ proposal
      elseif p === 3
        llc, prc, ψc = update_ψ!(llc, prc, ψc, nfos, L, ψ_prior)
      
      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)
        llc, ns, ne, nfos, L = 
                  update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ns, ne, nfos, L)
      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, ψc
end




"""
    mcmc_cobd(Ξ        ::Vector{sTfbd},
              idf      ::Array{iBffs,1},
              llc      ::Float64,
              prc      ::Float64,
              λc       ::Float64,
              μc       ::Float64,
              ψc       ::Float64,
              λ_prior  ::NTuple{2,Float64},
              μ_prior  ::NTuple{2,Float64},
              ψ_prior  ::NTuple{2,Float64},
              mc       ::Float64,
              th       ::Float64,
              stem     ::Bool,
              niter    ::Int64,
              nthin    ::Int64,
              pup      ::Array{Int64,1}, 
              prints   ::Int64)

MCMC da chain for constant fossilized birth-death using forward simulation.
"""
function mcmc_cobd(Ξ      ::Vector{sTfbd},
                   idf     ::Array{iBffs,1},
                   llc     ::Float64,
                   prc     ::Float64,
                   λc      ::Float64,
                   μc      ::Float64,
                   ψc      ::Float64,
                   λ_prior ::NTuple{2,Float64},
                   μ_prior ::NTuple{2,Float64},
                   ψ_prior  ::NTuple{2,Float64},
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Bool,
                   niter   ::Int64,
                   nthin   ::Int64,
                   pup     ::Array{Int64,1}, 
                   prints  ::Int64)

  el = lastindex(idf)
  ns = Float64(nnodesbifurcation(Ξ))
  ne = Float64(ntipsextinct(Ξ))
  nfos = Float64(nfossils(Ξ))
  L  = treelength(Ξ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 6)

  # make tree vector
  treev  = sTfbd[]

  pbar = Progress(niter, prints, "running mcmc...", 0)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1
        llc, prc, λc, mc = 
                update_λ!(llc, prc, λc, ns, L, μc, mc, th, stem, λ_prior)

      # μ proposal
      elseif p === 2
        llc, prc, μc, mc = 
                update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

      # ψ proposal
      elseif p === 3
        llc, prc, ψc = update_ψ!(llc, prc, ψc, nfos, L, ψ_prior)
      
      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)
        llc, ns, ne, nfos, L = 
                  update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ns, ne, nfos, L)
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
        push!(treev, couple(deepcopy(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev, λc, μc, ψc
end




#="""
    ref_posterior(Ξ        ::Vector{sTfbd},
                  idf      ::Array{iBffs,1},
                  λc       ::Float64,
                  μc       ::Float64,
                  μtn      ::Float64,
                  ψtn      ::Float64,
                  mc       ::Float64,
                  th       ::Float64,
                  stem     ::Bool,
                  λ_prior  ::NTuple{2,Float64},
                  μ_prior  ::NTuple{2,Float64},
                  ψ_prior  ::NTuple{2,Float64},
                  λ_rdist  ::NTuple{2,Float64},
                  μ_rdist  ::NTuple{2,Float64},
                  ψ_rdist  ::NTuple{2,Float64},
                  nitpp    ::Int64,
                  nthpp    ::Int64,
                  βs       ::Vector{Float64},
                  pup      ::Array{Int64,1})

MCMC da chain for constant fossilized birth-death using forward simulation.
"""
function ref_posterior(Ξ        ::Vector{sTfbd},
                       idf      ::Array{iBffs,1},
                       λc       ::Float64,
                       μc       ::Float64,
                       μtn      ::Float64,
                       ψtn      ::Float64,
                       mc       ::Float64,
                       th       ::Float64,
                       stem     ::Bool,
                       λ_prior  ::NTuple{2,Float64},
                       μ_prior  ::NTuple{2,Float64},
                       ψ_prior  ::NTuple{2,Float64},
                       λ_rdist  ::NTuple{2,Float64},
                       μ_rdist  ::NTuple{2,Float64},
                       ψ_rdist  ::NTuple{2,Float64},
                       nitpp    ::Int64,
                       nthpp    ::Int64,
                       βs       ::Vector{Float64},
                       pup      ::Array{Int64,1})

  K = lastindex(βs)

  # make log-likelihood table per power
  nlg = fld(nitpp, nthpp)
  pp  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

  el = lastindex(idf)
  ns = Float64(nnodesinternal(Ξ))
  ne = Float64(ntipsextinct(Ξ))
  L  = treelength(Ξ)
  
  nsi = stem ? 0.0 : log(λc)

  llc = llik_cobd(Ξ, λc, μc, ψc) - nsi + log(mc) + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2]) + 
        logdgamma(μc, μ_prior[1], μ_prior[2]) + 
        logdgamma(ψc, ψ_prior[1], ψ_prior[2])

  for k in 2:K

    βi  = βs[k]
    rdc = logdgamma(λc, λ_rdist[1], λ_rdist[2]) + 
          logdtnorm(μc, μ_rdist[1], μ_rdist[2]) + 
          logdgamma(ψc, ψ_rdist[1], ψ_rdist[2])

    # logging
    lth, lit = 0, 0

    for it in Base.OneTo(nitpp)

      shuffle!(pup)

      for p in pup

        # λ proposal
        if p === 1

          llc, prc, rdc, λc, mc = 
            update_λ!(llc, prc, rdc, λc, ns, L, μc, mc, th, stem, 
              λ_prior, λ_rdist, βi)

        # forward simulation proposal proposal
        elseif p === 2 

          llc, prc, rdc, μc, mc = 
            update_μ!(llc, prc, rdc, μc, ne, L, μtn, λc, mc, th, stem, 
              μ_prior, μ_rdist, βi)

        else

          bix = ceil(Int64,rand()*el)
          llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L)

        end
      end

      # log log-likelihood
      lth += 1
      if lth === nthpp
        lit += 1
        pp[k][lit] = llc + prc - rdc
        lth = 0
      end
    end

    @info string(βi," power done")
  end

  return pp
end=#




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{sTfbd},
               idf  ::Vector{iBffs},
               llc  ::Float64,
               λ    ::Float64, 
               μ    ::Float64,
               ψ    ::Float64,
               ns   ::Float64,
               ne   ::Float64,
               nfos ::Float64,
               L    ::Float64)

Forward simulation proposal function for constant fossilized birth-death.
"""
function update_fs!(bix    ::Int64,
                    Ξ      ::Vector{sTfbd},
                    idf    ::Vector{iBffs},
                    llc    ::Float64,
                    λ      ::Float64, 
                    μ      ::Float64,
                    ψ      ::Float64,
                    ns     ::Float64,
                    ne     ::Float64,
                    nfos   ::Float64,
                    L      ::Float64)

  bi = idf[bix]

  # forward simulate an internal branch
  ξp, np, ntp = fsbi(bi, λ, μ, ψ, 100)

  # retained conditional on survival
  if ntp > 0

    itb = it(bi)   # is it terminal
    iψb = ifos(bi) # is it a fossil
    ρbi = ρi(bi)   # get branch sampling fraction
    nc  = ni(bi)   # current ni
    ntc = nt(bi)   # current nt

    # current tree
    ξc  = Ξ[bix]

    # if terminal non-fossil branch
    if itb && !iψb
      llr = log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc))
      acr = 0.0
    else
      np  -= 1
      llr = log((1.0 - ρbi)^(np - nc))
      acr = log(Float64(ntp)/Float64(ntc))
    end

    # MH ratio
    if -randexp() < llr + acr

      # update ns, ne, nfos & L
      ns +=   Float64(nnodesbifurcation(ξp) - nnodesbifurcation(ξc))
      ne +=   Float64(ntipsextinct(ξp)      - ntipsextinct(ξc))
      nfos += Float64(nfossils(ξp)          - nfossils(ξc))
      L  +=   treelength(ξp)                - treelength(ξc)

      # likelihood ratio
      llr += llik_cobd(ξp, λ, μ, ψ) - llik_cobd(ξc, λ, μ, ψ)

      Ξ[bix] = ξp     # set new decoupled tree
      llc += llr      # set new likelihood
      
      setni!(bi, np)  # set new ni
      setnt!(bi, ntp) # set new nt
    end
  end

  return llc, ns, ne, nfos, L
end





"""
    fsbi(bi::iBffs, λ::Float64, μ::Float64, ψ::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi(bi::iBffs, λ::Float64, μ::Float64, ψ::Float64, ntry::Int64)

  # times
  tfb = tf(bi)

  ext = 0
  # condition on non-extinction (helps in mixing)
  while ext < ntry

    # forward simulation during branch length
    t0, na, nfos = sim_cobd(e(bi), λ, μ, ψ, 0, 0)

    if iszero(nfos) # Exclude if any fossil is sampled
      nat = na
      
      if isone(na)
        # fix the only tip alive
        fixalive!(t0)
      
      elseif na > 1
        # fix a random tip
        fixrtip!(t0)

        if !it(bi) || ifos(bi)
          # add tips until the present
          tx, na, nfos = tip_sims!(t0, tfb, λ, μ, ψ, na, nfos)
          if !iszero(nfos) ext += 1 ; continue end
        end
      end

      if ifos(bi)
        # replace extant tip by a fossil
        fossilizefixedtip!(t0)

        # if the branch is a fossil tip, complete it with a simulated tree
        if it(bi)
          tx, na, nfos = fixedtip_sim!(t0, tfb, λ, μ, ψ, na, nfos)
          if !iszero(nfos) ext += 1 ; continue end
        end
      end
      
      return t0, na, nat
    end

    ext += 1
  end

  return sTfbd(), 0, 0
end




"""
    tip_sims!(tree::sTfbd, t::Float64, λ::Float64, μ::Float64, 
              ψ::Float64, na::Int64, nfos::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::sTfbd, t::Float64, λ::Float64, μ::Float64, 
                   ψ::Float64, na::Int64, nfos::Int64)
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if !defd1 && !defd2
    # tips
    if !isfix(tree) && isalive(tree)

      # simulate
      stree, na, nfos = sim_cobd(t, λ, μ, ψ, na-1, 0)

      if iszero(nfos)
        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    end
  else
    # bifurcations and sampled fossil ancestors
    if defd1 tree.d1, na, nfos = tip_sims!(tree.d1, t, λ, μ, ψ, na, nfos) end
    if !iszero(nfos) return tree, na, nfos end
    if defd2 tree.d2, na, nfos = tip_sims!(tree.d2, t, λ, μ, ψ, na, nfos) end
  end

  return tree, na, nfos
end




"""
    fixedtip_sim!(tree::sTfbd, t::Float64, λ::Float64, μ::Float64, 
                  ψ::Float64, na::Int64, nfos::Int64)

Continue simulation until time `t` for the fixed tip in `tree`. 
"""
function fixedtip_sim!(tree::sTfbd, t::Float64, λ::Float64, μ::Float64, 
                       ψ::Float64, na::Int64, nfos::Int64)
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  # tips
  if !defd1 && !defd2
    # simulate
    stree, na, nfos = sim_cobd(t, λ, μ, ψ, na-1, 0)
    if iszero(nfos)
      # merge to current tip
      tree.d1 = stree
    end
  
  # bifurcations and sampled fossil ancestors
  else
    if defd1 && isfix(tree.d1)
      tree.d1, na, nfos = fixedtip_sim!(tree.d1, t, λ, μ, ψ, na, nfos)
      if !iszero(nfos) return tree, na, nfos end
    end
    if defd2 && isfix(tree.d2)
      tree.d2, na, nfos = fixedtip_sim!(tree.d2, t, λ, μ, ψ, na, nfos)
    end
  end

  return tree, na, nfos
end




"""
    update_ψ!(llc    ::Float64,
              prc    ::Float64,
              ψc     ::Float64,
              nfos   ::Float64,
              L      ::Float64,
              ψ_prior::NTuple{2,Float64})

Gibbs sampling of `ψ` for constant fossilized birth-death.
"""
function update_ψ!(llc    ::Float64,
                   prc    ::Float64,
                   ψc     ::Float64,
                   nfos   ::Float64,
                   L      ::Float64,
                   ψ_prior::NTuple{2,Float64})

  ψp  = randgamma(ψ_prior[1] + nfos, ψ_prior[2] + L)

  llc += nfos * log(ψp/ψc) + L * (ψc - ψp)
  prc += llrdgamma(ψp, ψc, ψ_prior[1], ψ_prior[2])
  ψc   = ψp

  return llc, prc, ψc 
end


