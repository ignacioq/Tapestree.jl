#=

constant birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cbd(tree    ::sT_label, 
               out_file::String;
               λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               niter   ::Int64                 = 1_000,
               nthin   ::Int64                 = 10,
               nburn   ::Int64                 = 200,
               marginal::Bool                  = false,
               nitpp   ::Int64                 = 100, 
               nthpp   ::Int64                 = 10,
               K       ::Int64                 = 10,
               ϵi      ::Float64               = 0.4,
               λi      ::Float64               = NaN,
               μi      ::Float64               = NaN,
               pupdp   ::NTuple{3,Float64}     = (0.2,0.2,0.2),
               prints  ::Int64                 = 5,
               tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant birth-death.
"""
function insane_cbd(tree    ::sT_label, 
                    out_file::String;
                    λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    niter   ::Int64                 = 1_000,
                    nthin   ::Int64                 = 10,
                    nburn   ::Int64                 = 200,
                    marginal::Bool                  = false,
                    nitpp   ::Int64                 = 100, 
                    nthpp   ::Int64                 = 10,
                    K       ::Int64                 = 11,
                    ϵi      ::Float64               = 0.4,
                    λi      ::Float64               = NaN,
                    μi      ::Float64               = NaN,
                    pupdp   ::NTuple{3,Float64}     = (0.2,0.2,0.2),
                    prints  ::Int64                 = 5,
                    tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)
  stem = !iszero(e(tree))

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ)

  # starting parameters
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), ti(idf[1]), ϵi)
  else
    λc, μc = λi, μi
  end
  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 1_000, stem)

  # make a decoupled tree and fix it
  Ξ = sTbd[]
  sTbd!(Ξ, tree)

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant birth-death with forward simulation"

  # adaptive phase
  llc, prc, λc, μc, mc = 
      mcmc_burn_cbd(Ξ, idf, λ_prior, μ_prior, nburn, λc, μc, mc, th, stem,
        pup, prints)

  # mcmc
  r, treev, λc, μc, mc = mcmc_cbd(Ξ, idf, llc, prc, λc, μc, mc, th, stem,
    λ_prior, μ_prior, niter, nthin, pup, prints)

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2))

  write_ssr(r, pardic, out_file)

  if marginal

     # reference distribution
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
    pp = ref_posterior(Ξ, idf, λc, μc, v, mc, th, stem, 
      λ_prior, μ_prior, λ_rdist, μ_rdist, nitpp, nthpp, βs, pup)

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

    ml = gss(pp, βs)
  else
    ml = NaN
  end

  return r, treev, ml
end




"""
    mcmc_burn_cbd(Ξ      ::Vector{sTbd},
                  idf    ::Array{iBffs,1},
                  λ_prior::NTuple{2,Float64},
                  μ_prior::NTuple{2,Float64},
                  nburn  ::Int64,
                  λc     ::Float64,
                  μc     ::Float64,
                  mc     ::Float64,
                  th     ::Float64,
                  stem   ::Bool,
                  pup    ::Array{Int64,1}, 
                  prints ::Int64)

Adaptive MCMC phase for da chain for constant birth-death using forward
simulation.
"""
function mcmc_burn_cbd(Ξ      ::Vector{sTbd},
                       idf    ::Array{iBffs,1},
                       λ_prior::NTuple{2,Float64},
                       μ_prior::NTuple{2,Float64},
                       nburn  ::Int64,
                       λc     ::Float64,
                       μc     ::Float64,
                       mc     ::Float64,
                       th     ::Float64,
                       stem   ::Bool,
                       pup    ::Array{Int64,1}, 
                       prints ::Int64)

  el  = lastindex(idf)
  L   = treelength(Ξ)     # tree length
  nsi = stem ? 0.0 : log(λc)
  ns  = Float64(el-1)*0.5 # number of speciation events
  ne  = 0.0               # number of extinction events

  # likelihood
  llc = llik_cbd(Ξ, λc, μc) - nsi + log(mc) + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2]) + 
        logdgamma(μc, μ_prior[1], μ_prior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

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

      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L)
      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, mc
end




"""
    mcmc_cbd(Ξ      ::Vector{sTbd},
             idf    ::Array{iBffs,1},
             llc    ::Float64,
             prc    ::Float64,
             λc     ::Float64,
             μc     ::Float64,
             mc     ::Float64,
             th     ::Float64,
             stem   ::Bool,
             λ_prior::NTuple{2,Float64},
             μ_prior::NTuple{2,Float64},
             niter  ::Int64,
             nthin  ::Int64,
             pup    ::Array{Int64,1}, 
             prints ::Int64)

MCMC da chain for constant birth-death using forward simulation.
"""
function mcmc_cbd(Ξ      ::Vector{sTbd},
                  idf    ::Array{iBffs,1},
                  llc    ::Float64,
                  prc    ::Float64,
                  λc     ::Float64,
                  μc     ::Float64,
                  mc     ::Float64,
                  th     ::Float64,
                  stem   ::Bool,
                  λ_prior::NTuple{2,Float64},
                  μ_prior::NTuple{2,Float64},
                  niter  ::Int64,
                  nthin  ::Int64,
                  pup    ::Array{Int64,1}, 
                  prints ::Int64)

  el = lastindex(idf)
  ns = Float64(nnodesinternal(Ξ))
  ne = Float64(ntipsextinct(Ξ))
  L  = treelength(Ξ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 5)

  # make tree vector
  treev  = sTbd[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc = 
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, stem, λ_prior)

        # llci = llik_cbd(Ξ, λc, μc) - nsi + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc = 
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

        # llci = llik_cbd(Ξ, λc, μc) - nsi + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)
        llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L)

        # llci = llik_cbd(Ξ, λc, μc) - nsi + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
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
        push!(treev, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev, λc, μc, mc
end




"""
    ref_posterior(Ξ      ::Vector{sTbd},
                  idf    ::Array{iBffs,1},
                  λc     ::Float64,
                  μc     ::Float64,
                  μtn    ::Float64,
                  mc     ::Float64,
                  th     ::Float64,
                  stem   ::Bool,
                  λ_prior::NTuple{2,Float64},
                  μ_prior::NTuple{2,Float64},
                  λ_rdist::NTuple{2,Float64},
                  μ_rdist::NTuple{2,Float64},
                  nitpp  ::Int64,
                  nthpp  ::Int64,
                  βs     ::Vector{Float64},
                  pup    ::Array{Int64,1})

MCMC da chain for constant birth-death using forward simulation.
"""
function ref_posterior(Ξ      ::Vector{sTbd},
                       idf    ::Array{iBffs,1},
                       λc     ::Float64,
                       μc     ::Float64,
                       μtn    ::Float64,
                       mc     ::Float64,
                       th     ::Float64,
                       stem   ::Bool,
                       λ_prior::NTuple{2,Float64},
                       μ_prior::NTuple{2,Float64},
                       λ_rdist::NTuple{2,Float64},
                       μ_rdist::NTuple{2,Float64},
                       nitpp  ::Int64,
                       nthpp  ::Int64,
                       βs     ::Vector{Float64},
                       pup    ::Array{Int64,1})

  K = lastindex(βs)

  # make log-likelihood table per power
  nlg = fld(nitpp, nthpp)
  pp  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

  el = lastindex(idf)
  ns = Float64(nnodesinternal(Ξ))
  ne = Float64(ntipsextinct(Ξ))
  L  = treelength(Ξ)

  nsi = stem ? 0.0 : log(λc)

  llc = llik_cbd(Ξ, λc, μc) - nsi + log(mc) + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2]) + 
        logdgamma(μc, μ_prior[1], μ_prior[2])

  for k in 2:K

    βi  = βs[k]
    rdc = logdgamma(λc, λ_rdist[1], λ_rdist[2]) + 
          logdtnorm(μc, μ_rdist[1], μ_rdist[2])

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
end




"""
    update_fs!(bix::Int64,
               Ξ  ::Vector{sTbd},
               idf::Vector{iBffs},
               llc::Float64,
               λ  ::Float64, 
               μ  ::Float64,
               ns ::Float64,
               ne ::Float64,
               L  ::Float64)

Forward simulation proposal function for constant birth-death.
"""
function update_fs!(bix::Int64,
                    Ξ  ::Vector{sTbd},
                    idf::Vector{iBffs},
                    llc::Float64,
                    λ  ::Float64, 
                    μ  ::Float64,
                    ns ::Float64,
                    ne ::Float64,
                    L  ::Float64)

  bi = idf[bix]

  # forward simulate an internal branch
  ξp, nap, ntp = fsbi(bi, λ, μ, 100)

  # retained conditional on survival
  if ntp > 0

    itb = it(bi) # is it terminal
    ρbi = ρi(bi) # get branch sampling fraction
    nac = ni(bi) # current ni
    ntc = nt(bi) # current nt

    # current tree
    ξc  = Ξ[bix]

    # if terminal branch
    if itb
      llr = log(Float64(nap)/Float64(nac) * (1.0 - ρbi)^(nap - nac))
      acr = 0.0
    else
      nap -= 1
      llr  = log((1.0 - ρbi)^(nap - nac))
      acr  = log(Float64(ntp)/Float64(ntc))
    end

    # MH ratio
    if -randexp() < llr + acr

      # update ns, ne & L
      ns += Float64(nnodesinternal(ξp) - nnodesinternal(ξc))
      ne += Float64(ntipsextinct(ξp)   - ntipsextinct(ξc))
      L  += treelength(ξp)             - treelength(ξc)

      # likelihood ratio
      llr += llik_cbd(ξp, λ, μ) - llik_cbd(ξc, λ, μ)

      Ξ[bix] = ξp     # set new decoupled tree
      llc += llr      # set new likelihood
      setni!(bi, nap) # set new ni
      setnt!(bi, ntp) # set new nt
    end
  end

  return llc, ns, ne, L
end




"""
    fsbi(bi::iBffs, λ::Float64, μ::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi(bi::iBffs, λ::Float64, μ::Float64, ntry::Int64)

  # times
  tfb = tf(bi)

  ext = 0
  # condition on non-extinction (helps in mixing)
  while ext < ntry 
    ext += 1

    # forward simulation during branch length
    t0, na = sim_cbd(e(bi), λ, μ, 0)

    nat = na

    if isone(na)
      fixalive!(t0)

      return t0, na, nat
    elseif na > 1
      # fix random tip
      _fixrtip!(t0, na)

      if !it(bi)
        # add tips until the present
        tx, na = tip_sims!(t0, tfb, λ, μ, na)
      end

      return t0, na, nat
    end
  end

  return sTbd(0.0, false, false), 0, 0
end




"""
    tip_sims!(tree::sTbd, t::Float64, λ::Float64, μ::Float64, na::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::sTbd, t::Float64, λ::Float64, μ::Float64, na::Int64)

  if istip(tree) 
    if !isfix(tree) && isalive(tree)

      # simulate
      stree, na = sim_cbd(t, λ, μ, na-1)

      # merge to current tip
      sete!(tree, e(tree) + e(stree))
      setproperty!(tree, :iμ, isextinct(stree))
      if isdefined(stree, :d1)
        tree.d1 = stree.d1
        tree.d2 = stree.d2
      end
    end
  else
    tree.d1, na = tip_sims!(tree.d1, t, λ, μ, na)
    tree.d2, na = tip_sims!(tree.d2, t, λ, μ, na)
  end

  return tree, na
end




"""
    update_λ!(llc    ::Float64,
              prc    ::Float64,
              λc     ::Float64,
              ns     ::Float64,
              L      ::Float64,
              μc     ::Float64,
              mc     ::Float64,
              th     ::Float64,
              stem   ::Bool,
              λ_prior::NTuple{2,Float64})

Mixed HM-Gibbs sampling of `λ` for constant birth-death.
"""
function update_λ!(llc    ::Float64,
                   prc    ::Float64,
                   λc     ::Float64,
                   ns     ::Float64,
                   L      ::Float64,
                   μc     ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Bool,
                   λ_prior::NTuple{2,Float64})

  nsi = stem ? 0.0 : 1.0

  λp  = randgamma(λ_prior[1] + ns - nsi, λ_prior[2] + L)

  mp  = m_surv_cbd(th, λp, μc, 1_000, stem) 
  llr = log(mp/mc) 

  if -randexp() < llr
    llc += (ns - nsi) * log(λp/λc) + L * (λc - λp) + llr
    prc += llrdgamma(λp, λc, λ_prior[1], λ_prior[2])
    λc   = λp
    mc   = mp
  end

  return llc, prc, λc, mc
end




"""
    update_λ!(llc    ::Float64,
              prc    ::Float64,
              rdc    ::Float64,
              λc     ::Float64,
              ns     ::Float64,
              L      ::Float64,
              μc     ::Float64,
              mc     ::Float64,
              th     ::Float64,
              stem   ::Bool,
              λ_prior::NTuple{2,Float64},
              λ_rdist::NTuple{2,Float64},
              pow    ::Float64)

Mixed HM-Gibbs of `λ` for constant birth-death with reference distribution.
"""
function update_λ!(llc    ::Float64,
                   prc    ::Float64,
                   rdc    ::Float64,
                   λc     ::Float64,
                   ns     ::Float64,
                   L      ::Float64,
                   μc     ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Bool,
                   λ_prior::NTuple{2,Float64},
                   λ_rdist::NTuple{2,Float64},
                   pow    ::Float64)

  nsi = stem ? 0.0 : 1.0

  λp  = randgamma((λ_prior[1] + ns - nsi) * pow + λ_rdist[1] * (1.0 - pow),
                  (λ_prior[2] + L) * pow        + λ_rdist[2] * (1.0 - pow)) 
  mp  = m_surv_cbd(th, λp, μc, 1_000, stem) 
  llr = log(mp/mc)

  if -randexp() < (pow * llr)
    llc += (ns - nsi) * log(λp/λc) + L * (λc - λp) + llr
    prc += llrdgamma(λp, λc, λ_prior[1], λ_prior[2])
    rdc += llrdgamma(λp, λc, λ_rdist[1], λ_rdist[2])
    λc   = λp
    mc   = mp
  end

  return llc, prc, rdc, λc, mc
end




"""
    update_μ!(llc    ::Float64,
              prc    ::Float64,
              μc     ::Float64,
              ne     ::Float64,
              L      ::Float64,
              λc     ::Float64,
              mc     ::Float64,
              th     ::Float64,
              stem   ::Bool,
              μ_prior::NTuple{2,Float64})

Mixed HM-Gibbs of `μ` for constant birth-death.
"""
function update_μ!(llc    ::Float64,
                   prc    ::Float64,
                   μc     ::Float64,
                   ne     ::Float64,
                   L      ::Float64,
                   λc     ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Bool,
                   μ_prior::NTuple{2,Float64})

  μp  = randgamma(μ_prior[1] + ne, μ_prior[2] + L)

  mp  = m_surv_cbd(th, λc, μp, 1_000, stem)
  llr = log(mp/mc) 

  if -randexp() < llr
    llc += ne * log(μp/μc) + L * (μc - μp) + llr
    prc += llrdgamma(μp, μc, μ_prior[1], μ_prior[2])
    μc   = μp
    mc   = mp
  end

  return llc, prc, μc, mc
end




"""
    update_μ!(llc    ::Float64,
              prc    ::Float64,
              rdc    ::Float64,
              μc     ::Float64,
              ne     ::Float64,
              L      ::Float64,
              μtn    ::Float64,
              λc     ::Float64,
              mc     ::Float64,
              th     ::Float64,
              stem   ::Bool,
              μ_prior::NTuple{2,Float64},
              μ_rdist::NTuple{2,Float64},
              pow    ::Float64)

Mixed HM-Gibbs of `μ` for constant birth-death with reference distribution.
"""
function update_μ!(llc    ::Float64,
                   prc    ::Float64,
                   rdc    ::Float64,
                   μc     ::Float64,
                   ne     ::Float64,
                   L      ::Float64,
                   μtn    ::Float64,
                   λc     ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Bool,
                   μ_prior::NTuple{2,Float64},
                   μ_rdist::NTuple{2,Float64},
                   pow    ::Float64)

  μp  = mulupt(μc, μtn)::Float64
  mp  = m_surv_cbd(th, λc, μp, 1_000, stem)

  μr  = log(μp/μc)
  llr = ne * μr + L * (μc - μp) + log(mp/mc)
  prr = llrdgamma(μp, μc, μ_prior[1], μ_prior[2])
  rdr = llrdtnorm(μp, μc, μ_rdist[1], μ_rdist[2])

  if -randexp() < (pow * (llr + prr) + (1.0 - pow) * rdr + μr)
    llc += llr
    prc += prr
    rdc += rdr
    μc   = μp
    mc   = mp
  end

  return llc, prc, rdc, μc, mc
end


