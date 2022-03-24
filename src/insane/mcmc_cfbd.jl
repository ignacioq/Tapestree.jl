#=

constant fossilized birth-death MCMC using forward simulation

Jérémy Andréoletti
Adapted from birth-death MCMC by Ignacio Quintero Mächler

v(°-°v)

Created 07 10 2021
=#




"""
    insane_cfbd(tree    ::sTf_label,
                out_file::String;
                λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
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
                ψi      ::Float64               = NaN,
                pupdp   ::NTuple{4,Float64}     = (0.2,0.2,0.2,0.2),
                prints  ::Int64                 = 5,
                tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant fossilized birth-death.
"""
function insane_cfbd(tree    ::sTf_label,
                     out_file::String;
                     λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
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
                     ψi      ::Float64               = NaN,
                     pupdp   ::NTuple{4,Float64}     = (0.2,0.2,0.2,0.2),
                     prints  ::Int64                 = 5,
                     tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
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
      λc, μc = moments(Float64(n), th, ϵi)
    end
    # if no sampled fossil
    if iszero(nfossils(tree))
      ψc = prod(ψ_prior)
    else
      ψc = Float64(nfossils(tree))/Float64(treelength(tree))
    end
  else
    λc, μc, ψc = λi, μi, ψi
  end

  # define conditioning
  if ntipsalive(tree) > 0
    # if crown conditioning
    if def1(tree) && def2(tree) &&
       ntipsalive(tree.d1) > 0 && ntipsalive(tree.d2) > 0
      stem = 1
    # if stem conditioning
    else
      stem = 0
    end
  # no survival
  else
    stem = 2
  end

  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 500, stem)

  # make a decoupled tree and fix it
  Ξ = sTfbd[]
  sTfbd!(Ξ, tree)

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant fossilized birth-death"

  # adaptive phase
  llc, prc, λc, μc, ψc, mc =
     mcmc_burn_cfbd(Ξ, idf, λ_prior, μ_prior, ψ_prior, nburn,
        λc, μc, ψc, mc, th, stem, pup, prints)

  # mcmc
  r, treev, λc, μc, ψc =
    mcmc_cfbd(Ξ, idf, llc, prc, λc, μc, ψc, λ_prior, μ_prior,
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
    mcmc_burn_cfbd(Ξ      ::Vector{sTfbd},
                   idf    ::Array{iBffs,1},
                   λ_prior::NTuple{2,Float64},
                   μ_prior::NTuple{2,Float64},
                   ψ_prior::NTuple{2,Float64},
                   nburn  ::Int64,
                   λc     ::Float64,
                   μc     ::Float64,
                   ψc     ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Bool,
                   pup    ::Array{Int64,1},
                   prints ::Int64)

Adaptive MCMC phase for da chain for constant fossilized birth-death using
forward simulation.
"""
function mcmc_burn_cfbd(Ξ      ::Vector{sTfbd},
                        idf    ::Array{iBffs,1},
                        λ_prior::NTuple{2,Float64},
                        μ_prior::NTuple{2,Float64},
                        ψ_prior::NTuple{2,Float64},
                        nburn  ::Int64,
                        λc     ::Float64,
                        μc     ::Float64,
                        ψc     ::Float64,
                        mc     ::Float64,
                        th     ::Float64,
                        stem   ::Int64,
                        pup    ::Array{Int64,1},
                        prints ::Int64)

  el = lastindex(idf)                # number of branches
  L  = treelength(Ξ)                 # tree length
  nf = Float64(nfossils(Ξ))          # number of fossilization events
  ns = Float64(nnodesbifurcation(Ξ)) # number of speciation events
  ne = Float64(ntipsextinct(Ξ))      # number of extinction events

  # likelihood
  llc = llik_cfbd(Ξ, λc, μc, ψc) - stem*log(λc) + log(mc) + prob_ρ(idf)
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
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

      # ψ proposal
      elseif p === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)
        llc, ns, ne, L =
          update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ns, ne, L)

      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, ψc, mc
end




"""
    mcmc_cfbd(Ξ      ::Vector{sTfbd},
              idf    ::Array{iBffs,1},
              llc    ::Float64,
              prc    ::Float64,
              λc     ::Float64,
              μc     ::Float64,
              ψc     ::Float64,
              λ_prior::NTuple{2,Float64},
              μ_prior::NTuple{2,Float64},
              ψ_prior::NTuple{2,Float64},
              mc     ::Float64,
              th     ::Float64,
              stem   ::Bool,
              niter  ::Int64,
              nthin  ::Int64,
              pup    ::Array{Int64,1},
              prints ::Int64)

MCMC da chain for constant fossilized birth-death using forward simulation.
"""
function mcmc_cfbd(Ξ      ::Vector{sTfbd},
                   idf    ::Array{iBffs,1},
                   llc    ::Float64,
                   prc    ::Float64,
                   λc     ::Float64,
                   μc     ::Float64,
                   ψc     ::Float64,
                   λ_prior::NTuple{2,Float64},
                   μ_prior::NTuple{2,Float64},
                   ψ_prior::NTuple{2,Float64},
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Int64,
                   niter  ::Int64,
                   nthin  ::Int64,
                   pup    ::Array{Int64,1},
                   prints ::Int64)

  el = lastindex(idf)
  ns = Float64(nnodesbifurcation(Ξ))
  ne = Float64(ntipsextinct(Ξ))
  nf = Float64(nfossils(Ξ))
  L  = treelength(Ξ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 6)

  # make tree vector
  treev  = sTfbd[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, stem, λ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc) - stem * log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc) - stem * log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # ψ proposal
      elseif p === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc) - stem * log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)
        llc, ns, ne, L =
          update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ns, ne, L)

        # llci = llik_cfbd(Ξ, λc, μc, ψc) - stem * log(λc) + log(mc) + prob_ρ(idf)
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
        R[lit,6] = ψc
        push!(treev, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev, λc, μc, ψc
end




"""
    update_fs!(bix::Int64,
               Ξ  ::Vector{sTfbd},
               idf::Vector{iBffs},
               llc::Float64,
               λ  ::Float64,
               μ  ::Float64,
               ψ  ::Float64,
               ns ::Float64,
               ne ::Float64,
               L  ::Float64)

Forward simulation proposal function for constant fossilized birth-death.
"""
function update_fs!(bix::Int64,
                    Ξ  ::Vector{sTfbd},
                    idf::Vector{iBffs},
                    llc::Float64,
                    λ  ::Float64,
                    μ  ::Float64,
                    ψ  ::Float64,
                    ns ::Float64,
                    ne ::Float64,
                    L  ::Float64)

  bi = idf[bix]

  # if terminal alive
  if it(bi) && !isfossil(bi)
    ξp, llr = fsbi_t(bi, λ, μ, ψ)
  else
    ξp, llr = fsbi_i(bi, λ, μ, ψ)
  end

  if isfinite(llr)
    ξc  = Ξ[bix]

    # update llc, ns, ne & L
    llc += llik_cfbd(ξp, λ, μ, ψ)        - llik_cfbd(ξc, λ, μ, ψ) + llr
    ns  += Float64(nnodesbifurcation(ξp) - nnodesbifurcation(ξc))
    ne  += Float64(ntipsextinct(ξp)      - ntipsextinct(ξc))
    L   += treelength(ξp)                - treelength(ξc)

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, ns, ne, L
end





"""
    fsbi_t(bi::iBffs, λ::Float64, μ::Float64, ψ::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs, λ::Float64, μ::Float64, ψ::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_cfbd_t(e(bi), λ, μ, ψ, lc, lU, Iρi, 0, 1, 500)

  if na > 0 && isfinite(llr) 

    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, -Inf
  end
end




"""
    fsbi_i(bi::iBffs, λ::Float64, μ::Float64, ψ::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi::iBffs, λ::Float64, μ::Float64, ψ::Float64)


  # forward simulation during branch length
  t0, na, nf, nn = _sim_cfbd_i(e(bi), λ, μ, ψ, 0, 0, 1, 500)

  if na < 1 || nf > 0 || nn >= 500
    return t0, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if lU < acr

    _fixrtip!(t0, na)

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, acr, lU, Iρi, na, nn)
    end

    if lU < acr

      if isfossil(bi)
        # fossilize extant tip
        fossilizefixedtip!(t0)

        # if terminal fossil branch
        if it(bi)
          tx, na, nn, acr = 
            fossiltip_sim!(t0, tf(bi), λ, μ, ψ, acr, lU, Iρi, na, nn)
        end
      end

      if lU < acr
        na -= !it(bi)
        llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
        setnt!(bi, ntp)                # set new nt
        setni!(bi, na)                 # set new ni

        return t0, llr
      end
    end
  end

  return t0, NaN
end




"""
    tip_sims!(tree::sTfbd,
              t   ::Float64,
              λ   ::Float64,
              μ   ::Float64,
              ψ   ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nf  ::Int64
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTfbd,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

  if lU < lr && nn < 500

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nn, lr = 
          _sim_cfbd_it(t, λ, μ, ψ, lr, lU, Iρi, na-1, nn, 500)

        if isnan(lr) || nn >= 500
          return tree, na, nn, NaN
        end

        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        if def1(stree)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr = 
        tip_sims!(tree.d1, t, λ, μ, ψ, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, λ, μ, ψ, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::sTfbd,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

Continue simulation until time `t` for the fixed tip in `tree`.
"""
function fossiltip_sim!(tree::sTfbd,
                        t   ::Float64,
                        λ   ::Float64,
                        μ   ::Float64,
                        ψ   ::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        Iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64)

  if lU < lr && nn < 500
    if istip(tree)
      stree, na, nn, lr = 
        _sim_cfbd_it(t, λ, μ, ψ, lr, lU, Iρi, na-1, nn, 500)

      if isnan(lr) || nn >= 500
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr = 
        fossiltip_sim!(tree.d1, t, λ, μ, ψ, lr, lU, Iρi, na, nn)
    else
      tree.d2, na, nn, lr = 
        fossiltip_sim!(tree.d2, t, λ, μ, ψ, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_ψ!(llc    ::Float64,
              prc    ::Float64,
              ψc     ::Float64,
              nf     ::Float64,
              L      ::Float64,
              ψ_prior::NTuple{2,Float64})

Gibbs sampling of `ψ` for constant fossilized birth-death.
"""
function update_ψ!(llc    ::Float64,
                   prc    ::Float64,
                   ψc     ::Float64,
                   nf     ::Float64,
                   L      ::Float64,
                   ψ_prior::NTuple{2,Float64})

  ψp  = randgamma(ψ_prior[1] + nf, ψ_prior[2] + L)

  llc += nf * log(ψp/ψc) + L * (ψc - ψp)
  prc += llrdgamma(ψp, ψc, ψ_prior[1], ψ_prior[2])

  return llc, prc, ψp
end



