#=

constant birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cpe(tree    ::sT_label;
               xa      ::Dict{String, Float64};
               λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               σa_prior::NTuple{2,Float64}     = (1.0, 0.5),
               σk_prior::NTuple{2,Float64}     = (1.0, 0.5),
               niter   ::Int64                 = 1_000,
               nthin   ::Int64                 = 10,
               nburn   ::Int64                 = 200,
               nflush  ::Int64                 = nthin,
               ofile   ::String                = string(homedir(), "/pe"),
               ϵi      ::Float64               = 0.4,
               λi      ::Float64               = NaN,
               μi      ::Float64               = NaN,
               pupdp   ::NTuple{3,Float64}     = (0.2,0.2,0.2),
               prints  ::Int64                 = 5,
               survival::Bool                  = true,
               mxthf   ::Float64               = Inf,
               tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant birth-death punctuated equilibrium.
"""
function insane_cpe(tree    ::sT_label;
                    xa      ::Dict{String, Float64};
                    xs      ::Dict{String, Float64} = Dict{String,Float64}(),
                    λ_prior ::NTuple{2,Float64}      = (1.0, 1.0),
                    μ_prior ::NTuple{2,Float64}      = (1.0, 1.0),
                    σa_prior::NTuple{2,Float64}      = (0.05, 0.05),
                    σk_prior::NTuple{2,Float64}      = (0.05, 0.05),
                    niter   ::Int64                  = 1_000,
                    nthin   ::Int64                  = 10,
                    nburn   ::Int64                  = 200,
                    nflush  ::Int64                  = nthin,
                    ofile   ::String                 = string(homedir(), "/pe"),
                    ϵi      ::Float64                = 0.4,
                    λi      ::Float64                = NaN,
                    μi      ::Float64                = NaN,
                    rak     ::Float64               = 0.5,
                    pupdp   ::NTuple{4,Float64}      = (0.2, 0.2, 0.2, 0.2. 0.2),
                    prints  ::Int64                  = 5,
                    survival::Bool                   = true,
                    mxthf   ::Float64                = Inf,
                    tρ      ::Dict{String, Float64}  = Dict("" => 1.0))

  n  = ntips(tree)
  th = treeheight(tree)

  surv = 0   # condition on survival of 0, 1, or 2 starting lineages
  rmλ  = 0.0 # condition on first speciation event
  if survival 
    if iszero(e(tree)) 
      surv += 2
      rmλ  += 1.0
    else
      surv += 1
    end
  end

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # set trait uncertainty
  if isempty(xs)
    xs = Dict(tl[i] => 0.0 for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxi = make_idf(tree, tρ, xa, xs, th * mxthf)

  # starting parameters
  if isnan(λi) || isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end

  σai = σki = rak*σxi

  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 5_000, surv)

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, xr, σai, σki, sTpe)

  # get vector of internal edges
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant birth-death punctuated equilibrium"

  # adaptive phase
  llc, prc, λc, μc, mc, ns, L =
      mcmc_burn_cpe(Ξ, idf, λ_prior, μ_prior, nburn, λc, μc, mc, th, rmλ, surv,
        pup, prints)

  # mcmc
  r, treev, λc, μc, mc = mcmc_cpe(Ξ, idf, llc, prc, λc, μc, mc, ns, L, 
    th, rmλ, surv, λ_prior, μ_prior, pup, niter, nthin, nflush, ofile, prints)

  # if marginal

  #    # reference distribution
  #   βs = [range(0.0, 1.0, K)...]
  #   reverse!(βs)

  #   # make reference posterior for `λ`
  #   @views p = r[:,4]
  #   m     = mean(p)
  #   v     = var(p)
  #   λ_rdist = (m^2/v, m/v)

  #   # make reference posterior for `μ`
  #   @views p = r[:,5]
  #   m  = mean(p)
  #   sd = std(p)

  #   if sum(x -> x < 0.2, p) > sum(x -> 0.2 < x < 0.4, p)
  #     μ0 = 0.0
  #   else
  #     μ0 = m
  #   end

  #   σ0 = max(0.5, sd)

  #   x1 = run_newton(μ0, σ0, m, sd)

  #   μ_rdist = (x1[1], x1[2])

  #   # marginal likelihood
  #   pp = ref_posterior(Ξ, idf, λc, μc, v, mc, th, crown,
  #     λ_prior, μ_prior, λ_rdist, μ_rdist, nitpp, nthpp, βs, pup)

  #   # process with reference distribution the posterior
  #   p1 = Vector{Float64}(undef, size(r,1))
  #   for i in Base.OneTo(size(r,1))
  #     p1[i] = r[i,2] + r[i,3] -
  #             logdgamma(r[i,4], λ_rdist[1], λ_rdist[2]) -
  #             logdtnorm(r[i,5], μ_rdist[1], μ_rdist[2])
  #   end
  #   pp[1] = p1

  #   reverse!(pp)
  #   reverse!(βs)

  #   ml = gss(pp, βs)
  # else
  #   ml = NaN
  # end

  return r, treev
end




"""
    mcmc_burn_cpe(Ξ      ::Vector{sTpe},
                  idf    ::Array{iBffs,1},
                  λ_prior::NTuple{2,Float64},
                  μ_prior::NTuple{2,Float64},
                  nburn  ::Int64,
                  λc     ::Float64,
                  μc     ::Float64,
                  mc     ::Float64,
                  th     ::Float64,
                  rmλ    ::Float64,
                  surv   ::Int64,
                  pup    ::Array{Int64,1},
                  prints ::Int64)

Adaptive MCMC phase for da chain for constant birth-death using forward
simulation.
"""
function mcmc_burn_cpe(Ξ        ::Vector{sTpe},
                       idf     ::Array{iBffs,1},
                       λ_prior ::NTuple{2,Float64},
                       μ_prior ::NTuple{2,Float64},
                       σa_prior::NTuple{2,Float64},
                       σk_prior::NTuple{2,Float64},
                       nburn  ::Int64,
                       λc     ::Float64,
                       μc     ::Float64,
                       σac    ::Float64,
                       σkc    ::Float64,
                       mc     ::Float64,
                       th     ::Float64,
                       rmλ    ::Float64,
                       surv   ::Int64,
                       pup    ::Array{Int64,1},
                       prints ::Int64)

  el  = lastindex(idf)
  L   = treelength(Ξ)          # tree length
  ns  = nnodesbifurcation(idf) # number of speciation events
  nin = Int64(ns)              # number of internal nodes
  ne  = 0.0                    # number of extinction events

  # likelihood
  llc = llik_cpe(Ξ, idf, λc, μc, σac, σkc, ns) - rmλ * log(λc) + 
        log(mc) + prob_ρ(idf)

  # prior
  prc = logdgamma(λc, λ_prior[1], λ_prior[2])         +
        logdgamma(μc, μ_prior[1], μ_prior[2])         +
        logdinvgamma(σac^2, σa_prior[1], σa_prior[2]) +
        logdinvgamma(σkc^2, σk_prior[1], σk_prior[2])

  # tracked quantities
  sσa = [(xi(ξ) - xf(ξ))^2/e(ξ) for ξ in Ξ]
  filter!(!isnan, sσa)
  sσa = sum(sσa)
  sσk = 0.0

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, rmλ, surv, λ_prior)

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, surv, μ_prior)

      # σa (anagenetic) proposal
      elseif p === 3

        ll, prc, σac = 
          update_σ!(σac, 0.5*sσa, 2.0*ns + (1.0-rmλ), ll, prc, σa_prior)

      # σk (cladogenetic) proposal
      elseif p === 4

        ll, prc, σkc = update_σ!(σkc, 0.5*sσk, ns, ll, prc, σk_prior)

      # update inner nodes traits
      elseif p === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]
        ll, sσa, sσk = update_x!(bix, Ξ, idf, σa, σk, ll, sσa, sσk)

      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L)
      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, mc, ns, L
end




"""
    mcmc_cpe(Ξ      ::Vector{sTpe},
             idf    ::Array{iBffs,1},
             llc    ::Float64,
             prc    ::Float64,
             λc     ::Float64,
             μc     ::Float64,
             mc     ::Float64,
             ns     ::Float64,
             L      ::Float64,
             th     ::Float64,
             rmλ    ::Float64,
             surv   ::Int64,
             λ_prior::NTuple{2,Float64},
             μ_prior::NTuple{2,Float64},
             pup    ::Array{Int64,1},
             niter  ::Int64,
             nthin  ::Int64,
             nflush ::Int64,
             ofile  ::String,
             prints ::Int64)

MCMC da chain for constant birth-death using forward simulation.
"""
function mcmc_cpe(Ξ      ::Vector{sTpe},
                  idf    ::Array{iBffs,1},
                  llc    ::Float64,
                  prc    ::Float64,
                  λc     ::Float64,
                  μc     ::Float64,
                  mc     ::Float64,
                  ns     ::Float64,
                  L      ::Float64,
                  th     ::Float64,
                  rmλ    ::Float64,
                  surv   ::Int64,
                  λ_prior::NTuple{2,Float64},
                  μ_prior::NTuple{2,Float64},
                  pup    ::Array{Int64,1},
                  niter  ::Int64,
                  nthin  ::Int64,
                  nflush ::Int64,
                  ofile  ::String,
                  prints ::Int64)

  el = lastindex(idf)
  ne = Float64(ntipsextinct(Ξ))

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  r = Array{Float64,2}(undef, nlogs, 5)

  treev = sTpe[]     # make tree vector
  sthin = 0          # flush to file
  io    = IOBuffer() # buffer 

  open(ofile*".log", "w") do of
    write(of, "iteration\tlikelihood\tprior\tlambda\tmu\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      pbar = Progress(niter, prints, "running mcmc...", 20)

      for it in Base.OneTo(niter)

        shuffle!(pup)

        for p in pup

          # λ proposal
          if p === 1

            llc, prc, λc, mc =
              update_λ!(llc, prc, λc, ns, L, μc, mc, th, rmλ, surv, λ_prior)

            # llci = llik_cpe(Ξ, λc, μc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
            # if !isapprox(llci, llc, atol = 1e-6)
            #    @show llci, llc, it, p
            #    return
            # end

          # μ proposal
          elseif p === 2

            llc, prc, μc, mc =
              update_μ!(llc, prc, μc, ne, L, λc, mc, th, surv, μ_prior)

            # llci = llik_cpe(Ξ, λc, μc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
            # if !isapprox(llci, llc, atol = 1e-6)
            #    @show llci, llc, it, p
            #    return
            # end

          # forward simulation proposal proposal
          else

            bix = ceil(Int64,rand()*el)
            llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L)

            # llci = llik_cpe(Ξ, λc, μc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
            # if !isapprox(llci, llc, atol = 1e-6)
            #    @show llci, llc, it, p
            #    return
            # end
          end
        end

        # log parameters
        lthin += 1
        if lthin === nthin

          lit += 1
          @inbounds begin
            r[lit,1] = Float64(it)
            r[lit,2] = llc
            r[lit,3] = prc
            r[lit,4] = λc
            r[lit,5] = μc
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthin += 1
        if sthin === nflush
          print(of, Float64(it), '\t', llc, '\t', prc, '\t', λc,'\t', μc, '\n')
          flush(of)
          ibuffer(io, couple(Ξ, idf, 1))
          write(io, '\n')
          write(tf, take!(io))
          flush(tf)
          sthin = 0
        end

        next!(pbar)
      end
    end
  end

  return r, treev, λc, μc, mc
end




"""
    update_x!(bix ::Int64,
              Ξ   ::Vector{sTpe},
              idf ::Vector{iBffs},
              σa  ::Float64,
              σk  ::Float64,
              ll  ::Float64,
              sσa ::Float64,
              sσk ::Float64)

Perform a punkeep trait update for an internal branch and its descendants.
"""
function update_x!(bix ::Int64,
                   Ξ   ::Vector{sTpe},
                   idf ::Vector{iBffs},
                   σa  ::Float64,
                   σk  ::Float64,
                   ll  ::Float64,
                   sσa ::Float64,
                   sσk ::Float64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  b1   = idf[i1]
  i2   = d2(bi)
  ξ1   = Ξ[i1]
  root = iszero(pa(bi))

  # if mrca
  if root && iszero(e(ξi))
    # if crown
    ll, sσa, sσk = _crown_update!(ξi, ξ1, ξ2, σa, σk, ll, sσa, sσk)
  else
  # if stem
    if root
      ll, sσa = _stem_update!(ξi, σa, ll, sσa)
    end

    # updates within the parent branch
    ll, sσa, sσk = _update_node_x!(ξi, σa, σk, ll, sσa, sσk)

    # get fixed tip
    lξi = fixtip(ξi)

    isd = iszero(i2)
    # if duo
    if isd
      ll, sσa = _update_duo_x!(lξi, ξ1, σa, ll, sσa)
    # if triad
    else
      ξ2  = Ξ[i2]
      ll, sσa, sσk = _update_node_x!(ξi, ξ1, ξ2, σa, σk, ll, sσa, sσk)
    end

    ### update daughters
    ## D1
    # if leaf
    if iszero(d1(b1))
      if ifx(b1) 
          ll, sσa, sσk = 
            _update_leaf_x!(ξ1, xavg(b1), xstd(b1), σa, σk, ll, sσa, sσk)
      else
          ll, sσa, sσk = _update_leaf_x!(ξ1, σa, σk, ll, sσa, sσk)
      end
    # if not leaf
    else
      ll, sσa, sσk = _update_node_x!(ξ1, σa, σk, ll, sσa, sσk)
    end

    if !isd
      ## D2
      b2 = idf[i2]
      # if leaf
      if iszero(d1(b2))
        ξ2 = Ξ[i2]
        if ifx(b2)
          ll, sσa, sσk = 
            _update_leaf_x!(ξ2, xavg(b2), xstd(b2), σa, σk, ll, sσa, sσk)
        else
          ll, sσa, sσk = _update_leaf_x!(ξ2, σa, σk, ll, sσa, sσk)
        end
      # if not leaf
      else
        ll, sσa, sσk = _update_node_x!(ξ2, σa, σk, ll, sσa, sσk)
      end
    end

  return ll, sσa, sσk
end




"""
    update_fs!(bix::Int64,
               Ξ  ::Vector{sTpe},
               idf::Vector{iBffs},
               llc::Float64,
               λ  ::Float64,
               μ  ::Float64,
               σa ::Float64,
               σk ::Float64,
               ns ::Float64,
               ne ::Float64,
               L  ::Float64,
               sσa::Float64, 
               sσk::Float64)

Forward simulation proposal function for constant punkeek.
"""
function update_fs!(bix::Int64,
                    Ξ  ::Vector{sTpe},
                    idf::Vector{iBffs},
                    llc::Float64,
                    λ  ::Float64,
                    μ  ::Float64,
                    σa ::Float64,
                    σk ::Float64,
                    ns ::Float64,
                    ne ::Float64,
                    L  ::Float64,
                    sσa::Float64, 
                    sσk::Float64,
                    xis::Vector{Float64},
                    xfs::Vector{Float64},
                    es ::Vector{Float64})

  bi = idf[bix]
  ξc = Ξ[bix]

   # if terminal branch
  if iszero(d1(bi))

    xav = xsd = NaN
    if !isnothing(xavg(bi))
      xav, xsd = xavg(bi), xstd(bi)
    end

    ξp, llr = fsbi_t(bi, xav, xsd, ξc, λ, μ, σa, σk, xis, xfs, es)

  # if mid branch
  elseif iszero(d2(bi))

    """
    here
    """

    ξp, llr, sdXr = fsbi_m(bi, Ξ[d1(bi)], xi(ξc), λ, μ, σx)


  # if trio branch
  else

    ξp, llr, sdXr = fsbi_i(bi, Ξ[d1(bi)], Ξ[d2(bi)], xi(ξc), λ, μ, σx)
  end

  if isfinite(llr)
    ξc  = Ξ[bix]

    # update llc, ns, ne & L

    # Make one function that joins these tree traversals
    llik_cpe_trackers(ξp, λ, μ,)...

    llc += llik_cpe(ξp, λ, μ) - llik_cpe(ξc, λ, μ) + llr
    ns  += Float64(nnodesinternal(ξp) - nnodesinternal(ξc))
    ne  += Float64(ntipsextinct(ξp)   - ntipsextinct(ξc))
    L   += treelength(ξp)             - treelength(ξc)

    sσa +=
    sσk +=

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, ns, ne, L
end




"""
    fsbi_t(bi  ::iBffs,
           xav::Float64,
           xst::Float64,
           ξc  ::sTpe,
           λ   ::Float64,
           μ   ::Float64,
           σa  ::Float64,
           σk  ::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi  ::iBffs,
                xav::Float64,
                xst::Float64,
                ξc  ::sTpe,
                λ   ::Float64,
                μ   ::Float64,
                σa  ::Float64,
                σk  ::Float64,
                xis::Vector{Float64},
                xfs::Vector{Float64},
                es ::Vector{Float64})

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inverse branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  empty!(xis)
  empty!(xfs)
  empty!(es)

  t0, na, nn, llr =
    _sim_cpe_t(e(bi), λ, μ, xi(ξc), σa, σk, lc, lU, iρi, 0, 1, 500, 
               xis, xfs, es)

  if na < 1 || isnan(llr)
    return t0, NaN
  end

  # if fix node
  if ifx(bi)

    # sample tip
    wti = fIrand(na) + 1
    lξc = fixtip(ξc)

    xp = NaN
    if iszero(xst)
      xp  = xavg
      acr = logdnorm(xavg, xis[wti], es[wti]*σa^2) - 
            logdnorm(xavg, xi(lξc),     e(lξc)*σa^2)
    else
      xp  = xfs[wti]
      acr = duoldnorm(     xp, xis[wti], xavg, es[wti]*σa^2, xst) - 
            duoldnorm(xf(lξc),  xi(lξc), xavg,  e(lξc)*σa^2, xst)
    end

    if lU < acr + llr

      if wti <= div(na,2)
        fixtip1!(t0, wti, 0, xp)
      else
        fixtip2!(t0, na - wti + 1, 0, xp)
      end

      setni!(bi, na)    # set new ni
      return t0, llr
    end

  else
    if lU < llr
      _fixrtip!(t0, na)

      setni!(bi, na)    # set new ni
      return t0, llr
    end
  end

  return t0, NaN
end




"""
    fsbi_m(bi::iBffs,
           ξ1::sTpe,
           λ ::Float64,
           μ ::Float64,
           x0::Float64,
           σa::Float64,
           σk::Float64)

Forward simulation for internal branch.
"""
function fsbi_m(bi::iBffs,
                ξ1::sTpe,
                λ ::Float64,
                μ ::Float64,
                x0::Float64,
                σa::Float64,
                σk::Float64)

  t0, na, nn = _sim_cpe_i(e(bi), λ, μ, x0, σa, σk, 0, 1, 500)

  if na < 1 || nn > 999
    return t0, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # inverse branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # get fix `x`
  xp = fixrtip!(t0, na, NaN)

  # acceptance ration with respect to daughters
  llr = llrdnorm_x(xp, xi(ξ1), xf(ξ1), e(ξ1)*σa^2)

  if lU < acr + llr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr = tip_sims!(t0, tf(bi), λ, μ, acr, lU, iρi, na, nn)
    end

    if lU < acr
      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
      setnt!(bi, ntp)                # set new nt
      setni!(bi, na)                 # set new ni

      return t0, llr
    end
  end

  return t0, NaN
end



"""
    fsbi_i(bi::iBffs, λ::Float64, μ::Float64)

Forward simulation for internal branch.
"""
function fsbi_i(bi::iBffs, λ::Float64, μ::Float64)

  t0, na, nn = _sim_cpe_i(e(bi), λ, μ, 0, 1, 1_000)

  if na < 1 || nn > 999
    return t0, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  if lU < acr

    _fixrtip!(t0, na) # fix random tip

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr = tip_sims!(t0, tf(bi), λ, μ, acr, lU, iρi, na, nn)
    end

    if lU < acr
      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
      setnt!(bi, ntp)                # set new nt
      setni!(bi, na)                 # set new ni

      return t0, llr
    end
  end

  return t0, NaN
end




"""
    tip_sims!(tree::sTpe,
              t   ::Float64,
              λ   ::Float64,
              μ   ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTpe,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn ::Int64)

  if lU < lr && nn < 1_000

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nn, lr = 
          _sim_cpe_it(t, λ, μ, lr, lU, iρi, na-1, nn, 1_000)

        if isnan(lr) || nn > 999
          return tree, na, nn, NaN
        end

        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr = tip_sims!(tree.d1, t, λ, μ, lr, lU, iρi, na, nn)
      tree.d2, na, nn, lr = tip_sims!(tree.d2, t, λ, μ, lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end





