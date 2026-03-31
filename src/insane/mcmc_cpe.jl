#=

constant birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cpe(tree    ::sT_label,
               xa      ::Dict{String, Float64};
               xs      ::Dict{String, Float64} = Dict{String,Float64}(),
               λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               σa_prior::NTuple{2,Float64}     = (0.05, 0.05),
               σk_prior::NTuple{2,Float64}     = (0.05, 0.05),
               niter   ::Int64                 = 1_000,
               nthin   ::Int64                 = 10,
               nburn   ::Int64                 = 200,
               nflush  ::Int64                 = nthin,
               ofile   ::String                = string(homedir(), "/cpe"),
               ϵi      ::Float64               = 0.4,
               λi      ::Float64               = NaN,
               μi      ::Float64               = NaN,
               pupdp   ::NTuple{6,Float64}     = (1e-2, 1e-2, 1e-2, 1e-2, 0.1, 0.2),
               prints  ::Int64                 = 5,
               survival::Bool                  = true,
               mxthf   ::Float64               = 0.1,
               tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant birth-death punctuated equilibrium.
"""
function insane_cpe(tree    ::sT_label,
                    xa      ::Dict{String, Float64};
                    xs      ::Dict{String, Float64} = Dict{String,Float64}(),
                    λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    σa_prior::NTuple{2,Float64}     = (0.05, 0.05),
                    σk_prior::NTuple{2,Float64}     = (0.05, 0.05),
                    niter   ::Int64                 = 1_000,
                    nthin   ::Int64                 = 10,
                    nburn   ::Int64                 = 200,
                    nflush  ::Int64                 = nthin,
                    ofile   ::String                = string(homedir(), "/cpe"),
                    ϵi      ::Float64               = 0.4,
                    λi      ::Float64               = NaN,
                    μi      ::Float64               = NaN,
                    pupdp   ::NTuple{6,Float64}     = (1e-2, 1e-2, 1e-2, 1e-2, 0.1, 0.2),
                    prints  ::Int64                 = 5,
                    survival::Bool                  = true,
                    mxthf   ::Float64               = 0.1,
                    tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n  = ntips(tree)
  th = treeheight(tree)

  surv = 0   # condition on survival of 0, 1, or 2 starting lineages
  rmλ  = 0.0 # condition on first speciation event
  if iszero(e(tree)) 
    rmλ  += 1.0
    surv += survival ? 2 : 0
  else
    surv += survival ? 1 : 0
  end

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxi = make_idf(tree, tρ, xa, xs, th * mxthf)

  # starting parameters
  λc, μc = λi, μi
  if isnan(λi) || isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  end

  σac = σkc = σxi

  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 5_000, surv)

  # get vector of internal edges
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(6)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running constant punctuated equilibrium"

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, xr, σkc, sTpe)

  # adaptive phase
  llc, prc, λc, μc, σac, σkc, mc, ns, ne, L, sσa, sσk, nσs =
      mcmc_burn_cpe(Ξ, idf, λ_prior, μ_prior, σa_prior, σk_prior, nburn, 
        λc, μc, σac, σkc, mc, th, rmλ, inodes, surv, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cpe(Ξ, idf, llc, prc, λc, μc, σac, σkc, mc, ns, ne, L, sσa, sσk, nσs,
      th, rmλ, inodes, surv, λ_prior, μ_prior, σa_prior, σk_prior, pup, 
      niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_cpe(Ξ        ::Vector{sTpe},
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
                  inodes ::Vector{Int64},
                  surv   ::Int64,
                  pup    ::Array{Int64,1},
                  prints ::Int64)

Burn-in for constant birth-death punctuated equilibrium.
"""
function mcmc_burn_cpe(Ξ       ::Vector{sTpe},
                       idf     ::Array{iBffs,1},
                       λ_prior ::NTuple{2,Float64},
                       μ_prior ::NTuple{2,Float64},
                       σa_prior::NTuple{2,Float64},
                       σk_prior::NTuple{2,Float64},
                       nburn   ::Int64,
                       λc      ::Float64,
                       μc      ::Float64,
                       σac     ::Float64,
                       σkc     ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       rmλ     ::Float64,
                       inodes  ::Vector{Int64},
                       surv    ::Int64,
                       pup     ::Array{Int64,1},
                       prints  ::Int64)

  el  = lastindex(idf)
  L   = treelength(Ξ)          # tree length
  ns  = nnodesbifurcation(idf) # number of speciation events
  nin = lastindex(inodes)      # number of internal nodes
  ne  = 0.0                    # number of extinction events

  # likelihood
  llc = llik_cpe(Ξ, idf, λc, μc, σac, σkc, nnodesbifurcation(idf)) - rmλ * log(λc) + 
        log(mc) + prob_ρ(idf)

  # prior
  prc = logdgamma(λc, λ_prior[1], λ_prior[2])         +
        logdgamma(μc, μ_prior[1], μ_prior[2])         +
        logdinvgamma(σac^2, σa_prior[1], σa_prior[2]) +
        logdinvgamma(σkc^2, σk_prior[1], σk_prior[2])

  # tracked quantities
  sσa, sσk = gibbs_quanta(Ξ, idf)

  # n number to sum to ns for σa updates
  nσs = Float64(lastindex(idf)) - 2.0*nnodesbifurcation(idf) - rmλ

 # empty vector
  xis = Float64[]
  xfs = Float64[]
  es  = Float64[]
  pv  = Float64[]

  pbar = Progress(nburn, dt = prints, desc = "burn-in mcmc...", barlen = 20)

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

        llc, prc, σac = 
          update_σ!(σac, 0.5*sσa, 2.0*ns + nσs, llc, prc, σa_prior)

      # σk (cladogenetic) proposal
      elseif p === 4

        llc, prc, σkc = update_σ!(σkc, 0.5*sσk, ns, llc, prc, σk_prior)

      # update inner nodes traits
      elseif p === 5

        bix = inodes[fIrand(nin) + 1]
        llc, sσa, sσk = update_x!(bix, Ξ, idf, σac, σkc, llc, sσa, sσk)

      # forward simulation proposal proposal
      else

        bix = fIrand(el) + 1
        llc, ns, ne, L, sσa, sσk = 
          update_fs!(bix, Ξ, idf, llc, λc, μc, σac, σkc, ns, ne, L, sσa, sσk,
            xis, xfs, es, pv)

      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, σac, σkc, mc, ns, ne, L, sσa, sσk, nσs
end




"""
    mcmc_cpe(Ξ       ::Vector{sTpe},
             idf     ::Array{iBffs,1},
             llc     ::Float64,
             prc     ::Float64,
             λc      ::Float64,
             μc      ::Float64,
             σac     ::Float64,
             σkc     ::Float64,
             mc      ::Float64,
             ns      ::Float64,
             ne      ::Float64,
             L       ::Float64,
             sσa     ::Float64, 
             sσk     ::Float64,
             th      ::Float64,
             rmλ     ::Float64,
             inodes  ::Vector{Int64},
             surv    ::Int64,
             λ_prior ::NTuple{2,Float64},
             μ_prior ::NTuple{2,Float64},
             σa_prior::NTuple{2,Float64},
             σk_prior::NTuple{2,Float64},
             pup     ::Array{Int64,1},
             niter   ::Int64,
             nthin   ::Int64,
             nflush  ::Int64,
             ofile   ::String,
             prints  ::Int64)

Sampling for constant birth-death punctuated equilibrium.
"""
function mcmc_cpe(Ξ       ::Vector{sTpe},
                  idf     ::Array{iBffs,1},
                  llc     ::Float64,
                  prc     ::Float64,
                  λc      ::Float64,
                  μc      ::Float64,
                  σac     ::Float64,
                  σkc     ::Float64,
                  mc      ::Float64,
                  ns      ::Float64,
                  ne      ::Float64,
                  L       ::Float64,
                  sσa     ::Float64, 
                  sσk     ::Float64,
                  nσs     ::Float64,
                  th      ::Float64,
                  rmλ     ::Float64,
                  inodes  ::Vector{Int64},
                  surv    ::Int64,
                  λ_prior ::NTuple{2,Float64},
                  μ_prior ::NTuple{2,Float64},
                  σa_prior::NTuple{2,Float64},
                  σk_prior::NTuple{2,Float64},
                  pup     ::Array{Int64,1},
                  niter   ::Int64,
                  nthin   ::Int64,
                  nflush  ::Int64,
                  ofile   ::String,
                  prints  ::Int64)

  el  = lastindex(idf)
  nin = lastindex(inodes)

  # logging
  nlogs = fld(niter,nthin)
  lthin = lit = sthin = zero(Int64)

  # parameter results
  r = Array{Float64,2}(undef, nlogs, 8)

  # empty vectors
  xis = Float64[]
  xfs = Float64[]
  es  = Float64[]
  pv  = Float64[]

  treev = sTpe[]     # make tree vector
  io    = IOBuffer() # buffer 

  open(ofile*".log", "w") do of 
    write(of, "iteration\tlikelihood\tprior\tlambda\tmu\tx0\tsigma_a\tsigma_k\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, λc = λc, μc = μc, σac = σac, σkc = σkc, mc = mc, ns = ns, ne = ne, L = L, sσa = sσa, sσk = sσk, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for p in pup

             # λ proposal
            if p === 1

              llc, prc, λc, mc =
                update_λ!(llc, prc, λc, ns, L, μc, mc, th, rmλ, surv, λ_prior)

              # llci = llik_cpe(Ξ, idf, λc, μc, σac, σkc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # μ proposal
            elseif p === 2

              llc, prc, μc, mc =
                update_μ!(llc, prc, μc, ne, L, λc, mc, th, surv, μ_prior)

              # llci = llik_cpe(Ξ, idf, λc, μc, σac, σkc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # σa (anagenetic) proposal
            elseif p === 3

              llc, prc, σac = 
                update_σ!(σac, 0.5*sσa, 2.0*ns + nσs, llc, prc, σa_prior)

              # llci = llik_cpe(Ξ, idf, λc, μc, σac, σkc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # σk (cladogenetic) proposal
            elseif p === 4

              llc, prc, σkc = update_σ!(σkc, 0.5*sσk, ns, llc, prc, σk_prior)

              # llci = llik_cpe(Ξ, idf, λc, μc, σac, σkc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # update inner nodes traits
            elseif p === 5

              bix = inodes[fIrand(nin) + 1]

              llc, sσa, sσk = update_x!(bix, Ξ, idf, σac, σkc, llc, sσa, sσk)

              # llci = llik_cpe(Ξ, idf, λc, μc, σac, σkc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # forward simulation proposal proposal
            else

              bix = fIrand(el) + 1
              llc, ns, ne, L, sσa, sσk = 
                update_fs!(bix, Ξ, idf, llc, λc, μc, σac, σkc, ns, ne, L, 
                  sσa, sσk, xis, xfs, es, pv)

              # llci = llik_cpe(Ξ, idf, λc, μc, σac, σkc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            end
          end

          # log parameters
          lthin += one(Int64)
          if lthin === nthin

            lit += 1
            @inbounds begin
              r[lit,1] = Float64(it)
              r[lit,2] = llc
              r[lit,3] = prc
              r[lit,4] = λc
              r[lit,5] = μc
              r[lit,6] = xi(Ξ[1])
              r[lit,7] = σac
              r[lit,8] = σkc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', λc,'\t', μc, '\t', xi(Ξ[1]), '\t', σac, '\t', σkc, '\n')
            flush(of)
            ibuffer(io, couple(Ξ, idf, 1))
            write(io, '\n')
            write(tf, take!(io))
            flush(tf)
            sthin = zero(Int64)
          end

          next!(pbar)
        end
      end
    end
  end

  return r, treev
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

Perform a punkeek trait update for an internal branch and its descendants.
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

  ## update parent
  # if mrca
  if root && iszero(e(ξi))
    # if crown
    ξ2  = Ξ[i2]
    ll, sσa, sσk = _crown_update!(ξi, ξ1, ξ2, σa, σk, ll, sσa, sσk)
  else
  # if stem
    if root
      ll, sσa = _stem_update!(ξi, σa, ll, sσa)
    end

    # updates within the parent branch
    ll, sσa, sσk = 
      _update_node!(ξi, NaN, NaN, σa, σk, ll, sσa, sσk, false)

    # get fixed tip
    lξi = fixtip(ξi)

    isd = iszero(i2)
    # if duo
    if isd
      ll, sσa = _update_duo!(lξi, ξ1, σa, ll, sσa)
    # if triad
    else
      ll, sσa, sσk = _update_quartet!(lξi, ξ1, Ξ[i2], σa, σk, ll, sσa, sσk)
    end

    ## update daughters
    b1 = idf[i1]
    ll, sσa, sσk = 
      _update_node!(ξ1, xavg(b1), xstd(b1), σa, σk, 
        ll, sσa, sσk, iszero(d1(b1)))

    # if triad
    if !isd
      b2 = idf[i2]
      ll, sσa, sσk = 
        _update_node!(Ξ[i2], xavg(b2), xstd(b2), σa, σk, 
          ll, sσa, sσk, iszero(d1(b2)))
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
               sσk::Float64,
               xis::Vector{Float64},
               xfs::Vector{Float64},
               es ::Vector{Float64},
               pv ::Vector{Float64})

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
                    es ::Vector{Float64},
                    pv ::Vector{Float64})

  bi = idf[bix]
  ξc = Ξ[bix]

  llr = NaN
  sσar = sσkr = 0.0
   # if terminal branch
  if iszero(d1(bi))

    xav = xsd = NaN
    if !isnothing(xavg(bi))
      xav, xsd = xavg(bi), xstd(bi)
    end

    ξp, llr = fsbi_t(bi, xav, xsd, ξc, λ, μ, σa, σk, xis, xfs, es, pv)
  # if mid branch
  elseif iszero(d2(bi))

    ξp, llr, sσar = fsbi_m(bi, ξc, Ξ[d1(bi)], λ, μ, σa, σk, xfs, pv)
  # if trio branch
  elseif e(bi) > 0.0

    ξp, llr, sσar, sσkr = 
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], λ, μ, σa, σk, xfs, pv)
  end

  if isfinite(llr)

    σa2, σk2 = σa^2, σk^2

    ll1, ns1, ne1, L1, sσa1, sσk1 = 
      llik_cpe_track(ξp, λ, μ, σa2, σk2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ll0, ns0, ne0, L0, sσa0, sσk0 = 
      llik_cpe_track(ξc, λ, μ, σa2, σk2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    llc += ll1  - ll0 + llr
    ns  += ns1  - ns0
    ne  += ne1  - ne0
    L   += L1   - L0
    sσa += sσa1 - sσa0 + sσar
    sσk += sσk1 - sσk0 + sσkr

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, ns, ne, L, sσa, sσk
end




"""
    fsbi_t(bi ::iBffs,
           xav::Float64,
           xsd::Float64,
           ξi ::sTpe,
           λ  ::Float64,
           μ  ::Float64,
           σa ::Float64,
           σk ::Float64,
           xis::Vector{Float64},
           xfs::Vector{Float64},
           es ::Vector{Float64},
           pv ::Vector{Float64})

Forward simulation for terminal branch.
"""
function fsbi_t(bi ::iBffs,
                xav::Float64,
                xsd::Float64,
                ξi ::sTpe,
                λ  ::Float64,
                μ  ::Float64,
                σa ::Float64,
                σk ::Float64,
                xis::Vector{Float64},
                xfs::Vector{Float64},
                es ::Vector{Float64},
                pv ::Vector{Float64})

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
    _sim_cpe_t(e(bi), λ, μ, xi(ξi), σa, σk, lc, iρi, 0, 1, 500, 
               xis, xfs, es)

  if na < 1 || isnan(llr)
    return t0, NaN
  end

  # if fix node
  if ifx(bi)

    # propose trait value (if no uncertainty, then xp = xav)
    xp = xav
    if xsd > 0.0
      xp = rnorm(xav, xsd)
    end
    wt, acr, xp  = wfix_t(ξi, e(bi), xp, 0.0, xis, es, σa, na, nac, pv)

    if lU < acr + llr

      if wt <= div(na,2)
        fixtip1!(t0, wt, 0, xp)
      else
        fixtip2!(t0, na - wt + 1, 0, xp)
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
    wfix_t(ξi ::sTpe,
           ei ::Float64,
           xav::Float64,
           acr::Float64,
           xis::Vector{Float64},
           es ::Vector{Float64},
           σa ::Float64,
           na ::Int64,
           nac::Int64,
           pv ::Vector{Float64})

Choose most likely simulated lineage to fix with respect to the
trait value **without uncertainty** of terminal branches.
"""
function wfix_t(ξi ::sTpe,
                ei ::Float64,
                xav::Float64,
                acr::Float64,
                xis::Vector{Float64},
                es ::Vector{Float64},
                σa ::Float64,
                na ::Int64,
                nac::Int64,
                pv ::Vector{Float64})

  wt, sp, pp = 0, 0.0, NaN
  if isone(na)
    # sample from proposal
    empty!(pv)
    for i in Base.OneTo(na)
      p = dnorm(xav, xis[i], sqrt(es[i])*σa)
      push!(pv, p)
      sp += p
    end

    if iszero(sp)
      return 0, NaN, NaN
    end

    wt = _samplefast(pv, sp, na)
    pp = pv[wt]
  else
    pp = sp = dnorm(xav, xis[1], sqrt(es[1])*σa)
    wt = 1
  end

  # extract current `xis` and estimate ratio
  sc, pc = 0.0, NaN
  if isone(nac)
    empty!(xis)
    empty!(es)
    nac, xic = _xatt!(ξi, ei, xis, es, 0.0, 0, NaN)

    for i in Base.OneTo(nac)
      p   = dnorm(xav, xis[i], sqrt(es[i])*σa)
      sc += p
      if xic === xis[i]
        pc = p
      end
    end
  else
    pc = sc = 1.0
  end

  acr += log(pp) - log(pp/sp) + log(pc/sc)

  return wt, acr, xav
end




"""
    fsbi_m(bi::iBffs,
           ξi::sTpe,
           ξ1::sTpe,
           λ ::Float64,
           μ ::Float64,
           σa::Float64,
           σk::Float64,
           xfs::Vector{Float64},
           pv ::Vector{Float64})

Forward simulation for internal branch.
"""
function fsbi_m(bi::iBffs,
                ξi::sTpe,
                ξ1::sTpe,
                λ ::Float64,
                μ ::Float64,
                σa::Float64,
                σk::Float64,
                xfs::Vector{Float64},
                pv ::Vector{Float64})

  # forward simulation during branch length
  empty!(xfs)

  t0, na, nn = _sim_cpe_i(e(bi), λ, μ, xi(ξi), σa, σk, 0, 1, 500, xfs)

  if na < 1 || nn > 499
    return t0, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # add sampling fraction
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # inverse branch sampling fraction
  acr = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  ## choose most likely lineage to fix
  xp, wt, pp, pc, acr = wfix_m(ξi, ξ1, e(bi), acr, xfs, σa, na, pv)

  if lU < acr

    # fix the tip
    if wt <= div(na, 2)
      fixtip1!(t0, wt, 0, xp)
    else
      fixtip2!(t0, na - wt + 1, 0, xp)
    end

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr = 
        tip_sims!(t0, tf(bi), λ, μ, σa, σk, acr, lU, iρi, na, nn, 500)
    end

    if lU < acr
      na  -= 1
      llr  = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc)
      setnt!(bi, ntp)  # set new nt
      setni!(bi, na)   # set new ni

      sσar = ((xp - xf(ξ1))^2 - (xi(ξ1) - xf(ξ1))^2)/e(ξ1)
      setxi!(ξ1, xp)   # set new xp for initial x

      return t0, llr, sσar
    end
  end

  return t0, NaN, NaN
end




"""
    wfix_m(ξi ::sTpe,
           ξ1 ::sTpe,
           ei ::Float64,
           acr::Float64,
           xfs::Vector{Float64},
           σa ::Float64,
           na ::Int64,
           pv ::Vector{Float64})

Choose most likely simulated lineage to fix with respect to daughter
for `mid` branches.
"""
function wfix_m(ξi ::sTpe,
                ξ1 ::sTpe,
                ei ::Float64,
                acr::Float64,
                xfs::Vector{Float64},
                σa ::Float64,
                na ::Int64,
                pv ::Vector{Float64})

  # select best from proposal
  xf1, sre1σa = xf(ξ1), sqrt(e(ξ1))*σa
  
  empty!(pv)
  sp = 0.0
  for xfi in xfs
    p   = dnorm(xf1, xfi, sre1σa)
    push!(pv, p)
    sp += p
  end

  if iszero(sp)
    return NaN, 0, NaN, NaN, NaN
  end

  wt = _samplefast(pv, sp, na)
  pp = pv[wt]
  xp = xfs[wt]

  # extract current xis and estimate ratio
  empty!(xfs)
  xc, shc = _xatt!(ξi, σa^2, ei, xfs, 0.0, NaN, false)

  sc, pc = 0.0, NaN
  for xfi in xfs
    p   = dnorm(xf1, xfi, sre1σa)
    sc += p
    if xc === xfi
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp/sc)

  return xp, wt, pp, pc, acr
end




"""
    fsbi_i(bi ::iBffs,
           ξi ::sTpe,
           ξ1 ::sTpe,
           ξ2 ::sTpe,
           λ  ::Float64,
           μ  ::Float64,
           σa ::Float64,
           σk ::Float64,
           xfs::Vector{Float64},
           pv ::Vector{Float64})

Forward simulation for internal branch.
"""
function fsbi_i(bi ::iBffs,
                ξi ::sTpe,
                ξ1 ::sTpe,
                ξ2 ::sTpe,
                λ  ::Float64,
                μ  ::Float64,
                σa ::Float64,
                σk ::Float64,
                xfs::Vector{Float64},
                pv ::Vector{Float64})

  # forward simulation during branch length
  empty!(xfs)

  t0, na, nn = _sim_cpe_i(e(bi), λ, μ, xi(ξi), σa, σk, 0, 1, 500, xfs)

  if na < 1 || nn >= 500
    return t0, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # add sampling fraction
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  ## choose most likely lineage to fix
  wt, xp, xkp, shp, pp, xc, shc, pc, acr =
    wfix_i(ξi, ξ1, ξ2, e(bi), acr, xfs, σa^2, σk^2, na, pv) 

  if lU < acr

    # fix the tip
    if wt <= div(na,2)
      fixtip1!(t0, wt, 0, shp)
    else
      fixtip2!(t0, na - wt + 1, 0, shp)
    end

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr = 
        tip_sims!(t0, tf(bi), λ, μ, σa, σk, acr, lU, iρi, na, nn, 500)
    end

    if lU < acr
      na  -= 1
      llr  = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + pp - pc
      setnt!(bi,  ntp)  # set new nt
      setni!(bi,  na)   # set new ni

      ξac, ξkc = if shc ξ2, ξ1 else ξ1, ξ2 end
      ξap, ξkp = if shp ξ2, ξ1 else ξ1, ξ2 end

      sσar = (xp  - xf(ξap))^2/e(ξap) - (xc      - xf(ξac))^2/e(ξac) +
             (xkp - xf(ξkp))^2/e(ξkp) - (xi(ξkc) - xf(ξkc))^2/e(ξkc)

      sσkr = (xp - xkp)^2 - (xc - xi(ξkc))^2
      setxi!(ξap, xp)   # set new xp for initial anagenetic daughter
      setxi!(ξkp, xkp)  # set new xkp for initial cladogenetic daughter

      return t0, llr, sσar, sσkr
    end
  end

  return t0, NaN, NaN, NaN
end




# """
#     wfix_i(ξi ::sTpe,
#            ξ1 ::sTpe,
#            ξ2 ::sTpe,
#            ei ::Float64,
#            acr::Float64,
#            xfs::Vector{Float64},
#            σa2::Float64,
#            σk2::Float64)

# Choose most likely simulated lineage to fix with respect to daughter
# for bifurcating `i` branches.
# """
# function wfix_i(ξi ::sTpe,
#                 ξ1 ::sTpe,
#                 ξ2 ::sTpe,
#                 ei ::Float64,
#                 acr::Float64,
#                 xfs::Vector{Float64},
#                 σa2::Float64,
#                 σk2::Float64)

#   xi1, xi2, xf1, xf2, e1, e2 = xi(ξ1), xi(ξ2), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2)

#   ntp = lastindex(xfs)
#   # select one tip at random
#   wt  = fIrand(ntp) + 1
#   xfp = xfs[wt]
#   pk1 = llik_cpe_dyad(xfp, xf2, xf1, e2, e1, σa2, σk2)
#   pk2 = llik_cpe_dyad(xfp, xf1, xf2, e1, e2, σa2, σk2)
#   o12 = exp(pk1 - pk2)  # odds
#   p1  = o12/(1.0 + o12) # probability of d1 cladogenetic
#   shp = rand() < p1
#   lpp = shp ? pk1 - log(p1) : pk2 - log(1.0 - p1)

#   # proposal cladogenetic and likelihood
#   xkp, ll3p = NaN, NaN
#   if shp
#     xkp   = duoprop(xfp, xf1, σk2, e1*σa2)
#     ll3p  = llik_cpe_trio(xfp, xkp, xf2, xf1, e2, e1, σa2, σk2)
#   else
#     xkp  = duoprop(xfp, xf2, σk2, e2*σa2)
#     ll3p = llik_cpe_trio(xfp, xkp, xf1, xf2, e1, e2, σa2, σk2)
#   end

#   # extract current xis and estimate ratio
#   empty!(xfs)
#   xfc, shc = _xatt!(ξi, σa2, ei, xfs, 0.0, NaN, false)
#   ntc = Float64(lastindex(xfs))
#   pk1 = llik_cpe_dyad(xfc, xf2, xf1, e2, e1, σa2, σk2)
#   pk2 = llik_cpe_dyad(xfc, xf1, xf2, e1, e2, σa2, σk2)
#   o12 = exp(pk1 - pk2)  # odds
#   p1  = o12/(1.0 + o12) # probability of d1 cladogenetic

#   lpc, ll3c = NaN, NaN
#   if shc
#     lpc  = pk1 - log(p1)
#     ll3c = llik_cpe_trio(xfc, xi1, xf2, xf1, e2, e1, σa2, σk2)
#   else
#     lpc  = pk2 - log(1.0 - p1)
#     ll3c = llik_cpe_trio(xfc, xi2, xf1, xf2, e1, e2, σa2, σk2)
#   end

#   # add to acceptance ratio
#   acr += lpp - lpc + log(Float64(ntp)/ntc)

#   return wt, xfp, xkp, shp, ll3p, xfc, shc, ll3c, acr
# end





"""
    wfix_i(ξi ::sTpe,
           ξ1 ::sTpe,
           ξ2 ::sTpe,
           ei ::Float64,
           acr::Float64,
           xfs::Vector{Float64},
           σa2::Float64,
           σk2::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for bifurcating `i` branches.
"""
function wfix_i(ξi ::sTpe,
                ξ1 ::sTpe,
                ξ2 ::sTpe,
                ei ::Float64,
                acr::Float64,
                xfs::Vector{Float64},
                σa2::Float64,
                σk2::Float64,
                na ::Int64,
                pv ::Vector{Float64})

  xi1, xi2, xf1, xf2, e1, e2 = xi(ξ1), xi(ξ2), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2)

  # select best from proposal
  empty!(pv)
  sp = 0.0
  for xfi in xfs
    pk1 = exp(llik_cpe_dyad(xfi, xf2, xf1, e2, e1, σa2, σk2))
    pk2 = exp(llik_cpe_dyad(xfi, xf1, xf2, e1, e2, σa2, σk2))
    push!(pv, pk1, pk2)
    sp += pk1 + pk2
  end

  if iszero(sp)
    return 0, NaN, NaN, false, NaN, NaN, false, NaN, NaN
  end

  wi  = _samplefast(pv, sp, 2*na)
  pkp = pv[wi]
  pp  = pkp/sp
  shp = isodd(wi)
  wt  = div(wi, 2) + (shp ? 1 : 0)
  xp  = xfs[wt]

  # proposal cladogenetic and likelihood
  xkp, ll3p = NaN, NaN
  if shp
    xkp = duoprop(xp, xf1, σk2, e1*σa2)
    ll3p  = llik_cpe_trio(xp, xkp, xf2, xf1, e2, e1, σa2, σk2)
  else
    xkp = duoprop(xp, xf2, σk2, e2*σa2)
    ll3p  = llik_cpe_trio(xp, xkp, xf1, xf2, e1, e2, σa2, σk2)
  end

  # extract current xis and estimate ratio
  empty!(xfs)
  xc, shc = _xatt!(ξi, σa2, ei, xfs, 0.0, NaN, false)

  sc, ll3c = 0.0, NaN
  for xfi in xfs
    pk1 = exp(llik_cpe_dyad(xfi, xf2, xf1, e2, e1, σa2, σk2))
    pk2 = exp(llik_cpe_dyad(xfi, xf1, xf2, e1, e2, σa2, σk2))
    sc += pk1 + pk2
    if xc === xfi
      if shc
        pc  = pk1
        ll3c = llik_cpe_trio(xfi, xi1, xf2, xf1, e2, e1, σa2, σk2)
      else
        pc  = pk2
        ll3c = llik_cpe_trio(xfi, xi2, xf1, xf2, e1, e2, σa2, σk2)
      end
    end
  end
  pkc = pc
  pc /= sc

  # likelihood ratio and acceptance
  acr += log(pkp/pp) - log(pkc/pc)

  return wt, xp, xkp, shp, ll3p, xc, shc, ll3c, acr
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
                   σa  ::Float64,
                   σk  ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64)

  if lU < lr && nn < nlim

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nn, lr = 
          _sim_cpe_it(t, λ, μ, xf(tree), σa, σk, lr, lU, iρi, na-1, nn, nlim)

        if isnan(lr) || nn >= nlim
          return tree, na, nn, NaN
        end

        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        setxf!(tree, xf(stree))
        setsh!(tree, sh(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr = 
        tip_sims!(tree.d1, t, λ, μ, σa, σk, lr, lU, iρi, na, nn, nlim)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, λ, μ, σa, σk, lr, lU, iρi, na, nn, nlim)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end





