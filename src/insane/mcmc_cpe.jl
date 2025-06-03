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
               pupdp   ::NTuple{6,Float64}     = (0.2, 0.2, 0.2, 0.2, 0.2, 0.8),
               prints  ::Int64                 = 5,
               survival::Bool                  = true,
               mxthf   ::Float64               = Inf,
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
                    pupdp   ::NTuple{6,Float64}     = (0.2, 0.2, 0.2, 0.2, 0.2, 0.8),
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

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, xr, σac, σkc, sTpe)

  # get vector of internal edges
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(6)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running constant punctuated equilibrium"

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
                       inodes ::Vector{Int64},
                       surv   ::Int64,
                       pup    ::Array{Int64,1},
                       prints ::Int64)

  el  = lastindex(idf)
  L   = treelength(Ξ)          # tree length
  ns  = nnodesbifurcation(idf) # number of speciation events
  nin = lastindex(inodes)      # number of internal nodes
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
  sσa, sσk = ssσak(Ξ, idf)

  # n number to sum to ns for σa updates
  nσs = Float64(lastindex(idf)) - 2.0*ns - rmλ

 # empty vector
  xis = Float64[]
  xfs = Float64[]
  es  = Float64[]

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

        llc, prc, σac = 
          update_σ!(σac, 0.5*sσa, 2.0*ns + nσs, llc, prc, σa_prior)

      # σk (cladogenetic) proposal
      elseif p === 4

        llc, prc, σkc = update_σ!(σkc, 0.5*sσk, ns, llc, prc, σk_prior)

      # update inner nodes traits
      elseif p === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]
        llc, sσa, sσk = update_x!(bix, Ξ, idf, σac, σkc, llc, sσa, sσk)

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)
        llc, ns, ne, L, sσa, sσk = 
          update_fs!(bix, Ξ, idf, llc, λc, μc, σac, σkc, ns, ne, L, sσa, sσk,
            xis, xfs, es)

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

  # empty vector
  xis = Float64[]
  xfs = Float64[]
  es  = Float64[]

  treev = sTpe[]     # make tree vector
  io    = IOBuffer() # buffer 

  open(ofile*".log", "w") do of 
    write(of, "iteration\tlikelihood\tprior\tlambda\tmu\tx0\tsigma_a\tsigma_k\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, λc = λc, μc = μc, σac = σac, σkc = σkc, mc = mc, ns = ns, ne = ne, L = L, sσa = sσa, sσk = sσk, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, prints, "running mcmc...", 20)

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

              nix = ceil(Int64,rand()*nin)
              bix = inodes[nix]
              llc, sσa, sσk = update_x!(bix, Ξ, idf, σac, σkc, llc, sσa, sσk)

              # llci = llik_cpe(Ξ, idf, λc, μc, σac, σkc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # forward simulation proposal proposal
            else

              bix = ceil(Int64,rand()*el)
              llc, ns, ne, L, sσa, sσk = 
                update_fs!(bix, Ξ, idf, llc, λc, μc, σac, σkc, ns, ne, L, 
                  sσa, sσk, xis, xfs, es)

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
      ll, sσa, sσk = _update_node_x!(lξi, ξ1, ξ2, σa, σk, ll, sσa, sσk)
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

  llr = NaN
  sσar = sσkr = 0.0
   # if terminal branch
  if iszero(d1(bi))

    xav = xsd = NaN
    if !isnothing(xavg(bi))
      xav, xsd = xavg(bi), xstd(bi)
    end

    ξp, llr = fsbi_t(bi, xav, xsd, ξc, λ, μ, σa, σk, xis, xfs, es)
  # if mid branch
  elseif iszero(d2(bi))

    ξp, llr, sσar = fsbi_m(bi, ξc, Ξ[d1(bi)], λ, μ, σa, σk, xis, xfs)
  # if trio branch
  elseif e(bi) > 0.0

    ξp, llr, sσar, sσkr = 
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], λ, μ, σa, σk, xis, xfs)
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
function fsbi_t(bi ::iBffs,
                xav::Float64,
                xst::Float64,
                ξc ::sTpe,
                λ  ::Float64,
                μ  ::Float64,
                σa ::Float64,
                σk ::Float64,
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
    _sim_cpe_t(e(bi), λ, μ, xi(ξc), σa, σk, lc, iρi, 0, 1, 500, 
               xis, xfs, es)

  if na < 1 || isnan(llr)
    return t0, NaN
  end

  # if fix node
  if ifx(bi)

    # if no uncertainty around trait value
    if iszero(xst)
       wti, acr, xp  = wfix_t(ξc, e(bi), xav, 0.0, xis, es, σa, na)

    # if uncertainty around trait value
    else
       wti, acr, xp  = wfix_t(ξc, e(bi), xav, xst, 0.0, xis, xfs, es, σa, na)
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
    wfix_t(ξi ::T, 
           ξ1 ::T,
           ei ::Float64,
           acr::Float64,
           xfs::Vector{Float64}, 
           xcs::Vector{Float64}, 
           σa2::Float64) where {T <: Tpe}

Choose most likely simulated lineage to fix with respect to the
trait value **without uncertainty** of terminal branches.
"""
function wfix_t(ξi ::T,
                ei ::Float64,
                xav::Float64,
                acr::Float64,
                xis::Vector{Float64},
                es ::Vector{Float64},
                σa ::Float64,
                na ::Int64) where {T <: Tpe}

  # select best from proposal
  sp, wt, pp = 0.0, 0, -Inf
  for i in Base.OneTo(na)
    p   = dnorm(xav, xis[i], sqrt(es[i])*σa)
    sp += p
    if p > pp
      pp = p
      wt = i
    end
  end

  # extract current xcs and estimate ratio
  empty!(xis)
  empty!(es)
  nac, xic = _xisatt!(ξi, ei, xis, es, 0.0, 0, NaN)

  sc = 0.0
  for i in Base.OneTo(nac)
    sc += dnorm(xav, xis[i], sqrt(es[i])*σa)
  end

  # likelihood ratio and acceptance
  acr += log(sp) - log(sc)

  return wt, acr, xav
end






"""
    wfix_t(ξi ::T,
           ei ::Float64,
           xav::Float64,
           xst::Float64,
           acr::Float64,
           xis::Vector{Float64},
           xfs::Vector{Float64},
           es ::Vector{Float64},
           σa ::Float64,
           na ::Int64) where {T <: Tpe}

Choose most likely simulated lineage to fix with respect to the
trait value **with uncertainty** of terminal branches.
"""
function wfix_t(ξi ::T,
                ei ::Float64,
                xav::Float64,
                xst::Float64,
                acr::Float64,
                xis::Vector{Float64},
                xfs::Vector{Float64},
                es ::Vector{Float64},
                σa ::Float64,
                na ::Int64) where {T <: Tpe}

  # select best from proposal
  sp, wt, xp, pp = 0.0, 0, NaN, -Inf
  for i in Base.OneTo(na)
    p   = duodnorm(xfs[i], xis[i], xav, sqrt(es[i])*σa, xst)
    sp += p
    if p > pp
      pp = p
      xp = xfs[i]
      wt = i
    end
  end

  # extract current xcs and estimate ratio
  empty!(xix)
  empty!(xfx)
  empty!(es)
  nac, xc, xic = _xisatt!(ξi, ei, xcs, es, 0.0, 0, NaN, NaN)

  sc, pc = 0.0, NaN
  for i in Base.OneTo(nac)
    p   = duodnorm(xc, xis[i], xav, sqrt(es[i])*σa, xst)
    sc += p
    if xis[i] === xic
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp) - log(sc)

  return wt, acr, xp
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
           xcs::Vector{Float64})

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
                xcs::Vector{Float64})

  # forward simulation during branch length
  empty!(xfs)

  t0, na, nn = _sim_cpe_i(e(bi), λ, μ, xi(ξi), σa, σk, 0, 1, 500, xfs)

  if na < 1 || nn >= 500
    return t0, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # add sampling fraction
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # inverse branch sampling fraction
  acr = Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  ## choose most likely lineage to fix
  xp, wt, pp, pc, acr = wfix_m(ξi, ξ1, e(bi), acr, xfs, xcs, σa)

  if lU < acr

    # fix the tip
    if wt <= div(na,2)
      fixtip1!(t0, wt, 0)
    else
      fixtip2!(t0, na - wt + 1, 0)
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
    wfix_m(ξi ::T, 
           ξ1 ::T,
           ei ::Float64,
           acr::Float64,
           xfs::Vector{Float64}, 
           xcs::Vector{Float64}, 
           σa2::Float64) where {T <: Tpe}

Choose most likely simulated lineage to fix with respect to daughter
for `mid` branches.
"""
function wfix_m(ξi ::T,
                ξ1 ::T,
                ei ::Float64,
                acr::Float64,
                xfs::Vector{Float64},
                xcs::Vector{Float64},
                σa ::Float64) where {T <: Tpe}

  # select best from proposal
  xf1, sre1 = xf(ξ1), sqrt(e(ξ1))
  sp, i, wt, xp, pp = 0.0, 0, 0, NaN, -Inf
  for xfi in xfs
    p   = dnorm(xfi, xf1, sre1*σa)
    sp += p
    i  += 1
    if p > pp
      pp  = p
      xp  = xfi
      wt  = i
    end
  end

  # extract current xcs and estimate ratio
  empty!(xcs)
  xc, shc = _xatt!(ξi, ei, xcs, 0.0, NaN, false)

  sc, pc = 0.0, NaN
  for xci in xcs
    p   = dnorm(xci, xf1, sre1*σa)
    sc += p
    if xc === xci
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp) - log(sc)

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
           xfs::Vector{Float64})

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
                xcs::Vector{Float64})

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
  acr = Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  ## choose most likely lineage to fix
  wt, xp, shp, pp, xc, shc, pc, acr = 
    wfix_i(ξi, ξ1, ξ2, e(bi), acr, xfs, xcs, σa^2, σk^2) 

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

      sσar = (xp - xf(ξap))^2/e(ξap)      - (xc - xf(ξac))^2/e(ξac)      +
             (xi(ξkp) - xf(ξkp))^2/e(ξkp) - (xi(ξkc) - xf(ξkc))^2/e(ξkc)
      sσkr = (xp - xi(ξkp))^2 - (xc - xi(ξkc))^2
      setxi!(ξap, xp)   # set new xp for initial anagenetic daughter

      return t0, llr, sσar, sσkr
    end
  end

  return t0, NaN, NaN, NaN
end




"""
    wfix_i(ξi ::T,
           ξ1 ::T,
           ξ2 ::T,
           ei ::Float64,
           acr::Float64,
           xfs::Vector{Float64},
           xcs::Vector{Float64},
           σa2::Float64,
           σk2::Float64) where {T <: Tpe}

Choose most likely simulated lineage to fix with respect to daughter
for bifurcating `i` branches.
"""
function wfix_i(ξi ::T,
                ξ1 ::T,
                ξ2 ::T,
                ei ::Float64,
                acr::Float64,
                xfs::Vector{Float64},
                xcs::Vector{Float64},
                σa2::Float64,
                σk2::Float64) where {T <: Tpe}

  # select best from proposal
  sp, i, wt, xp, pp, shp = 0.0, 0, 0, NaN, -Inf, false
  for xfi in xfs
    i  += 1
    pk1 = llik_trio(xfi, xi(ξ1), xf(ξ2), xf(ξ1), e(ξ2), e(ξ1), σa2, σk2)
    pk2 = llik_trio(xfi, xi(ξ2), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2), σa2, σk2)
    sp += exp(pk1) + exp(pk2)
    pfi = max(pk1, pk2)

    if pfi > pp
      pp  = pfi
      xp  = xfi
      shp = pk1 > pk2
      wt  = i
    end
  end

  # extract current xcs and estimate ratio
  empty!(xcs)
  xc, shc = _xatt!(ξi, ei, xcs, 0.0, NaN, false)

  sc, pc = 0.0, NaN
  for xci in xcs
    pk1 = llik_trio(xci, xi(ξ1), xf(ξ2), xf(ξ1), e(ξ2), e(ξ1), σa2, σk2)
    pk2 = llik_trio(xci, xi(ξ2), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2), σa2, σk2)
    sc += exp(pk1) + exp(pk2)

    if xc === xci
      pc = shc ? pk1 : pk2
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp) - log(sc)

  return wt, xp, shp, pp, xc, shc, pc, acr
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





