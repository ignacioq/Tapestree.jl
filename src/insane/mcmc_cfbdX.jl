#=

constant birth-death MCMC using forward simulation with traits

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cfbd(tree    ::sT_label,
               X       ::Dict{String, Float64},
               out_file::String;
               λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               σx_prior::NTuple{2,Float64}     = (0.05, 0.05),
               x0_prior::NTuple{2,Float64}     = (0.0, 10.0),
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
function insane_cfbd(tree    ::sTf_label,
                     X       ::Dict{String, Float64},
                     out_file::String;
                     λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     σx_prior::NTuple{2,Float64}     = (0.05, 0.05),
                     x0_prior::NTuple{2,Float64}     = (0.0, 10.0),
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
                     pupdp   ::NTuple{6,Float64}     = (0.2,0.2,0.2,0.2,0.2,0.2),
                     prints  ::Int64                 = 5,
                     tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxc = make_idf(tree, tρ, X)

  # other starting parameters
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
      crown = 1
    # if stem conditioning
    else
      crown = 0
    end
  # no survival
  else
    crown = 2
  end
  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 500, crown)

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, xr, sTfbdX)

  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(Ξ))
    if !it(idf[i]) || isfossil(idf[i])
      push!(inodes, i)
    end
  end

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(6)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant birth-death and trait evolution"

  # adaptive phase
  llc, prc, λc, μc, ψc, σxc, mc =
      mcmc_burn_cfbd(Ξ, idf, λ_prior, μ_prior, ψ_prior, σx_prior, x0_prior, nburn,
        λc, μc, ψc, σxc, mc, th, crown, inodes, pup, prints)

  # mcmc
  r, treev =
    mcmc_cfbd(Ξ, idf, llc, prc, λc, μc, ψc, σxc, mc, th, crown,
      λ_prior, μ_prior, ψ_prior, σx_prior, x0_prior, niter, nthin, inodes, pup, prints)

  pardic = Dict("lambda"  => 1,
                "mu"      => 2,
                "x0"      => 3,
                "sigma_x" => 4)

  write_ssr(r, pardic, out_file)

  return r, treev, NaN
end




"""
    mcmc_burn_cfbd(Ξ      ::Vector{sTbd},
                  idf    ::Array{iBffs,1},
                  λ_prior::NTuple{2,Float64},
                  μ_prior::NTuple{2,Float64},
                  nburn  ::Int64,
                  λc     ::Float64,
                  μc     ::Float64,
                  mc     ::Float64,
                  th     ::Float64,
                  crown   ::Bool,
                  pup    ::Array{Int64,1},
                  prints ::Int64)

Adaptive MCMC phase for da chain for constant birth-death using forward
simulation.
"""
function mcmc_burn_cfbd(Ξ       ::Vector{sTfbdX},
                        idf     ::Array{iBffs,1},
                        λ_prior ::NTuple{2,Float64},
                        μ_prior ::NTuple{2,Float64},
                        ψ_prior ::NTuple{2,Float64},
                        σx_prior::NTuple{2,Float64},
                        x0_prior::NTuple{2,Float64},
                        nburn   ::Int64,
                        λc      ::Float64,
                        μc      ::Float64,
                        ψc      ::Float64,
                        σxc     ::Float64,
                        mc      ::Float64,
                        th      ::Float64,
                        crown    ::Int64,
                        inodes  ::Array{Int64,1},
                        pup     ::Array{Int64,1},
                        prints  ::Int64)

  el      = lastindex(idf)
  nin     = lastindex(inodes)
  L       = treelength(Ξ)        # tree length
  nf      = Float64(nfossils(Ξ)) # number of fossilization events
  ns      = Float64(nnodesbifurcation(Ξ)) # number of speciation events
  ne      = Float64(ntipsextinct(Ξ))      # number of extinction events
  sdX, nX = sdeltaX(Ξ)           # standardized trait differences

  # likelihood
  llc = llik_cfbd(Ξ, λc, μc, ψc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
  prc = logdgamma(λc,       λ_prior[1],  λ_prior[2])    +
        logdgamma(μc,       μ_prior[1],  μ_prior[2])    +
        logdgamma(ψc,        ψ_prior[1], ψ_prior[2])    +
        logdinvgamma(σxc^2, σx_prior[1], σx_prior[2])   +
        logdnorm(xi(Ξ[1]),  x0_prior[1], x0_prior[2]^2)

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, crown, λ_prior)

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, crown, μ_prior)

       # psi update
      elseif p === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

       # sigma_x update
      elseif p === 4

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

      # X ancestors update
      elseif p === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX =
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, x0_prior)

      # forward simulation proposal update
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L, sdX, nX =
          update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, σxc, ns, ne, L, sdX, nX)

      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, ψc, σxc, mc
end




"""
    mcmc_cfbd(Ξ      ::Vector{sTbd},
             idf    ::Array{iBffs,1},
             llc    ::Float64,
             prc    ::Float64,
             λc     ::Float64,
             μc     ::Float64,
             mc     ::Float64,
             th     ::Float64,
             crown   ::Bool,
             λ_prior::NTuple{2,Float64},
             μ_prior::NTuple{2,Float64},
             niter  ::Int64,
             nthin  ::Int64,
             pup    ::Array{Int64,1},
             prints ::Int64)

MCMC da chain for constant birth-death using forward simulation.
"""
function mcmc_cfbd(Ξ      ::Vector{sTfbdX},
                   idf     ::Array{iBffs,1},
                   llc     ::Float64,
                   prc     ::Float64,
                   λc      ::Float64,
                   μc      ::Float64,
                   ψc      ::Float64,
                   σxc     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   crown    ::Int64,
                   λ_prior ::NTuple{2,Float64},
                   μ_prior ::NTuple{2,Float64},
                   ψ_prior ::NTuple{2,Float64},
                   σx_prior::NTuple{2,Float64},
                   x0_prior::NTuple{2,Float64},
                   niter   ::Int64,
                   nthin   ::Int64,
                   inodes  ::Array{Int64,1},
                   pup     ::Array{Int64,1},
                   prints  ::Int64)

  el      = lastindex(idf)
  nin     = lastindex(inodes)
  ns      = Float64(nnodesbifurcation(Ξ))
  ne      = Float64(ntipsextinct(Ξ))
  nf      = Float64(nfossils(Ξ))
  L       = treelength(Ξ)
  sdX, nX = sdeltaX(Ξ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 7)

  # make tree vector
  treev  = sTfbdX[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, crown, λ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, crown, μ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

       # psi update
      elseif p === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

       # sigma_x update
      elseif p === 4

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # X ancestors update
      elseif p === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX =
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, x0_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L, sdX, nX =
          update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, σxc, ns, ne, L, sdX, nX)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
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
        R[lit,6] = xi(Ξ[1])
        R[lit,7] = σxc
        push!(treev, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev
end




"""
    update_fs!(bix::Int64,
               Ξ  ::Vector{sTfbdX},
               idf::Vector{iBffs},
               llc::Float64,
               λ  ::Float64,
               μ  ::Float64,
               ψ  ::Float64,
               σx ::Float64,
               ns ::Float64,
               ne ::Float64,
               L  ::Float64,
               sdX::Float64,
               nX ::Float64)

Forward simulation proposal function for constant fossilized birth-death.
"""
function update_fs!(bix::Int64,
                    Ξ  ::Vector{sTfbdX},
                    idf::Vector{iBffs},
                    llc::Float64,
                    λ  ::Float64,
                    μ  ::Float64,
                    ψ  ::Float64,
                    σx ::Float64,
                    ns ::Float64,
                    ne ::Float64,
                    L  ::Float64,
                    sdX::Float64,
                    nX ::Float64)

  bi = idf[bix]
  ξc = Ξ[bix]

  if it(bi)
    if isfossil(bi)
      ξp, llr = fsbi_ft(bi, ξc, λ, μ, ψ, σx)
    else
      ξp, llr = fsbi_t( bi, ξc, λ, μ, ψ, σx)
    end
    sdXr = 0.0
  else
    if isfossil(bi)
      ξp, llr, sdXr = fsbi_fi(bi, ξc, Ξ[d1(bi)], λ, μ, ψ, σx)
    else
      ξp, llr, sdXr = fsbi_i( bi, Ξ[d1(bi)], Ξ[d2(bi)], xi(ξc), λ, μ, ψ, σx)
    end
  end

  if isfinite(llr)

    # update llc, ns, ne & L
    llc += llik_cfbd(ξp, λ, μ, ψ, σx)    - llik_cfbd(ξc, λ, μ, ψ, σx) + llr
    ns  += Float64(nnodesbifurcation(ξp) - nnodesbifurcation(ξc))
    ne  += Float64(ntipsextinct(ξp)      - ntipsextinct(ξc))
    L   += treelength(ξp)                - treelength(ξc)
    sdXp, nXp = _sdeltaX(ξp, 0.0, 0.0)
    sdXc, nXc = _sdeltaX(ξc, 0.0, 0.0)
    sdX += sdXp - sdXc + sdXr
    nX  += nXp  - nXc

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, ns, ne, L, sdX, nX
end




"""
    fsbi_t(bi::iBffs,
           ξc::sTfbdX,
           λ ::Float64,
           μ ::Float64,
           ψ ::Float64,
           σx::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs,
                ξc::sTfbdX,
                λ ::Float64,
                μ ::Float64,
                ψ ::Float64,
                σx::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  xist = Float64[]
  xfst = Float64[]
  est  = Float64[]

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_cfbd_t(e(bi), λ, μ, ψ, xi(ξc), σx, lc, lU, Iρi, 0, 1, 500,
      xist, xfst, est)

  if na < 1 || isnan(llr) || nn >= 500
    return t0, NaN
  end

  # if fix node
  if fx(bi)
    #get fix `x` and edge
    xc  = fixed_xt(ξc)

    # sample tip according to fix `x` value
    acr = 0.0
    wp = Float64[]
    @simd for i in Base.OneTo(na)
      srt = sqrt(est[i])
      wi  = dnorm_bm(xfst[i], xist[i], srt*σx)/dnorm_bm(xc, xist[i], srt*σx)
      push!(wp, wi)
      acr += wi
    end
    acr = log(acr)

    if isfinite(acr) && lU <  acr + llr
      # sample tip
      wti = sample(wp)

      if wti <= div(na,2)
        fixtip1!(t0, wti, 0, xc)
      else
        fixtip2!(t0, na - wti + 1, 0, xc)
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
    fsbi_ft(bi::iBffs,
            ξc::sTfbdX,
            λ ::Float64,
            μ ::Float64,
            ψ ::Float64,
            σx::Float64)

Forward simulation for fossil terminal branch.
"""
function fsbi_ft(bi::iBffs,
                 ξc::sTfbdX,
                 λ ::Float64,
                 μ ::Float64,
                 ψ ::Float64,
                 σx::Float64)

  lU = -randexp() # log-probability

  ntp = na

  # add sampling fraction
  nac = ni(bi)                # current ni
  Iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = - Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if fx(bi)

    xist = Float64[]
    xfst = Float64[]
    est  = Float64[]

    # forward simulation during branch length
    t0, na, nf, nn  =
      _sim_cfbd_i(e(bi), λ, μ, ψ, xi(ξc), σx, 0, 0, 1, 500, xist, xfst, est)

    if na < 1 || nf > 0 || nn >= 500
      return t0, NaN
    end

    xc = fossil_xt(ξc)

    # sample tip according to fix `x` value
    wp = Float64[]
    ap = 0.0
    @simd for i in Base.OneTo(na)
      srt = sqrt(est[i])
      wi  = dnorm_bm(xfst[i], xist[i], srt*σx)/dnorm_bm(xc, xist[i], srt*σx)
      push!(wp, wi)
      ap += wi
    end
    acr = log(ap)

    if isfinite(acr) && lU < acr
      wti = sample(wp)

      if wti <= div(na,2)
        fixtip1!(t0, wti, 0, xc)
      else
        fixtip2!(t0, na - wti + 1, 0, xc)
      end
    else
      return t0, NaN
    end
  # if not a fix node
  else
    t0, na, nf, nn  = _sim_cfbd_i(e(bi), λ, μ, ψ, xi(ξc), σx, 0, 0, 1, 500)

    if na < 1 || nf > 0 || nn >= 500
      return t0, NaN
    end

    acr = 0.0
    _fixrtip!(t0, na)
  end


  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, σx, acr, lU, Iρi, na, nn)
    end

    if lU < acr

      # fossilize extant tip
      fossilizefixedtip!(t0)

      tx, na, nn, acr =
        fossiltip_sim!(t0, tf(bi), λ, μ, ψ, σx, acr, lU, Iρi, na, nn)

      if lU < acr
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
    fsbi_fi(bi::iBffs,
            ξ1::sTfbdX,
            x0::Float64,
            λ ::Float64,
            μ ::Float64,
            ψ ::Float64,
            σx::Float64)

Forward simulation for fossil internal branch.
"""
function fsbi_fi(bi::iBffs,
                 ξc::sTfbdX,
                 ξ1::sTfbdX,
                 λ ::Float64,
                 μ ::Float64,
                 ψ ::Float64,
                 σx::Float64)

  if fx(bi)

    xist = Float64[]
    xfst = Float64[]
    est  = Float64[]

    # forward simulation during branch length
    t0, na, nf, nn  =
      _sim_cfbd_i(e(bi), λ, μ, ψ, xi(ξc), σx, 0, 0, 1, 500, xist, xfst, est)

    if na < 1 || nf > 0 || nn >= 500
      return t0, NaN, NaN
    end

    xc = fossil_xt(ξc)

    # sample tip according to fix `x` value
    wp = Float64[]
    ap = 0.0
    @simd for i in Base.OneTo(na)
      srt = sqrt(est[i])
      wi  = dnorm_bm(xfst[i], xist[i], srt*σx)/dnorm_bm(xc, xist[i], srt*σx)
      push!(wp, wi)
      ap += wi
    end
    acr = log(ap)

    if isfinite(acr)
      wti = sample(wp)

      if wti <= div(na,2)
        fixtip1!(t0, wti, 0, xc)
      else
        fixtip2!(t0, na - wti + 1, 0, xc)
      end
    else
      return t0, NaN, NaN
    end
    xp  = xc
    llr = 0.0
  # if not a fix node
  else
    t0, na, nf, nn = _sim_cfbd_i(e(bi), λ, μ, ψ, xi(ξc), σx, 0, 0, 1, 500)

    if na < 1 || nf > 0 || nn >= 500
      return t0, NaN, NaN
    end
    xp = fixrtip!(t0, na, NaN)
    # acceptance ration with respect to daughter
    llr = lrdnorm_bm_x(xp, xi(ξ1), xf(ξ1), sqrt(e(ξ1))*σx)

    acr = 0.0
  end

  lU = -randexp() # log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if lU < acr + llr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, σx, acr, lU, Iρi, na, nn)
    end

    if lU < acr + llr
      # fossilize extant tip
      fossilizefixedtip!(t0)
      na -= 1

      llr += (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      sdX = ((xp - xf(ξ1))^2 - (xi(ξ1) - xf(ξ1))^2)/(2.0*e(ξ1))
      setnt!(bi, ntp)                # set new nt
      setni!(bi, na)                 # set new ni
      setxi!(ξ1, xp)                 # set new xp for initial x

      return t0, llr, sdX
    end
  end

  return t0, NaN, NaN
end




"""
    fsbi_i(bi::iBffs,
           ξ1::sTfbdX,
           ξ2::sTfbdX,
           x0::Float64,
           λ ::Float64,
           μ ::Float64,
           ψ ::Float64,
           σx::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi::iBffs,
                ξ1::sTfbdX,
                ξ2::sTfbdX,
                x0::Float64,
                λ ::Float64,
                μ ::Float64,
                ψ ::Float64,
                σx::Float64)

  # forward simulation during branch length
  t0, na, nf, nn = _sim_cfbd_i(e(bi), λ, μ, ψ, x0, σx, 0, 0, 1, 500)

  if na < 1 || nf > 0 || nn >= 500
    return t0, NaN, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # get fix `x`
  xp = fixrtip!(t0, na, NaN)

  # acceptance ration with respect to daughters
  llr = duoldnorm(xp,     xf(ξ1), xf(ξ2), e(ξ1), e(ξ2), σx) -
        duoldnorm(xi(ξ1), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2), σx)

  if lU < acr + llr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, σx, acr, lU, Iρi, na, nn)
    end

    if lU < acr + llr
      na  -= 1

      llr += (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      sdX = ((xp - xf(ξ1))^2 - (xi(ξ1) - xf(ξ1))^2)/(2.0*e(ξ1)) +
            ((xp - xf(ξ2))^2 - (xi(ξ2) - xf(ξ2))^2)/(2.0*e(ξ2))
      setnt!(bi, ntp)                # set new nt
      setni!(bi, na)                 # set new ni
      setxi!(ξ1, xp)                 # set new xp for initial x
      setxi!(ξ2, xp)                 # set new xp for initial x

      return t0, llr, sdX
    end
  end

  return t0, NaN, NaN
end




"""
    tip_sims!(tree::sTfbdX,
              t   ::Float64,
              λ   ::Float64,
              μ   ::Float64,
              ψ   ::Float64,
              σx  ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTfbdX,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   σx  ::Float64,
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
          _sim_cfbd_it(t, λ, μ, ψ, xf(tree), σx, lr, lU, Iρi, na-1, nn, 500)

        if isnan(lr) || nn >= 500
          return tree, na, nn, NaN
        end

        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        setxf!(tree, xf(stree))
        if def1(stree)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr =
        tip_sims!(tree.d1, t, λ, μ, ψ, σx, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr =
        tip_sims!(tree.d2, t, λ, μ, ψ, σx, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::sTfbdX,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   σx  ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

Continue simulation until time `t` for the fixed tip in `tree`.
"""
function fossiltip_sim!(tree::sTfbdX,
                        t   ::Float64,
                        λ   ::Float64,
                        μ   ::Float64,
                        ψ   ::Float64,
                        σx  ::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        Iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64)

  if lU < lr && nn < 500
    if istip(tree)

      stree, na, nn, lr =
        _sim_cfbd_it(t, λ, μ, ψ, xf(tree), σx, lr, lU, Iρi, na-1, nn, 500)

      if isnan(lr) || nn >= 500
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr =
        fossiltip_sim!(tree.d1, t, λ, μ, ψ, σx, lr, lU, Iρi, na, nn)
    else
      tree.d2, na, nn, lr =
        fossiltip_sim!(tree.d2, t, λ, μ, ψ, σx, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_x!(bix     ::Int64,
              Ξ       ::Vector{sTfbdX},
              idf     ::Vector{iBffs},
              σx      ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              sdX     ::Float64,
              x0_prior::NTuple{2, Float64})

Make a `gbm` update for an internal branch and its descendants.
"""
function update_x!(bix     ::Int64,
                   Ξ       ::Vector{sTfbdX},
                   idf     ::Vector{iBffs},
                   σx      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   sdX     ::Float64,
                   x0_prior::NTuple{2, Float64})

  ξi  = Ξ[bix]
  bi  = idf[bix]

  if !it(bi)
    id1 = d1(bi)
    id2 = d2(bi)
    ξ1  = Ξ[id1]
    if !isfossil(bi)
      ξ2  = Ξ[id2]
    end
  end

  root = iszero(pa(bi))
  # if crown root
  if root && iszero(e(bi))
    if !fx(bi)
      llc, prc, sdX =
         _crown_update_x!(ξi, ξ1, ξ2, σx, llc, prc, sdX, x0_prior)
    end
  else
    # if crown
    if root
      if !fx(bi)
        llc, prc, sdX = _stem_update_x!(ξi, σx, llc, prc, sdX, x0_prior)
      end
    end

    # updates within the parent branch
    llc, sdX = _update_x!(ξi, σx, llc, sdX, !fx(bi) && it(bi))

    if !it(bi)
    # get fixed tip
      lξi = fixtip(ξi)

      # make between decoupled trees node update
      if !fx(bi)
        if isfossil(bi)
          llc, sdX = _update_duo_x!(  lξi, ξ1, σx, llc, sdX)
        else
          llc, sdX = _update_triad_x!(lξi, ξ1, ξ2, σx, llc, sdX)
        end
      end
    end
  end

  if !it(bi)
    # carry on updates in the daughters
    llc, sdX = _update_x!(ξ1, σx, llc, sdX, !fx(idf[id1]) && it(idf[id1]))
    if !isfossil(bi)
      llc, sdX = _update_x!(ξ2, σx, llc, sdX, !fx(idf[id2]) && it(idf[id2]))
    end
  end

  return llc, prc, sdX
end




"""
    _update_x!(tree::T,
               σx  ::Float64,
               llc ::Float64,
               sdX ::Float64,
               ufx ::Bool) where {T <: sTX}

Do gbm updates on a decoupled tree recursively.
"""
function _update_x!(tree::sTfbdX,
                    σx  ::Float64,
                    llc ::Float64,
                    sdX ::Float64,
                    ufx ::Bool)

  if def1(tree)
    if def2(tree)
      llc, sdX = _update_triad_x!(tree, tree.d1, tree.d2, σx, llc, sdX)
      llc, sdX = _update_x!(tree.d1, σx, llc, sdX, ufx)
      llc, sdX = _update_x!(tree.d2, σx, llc, sdX, ufx)
    else
      if ufx
        llc, sdX = _update_duo_x!(tree, tree.d1, σx, llc, sdX)
      end
      llc, sdX = _update_x!(tree.d1, σx, llc, sdX, ufx)
    end
  elseif isfix(tree)
    if ufx
      llc, sdX = _update_tip_x!(tree, σx, llc, sdX)
    end
  else
    llc, sdX = _update_tip_x!(tree, σx, llc, sdX)
  end

  return llc, sdX
end




"""
    _update_duo_x!(tree::T,
                   tre1::T,
                   σx  ::Float64,
                   llc ::Float64,
                   sdX ::Float64) where {T <: sTX}

Make gibbs node update for trait.
"""
function _update_duo_x!(tree::T,
                        tre1::T,
                        σx  ::Float64,
                        llc ::Float64,
                        sdX ::Float64) where {T <: sTX}

  xa = xi(tree)
  xo = xf(tree)
  x1 = xf(tre1)
  ea = e(tree)
  e1 = e(tre1)

  # gibbs sampling
  xn = duoprop(xa, x1, ea, e1, σx)

  setxf!(tree, xn)
  setxi!(tre1, xn)

  # update llc, prc and sdX
  llc += duoldnorm(xn, xa, x1, ea, e1, σx) -
         duoldnorm(xo, xa, x1, ea, e1, σx)

  sdX += ((xa - xn)^2 - (xa - xo)^2)/(2.0*ea) +
         ((xn - x1)^2 - (xo - x1)^2)/(2.0*e1)

  return llc, sdX
end





