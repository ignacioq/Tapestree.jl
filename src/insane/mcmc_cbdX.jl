#=

constant birth-death MCMC using forward simulation with traits

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cbd(tree    ::sT_label,
               X       ::Dict{String, Float64},
               out_file::String;
               λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
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
function insane_cbd(tree    ::sT_label,
                    X       ::Dict{String, Float64},
                    out_file::String;
                    λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
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
                    pupdp   ::NTuple{5,Float64}     = (0.2,0.2,0.2,0.2,0.2),
                    prints  ::Int64                 = 5,
                    tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)
  crown = Int64(iszero(e(tree)))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxc = make_idf(tree, tρ, X)

  # starting parameters
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), ti(idf[1]), ϵi)
  else
    λc, μc = λi, μi
  end
  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 1_000, crown)

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, xr, sTbdX)

  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(Ξ))
    !it(idf[i]) && push!(inodes, i)
  end

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant birth-death and trait evolution"

  # adaptive phase
  llc, prc, λc, μc, σxc, mc =
      mcmc_burn_cbd(Ξ, idf, λ_prior, μ_prior, σx_prior, x0_prior, nburn,
        λc, μc, σxc, mc, th, crown, inodes, pup, prints)

  # mcmc
  r, treev =
    mcmc_cbd(Ξ, idf, llc, prc, λc, μc, σxc, mc, th, crown,
      λ_prior, μ_prior, σx_prior, x0_prior, niter, nthin, inodes, pup, prints)

  pardic = Dict("lambda"  => 1,
                "mu"      => 2,
                "x0"      => 3,
                "sigma_x" => 4)

  write_ssr(r, pardic, out_file)

  return r, treev, NaN
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
                  crown   ::Bool,
                  pup    ::Array{Int64,1},
                  prints ::Int64)

Adaptive MCMC phase for da chain for constant birth-death using forward
simulation.
"""
function mcmc_burn_cbd(Ξ       ::Vector{sTbdX},
                       idf     ::Array{iBffs,1},
                       λ_prior ::NTuple{2,Float64},
                       μ_prior ::NTuple{2,Float64},
                       σx_prior::NTuple{2,Float64},
                       x0_prior::NTuple{2,Float64},
                       nburn   ::Int64,
                       λc      ::Float64,
                       μc      ::Float64,
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
  ns      = Float64(el-1)/2.0    # number of speciation events
  ne      = 0.0                  # number of extinction events
  sdX, nX = sdeltaX(Ξ)           # standardized trait differences

  # likelihood
  llc = llik_cbd(Ξ, λc, μc, σxc) - Float64(crown) * log(λc) + log(mc) + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2])           +
        logdgamma(μc, μ_prior[1], μ_prior[2])           +
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

       # sigma_x update
      elseif p === 3

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

      # X ancestors update
      elseif p === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX =
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, x0_prior)

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L, sdX, nX =
          update_fs!(bix, Ξ, idf, llc, λc, μc, σxc, ns, ne, L, sdX, nX)

      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, σxc, mc
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
             crown   ::Bool,
             λ_prior::NTuple{2,Float64},
             μ_prior::NTuple{2,Float64},
             niter  ::Int64,
             nthin  ::Int64,
             pup    ::Array{Int64,1},
             prints ::Int64)

MCMC da chain for constant birth-death using forward simulation.
"""
function mcmc_cbd(Ξ       ::Vector{sTbdX},
                  idf     ::Array{iBffs,1},
                  llc     ::Float64,
                  prc     ::Float64,
                  λc      ::Float64,
                  μc      ::Float64,
                  σxc     ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  crown    ::Int64,
                  λ_prior ::NTuple{2,Float64},
                  μ_prior ::NTuple{2,Float64},
                  σx_prior::NTuple{2,Float64},
                  x0_prior::NTuple{2,Float64},
                  niter   ::Int64,
                  nthin   ::Int64,
                  inodes  ::Array{Int64,1},
                  pup     ::Array{Int64,1},
                  prints  ::Int64)

  el      = lastindex(idf)
  nin     = lastindex(inodes)
  ns      = Float64(nnodesinternal(Ξ))
  ne      = Float64(ntipsextinct(Ξ))
  L       = treelength(Ξ)
  sdX, nX = sdeltaX(Ξ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 7)

  # make tree vector
  treev  = sTbdX[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, crown, λ_prior)

        llci = llik_cbd(Ξ, λc, μc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        if !isapprox(llci, llc, atol = 1e-6)
           @show llci, llc, it, p
           return
        end

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, crown, μ_prior)

        llci = llik_cbd(Ξ, λc, μc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        if !isapprox(llci, llc, atol = 1e-6)
           @show llci, llc, it, p
           return
        end

       # sigma_x update
      elseif p === 3

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

        llci = llik_cbd(Ξ, λc, μc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        if !isapprox(llci, llc, atol = 1e-6)
           @show llci, llc, it, p
           return
        end

      # X ancestors update
      elseif p === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX =
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, x0_prior)

        llci = llik_cbd(Ξ, λc, μc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        if !isapprox(llci, llc, atol = 1e-6)
           @show llci, llc, it, p
           return
        end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L, sdX, nX =
          update_fs!(bix, Ξ, idf, llc, λc, μc, σxc, ns, ne, L, sdX, nX)

        llci = llik_cbd(Ξ, λc, μc, σxc) - Float64(crown)*log(λc) + log(mc) + prob_ρ(idf)
        if !isapprox(llci, llc, atol = 1e-6)
           @show llci, llc, it, p
           return
        end

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
               Ξ  ::Vector{sTbdX},
               idf::Vector{iBffs},
               llc::Float64,
               λ  ::Float64,
               σx ::Float64,
               ns ::Float64,
               μ  ::Float64,
               ne ::Float64,
               L  ::Float64,
               sdX::Float64,
               nX ::Float64)

Forward simulation proposal function for constant pure-birth.
"""
function update_fs!(bix::Int64,
                    Ξ  ::Vector{sTbdX},
                    idf::Vector{iBffs},
                    llc::Float64,
                    λ  ::Float64,
                    μ  ::Float64,
                    σx ::Float64,
                    ns ::Float64,
                    ne ::Float64,
                    L  ::Float64,
                    sdX::Float64,
                    nX ::Float64)

  bi = idf[bix]
  ξc  = Ξ[bix]

  if it(bi) # is it terminal
    ξp, llr = fsbi_t(bi, ξc, λ, μ, σx)
    sdXr = 0.0
  else
    ξp, llr, sdXr = 
      fsbi_i(bi, Ξ[d1(bi)], Ξ[d2(bi)], xi(ξc), λ, μ, σx)
  end

  if isfinite(llr)

    # update llc, ns & L
    llc += llr + llik_cbd(ξp, λ, μ, σx) - llik_cbd(ξc, λ, μ, σx)
    ns  += Float64(nnodesinternal(ξp) - nnodesinternal(ξc))
    ne  += Float64(ntipsextinct(ξp)   - ntipsextinct(ξc))
    L   += treelength(ξp)             - treelength(ξc)
    sdXp, nXp = _sdeltaX(ξp, 0.0, 0.0)
    sdXc, nXc = _sdeltaX(ξc, 0.0, 0.0)
    sdX += sdXp - sdXc + sdXr
    nX  += nXp  - nXc

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ns, ne, L, sdX, nX
end




"""
    fsbi_t(bi::iBffs,
           ξc::sTbdX,
           λ ::Float64,
           μ ::Float64,
           σx::Float64)

Forward simulation for terminal branch `bi`.
"""
function fsbi_t(bi::iBffs,
                ξc::sTbdX,
                λ ::Float64,
                μ ::Float64,
                σx::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  xist = Float64[]
  xfst = Float64[]
  est  = Float64[]

  t0, na, nn, llr =
    _sim_cbd_t(e(bi), λ, μ, xi(ξc), σx, lc, lU, Iρi, 0, 1, 500, xist, xfst, est)

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

    if lU <  acr + llr
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
  end

  return t0, NaN
end




"""
    fsbi_i(bi::iBffs,
           ξc::sTbdX,
           λ ::Float64,
           σx::Float64)

Forward simulation for terminal branch `bi`.
"""
function fsbi_i(bi::iBffs,
                ξ1::sTbdX,
                ξ2::sTbdX,
                x0::Float64,
                λ ::Float64,
                μ ::Float64,
                σx::Float64)

  t0, na, nn = _sim_cbd_i(e(bi), λ, μ, x0, σx, 0, 1, 500)

  if na < 1 || nn >= 500
    return t0, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))

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

    if na > 1
      # simulated remaining tips until the present
      tx, na, nn, acr = tip_sims!(t0, tf(bi), λ, μ, σx, acr, lU, Iρi, na, nn)
    end

    if lU < acr + llr
      na -= 1
      llr += (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setni!(bi, na)  # set new ni
      setnt!(bi, ntp) # set new nt

      sdX = ((xp - xf(ξ1))^2 - (xi(ξ1) - xf(ξ1))^2)/(2.0*e(ξ1)) +
            ((xp - xf(ξ2))^2 - (xi(ξ2) - xf(ξ2))^2)/(2.0*e(ξ2))
      setxi!(ξ1, xp) # set new xp for initial x
      setxi!(ξ2, xp) # set new xp for initial x

      return t0, llr, sdX
    end

  end

  return t0, NaN, NaN
end




"""
    tip_sims!(tree::sTbdX,
              t   ::Float64,
              λ   ::Float64,
              μ   ::Float64,
              σx  ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTbdX,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   σx  ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

  if na < 500 && lU < lr

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nn, lr = 
          _sim_cbd_it(t, λ, μ, xf(tree), σx, lr, lU, Iρi, na-1, nn, 500)

        if isnan(lr) || nn >= 500
          return tree, na, nn, NaN
        end

        # merge to current tip
        sete!( tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        setxf!(tree, xf(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr = tip_sims!(tree.d1, t, λ, μ, σx, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr = tip_sims!(tree.d2, t, λ, μ, σx, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end



