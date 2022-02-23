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
  stem = !iszero(e(tree))

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
  mc = m_surv_cbd(th, λc, μc, 1_000, stem)

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

  @info "Running constant birth-death and trait evolution with forward simulation"

  # adaptive phase
  llc, prc, λc, μc, σxc, mc = 
      mcmc_burn_cbd(Ξ, idf, λ_prior, μ_prior, σx_prior, x0_prior, nburn, 
        λc, μc, σxc, mc, th, stem, inodes, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cbd(Ξ, idf, llc, prc, λc, μc, σxc, mc, th, stem,
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
                  stem   ::Bool,
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
                       stem    ::Bool,
                       inodes  ::Array{Int64,1},
                       pup     ::Array{Int64,1}, 
                       prints  ::Int64)

  el      = lastindex(idf)
  nin = lastindex(inodes)
  L       = treelength(Ξ)        # tree length
  ns      = Float64(el-1)/2.0    # number of speciation events
  ne      = 0.0                  # number of extinction events
  sdX, nX = sdeltaX(Ξ)           # standardized trait differences
  nsi     = stem ? 0.0 : log(λc) # if stem or crown

  # likelihood
  llc = llik_cbd(Ξ, λc, μc, σxc) - nsi + log(mc) + prob_ρ(idf)
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
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, stem, λ_prior)

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc = 
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

       # sigma_x update
      elseif p === 3

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

      # X ancestors update
      elseif p === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX =
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, stem, x0_prior)

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
             stem   ::Bool,
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
                  stem    ::Bool,
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
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, stem, λ_prior)

        # llci = llik_cbd(Ξ, λc, μc, σxc) - log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc = 
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

        # llci = llik_cbd(Ξ, λc, μc, σxc) - log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

       # sigma_x update
      elseif p === 3

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

        # llci = llik_cbd(Ξ, λc, μc, σxc) - log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # X ancestors update
      elseif p === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX =
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, stem, x0_prior)

        # llci = llik_cbd(Ξ, λc, μc, σxc) - log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end


      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L, sdX, nX = 
          update_fs!(bix, Ξ, idf, llc, λc, μc, σxc, ns, ne, L, sdX, nX)

        # llci = llik_cbd(Ξ, λc, μc, σxc) - log(λc) + log(mc) + prob_ρ(idf)
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
               Ξ  ::Vector{sTbdX},
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

  # forward simulate an internal branch
  ξp, np, ntp, xt = fsbi_cbdX(bi, λ, μ, xi(Ξ[bix]), σx, 100)

  # retained conditional on survival
  if ntp > 0

    itb = it(bi) # is it terminal
    ρbi = ρi(bi) # get branch sampling fraction
    nc  = ni(bi) # current ni
    ntc = nt(bi) # current nt

    # current tree
    ξc  = Ξ[bix]

    # if terminal branch
    if itb
      llr = log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc))
      xt  = fixed_xt(ξc)       # get previous `x` of fixed tip
      acr = _match_tip_x!(ξp, xt, σx)
    else
      np -= 1
      ξ1  = Ξ[d1(bi)]
      ξ2  = Ξ[d2(bi)]
      llr = log((1.0 - ρbi)^(np - nc))                          +
            duoldnorm(xt,     xf(ξ1), xf(ξ2), e(ξ1), e(ξ2), σx) -
            duoldnorm(xi(ξ1), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2), σx)
      acr = log(Float64(ntp)/Float64(ntc))
    end

    # MH ratio
    if -randexp() < llr + acr

      # update ns, ne & L
      ns += Float64(nnodesinternal(ξp) - nnodesinternal(ξc))
      ne += Float64(ntipsextinct(ξp)   - ntipsextinct(ξc))
      L  += treelength(ξp)             - treelength(ξc)
      sdXp, nXp = _sdeltaX(ξp, 0.0, 0.0) 
      sdXc, nXc = _sdeltaX(ξc, 0.0, 0.0) 
      sdX += sdXp - sdXc
      nX  += nXp  - nXc

      # likelihood ratio
      llr += llik_cbd(ξp, λ, μ, σx) - llik_cbd(ξc, λ, μ, σx)

      Ξ[bix] = ξp     # set new decoupled tree
      llc += llr      # set new likelihood
      setni!(bi, np)  # set new ni
      setnt!(bi, ntp) # set new nt

      if !itb
        sdX += ((xt - xf(ξ1))^2 - (xi(ξ1) - xf(ξ1))^2)/(2.0*e(ξ1)) +
               ((xt - xf(ξ2))^2 - (xi(ξ1) - xf(ξ2))^2)/(2.0*e(ξ2))
        setxi!(ξ1, xt) # set new xt
        setxi!(ξ2, xt) # set new xt
      end

    end
  end

  return llc, ns, ne, L, sdX, nX
end




"""
    fsbi(bi::iBffs, λ::Float64, μ::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi_cbdX(bi::iBffs, 
                   λ ::Float64, 
                   μ ::Float64, 
                   x0::Float64, 
                   σx::Float64, 
                   ntry::Int64)

  # times
  tfb = tf(bi)

  ext = 0
  # condition on non-extinction (helps in mixing)
  while ext < ntry 
    ext += 1

    # forward simulation during branch length
    t0, na = sim_cbd(e(bi), λ, μ, x0, σx, 0)

    nat = na

    if isone(na)
     f, xt = fixalive!(t0, NaN)

      return t0, na, nat, xt
    elseif na > 1
      # fix random tip
      xt = fixrtip!(t0, na, NaN)

      if !it(bi)
        # add tips until the present
        tx, na = tip_sims!(t0, tfb, λ, μ, σx, na)
      end

      return t0, na, nat, xt
    end
  end

  return sTbdX(0.0, false, false, 0.0, 0.0), 0, 0, NaN
end




"""
    tip_sims!(tree::sTbdX, t::Float64, λ::Float64, μ::Float64, na::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::sTbdX, 
                   t   ::Float64, 
                   λ   ::Float64, 
                   μ   ::Float64, 
                   σx  ::Float64, 
                   na  ::Int64)

  if istip(tree) 
    if !isfix(tree) && isalive(tree)

      # simulate
      stree, na = sim_cbd(t, λ, μ, xf(tree), σx, na-1)

      # merge to current tip
      sete!(tree, e(tree) + e(stree))
      setproperty!(tree, :iμ, isextinct(stree))
      setxf!(tree, xf(stree))
      if isdefined(stree, :d1)
        tree.d1 = stree.d1
        tree.d2 = stree.d2
      end
    end
  else
    tree.d1, na = tip_sims!(tree.d1, t, λ, μ, σx, na)
    tree.d2, na = tip_sims!(tree.d2, t, λ, μ, σx, na)
  end

  return tree, na
end




