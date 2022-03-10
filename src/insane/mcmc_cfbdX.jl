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
  stem = !iszero(e(tree))

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
  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 500, stem)

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
        λc, μc, ψc, σxc, mc, th, stem, inodes, pup, prints)

  # mcmc
  r, treev =
    mcmc_cfbd(Ξ, idf, llc, prc, λc, μc, ψc, σxc, mc, th, stem,
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
                  stem   ::Bool,
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
                        stem    ::Bool,
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
  llc = llik_cfbd(Ξ, λc, μc, ψc, σxc) - !stem*log(λc) + log(mc) + prob_ρ(idf)
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
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, stem, λ_prior)

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

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
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, stem, x0_prior)

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
             stem   ::Bool,
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
                   stem    ::Bool,
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
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, stem, λ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - !stem*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, stem, μ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - !stem*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

       # psi update
      elseif p === 3

        llc, prc, ψc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - !stem*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

       # sigma_x update
      elseif p === 4

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - !stem*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # X ancestors update
      elseif p === 5

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX =
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, stem, x0_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - !stem*log(λc) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L, sdX, nX =
          update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, σxc, ns, ne, L, sdX, nX)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, σxc) - !stem*log(λc) + log(mc) + prob_ρ(idf)
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

Forward simulation proposal function for constant birth-death.
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

  # forward simulate an internal branch
  ξp, np, ntp, xt = fsbi_cfbdX(bi, λ, μ, ψ, xi(Ξ[bix]), σx, 100)

  # retained conditional on survival
  if ntp > 0

    itb = it(bi) # is it terminal
    iψb = isfossil(bi) # is it a fossil
    ρbi = ρi(bi) # get branch sampling fraction
    nc  = ni(bi) # current ni
    ntc = nt(bi) # current nt
    fxi = fx(bi) # is a fix X node

    # current tree
    ξc  = Ξ[bix]

    # if terminal non-fossil branch
    if itb && !iψb
      llr = log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc))
      xt  = fixed_xt(ξc)       # get previous `x` of fixed tip
      if fxi
        acr = _match_tip_x!(ξp, xt, σx)
      else
        acr = 0.0
      end
    else
      np -= !iψb
      llr = log((1.0 - ρbi)^(np - nc))
      acr = log(Float64(ntp)/Float64(ntc))

      if itb
        if fxi
          xt   = fossil_xt(ξc)       # get previous `x` of fixed tip
          acr += _match_fossil_x!(ξp, xt, σx)
        end
      else
        if fxi
          acr += _match_tip_x!(ξp, xt, σx)
        else
          ξ1  = Ξ[d1(bi)]
          if iψb
            llr += lrdnorm_bm_x(xt, xi(ξ1), xf(ξ1), sqrt(e(ξ1))*σx)
          else
            ξ2  = Ξ[d2(bi)]
            llr += duoldnorm(xt,     xf(ξ1), xf(ξ2), e(ξ1), e(ξ2), σx) -
                   duoldnorm(xi(ξ1), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2), σx)
          end
        end
      end
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
      llr += llik_cfbd(ξp, λ, μ, ψ, σx) - llik_cfbd(ξc, λ, μ, ψ, σx)

      Ξ[bix] = ξp     # set new decoupled tree
      llc += llr      # set new likelihood
      setni!(bi, np)  # set new ni
      setnt!(bi, ntp) # set new nt

      if !fxi && !itb
        sdX += ((xt - xf(ξ1))^2 - (xi(ξ1) - xf(ξ1))^2)/(2.0*e(ξ1))
        setxi!(ξ1, xt) # set new xt for initial x
        if !iψb
          sdX += ((xt - xf(ξ2))^2 - (xi(ξ2) - xf(ξ2))^2)/(2.0*e(ξ2))
          setxi!(ξ2, xt) # set new xt for initial x
        end
      end
    end
  end

  return llc, ns, ne, L, sdX, nX
end




"""
    fsbi_cfbdX(bi  ::iBffs,
               λ   ::Float64,
               μ   ::Float64,
               ψ   ::Float64,
               x0  ::Float64,
               σx  ::Float64,
               ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi_cfbdX(bi  ::iBffs,
                    λ   ::Float64,
                    μ   ::Float64,
                    ψ   ::Float64,
                    x0  ::Float64,
                    σx  ::Float64,
                    ntry::Int64)

  # times
  tfb = tf(bi)

  ext = 0
  # condition on non-extinction (helps in mixing)
  while ext < ntry
    ext += 1

    # forward simulation during branch length
    t0, na, nf = sim_cfbd(e(bi), λ, μ, ψ, x0, σx, 0, 0)

    if na > 0 && iszero(nf) # exclude if any fossil is sampled
      nat = na

      if isone(na)
        # fix the only tip alive
        f, xt = fixalive!(t0, NaN)

      elseif na > 1
        # fix random tip
        xt = fixrtip!(t0, na, NaN)

        if !it(bi) || isfossil(bi)
          # add tips until the present
          tx, na, nf = tip_sims!(t0, tfb, λ, μ, ψ, σx, na, nf)
          if !iszero(nf)
            ext += 1
            continue
          end
        end
      end

      if isfossil(bi)
        # replace extant tip by a fossil
        fossilizefixedtip!(t0)

        # if terminal fossil branch
        if it(bi)
          tx, na, nf = fossiltip_sim!(t0, tfb, λ, μ, ψ, σx, na, nf)
          if !iszero(nf)
            ext += 1
            continue
          end
        end
      end

      return t0, na, nat, xt
    end
  end

  return sTfbdX(0.0, false, false, false, 0.0, 0.0), 0, 0, NaN
end





"""
    tip_sims!(tree::sTfbdX,
              t   ::Float64,
              λ   ::Float64,
              μ   ::Float64,
              ψ   ::Float64,
              σx  ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTfbdX,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   σx  ::Float64,
                   na  ::Int64,
                   nf  ::Int64)

  if iszero(nf)
    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nf = sim_cfbd(t, λ, μ, ψ, xf(tree), σx, na-1, nf)

        if iszero(nf)
          # merge to current tip
          sete!(tree, e(tree) + e(stree))
          setproperty!(tree, :iμ, isextinct(stree))
          setxf!(tree, xf(stree))
          if def1(stree)
            tree.d1 = stree.d1
            tree.d2 = stree.d2
          end
        end
      end
    else
      tree.d1, na, nf = tip_sims!(tree.d1, t, λ, μ, ψ, σx, na, nf)
      tree.d2, na, nf = tip_sims!(tree.d2, t, λ, μ, ψ, σx, na, nf)
    end
  end

  return tree, na, nf
end




"""
    fossiltip_sim!(tree::sTfbdX,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   σx  ::Float64,
                   na  ::Int64,
                   nf  ::Int64)

Continue simulation until time `t` for the fixed tip in `tree`.
"""
function fossiltip_sim!(tree::sTfbdX,
                        t   ::Float64,
                        λ   ::Float64,
                        μ   ::Float64,
                        ψ   ::Float64,
                        σx  ::Float64,
                        na  ::Int64,
                        nf  ::Int64)

  if iszero(nf)
    if istip(tree)
      stree, na, nf = sim_cfbd(t, λ, μ, ψ, xf(tree), σx, na-1, 0)
      if iszero(nf)
        # merge to current tip
        tree.d1 = stree
      end
    elseif isfix(tree.d1)
      tree.d1, na, nf = fossiltip_sim!(tree.d1, t, λ, μ, ψ, σx, na, nf)
    else
      tree.d2, na, nf = fossiltip_sim!(tree.d2, t, λ, μ, ψ, σx, na, nf)
    end
  end

  return tree, na, nf
end




"""
    _match_fossil_x!(tree::T,
                     xt  ::Float64,
                     σx  ::Float64) where {T <: sTX}

Make joint proposal to match simulation with tip fixed `x` value.
"""
function _match_fossil_x!(tree::T,
                          xt  ::Float64,
                          σx  ::Float64) where {T <: sTX}

  if isfossil(tree)

    xa = xi(tree)
    xn = xf(tree)
    x1 = xf(tree.d1)
    ea = e(tree)
    e1 = e(tree.d1)
    acr = duoldnorm(xt, xa, x1, ea, e1, σx) -
          duoldnorm(xn, xa, x1, ea, e1, σx)
    setxf!(tree,    xt)
    setxi!(tree.d1, xt)

    return acr
  else
    if isfix(tree.d1)
      acr = _match_fossil_x!(tree.d1, xt, σx)
    else
      acr = _match_fossil_x!(tree.d2, xt, σx)
    end
  end

  return acr
end




"""
    update_x!(bix     ::Int64,
              Ξ       ::Vector{sTfbdX},
              idf     ::Vector{iBffs},
              σx      ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              sdX     ::Float64,
              stem    ::Bool,
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
                   stem    ::Bool,
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
  if root && !stem
    if !fx(bi)
      llc, prc, sdX =
         _crown_update_x!(ξi, ξ1, ξ2, σx, llc, prc, sdX, x0_prior)
    end
  else
    # if stem
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




"""
    _update_tip_x!(tree::T,
                   σx  ::Float64,
                   llc ::Float64,
                   sdX ::Float64) where {T <: sTX}

Make gibbs node update for trait.
"""
function _update_tip_x!(tree::T,
                        σx  ::Float64,
                        llc ::Float64,
                        sdX ::Float64) where {T <: sTX}

  xa = xi(tree)
  xo = xf(tree)
  ea = e(tree)

  # gibbs sampling
  s = sqrt(ea)*σx
  xn = rnorm(xa, s)
  setxf!(tree, xn)

  # update llc and sdX
  llc += llrdnorm_x(xn, xo, xa, s^2)
  sdX += ((xn - xa)^2 - (xo - xa)^2)/(2.0*ea)

  return llc, sdX
end



