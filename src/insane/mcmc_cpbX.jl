#=

pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    insane_cpb(tree    ::sT_label, 
               X       ::Dict{String, Float64},
               out_file::String;
               λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               niter   ::Int64                 = 1_000,
               nthin   ::Int64                 = 10,
               nburn   ::Int64                 = 200,
               tune_int::Int64                 = 100,
               marginal::Bool                  = false,
               nitpp   ::Int64                 = 100, 
               nthpp   ::Int64                 = 10,
               K       ::Int64                 = 10,
               λi      ::Float64               = NaN,
               pupdp   ::NTuple{2,Float64}     = (0.2, 0.2),
               prints  ::Int64                 = 5,
               tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant pure-birth with Brownian motion trait evolution.
"""
function insane_cpb(tree    ::sT_label, 
                    X       ::Dict{String, Float64},
                    out_file::String;
                    λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    σx_prior::NTuple{2,Float64}     = (0.05, 0.05),
                    x0_prior::NTuple{2,Float64}     = (0.0, 10.0),
                    niter   ::Int64                 = 1_000,
                    nthin   ::Int64                 = 10,
                    nburn   ::Int64                 = 200,
                    tune_int::Int64                 = 100,
                    marginal::Bool                  = false,
                    nitpp   ::Int64                 = 100, 
                    nthpp   ::Int64                 = 10,
                    K       ::Int64                 = 10,
                    λi      ::Float64               = NaN,
                    pupdp   ::NTuple{4,Float64}     = (0.2, 0.2, 0.3, 0.3),
                    prints  ::Int64                 = 5,
                    tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  stem = !iszero(e(tree))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory and pic reconstruction
  idf, xr, σxc = make_idf(tree, tρ, X)

  # starting parameters
  if isnan(λi)
    λc = Float64(n-2)/treelength(tree)
  else
    λc = λi
  end

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, xr, sTpbX)

  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(Ξ))
    !it(idf[i]) && push!(inodes, i)
  end

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant pure-birth and BM trait evolution with forward simulation"

  # adaptive phase
  llc, prc, λc, σxc, sdX, nX = 
    mcmc_burn_cpb(Ξ, idf, λ_prior, σx_prior, x0_prior, 
      nburn, λc, σxc, pup, inodes, prints, stem)

  # mcmc
  r, treev, λc, σxc = 
    mcmc_cpb(Ξ, idf, llc, prc, λc, σxc, sdX, nX, λ_prior, σx_prior, x0_prior, 
     niter, nthin, pup, inodes, prints, stem)

  pardic = Dict(("lambda"  => 1,
                 "x0"      => 2,
                 "sigma_x" => 3))

  write_ssr(r, pardic, out_file)

  if marginal
    # reference distribution
    βs = [range(0.0, 1.0, K)...]::Vector{Float64}
    reverse!(βs)

    @views p = r[:,4]

    # make reference posterior
    m     = mean(p)
    v     = var(p)
    λ_refd = (m^2/v, m/v)

    # marginal likelihood
    pp = ref_posterior(Ξ, idf, λc, λ_prior, λ_refd, nitpp, nthpp, βs, pup, stem)

    # process with reference distribution the posterior
    p1 = Vector{Float64}(undef, size(r,1))
    for i in Base.OneTo(size(r,1))
      p1[i] = r[i,2] + r[i,3] - logdgamma(r[i,4], λ_refd[1], λ_refd[2])
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
    mcmc_burn_cpb(Ξ      ::Vector{sTpbX}, 
                  idf    ::Array{iBffs,1},
                  λ_prior::NTuple{2,Float64},
                  nburn  ::Int64,
                  λc     ::Float64,
                  pup    ::Array{Int64,1}, 
                  prints ::Int64,
                  stem   ::Bool)

MCMC chain for constant pure-birth.
"""
function mcmc_burn_cpb(Ξ       ::Vector{sTpbX}, 
                       idf     ::Array{iBffs,1},
                       λ_prior ::NTuple{2,Float64},
                       σx_prior::NTuple{2,Float64},
                       x0_prior::NTuple{2,Float64},
                       nburn   ::Int64,
                       λc      ::Float64,
                       σxc     ::Float64,
                       pup     ::Array{Int64,1}, 
                       inodes  ::Array{Int64,1},
                       prints  ::Int64,
                       stem    ::Bool)

  el      = lastindex(idf)
  nin     = lastindex(inodes)
  L       = treelength(Ξ)        # tree length
  ns      = Float64(el-1)*0.5    # number of speciation events
  sdX, nX = sdeltaX(Ξ)           # standardized trait differences
  nsi     = stem ? 0.0 : log(λc) # if stem or crown

  #likelihood
  llc = llik_cpb(Ξ, λc, σxc) - nsi + prob_ρ(idf)
  prc = logdgamma(λc,        λ_prior[1], λ_prior[2])    + 
        logdinvgamma(σxc^2, σx_prior[1], σx_prior[2])   +
        logdnorm(xi(Ξ[1]),  x0_prior[1], x0_prior[2]^2)

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc = update_λ!(llc, prc, λc, ns, L, stem, λ_prior)

      # sigma_x update
      elseif p === 2

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

      # X ancestors update
      elseif p === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX =
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, stem, x0_prior)

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, L, sdX, nX = 
          update_fs!(bix, Ξ, idf, llc, λc, σxc, ns, L, sdX, nX)

      end
    end

    next!(pbar)
  end

  return llc, prc, λc, σxc, sdX, nX
end




"""
    mcmc_cpb(Ξ       ::Vector{sTpbX},
             idf     ::Array{iBffs,1},
             llc     ::Float64,
             prc     ::Float64,
             λc      ::Float64,
             σxc     ::Float64,
             sdX     ::Float64, 
             nX      ::Float64,
             λ_prior ::NTuple{2,Float64},
             σx_prior::NTuple{2,Float64},
             x0_prior::NTuple{2,Float64},
             niter   ::Int64,
             nthin   ::Int64,
             pup     ::Array{Int64,1}, 
             prints  ::Int64,
             stem    ::Bool)

MCMC chain for constant pure-birth.
"""
function mcmc_cpb(Ξ       ::Vector{sTpbX},
                  idf     ::Array{iBffs,1},
                  llc     ::Float64,
                  prc     ::Float64,
                  λc      ::Float64,
                  σxc     ::Float64,
                  sdX     ::Float64, 
                  nX      ::Float64,
                  λ_prior ::NTuple{2,Float64},
                  σx_prior::NTuple{2,Float64},
                  x0_prior::NTuple{2,Float64},
                  niter   ::Int64,
                  nthin   ::Int64,
                  pup     ::Array{Int64,1}, 
                  inodes  ::Array{Int64,1},
                  prints  ::Int64,
                  stem    ::Bool)

  el  = lastindex(idf)
  nin = lastindex(inodes)
  ns  = Float64(nnodesinternal(Ξ))
  L   = treelength(Ξ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  r = Array{Float64,2}(undef, nlogs, 6)

  # make tree vector
  treev  = sTpbX[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc = update_λ!(llc, prc, λc, ns, L, stem, λ_prior)

        # llci = llik_cpb(Ξ, λc, σxc) - log(λc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # sigma_x update
      elseif p === 2

        llc, prc, σxc = update_σx!(σxc, sdX, nX, llc, prc, σx_prior)

        # llci = llik_cpb(Ξ, λc, σxc) - log(λc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # X ancestors update
      elseif p === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, sdX = 
          update_x!(bix, Ξ, idf, σxc, llc, prc, sdX, stem, x0_prior)

        # llci = llik_cpb(Ξ, λc, σxc) - log(λc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, L, sdX, nX = 
          update_fs!(bix, Ξ, idf, llc, λc, σxc, ns, L, sdX, nX)

        # llci = llik_cpb(Ξ, λc, σxc) - log(λc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      end
    end

    lthin += 1
    if lthin == nthin
      lit += 1
      @inbounds begin
        r[lit,1] = Float64(lit)
        r[lit,2] = llc
        r[lit,3] = prc
        r[lit,4] = λc
        r[lit,5] = xi(Ξ[1])
        r[lit,6] = σxc
        push!(treev, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return r, treev, λc, σxc
end






"""
    update_fs!(bix::Int64,
               Ξ  ::Vector{sTpbX},
               idf::Vector{iBffs},
               llc::Float64,
               λ  ::Float64,
               σx ::Float64,
               ns ::Float64,
               L  ::Float64,
               sdX::Float64,
               nX ::Float64)

Forward simulation proposal function for constant pure-birth.
"""
function update_fs!(bix::Int64,
                    Ξ  ::Vector{sTpbX},
                    idf::Vector{iBffs},
                    llc::Float64,
                    λ  ::Float64,
                    σx ::Float64,
                    ns ::Float64,
                    L  ::Float64,
                    sdX::Float64,
                    nX ::Float64)

  bi = idf[bix]

  # forward simulate an internal branch
  ξp, np, ntp, xt = fsbi_pbX(bi, λ, xi(Ξ[bix]), σx, 1_000)

  itb = it(bi) # is it terminal
  ρbi = ρi(bi) # get branch sampling fraction
  nc  = ni(bi) # current ni
  ntc = nt(bi) # current nt

  if ntp > 0

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

      # update ns, L & sdX
      ns  += Float64(nnodesinternal(ξp) - nnodesinternal(ξc))
      L   += treelength(ξp)             - treelength(ξc)
      sdXp, nXp = _sdeltaX(ξp, 0.0, 0.0) 
      sdXc, nXc = _sdeltaX(ξc, 0.0, 0.0) 
      sdX += sdXp - sdXc
      nX  += nXp  - nXc

      # add likelihood ratio
      llr += llik_cpb(ξp, λ, σx) - llik_cpb(ξc, λ, σx)

      Ξ[bix] = ξp     # set new decoupled tree
      llc   += llr    # set new likelihood
      setni!(bi, np)  # set new ni
      setnt!(bi, ntp) # set new nt

      if !it(bi)
        sdX += ((xt - xf(ξ1))^2 - (xi(ξ1) - xf(ξ1))^2)/(2.0*e(ξ1)) +
               ((xt - xf(ξ2))^2 - (xi(ξ1) - xf(ξ2))^2)/(2.0*e(ξ2))
        setxi!(ξ1, xt) # set new xt
        setxi!(ξ2, xt) # set new xt
      end

    end
  end

  return llc, ns, L, sdX, nX
end




"""
    fsbi_pbX(bi::iBffs, λ::Float64, x0::Float64, σx::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi_pbX(bi::iBffs, λ::Float64, x0::Float64, σx::Float64, ntry::Int64)

  # times
  tfb = tf(bi)

  ext = 0
  # condition on non-extinction (helps in mixing)
  while ext < ntry 
    ext += 1

    # forward simulation during branch length
    t0, na = sim_cpb(e(bi), λ, x0, σx, 0)

    nat = na

    if isone(na)
      f, xt = fixalive!(t0, NaN)

      return t0, na, nat, xt
    elseif na > 1
      # fix random tip
      xt = fixrtip!(t0, 0, NaN)

      if !it(bi)
        # add tips until the present
        tx, na = tip_sims!(t0, tfb, λ, σx, na)
      end

      return t0, na, nat, xt
    end
  end

  return sTpbX(0.0, false, 0.0, 0.0), 0, 0, NaN
end




"""
    tip_sims!(tree::sTpbX, t::Float64, λ::Float64, σx::Float64, na::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::sTpbX, t::Float64, λ::Float64, σx::Float64, na::Int64)

  if istip(tree) 
    if !isfix(tree)

      # simulate
      stree, na = sim_cpb(t, λ, xf(tree), σx, na-1)

      # merge to current tip
      sete!( tree, e(tree) + e(stree))
      setxf!(tree, xf(stree))
      if isdefined(stree, :d1)
        tree.d1 = stree.d1
        tree.d2 = stree.d2
      end
    end
  else
    tree.d1, na = tip_sims!(tree.d1, t, λ, σx, na)
    tree.d2, na = tip_sims!(tree.d2, t, λ, σx, na)
  end

  return tree, na
end




"""
    _match_tip_x!(tree::T,
                 xt  ::Float64,
                 σx  ::Float64) where {T <: iTreeX}

Make joint proposal to match simulation with tip fixed `x` value.
"""
function _match_tip_x!(tree::T,
                       xt  ::Float64,
                       σx  ::Float64) where {T <: iTreeX}

  if istip(tree)
    xa   = xi(tree)
    xn   = xf(tree)
    el   = e(tree)
    σsrt = sqrt(el) * σx
    # acceptance rate
    acr  = ldnorm_bm(xt, xa, σsrt) - ldnorm_bm(xn, xa, σsrt)
    setxf!(tree, xt)
  else
    if isfix(tree.d1)
      acr = _match_tip_x!(tree.d1, xt, σx)
    else
      acr = _match_tip_x!(tree.d2, xt, σx)
    end
  end

  return acr
end




"""
    update_x!(bix     ::Int64,
              Ξ       ::Vector{sTpbX},
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
                   Ξ       ::Vector{sTpbX},
                   idf     ::Vector{iBffs},
                   σx      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   sdX     ::Float64,
                   stem    ::Bool,
                   x0_prior::NTuple{2, Float64})

  ξi = Ξ[bix]
  bi = idf[bix]
  ξ1 = Ξ[d1(bi)]
  ξ2 = Ξ[d2(bi)]

  root = iszero(pa(bi))
  # if crown root
  if root && !stem
    llc, prc, sdX = 
       _crown_update_x!(ξi, ξ1, ξ2, σx, llc, prc, sdX, x0_prior)
  else
    # if stem
    if root
      llc, prc, sdX = _stem_update_x!(ξi, σx, llc, prc, sdX, x0_prior)
    end

    # updates within the parent branch
    llc, sdX = _update_x!(ξi, σx, llc, sdX)

    # get fixed tip 
    lξi = fixtip(ξi) 

    # make between decoupled trees node update
    llc, sdX = _update_triad_x!(lξi, ξ1, ξ2, σx, llc, sdX)
  end

  # carry on updates in the daughters
  llc, sdX = _update_x!(ξ1, σx, llc, sdX)
  llc, sdX = _update_x!(ξ2, σx, llc, sdX)

  return llc, prc, sdX
end




"""
    _stem_update_x!(ξi      ::T,
                    σx      ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    sdX     ::Float64,
                    x0_prior::NTuple{2,Float64}) where {T <: iTreeX}

Make crown update for trait.
"""
function _stem_update_x!(ξi      ::T,
                         σx      ::Float64,
                         llc     ::Float64,
                         prc     ::Float64,
                         sdX     ::Float64,
                         x0_prior::NTuple{2,Float64}) where {T <: iTreeX}

  m0, σx0 = x0_prior
  σx02    = σx0^2
  xo      = xi(ξi)
  m       = xf(ξi)
  el      = e(ξi)
  s2      = el * σx^2

  # gibbs sampling
  xn = rnorm((m * σx02 + m0 * s2) / (s2 + σx02), sqrt(s2 * σx02 / (s2 + σx02)))

  setxi!(ξi, xn)

  # update llc, prc and sdX
  llc += llrdnorm_μ(m, xn, xo, s2)
  prc += llrdnorm_x(xn, xo, m0, σx02)

  sdX += ((xn - m)^2 - (xo - m)^2)/(2.0*el)

  return llc, prc, sdX
end




"""
    _crown_update_x!(ξi      ::T,
                     ξ1      ::T,
                     ξ2      ::T,
                     σx      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     sdX     ::Float64,
                     x0_prior::NTuple{2,Float64}) where {T <: iTreeX}

Make crown update for trait.
"""
function _crown_update_x!(ξi      ::T,
                          ξ1      ::T,
                          ξ2      ::T,
                          σx      ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          sdX     ::Float64,
                          x0_prior::NTuple{2,Float64}) where {T <: iTreeX}

  m0, σx0 = x0_prior
  σx02    = σx0^2
  xo      = xi(ξi)
  x1      = xf(ξ1)
  x2      = xf(ξ2)
  e1      = e(ξ1)
  e2      = e(ξ2)

  inve   = 1.0/(e1 + e2)
  m      = (e2 * inve * x1 + e1 * inve * x2)
  sigma2 = e1 * e2 * inve * σx^2

  # gibbs sampling
  xn = rnorm((m * σx02 + m0 * sigma2) / (sigma2 + σx02),
             sqrt(sigma2 * σx02 / (sigma2 + σx02)))

  setxi!(ξi, xn)
  setxf!(ξi, xn)
  setxi!(ξ1, xn)
  setxi!(ξ2, xn)

  # update llc, prc and sdX
  llc += duoldnorm(xn, x1, x2, e1, e2, σx) - 
         duoldnorm(xo, x1, x2, e1, e2, σx) 

  prc += llrdnorm_x(xn, xo, m0, σx02)

  sdX += ((xn - x1)^2 - (xo - x1)^2)/(2.0*e1) +
         ((xn - x2)^2 - (xo - x2)^2)/(2.0*e2)

  return llc, prc, sdX
end




"""
    _update_x!(tree::T,
               σx  ::Float64,
               llc ::Float64,
               sdX ::Float64) where {T <:i TreeX}

Do gbm updates on a decoupled tree recursively.
"""
function _update_x!(tree::T,
                    σx  ::Float64,
                    llc ::Float64,
                    sdX ::Float64) where {T <: iTreeX}

  if isdefined(tree, :d1)
    llc, sdX = _update_triad_x!(tree, tree.d1, tree.d2, σx, llc, sdX)
    llc, sdX = _update_x!(tree.d1, σx, llc, sdX)
    llc, sdX = _update_x!(tree.d2, σx, llc, sdX)
  end

  return llc, sdX
end




"""
    _update_triad_x!(tree::T,
                     σx  ::Float64,
                     llc ::Float64,
                     sdX ::Float64) where {T <: iTreeX}

Make gibbs node update for trait.
"""
function _update_triad_x!(tree::T,
                          tre1::T,
                          tre2::T,
                          σx  ::Float64,
                          llc ::Float64,
                          sdX ::Float64) where {T <: iTreeX}

  xa = xi(tree)
  xo = xf(tree)
  x1 = xf(tre1)
  x2 = xf(tre2)
  ea = e(tree)
  e1 = e(tre1)
  e2 = e(tre2)

  # gibbs sampling
  xn = trioprop(xa, x1, x2, ea, e1, e2, σx)

  setxf!(tree, xn)
  setxi!(tre1, xn)
  setxi!(tre2, xn)

  # update llc, prc and sdX
  llc += trioldnorm(xn, xa, x1, x2, ea, e1, e2, σx) - 
         trioldnorm(xo, xa, x1, x2, ea, e1, e2, σx) 

  sdX += ((xa - xn)^2 - (xa - xo)^2)/(2.0*ea) +
         ((xn - x1)^2 - (xo - x1)^2)/(2.0*e1) +
         ((xn - x2)^2 - (xo - x2)^2)/(2.0*e2)

  return llc, sdX
end




"""
    update_σx!(σxc     ::Float64,
               sdX     ::Float64,
               nX      ::Float64,
               llc     ::Float64,
               prc     ::Float64,
               σx_prior::NTuple{2,Float64})

Gibbs update for `σx`.
"""
function update_σx!(σxc     ::Float64,
                    sdX     ::Float64,
                    nX      ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    σx_prior::NTuple{2,Float64})

  σx_p1 = σx_prior[1]
  σx_p2 = σx_prior[2]

  # Gibbs update for σ
  σxp2 = randinvgamma(σx_p1 + 0.5 * nX, σx_p2 + sdX)
  σxp  = sqrt(σxp2)

  # update likelihood and prior
  llc += sdX*(1.0/σxc^2 - 1.0/σxp2) - nX*(log(σxp/σxc))

  prc += llrdinvgamma(σxp2, σxc^2, σx_p1, σx_p2)

  return llc, prc, σxp
end

