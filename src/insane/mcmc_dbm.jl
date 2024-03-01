#=

Diffused Brownian motion MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#




"""
    insane_dbm(tree   ::Tlabel,
               xa     ::Dict{String, Float64};
               xs     ::Dict{String, Float64} = Dict("" => 0.0),
               γ_prior::NTuple{2,Float64}     = (3.0, 0.5),
               niter  ::Int64                 = 1_000,
               nthin  ::Int64                 = 10,
               nburn  ::Int64                 = 200,
               nflush ::Int64                 = nthin,
               ofile  ::String                = string(homedir(), "/ipb"),
               γi     ::Float64                = 0.1,
               pupdp  ::NTuple{2,Float64}     = (0.1, 0.9),
               δt     ::Float64               = 1e-3,
               prints ::Int64                 = 5)

Run diffused Brownian motion trait evolution model.
"""
function insane_dbm(tree   ::Tlabel,
                    xa     ::Dict{String, Float64};
                    xs     ::Dict{String, Float64} = Dict{String,Float64}(),
                    γ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                    niter  ::Int64                 = 1_000,
                    nthin  ::Int64                 = 10,
                    nburn  ::Int64                 = 200,
                    nflush ::Int64                 = nthin,
                    ofile  ::String                = string(homedir(), "/dbm"),
                    γi     ::Float64               = 0.1,
                    pupdp  ::NTuple{3,Float64}     = (0.1, 0.05, 0.9),
                    δt     ::Float64               = 1e-3,
                    stn    ::Float64               = 0.1,
                    prints ::Int64                 = 5)

  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # set tips sampling fraction
  tl = labels(tree)
  tρ = Dict(tl[i] => 1.0 for i in 1:n)

  if iszero(length(xs))
    xs = Dict(tl[i] => 0.0 for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxi = make_idf(tree, tρ, xa, xs, Inf)

  # make a decoupled tree
  Ξ = make_Ξ(idf, xr, log(σxi), γi, δt, srδt, sTxs)

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # parameter updates (1: γ, 2: scale 3: gbm)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running diffused Brownian motion"

  # burn-in phase
  Ξ, idf, ll, ssσ, nσ, prc, γc, stn =
    mcmc_burn_dbm(Ξ, idf, γ_prior, nburn, γi, stn, δt, srδt, 
      inodes, pup, prints)

  # mcmc
  r, treev = mcmc_dbm(Ξ, idf, ll, ssσ, nσ, prc, γc, stn, γ_prior, δt, srδt, 
               inodes, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_dbm(Ξ       ::Vector{sTxs},
                  idf     ::Vector{iBffs},
                  γ_prior::NTuple{2,Float64},
                  nburn   ::Int64,
                  γc      ::Float64,
                  δt      ::Float64,
                  srδt    ::Float64,
                  inodes  ::Array{Int64,1},
                  pup     ::Array{Int64,1},
                  prints  ::Int64)

MCMC burn-in chain for diffused Brownian motion.
"""
function mcmc_burn_dbm(Ξ       ::Vector{sTxs},
                       idf     ::Vector{iBffs},
                       γ_prior::NTuple{2,Float64},
                       nburn   ::Int64,
                       γc      ::Float64,
                       stn     ::Float64,
                       δt      ::Float64,
                       srδt    ::Float64,
                       inodes  ::Array{Int64,1},
                       pup     ::Array{Int64,1},
                       prints  ::Int64)

  # starting likelihood and prior
  ll  = zeros(lastindex(Ξ)) 
  llik_dbm_v!(ll, Ξ, γc, δt)
  prc = logdinvgamma(γc^2, γ_prior[1], γ_prior[2])

  # sum squares in log-σ(t)
  ssσ, nσ = sss_v(Ξ, lσ)
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  # for scale tuning
  ltn = 0
  lup = lac = 0.0

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update rate drift `γ`
      if pupi === 1

        prc, γc = update_γ!(γc, ssσ, nσ, ll, prc, γ_prior)

      elseif pupi === 2

        lac += update_scale!(Ξ, ll, stn, δt)
        lup += 1.0

      # update traits and rates
      else

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]
        update_x!(bix, Ξ, idf, γc, ll, ssσ, δt, srδt)

      end
    end

    ltn += 1
    if ltn === 100
      stn = tune(stn, lac/lup)
      ltn = 0
    end

    next!(pbar)
  end

  return Ξ, idf, ll, ssσ, nσ, prc, γc, stn
end




"""
    mcmc_dbm(Ξ       ::Vector{sTxs},
             idf     ::Vector{iBffs},
             ll      ::Vector{Float64},
             ssσ     ::Vector{Float64},
             nσ      ::Vector{Float64},
             prc     ::Float64,
             γc     ::Float64,
             γ_prior ::NTuple{2,Float64},
             δt      ::Float64,
             srδt    ::Float64,
             inodes  ::Array{Int64,1},
             pup     ::Vector{Int64},
             niter   ::Int64,
             nthin   ::Int64,
             nflush  ::Int64,
             ofile   ::String,
             prints  ::Int64)

MCMC chain for diffused Brownian motion.
"""
function mcmc_dbm(Ξ       ::Vector{sTxs},
                  idf     ::Vector{iBffs},
                  ll      ::Vector{Float64},
                  ssσ     ::Vector{Float64},
                  nσ      ::Vector{Float64},
                  prc     ::Float64,
                  γc      ::Float64,
                  stn     ::Float64,
                  γ_prior ::NTuple{2,Float64},
                  δt      ::Float64,
                  srδt    ::Float64,
                  inodes  ::Array{Int64,1},
                  pup     ::Vector{Int64},
                  niter   ::Int64,
                  nthin   ::Int64,
                  nflush  ::Int64,
                  ofile   ::String,
                  prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  r = Array{Float64,2}(undef, nlogs, 6)

  # make Ξ vector
  treev   = sTxs[]
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  # flush to file
  sthin = 0

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tx_root\tsigma_root\tgamma\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      pbar = Progress(niter, prints, "running mcmc...", 20)

      for it in Base.OneTo(niter)

        shuffle!(pup)

        for pupi in pup

          ## parameter updates
          # update rate drift `γ`
          if pupi === 1

            prc, γc = update_γ!(γc, ssσ, nσ, ll, prc, γ_prior)

            # ll0 = llik_dbm(Ξ, γc, δt)
            # if !isapprox(ll0, sum(ll), atol = 1e-4)
            #    @show ll0, llc, it, pupi
            #    return
            # end

          elseif pupi === 2

            lac = update_scale!(Ξ, ll, stn, δt)

            # ll0 = llik_dbm(Ξ, γc, δt)
            # if !isapprox(ll0, sum(ll), atol = 1e-4)
            #    @show ll0, llc, it, pupi
            #    return
            # end

          # update traits and rates
          else

            nix = ceil(Int64,rand()*nin)
            bix = inodes[nix]
            update_x!(bix, Ξ, idf, γc, ll, ssσ, δt, srδt)

            # ll0 = llik_dbm(Ξ, γc, δt)
            # if !isapprox(ll0, sum(ll), atol = 1e-4)
            #    @show ll0, llc, it, pupi
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
            r[lit,2] = sum(ll)
            r[lit,3] = prc
            r[lit,4] = xv(Ξ[1])[1]
            r[lit,5] = exp(lσ(Ξ[1])[1])
            r[lit,6] = γc
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthin += 1
        if sthin === nflush
          write(of, 
            string(Float64(it), "\t", sum(ll), "\t", prc, "\t", 
              xv(Ξ[1])[1],"\t", exp(lσ(Ξ[1])[1]), "\t", γc, "\n"))
          flush(of)
          write(tf, 
            string(istring(couple(Ξ, idf, 1)), "\n"))
          flush(tf)
          sthin = 0
        end

        next!(pbar)
      end
    end
  end

  return r, treev
end






"""
    update_γ!(γc     ::Float64,
              ssσ    ::Vector{Float64},
              ns     ::Vector{Float64},
              ll     ::Vector{Float64},
              prc    ::Float64,
              γ_prior::NTuple{2,Float64})

Gibbs update for `γ`.
"""
function update_γ!(γc     ::Float64,
                   ssσ    ::Vector{Float64},
                   ns     ::Vector{Float64},
                   ll     ::Vector{Float64},
                   prc    ::Float64,
                   γ_prior::NTuple{2,Float64})

  γ_p1, γ_p2 = γ_prior

  # Gibbs update for σ
  γp2 = randinvgamma(γ_p1 + 0.5 * sum(ns), γ_p2 + sum(ssσ))

  # update prior
  prc += llrdinvgamma(γp2, γc^2, γ_p1, γ_p2)

  γp = sqrt(γp2)

  # update likelihoods
  for i in Base.OneTo(lastindex(ll))
    ll[i] += ssσ[i]*(1.0/γc^2 - 1.0/γp2) - ns[i]*(log(γp/γc))
  end

  return prc, γp
end




"""
    update_scale!(Ξ  ::Vector{T},
                  idf::Vector{iBffs},
                  llc::Float64,
                  ir ::Float64,
                  ns ::Float64,
                  stn::Float64) where {T <: iTree}

Update scale for speciation.
"""
function update_scale!(Ξ   ::Vector{sTxs},
                       ll  ::Vector{Float64},
                       stn ::Float64,
                       δt  ::Float64)

  # sample log(scaling factor)
  s = randn()*stn

  # likelihood ratio
  llr = llr_scale(Ξ, s, δt)

  acc = 0.0

  if -randexp() < sum(llr)
    acc += 1.0
    scale_rate!(Ξ, lσ, s)
    @turbo for i in Base.OneTo(lastindex(ll))
      ll[i] += llr[i]
    end
  end

  return acc
end




"""
    update_x!(bix  ::Int64,
              Ξ    ::Vector{sTxs},
              idf  ::Vector{iBffs},
              γ    ::Float64,
              ll   ::Vector{Float64},
              ssσ  ::Vector{Float64},
              δt   ::Float64,
              srδt ::Float64)

Make a `dbm` update for an internal branch and its descendants.
"""
function update_x!(bix  ::Int64,
                   Ξ    ::Vector{sTxs},
                   idf  ::Vector{iBffs},
                   γ    ::Float64,
                   ll   ::Vector{Float64},
                   ssσ  ::Vector{Float64},
                   δt   ::Float64,
                   srδt ::Float64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  b1   = idf[i1]
  i2   = d2(bi)
  ξ1   = Ξ[i1]
  root = iszero(pa(bi))

  # if crown root
  if root && iszero(e(ξi))
    ll[i1], ll[i2], ssσ[i1], ssσ[i2] =
      _crown_update!(ξi, ξ1, Ξ[i2], γ, δt, srδt)
  # if stem
  elseif root
    ll[bix], ssσ[bix] = _stem_update!(ξi, γ, δt, srδt)
  # if duo
  elseif iszero(i2)
    if ifx(bi)
      ll[bix], ll[i1], ssσ[bix], ssσ[i1] = 
        _update_duo_x!(ξi, ξ1, xavg(bi), xstd(bi), γ, δt, srδt)
    else
      ll[bix], ll[i1], ssσ[bix], ssσ[i1] = _update_duo_x!(ξi, ξ1, γ, δt, srδt)
    end
  # if triad
  else
    ξ2 = Ξ[i2]
    ll[bix], ll[i1], ll[i2], ssσ[bix], ssσ[i1], ssσ[i2] =
      _update_triad_x!(ξi, ξ1, ξ2, γ, δt, srδt)
  end

  # update daughters
  if iszero(d1(b1))
    if ifx(b1) 
      ll[i1], ssσ[i1] = _update_leaf_x!(ξ1, xavg(b1), xstd(b1), γ, δt, srδt)
    else
      ll[i1], ssσ[i1] = _update_leaf_x!(ξ1, γ, δt, srδt)
    end
  end

  if i2 > 0
    b2 = idf[i2]
    if iszero(d1(b2))
      ξ2 = Ξ[i2]
      if ifx(b2)
        ll[i2], ssσ[i2] = _update_leaf_x!(ξ2, xavg(b2), xstd(b2), γ, δt, srδt)
      else
        ll[i2], ssσ[i2] = _update_leaf_x!(ξ2, γ, δt, srδt)
      end
    end
  end

  return nothing
end



