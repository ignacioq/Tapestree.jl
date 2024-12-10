#=

Diffused Brownian motion MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#




"""
    insane_dbm(tree     ::Tlabel,
               xa       ::Dict{String, Float64};
               xs       ::Dict{String, Float64} = Dict{String,Float64}(),
               αx_prior ::NTuple{2,Float64} = (0.0, 10.0),
               ασ_prior ::NTuple{2,Float64} = (0.0, 10.0),
               γ_prior  ::NTuple{2,Float64} = (0.05, 0.05),
               niter    ::Int64             = 1_000,
               nthin    ::Int64             = 10,
               nburn    ::Int64             = 200,
               nflush   ::Int64             = nthin,
               ofile    ::String            = string(homedir(), "/dbm"),
               αi       ::Float64           = 0.0,
               γi       ::Float64           = 1e-3,
               pupdp    ::NTuple{4,Float64} = (0.1, 0.1, 0.05, 0.9),
               δt       ::Float64           = 1e-3,
               stn      ::Float64           = 0.1,
               mxthf    ::Float64           = Inf,
               prints   ::Int64             = 5)

Run diffused Brownian motion trait evolution model.
"""
function insane_dbm(tree     ::Tlabel,
                    xa       ::Dict{String, Float64};
                    xs       ::Dict{String, Float64} = Dict{String,Float64}(),
                    αx_prior ::NTuple{2,Float64} = (0.0, 10.0),
                    ασ_prior ::NTuple{2,Float64} = (0.0, 10.0),
                    γ_prior  ::NTuple{2,Float64} = (0.05, 0.05),
                    niter    ::Int64             = 1_000,
                    nthin    ::Int64             = 10,
                    nburn    ::Int64             = 200,
                    nflush   ::Int64             = nthin,
                    ofile    ::String            = string(homedir(), "/dbm"),
                    αxi      ::Float64           = 0.0,
                    ασi      ::Float64           = 0.0,
                    γi       ::Float64           = 1e-3,
                    pupdp    ::NTuple{5,Float64} = (0.1, 0.1, 0.1, 0.05, 0.9),
                    δt       ::Float64           = 1e-3,
                    stn      ::Float64           = 0.1,
                    mxthf    ::Float64           = Inf,
                    prints   ::Int64             = 5)

  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # set tips sampling fraction
  tl = labels(tree)
  tρ = Dict(tl[i] => 1.0 for i in 1:n)

  # set errors
  if iszero(length(xs))
    xs = Dict(tl[i] => 0.0 for i in 1:n)
  end

  # estimate branch split (multiple of δt)
  maxt = δt * floor(th * mxthf/δt)

  # make fix tree directory
  idf, xr, σxi = make_idf(tree, tρ, xa, xs, maxt)

  # make a decoupled tree
  Ξ = make_Ξ(idf, xr, log(σxi), γi, δt, srδt, sTxs)

  # get vector of internal edges
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # parameter updates (1: αx, 2: ασ, 3:γ, 4: scale 5: gbm)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running diffused Brownian motion"

  # burn-in phase
  Ξ, idf, ll, Ls, Xs, ddσ, ssσ, nσ, prc, αxc, ασc, γc, stn =
    mcmc_burn_dbm(Ξ, idf, αx_prior, ασ_prior, γ_prior, nburn, αxi, ασi, γi, stn, 
      δt, srδt, inodes, pup, prints)

  # mcmc
  r, treev = mcmc_dbm(Ξ, idf, ll, Ls, Xs, ddσ, ssσ, nσ, prc, αxc, ασc, γc, stn, 
              αx_prior, ασ_prior, γ_prior, δt, srδt, inodes, pup, 
              niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_dbm(Ξ        ::Vector{sTxs},
                  idf      ::Vector{iBffs},
                  αx_prior ::NTuple{2,Float64},
                  ασ_prior ::NTuple{2,Float64},
                  γ_prior  ::NTuple{2,Float64},
                  nburn    ::Int64,
                  αxc      ::Float64,
                  ασc      ::Float64,
                  γc       ::Float64,
                  stn      ::Float64,
                  δt       ::Float64,
                  srδt     ::Float64,
                  inodes   ::Array{Int64,1},
                  pup      ::Array{Int64,1},
                  prints   ::Int64)

MCMC burn-in chain for diffused Brownian motion.
"""
function mcmc_burn_dbm(Ξ        ::Vector{sTxs},
                       idf      ::Vector{iBffs},
                       αx_prior ::NTuple{2,Float64},
                       ασ_prior ::NTuple{2,Float64},
                       γ_prior  ::NTuple{2,Float64},
                       nburn    ::Int64,
                       αxc      ::Float64,
                       ασc      ::Float64,
                       γc       ::Float64,
                       stn      ::Float64,
                       δt       ::Float64,
                       srδt     ::Float64,
                       inodes   ::Array{Int64,1},
                       pup      ::Array{Int64,1},
                       prints   ::Int64)

  # starting likelihood and prior
  ll = zeros(lastindex(Ξ))
  llik_dbm_v!(ll, Ξ, αxc, ασc, γc, δt)
  prc = logdnorm(αxc,      αx_prior[1], αx_prior[2]^2) + 
        logdnorm(ασc,      ασ_prior[1], ασ_prior[2]^2) + 
        logdinvgamma(γc^2, γ_prior[1],  γ_prior[2])

  L   = [e(ξ) for ξ in Ξ]  # edge lengths
  nin = lastindex(inodes)  # number of internal nodes
  el  = lastindex(idf)     # number of edges

  # delta change, sum squares, path length in log-σ(t)
  Ls, Xs, ddσ, ssσ, nσ = sss_v(Ξ, xv, lσ2, ασc)

  # for scale tuning
  ltn = 0
  lup = lac = 0.0

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update trait drift `αx`
      if pupi === 1

        prc, αxc = update_αx!(αxc, Ls, Xs, ll, prc, αx_prior)

      # update rate drift `ασ`
      elseif pupi === 2

        prc, ασc = update_ασ!(ασc, γc, L, ddσ, ll, prc, ασ_prior)

        _ss!(ssσ, Ξ, lσ2, ασc)

      # update rate diffusion `γ`
      elseif pupi === 3

        prc, γc = update_γ!(γc, ssσ, nσ, ll, prc, γ_prior)

      elseif pupi === 4

        lac += update_scale!(Ξ, αxc, ll, Ls, Xs, stn, δt)
        lup += 1.0

      # update traits and rates
      else

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]
        update_xs!(bix, Ξ, idf, αxc, ασc, γc, ll, Ls, Xs, ddσ, ssσ, δt, srδt)

      end
    end

    ltn += 1
    if ltn === 100
      stn = tune(stn, lac/lup)
      ltn = 0
    end

    next!(pbar)
  end

  return Ξ, idf, ll, Ls, Xs, ddσ, ssσ, nσ, prc, αxc, ασc, γc, stn
end




"""
    mcmc_dbm(Ξ        ::Vector{sTxs},
             idf      ::Vector{iBffs},
             ll       ::Vector{Float64},
             Ls       ::Vector{Float64},
             Xs       ::Vector{Float64},
             ddσ      ::Vector{Float64},
             ssσ      ::Vector{Float64},
             nσ       ::Vector{Float64},
             prc      ::Float64,
             αxc      ::Float64,
             ασc      ::Float64,
             γc       ::Float64,
             stn      ::Float64,
             ασ_prior ::NTuple{2,Float64},
             αx_prior ::NTuple{2,Float64},
             γ_prior  ::NTuple{2,Float64},
             δt       ::Float64,
             srδt     ::Float64,
             inodes   ::Array{Int64,1},
             pup      ::Vector{Int64},
             niter    ::Int64,
             nthin    ::Int64,
             nflush   ::Int64,
             ofile    ::String,
             prints   ::Int64)

MCMC chain for diffused Brownian motion.
"""
function mcmc_dbm(Ξ        ::Vector{sTxs},
                  idf      ::Vector{iBffs},
                  ll       ::Vector{Float64},
                  Ls       ::Vector{Float64},
                  Xs       ::Vector{Float64},
                  ddσ      ::Vector{Float64},
                  ssσ      ::Vector{Float64},
                  nσ       ::Vector{Float64},
                  prc      ::Float64,
                  αxc      ::Float64,
                  ασc      ::Float64,
                  γc       ::Float64,
                  stn      ::Float64,
                  ασ_prior ::NTuple{2,Float64},
                  αx_prior ::NTuple{2,Float64},
                  γ_prior  ::NTuple{2,Float64},
                  δt       ::Float64,
                  srδt     ::Float64,
                  inodes   ::Array{Int64,1},
                  pup      ::Vector{Int64},
                  niter    ::Int64,
                  nthin    ::Int64,
                  nflush   ::Int64,
                  ofile    ::String,
                  prints   ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin = lit = sthin = zero(Int64)

  r = Array{Float64,2}(undef, nlogs, 8)

  # make Ξ vector
  treev   = sTxs[]
  L       = [e(ξ) for ξ in Ξ]  # edge lengths
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  # flush to file

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tx_root\tsigma2_root\talpha_x\talpha_sigma\tgamma\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let prc = prc, αxc = αxc, ασc = ασc, γc = γc, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, prints, "running mcmc...", 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update trait drift `αx`
            if pupi === 1

              prc, αxc = update_αx!(αxc, Ls, Xs, ll, prc, αx_prior)

              # ll0 = llik_dbm(Ξ, αxc, ασc, γc, δt)
              # if !isapprox(ll0, sum(ll))
              #    @show ll0, sum(ll), it, pupi
              #    return
              # end

            # update rate drift `ασ`
            elseif pupi === 2


              prc, ασc = update_ασ!(ασc, γc, L, ddσ, ll, prc, ασ_prior)

              _ss!(ssσ, Ξ, lσ2, ασc)

              # ll0 = llik_dbm(Ξ, αxc, ασc, γc, δt)
              # if !isapprox(ll0, sum(ll))
              #    @show ll0, sum(ll), it, pupi
              #    return
              # end

            # update rate diffusion `γ`
            elseif pupi === 3

              prc, γc = update_γ!(γc, ssσ, nσ, ll, prc, γ_prior)

              # ll0 = llik_dbm(Ξ, αxc, ασc, γc, δt)
              # if !isapprox(ll0, sum(ll))
              #    @show ll0, sum(ll), it, pupi
              #    return
              # end

            elseif pupi === 4

              lac = update_scale!(Ξ, αxc, ll, Ls, Xs, stn, δt)

              # ll0 = llik_dbm(Ξ, αxc, ασc, γc, δt)
              # if !isapprox(ll0, sum(ll))
              #    @show ll0, sum(ll), it, pupi
              #    return
              # end

            # update traits and rates
            else

              nix = ceil(Int64,rand()*nin)
              bix = inodes[nix]
              update_xs!(bix, Ξ, idf, αxc, ασc, γc, ll, Ls, Xs, ddσ, ssσ, δt, srδt)

              # ll0 = llik_dbm(Ξ, αxc, ασc, γc, δt)
              # if !isapprox(ll0, sum(ll))
              #    @show ll0, sum(ll), it, pupi
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
              r[lit,5] = exp(lσ2(Ξ[1])[1])
              r[lit,6] = αxc
              r[lit,7] = ασc
              r[lit,8] = γc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            write(of, 
              string(Float64(it), "\t", sum(ll), "\t", prc, "\t", 
                xv(Ξ[1])[1],"\t", exp(lσ2(Ξ[1])[1]), "\t", αxc, "\t", ασc, "\t",
                 γc, "\n"))
            flush(of)
            write(tf, 
              string(istring(couple(Ξ, idf, 1)), "\n"))
            flush(tf)
            sthin = zero(Int64)
          end

          next!(pbar)
        end
        
        return r, treev
      end
    end
  end
end




"""
    update_αx!(αxc     ::Float64,
               Ls      ::Vector{Float64},
               Xs      ::Vector{Float64}, 
               ll      ::Vector{Float64},
               prc     ::Float64,
               αx_prior::NTuple{2,Float64})

Gibbs update for `αx`.
"""
function update_αx!(αxc     ::Float64,
                    Ls      ::Vector{Float64},
                    Xs      ::Vector{Float64}, 
                    ll      ::Vector{Float64},
                    prc     ::Float64,
                    αx_prior::NTuple{2,Float64})

  # ratio
  ν   = αx_prior[1]
  τ2  = αx_prior[2]^2
  idn = 1.0/(sum(Ls)*τ2 + 1.0)

  # gibbs update for σ
  αxp = rnorm((sum(Xs)*τ2 + ν)*idn, sqrt(τ2*idn))

  # update prior
  prc += llrdnorm_x(αxp, αxc, ν, τ2)

  # update likelihoods
  for i in Base.OneTo(lastindex(ll))
    iszero(Ls[i]) && continue
    ll[i] += 0.5*Ls[i]*(αxc^2 - αxp^2 + 2.0*Xs[i]*(αxp - αxc)/Ls[i])
  end

  return prc, αxp
end




"""
    update_ασ!(ασc     ::Float64,
               γ       ::Float64,
               L       ::Vector{Float64},
               ddσ     ::Vector{Float64},
               ll      ::Vector{Float64},
               prc     ::Float64,
               ασ_prior::NTuple{2,Float64})

Gibbs update for `ασ`.
"""
function update_ασ!(ασc     ::Float64,
                    γ       ::Float64,
                    L       ::Vector{Float64},
                    ddσ     ::Vector{Float64},
                    ll      ::Vector{Float64},
                    prc     ::Float64,
                    ασ_prior::NTuple{2,Float64})

  # ratio
  ν   = ασ_prior[1]
  τ2  = ασ_prior[2]^2
  γ2  = γ^2
  rs  = γ2/τ2
  idn = 1.0/(rs + sum(L))

  # gibbs update for αp
  ασp = rnorm((sum(ddσ) + rs*ν)*idn, sqrt(γ2*idn))

  # update prior
  prc += llrdnorm_x(ασp, ασc, ν, τ2)

  # update likelihoods
  for i in Base.OneTo(lastindex(ll))
    iszero(L[i]) && continue
    ll[i] += 0.5*L[i]/γ2*(ασc^2 - ασp^2 + 2.0*ddσ[i]*(ασp - ασc)/L[i])
  end

  return prc, ασp
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
    update_scale!(Ξ   ::Vector{sTxs},
                  αx  ::Float64,
                  ll  ::Vector{Float64},
                  Ls  ::Vector{Float64},
                  Xs  ::Vector{Float64},
                  stn ::Float64,
                  δt  ::Float64)

Update scale for `σ`.
"""
function update_scale!(Ξ   ::Vector{sTxs},
                       αx  ::Float64,
                       ll  ::Vector{Float64},
                       Ls  ::Vector{Float64},
                       Xs  ::Vector{Float64},
                       stn ::Float64,
                       δt  ::Float64)

  # sample log(scaling factor)
  s = randn()*stn

  # likelihood ratio
  llr, Lsr, Xsr = llr_scale(Ξ, αx, s, δt)

  acc = 0.0
  if -randexp() < sum(llr)
    acc += 1.0
    scale_rate!(Ξ, lσ2, s)
    @turbo for i in Base.OneTo(lastindex(ll))
      ll[i] += llr[i]
      Ls[i] += Lsr[i]
      Xs[i] += Xsr[i]
    end
  end

  return acc
end




"""
    update_xs!(bix ::Int64,
               Ξ   ::Vector{sTxs},
               idf ::Vector{iBffs},
               αx  ::Float64,
               ασ  ::Float64,
               γ   ::Float64,
               ll  ::Vector{Float64},
               Ls  ::Vector{Float64},
               Xs  ::Vector{Float64},
               ddσ ::Vector{Float64},
               ssσ ::Vector{Float64},
               δt  ::Float64,
               srδt::Float64)

Make a `dbm` update for an internal branch and its descendants.
"""
function update_xs!(bix ::Int64,
                    Ξ   ::Vector{sTxs},
                    idf ::Vector{iBffs},
                    αx  ::Float64,
                    ασ  ::Float64,
                    γ   ::Float64,
                    ll  ::Vector{Float64},
                    Ls  ::Vector{Float64},
                    Xs  ::Vector{Float64},
                    ddσ ::Vector{Float64},
                    ssσ ::Vector{Float64},
                    δt  ::Float64,
                    srδt::Float64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  b1   = idf[i1]
  i2   = d2(bi)
  ξ1   = Ξ[i1]
  root = iszero(pa(bi))

  # if mrca
  if root && iszero(e(ξi))
    #if stem fossil
    if isfossil(bi)
      ll[i1], Ls[i1], Xs[i1], ddσ[i1], ssσ[i1] = 
        _fstem_update!(ξi, ξ1, αx, ασ, γ, δt, srδt)
    # if crown
    else
      ll[i1], ll[i2], Ls[i1], Ls[i2], Xs[i1], Xs[i2], ddσ[i1], ddσ[i2], 
      ssσ[i1], ssσ[i2] = _crown_update!(ξi, ξ1, Ξ[i2], αx, ασ, γ, δt, srδt)
    end
  # if stem
  elseif root
    ll[bix],  Ls[bix], Xs[bix], ddσ[bix], ssσ[bix] = 
      _stem_update!(ξi, αx, ασ, γ, δt, srδt)
  # if duo
  elseif iszero(i2)
    if ifx(bi)
      ll[bix], ll[i1], Ls[bix], Ls[i1], Xs[bix], Xs[i1], 
      ddσ[bix], ddσ[i1], ssσ[bix], ssσ[i1] = 
        _update_duo_x!(ξi, ξ1, xavg(bi), xstd(bi), αx, ασ, γ, δt, srδt)
    else
      ll[bix], ll[i1], Ls[bix], Ls[i1], Xs[bix], Xs[i1], 
      ddσ[bix], ddσ[i1], ssσ[bix], ssσ[i1] = 
        _update_duo_x!(ξi, ξ1, αx, ασ, γ, δt, srδt)
    end
  # if triad
  else
    ξ2 = Ξ[i2]
    ll[bix],  ll[i1],   ll[i2], 
    Ls[bix],  Ls[i1],   Ls[i2], 
    Xs[bix],  Xs[i1],   Xs[i2], 
    ddσ[bix], ddσ[i1], ddσ[i2], 
    ssσ[bix], ssσ[i1], ssσ[i2] =
      _update_triad_x!(ξi, ξ1, ξ2, αx, ασ, γ, δt, srδt)
  end

  # update daughters
  if iszero(d1(b1))
    if ifx(b1) 
      ll[i1], Ls[i1], Xs[i1], ddσ[i1], ssσ[i1] = 
        _update_leaf_x!(ξ1, xavg(b1), xstd(b1), αx, ασ, γ, δt, srδt)
    else
      ll[i1], Ls[i1], Xs[i1], ddσ[i1], ssσ[i1] = 
        _update_leaf_x!(ξ1, αx, ασ, γ, δt, srδt)
    end
  end

  if i2 > 0
    b2 = idf[i2]
    if iszero(d1(b2))
      ξ2 = Ξ[i2]
      if ifx(b2)
        ll[i2], Ls[i2], Xs[i2], ddσ[i2], ssσ[i2] = 
          _update_leaf_x!(ξ2, xavg(b2), xstd(b2), αx, ασ, γ, δt, srδt)
      else
        ll[i2], Ls[i2], Xs[i2], ddσ[i2], ssσ[i2] = 
          _update_leaf_x!(ξ2, αx, ασ, γ, δt, srδt)
      end
    end
  end

  return nothing
end



