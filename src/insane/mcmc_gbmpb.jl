#=

GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#




"""
    insane_gbmpb(tree    ::sT_label;
                 α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64}     = (3.0, 0.5),
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 nflush  ::Int64                 = nthin,
                 ofile   ::String                = string(homedir(), "/ipb"),
                 αi      ::Float64               = 0.0,
                 σλi     ::Float64               = 0.1,
                 pupdp   ::NTuple{5,Float64}     = (0.01, 0.01, 0.01, 0.1, 0.2),
                 δt      ::Float64               = 1e-3,
                 prints  ::Int64                 = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `pbd`.
"""
function insane_gbmpb(tree    ::sT_label;
                      α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64}     = (3.0, 0.5),
                      niter   ::Int64                 = 1_000,
                      nthin   ::Int64                 = 10,
                      nburn   ::Int64                 = 200,
                      nflush  ::Int64                 = nthin,
                      ofile   ::String                = string(homedir(), "/ipb"),
                      αi      ::Float64               = 0.0,
                      σλi     ::Float64               = 0.1,
                      pupdp   ::NTuple{5,Float64}     = (0.01, 0.01, 0.01, 0.1, 0.2),
                      δt      ::Float64               = 1e-3,
                      prints  ::Int64                 = 5,
                      stn     ::Float64               = 0.5,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ, Inf)

  # make a decoupled tree
  Ξ = make_Ξ(idf, λmle_cpb(tree), αi, σλi, δt, srδt, iTpb)

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf))  if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σ, 3: scale, 4: gbm, 5: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running pure-birth gbm"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, ns, stn =
    mcmc_burn_gbmpb(Ξ, idf, α_prior, σλ_prior, nburn, αi, σλi, stn,
      δt, srδt, inodes, pup, prints)

  # mcmc
  r, treev = mcmc_gbmpb(Ξ, idf, llc, prc, αc, σλc, ns, stn, α_prior, σλ_prior,
              δt, srδt, inodes, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end



"""
    mcmc_burn_gbmpb(Ξ       ::Vector{iTpb},
                    idf     ::Vector{iBffs},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    stn     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

MCMC burn-in chain for `pbd`.
"""
function mcmc_burn_gbmpb(Ξ       ::Vector{iTpb},
                         idf     ::Vector{iBffs},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         stn     ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  nsi = Float64(iszero(e(Ξ[1])))

  # starting likelihood and prior
  llc = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) - nsi*lλ(Ξ[1])[1] + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2]) +
        logdnorm(αc,         α_prior[1],  α_prior[2]^2)

  L   = treelength(Ξ)      # tree length
  nin = lastindex(inodes)  # number of internal nodes
  el  = lastindex(idf)     # number of branches
  ns  = Float64(nin) - nsi # number of speciation events in likelihood

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, nλ, irλ = 
    _ss_ir_dd(Ξ, lλ, αc)

  # for scale tuning
  ltn = 0
  lup = lac = 0.0

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update drift
      if pupi === 1

        llc, prc, αc = update_α!(αc, σλc, L, ddλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ = _ss(Ξ, lλ, αc)

      # update diffusion
      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, ssλ, nλ, llc, prc, σλ_prior)

      # update scale
      elseif pupi === 3

        llc, irλ, acc = update_scale!(Ξ, idf, llc, irλ, ns, stn)

        lac += acc
        lup += 1.0
      # update gbm
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, ddλ, ssλ, irλ =
          update_gbm!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, irλ, δt, srδt)

      # forward simulation
      else

        bix = ceil(Int64,rand()*el)

        llc, ddλ, ssλ, nλ, irλ, ns, L =
          update_fs!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, nλ, irλ, ns, L, 
            δt, srδt)

      end
    end

    ltn += 1
    if ltn === 100
      stn = tune(stn, lac/lup)
      ltn = 0
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, ns, stn
end




"""
    mcmc_gbmpb(Ξ       ::Vector{iTpb},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               ns      ::Float64,
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Vector{Int64},
               niter   ::Int64,
               nthin   ::Int64,
               nflush  ::Int64,
               ofile   ::String,
               prints  ::Int64)

MCMC chain for pure-birth diffusion.
"""
function mcmc_gbmpb(Ξ       ::Vector{iTpb},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    ns      ::Float64,
                    stn     ::Float64,
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
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
  treev = iTpb[]

  L   = treelength(Ξ)      # tree length
  nin = lastindex(inodes)  # number of internal nodes
  el  = lastindex(idf)     # number of branches

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, nλ, irλ = 
    _ss_ir_dd(Ξ, lλ, αc)

  # flush to file
  sthin = 0

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\talpha\tsigma_lambda\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      pbar = Progress(niter, prints, "running mcmc...", 20)

      for it in Base.OneTo(niter)

        shuffle!(pup)

        for pupi in pup

          ## parameter updates
          # update drift
          if pupi === 1

            llc, prc, αc = update_α!(αc, σλc, L, ddλ, llc, prc, α_prior)

            # update ssλ with new drift `α`
            ssλ = _ss(Ξ, lλ, αc)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
            # if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi
            #    return
            # end

          # update diffusion rate
          elseif pupi === 2

            llc, prc, σλc = update_σ!(σλc, ssλ, nλ, llc, prc, σλ_prior)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
            # if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi
            #    return
            # end

          # update scale
          elseif pupi === 3

            llc, irλ, acc = update_scale!(Ξ, idf, llc, irλ, ns, stn)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
            # if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi
            #    return
            # end

          # update gbm
          elseif pupi === 4

            nix = ceil(Int64,rand()*nin)
            bix = inodes[nix]

            llc, ddλ, ssλ, irλ =
              update_gbm!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, irλ, δt, srδt)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
            # if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi
            #    return
            # end

          # update by forward simulation
          else

            bix = ceil(Int64,rand()*el)

            llc, ddλ, ssλ, nλ, irλ, ns, L =
              update_fs!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, nλ, irλ, ns, L, 
                δt, srδt)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
            # if !isapprox(ll0, llc, atol = 1e-4)
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
            r[lit,2] = llc
            r[lit,3] = prc
            r[lit,4] = exp(lλ(Ξ[1])[1])
            r[lit,5] = αc
            r[lit,6] = σλc
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthin += 1
        if sthin === nflush
          write(of, 
            string(Float64(it), "\t", llc, "\t", prc, "\t", 
              exp(lλ(Ξ[1])[1]),"\t", αc, "\t", σλc, "\n"))
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
    update_α!(αc     ::Float64,
              σλ     ::Float64,
              L      ::Float64,
              ddλ    ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              α_prior::NTuple{2,Float64})

Gibbs update for `α`.
"""
function update_α!(αc     ::Float64,
                   σλ     ::Float64,
                   L      ::Float64,
                   ddλ    ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   α_prior::NTuple{2,Float64})

  # ratio
  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2

  # gibbs update for σ
  αp = rnorm((ddλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  # update prior
  prc += llrdnorm_x(αp, αc, ν, τ2)

  # update likelihood
  llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*ddλ*(αp - αc)/L)

  return llc, prc, αp
end




"""
    update_σ!(σλc     ::Float64,
              ssλ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64})

Gibbs update for `σλ`.
"""
function update_σ!(σλc     ::Float64,
                   ssλ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   σλ_prior::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)

  # update prior
  prc += llrdinvgamma(σλp2, σλc^2, σλ_p1, σλ_p2)

  σλp = sqrt(σλp2)

  # update likelihood
  llc += ssλ*(1.0/σλc^2 - 1.0/σλp2) - n*(log(σλp/σλc))

  return llc, prc, σλp
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
function update_scale!(Ξ  ::Vector{T},
                       idf::Vector{iBffs},
                       llc::Float64,
                       ir ::Float64,
                       ns ::Float64,
                       stn::Float64) where {T <: iTree}

  # sample log(scaling factor)
  s = randn()*stn

  # likelihood ratio
  iri = (1.0 - exp(s)) * ir
  llr = ns * s + iri

  acc = 0.0

  if -randexp() < llr
    acc += 1.0
    llc += llr
    ir  -= iri
    scale_rate!(Ξ, lλ, s)
    scale_rate!(idf, s)
  end

  return llc, ir, acc
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTpb},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                llc  ::Float64,
                ddλ  ::Float64,
                ssλ  ::Float64,
                irλ  ::Float64,
                δt   ::Float64,
                srδt ::Float64)

Make a `gbm` update for an interna branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTpb},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     llc  ::Float64,
                     ddλ  ::Float64,
                     ssλ  ::Float64,
                     irλ  ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  i2   = d2(bi)
  ξ1   = Ξ[i1]
  ξ2   = Ξ[i2]
  root = iszero(pa(bi))

  # if crown root
  if root && iszero(e(ξi))
    llc, ddλ, ssλ, irλ =
      _crown_update!(ξi, ξ1, ξ2, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt)
    setλt!(bi, lλ(ξi)[1])
  else
    # if stem
    if root
      llc, ddλ, ssλ, irλ = 
        _stem_update!(ξi, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt)
    end

    # updates within the parent branch
    llc, ddλ, ssλ, irλ = 
      _update_gbm!(ξi, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt, false)

    # get fixed tip
    lξi = fixtip(ξi)

    # make between decoupled trees node update
    llc, ddλ, ssλ, irλ = 
      update_triad_pb!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2),
        fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, llc, ddλ, ssλ, irλ, δt, srδt)

    # set fixed `λ(t)` in branch
    setλt!(bi, lλ(lξi)[end])
  end

  # # carry on updates in the daughters
  llc, ddλ, ssλ, irλ = 
    _update_gbm!(ξ1, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt, iszero(d1(idf[i1])))
  llc, ddλ, ssλ, irλ = 
    _update_gbm!(ξ2, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt, iszero(d1(idf[i2])))

  return llc, ddλ, ssλ, irλ
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{iTpb},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               llc  ::Float64,
               ddλ  ::Float64,
               ssλ  ::Float64,
               nλ   ::Float64,
               irλ  ::Float64,
               ns   ::Float64,
               L    ::Float64,
               δt   ::Float64,
               srδt ::Float64)

Forward simulation proposal function for pure birth diffusion.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{iTpb},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    llc  ::Float64,
                    ddλ  ::Float64,
                    ssλ  ::Float64,
                    nλ   ::Float64,
                    irλ  ::Float64,
                    ns   ::Float64,
                    L    ::Float64,
                    δt   ::Float64,
                    srδt ::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]

  # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, ξc, α, σλ, δt, srδt)
    ddrλ = ssrλ = irrλ = 0.0
  # if internal
  else
    ξp, llr, ddrλ, ssrλ, irrλ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, δt, srδt)
  end

  # if accepted
  if isfinite(llr)
    ll1, ddλ1, ssλ1, nλ1, irλ1, ns1 = llik_gbm_ssλ(ξp, α, σλ, δt, srδt, 0.0)
    ll0, ddλ0, ssλ0, nλ0, irλ0, ns0 = llik_gbm_ssλ(ξc, α, σλ, δt, srδt, 0.0)

    # update llr, ssλ, nλ, L
    llc += ll1  - ll0 + llr
    ddλ += ddλ1 - ddλ0 + ddrλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    nλ  += nλ1  - nλ0
    irλ += irλ1 - irλ0 + irrλ
    ns  += ns1  - ns0
    L   += treelength(ξp) - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ssλ, nλ, irλ, ns, L
end




"""
    fsbi_t(bi  ::iBffs,
           ξc  ::iTpb,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTpb,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, nap, nn, llr =
    _sim_gbmpb_t(e(bi), lλ(ξc)[1], α, σλ, δt, srδt, lc, lU, Iρi, 0, 1, 500)

  if isfinite(llr)
    _fixrtip!(t0, nap) # fix random tip
    setni!(bi, nap)    # set new ni

    return t0, llr
  else
    return t0, NaN
  end
end




"""
    fsbi_i(bi  ::iBffs,
           ξ1  ::iTpb,
           ξ2  ::iTpb,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTpb,
                ξ1  ::iTpb,
                ξ2  ::iTpb,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  t0, na = _sim_gbmpb(e(bi), lλ(ξc)[1], α, σλ, δt, srδt, 1, 1_000)

  if na >= 1_000
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

 # fix random tip
  λf = fixrtip!(t0, na, NaN)

  llrd, acrd, drλ, ssrλ, irrλ, λ1p, λ2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulated remaining tips until the present
    t0, na, acr =
      tip_sims!(t0, tf(bi), α, σλ, δt, srδt, acr, lU, Iρi, na)

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      l1  = lastindex(λ1p)
      l2  = lastindex(λ2p)
      setnt!(bi, ntp)                    # set new nt
      setni!(bi, na)                     # set new ni
      setλt!(bi, λf)                     # set new λt
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 2 λ vector

      return t0, llr, drλ, ssrλ, irrλ
    else
      return t0, NaN, NaN, NaN, NaN
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTpb,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTpb,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64)

 if lU < lr && na < 1_000

    if istip(tree)
      if !isfix(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)

        # simulate
        stree, na, lr =
          _sim_gbmpb_it(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, δt, srδt,
            lr, lU, Iρi, na, 1_000)

        if isnan(lr) || na >= 1_000
          return tree, na, NaN
        end

        sete!(tree, e(tree) + e(stree))

        lλs = lλ(stree)

        if lastindex(lλs) === 2
          setfdt!(tree, fdt(tree) + fdt(stree))
        else
          setfdt!(tree, fdt(stree))
        end

        pop!(lλ0)
        popfirst!(lλs)
        append!(lλ0, lλs)

        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, lr = tip_sims!(tree.d1, t, α, σλ, δt, srδt, lr, lU, Iρi, na)
      tree.d2, na, lr = tip_sims!(tree.d2, t, α, σλ, δt, srδt, lr, lU, Iρi, na)
    end

    return tree, na, lr
  end

  return tree, na, NaN
end




"""
    tune(window::Float64, acc_rate::Float64)

Tune proposal based on acceptance rate.
"""
function tune(window::Float64, acc_rate::Float64)

  if acc_rate > 0.234
    window *= (1.0 + (acc_rate - 0.234) * 1.3054830287206267)
  else
    window /= (2.0 - acc_rate * 4.273504273504273)
  end

  return window
end




# """
#     update_α!(αc     ::Float64,
#               σλ     ::Float64,
#               L      ::Float64,
#               dλ     ::Float64,
#               llc    ::Float64,
#               prc    ::Float64,
#               rdc    ::Float64,
#               α_prior::NTuple{2,Float64},
#               α_rdist ::NTuple{2,Float64},
#               pow    ::Float64)


# Gibbs update for `α` given reference distribution.
# """
# function update_α!(αc     ::Float64,
#                    σλ     ::Float64,
#                    L      ::Float64,
#                    dλ     ::Float64,
#                    llc    ::Float64,
#                    prc    ::Float64,
#                    rdc    ::Float64,
#                    α_prior::NTuple{2,Float64},
#                    α_rdist ::NTuple{2,Float64},
#                    pow    ::Float64)

#   # ratio
#   ν   = α_prior[1]
#   τ2  = α_prior[2]^2
#   σλ2 = σλ^2
#   rs  = σλ2/τ2

#   cpow = (1.0 - pow)

#   # gibbs update for α
#   m   = (dλ + rs*ν)/(rs + L)
#   s2  = σλ2/(rs + L)
#   m0  = α_rdist[1]
#   s02 = α_rdist[2]^2
#   αp  = rnorm((m0 * s2 * cpow + m * s02 * pow) / (pow * s02 + s2 * cpow),
#               sqrt( s2 * s02 /  (pow * s02 + s2 * cpow)) )

#   # update likelihood, prior and reference
#   llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/L)
#   prc += llrdnorm_x(αp, αc, ν, τ2)
#   rdc += llrdnorm_x(αp, αc, m0, s02)

#   return llc, prc, rdc, αp
# end





# """
#     update_σ!(σc     ::Float64,
#               ss     ::Float64,
#               n       ::Float64,
#               llc     ::Float64,
#               prc     ::Float64,
#               σλ_prior::NTuple{2,Float64})

# Gibbs update for `σλ` given reference distribution.
# """
# function update_σ!(σc     ::Float64,
#                    ss     ::Float64,
#                    n       ::Float64,
#                    llc     ::Float64,
#                    prc     ::Float64,
#                    rdc     ::Float64,
#                    σλ_prior::NTuple{2,Float64},
#                    σλ_rdist ::NTuple{2,Float64},
#                    pow     ::Float64)

#   σλ_p1 = σλ_prior[1]
#   σλ_p2 = σλ_prior[2]

#   # Gibbs update for σ
#   σλp2 = randinvgamma((σλ_p1 + 0.5 * n) * pow + σλ_rdist[1] * (1.0 - pow),
#                       (σλ_p2 + ss) * pow     + σλ_rdist[2] * (1.0 - pow))

#   # update likelihood, prior and reference
#   σλp = sqrt(σλp2)
#   llc += ss*(1.0/σc^2 - 1.0/σλp2) - n*(log(σλp/σc))
#   prc += llrdinvgamma(σλp2, σλc^2, σλ_p1, σλ_p2)
#   rdc += llrdinvgamma(σλp2, σλc^2, σλ_rdist[1], σλ_rdist[2])

#   return llc, prc, rdc, σλp
# end






# """
#     ref_posterior(Ξ       ::Vector{iTpb},
#                   idf     ::Vector{iBffs},
#                   llc     ::Float64,
#                   prc     ::Float64,
#                   αc      ::Float64,
#                   σλc     ::Float64,
#                   λ0_prior::NTuple{2,Float64},
#                   α_prior ::NTuple{2,Float64},
#                   σλ_prior::NTuple{2,Float64},
#                   λ0rdist::NTuple{2,Float64},
#                   α_rdist ::NTuple{2,Float64},
#                   σλ_rdist::NTuple{2,Float64},
#                   nitpp   ::Int64,
#                   nthpp   ::Int64,
#                   βs      ::Vector{Float64},
#                   δt      ::Float64,
#                   srδt    ::Float64,
#                   inodes  ::Array{Int64,1},
#                   pup     ::Array{Int64,1},
#                   prints  ::Int64)

# MCMC chain for GBM pure-birth.
# """
# function ref_posterior(Ξ       ::Vector{iTpb},
#                        idf     ::Vector{iBffs},
#                        llc     ::Float64,
#                        prc     ::Float64,
#                        αc      ::Float64,
#                        σλc     ::Float64,
#                        α_prior ::NTuple{2,Float64},
#                        σλ_prior::NTuple{2,Float64},
#                        α_rdist ::NTuple{2,Float64},
#                        σλ_rdist::NTuple{2,Float64},
#                        nitpp   ::Int64,
#                        nthpp   ::Int64,
#                        βs      ::Vector{Float64},
#                        δt      ::Float64,
#                        srδt    ::Float64,
#                        inodes  ::Array{Int64,1},
#                        pup     ::Array{Int64,1})

#   # starting likelihood and prior
#   llc = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
#   prc = logdnorm(αc,         α_prior[1], α_prior[2]^2) +
#         logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])

#   K = lastindex(βs)

#   # make log-likelihood table per power
#   nlg = fld(nitpp, nthpp)
#   pp  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

#   L       = treelength(Ξ)      # tree length
#   dλ      = deltaλ(Ξ)          # delta change in λ
#   ssλ, nλ = sss_gbm(Ξ, αc)     # sum squares in λ
#   nin     = lastindex(inodes)  # number of internal nodes
#   el      = lastindex(idf)     # number of branches

#   for k in 2:K

#     βi  = βs[k]
#     rdc = logdnorm(       αc,  α_rdist[1], α_rdist[2]^2) +
#           logdinvgamma(σλc^2, σλ_rdist[1], σλ_rdist[2])

#     # logging
#     lth, lit = 0, 0

#     for it in Base.OneTo(nitpp)

#       shuffle!(pup)

#       for pupi in pup

#         ## parameter updates
#         # update drift
#         if pupi === 1

#           llc, prc, rdc, αc = update_α!(αc, σλc, L, dλ, llc, prc, rdc,
#             α_prior, α_rdist, βi)

#           # update ssλ with new drift `α`
#           ssλ, nλ = sss_gbm(Ξ, αc)

#         # update diffusion rate
#         elseif pupi === 2

#           llc, prc, rdc, σλc = update_σ!(σλc, ssλ, nλ, llc, prc, rdc,
#             σλ_prior, σλ_rdist, βi)

#         # update gbm
#         elseif pupi === 3

#           nix = ceil(Int64,rand()*nin)
#           bix = inodes[nix]

#           llc, dλ, ssλ =
#             update_gbm!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, δt, srδt)

#         # update by forward simulation
#         else
#           bix = ceil(Int64,rand()*el)

#           llc, dλ, ssλ, nλ, L =
#             update_fs!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, nλ, L, δt, srδt)

#         end
#       end

#       # log log-likelihood
#       lth += 1
#       if lth === nthpp
#         lit += 1
#         pp[k][lit] = llc + prc - rdc
#         lth = 0
#       end
#     end

#     @info string(βi," power done")
#   end

#   return pp
# end


