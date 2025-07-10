#=

GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#




"""
    insane_cladspb(tree    ::sT_label;
                 λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                 α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                 σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 nflush  ::Int64                 = nthin,
                 ofile   ::String                = string(homedir(), "/ipb"),
                 αi      ::Float64               = 0.0,
                 σλi     ::Float64               = 0.1,
                 pupdp   ::NTuple{5,Float64}     = (1e-3, 1e-3, 1e-3, 0.2, 0.2),
                 δt      ::Float64               = 1e-3,
                 prints  ::Int64                 = 5,
                 stn     ::Float64               = 0.5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for clads `pd`.
"""
function insane_cladspb(tree    ::sT_label;
                        λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                        α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                        σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                        niter   ::Int64                 = 1_000,
                        nthin   ::Int64                 = 10,
                        nburn   ::Int64                 = 200,
                        nflush  ::Int64                 = nthin,
                        ofile   ::String                = string(homedir(), "/cladspb"),
                        αi      ::Float64               = 0.0,
                        σλi     ::Float64               = 0.1,
                        pupdp   ::NTuple{5,Float64}     = (1e-3, 1e-3, 1e-3, 0.1, 0.2),
                        prints  ::Int64                 = 5,
                        stn     ::Float64               = 0.5,
                        tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)

  # turn to logarithmic terms
  λ0_prior = (log(λ0_prior[1]), 2*log(λ0_prior[2]))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ, Inf)

  # make a decoupled tree
  Ξ = make_Ξ(idf, λmle_cpb(tree), αi, σλi, cTpb)

  # if rm first speciation event (condition on observing the tree)
  rmλ = Float64(iszero(e(Ξ[1])))

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σ, 3: scale, 4: gbm, 5: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running pure-birth clads"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, ns, stn =
    mcmc_burn_cladspb(Ξ, idf, λ0_prior, α_prior, σλ_prior, nburn, αi, σλi, stn,
      rmλ, inodes, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cladspb(Ξ, idf, llc, prc, αc, σλc, ns, stn, λ0_prior, α_prior, σλ_prior,
      rmλ, inodes, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end



"""
    mcmc_burn_cladspb(Ξ       ::Vector{cTpb},
                    idf     ::Vector{iBffs},
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    stn     ::Float64,
                    rmλ     ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

MCMC burn-in chain for `pbd`.
"""
function mcmc_burn_cladspb(Ξ       ::Vector{cTpb},
                           idf     ::Vector{iBffs},
                           λ0_prior::NTuple{2,Float64},
                           α_prior ::NTuple{2,Float64},
                           σλ_prior::NTuple{2,Float64},
                           nburn   ::Int64,
                           αc      ::Float64,
                           σλc     ::Float64,
                           stn     ::Float64,
                           rmλ     ::Float64,
                           inodes  ::Array{Int64,1},
                           pup     ::Array{Int64,1},
                           prints  ::Int64)

  # starting likelihood and prior
  lλ0 = lλ(Ξ[1])
  llc = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ0 + prob_ρ(idf)
  prc = logdnorm(lλ0,       λ0_prior[1], λ0_prior[2])   +
        logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])   +
        logdnorm(αc,         α_prior[1],  α_prior[2]^2)

  L   = treelength(Ξ)                           # tree length
  nin = lastindex(inodes)                       # number of internal nodes
  el  = lastindex(idf)                          # number of branches
  ns  = sum(x -> Float64(d2(x) > 0), idf) - rmλ # number of speciation events in likelihood

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, nλ, irλ = _ss_ir_dd(Ξ, idf, lλ, αc)

  # for scale tuning
  ltn = 0
  lup = lac = 0.0

  pbar = Progress(nburn, dt = prints, desc = "burning mcmc...", barlen = 20)

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

        llc, prc, irλ, acc = 
          update_scale!(Ξ, idf, llc, prc, irλ, ns, stn, λ0_prior)

        lac += acc
        lup += 1.0

      # update gbm
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ddλ, ssλ, irλ =
          update_gbm!(bix, Ξ, idf, αc, σλc, llc, prc, ddλ, ssλ, irλ, 
            δt, srδt, λ0_prior)

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
    mcmc_cladspb(Ξ       ::Vector{cTpb},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               ns      ::Float64,
               λ0_prior::NTuple{2,Float64},
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
function mcmc_cladspb(Ξ       ::Vector{cTpb},
                      idf     ::Vector{iBffs},
                      llc     ::Float64,
                      prc     ::Float64,
                      αc      ::Float64,
                      σλc     ::Float64,
                      ns      ::Float64,
                      stn     ::Float64,
                      λ0_prior::NTuple{2,Float64},
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
  lthin = lit = sthin = zero(Int64)

  r = Array{Float64,2}(undef, nlogs, 6)

  L   = treelength(Ξ)      # tree length
  nin = lastindex(inodes)  # number of internal nodes
  el  = lastindex(idf)     # number of branches

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, nλ, irλ = _ss_ir_dd(Ξ, idf, lλ, αc)

  treev = cTpb[]  # make Ξ vector
  io = IOBuffer() # buffer 

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\talpha\tsigma_lambda\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αc = αc, σλc = σλc, ns = ns, nλ = nλ, ssλ = ssλ, ddλ = ddλ, irλ = irλ, L = L, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update drift
            if pupi === 1

              llc, prc, αc = update_α!(αc, σλc, nλ, ddλ, llc, prc, α_prior)

              # update ssλ with new drift `α`
              ssλ = _ss(Ξ, idf, lλ, αc)

              ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              if !isapprox(ll0, llc, atol = 1e-4)
                 @show ll0, llc, it, pupi
                 return
              end

            # update diffusion rate
            elseif pupi === 2

              llc, prc, σλc = update_σ!(σλc, ssλ, nλ, llc, prc, σλ_prior)

              ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              if !isapprox(ll0, llc, atol = 1e-4)
                 @show ll0, llc, it, pupi
                 return
              end

            # update scale
            elseif pupi === 3

              llc, prc, irλ, acc = 
                update_scale!(Ξ, idf, llc, prc, irλ, ns, stn, λ0_prior)

              ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              if !isapprox(ll0, llc, atol = 1e-4)
                 @show ll0, llc, it, pupi
                 return
              end

            # update internal λ
            elseif pupi === 4

              bix = inodes[fIrand(nin) + 1]

              llc, prc, ddλ, ssλ, irλ =
                update_internal!(bix, Ξ, idf, αc, σλc, llc, prc, ddλ, ssλ, irλ, 
                  λ0_prior)

              ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              if !isapprox(ll0, llc, atol = 1e-4)
                 @show ll0, llc, it, pupi
                 return
              end

            # update by forward simulation
            else

              bix = fIrand(el) + 1

              llc, ddλ, ssλ, nλ, irλ, ns, L =
                update_fs!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, nλ, irλ, ns, L, 
                  δt, srδt)

              ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              if !isapprox(ll0, llc, atol = 1e-4)
                 @show ll0, llc, it, pupi
                 return
              end

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
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
                  exp(lλ(Ξ[1])[1]),'\t', αc, '\t', σλc, '\n')
            flush(of)
            ibuffer(io, couple(Ξ, idf, 1))
            write(io, '\n')
            write(tf, take!(io))
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
    update_scale!(Ξ       ::Vector{T},
                  idf     ::Vector{iBffs},
                  llc     ::Float64,
                  prc     ::Float64,
                  ir      ::Float64,
                  ns      ::Float64,
                  stn     ::Float64,
                  λ0_prior::NTuple{2,Float64}) where {T <: iTree}

Update scale for speciation.
"""
function update_scale!(Ξ       ::Vector{cTpb},
                       idf     ::Vector{iBffs},
                       llc     ::Float64,
                       prc     ::Float64,
                       ir      ::Float64,
                       ns      ::Float64,
                       stn     ::Float64,
                       λ0_prior::NTuple{2,Float64})

  # sample log(scaling factor)
  s = randn()*stn

  # likelihood ratio
  iri = (1.0 - exp(s)) * ir
  llr = ns * s + iri

  lλ0 = lλ(Ξ[1])

  # prior ratio
  prr = llrdnorm_x(lλ0 + s, lλ0, λ0_prior[1], λ0_prior[2]) 

  acc = 0.0

  if -randexp() < llr + prr
    acc += 1.0
    llc += llr
    prc += prr
    ir  -= iri
    scale_rateλ!(Ξ, s)
    scale_rate!(idf, s)
  end

  return llc, prc, ir, acc
end




"""
    update_internal!(bix     ::Int64,
                     Ξ       ::Vector{cTpb},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ddλ     ::Float64,
                     ssλ     ::Float64,
                     irλ     ::Float64,
                     λ0_prior::NTuple{2,Float64})

Make a `gbm` update for an internal branch and its descendants.
"""
function update_internal!(bix     ::Int64,
                          Ξ       ::Vector{cTpb},
                          idf     ::Vector{iBffs},
                          α       ::Float64,
                          σλ      ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          ddλ     ::Float64,
                          ssλ     ::Float64,
                          irλ     ::Float64,
                          λ0_prior::NTuple{2,Float64})

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  i2   = d2(bi)
  ia   = pa(bi)
  ξ1   = Ξ[i1]
  ξ2   = Ξ[i2]
  root = iszero(ia)
  λa   = NaN  # ancestral speciation

  # if crown root
  if root && iszero(e(ξi))
    llc, prc, ddλ, ssλ =
      _crown_update!(ξi, ξ1, ξ2, α, σλ, llc, prc, ddλ, ssλ, λ0_prior)
    setλt!(bi, lλ(ξi))
    λa = lλ(ξi)
  else
    # if stem
    if root
      if istip(ξi)
        llc, prc, ddλ, ssλ, irλ = 
          _stem_update!(ξi, lλ(ξ1), lλ(ξ2), 
            α, σλ, llc, prc, ddλ, ssλ, irλ, λ0_prior)
        λa = lλ(ξi)
      else
        llc, prc, ddλ, ssλ, irλ = 
          _stem_update!(ξi, lλ(ξ1.d1), lλ(ξ2.d2),
            α, σλ, llc, prc, ddλ, ssλ, irλ, λ0_prior)

        # updates within the stem daughter branches
        llc, ddλ, ssλ, irλ = 
          _update_internal!(ξi.d1, lλ(ξi), α, σλ, llc, ddλ, ssλ, irλ, false)
        llc, ddλ, ssλ, irλ = 
          _update_internal!(ξi.d2, lλ(ξi), α, σλ, llc, ddλ, ssλ, irλ, false)

        # get fixed tip and ancestral rate
        lξi, λa = fixtip(ξi, λa)

        # make fixed branch update
        llc, ddλ, ssλ, irλ, λa = 
          update_triad!(lξi, ξ1, ξ2, λa, α, σλ, llc, ddλ, ssλ, irλ)
      end
    else
      λa = λt(idf[ia])

      # updates within the parent branch
      llc, ddλ, ssλ, irλ = 
        _update_internal!(ξi, λa, α, σλ, llc, ddλ, ssλ, irλ, false)

      # get fixed tip and ancestral rate
      lξi, λa = fixtip(ξi, λa)

      # make fixed branch update
      llc, ddλ, ssλ, irλ, λa = 
        update_triad!(lξi, ξ1, ξ2, λa, α, σλ, llc, ddλ, ssλ, irλ)
    end

    # set fixed `λ(t)` in branch
    setλt!(bi, λa)
  end

  # # carry on updates in the daughters
  llc, ddλ, ssλ, irλ = 
    _update_internal!(ξ1, λa, α, σλ, llc, ddλ, ssλ, irλ, iszero(d1(idf[i1])))
  llc, ddλ, ssλ, irλ = 
    _update_internal!(ξ2, λa, α, σλ, llc, ddλ, ssλ, irλ, iszero(d1(idf[i2])))

  return llc, prc, ddλ, ssλ, irλ
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{cTpb},
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
                    Ξ    ::Vector{cTpb},
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

    # update quantities
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
           ξc  ::cTpb,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::cTpb,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  t0, nap, nn, llr =
    _sim_cladspb_t(e(bi), lλ(ξc)[1], α, σλ, δt, srδt, lc, lU, iρi, 0, 1, 500)

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
           ξ1  ::cTpb,
           ξ2  ::cTpb,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::cTpb,
                ξ1  ::cTpb,
                ξ2  ::cTpb,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  t0, na = _sim_cladspb(e(bi), lλ(ξc)[1], α, σλ, δt, srδt, 1, 1_000)

  if na >= 1_000
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

 # fix random tip
  λf = fixrtip!(t0, na, NaN)

  llrd, acrd, drλ, ssrλ, irrλ, λ1p, λ2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulated remaining tips until the present
    t0, na, acr =
      tip_sims!(t0, tf(bi), α, σλ, δt, srδt, acr, lU, iρi, na)

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
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
    tip_sims!(tree::cTpb,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::cTpb,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64)

 if lU < lr && na < 1_000

    if istip(tree)
      if !isfix(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)

        # simulate
        stree, na, lr =
          _sim_cladspb_it(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, δt, srδt,
            lr, lU, iρi, na, 1_000)

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
      tree.d1, na, lr = tip_sims!(tree.d1, t, α, σλ, δt, srδt, lr, lU, iρi, na)
      tree.d2, na, lr = tip_sims!(tree.d2, t, α, σλ, δt, srδt, lr, lU, iρi, na)
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





