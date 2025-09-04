#=

GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#




"""
    insane_cladsb(tree    ::sT_label;
                  λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                  α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                  σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                  niter   ::Int64                 = 1_000,
                  nthin   ::Int64                 = 10,
                  nburn   ::Int64                 = 200,
                  nflush  ::Int64                 = nthin,
                  ofile   ::String                = string(homedir(), "/cladsb"),
                  αi      ::Float64               = 0.0,
                  σλi     ::Float64               = 0.1,
                  pupdp   ::NTuple{5,Float64}     = (1e-3, 1e-3, 1e-4, 0.1, 0.2),
                  prints  ::Int64                 = 5,
                  stn     ::Float64               = 0.5,
                  tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for clads `pd`.
"""
function insane_cladsb(tree    ::sT_label;
                       λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                       α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                       σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                       niter   ::Int64                 = 1_000,
                       nthin   ::Int64                 = 10,
                       nburn   ::Int64                 = 200,
                       nflush  ::Int64                 = nthin,
                       ofile   ::String                = string(homedir(), "/cladsb"),
                       αi      ::Float64               = 0.0,
                       σλi     ::Float64               = 0.1,
                       pupdp   ::NTuple{5,Float64}     = (1e-3, 1e-3, 1e-4, 0.1, 0.2),
                       prints  ::Int64                 = 5,
                       stn     ::Float64               = 0.5,
                       tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)

  # turn to logarithmic terms
  λ0_prior = (log(λ0_prior[1]), 2.0*log(λ0_prior[2]))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ, Inf)

  # make a decoupled tree
  Ξ = make_Ξ(idf, λmle_cb(tree), cTb)

  # if rm first speciation event (condition on observing the tree)
  rmλ = Float64(iszero(e(Ξ[1])))

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σ, 3: scale, 4: internal, 5: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running pure-birth clads (μ(t) = 0)"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, ns, stn =
    mcmc_burn_cladsb(Ξ, idf, λ0_prior, α_prior, σλ_prior, nburn, αi, σλi, stn,
      rmλ, inodes, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cladsb(Ξ, idf, llc, prc, αc, σλc, ns, stn, rmλ, λ0_prior, α_prior, 
      σλ_prior, inodes, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end



"""
    mcmc_burn_cladsb(Ξ       ::Vector{cTb},
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
function mcmc_burn_cladsb(Ξ       ::Vector{cTb},
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

  nin = lastindex(inodes)                       # number of internal nodes
  el  = lastindex(idf)                          # number of branches
  ns  = sum(x -> Float64(d2(x) > 0), idf) - rmλ # number of speciation events in likelihood
  λfs = Float64[]

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ = _dd_ss(Ξ, idf, αc)

  # for scale tuning
  ltn = zero(Int64)
  lup = lac = zero(Float64)

  pbar = Progress(nburn, dt = prints, desc = "burning mcmc...", barlen = 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update drift
      if pupi === 1

        llc, prc, αc = 
          update_α!(αc, σλc, 2.0*(ns + rmλ), ddλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ = _ss(Ξ, idf, αc)

      # update diffusion
      elseif pupi === 2

        llc, prc, σλc = 
          update_σ!(σλc, ssλ, 2.0*(ns + rmλ), llc, prc, σλ_prior)

      # update scale
      elseif pupi === 3

        llc, prc, acc = 
          update_scale!(Ξ, idf, llc, prc, ns, stn, λ0_prior)

        lac += acc
        lup += 1.0

      # update gbm
      elseif pupi === 4

        bix = inodes[fIrand(nin) + 1]

        llc, prc, ddλ, ssλ =
          update_internal!(bix, Ξ, idf, αc, σλc, llc, prc, ddλ, ssλ, 
            λ0_prior)

      # forward simulation
      else

        bix = fIrand(el) + 1

        llc, ddλ, ssλ, ns =
          update_fs!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, ns, λfs)
      end
    end

    ltn += one(Int64)
    if ltn === 100
      stn = tune(stn, lac/lup)
      ltn = zero(Int64)
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, ns, stn
end




"""
    mcmc_cladsb(Ξ       ::Vector{cTb},
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
function mcmc_cladsb(Ξ       ::Vector{cTb},
                      idf     ::Vector{iBffs},
                      llc     ::Float64,
                      prc     ::Float64,
                      αc      ::Float64,
                      σλc     ::Float64,
                      ns      ::Float64,
                      stn     ::Float64,
                      rmλ     ::Float64,
                      λ0_prior::NTuple{2,Float64},
                      α_prior ::NTuple{2,Float64},
                      σλ_prior::NTuple{2,Float64},
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

  r   = Array{Float64,2}(undef, nlogs, 6)

  nin = lastindex(inodes)  # number of internal nodes
  el  = lastindex(idf)     # number of branches

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ = _dd_ss(Ξ, idf, αc)

  λfs   = Float64[]
  treev = cTb[]  # make Ξ vector
  io = IOBuffer() # buffer 

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\talpha\tsigma_lambda\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αc = αc, σλc = σλc, ns = ns, ssλ = ssλ, ddλ = ddλ, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update drift
            if pupi === 1

              llc, prc, αc = 
                update_α!(αc, σλc, 2.0*(ns + rmλ), ddλ, llc, prc, α_prior)

              # update ssλ with new drift `α`
              ssλ = _ss(Ξ, idf, αc)

              # ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update diffusion rate
            elseif pupi === 2

              llc, prc, σλc = 
                update_σ!(σλc, ssλ, 2.0*(ns + rmλ), llc, prc, σλ_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update scale
            elseif pupi === 3

              llc, prc, acc = 
                update_scale!(Ξ, idf, llc, prc, ns, stn, λ0_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update internal λ
            elseif pupi === 4

              bix = inodes[fIrand(nin) + 1]

              llc, prc, ddλ, ssλ =
                update_internal!(bix, Ξ, idf, αc, σλc, llc, prc, ddλ, ssλ, 
                  λ0_prior)

              # ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update by forward simulation
            else

              bix = fIrand(el) + 1

              llc, ddλ, ssλ, ns =
                update_fs!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, ns, λfs)

              # ll0 = llik_clads(Ξ, idf, αc, σλc) - rmλ*lλ(Ξ[1]) + prob_ρ(idf)
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
              r[lit,4] = exp(lλ(Ξ[1]))
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
                  exp(lλ(Ξ[1])),'\t', αc, '\t', σλc, '\n')
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
function update_scale!(Ξ       ::Vector{cTb},
                       idf     ::Vector{iBffs},
                       llc     ::Float64,
                       prc     ::Float64,
                       ns      ::Float64,
                       stn     ::Float64,
                       λ0_prior::NTuple{2,Float64})

  # sample log(scaling factor)
  s = randn()*stn

  # likelihood ratio
  ir  = _ir(Ξ)
  llr = ns * s + (1.0 - exp(s)) * ir

  lλ0 = lλ(Ξ[1])

  # prior ratio
  prr = llrdnorm_x(lλ0 + s, lλ0, λ0_prior[1], λ0_prior[2]) 

  acc = 0.0

  if -randexp() < llr + prr
    acc += 1.0
    llc += llr
    prc += prr
    scale_rate!(Ξ,   addlλ!, s)
    scale_rate!(idf, addlλ!, s)
  end

  return llc, prc, acc
end




"""
    update_internal!(bix     ::Int64,
                     Ξ       ::Vector{cTb},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ddλ     ::Float64,
                     ssλ     ::Float64,
                     λ0_prior::NTuple{2,Float64})

Make a `gbm` update for an internal branch and its descendants.
"""
function update_internal!(bix     ::Int64,
                          Ξ       ::Vector{cTb},
                          idf     ::Vector{iBffs},
                          α       ::Float64,
                          σλ      ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          ddλ     ::Float64,
                          ssλ     ::Float64,
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
    λa = lλ(ξi)
    setλt!(bi, λa)
  else
    # if stem
    if root
      if istip(ξi)
        llc, prc, ddλ, ssλ = 
          _stem_update!(ξi, lλ(ξ1), lλ(ξ2), 
            α, σλ, llc, prc, ddλ, ssλ, λ0_prior)
        λa = lλ(ξi)
      else
        llc, prc, ddλ, ssλ = 
          _stem_update!(ξi, lλ(ξ1.d1), lλ(ξ2.d2),
            α, σλ, llc, prc, ddλ, ssλ, λ0_prior)

        # updates within the stem daughter branches
        llc, ddλ, ssλ = 
          _update_internal!(ξi.d1, lλ(ξi), α, σλ, llc, ddλ, ssλ, false)
        llc, ddλ, ssλ = 
          _update_internal!(ξi.d2, lλ(ξi), α, σλ, llc, ddλ, ssλ, false)

        # get fixed tip and ancestral rate
        lξi, λa = fixtip(ξi, λa)

        # make fixed branch update
        llc, ddλ, ssλ, λa = 
          update_triad!(lξi, ξ1, ξ2, λa, α, σλ, llc, ddλ, ssλ)
      end
    else
      λa = λt(idf[ia])

      # updates within the parent branch
      llc, ddλ, ssλ, λx = 
        _update_internal!(ξi, λa, α, σλ, llc, ddλ, ssλ, false)

      # get fixed tip and ancestral rate
      lξi, λa = fixtip(ξi, λa)

      # make fixed branch update
      llc, ddλ, ssλ, λa = 
        update_triad!(lξi, ξ1, ξ2, λa, α, σλ, llc, ddλ, ssλ)
    end

    # set fixed `λ(t)` in branch
    setλt!(bi, λa)
  end

  # # carry on updates in the daughters
  llc, ddλ, ssλ, λx = 
    _update_internal!(ξ1, λa, α, σλ, llc, ddλ, ssλ, iszero(d1(idf[i1])))
  llc, ddλ, ssλ, λx = 
    _update_internal!(ξ2, λa, α, σλ, llc, ddλ, ssλ, iszero(d1(idf[i2])))

  return llc, prc, ddλ, ssλ
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{cTb},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               llc  ::Float64,
               ddλ  ::Float64,
               ssλ  ::Float64,
               ns   ::Float64,
               λfs  ::Vector{Float64})

Forward simulation proposal function for pure birth diffusion.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{cTb},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    llc  ::Float64,
                    ddλ  ::Float64,
                    ssλ  ::Float64,
                    ns   ::Float64,
                    λfs  ::Vector{Float64})

  bi = idf[bix]
  ξc = Ξ[bix]
  ia = pa(bi)

  λa = NaN
  if ia > 0
    λa = λt(idf[ia])
  end

  # if terminal
  ssλr = ddλr = zero(Float64)
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, λa, α, σλ)
  # if internal
  else
    ξp, llr, ddλr, ssλr = 
      fsbi_i(bi, ξc, λa, lλ(Ξ[d1(bi)]), lλ(Ξ[d2(bi)]), α, σλ, λfs)
  end

  # if accepted
  if isfinite(llr)

    llc, ddλ, ssλ, ns = 
      llik_cladsb_track!(ξc, α, σλ, llc, ddλ, ssλ, ns, -)
    llc, ddλ, ssλ, ns = 
      llik_cladsb_track!(ξp, α, σλ, llc, ddλ, ssλ, ns, +)

    # first change from ancestor
    if ia > 0
      λp, λc = lλ(ξp), lλ(ξc)
      llc += llrdnorm_x(λp, λc, λa + α, σλ^2)
      ddλ += λp - λc
      ssλ += 0.5*((λp - λa - α)^2 - (λc - λa - α)^2)
    end

    # update quantities
    ddλ += ddλr
    ssλ += ssλr
    llc += llr

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ssλ, ns
end




"""
    fsbi_t(bi  ::iBffs,
           ξc  ::cTb,
           λa  ::Float64,
           α   ::Float64,
           σλ  ::Float64)

Forward simulation for terminal branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                λa  ::Float64,
                α   ::Float64,
                σλ  ::Float64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  λi = rnorm(λa + α, σλ)

  t0, nap, nn, llr =
    _sim_cladsb_t(e(bi), λi, α, σλ, lc, lU, iρi, 0, 1, 500)

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
           ξc  ::cTb,
           λa  ::Float64,
           λ1  ::Float64,
           λ2  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           λfs ::Vector{Float64})

Forward simulation for internal branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::cTb,
                λa  ::Float64,
                λ1  ::Float64,
                λ2  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                λfs ::Vector{Float64})

  if isnan(λa)
    λi = lλ(ξc)
  else
    λi = rnorm(λa + α, σλ)
  end

  empty!(λfs)

  # forward simulation during branch length
  t0, na = _sim_cladsb_i(e(bi), λi, α, σλ, 1, 500, λfs)

  if na > 499
    return t0, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr  = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # choose most likely lineage to fix
  wt, λp, pp, λc, pc, acr, ddr, ssr = 
    wfix_i(ξc, e(bi), λfs, λ1, λ2, α, σλ, acr)

  if lU < acr

    # fix the tip
    if wt <= div(na,2)
      fixtip1!(t0, wt, 0)
    else
      fixtip2!(t0, na - wt + 1, 0)
    end

    # simulated remaining tips until the present
    t0, na, acr =
      tip_sims!(t0, tf(bi), α, σλ, acr, lU, iρi, na)

    if lU < acr
      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc) + λp - λc
      setni!(bi, na)                     # set new ni
      setλt!(bi, λp)                     # set new λt

      return t0, llr, ddr, ssr
    end
  end

  return t0, NaN, NaN, NaN
end




"""
    wfix_i(ξi ::cTb,
           ei ::Float64,
           λfs::Vector{Float64},
           λ1 ::Float64,
           λ2 ::Float64,
           α  ::Float64,
           σλ ::Float64,
           ss ::Float64,
           dd ::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for bifurcating `i` branches.
"""
function wfix_i(ξi ::cTb,
                ei ::Float64,
                λfs::Vector{Float64},
                λ1 ::Float64,
                λ2 ::Float64,
                α  ::Float64,
                σλ ::Float64,
                acr::Float64)

  # select best from proposal
  sp, i, wt, λp, pp = 0.0, 0, 0, NaN, -Inf
  for λfi in λfs
    i += 1
    p   = dnorm2(λ1, λ2, λfi + α, σλ)
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      wt  = i
    end
  end

  # extract current xis and estimate ratio
  empty!(λfs)
  λc = _λat!(ξi, ei, λfs, 0.0, NaN)

  sc, pc = 0.0, NaN
  for λfi in λfs
    p   = dnorm2(λ1, λ2, λfi + α, σλ)
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood and acceptance ratio
  acr += log(sp/sc) + λp - λc
  ddr  = 2.0*(λc - λp)
  ssr  = 0.5*((λ1 - λp - α)^2 + (λ2 - λp - α)^2 - 
              (λ1 - λc - α)^2 - (λ2 - λc - α)^2)

  return wt, λp, pp, λc, pc, acr, ddr, ssr
end




"""
    tip_sims!(tree::cTb,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::cTb,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64)

 if lU < lr && na < 500

    if istip(tree)
      if !isfix(tree)

        # simulate
        stree, na, lr =
          _sim_cladsb_it(t, lλ(tree), α, σλ, lr, lU, iρi, na, 500)

        if isnan(lr) || na > 499
          return tree, na, NaN
        end

        sete!(tree, e(tree) + e(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, lr = tip_sims!(tree.d1, t, α, σλ, lr, lU, iρi, na)
      tree.d2, na, lr = tip_sims!(tree.d2, t, α, σλ, lr, lU, iρi, na)
    end

    return tree, na, lr
  end

  return tree, na, NaN
end





