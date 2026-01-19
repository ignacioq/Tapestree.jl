#=

Anagenetic `gbmbd` MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    insane_gbmbd(tree    ::sT_label,
                 tv      ::Vector{Vector{Float64}},
                 ev      ::Vector{Vector{Float64}};
                 λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                 μ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                 α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                 σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                 σμ_prior::NTuple{2,Float64}     = (3.0, 0.1),
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 nflush  ::Int64                 = nthin,
                 ofile   ::String                = string(homedir(), "/ibd_efx"),
                 ϵi      ::Float64               = 0.2,
                 λi      ::Float64               = NaN,
                 αi      ::Float64               = 0.0,
                 σλi     ::Float64               = 0.01,
                 σμi     ::Float64               = 0.01,
                 pupdp   ::NTuple{4,Float64}     = (1e-3, 1e-3, 0.1, 0.2),
                 δt      ::Float64               = 1e-3,
                 prints  ::Int64                 = 5,
                 survival::Bool                  = true,
                 mxthf   ::Float64               = 0.1,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `gbm-bd` with fixed extinction.
"""
function insane_gbmbd(tree    ::sT_label,
                      tv      ::Vector{Vector{Float64}},
                      ev      ::Vector{Vector{Float64}};
                      λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                      μ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                      α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                      σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                      σμ_prior::NTuple{2,Float64}     = (3.0, 0.1),
                      niter   ::Int64                 = 1_000,
                      nthin   ::Int64                 = 10,
                      nburn   ::Int64                 = 200,
                      nflush  ::Int64                 = nthin,
                      ofile   ::String                = string(homedir(), "/ibd_efx"),
                      ϵi      ::Float64               = 0.2,
                      λi      ::Float64               = NaN,
                      αi      ::Float64               = 0.0,
                      σλi     ::Float64               = 0.01,
                      σμi     ::Float64               = 0.01,
                      pupdp   ::NTuple{4,Float64}     = (1e-3, 1e-3, 0.1, 0.2),
                      δt      ::Float64               = 1e-3,
                      prints  ::Int64                 = 5,
                      survival::Bool                  = true,
                      mxthf   ::Float64               = 0.1,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n     = ntips(tree)
  th    = treeheight(tree)
  δt   *= max(0.1,round(th, RoundDown, digits = 2))
  srδt  = sqrt(δt)

  surv = 0   # condition on survival of 0, 1, or 2 starting lineages
  rmλ  = 0.0 # condition on first speciation event
  if iszero(e(tree)) 
    rmλ  += 1.0
    surv += survival ? 2 : 0
  else
    surv += survival ? 1 : 0
  end

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # sort extinction vectors 
  ev = deepcopy(ev)
  tv = deepcopy(tv)
  for i in Base.OneTo(lastindex(tv))
    io = sortperm(tv[i], rev = true)
    tv[i] = tv[i][io]
    ev[i] = log.(ev[i][io])
  end

  # estimate branch split (multiple of δt)
  ndts = floor(th * mxthf/δt)
  maxt = δt * ndts

  # make fix tree directory
  idf = make_idf(tree, tρ, maxt)

  # make one per unit in `idf`
  ixi = nnodesinternal(tree)*2 + 2
  ixs = Int64[]
  for bi in reverse(idf)
    ixi -= ismid(bi) ? 0 : 1 
    push!(ixs, ixi)
  end
  reverse!(ixs)

  tv = tv[ixs]
  ev = ev[ixs]

  # find initial extinction rate
  ix = findfirst(x -> x < th, tv[1]) - 1
  μi = exp(linpred(th, tv[1][ix], tv[1][ix+1], ev[1][ix], ev[1][ix+1]))

  # starting parameters (using method of moments)
  if isnan(λi)
    λc, μc = μi/ϵi, μi
  else
    λc, μc = λi, μi
  end

  # make a decoupled tree
  Ξ, ixiv, ixfv = make_Ξ(idf, λc, αi, σλi, tv, ev, δt, srδt, iTbd)

  # survival
  mc = m_surv_gbmbd(th, log(λc), log(μc), αi, σλi, σμi, δt, srδt, 1_000, surv)

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf))  if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σλ, 3: σμ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death diffusion with fixed extinction function"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, σμc, mc =
    mcmc_burn_gbmbd(Ξ, idf, λ0_prior, μ0_prior, α_prior, σλ_prior, σμ_prior,
      nburn, αi, σλi, σμi, mc, tv, ev, ixiv, ixfv, th, rmλ, surv, 
      δt, srδt, inodes, pup, prints)

  # mcmc
  r, treev =
    mcmc_gbmbd(Ξ, idf, llc, prc, αc, σλc, σμc, mc, tv, ev, ixiv, ixfv, 
      th, rmλ, surv, λ0_prior, α_prior, σλ_prior, σμ_prior, δt, srδt, 
      inodes, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_gbmbd(Ξ       ::Vector{iTbd},
                    idf     ::Vector{iBffs},
                    λ0_prior::NTuple{2,Float64},
                    μ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    mc      ::Float64,
                    tv      ::Vector{Vector{Float64}},
                    ev      ::Vector{Vector{Float64}},
                    ixiv    ::Array{Int64,1},
                    ixfv    ::Array{Int64,1},
                    th      ::Float64,
                    crown   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

MCMC burn-in chain for `gbmbd`.
"""
function mcmc_burn_gbmbd(Ξ       ::Vector{iTbd},
                         idf     ::Vector{iBffs},
                         λ0_prior::NTuple{2,Float64},
                         μ0_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         σμ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         σμc     ::Float64,
                         mc      ::Float64,
                         tv      ::Vector{Vector{Float64}},
                         ev      ::Vector{Vector{Float64}},
                         ixiv    ::Array{Int64,1},
                         ixfv    ::Array{Int64,1},
                         th      ::Float64,
                         rmλ     ::Float64,
                         surv    ::Int64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  lλ0 = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - rmλ * lλ0 +
        log(mc) + prob_ρ(idf)
  prc = logdnorm(lλ0,         λ0_prior[1], λ0_prior[2])   +
        logdnorm(lμ(Ξ[1])[1], μ0_prior[1], μ0_prior[2])   +
        logdinvgamma(σλc^2,   σλ_prior[1], σλ_prior[2])   +
        logdinvgamma(σμc^2,   σμ_prior[1], σμ_prior[2])   +
        logdnorm(αc,           α_prior[1],  α_prior[2]^2)

  L   = treelength(Ξ)      # tree length
  nin = lastindex(inodes) # number of internal nodes
  el  = lastindex(idf)    # number of branches

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, ssμ, nλ = _dd_ss(Ξ, αc)

  pbar = Progress(nburn, dt = prints, desc = "burning mcmc...", barlen = 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update α
      if pupi === 1

        llc, prc, αc, mc  =
          update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, L, ddλ, llc, prc,
            mc, th, surv, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ = _ss(Ξ, αc)

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nλ,
            llc, prc, mc, th, surv, δt, srδt, σλ_prior, σμ_prior)

      # gbm update
      elseif pupi === 3

        bix = inodes[fIrand(nin) + 1]

        llc, prc, ddλ, ssλ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, prc, ddλ, ssλ,
            mc, th, δt, srδt, λ0_prior, surv)

      # forward simulation update
      else

        bix = fIrand(el) + 1

        llc, ddλ, ssλ, ssμ, nλ, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, llc, ddλ, ssλ, ssμ, nλ, L,
            ixiv, ixfv, tv[bix], ev[bix], δt, srδt)

      end
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, σμc, mc
end




"""
    mcmc_gbmbd(Ξ       ::Vector{iTbd},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               σμc     ::Float64,
               mc      ::Float64,
               tv      ::Vector{Vector{Float64}},
               ev      ::Vector{Vector{Float64}},
               ixiv    ::Array{Int64,1},
               ixfv    ::Array{Int64,1},
               th      ::Float64,
               rmλ     ::Float64,
               surv    ::Int64,
               λ0_prior::NTuple{2,Float64},
               μ0_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               σμ_prior::NTuple{2,Float64},
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Vector{Int64},
               niter   ::Int64,
               nthin   ::Int64,
               nflush  ::Int64,
               ofile   ::String,
               prints  ::Int64)

MCMC chain for `gbmbd`.
"""
function mcmc_gbmbd(Ξ       ::Vector{iTbd},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    mc      ::Float64,
                    tv      ::Vector{Vector{Float64}},
                    ev      ::Vector{Vector{Float64}},
                    ixiv    ::Array{Int64,1},
                    ixfv    ::Array{Int64,1},
                    th      ::Float64,
                    rmλ     ::Float64,
                    surv    ::Int64,
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
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

  L   = treelength(Ξ)     # tree length
  nin = lastindex(inodes) # number of internal nodes
  el  = lastindex(idf)    # number of branches

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, ssμ, nλ = _dd_ss(Ξ, αc)

  # parameter results
  r = Array{Float64,2}(undef, nlogs, 8)

  treev = iTbd[]          # make tree vector
  io    = IOBuffer()      # buffer 

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\tmu_root\talpha\tsigma_lambda\tsigma_mu\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αc = αc, σλc = σλc, σμc = σμc, mc = mc, nλ = nλ, ssλ = ssλ, ssμ = ssμ, ddλ = ddλ, L = L, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          # parameter updates
          for pupi in pup

            # update α
            if pupi === 1

              llc, prc, αc, mc  =
                update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, L, ddλ, llc, prc,
                  mc, th, surv, δt, srδt, α_prior)

              # update ssλ with new drift `α`
              ssλ, ssμ = _ss(Ξ, αc)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - Float64(surv > 1) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              #  if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, i, pupi, Ξ
              #    return
              # end

            # σλ & σμ update
            elseif pupi === 2

              llc, prc, σλc, σμc, mc =
                update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αc, ssλ, ssμ, nλ,
                  llc, prc, mc, th, surv, δt, srδt, σλ_prior, σμ_prior)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - Float64(surv > 1) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              #  if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, i, pupi, Ξ
              #    return
              # end

            # gbm update
            elseif pupi === 3

              bix = inodes[fIrand(nin) + 1]

              llc, prc, ddλ, ssλ, mc =
                update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, prc, ddλ, ssλ,
                  mc, th, δt, srδt, λ0_prior, surv)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - Float64(surv > 1) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              #  if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, i, pupi, Ξ
              #    return
              # end

            # forward simulation update
            else

              bix = fIrand(el) + 1

              llc, ddλ, ssλ, ssμ, nλ, L =
                update_fs!(bix, Ξ, idf, αc, σλc, σμc, llc, ddλ, ssλ, ssμ, nλ, L,
                  ixiv, ixfv, tv[bix], ev[bix], δt, srδt)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, δt, srδt) - Float64(surv > 1) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              #  if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, i, pupi, Ξ
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
              r[lit,5] = exp(lμ(Ξ[1])[1])
              r[lit,6] = αc
              r[lit,7] = σλc
              r[lit,8] = σμc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
                 exp(lλ(Ξ[1])[1]),'\t', exp(lμ(Ξ[1])[1]), '\t', αc, '\t',
                 σλc, '\t', σμc,'\n')
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
    update_fs!(bix    ::Int64,
               Ξ      ::Vector{iTbd},
               idf    ::Vector{iBffs},
               α      ::Float64,
               σλ     ::Float64,
               σμ     ::Float64,
               llc    ::Float64,
               ddλ    ::Float64,
               ssλ    ::Float64,
               ssμ    ::Float64,
               nλ     ::Float64,
               L      ::Float64,
               ixiv   ::Vector{Int64},
               ixfv   ::Vector{Int64},
               tv     ::Vector{Float64},
               ev     ::Vector{Float64},
               δt     ::Float64,
               srδt   ::Float64)

Forward simulation proposal function for `gbmbd`.
"""
function update_fs!(bix    ::Int64,
                    Ξ      ::Vector{iTbd},
                    idf    ::Vector{iBffs},
                    α      ::Float64,
                    σλ     ::Float64,
                    σμ     ::Float64,
                    llc    ::Float64,
                    ddλ    ::Float64,
                    ssλ    ::Float64,
                    ssμ    ::Float64,
                    nλ     ::Float64,
                    L      ::Float64,
                    ixiv   ::Vector{Int64},
                    ixfv   ::Vector{Int64},
                    tv     ::Vector{Float64},
                    ev     ::Vector{Float64},
                    δt     ::Float64,
                    srδt   ::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]
  ixi = ixiv[bix]
  ixf = ixfv[bix]

  # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, ixi, ξc, α, σλ, tv, ev, δt, srδt)
    drλ = ssrλ = 0.0

  # if mid
  elseif iszero(d2(bi))
    ξp, llr, drλ, ssrλ =
      fsbi_m(bi, ixi, ixf, ξc, Ξ[d1(bi)], α, σλ, tv, ev, δt, srδt)

  # if internal
  else
    ξp, llr, drλ, ssrλ  =
      fsbi_i(bi, ixi, ixf, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, tv, ev, δt, srδt)
  end

  # if accepted
  if isfinite(llr)
    ll1, ddλ1, ssλ1, ssμ1, nλ1 = llik_gbm_ss(ξp, α, σλ, σμ, δt, srδt)
    ll0, ddλ0, ssλ0, ssμ0, nλ0 = llik_gbm_ss(ξc, α, σλ, σμ, δt, srδt)

    # update llr, ssλ, nλ, sns, ne, L,
    llc += ll1  - ll0  + llr
    ddλ  += ddλ1  - ddλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    ssμ += ssμ1 - ssμ0
    nλ  += nλ1  - nλ0
    L   += treelength(ξp)   - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ssλ, ssμ, nλ, L
end




"""
    fsbi_t(bi  ::iBffs,
           ix  ::Int64,
           ξc  ::iTbd,
           α   ::Float64,
           σλ  ::Float64,
           tv  ::Vector{Float64},
           ev  ::Vector{Float64},
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_t(bi  ::iBffs,
                ix  ::Int64,
                ξc  ::iTbd,
                α   ::Float64,
                σλ  ::Float64,
                tv  ::Vector{Float64},
                ev  ::Vector{Float64},
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_gbmbd_t(e(bi), lλ(ξc)[1], α, σλ, ix, tv, ev, δt, srδt, 
      lc, lU, Iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr)
    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, -Inf
  end
end



"""
    fsbi_m(bi  ::iBffs,
           ix  ::Int64,
           ξc  ::iTbd,
           ξ1  ::iTbd,
           α   ::Float64,
           σλ  ::Float64,
           tv  ::Vector{Float64},
           ev  ::Vector{Float64},
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_m(bi  ::iBffs,
                ixi ::Int64,
                ixf ::Int64,
                ξc  ::iTbd,
                ξ1  ::iTbd,
                α   ::Float64,
                σλ  ::Float64,
                tv  ::Vector{Float64},
                ev  ::Vector{Float64},
                δt  ::Float64,
                srδt::Float64)

  t0, na, nn =
    _sim_gbmbd_i(ti(bi), tf(bi), lλ(ξc)[1], α, σλ, ixi, tv, ev, 
      δt, srδt, 0, 1, 1_000)

  if na < 1 || nn > 999
    return t0, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample and fix random  tip
  λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, λ1p =
    _daughter_update!(ξ1, λf, α, σλ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, ixf, tv, ev, 
          δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      l1  = lastindex(λ1p)
      setnt!(bi, ntp)                       # set new nt
      setni!(bi, na)                        # set new ni
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector

      return t0, llr, drλ, ssrλ
    end
  end

  return t0, NaN, NaN, NaN
end






"""
    fsbi_i(bi  ::iBffs,
           ix  ::Int64,
           ξc  ::iTbd,
           ξ1  ::iTbd,
           ξ2  ::iTbd,
           α   ::Float64,
           σλ  ::Float64,
           tv  ::Vector{Float64},
           ev  ::Vector{Float64},
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ixi ::Int64,
                ixf ::Int64,
                ξc  ::iTbd,
                ξ1  ::iTbd,
                ξ2  ::iTbd,
                α   ::Float64,
                σλ  ::Float64,
                tv  ::Vector{Float64},
                ev  ::Vector{Float64},
                δt  ::Float64,
                srδt::Float64)

  t0, na, nn =
    _sim_gbmbd_i(ti(bi), tf(bi), lλ(ξc)[1], α, σλ, ixi, tv, ev, 
      δt, srδt, 0, 1, 1_000)

  if na < 1 || nn > 999
    return t0, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample and fix random  tip
  λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, λ1p, λ2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, ixf, tv, ev, 
          δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      l1  = lastindex(λ1p)
      l2  = lastindex(λ2p)
      setnt!(bi, ntp)                       # set new nt
      setni!(bi, na)                        # set new ni
      setλt!(bi, λf)                        # set new λt
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 2 λ vector

      return t0, llr, drλ, ssrλ
    end
  end

  return t0, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTbd,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTbd,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   ix  ::Int64,
                   tv  ::Vector{Float64},
                   ev  ::Vector{Float64},
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

  if lU < lr && nn < 1_000

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)
        lμ0  = lμ(tree)
        l    = lastindex(lλ0)

        # simulate
        stree, na, nn, lr =
          _sim_gbmbd_it(max(δt-fdti, 0.0), t, lλ0[l], α, σλ, ix, tv, ev, 
            δt, srδt, lr, lU, Iρi, na-1, nn, 1_000)

        if isnan(lr) || nn > 999
          return tree, na, nn, NaN
        end

        setproperty!(tree, :iμ, isextinct(stree))
        sete!(tree, e(tree) + e(stree))

        lλs = lλ(stree)
        lμs = lμ(stree)

        if lastindex(lλs) === 2
          setfdt!(tree, fdt(tree) + fdt(stree))
        else
          setfdt!(tree, fdt(stree))
        end

        pop!(lλ0)
        pop!(lμ0)
        popfirst!(lλs)
        popfirst!(lμs)
        append!(lλ0, lλs)
        append!(lμ0, lμs)

        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else

      tree.d1, na, nn, lr =
        tip_sims!(tree.d1, t, α, σλ, ix, tv, ev, δt, srδt, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr =
        tip_sims!(tree.d2, t, α, σλ, ix, tv, ev, δt, srδt, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_gbm!(bix     ::Int64,
                Ξ       ::Vector{iTbd},
                idf     ::Vector{iBffs},
                α       ::Float64,
                σλ      ::Float64,
                σμ      ::Float64,
                llc     ::Float64,
                prc     ::Float64,
                ddλ     ::Float64,
                ssλ     ::Float64,
                mc      ::Float64,
                th      ::Float64,
                δt      ::Float64,
                srδt    ::Float64,
                λ0_prior::NTuple{2,Float64}, 
                surv    ::Int64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix     ::Int64,
                     Ξ       ::Vector{iTbd},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     σμ      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ddλ     ::Float64,
                     ssλ     ::Float64,
                     mc      ::Float64,
                     th      ::Float64,
                     δt      ::Float64,
                     srδt    ::Float64,
                     λ0_prior::NTuple{2,Float64}, 
                     surv    ::Int64)

  # here
  @inbounds begin

    ξi   = Ξ[bix]
    bi   = idf[bix]
    i1   = d1(bi)
    i2   = d2(bi)
    ξ1   = Ξ[i1]
    root = iszero(pa(bi))

    # if crown
    if root && iszero(e(bi))
      ξ2 = Ξ[i2]
      llc, prc, ddλ, ssλ, mc =
        _crown_update!(ξi, ξ1, ξ2, α, σλ, σμ, llc, prc, ddλ, ssλ, mc, th,
          δt, srδt, λ0_prior, surv)
      setλt!(bi, lλ(ξi)[1])
    else
      # if stem
      if root
        llc, prc, ddλ, ssλ, mc =
          _stem_update!(ξi, α, σλ, σμ, llc, prc, ddλ, ssλ, mc, th, 
            δt, srδt, λ0_prior, surv)
      end

      # updates within the parent branch
      llc, ddλ, ssλ =
        _update_gbm!(ξi, α, σλ, llc, ddλ, ssλ, δt, srδt, false)

      # get fixed tip
      lξi = fixtip(ξi)

      if iszero(i2)

        llc, ssλ =
          update_duo!(lλ(lξi), lλ(ξ1), e(lξi), e(ξ1),
            fdt(lξi), fdt(ξ1), α, σλ, llc, ssλ, δt, srδt)

      # if internal branch
      else
        ξ2 = Ξ[i2]
        # make between decoupled trees node update
        llc, ddλ, ssλ =
          update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2),
            e(lξi), e(ξ1), e(ξ2), fdt(lξi), fdt(ξ1), fdt(ξ2),
            α, σλ, llc, ddλ, ssλ, δt, srδt)

        # set fixed `λ(t)` in branch
        setλt!(bi, lλ(ξ1)[1])
      end
    end

    # # carry on updates in the daughters
    llc, ddλ, ssλ =
      _update_gbm!(ξ1, α, σλ, llc, ddλ, ssλ, δt, srδt, iszero(d1(idf[i1])))
    if i2 > 0
      llc, ddλ, ssλ =
        _update_gbm!(Ξ[i2], α, σλ, llc, ddλ, ssλ, δt, srδt, iszero(d1(idf[i2])))
    end
  end

  return llc, prc, ddλ, ssλ, mc
end


