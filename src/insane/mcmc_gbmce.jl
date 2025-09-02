#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmce(tree    ::sT_label;
                 λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                 α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                 σλ_prior::NTuple{2,Float64}     = (3.0, 0.05),
                 μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 nflush  ::Int64                 = nthin,
                 ofile   ::String                = string(homedir(), "/ice"),
                 λi      ::Float64               = NaN,
                 αi      ::Float64               = 0.0,
                 σλi     ::Float64               = 0.01,
                 μi      ::Float64               = NaN,
                 ϵi      ::Float64               = 0.2,
                 pupdp   ::NTuple{5,Float64}     = (1e-3, 1e-3, 1e-4, 0.1, 0.2),
                 δt      ::Float64               = 1e-3,
                 prints  ::Int64                 = 5,
                 survival::Bool                  = true,
                 mxthf   ::Float64               = 0.1,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `gbm-ce`.
"""
function insane_gbmce(tree    ::sT_label;
                      λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                      α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                      σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                      μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                      niter   ::Int64                 = 1_000,
                      nthin   ::Int64                 = 10,
                      nburn   ::Int64                 = 200,
                      nflush  ::Int64                 = nthin,
                      ofile   ::String                = string(homedir(), "/ice"),
                      λi      ::Float64               = NaN,
                      αi      ::Float64               = 0.0,
                      σλi     ::Float64               = 0.01,
                      μi      ::Float64               = NaN,
                      ϵi      ::Float64               = 0.2,
                      pupdp   ::NTuple{5,Float64}     = (1e-3, 1e-3, 1e-4, 0.1, 0.2),
                      δt      ::Float64               = 1e-3,
                      prints  ::Int64                 = 5,
                      survival::Bool                  = true,
                      mxthf   ::Float64               = 0.1,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # turn to logarithmic terms
  λ0_prior = (log(λ0_prior[1]), 2*log(λ0_prior[2]))

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
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # estimate branch split (multiple of δt)
  ndts = floor(th * mxthf/δt)
  maxt = δt * ndts

  # make fix tree directory
  idf = make_idf(tree, tρ, maxt)

   # starting parameters (using method of moments)
  λc, μc = λi, μi
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  end

  # make a decoupled tree
  Ξ = make_Ξ(idf, λc, αi, σλi, δt, srδt, iTce)

  #survival
  mc = m_surv_gbmce(th, log(λc), αi, σλi, μc, δt, srδt, 1_000, surv)

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf))  if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σλ, 3: μ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death gbm with constant μ"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, μc, mc =
    mcmc_burn_gbmce(Ξ, idf, λ0_prior, α_prior, σλ_prior, μ_prior,
      nburn, αi, σλi, μc, mc, th, rmλ, surv, δt, srδt, inodes, pup, prints)

  # mcmc
  r, treev =
    mcmc_gbmce(Ξ, idf, llc, prc, αc, σλc, μc, mc, th, surv,
      λ0_prior, α_prior, σλ_prior, μ_prior, δt, srδt, inodes, pup, 
      niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_gbmce(Ξ       ::Vector{iTce},
                    idf     ::Vector{iBffs},
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    μ_prior ::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    μc      ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    rmλ     ::Float64,
                    surv    ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Vector{Int64},
                    pup     ::Vector{Int64},
                    prints  ::Int64)

MCMC burn-in chain for `ice`.
"""
function mcmc_burn_gbmce(Ξ       ::Vector{iTce},
                         idf     ::Vector{iBffs},
                         λ0_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         μ_prior ::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         μc      ::Float64,
                         mc      ::Float64,
                         th      ::Float64,
                         rmλ     ::Float64,
                         surv    ::Int64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Vector{Int64},
                         pup     ::Vector{Int64},
                         prints  ::Int64)

  lλ0 = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) - rmλ * lλ0 + 
        log(mc) + prob_ρ(idf)
  prc = logdnorm(lλ0,       λ0_prior[1], λ0_prior[2])   +
        logdnorm(αc,         α_prior[1],  α_prior[2]^2) +
        logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])   +
        logdgamma(μc,        μ_prior[1],  μ_prior[2])

  L       = treelength(Ξ)      # tree length
  ne      = 0.0                # number of extinction events
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, nλ = _ss_dd(Ξ, αc)

  pbar = Progress(nburn, dt = prints, desc = "burning mcmc...", barlen = 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update drift
      if pupi === 1

        llc, prc, αc, mc =
          update_α!(αc, lλ(Ξ[1])[1], σλc, μc, L, ddλ, llc, prc, mc, th, surv,
            δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ = _ss(Ξ, lλ, αc)

      # update sigma
      elseif pupi === 2

        llc, prc, σλc, mc =
          update_σ!(σλc, lλ(Ξ[1])[1], αc, μc, ssλ, nλ, llc, prc, mc, th, surv,
            δt, srδt, σλ_prior)

      # update extinction
      elseif pupi === 3

        llc, prc, μc, mc =
          update_μ!(μc, lλ(Ξ[1])[1], αc, σλc, llc, prc, ne, L, mc, th, surv,
            δt, srδt, μ_prior)

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ddλ, ssλ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, μc, llc, prc, ddλ, ssλ, mc, th,
            δt, srδt, λ0_prior, surv)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, ddλ, ssλ, nλ, ne, L =
          update_fs!(bix, Ξ, idf, αc, σλc, μc, llc, ddλ, ssλ, nλ, ne, L, 
            δt, srδt)
      end
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, μc, mc
end




"""
    mcmc_gbmce(Ξ       ::Vector{iTce},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               μc      ::Float64,
               mc      ::Float64,
               th      ::Float64,
               surv    ::Int64,
               λ0_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               μ_prior ::NTuple{2,Float64},
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Vector{Int64},
               niter   ::Int64,
               nthin   ::Int64,
               nflush  ::Int64,
               ofile   ::String,
               prints  ::Int64)

MCMC chain for `gbmce`.
"""
function mcmc_gbmce(Ξ       ::Vector{iTce},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    μc      ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    surv    ::Int64,
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    μ_prior ::NTuple{2,Float64},
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

  L       = treelength(Ξ)            # tree length
  ne      = Float64(ntipsextinct(Ξ)) # number of extinction events
  nin     = lastindex(inodes)        # number of internal nodes
  el      = lastindex(idf)           # number of branches

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, nλ = _ss_dd(Ξ, αc)

  # parameter results
  r = Array{Float64,2}(undef, nlogs, 7)

  treev = iTce[]     # make tree vector
  io    = IOBuffer() # buffer 

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\talpha\tsigma_lambda\tmu\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αc = αc, σλc = σλc, μc = μc, mc = mc, nλ = nλ, ssλ = ssλ, ddλ = ddλ, L = L, ne = ne, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          # parameter updates
          for pupi in pup

            # update α
            if pupi === 1

              llc, prc, αc, mc =
                update_α!(αc, lλ(Ξ[1])[1], σλc, μc, L, ddλ, llc, prc, mc, th, surv,
                  δt, srδt, α_prior)

              # update ssλ with new drift `α`
              ssλ = _ss(Ξ, lλ, αc)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) -Float64(surv > 0) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-5)
              #    @show ll0, llc, pupi, i
              #    return
              # end

            # update σλ
            elseif pupi === 2

              llc, prc, σλc, mc =
                update_σ!(σλc, lλ(Ξ[1])[1], αc, μc, ssλ, nλ, llc, prc, mc, th, surv,
                  δt, srδt, σλ_prior)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) -Float64(surv > 0) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-5)
              #    @show ll0, llc, pupi, i
              #    return
              # end

            elseif pupi === 3

              llc, prc, μc, mc =
                update_μ!(μc, lλ(Ξ[1])[1], αc, σλc, llc, prc, ne, L, mc, th, surv,
                  δt, srδt, μ_prior)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) -Float64(surv > 0) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-5)
              #    @show ll0, llc, pupi, i
              #    return
              # end

            # gbm update
            elseif pupi === 4

              nix = ceil(Int64,rand()*nin)
              bix = inodes[nix]

              llc, prc, ddλ, ssλ, mc =
                update_gbm!(bix, Ξ, idf, αc, σλc, μc, llc, prc, ddλ, ssλ, mc, th,
                  δt, srδt, λ0_prior, surv)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) -Float64(surv > 0) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-5)
              #    @show ll0, llc, pupi, i
              #    return
              # end

            # forward simulation update
            else

              bix = ceil(Int64,rand()*el)

              llc, ddλ, ssλ, nλ, ne, L =
                update_fs!(bix, Ξ, idf, αc, σλc, μc, llc, ddλ, ssλ, nλ, ne, L,
                  δt, srδt)

              # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) -Float64(surv > 0) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-5)
              #    @show ll0, llc, pupi, i
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
              r[lit,7] = μc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
                  exp(lλ(Ξ[1])[1]),'\t',  αc, '\t', σλc, '\t', μc,'\n')
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
    update_fs!(bix    ::Int64,
               Ξ      ::Vector{iTce},
               idf    ::Vector{iBffs},
               α      ::Float64,
               σλ     ::Float64,
               μ      ::Float64,
               llc    ::Float64,
               ddλ     ::Float64,
               ssλ    ::Float64,
               nλ     ::Float64,
               ne     ::Float64,
               L      ::Float64,
               δt     ::Float64,
               srδt   ::Float64)

Forward simulation proposal function for `gbmce`.
"""
function update_fs!(bix    ::Int64,
                    Ξ      ::Vector{iTce},
                    idf    ::Vector{iBffs},
                    α      ::Float64,
                    σλ     ::Float64,
                    μ      ::Float64,
                    llc    ::Float64,
                    ddλ     ::Float64,
                    ssλ    ::Float64,
                    nλ     ::Float64,
                    ne     ::Float64,
                    L      ::Float64,
                    δt     ::Float64,
                    srδt   ::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]

  # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, ξc, α, σλ, μ, δt, srδt)
    drλ = ssrλ = 0.0

  # if mid
  elseif iszero(d2(bi))
    ξp, llr, drλ, ssrλ =
      fsbi_m(bi, ξc, Ξ[d1(bi)], α, σλ, μ, δt, srδt)

  # if internal
  else
    ξp, llr, drλ, ssrλ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, μ, δt, srδt)
  end

  # if accepted
  if isfinite(llr)
    ll1, ddλ1, ssλ1, nλ1 = llik_gbm_ssλ(ξp, α, σλ, μ, δt, srδt)
    ll0, ddλ0, ssλ0, nλ0 = llik_gbm_ssλ(ξc, α, σλ, μ, δt, srδt)

    # update llr, ssλ, nλ, ne, L
    llc += llr  + ll1  - ll0
    ddλ  += ddλ1  - ddλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    nλ  += nλ1  - nλ0
    ne  += ntipsextinct(ξp) - ntipsextinct(ξc)
    L   += treelength(ξp)   - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ssλ, nλ, ne, L
end




"""
    fsbi_t(bi  ::iBffs,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           μ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTce,
                α   ::Float64,
                σλ  ::Float64,
                μ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_gbmce_t(e(bi), lλ(ξc)[1], α, σλ, μ, δt, srδt, lc, lU, iρi, 0, 1, 1_000)

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
           ξc  ::iTce,
           ξ1  ::iTce,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           μ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_m(bi  ::iBffs,
                ξc  ::iTce,
                ξ1  ::iTce,
                α   ::Float64,
                σλ  ::Float64,
                μ   ::Float64,
                δt  ::Float64,
                srδt::Float64)


  t0, na, nn =
    _sim_gbmce(e(bi), lλ(ξc)[1], α, σλ, μ, δt, srδt, 0, 1, 1_000)

  if na < 1 || nn > 999
    return t0, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  λf = fixrtip!(t0, na, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, λ1p =
    _daughter_update!(ξ1, λf, α, σλ, μ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr = 
        tip_sims!(t0, tf(bi), α, σλ, μ, δt, srδt, acr, lU, iρi, na, nn)
    end

    if lU < acr
      na -= 1
      llr = llrd + (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
      setnt!(bi, ntp)                                   # set new nt
      setni!(bi, na)                                    # set new ni
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, lastindex(λ1p)) # set new daughter 1 λ vector

      return t0, llr, drλ, ssrλ
    end
  end

  return t0, NaN, NaN, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTce,
           ξ1  ::iTce,
           ξ2  ::iTce,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           μ   ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTce,
                ξ1  ::iTce,
                ξ2  ::iTce,
                α   ::Float64,
                σλ  ::Float64,
                μ   ::Float64,
                δt  ::Float64,
                srδt::Float64)


  t0, na, nn =
    _sim_gbmce(e(bi), lλ(ξc)[1], α, σλ, μ, δt, srδt, 0, 1, 1_000)

  if na < 1 || nn > 999
    return t0, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  λf = fixrtip!(t0, na, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, λ1p, λ2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, μ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr = 
        tip_sims!(t0, tf(bi), α, σλ, μ, δt, srδt, acr, lU, iρi, na, nn)
    end

    if lU < acr
      na -= 1
      llr = llrd + (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
      setnt!(bi, ntp)                                   # set new nt
      setni!(bi, na)                                    # set new ni
      setλt!(bi, λf)                                    # set new λt
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, lastindex(λ1p)) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, lastindex(λ2p)) # set new daughter 2 λ vector

      return t0, llr, drλ, ssrλ
    end
  end

  return t0, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTce,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTce,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

  if lU < lr && nn < 1_000

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)

        # simulate
        stree, na, nn, lr =
          _sim_gbmce_it(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, μ, δt, srδt,
                     lr, lU, iρi, na-1, nn, 1_000)

        if isnan(lr) || nn > 999
          return tree, na, nn, NaN
        end

        setproperty!(tree, :iμ, isextinct(stree))
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
      tree.d1, na, nn, lr = 
        tip_sims!(tree.d1, t, α, σλ, μ, δt, srδt, lr, lU, iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, α, σλ, μ, δt, srδt, lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_gbm!(bix     ::Int64,
                Ξ       ::Vector{iTce},
                idf     ::Vector{iBffs},
                α       ::Float64,
                σλ      ::Float64,
                μ       ::Float64,
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
                     Ξ       ::Vector{iTce},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     μ       ::Float64,
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
        _crown_update!(ξi, ξ1, ξ2, α, σλ, μ, llc, prc, ddλ, ssλ, mc, th, 
          δt, srδt, λ0_prior, surv)
      setλt!(bi, lλ(ξi)[1])
    else
      # if stem
      if root
        llc, prc, ddλ, ssλ, mc =
          _stem_update!(ξi, α, σλ, μ, llc, prc, ddλ, ssλ, mc, th, 
            δt, srδt, λ0_prior, surv)
      end

      # parent branch update
      llc, ddλ, ssλ = _update_gbm!(ξi, α, σλ, μ, llc, ddλ, ssλ, δt, srδt, false)

      # get fixed tip
      lξi = fixtip(ξi)

      if iszero(i2)

        llc, ssλ =
          update_duo!(lλ(lξi), lλ(ξ1), e(lξi), e(ξ1),
            fdt(lξi), fdt(ξ1), α, σλ, μ, llc, ssλ, δt, srδt)

      else
        ξ2 = Ξ[i2]
        # make between decoupled trees node update
        llc, ddλ, ssλ =
          update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2),
            fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, μ, llc, ddλ, ssλ, δt, srδt)

        # set fixed `λ(t)` in branch
        setλt!(bi, lλ(ξ1)[1])
      end
    end

    # carry on updates in the daughters
    llc, ddλ, ssλ = 
      _update_gbm!(ξ1, α, σλ, μ, llc, ddλ, ssλ, δt, srδt, iszero(d1(idf[i1])))

    if i2 > 0
      ξ2 = Ξ[i2]
      llc, ddλ, ssλ = 
        _update_gbm!(ξ2, α, σλ, μ, llc, ddλ, ssλ, δt, srδt, iszero(d1(idf[i2])))
    end
  end

  return llc, prc, ddλ, ssλ, mc
end




"""
    update_α!(αc     ::Float64,
              λ0     ::Float64,
              σλ     ::Float64,
              μ      ::Float64,
              L      ::Float64,
              ddλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              mc     ::Float64,
              th     ::Float64,
              crown  ::Int64,
              δt     ::Float64,
              srδt   ::Float64,
              α_prior::NTuple{2,Float64})

Gibbs update for `α`.
"""
function update_α!(αc     ::Float64,
                   λ0     ::Float64,
                   σλ     ::Float64,
                   μ      ::Float64,
                   L      ::Float64,
                   ddλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   surv  ::Int64,
                   δt     ::Float64,
                   srδt   ::Float64,
                   α_prior::NTuple{2,Float64})

  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2
  αp  = rnorm((ddλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  mp  = m_surv_gbmce(th, λ0, αp, σλ, μ, δt, srδt, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*ddλ*(αp - αc)/L) + llr
    prc += llrdnorm_x(αp, αc, ν, τ2)
    αc   = αp
    mc   = mp
  end

  return llc, prc, αc, mc
end




"""

    update_σ!(σλc     ::Float64,
              λ0      ::Float64,
              α       ::Float64,
              μ       ::Float64,
              ssλ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              mc      ::Float64,
              th      ::Float64,
              crown   ::Int64,
              δt      ::Float64,
              srδt    ::Float64,
              σλ_prior::NTuple{2,Float64})

Gibbs update for `σλ`.
"""
function update_σ!(σλc     ::Float64,
                   λ0      ::Float64,
                   α       ::Float64,
                   μ       ::Float64,
                   ssλ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   surv   ::Int64,
                   δt      ::Float64,
                   srδt    ::Float64,
                   σλ_prior::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)
  σλp  = sqrt(σλp2)

  mp  = m_surv_gbmce(th, λ0, α, σλp, μ, δt, srδt, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += ssλ*(1.0/σλc^2 - 1.0/σλp2) - n*(log(σλp/σλc)) + llr
    prc += llrdinvgamma(σλp2, σλc^2, σλ_p1, σλ_p2)
    σλc  = σλp
    mc   = mp
  end

  return llc, prc, σλc, mc
end




"""
    update_μ!(μc     ::Float64,
              λ0     ::Float64,
              α      ::Float64,
              σλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              ne     ::Float64,
              L      ::Float64,
              mc     ::Float64,
              th     ::Float64,
              surv  ::Int64,
              δt     ::Float64,
              srδt   ::Float64,
              μ_prior::NTuple{2,Float64})

Gibbs-MH update for `μ`.
"""
function update_μ!(μc     ::Float64,
                   λ0     ::Float64,
                   α      ::Float64,
                   σλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   ne     ::Float64,
                   L      ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   surv  ::Int64,
                   δt     ::Float64,
                   srδt   ::Float64,
                   μ_prior::NTuple{2,Float64})

  μp  = randgamma(μ_prior[1] + ne, μ_prior[2] + L)

  mp  = m_surv_gbmce(th, λ0, α, σλ, μp, δt, srδt, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += ne * log(μp/μc) + L * (μc - μp) + llr
    prc += llrdgamma(μp, μc, μ_prior[1], μ_prior[2])
    μc   = μp
    mc   = mp
  end

  return llc, prc, μc, mc
end




# """
#     update_μ!(psi   ::Vector{iTce},
#               llc   ::Float64,
#               prc   ::Float64,
#               rdc   ::Float64,
#               μc    ::Float64,
#               μtn   ::Float64,
#               ne    ::Float64,
#               L     ::Float64,
#               sns   ::NTuple{3,BitVector},
#               μ_prior::Float64,
#               μ_refd ::NTuple{2,Float64},
#               scond ::Function,
#               pow   ::Float64)

# MCMC update for `μ`.
# """
# function update_μ!(psi   ::Vector{iTce},
#                    llc   ::Float64,
#                    prc   ::Float64,
#                    rdc   ::Float64,
#                    μc    ::Float64,
#                    μtn   ::Float64,
#                    ne    ::Float64,
#                    L     ::Float64,
#                    sns   ::NTuple{3,BitVector},
#                    μ_prior::NTuple{2,Float64},
#                    μ_refd ::NTuple{2,Float64},
#                    scond ::Function,
#                    pow   ::Float64)

#   # parameter proposal
#   μp = mulupt(μc, μtn)::Float64

#   # log likelihood and prior ratio
#   μr   = log(μp/μc)
#   llr  = ne*μr + L*(μc - μp) + scond(psi, μp, sns) - scond(psi, μc, sns)
#   prr = llrdgamma(μp, μc, μ_prior[1], μ_prior[2])
#   rdr = llrdtnorm(μp, μc, μ_refd[1],  μ_refd[2])


#   if -randexp() < (pow * (llr + prr) + (1.0 - pow) * rdr + μr)
#     llc += llr
#     prc += prr
#     rdc += rdr
#     μc   = μp
#   end

#   return llc, prc, rdc, μc
# end

