#=

Anagenetic GBM birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmce(tree    ::sT_label,
                 out_file::String;
                 λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                 α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                 μ_prior ::NTuple{2,Float64} = (1.0, 1.0),
                 niter   ::Int64             = 1_000,
                 nthin   ::Int64             = 10,
                 nburn   ::Int64             = 200,
                 marginal::Bool              = false,
                 nitpp   ::Int64             = 100,
                 nthpp   ::Int64             = 10,
                 K       ::Int64             = 11,
                 λi      ::Float64           = NaN,
                 αi      ::Float64           = 0.0,
                 σλi     ::Float64           = 0.01,
                 μi      ::Float64           = NaN,
                 ϵi      ::Float64           = 0.2,
                 pupdp   ::NTuple{5,Float64} = (0.1,0.1,0.1,0.2,0.2),
                 nlim    ::Int64             = 500,
                 δt      ::Float64           = 1e-2,
                 prints  ::Int64             = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `gbm-ce`.
"""
function insane_gbmce(tree    ::sT_label,
                      out_file::String;
                      λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                      α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                      μ_prior ::NTuple{2,Float64} = (1.0, 1.0),
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      marginal::Bool              = false,
                      nitpp   ::Int64             = 100,
                      nthpp   ::Int64             = 10,
                      K       ::Int64             = 11,
                      λi      ::Float64           = NaN,
                      αi      ::Float64           = 0.0,
                      σλi     ::Float64           = 0.01,
                      μi      ::Float64           = NaN,
                      ϵi      ::Float64           = 0.2,
                      pupdp   ::NTuple{5,Float64} = (0.1,0.1,0.1,0.2,0.2),
                      nlim    ::Int64             = 500,
                      δt      ::Float64           = 1e-2,
                      prints  ::Int64             = 5,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)
  stem = !iszero(e(tree))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ)

   # starting parameters (using method of moments)
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), ti(idf[1]), ϵi)
  else
    λc, μc = λi, μi
  end
  mc = m_surv_gbmce(th, log(λc), αi, σλi, μc, δt, srδt, 500, stem)

  # make a decoupled tree
  Ξ = make_Ξ(idf, log(λc), αi, σλi, δt, srδt, iTce)

  # set end of fix branch speciation times and
  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    setλt!(bi, lλ(Ξ[i])[end])
    if !it(bi)
      push!(inodes, i)
    end
  end

  # parameter updates (1: α, 2: σλ, 3: μ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running birth-death gbm with constant μ"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, μc, mc =
    mcmc_burn_gbmce(Ξ, idf, λa_prior, α_prior, σλ_prior, μ_prior,
      nburn, αi, σλi, μc, mc, th, stem, δt, srδt, inodes, pup, prints)

  # mcmc
  R, Ξv =
    mcmc_gbmce(Ξ, idf, llc, prc, αc, σλc, μc, mc, th, stem,
      λa_prior, α_prior, σλ_prior, μ_prior, niter, nthin, δt, srδt,
      inodes, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "alpha"        => 2,
                 "sigma_lambda" => 3,
                 "mu"           => 4))

  write_ssr(R, pardic, out_file)

  return R, Ξv
end




"""
    mcmc_burn_gbmce(Ξ       ::Vector{iTce},
                    idf     ::Vector{iBffs},
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    μ_prior ::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    μc      ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    stem    ::Bool,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Vector{Int64},
                    pup     ::Vector{Int64},
                    prints  ::Int64)

MCMC burn-in chain for `gbmce`.
"""
function mcmc_burn_gbmce(Ξ       ::Vector{iTce},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         μ_prior ::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         μc      ::Float64,
                         mc      ::Float64,
                         th      ::Float64,
                         stem    ::Bool,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Vector{Int64},
                         pup     ::Vector{Int64},
                         prints  ::Int64)

  λ0  = lλ(Ξ[1])[1]
  nsi = stem ? 0.0 : λ0

  llc = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) - nsi + log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2]) +
        logdunif(exp(λ0), λa_prior[1], λa_prior[2])   +
        logdnorm(αc,  α_prior[1], α_prior[2]^2)       +
        logdgamma(μc, μ_prior[1], μ_prior[2])

  # maximum bounds according to unfiorm priors
  lλxpr = log(λa_prior[2])

  L       = treelength(Ξ)      # tree length
  dλ      = deltaλ(Ξ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)     # sum squares in λ
  ne      = 0.0                # number of extinction events
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update drift
      if pupi === 1

        llc, prc, αc, mc =
          update_α!(αc, lλ(Ξ[1])[1], σλc, μc, L, dλ, llc, prc, mc, th, stem,
            δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

      # update sigma
      elseif pupi === 2

        llc, prc, σλc, mc =
          update_σ!(σλc, lλ(Ξ[1])[1], αc, μc, ssλ, nλ, llc, prc, mc, th, stem,
            δt, srδt, σλ_prior)

      # update extinction
      elseif pupi === 3

        llc, prc, μc, mc =
          update_μ!(μc, lλ(Ξ[1])[1], αc, σλc, llc, prc, ne, L, mc, th, stem,
            δt, srδt, μ_prior)

      # gbm update
      elseif pupi === 4

        # nix = ceil(Int64,rand()*nin)
        # bix = inodes[nix]
        bix = 1

        llc, dλ, ssλ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, μc, llc, dλ, ssλ, mc, th, stem,
            δt, srδt, lλxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, nλ, ne, L =
          update_fs!(bix, Ξ, idf, αc, σλc, μc, llc, dλ, ssλ, nλ, ne, L,
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
               stem    ::Bool,
               λa_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               μ_prior ::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Array{Int64,1},
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
                    stem    ::Bool,
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    μ_prior ::NTuple{2,Float64},
                    niter   ::Int64,
                    nthin   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # maximum bounds according to uniform priors
  lλxpr = log(λa_prior[2])

  L       = treelength(Ξ)            # tree length
  dλ      = deltaλ(Ξ)                # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)           # sum squares in λ
  ne      = Float64(ntipsextinct(Ξ)) # number of extinction events
  nin     = lastindex(inodes)        # number of internal nodes
  el      = lastindex(idf)           # number of branches

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 7)

  # make Ξ vector
  Ξv = iTce[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # check for extinct

      # update σλ or σμ
      if pupi === 1

        llc, prc, αc, mc =
          update_α!(αc, lλ(Ξ[1])[1], σλc, μc, L, dλ, llc, prc, mc, th, stem,
            δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) - !stem*lλ(Ξ[1])[1]  + log(mc) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      elseif pupi === 2

        llc, prc, σλc, mc =
          update_σ!(σλc, lλ(Ξ[1])[1], αc, μc, ssλ, nλ, llc, prc, mc, th, stem,
            δt, srδt, σλ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) - !stem*lλ(Ξ[1])[1]  + log(mc) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      elseif pupi === 3

         llc, prc, μc, mc =
            update_μ!(μc, lλ(Ξ[1])[1], αc, σλc, llc, prc, ne, L, mc, th, stem,
              δt, srδt, μ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) - !stem*lλ(Ξ[1])[1]  + log(mc) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      # gbm update
      elseif pupi === 4

        # nix = ceil(Int64,rand()*nin)
        # bix = inodes[nix]
        bix = 1

        llc, dλ, ssλ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, μc, llc, dλ, ssλ, mc, th, stem,
            δt, srδt, lλxpr)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) - !stem*lλ(Ξ[1])[1]  + log(mc) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, nλ, ne, L =
          update_fs!(bix, Ξ, idf, αc, σλc, μc, llc, dλ, ssλ, nλ, ne, L,
            δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, μc, δt, srδt) - !stem*lλ(Ξ[1])[1]  + log(mc) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-5)
        #    @show ll0, llc, pupi, i, Ξ
        #    return
        # end
      end
    end

    # log parameters
    lthin += 1
    if lthin === nthin
      lit += 1
      @inbounds begin
        R[lit,1] = Float64(lit)
        R[lit,2] = llc
        R[lit,3] = prc
        R[lit,4] = exp(lλ(Ξ[1])[1])
        R[lit,5] = αc
        R[lit,6] = σλc
        R[lit,7] = μc
        push!(Ξv, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ξv
end




"""
    update_fs!(bix    ::Int64,
               Ξ      ::Vector{iTce},
               idf    ::Vector{iBffs},
               α      ::Float64,
               σλ     ::Float64,
               μ      ::Float64,
               llc    ::Float64,
               dλ     ::Float64,
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
                    dλ     ::Float64,
                    ssλ    ::Float64,
                    nλ     ::Float64,
                    ne     ::Float64,
                    L      ::Float64,
                    δt     ::Float64,
                    srδt   ::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]

  # if terminal
  if it(bi)
    ξp, llr = fsbi_t(bi, lλ(ξc)[1], α, σλ, μ, δt, srδt)
    drλ  = 0.0
    ssrλ = 0.0
  # if internal
  else
    ξp, llr, drλ, ssrλ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], lλ(ξc)[1], α, σλ, μ, δt, srδt)
  end

  # if accepted
  if isfinite(llr)
    ll1, dλ1, ssλ1, nλ1 = llik_gbm_ssλ(ξp, α, σλ, μ, δt, srδt)
    ll0, dλ0, ssλ0, nλ0 = llik_gbm_ssλ(ξc, α, σλ, μ, δt, srδt)

    # update llr, ssλ, nλ, ne, L
    llc += llr  + ll1  - ll0
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    nλ  += nλ1  - nλ0
    ne  += ntipsextinct(ξp) - ntipsextinct(ξc)
    L   += treelength(ξp)   - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, dλ, ssλ, nλ, ne, L
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
                λ0  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                μ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_gbmce_t(e(bi), λ0, α, σλ, μ, δt, srδt, lc, lU, Iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr)
    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, -Inf
  end
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
                λ0  ::Float64,
                α   ::Float64,
                σλ  ::Float64,
                μ   ::Float64,
                δt  ::Float64,
                srδt::Float64)

  λsp = Float64[]

  t0, na, nn =
    _sim_gbmce_i(e(bi), λ0, α, σλ, μ, δt, srδt, 0, 1, 1_000, λsp)

  if na < 1 || nn >= 1_000
    return t0, NaN, NaN, NaN
  end

  # get current speciation rates at branch time
  λsc = λst(bi)

  e1  = e(ξ1)
  e2  = e(ξ2)
  λ1c = lλ(ξ1)
  λ2c = lλ(ξ2)
  l1  = lastindex(λ1c)
  l2  = lastindex(λ2c)
  λ1  = λ1c[l1]
  λ2  = λ2c[l2]

  # current acceptance ratio
  ac = 0.0
  for λi in λsc
    ac += exp(λi) * dnorm_bm(λi, λ1 - α*e1, sqrt(e1)*σλ) *
                    dnorm_bm(λi, λ2 - α*e2, sqrt(e2)*σλ)
  end
  ac = log(ac)

  # proposed acceptance ratio
  wp = Float64[]
  ap = 0.0
  for λi in λsp
    wi  = exp(λi) * dnorm_bm(λi, λ1 - α*e1, sqrt(e1)*σλ) *
                    dnorm_bm(λi, λ2 - α*e2, sqrt(e2)*σλ)
    ap += wi
    push!(wp, wi)
  end
  ap = log(ap)

  if isinf(ap)
    return t0, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = ap - ac

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample tip
  wti = sample(wp)
  λf  = λsp[wti]

  llrd, acrd, drλ, ssrλ, λ1p, λ2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, μ, δt, srδt)

  acr += acrd

  if lU < acr

     # fix sampled tip
    lw = lastindex(wp)

    if wti <= div(lw,2)
      fixtip1!(t0, wti, 0)
    else
      fixtip2!(t0, lw - wti + 1, 0)
    end

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr = 
        tip_sims!(t0, tf(bi), α, σλ, μ, δt, srδt, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1
      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setni!( bi, na)                       # set new ni
      setλt!( bi, λf)                       # set new λt
      setλst!(bi, λsp)                      # set new λst
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 2 λ vector

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
              Iρi ::Float64,
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
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

  if lU < lr && nn < 500

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)

        # simulate
        stree, na, nn, lr =
          _sim_gbmce_it(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, μ, δt, srδt,
                     lr, lU, Iρi, na-1, nn, 500)

        if isnan(lr) || nn >= 500
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
        tip_sims!(tree.d1, t, α, σλ, μ, δt, srδt, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, α, σλ, μ, δt, srδt, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTce},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                μ    ::Float64,
                llc  ::Float64,
                dλ   ::Float64,
                ssλ  ::Float64,
                mc   ::Float64,
                th   ::Float64,
                stem ::Bool,
                δt   ::Float64,
                srδt ::Float64,
                lλxpr::Float64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTce},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     μ    ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     mc   ::Float64,
                     th   ::Float64,
                     stem ::Bool,
                     δt   ::Float64,
                     srδt ::Float64,
                     lλxpr::Float64)

  @inbounds begin
    ξi   = Ξ[bix]
    bi   = idf[bix]
    ξ1   = Ξ[d1(bi)]
    ξ2   = Ξ[d2(bi)]
    ter1 = it(idf[d1(bi)])
    ter2 = it(idf[d2(bi)])


    root = iszero(pa(bi))
    # if crown
    if root && !stem
      llc, dλ, ssλ, mc =
        _crown_update!(ξi, ξ1, ξ2, α, σλ, μ, llc, dλ, ssλ, mc, th,
          δt, srδt, lλxpr)
      setλt!(bi, lλ(ξi)[1])
    else
      # if stem
      if root
        llc, dλ, ssλ, mc =
          _stem_update!(ξi, α, σλ, μ, llc, dλ, ssλ, mc, th, δt, srδt, lλxpr)
      end

      # parent branch update
      llc, dλ, ssλ =
        _update_gbm!(ξi, α, σλ, μ, llc, dλ, ssλ, δt, srδt, false)

      # get fixed tip
      lξi = fixtip(ξi)

      # make between decoupled trees node update
      llc, dλ, ssλ =
        update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2),
          fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, μ, llc, dλ, ssλ, δt, srδt)

      # set fixed `λ(t)` in branch
      setλt!(bi, lλ(lξi)[end])
    end

    # carry on updates in the daughters
    llc, dλ, ssλ = _update_gbm!(ξ1, α, σλ, μ, llc, dλ, ssλ, δt, srδt, ter1)
    llc, dλ, ssλ = _update_gbm!(ξ2, α, σλ, μ, llc, dλ, ssλ, δt, srδt, ter2)
  end

  return llc, dλ, ssλ, mc
end




"""
    update_α!(αc     ::Float64,
              λ0     ::Float64,
              σλ     ::Float64,
              μ      ::Float64,
              L      ::Float64,
              dλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              mc     ::Float64,
              th     ::Float64,
              stem   ::Bool,
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
                   dλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   stem   ::Bool,
                   δt     ::Float64,
                   srδt   ::Float64,
                   α_prior::NTuple{2,Float64})

  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2
  αp  = rnorm((dλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  mp  = m_surv_gbmce(th, λ0, αp, σλ, μ, δt, srδt, 500, stem)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/L) + llr
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
              stem    ::Bool,
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
                   stem    ::Bool,
                   δt      ::Float64,
                   srδt    ::Float64,
                   σλ_prior::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)
  σλp  = sqrt(σλp2)

  mp  = m_surv_gbmce(th, λ0, α, σλp, μ, δt, srδt, 500, stem)
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
              stem   ::Bool,
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
                   stem   ::Bool,
                   δt     ::Float64,
                   srδt   ::Float64,
                   μ_prior::NTuple{2,Float64})

  μp  = randgamma(μ_prior[1] + ne, μ_prior[2] + L)

  mp  = m_surv_gbmce(th, λ0, α, σλ, μp, δt, srδt, 500, stem)
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

