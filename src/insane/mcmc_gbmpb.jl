#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#





"""
    insane_gbmpb(tree    ::sT_label, 
                 out_file::String;
                 δt      ::Float64               = 1e-2,
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 marginal::Bool                  = false,
                 nitpp   ::Int64                 = 100, 
                 nthpp   ::Int64                 = 10,
                 K       ::Int64                 = 11,
                 σλi     ::Float64               = 0.1,
                 αi      ::Float64               = 0.0,
                 prints  ::Int64                 = 5,
                 pupdp   ::NTuple{4,Float64}     = (0.2, 0.2, 0.3, 0.3),
                 α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for GBM pure-birth.
"""
function insane_gbmpb(tree    ::sT_label, 
                      out_file::String;
                      δt      ::Float64               = 1e-2,
                      niter   ::Int64                 = 1_000,
                      nthin   ::Int64                 = 10,
                      nburn   ::Int64                 = 200,
                      marginal::Bool                  = false,
                      nitpp   ::Int64                 = 100, 
                      nthpp   ::Int64                 = 10,
                      K       ::Int64                 = 11,
                      σλi     ::Float64               = 0.1,
                      αi      ::Float64               = 0.0,
                      prints  ::Int64                 = 5,
                      pupdp   ::NTuple{4,Float64}     = (0.2, 0.2, 0.3, 0.3),
                      α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
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

  # lλ root node
  lλa = log(λmle_cpb(tree))

  # make fix tree directory
  idf = make_idf(tree, tρ)

  # make a decoupled tree
  Ξ = iTgbmpb[]
  iTgbmpb!(Ξ, tree, δt, srδt, lλa, αi, σλi)

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

  # parameter updates (1: α, 2: σ, 3: gbm, 4: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running pure-birth gbm"

  # burn-in phase
  llc, prc, αc, σλc =
    mcmc_burn_gbmpb(Ξ, idf, α_prior, σλ_prior, nburn, αi, σλi, 
      δt, srδt, inodes, pup, prints)

  # mcmc
  r, Ξv, αc, σλc = mcmc_gbmpb(Ξ, idf, llc, prc, αc, σλc, α_prior, σλ_prior, 
        niter, nthin, δt, srδt, inodes, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "alpha"        => 2,
                 "sigma_lambda" => 3))

  write_ssr(r, pardic, out_file)

  if marginal

     # reference distribution
    βs = [range(0.0, 1.0, K)...]
    reverse!(βs)

    # make reference Gaussian posterior for `α`
    @views p = r[:,5]
    m      = median(p)
    sd     = std(p)
    α_rdist = (m, sd)

    # make reference IGamma posterior for `σλ`
    @views p = r[:,6]
    m       = mean(p)
    v       = var(p)
    σλ_rdist = (m^2/v, m/v)

    # marginal likelihood
    pp = ref_posterior(Ξ, idf, llc, prc, αc, σλc,
      α_prior, σλ_prior, α_rdist, σλ_rdist,
      nitpp, nthpp, βs, δt, srδt, inodes, pup)

    # process with reference distribution the posterior
    p1 = Vector{Float64}(undef, size(r,1))
    for i in Base.OneTo(size(r,1))
      p1[i] = r[i,2] + r[i,3] - 
              logdnorm(    r[i,5],    α_rdist[1], α_rdist[2]^2) +
              logdinvgamma(r[i,6]^2, σλ_rdist[1], σλ_rdist[2])
    end
    pp[1] = p1

    reverse!(pp)
    reverse!(βs)

    ml = gss(pp, βs)
  else
    ml = NaN
  end

  return r, Ξv, ml
end



"""
    mcmc_burn_gbmpb(Ξ       ::Vector{iTgbmpb},
                    idf     ::Vector{iBffs},
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)


MCMC burn-in chain for GBM pure-birth.
"""
function mcmc_burn_gbmpb(Ξ       ::Vector{iTgbmpb},
                         idf     ::Vector{iBffs},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  # starting likelihood and prior
  llc = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2]) +
        logdnorm(αc, α_prior[1], α_prior[2]^2)

  L       = treelength(Ξ)      # tree length
  dλ      = deltaλ(Ξ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)     # sum squares in λ
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update drift
      if pupi === 1

        llc, prc, αc = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

      # update diffusion
      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior)

      # update gbm
      elseif pupi === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ = 
          update_gbm!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, δt, srδt)

      # forward simulation
      else
        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, nλ, L = 
          update_fs!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, nλ, L, δt, srδt)
      end

    end

    next!(pbar)
  end

  return llc, prc, αc, σλc
end




"""
    mcmc_gbmpb(Ξ       ::Vector{iTgbmpb},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               λ0_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Array{Int64,1},
               prints  ::Int64)

MCMC chain for GBM pure-birth.
"""
function mcmc_gbmpb(Ξ       ::Vector{iTgbmpb},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
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

  r = Array{Float64,2}(undef, nlogs, 6)

  # make Ξ vector
  Ξv = iTgbmpb[]

  L       = treelength(Ξ)      # tree length
  dλ      = deltaλ(Ξ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)     # sum squares in λ
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update drift
      if pupi === 1

        llc, prc, αc = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, it, pupi
        #    return 
        # end

      # update diffusion rate
      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, it, pupi
        #    return 
        # end

      # update gbm
      elseif pupi === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ = 
          update_gbm!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
        # if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, it, pupi
        #    return 
        # end

      # update by forward simulation
      else
        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, nλ, L = 
          update_fs!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, nλ, L, δt, srδt)

        # ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
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
        r[lit,1] = Float64(lit)
        r[lit,2] = llc
        r[lit,3] = prc
        r[lit,4] = exp(lλ(Ξ[1])[1])
        r[lit,5] = αc
        r[lit,6] = σλc
        push!(Ξv, couple(deepcopy(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return r, Ξv, αc, σλc
end





"""
    ref_posterior(Ξ       ::Vector{iTgbmpb},
                  idf     ::Vector{iBffs},
                  llc     ::Float64,
                  prc     ::Float64,
                  αc      ::Float64,
                  σλc     ::Float64,
                  λ0_prior::NTuple{2,Float64},
                  α_prior ::NTuple{2,Float64},
                  σλ_prior::NTuple{2,Float64},
                  λ0rdist::NTuple{2,Float64},
                  α_rdist ::NTuple{2,Float64},
                  σλ_rdist::NTuple{2,Float64},
                  nitpp   ::Int64,
                  nthpp   ::Int64,
                  βs      ::Vector{Float64},
                  δt      ::Float64,
                  srδt    ::Float64,
                  inodes  ::Array{Int64,1},
                  pup     ::Array{Int64,1},
                  prints  ::Int64)

MCMC chain for GBM pure-birth.
"""
function ref_posterior(Ξ       ::Vector{iTgbmpb},
                       idf     ::Vector{iBffs},
                       llc     ::Float64,
                       prc     ::Float64,
                       αc      ::Float64,
                       σλc     ::Float64,
                       α_prior ::NTuple{2,Float64},
                       σλ_prior::NTuple{2,Float64},
                       α_rdist ::NTuple{2,Float64},
                       σλ_rdist::NTuple{2,Float64},
                       nitpp   ::Int64,
                       nthpp   ::Int64,
                       βs      ::Vector{Float64},
                       δt      ::Float64,
                       srδt    ::Float64,
                       inodes  ::Array{Int64,1},
                       pup     ::Array{Int64,1})

  # starting likelihood and prior
  llc = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
  prc = logdnorm(αc,         α_prior[1], α_prior[2]^2) + 
        logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])

  K = lastindex(βs)

  # make log-likelihood table per power
  nlg = fld(nitpp, nthpp)
  pp  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

  L       = treelength(Ξ)      # tree length
  dλ      = deltaλ(Ξ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)     # sum squares in λ
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  for k in 2:K

    βi  = βs[k]
    rdc = logdnorm(       αc,  α_rdist[1], α_rdist[2]^2) +
          logdinvgamma(σλc^2, σλ_rdist[1], σλ_rdist[2])

    # logging
    lth, lit = 0, 0

    for it in Base.OneTo(nitpp)

      shuffle!(pup)

      for pupi in pup

        ## parameter updates
        # update drift
        if pupi === 1

          llc, prc, rdc, αc = update_α!(αc, σλc, L, dλ, llc, prc, rdc, 
            α_prior, α_rdist, βi)

          # update ssλ with new drift `α`
          ssλ, nλ = sss_gbm(Ξ, αc)

        # update diffusion rate
        elseif pupi === 2

          llc, prc, rdc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, rdc,
            σλ_prior, σλ_rdist, βi)

        # update gbm
        elseif pupi === 3

          nix = ceil(Int64,rand()*nin)
          bix = inodes[nix]

          llc, dλ, ssλ = 
            update_gbm!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, δt, srδt)

        # update by forward simulation
        else
          bix = ceil(Int64,rand()*el)

          llc, dλ, ssλ, nλ, L = 
            update_fs!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, nλ, L, δt, srδt)

        end
      end

      # log log-likelihood
      lth += 1
      if lth === nthpp
        lit += 1
        pp[k][lit] = llc + prc - rdc
        lth = 0
      end
    end

    @info string(βi," power done")
  end

  return pp
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{iTgbmpb},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               llc  ::Float64,
               dλ   ::Float64,
               ssλ  ::Float64,
               nλ   ::Float64,
               L    ::Float64,
               δt   ::Float64,
               srδt ::Float64)

Forward simulation proposal function for constant birth-death.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{iTgbmpb},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    llc  ::Float64,
                    dλ   ::Float64,
                    ssλ  ::Float64,
                    nλ   ::Float64,
                    L    ::Float64,
                    δt   ::Float64,
                    srδt ::Float64)

  bi  = idf[bix]
  ρbi = ρi(bi) # get branch sampling fraction
  nc  = ni(bi) # current ni
  ntc = nt(bi) # current nt
  itb = it(bi) # if terminal

  ξc  = Ξ[bix]
  if !itb
    ξ1  = Ξ[d1(bi)]
    ξ2  = Ξ[d2(bi)]
  end

  # forward simulate an internal branch
  ξp, np, ntp, λf = fsbi(bi, lλ(ξc)[1], α, σλ, δt, srδt)

  # check for non-exploding simulation
  if np >= 1000
    return llc, dλ, ssλ, nλ, L
  end

  # if terminal branch
  if itb
    llr  = log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc))
    acr  = llr
    drλ  = 0.0
    ssrλ = 0.0
  else
    np -= 1
    llr = log((1.0 - ρbi)^(np - nc))
    acr = llr + log(Float64(ntp)/Float64(ntc))
    # change daughters
    if isfinite(acr)

      llrd, acrd, drλ, ssrλ, λ1p, λ2p = 
        _daughters_update!(ξ1, ξ2, λf, α, σλ, δt, srδt)

      llr += llrd
      acr += acrd
    else
      acr = -Inf
    end
  end

  # MH ratio
  if -randexp() < acr

    ll1, dλ1, ssλ1, nλ1 = llik_gbm_ssλ(ξp, α, σλ, δt, srδt)
    ll0, dλ0, ssλ0, nλ0 = llik_gbm_ssλ(ξc, α, σλ, δt, srδt)

    # update llr, ssλ, nλ, L
    llr += ll1  - ll0
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    nλ  += nλ1  - nλ0
    L   += treelength(ξp) - treelength(ξc)

    Ξ[bix] = ξp          # set new tree
    llc += llr           # set new likelihood
    setni!(bi, np)       # set new ni
    setnt!(bi, ntp)      # set new nt
    setλt!(bi, λf)       # set new λt
    if !itb
      copyto!(lλ(ξ1), λ1p) # set new daughter 1 λ vector
      copyto!(lλ(ξ2), λ2p) # set new daughter 2 λ vector
    end
  end

  return llc, dλ, ssλ, nλ, L
end




"""
    fsbi(bi  ::iBffs, 
         λ0  ::Float64, 
         α   ::Float64, 
         σλ  ::Float64, 
         δt  ::Float64, 
         srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi(bi  ::iBffs,
              λ0  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64)

  # times
  tfb = tf(bi)

  # forward simulation during branch length
  t0, na = _sim_gbmpb(e(bi), λ0, α, σλ, δt, srδt, 1, 1_000)

  if na >= 1_000
    return iTgbmpb(), 1_000, 0, 0.0
  end

  nat = na

  if isone(na)
    f, λf = fixalive!(t0, NaN)

    return t0, na, nat, λf
  elseif na > 1
    # fix random tip
    λf = fixrtip!(t0, na, NaN)

    if !it(bi)
      # add tips until the present
      tx, na = tip_sims!(t0, tfb, α, σλ, δt, srδt, na)
    end

    return t0, na, nat, λf
  end

  return iTgbmpb(), 0, 0, 0.0
end




"""
    tip_sims!(tree::iTgbmpb, 
              t   ::Float64, 
              α   ::Float64, 
              σλ  ::Float64, 
              δt  ::Float64,
              srδt::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::iTgbmpb, 
                   t   ::Float64, 
                   α   ::Float64, 
                   σλ  ::Float64, 
                   δt  ::Float64,
                   srδt::Float64,
                   na  ::Int64)

  if istip(tree) 
    if !isfix(tree)

      fdti = fdt(tree)
      lλ0  = lλ(tree)

      # simulate
      stree, na = 
        _sim_gbmpb(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, δt, srδt, na, 1_000)

      if !isdefined(stree, :lλ)
        return tree, 1_000
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
    tree.d1, na = tip_sims!(tree.d1, t, α, σλ, δt, srδt, na)
    tree.d2, na = tip_sims!(tree.d2, t, α, σλ, δt, srδt, na)
  end

  return tree, na
end




"""
    update_gbm!(bix  ::Int64,
                xi  ::Vector{iTgbmpb},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                llc  ::Float64,
                dλ   ::Float64,
                ssλ  ::Float64,
                δt   ::Float64,
                srδt ::Float64)

Make a `gbm` update for an interna branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTgbmpb},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     δt   ::Float64,
                     srδt ::Float64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  ξ1   = Ξ[d1(bi)]
  ξ2   = Ξ[d2(bi)]
  ter1 = it(idf[d1(bi)]) 
  ter2 = it(idf[d2(bi)])

  # if crown root
  if iszero(pa(bi)) && iszero(e(ξi))
    llc, dλ, ssλ = 
      _crown_update!(ξi, ξ1, ξ2, α, σλ, llc, dλ, ssλ, δt, srδt)
    setλt!(bi, lλ(ξi)[1])
  else
    # if stem
    if iszero(pa(bi))
      llc, dλ, ssλ = _stem_update!(ξi, α, σλ, llc, dλ, ssλ, δt, srδt)
    end

    # updates within the parent branch
    llc, dλ, ssλ = _update_gbm!(ξi, α, σλ, llc, dλ, ssλ, δt, srδt, false)

    # get fixed tip 
    lξi = fixtip(ξi) 

    # make between decoupled trees node update
    llc, dλ, ssλ = update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2), 
      fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, llc, dλ, ssλ, δt, srδt)

    # set fixed `λ(t)` in branch
    setλt!(bi, lλ(lξi)[end])
  end

  # carry on updates in the daughters
  llc, dλ, ssλ = _update_gbm!(ξ1, α, σλ, llc, dλ, ssλ, δt, srδt, ter1)
  llc, dλ, ssλ = _update_gbm!(ξ2, α, σλ, llc, dλ, ssλ, δt, srδt, ter2)

  return llc, dλ, ssλ
end





"""
    update_α!(αc     ::Float64,
              σλ     ::Float64,
              L      ::Float64,
              dλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              α_prior::NTuple{2,Float64})

Gibbs update for `α`.
"""
function update_α!(αc     ::Float64,
                   σλ     ::Float64,
                   L      ::Float64,
                   dλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   α_prior::NTuple{2,Float64})

  # ratio 
  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2

  # gibbs update for σ
  αp = rnorm((dλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  # update prior
  prc += llrdnorm_x(αp, αc, ν, τ2)

  # update likelihood
  llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/L)

  return llc, prc, αp
end




"""
    update_α!(αc     ::Float64,
              σλ     ::Float64,
              L      ::Float64,
              dλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              rdc    ::Float64,
              α_prior::NTuple{2,Float64},
              α_rdist ::NTuple{2,Float64},
              pow    ::Float64)


Gibbs update for `α` given reference distribution.
"""
function update_α!(αc     ::Float64,
                   σλ     ::Float64,
                   L      ::Float64,
                   dλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   rdc    ::Float64,
                   α_prior::NTuple{2,Float64},
                   α_rdist ::NTuple{2,Float64},
                   pow    ::Float64)

  # ratio 
  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2

  cpow = (1.0 - pow)

  # gibbs update for α
  m   = (dλ + rs*ν)/(rs + L)
  s2  = σλ2/(rs + L)
  m0  = α_rdist[1]
  s02 = α_rdist[2]^2
  αp  = rnorm((m0 * s2 * cpow + m * s02 * pow) / (pow * s02 + s2 * cpow),
              sqrt( s2 * s02 /  (pow * s02 + s2 * cpow)) )

  # update likelihood, prior and reference
  llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/L)
  prc += llrdnorm_x(αp, αc, ν, τ2)
  rdc += llrdnorm_x(αp, αc, m0, s02)

  return llc, prc, rdc, αp
end




"""
    update_σ!(σλc     ::Float64,
              α       ::Float64,
              ssλ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64}) 

Gibbs update for `σλ`.
"""
function update_σ!(σλc     ::Float64,
                   α       ::Float64,
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
    update_σ!(σλc     ::Float64,
              α       ::Float64,
              ssλ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64}) 

Gibbs update for `σλ` given reference distribution.
"""
function update_σ!(σλc     ::Float64,
                   α       ::Float64,
                   ssλ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   rdc     ::Float64,
                   σλ_prior::NTuple{2,Float64},
                   σλ_rdist ::NTuple{2,Float64},
                   pow     ::Float64)

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma((σλ_p1 + 0.5 * n) * pow + σλ_rdist[1] * (1.0 - pow), 
                      (σλ_p2 + ssλ) * pow     + σλ_rdist[2] * (1.0 - pow)) 

  # update likelihood, prior and reference
  σλp = sqrt(σλp2)
  llc += ssλ*(1.0/σλc^2 - 1.0/σλp2) - n*(log(σλp/σλc))
  prc += llrdinvgamma(σλp2, σλc^2, σλ_p1, σλ_p2)
  rdc += llrdinvgamma(σλp2, σλc^2, σλ_rdist[1], σλ_rdist[2])

  return llc, prc, rdc, σλp
end


