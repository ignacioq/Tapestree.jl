#=

Anagenetic `gbmbd` MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmbd(tree    ::sT_label, 
                 out_file::String;
                 λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                 μa_prior::NTuple{2,Float64} = (0.0, 100.0),
                 α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                 σμ_prior::NTuple{2,Float64} = (0.05, 0.05),
                 niter   ::Int64             = 1_000,
                 nthin   ::Int64             = 10,
                 nburn   ::Int64             = 200,
                 ϵi      ::Float64           = 0.2,
                 λi      ::Float64           = NaN,
                 μi      ::Float64           = NaN,
                 αi      ::Float64           = 0.0,
                 σλi     ::Float64           = 0.01, 
                 σμi     ::Float64           = 0.01,
                 pupdp   ::NTuple{4,Float64} = (0.1, 0.1, 0.2, 0.2),
                 δt      ::Float64           = 1e-2,
                 prints  ::Int64             = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `gbm-bd`.
"""
function insane_gbmbd(tree    ::sT_label, 
                      out_file::String;
                      λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                      μa_prior::NTuple{2,Float64} = (0.0, 100.0),
                      α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64} = (0.05, 0.05),
                      σμ_prior::NTuple{2,Float64} = (0.05, 0.05),
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      ϵi      ::Float64           = 0.2,
                      λi      ::Float64           = NaN,
                      μi      ::Float64           = NaN,
                      αi      ::Float64           = 0.0,
                      σλi     ::Float64           = 0.01, 
                      σμi     ::Float64           = 0.01,
                      pupdp   ::NTuple{4,Float64} = (0.0, 0.1, 0.2, 0.2),
                      δt      ::Float64           = 1e-2,
                      prints  ::Int64             = 5,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1,round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

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
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end

  # make a decoupled tree
  Ψ = iTgbmbd[]
  iTgbmbd!(Ψ, tree, δt, srδt, log(λc), log(μc), αi, σλi, σμi)

 # set end of fix branch speciation times and
  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    setλt!(bi, lλ(Ψ[i])[end])
    if !it(bi)
      push!(inodes, i)
    end
  end

  # parameter updates (1: α, 2: σλ, 3: σμ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # conditioning functions
  sns = (BitVector(), BitVector(), BitVector())
  snodes! = make_snodes(idf, !iszero(e(tree)), iTgbmbd)
  snodes!(Ψ, sns)
  scond, scond0 = make_scond(idf, !iszero(e(tree)), iTgbmbd)

  @info "running birth-death gbm"




  # burn-in phase
  Ψp, Ψc, llc, prc, αc, σλc, σμc =
    mcmc_burn_gbmbd(Ψp, Ψc, λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, 
      nburn, αi, σλi, σμi, δt, srδt, idf, pup, nlim, prints, svf)

  # mcmc
  R, Ψv =
    mcmc_gbmbd(Ψp, Ψc, llc, prc, αc, σλc, σμc,
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, niter, nthin, δt, srδt, 
      idf, pup, nlim, prints, svf)

  pardic = Dict(("lambda_root"  => 1,
                 "mu_root"      => 2,
                 "alpha"        => 3,
                 "sigma_lambda" => 4,
                 "sigma_mu"     => 5,
                 "n_extinct"    => 6))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmbd(Ψp      ::iTgbmbd,
                    Ψc      ::iTgbmbd,
                    λa_prior::NTuple{2,Float64},
                    μa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc     ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idf     ::Array{iBffs,1},
                    pup     ::Array{Int64,1},
                    nlim    ::Int64,
                    prints  ::Int64,
                    svf     ::Function)

MCMC burn-in chain for `gbmbd`.
"""
function mcmc_burn_gbmbd(Ψ       ::Vector{iTgbmbd},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         μa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         σμ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc     ::Float64,
                         σλc     ::Float64,
                         σμc     ::Float64,
                         sns     ::NTuple{3,BitVector},
                         δt      ::Float64,
                         srδt    ::Float64,
                         pup     ::Array{Int64,1},
                         prints  ::Int64,
                         snodes! ::Function,
                         scond   ::Function,
                         scond0  ::Function)

  llc = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) + scond(Ψ, sns) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])      + 
        logdinvgamma(σμc^2, σμ_prior[1], σμ_prior[2])      + 
        logdnorm(αc, α_prior[1], α_prior[2]^2)             +
        logdunif(exp(lλ(Ψ[1])[1]), λa_prior[1], λa_prior[2]) +
        logdunif(exp(lμ(Ψ[1])[1]), μa_prior[1], μa_prior[2])

  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L            = treelength(Ψ)      # tree length
  dλ           = deltaλ(Ψ)         # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ψ, αc)    # sum squares in λ and μ
  nin          = lastindex(inodes) # number of internal nodes
  el           = lastindex(idf)    # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update α
      if pupi === 1

        llc, prc, αc  = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nλ = sss_gbm(Ψ, αc)

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc = 
          update_σ!(σλc, σμc, αc, ssλ, ssμ, nλ, llc, prc, σλ_prior, σμ_prior)

      # gbm update
      elseif pupi === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, ssμ = 
          update_gbm!(bix, Ψ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, sns, δt, 
            srδt, lλxpr, lμxpr)

      # forward simulation update
      else

        bi  = idf[ceil(Int64,rand()*nbr)]

        if iszero(sc(bi)) 
          llc = root_update!(Ψp, Ψc, αc, σλc, σμc, llc, δt, srδt, 
            lλmxpr, lμmxpr, icr)
        end

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, αc, σλc, σμc, δt, srδt, nlim, icr)
      end
    end

    next!(pbar)
  end

  return Ψp, Ψc, llc, prc, αc, σλc, σμc
end




"""
    mcmc_gbmbd(Ψp      ::iTgbmbd,
               Ψc      ::iTgbmbd,
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               σμc     ::Float64,
               λa_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               μa_prior::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               σμ_prior::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               idf     ::Array{iBffs,1},
               pup     ::Array{Int64,1},
               nlim    ::Int64,
               prints  ::Int64,
               svf     ::Function)

MCMC chain for `gbmbd`.
"""
function mcmc_gbmbd(Ψp      ::iTgbmbd,
                    Ψc      ::iTgbmbd,
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    σμc     ::Float64,
                    λa_prior::NTuple{2,Float64},
                    μa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    σμ_prior::NTuple{2,Float64},
                    niter   ::Int64,
                    nthin   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idf     ::Array{iBffs,1},
                    pup     ::Array{Int64,1},
                    nlim    ::Int64,
                    prints  ::Int64,
                    svf     ::Function)

  # crown or stem conditioning
  icr = iszero(e(Ψc))

  lλmxpr = log(λa_prior[2])
  lμmxpr = log(μa_prior[2])

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 9)

  # make Ψ vector
  Ψv = iTgbmbd[]

  # number of branches and of triads
  nbr  = lastindex(idf)

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for i in Base.OneTo(niter)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update α
      if pupi === 1

        llc, prc, αc  = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nλ = sss_gbm(Ψ, αc)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) + scond(Ψ, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, pupi
        #    return 
        # end

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc = 
          update_σ!(σλc, σμc, αc, ssλ, ssμ, nλ, llc, prc, σλ_prior, σμ_prior)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) + scond(Ψ, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, pupi
        #    return 
        # end

      # gbm update
      elseif pupi === 3

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, ssμ = 
          update_gbm!(bix, Ψ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, sns, δt, 
            srδt, lλxpr, lμxpr)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) + scond(Ψ, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, pupi
        #    return 
        # end

      # forward simulation update
      else

        bi  = idf[ceil(Int64,rand()*nbr)]

        if iszero(sc(bi)) 
          llc = root_update!(Ψp, Ψc, αc, σλc, σμc, llc, δt, srδt, 
            lλmxpr, lμmxpr, icr)
        end

        Ψp, Ψc, llc = 
          fsp(Ψp, Ψc, bi, llc, αc, σλc, σμc, δt, srδt, nlim, icr)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, σμc, δt, srδt) + scond(Ψ, sns) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, pupi
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
        R[lit,4] = exp(lλ(Ψc)[1])
        R[lit,5] = exp(lμ(Ψc)[1])
        R[lit,6] = αc
        R[lit,7] = σλc
        R[lit,8] = σμc
        R[lit,9] = ntipsextinct(Ψc)
        push!(Ψv, deepcopy(Ψc))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ψv
end







"""
    update_fs!(bix    ::Int64,
               Ψ      ::Vector{iTgbmbd},
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
               sns    ::NTuple{3,BitVector},
               δt     ::Float64,
               srδt   ::Float64,
               snodes!::Function, 
               scond0 ::Function)

Forward simulation proposal function for `gbmbd`.
"""
function update_fs!(bix    ::Int64,
                    Ψ      ::Vector{iTgbmbd},
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
                    sns    ::NTuple{3,BitVector},
                    δt     ::Float64,
                    srδt   ::Float64,
                    snodes!::Function, 
                    scond0 ::Function)

  bi  = idf[bix]
  itb = it(bi) # if is terminal

  ψc  = Ψ[bix]
  if !itb
    ψ1  = Ψ[d1(bi)]
    ψ2  = Ψ[d2(bi)]
  end

  # forward simulate an internal branch
  ψp, np, ntp, λf = fsbi_ce(bi, lλ(ψc)[1], α, σλ, σμ, δt, srδt)

  # check for survival or non-exploding simulation
  if np > 0

    ρbi = ρi(bi) # get branch sampling fraction
    nc  = ni(bi) # current ni
    ntc = nt(bi) # current nt

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
          _daughters_update!(ψ1, ψ2, λf, α, σλ, σμ, δt, srδt)

        llr += llrd
        acr += acrd
      else
        return llc, dλ, ssλ, ssμ, nλ, ne, L
      end
    end

    # MH ratio
    if -randexp() < acr

      ll1, dλ1, ssλ1, nλ1 = llik_gbm_ssλ(ψp, α, σλ, σμ, δt, srδt)
      ll0, dλ0, ssλ0, nλ0 = llik_gbm_ssλ(ψc, α, σλ, σμ, δt, srδt)

      # if stem or crown conditioned
      scn = (iszero(pa(bi)) && e(bi) > 0.0) || 
             (isone(pa(bi)) && iszero(e(Ψ[1])))
      if scn
        llr += scond0(ψp, σμ, itb) - scond0(ψc, σμ, itb)
      end

      # update llr, ssλ, nλ, sns, ne, L,
      llr += ll1  - ll0
      dλ  += dλ1  - dλ0  + drλ
      ssλ += ssλ1 - ssλ0 + ssrλ
      nλ  += nλ1  - nλ0
      ne  += ntipsextinct(ψp) - ntipsextinct(ψc)
      L   += treelength(ψp)   - treelength(ψc)

      Ψ[bix] = ψp          # set new tree
      llc += llr           # set new likelihood
      if scn
        snodes!(Ψ, sns)    # set new sns
      end
      setni!(bi, np)       # set new ni
      setnt!(bi, ntp)      # set new nt
      setλt!(bi, λf)       # set new λt
      if !itb
        copyto!(lλ(ψ1), λ1p) # set new daughter 1 λ vector
        copyto!(lλ(ψ2), λ2p) # set new daughter 2 λ vector
      end
    end
  end

  return llc, dλ, ssλ, ssμ, nλ, ne, L
end




"""
    fsbi_ce(bi  ::iBffs,
            λ0  ::Float64,
            α   ::Float64,
            σλ  ::Float64,
            μ   ::Float64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_ce(bi  ::iBffs,
                 λ0  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 μ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

  # times
  tfb = tf(bi)

  # forward simulation during branch length
  t0, na, nsp = _sim_gbmbd(e(bi), λ0, α, σλ, σμ, δt, srδt, 0, 1, 1_000)

  if nsp >= 1_000
    return iTgbmbd(), 0, 0, 0.0
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
      tx, na = tip_sims!(t0, tfb, α, σλ, σμ, δt, srδt, na)
    end

    return t0, na, nat, λf
  end

  return iTgbmbd(), 0, 0, 0.0
end




"""
    tip_sims!(tree::iTgbmbd,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::iTgbmbd,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   na  ::Int64)

  if istip(tree) 
    if !isfix(tree) && isalive(tree)

      fdti = fdt(tree)
      lλ0  = lλ(tree)

      # simulate
      stree, na, nsp = 
        _sim_gbmbd(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, σμ, δt, srδt, 
                   na - 1, 1, 1_000)

      if !isdefined(stree, :lλ)
        return tree, 1_000
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
    tree.d1, na = tip_sims!(tree.d1, t, α, σλ, σμ, δt, srδt, na)
    tree.d2, na = tip_sims!(tree.d2, t, α, σλ, σμ, δt, srδt, na)
  end

  return tree, na
end




"""
    update_gbm!(bix  ::Int64,
                Ψ    ::Vector{iTgbmbd},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                σμ   ::Float64,
                llc  ::Float64,
                dλ   ::Float64,
                ssλ  ::Float64,
                ssμ  ::Float64,
                sns  ::NTuple{3,BitVector},
                δt   ::Float64,
                srδt ::Float64,
                lλxpr::Float64,
                lμxpr::Float64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ψ    ::Vector{iTgbmbd},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     σμ   ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     ssμ  ::Float64,
                     sns  ::NTuple{3,BitVector},
                     δt   ::Float64,
                     srδt ::Float64,
                     lλxpr::Float64,
                     lμxpr::Float64)

  @inbounds begin
    ψi   = Ψ[bix]
    bi   = idf[bix]
    ψ1   = Ψ[d1(bi)]
    ψ2   = Ψ[d2(bi)]
    ter1 = it(idf[d1(bi)]) 
    ter2 = it(idf[d2(bi)])

    cn = false
    # if crown root
    if iszero(pa(bi)) && iszero(e(bi))
      llc, dλ, ssλ, ssμ = 
        _crown_update!(ψi, ψ1, ψ2, α, σλ, σμ, llc, dλ, ssλ, ssμ, 
          δt, srδt, lλxpr, lμxpr)
      setλt!(bi, lλ(ψi)[1])

      # carry on updates in the crown daughter branches
      llc, dλ, ssλ, ssμ = 
        _update_gbm!(ψ1, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter1, 
          sns[2], 1)
      llc, dλ, ssλ, ssμ = 
        _update_gbm!(ψ2, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter2, 
          sns[3], 1)
    else
      # if stem branch
      if iszero(pa(bi))
        llc, dλ, ssλ, ssμ = 
          _stem_update!(ψi, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, 
            lλxpr, lμxpr)

        # updates within the stem branch in stem conditioning
        llc, dλ, ssλ, ssμ = 
          _update_gbm!(ψi, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, false, 
            sns[1], 1)

        # if observed node should be conditioned
        cn = sns[1][end]

      # if crown branch
      elseif isone(pa(bi)) && iszero(e(Ψ[1]))
        wsn = bix === d1(idf[pa(bi)]) ? 2 : 3
        sni = sns[wsn]
        # updates within the crown branch with crown conditioning
        llc, dλ, ssλ, ssμ = 
          _update_gbm!(ψi, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, false, 
            sni, 1)

        # if observed node should be conditioned
        if lastindex(sni) > 0
          cn = sni[end]
        end
      else
        # updates within the parent branch
        llc, dλ, ssλ, ssμ = _update_gbm!(ψi, α, σλ, σμ, llc, dλ, ssλ, ssμ, 
          δt, srδt, false)
      end

      # get fixed tip 
      lψi = fixtip(ψi) 

      # make between decoupled trees node update
      llc, dλ, ssλ, ssμ = update_triad!(lλ(lψi), lλ(ψ1), lλ(ψ2), 
        lμ(lψi), lμ(ψ1), lμ(ψ2), e(lψi), e(ψ1), e(ψ2), 
        fdt(lψi), fdt(ψ1), fdt(ψ2), α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, cn)

      # set fixed `λ(t)` in branch
      setλt!(bi, lλ(lψi)[end])

      # carry on updates in the daughters
      llc, dλ, ssλ, ssμ = 
        _update_gbm!(ψ1, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter1)
      llc, dλ, ssλ, ssμ = 
        _update_gbm!(ψ2, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, ter2)
    end
  end

  return llc, dλ, ssλ, ssμ
end




"""
    update_σ!(σλc     ::Float64,
              σμc     ::Float64,
              α       ::Float64,
              ssλ     ::Float64,
              ssμ     ::Float64,
              nλ      ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64},
              σμ_prior::NTuple{2,Float64})

Gibbs update for `σλ` and `σμ`.
"""
function update_σ!(σλc     ::Float64,
                   σμc     ::Float64,
                   α       ::Float64,
                   ssλ     ::Float64,
                   ssμ     ::Float64,
                   nλ      ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   σλ_prior::NTuple{2,Float64},
                   σμ_prior::NTuple{2,Float64})

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_prior[1] + 0.5 * nλ, σλ_prior[2] + ssλ)
  σμp2 = randinvgamma(σμ_prior[1] + 0.5 * nλ, σμ_prior[2] + ssμ)

  # update prior
  prc += llrdinvgamma(σλp2, σλc^2, σλ_prior[1], σλ_prior[2]) + 
         llrdinvgamma(σμp2, σμc^2, σμ_prior[1], σμ_prior[2])

  σλp = sqrt(σλp2)
  σμp = sqrt(σμp2)

  # update likelihood
  llc += ssλ*(1.0/σλc^2 - 1.0/σλp^2) - nλ*(log(σλp/σλc)) + 
         ssμ*(1.0/σμc^2 - 1.0/σμp^2) - nλ*(log(σμp/σμc))

  return llc, prc, σλp, σμp
end







































# """
#     gbm!(Ψp   ::iTgbmbd,
#          Ψc   ::iTgbmbd,
#          bi   ::iBffs,
#          llc  ::Float64,
#          α    ::Float64,
#          σλ   ::Float64,
#          μ    ::Float64,
#          δt   ::Float64,
#          srδt ::Float64)

# Update the gbm for branch `bi`.
# """
# function gbm!(Ψp   ::iTgbmbd,
#               Ψc   ::iTgbmbd,
#               bi   ::iBffs,
#               llc  ::Float64,
#               α    ::Float64,
#               σλ   ::Float64,
#               σμ    ::Float64,
#               δt   ::Float64,
#               srδt ::Float64)

#   # get branch information
#   dri = dr(bi)
#   ldr = lastindex(dri)
#   itb = it(bi)

#   # go to branch to be updated
#   treec, treep = drtree(Ψc, Ψp, dri, ldr, 0)

#   llc = gbm_update!(treep, treec, α, σλ, σμ, llc, δt, srδt)

#   return llc
# end




# """
#     fsp(Ψp   ::iTgbmbd,
#         Ψc   ::iTgbmbd,
#         bi   ::iBffs,
#         llc  ::Float64,
#         α    ::Float64, 
#         σλ   ::Float64, 
#         σμ   ::Float64,
#         tsv  ::Array{Array{Float64,1},1},
#         bbλp ::Array{Array{Float64,1},1}, 
#         bbμp ::Array{Array{Float64,1},1}, 
#         bbλc ::Array{Array{Float64,1},1}, 
#         bbμc ::Array{Array{Float64,1},1}, 
#         bix  ::Int64,
#         triad::Array{Int64,1},
#         ter  ::BitArray{1},
#         δt   ::Float64, 
#         srδt ::Float64,
#         nlim ::Int64,
#         icr  ::Bool, 
#         wbc  ::Int64)

# Forward simulation proposal function for `gbmbd`.
# """
# function fsp(Ψp   ::iTgbmbd,
#              Ψc   ::iTgbmbd,
#              bi   ::iBffs,
#              llc  ::Float64,
#              α    ::Float64, 
#              σλ   ::Float64, 
#              σμ   ::Float64,
#              δt   ::Float64, 
#              srδt ::Float64,
#              nlim ::Int64,
#              icr  ::Bool)

#   # get branch information
#   dri = dr(bi)
#   ldr = lastindex(dri)
#   itb = it(bi)

#   # go to branch to be updated
#   treec, treep = drtree(Ψc, Ψp, dri, ldr, 0)

#   # forward simulation
#   t0, ret, λf, μf, dft0 = 
#     fsbi(bi, lλ(treec)[1], lμ(treec)[1], α, σλ, σμ, δt, srδt, nlim)

#   # if retain simulation
#   if ret

#     if !itb
#       # make daughter proposal to be concordant with `t0`
#       llr, acr = daughters_lprop!(treep, treec, λf, μf, α, σλ, σμ, δt, srδt)

#       # change last event by speciation for llr
#       iλ = λf

#       # get previous ending λ
#       treecd1, treecd2 = fixds(treec)

#       # acceptance ratio
#       acr += λf - lλ(treecd1)[1]

#     else
#       iλ  = 0.0
#       llr = 0.0
#       acr = 0.0
#     end

#     # mh ratio
#     if -randexp() < acr 

#       llr += llik_gbm(     t0, α, σλ, σμ, δt, srδt) + iλ - 
#              llik_gbm_f(treec, α, σλ, σμ, δt, srδt)

#       if icr && isone(sc(bi))
#         css = itb ? cond_surv_stem : cond_surv_stem_p
#         if dri[1]
#           llr += css(t0) - cond_surv_stem(treec)
#         else
#           llr += css(t0) - cond_surv_stem(treec)
#         end
#       elseif iszero(sc(bi))
#         llr += cond_surv_stem_p(t0) -
#                cond_surv_stem(treec) 
#       end

#       llc += llr

#       # copy daughters vectors
#       if !itb
#         treepd1, treepd2 = fixds(treep)

#         copyto!(lλ(treecd1), lλ(treepd1))
#         copyto!(lμ(treecd1), lμ(treepd1))
#         copyto!(lλ(treecd2), lλ(treepd2))
#         copyto!(lμ(treecd2), lμ(treepd2))
#       end

#       # make combined swap branch
#       Ψp, Ψc = swapbranch!(Ψp, Ψc, t0, dri, ldr, itb, 0)
#     end
#   end

#   return Ψp, Ψc, llc
# end




# """
#     fsbi(bi  ::iBffs, 
#          iλ  ::Float64, 
#          iμ  ::Float64, 
#          α   ::Float64, 
#          σλ  ::Float64, 
#          σμ  ::Float64, 
#          δt  ::Float64, 
#          srδt::Float64,
#          nlim::Int64)

# Forward `gbmbd` simulation for branch `bi`.
# """
# function fsbi(bi  ::iBffs, 
#               iλ  ::Float64, 
#               iμ  ::Float64, 
#               α   ::Float64, 
#               σλ  ::Float64, 
#               σμ  ::Float64, 
#               δt  ::Float64, 
#               srδt::Float64,
#               nlim::Int64)

#   # retain the simulation?
#   ret = true

#   # times
#   tfb = tf(bi)

#   # simulate tree
#   t0, nsp = _sim_gbmbd(ti(bi) - tfb, iλ, iμ, α, σλ, σμ, δt, srδt, 1, nlim)

#   na = ntipsalive(t0)

#   λf, μf, dft0 = NaN, NaN, NaN

#   # if simulation goes extinct or maximum number of species reached
#   if iszero(na) || nsp === nlim
#     ret = false
#   # if one surviving lineage
#   elseif isone(na)
#     f, λf, μf, dft0 = fixalive!(t0, NaN, NaN, NaN)
#   elseif na > 1
#     # if terminal branch
#     if it(bi)
#       ret = false
#     # if continue the simulation
#     else
#       # fix random tip and return end λ(t) and μ(t) 
#       λf, μf, dft0 = fixrtip!(t0, na, NaN, NaN, NaN)

#       for j in Base.OneTo(na - 1)
#         # get their final λ and μ to continue forward simulation
#         ix, λt, μt, fdti = fλμ1(t0, NaN, NaN, NaN, false)

#         for i in Base.OneTo(2)
#           st0, nsp = 
#             _sim_gbmbd(max(δt - fdti, 0.0), tfb, λt, μt, α, σλ, σμ, 
#                       δt, srδt, 1, nlim)
#           # if maximum number of species reached.
#           if nsp === nlim
#             if i === 2
#               ret = false
#             end
#             continue
#           end
#           # if goes extinct before the present
#           if iszero(ntipsalive(st0))
#             # graft to tip
#             addtotip(t0, st0, false)
#             break
#           end
#           # if not succeeded after 2 tries.
#           if i === 2
#             ret = false
#           end
#         end
#         !ret && break
#       end
#     end
#   end

#   # speciates at time `tfb`
#   if iszero(dft0)
#     ret = false
#   end

#   return t0, ret, λf, μf, dft0
# end




# """
#     add1(tree::iTgbmbd, stree::iTgbmbd, it::Int64, ix::Int64)

# Add `stree` to tip in `tree` given by `it` in `tree.d1` order.
# """
# function addtotip(tree::iTgbmbd, stree::iTgbmbd, ix::Bool) 

#   if istip(tree) 
#     if isalive(tree) && !isfix(tree)

#       sete!(tree, e(tree) + e(stree))
#       setproperty!(tree, :iμ, isextinct(stree))

#       lλ0 = lλ(tree)
#       lμ0 = lμ(tree)
#       lλs = lλ(stree)
#       lμs = lμ(stree)

#       if lastindex(lλs) === 2
#         setfdt!(tree, fdt(tree) + fdt(stree))
#       else
#         setfdt!(tree, fdt(stree))
#       end

#       pop!(lλ0)
#       pop!(lμ0)
#       popfirst!(lλs)
#       popfirst!(lμs)
#       append!(lλ0, lλs)
#       append!(lμ0, lμs)

#       if isdefined(stree, :d1)
#         tree.d1 = stree.d1
#         tree.d2 = stree.d2
#       end

#       ix = true
#     end

#     return ix 
#   end

#   if !ix
#     ix = addtotip(tree.d1, stree, ix)
#   end
#   if !ix
#     ix = addtotip(tree.d2, stree, ix)
#   end

#   return ix
# end




# """
#     fixrtip!(tree::iTgbmbd, 
#              na  ::Int64, 
#              λf  ::Float64, 
#              μf  ::Float64,
#              dft0::Float64) 

# Fixes the the path for a random non extinct tip.
# """
# function fixrtip!(tree::iTgbmbd, 
#                   na  ::Int64, 
#                   λf  ::Float64, 
#                   μf  ::Float64,
#                   dft0::Float64) 

#   fix!(tree)

#   if isdefined(tree, :d1)
#     if isextinct(tree.d1)
#       λf, μf, dft0 = 
#         fixrtip!(tree.d2, na, λf, μf, dft0)
#     elseif isextinct(tree.d2)
#       λf, μf, dft0 = 
#         fixrtip!(tree.d1, na, λf, μf, dft0)
#     else
#       na1 = ntipsalive(tree.d1)
#       # probability proportional to number of lineages
#       if (fIrand(na) + 1) > na1
#         λf, μf, dft0 = 
#           fixrtip!(tree.d2, na - na1, λf, μf, dft0)
#       else
#         λf, μf, dft0 = 
#           fixrtip!(tree.d1, na1, λf, μf, dft0)
#       end
#     end
#   else

#     dft0 = fdt(tree)
#     λv   = lλ(tree)
#     μv   = lμ(tree)
#     l    = lastindex(λv)
#     λf   = λv[l]
#     μf   = μv[l]
#   end

#   return λf, μf, dft0
# end




# """
#     fixalive!(tree::iTgbmbd,
#               λf  ::Float64,
#               μf  ::Float64,
#               dft0::Float64)

# Fixes the the path from root to the only species alive.
# """
# function fixalive!(tree::iTgbmbd,
#                    λf  ::Float64,
#                    μf  ::Float64,
#                    dft0::Float64)

#   if istip(tree) 
#     if isalive(tree)
#       fix!(tree)
#       dft0 = fdt(tree)
#       λv   = lλ(tree)
#       μv   = lμ(tree)
#       l    = lastindex(λv)
#       λf   = λv[l]
#       μf   = μv[l]

#       return true, λf, μf, dft0
#     end
#   else
#     f, λf, μf, dft0 = 
#       fixalive!(tree.d2, λf, μf, dft0)
#     if f 
#       fix!(tree)
#       return true, λf, μf, dft0
#     end
#     f, λf, μf, dft0 = 
#       fixalive!(tree.d1, λf, μf, dft0)
#     if f 
#       fix!(tree)
#       return true, λf, μf, dft0
#     end
#   end

#   return false, λf, μf, dft0
# end




# """
#     fλμ1(tree::iTgbmbd, 
#          λt  ::Float64, 
#          μt  ::Float64, 
#          fdti::Float64,
#          ix  ::Bool)

# Get end `λ` and `μ` for a tip in `tree` given in `tree.d1` order
# not taking into account the fixed tip.
# """
# function fλμ1(tree::iTgbmbd, 
#               λt  ::Float64, 
#               μt  ::Float64, 
#               fdti::Float64,
#               ix  ::Bool)

#   if istip(tree) 
#     if isalive(tree) &&!isfix(tree)
#       @inbounds begin
#         lλv  = lλ(tree)
#         l    = lastindex(lλv)
#         λt   = lλv[l]
#         μt   = lμ(tree)[l]
#         fdti = fdt(tree)
#       end
#       ix = true
#     end

#     return ix, λt, μt, fdti
#   end

#   if !ix
#     ix, λt, μt, fdti = fλμ1(tree.d1, λt, μt, fdti, ix)
#   end
#   if !ix
#     ix, λt, μt, fdti = fλμ1(tree.d2, λt, μt, fdti, ix)
#   end

#   return ix, λt, μt, fdti
# end






