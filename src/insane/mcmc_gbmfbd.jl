#=

Anagenetic `gbmfbd` MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    insane_gbmfbd(tree    ::sTf_label;
                  λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                  μa_prior::NTuple{2,Float64} = (0.0, 100.0),
                  α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                  σλ_prior::NTuple{2,Float64} = (3.0, 0.5),
                  σμ_prior::NTuple{2,Float64} = (3.0, 0.5),
                  ψ_prior ::NTuple{2,Float64} = (1.0, 1.0),
                  ψ_epoch ::Vector{Float64}   = Float64[],
                  f_epoch ::Vector{Int64}     = Int64[0],
                  niter   ::Int64             = 1_000,
                  nthin   ::Int64             = 10,
                  nburn   ::Int64             = 200,
                  nflush  ::Int64             = nthin,
                  ofile   ::String            = string(homedir(), "/ifbd"),
                  ϵi      ::Float64           = 0.2,
                  λi      ::Float64           = NaN,
                  μi      ::Float64           = NaN,
                  ψi      ::Float64           = NaN,
                  αi      ::Float64           = 0.0,
                  σλi     ::Float64           = 0.1,
                  σμi     ::Float64           = 0.1,
                  pupdp   ::NTuple{5,Float64} = (0.01, 0.01, 0.01, 0.1, 0.2),
                  δt      ::Float64           = 1e-3,
                  survival::Bool              = true,
                  mxthf   ::Float64           = Inf,
                  prints  ::Int64             = 5,
                  tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for fossilized birth-death diffusion `fbdd`.
"""
function insane_gbmfbd(tree    ::sTf_label;
                       λa_prior::NTuple{2,Float64} = (0.0, 100.0),
                       μa_prior::NTuple{2,Float64} = (0.0, 100.0),
                       α_prior ::NTuple{2,Float64} = (0.0, 10.0),
                       σλ_prior::NTuple{2,Float64} = (3.0, 0.5),
                       σμ_prior::NTuple{2,Float64} = (3.0, 0.5),
                       ψ_prior ::NTuple{2,Float64} = (1.0, 1.0),
                       ψ_epoch ::Vector{Float64}   = Float64[],
                       f_epoch ::Vector{Int64}     = Int64[0],
                       niter   ::Int64             = 1_000,
                       nthin   ::Int64             = 10,
                       nburn   ::Int64             = 200,
                       nflush  ::Int64             = nthin,
                       ofile   ::String            = string(homedir(), "/ifbd"),
                       ϵi      ::Float64           = 0.2,
                       λi      ::Float64           = NaN,
                       μi      ::Float64           = NaN,
                       ψi      ::Float64           = NaN,
                       αi      ::Float64           = 0.0,
                       σλi     ::Float64           = 0.1,
                       σμi     ::Float64           = 0.1,
                       pupdp   ::NTuple{5,Float64} = (0.01, 0.01, 0.01, 0.1, 0.2),
                       δt      ::Float64           = 1e-3,
                       survival::Bool              = true,
                       mxthf   ::Float64           = Inf,
                       prints  ::Int64             = 5,
                       tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  # `n` tips, `th` treeheight define δt
  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1,round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # only include epochs where the tree occurs
  sort!(ψ_epoch, rev = true)
  tix = findfirst(x -> x < th, ψ_epoch)
  if !isnothing(tix)
    ψ_epoch = ψ_epoch[tix:end]
  end
  nep  = lastindex(ψ_epoch) + 1

  # make initial fossils per epoch vector
  if lastindex(f_epoch) !== nep
    f_epoch = fill(0, nep)
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

  # starting parameters
  if isnan(λi) && isnan(μi) && isnan(ψi)
    # if only one tip
    if isone(n)
      λc = prod(λ_prior)
      μc = prod(μ_prior)
    else
      λc, μc = moments(Float64(n), th, ϵi)
    end
    # if no sampled fossil
    nf = nfossils(tree)
    if iszero(nf)
      ψc = prod(ψ_prior)
    else
      ψc = Float64(nf)/Float64(treelength(tree))
    end
  else
    λc, μc, ψc = λi, μi, ψi
  end

  # make ψ vector
  ψc = fill(ψc, nep)

  # condition on survival of 0, 1, or 2 starting lineages
  surv = 0
  if survival 
    if iszero(e(tree)) 
      if def1(tree)
        surv += (ntipsalive(tree.d1) > 0)
        if def2(tree)
          surv += (ntipsalive(tree.d2) > 0)
        end
      end
    else
      surv += (ntipsalive(tree) > 0)
    end
  end

  # M attempts of survival
  mc = m_surv_gbmbd(th, log(λc), log(μc), αi, σλi, σμi, δt, srδt, 1_000, surv)

  # make a decoupled tree
  Ξ = make_Ξ(idf, λc, μc, αi, σλi, σμi, δt, srδt, iTfbd)

  # set end of fix branch speciation times and get vector of internal branches
  # and make epoch start vectors and indices for each `ξ`
  inodes = Int64[]
  eixi   = Int64[]
  eixf   = Int64[]
  bst    = Float64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    if d1(bi) > 0
      push!(inodes, i)
    end
    tib = ti(bi)
    ei  = findfirst(x -> x < tib, ψ_epoch)
    ei  = isnothing(ei) ? nep : ei
    ef  = findfirst(x -> x < tf(bi), ψ_epoch)
    ef  = isnothing(ef) ? nep : ef
    push!(bst, tib)
    push!(eixi, ei)
    push!(eixf, ef)
  end

  # parameter updates (1: α, 2: σλ & σμ, 3: ψ, 4: gbm, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running fossilized birth-death diffusion"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, σμc, ψc, mc  =
    mcmc_burn_gbmfbd(Ξ, idf, λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior,
      ψ_prior, ψ_epoch, f_epoch, nburn, αi, σλi, σμi, ψc, mc, th, surv, 
      δt, srδt, bst, eixi, eixf, inodes, pup, prints)

  # mcmc
  r, treev =
    mcmc_gbmfbd(Ξ, idf, llc, prc, αc, σλc, σμc, ψc, mc, th, surv,
      λa_prior, μa_prior, α_prior, σλ_prior, σμ_prior, 
      ψ_prior, ψ_epoch, f_epoch, δt, srδt, bst, eixi, eixf, inodes, pup, 
      niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_gbmfbd(Ξ       ::Vector{iTfbd},
                     idf     ::Vector{iBffs},
                     λa_prior::NTuple{2,Float64},
                     μa_prior::NTuple{2,Float64},
                     α_prior ::NTuple{2,Float64},
                     σλ_prior::NTuple{2,Float64},
                     σμ_prior::NTuple{2,Float64},
                     ψ_prior ::NTuple{2,Float64},
                     ψ_epoch ::Vector{Float64},
                     f_epoch ::Vector{Int64},
                     nburn   ::Int64,
                     αc      ::Float64,
                     σλc     ::Float64,
                     σμc     ::Float64,
                     ψc      ::Vector{Float64},
                     mc      ::Float64,
                     th      ::Float64,
                     surv   ::Int64,
                     δt      ::Float64,
                     srδt    ::Float64,
                     bst     ::Vector{Float64},
                     eixi    ::Vector{Int64},
                     eixf    ::Vector{Int64},
                     inodes  ::Array{Int64,1},
                     pup     ::Array{Int64,1},
                     prints  ::Int64)


MCMC burn-in chain for `fbdd`.
"""
function mcmc_burn_gbmfbd(Ξ       ::Vector{iTfbd},
                         idf     ::Vector{iBffs},
                         λa_prior::NTuple{2,Float64},
                         μa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         σμ_prior::NTuple{2,Float64},
                         ψ_prior ::NTuple{2,Float64},
                         ψ_epoch ::Vector{Float64},
                         f_epoch ::Vector{Int64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         σμc     ::Float64,
                         ψc      ::Vector{Float64},
                         mc      ::Float64,
                         th      ::Float64,
                         surv    ::Int64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         bst     ::Vector{Float64},
                         eixi    ::Vector{Int64},
                         eixf    ::Vector{Int64},
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  λ0  = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) -
        (iszero(e(Ξ[1])) && !isfossil(idf[1])) * λ0 + log(mc) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2,        σλ_prior[1], σλ_prior[2])  +
        logdinvgamma(σμc^2,        σμ_prior[1], σμ_prior[2])  +
        logdnorm(αc,               α_prior[1],  α_prior[2]^2) +
        logdunif(exp(λ0),          λa_prior[1], λa_prior[2])  +
        logdunif(exp(lμ(Ξ[1])[1]), μa_prior[1], μa_prior[2])  +
        sum(logdgamma.(ψc, ψ_prior[1], ψ_prior[2]))

  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L            = treelength(Ξ, ψ_epoch, bst, eixi) # tree length
  nf           = nfossils(idf, ψ_epoch, f_epoch)   # number of fossilization events per epoch
  dλ           = deltaλ(Ξ)                         # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ξ, αc)                    # sum squares in λ and μ
  nin          = lastindex(inodes)                 # number of internal nodes
  el           = lastindex(idf)                    # number of branches
  nep          = lastindex(ψc)

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for i in Base.OneTo(nburn)

    shuffle!(pup)

    # parameter updates
    for pupi in pup

      # update α
      if pupi === 1

        llc, prc, αc, mc  =
          update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, sum(L), dλ, llc, prc,
            mc, th, surv, δt, srδt, α_prior)

        # update ssλ with new drift `α`
        ssλ, ssμ, nλ = sss_gbm(Ξ, αc)

      # σλ & σμ update
      elseif pupi === 2

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], ssλ, ssμ, nλ,
            llc, prc, mc, th, surv, δt, srδt, σλ_prior, σμ_prior)

      # psi update
      elseif pupi === 3

        llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # gbm update
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, dλ, ssλ, ssμ, mc =
          update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
            surv, δt, srδt, lλxpr, lμxpr)

      # forward simulation update
      else

        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, ssμ, nλ, L =
          update_fs!(bix, Ξ, idf, αc, σλc, σμc, ψc, llc, dλ, ssλ, ssμ, nλ, L,
            ψ_epoch, δt, srδt, eixi, eixf)

      end
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, σμc, ψc, mc
end




"""
    mcmc_gbmfbd(Ξ       ::Vector{iTfbd},
                idf     ::Vector{iBffs},
                llc     ::Float64,
                prc     ::Float64,
                αc      ::Float64,
                σλc     ::Float64,
                σμc     ::Float64,
                ψc      ::Vector{Float64},
                mc      ::Float64,
                th      ::Float64,
                surv   ::Int64,
                λa_prior::NTuple{2,Float64},
                μa_prior::NTuple{2,Float64},
                α_prior ::NTuple{2,Float64},
                σλ_prior::NTuple{2,Float64},
                σμ_prior::NTuple{2,Float64},
                ψ_prior ::NTuple{2,Float64},
                ψ_epoch ::Vector{Float64},
                f_epoch ::Vector{Int64},
                δt      ::Float64,
                srδt    ::Float64,
                bst     ::Vector{Float64},
                eixi    ::Vector{Int64},
                eixf    ::Vector{Int64},
                inodes  ::Array{Int64,1},
                pup     ::Vector{Int64},
                niter   ::Int64,
                nthin   ::Int64,
                nflush  ::Int64,
                ofile   ::String,
                prints  ::Int64)

MCMC chain for `gbmbd`.
"""
function mcmc_gbmfbd(Ξ       ::Vector{iTfbd},
                     idf     ::Vector{iBffs},
                     llc     ::Float64,
                     prc     ::Float64,
                     αc      ::Float64,
                     σλc     ::Float64,
                     σμc     ::Float64,
                     ψc      ::Vector{Float64},
                     mc      ::Float64,
                     th      ::Float64,
                     surv    ::Int64,
                     λa_prior::NTuple{2,Float64},
                     μa_prior::NTuple{2,Float64},
                     α_prior ::NTuple{2,Float64},
                     σλ_prior::NTuple{2,Float64},
                     σμ_prior::NTuple{2,Float64},
                     ψ_prior ::NTuple{2,Float64},
                     ψ_epoch ::Vector{Float64},
                     f_epoch ::Vector{Int64},
                     δt      ::Float64,
                     srδt    ::Float64,
                     bst     ::Vector{Float64},
                     eixi    ::Vector{Int64},
                     eixf    ::Vector{Int64},
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

  # crown or crown conditioning
  lλxpr = log(λa_prior[2])
  lμxpr = log(μa_prior[2])

  L            = treelength(Ξ, ψ_epoch, bst, eixi) # tree length
  nf           = nfossils(idf, ψ_epoch, f_epoch)   # number of fossilization events per epoch
  dλ           = deltaλ(Ξ)                         # delta change in λ
  ssλ, ssμ, nλ = sss_gbm(Ξ, αc)                    # sum squares in λ and μ
  nin          = lastindex(inodes)                 # number of internal nodes
  el           = lastindex(idf)                    # number of branches
  nep          = lastindex(ψc)

  # parameter results
  nep = lastindex(ψc)
  r = Array{Float64,2}(undef, nlogs, 8 + nep)

  # make Ξ vector
  treev = iTfbd[]

  # number of branches and of triads
  nbr  = lastindex(idf)

  # flush to file
  sthin = 0

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\tmu_root\talpha\tsigma_lambda\tsigma_mu\t"*join(["psi"*(isone(nep) ? "" : string("_",i)) for i in 1:nep], "\t")*"\n")
    flush(of)

    open(ofile*".txt", "w") do tf


      pbar = Progress(niter, prints, "running mcmc...", 20)

      for it in Base.OneTo(niter)

        shuffle!(pup)

        # parameter updates
        for pupi in pup

          # update α
          if pupi === 1

            llc, prc, αc, mc  =
              update_α!(αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], σλc, σμc, sum(L), dλ, llc, prc,
                mc, th, surv, δt, srδt, α_prior)

            # update ssλ with new drift `α`
            ssλ, ssμ, nλ = sss_gbm(Ξ, αc)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) - (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
            #  if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi, Ξ
            #    @show bix
            #    return
            # end

          # σλ & σμ update
          elseif pupi === 2

            llc, prc, σλc, σμc, mc =
              update_σ!(σλc, σμc, αc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], ssλ, ssμ, nλ,
                llc, prc, mc, th, surv, δt, srδt, σλ_prior, σμ_prior)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) - (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
            #  if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi, Ξ
            #    @show bix
            #    return
            # end

          # psi update
          elseif pupi === 3

            llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) - (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
            #  if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi, Ξ
            #    @show bix
            #    return
            # end

          # gbm update
          elseif pupi === 4

            nix = ceil(Int64,rand()*nin)
            bix = inodes[nix]

            llc, dλ, ssλ, ssμ, mc =
              update_gbm!(bix, Ξ, idf, αc, σλc, σμc, llc, dλ, ssλ, ssμ, mc, th,
                surv, δt, srδt, lλxpr, lμxpr)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) - (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
            #  if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi, Ξ
            #    @show bix
            #    return
            # end

          # forward simulation update
          else

            bix = ceil(Int64,rand()*el)

            llc, dλ, ssλ, ssμ, nλ, L =
              update_fs!(bix, Ξ, idf, αc, σλc, σμc, ψc, llc, dλ, ssλ, ssμ, nλ, L,
                ψ_epoch, δt, srδt, eixi, eixf)

            # ll0 = llik_gbm(Ξ, idf, αc, σλc, σμc, ψc, ψ_epoch, bst, eixi, δt, srδt) - (iszero(e(Ξ[1])) && !isfossil(idf[1])) * lλ(Ξ[1])[1] + log(mc) + prob_ρ(idf)
            #  if !isapprox(ll0, llc, atol = 1e-4)
            #    @show ll0, llc, it, pupi, Ξ
            #    @show bix
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
            @turbo for i in Base.OneTo(nep)
              r[lit,8 + i] = ψc[i]
            end
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthin += 1
        if sthin === nflush
          write(of, 
            string(Float64(it), "\t", llc, "\t", prc, "\t", 
              exp(lλ(Ξ[1])[1]),"\t", exp(lμ(Ξ[1])[1]), "\t", αc, "\t",
               σλc, "\t", σμc, "\t", join(ψc, "\t"), "\n"))
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
    update_fs!(bix ::Int64,
               Ξ   ::Vector{iTfbd},
               idf ::Vector{iBffs},
               α   ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               ψ   ::Vector{Float64},
               llc ::Float64,
               dλ  ::Float64,
               ssλ ::Float64,
               ssμ ::Float64,
               nλ  ::Float64,
               L   ::Vector{Float64},
               ψts ::Vector{Float64},
               δt  ::Float64,
               srδt::Float64,
               eixi::Vector{Int64},
               eixf::Vector{Int64})

Forward simulation proposal function for `gbmfbd`.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{iTfbd},
                    idf ::Vector{iBffs},
                    α   ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Vector{Float64},
                    llc ::Float64,
                    dλ  ::Float64,
                    ssλ ::Float64,
                    ssμ ::Float64,
                    nλ  ::Float64,
                    L   ::Vector{Float64},
                    ψts ::Vector{Float64},
                    δt  ::Float64,
                    srδt::Float64,
                    eixi::Vector{Int64},
                    eixf::Vector{Int64})

  bi  = idf[bix]
  ξc  = Ξ[bix]
  ixi = eixi[bix]

  # terminal branch
  if iszero(d1(bi))

    # fossil terminal branch
    if isfossil(bi)
      ixf = eixf[bix]

      ξp, llr = fsbi_t(bi, ξc, α, σλ, σμ, ψ, ψts, ixi, ixf, δt, srδt)

      # if terminal but not successful proposal, update extinct
      if !isfinite(llr)
        ξp, llr = fsbi_et(iTfbd_wofe(ξc), bi, α, σλ, σμ, ψ, ψts, ixf,
          δt, srδt)
      end

    # non-fossil terminal branch
    else
      ξp, llr = fsbi_t(bi, ξc, α, σλ, σμ, ψ, ψts, ixi, δt, srδt)
    end
    drλ = ssrλ = ssrμ = 0.0

  # internal non-bifurcating branch
  elseif iszero(d2(bi))

    ξp, llr, drλ, ssrλ, ssrμ =
      fsbi_m(bi, ξc, Ξ[d1(bi)], α, σλ, σμ, ψ, ψts, ixi, eixf[bix], δt, srδt)

  # internal bifurcating branch
  else

    ξp, llr, drλ, ssrλ, ssrμ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, σμ, ψ, ψts, 
        ixi, eixf[bix], δt, srδt)
  end

  if isfinite(llr)
    tii = ti(bi)

    nep = lastindex(ψts) + 1

    ll1, ixd, dλ1, ssλ1, ssμ1, nλ1 =
      llik_gbm_ss(ξp, α, σλ, σμ, ψ, tii, ψts, ixi, δt, srδt, nep)
    ll0, ixd, dλ0, ssλ0, ssμ0, nλ0 =
      llik_gbm_ss(ξc, α, σλ, σμ, ψ, tii, ψts, ixi, δt, srδt, nep)

    # update llr, ssλ, ssμ, nλ, sns
    llc += ll1  - ll0  + llr
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    ssμ += ssμ1 - ssμ0 + ssrμ
    nλ  += nλ1  - nλ0
    # update tree lengths
    Lc = zeros(Float64, nep)
    _treelength!(ξc, tii, Lc, ψts, ixi, nep)
    _treelength!(ξp, tii, L,  ψts, ixi, nep)
    @turbo for i in Base.OneTo(nep)
      L[i] -= Lc[i]
    end

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, dλ, ssλ, ssμ, nλ, L
end




"""
    fsbi_t(bi::iBffs,
           ξc  ::iTfbd,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Vector{Float64},
           ψts ::Vector{Float64},
           ix  ::Int64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs,
                ξc  ::iTfbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ix  ::Int64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nn, llr =
    _sim_gbmfbd_t(e(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ, ψts, ix, nep,
      δt, srδt, lc, lU, Iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr)
    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, -Inf
  end
end




"""
     fsbi_t(bi  ::iBffs,
            ξc  ::iTfbd,
            α   ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            ψ   ::Vector{Float64},
            ψts ::Vector{Float64},
            ix  ::Int64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTfbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ixi ::Int64,
                ixf ::Int64,
                δt  ::Float64,
                srδt::Float64)

   # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ,
      ψts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return t0, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if lU < acr

    _fixrtip!(t0, na) # fix random tip

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr

      # fossilize extant tip
      fossilizefixedtip!(t0)

      # if terminal fossil branch
      tx, na, nn, acr =
        fossiltip_sim!(t0, tf(bi), α, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)

      if lU < acr

        llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
        setnt!(bi, ntp)      # set new nt
        setni!(bi, na)       # set new ni

        return t0, llr
      end
    end
  end

  return t0, NaN
end




"""
    fsbi_et(t0  ::iTfbd,
            bi  ::iBffs,
            α   ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            ψ   ::Vector{Float64},
            ψts ::Vector{Float64},
            ixf ::Int64,
            δt  ::Float64,
            srδt::Float64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_et(t0  ::iTfbd,
                 bi  ::iBffs,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 ψ   ::Vector{Float64},
                 ψts ::Vector{Float64},
                 ixf ::Int64,
                 δt  ::Float64,
                 srδt::Float64)

  nep = lastindex(ψts) + 1
  lU  = -randexp()            # log-probability
  nac = ni(bi)                # current ni
  Iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # if terminal fossil branch
  tx, na, nn, acr =
    fossiltip_sim!(t0, tf(bi), α, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
      acr, lU, Iρi, 1, 1)

  if lU < acr

    llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
    setni!(bi, na)       # set new ni

    return t0, llr
  end

  return t0, NaN
end



"""
    fsbi_m(bi  ::iBffs,
           ξc  ::iTfbd,
           ξ1  ::iTfbd,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Vector{Float64},
           ψts ::Vector{Float64},
           ixi ::Int64,
           ixf ::Int64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for fossil internal branch `bi`.
"""
function fsbi_m(bi  ::iBffs,
                ξc  ::iTfbd,
                ξ1  ::iTfbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ixi ::Int64,
                ixf ::Int64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ,
      ψts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample and fix random  tip
  λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, ssrμ, λ1p, μ1p, =
    _daughter_update!(ξ1, λf, μf, α, σλ, σμ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      # fossilize extant tip
      isfossil(bi) && fossilizefixedtip!(t0)

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setnt!(bi, ntp)                       # set new nt
      setni!(bi, na)                        # set new ni
      l1 = lastindex(λ1p)
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector

      return t0, llr, drλ, ssrλ, ssrμ
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    fsbi_i(bi  ::iBffs,
           ξc  ::iTfbd,
           ξ1  ::iTfbd,
           ξ2  ::iTfbd,
           α   ::Float64,
           σλ  ::Float64,
           σμ  ::Float64,
           ψ   ::Vector{Float64},
           ψts ::Vector{Float64},
           ixi ::Int64,
           ixf ::Int64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for internal branch `bi`.
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTfbd,
                ξ1  ::iTfbd,
                ξ2  ::iTfbd,
                α   ::Float64,
                σλ  ::Float64,
                σμ  ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ixi ::Int64,
                ixf ::Int64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nf, nn =
    _sim_gbmfbd_i(ti(bi), tf(bi), lλ(ξc)[1], lμ(ξc)[1], α, σλ, σμ, ψ,
      ψts, ixi, nep, δt, srδt, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # sample and fix random  tip
  λf, μf = fixrtip!(t0, na, NaN, NaN) # fix random tip

  llrd, acrd, drλ, ssrλ, ssrμ, λ1p, λ2p, μ1p, μ2p =
    _daughters_update!(ξ1, ξ2, λf, μf, α, σλ, σμ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), α, σλ, σμ, ψ, ψts, ixf, nep, δt, srδt,
          acr, lU, Iρi, na, nn)
    end

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setnt!(bi, ntp)                       # set new nt
      setni!(bi,  na)                       # set new ni
      setλt!(bi,  λf)                       # set new λt
      l1 = lastindex(λ1p)
      l2 = lastindex(λ2p)
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 1 λ vector
      unsafe_copyto!(lμ(ξ1), 1, μ1p, 1, l1) # set new daughter 1 μ vector
      unsafe_copyto!(lμ(ξ2), 1, μ2p, 1, l2) # set new daughter 1 μ vector

      return t0, llr, drλ, ssrλ, ssrμ
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTfbd,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              ψ   ::Vector{Float64},
              ψts ::Vector{Float64},
              ix  ::Int64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTfbd,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   ψ   ::Vector{Float64},
                   ψts ::Vector{Float64},
                   ix  ::Int64,
                   nep ::Int64,
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
          _sim_gbmfbd_it(max(δt-fdti, 0.0), t, lλ0[l], lμ0[l], α, σλ, σμ, ψ,
            ψts, ix, nep, δt, srδt, lr, lU, Iρi, na-1, nn, 1_000)

        if !isfinite(lr) || nn > 999
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

        # merge to current tip
        if def1(stree)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr =
        tip_sims!(tree.d1, t, α, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr =
        tip_sims!(tree.d2, t, α, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::iTfbd,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   ψ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

Continue simulation until time `t` for the fixed tip in `tree`.
"""
function fossiltip_sim!(tree::iTfbd,
                        t   ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        ψ   ::Vector{Float64},
                        ψts ::Vector{Float64},
                        ix  ::Int64,
                        nep ::Int64,
                        δt  ::Float64,
                        srδt::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        Iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64)


  if lU < lr && nn < 1_000
    if istip(tree)

      nep = lastindex(ψts) + 1
      stree, na, nn, lr =
        _sim_gbmfbd_it(t, lλ(tree)[end], lμ(tree)[end], α, σλ, σμ, ψ,
          ψts, ix, nep, δt, srδt, lr, lU, Iρi, na-1, nn, 1_000)

      if !isfinite(lr) || nn > 999
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr =
        fossiltip_sim!(tree.d1, t, α, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, Iρi, na, nn)
    else
      tree.d2, na, nn, lr =
        fossiltip_sim!(tree.d2, t, α, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_gbm!(bix  ::Int64,
                Ξ    ::Vector{iTfbd},
                idf  ::Vector{iBffs},
                α    ::Float64,
                σλ   ::Float64,
                σμ   ::Float64,
                llc  ::Float64,
                dλ   ::Float64,
                ssλ  ::Float64,
                ssμ  ::Float64,
                mc   ::Float64,
                th   ::Float64,
                surv   ::Int64,
                δt   ::Float64,
                srδt ::Float64,
                lλxpr::Float64,
                lμxpr::Float64)

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix  ::Int64,
                     Ξ    ::Vector{iTfbd},
                     idf  ::Vector{iBffs},
                     α    ::Float64,
                     σλ   ::Float64,
                     σμ   ::Float64,
                     llc  ::Float64,
                     dλ   ::Float64,
                     ssλ  ::Float64,
                     ssμ  ::Float64,
                     mc   ::Float64,
                     th   ::Float64,
                     surv ::Int64,
                     δt   ::Float64,
                     srδt ::Float64,
                     lλxpr::Float64,
                     lμxpr::Float64)
  @inbounds begin

    ξi   = Ξ[bix]
    bi   = idf[bix]
    i1   = d1(bi)
    i2   = d2(bi)
    ξ1   = Ξ[i1]
    root = iszero(pa(bi))

    if root && iszero(e(bi))

      # if stem fossil
      if isfossil(bi)
        llc, dλ, ssλ, ssμ, mc =
          _fstem_update!(ξi, ξ1, α, σλ, σμ, llc, dλ, ssλ, ssμ, mc, th, 
            δt, srδt, lλxpr, lμxpr, surv)
      # if crown
      else
        llc, dλ, ssλ, ssμ, mc =
          _crown_update!(ξi, ξ1, Ξ[i2], α, σλ, σμ, llc, dλ, ssλ, ssμ, mc, th,
            δt, srδt, lλxpr, lμxpr, surv)
        setλt!(bi, lλ(ξi)[1])
      end
    else
      # if stem
      if root
        llc, dλ, ssλ, ssμ, mc =
          _stem_update!(ξi, α, σλ, σμ, llc, dλ, ssλ, ssμ, mc, th, 
            δt, srδt, lλxpr, lμxpr, surv)
      end

      # updates within the parent branch
      llc, dλ, ssλ, ssμ =
        _update_gbm!(ξi, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, false)

      # get fixed tip
      lξi = fixtip(ξi)

      # if mid branch
      if iszero(i2)

        # make between decoupled trees node update
        llc, ssλ, ssμ =
          update_duo!(lλ(lξi), lλ(ξ1), lμ(lξi), lμ(ξ1), e(lξi), e(ξ1),
            fdt(lξi), fdt(ξ1), α, σλ, σμ, llc, ssλ, ssμ, δt, srδt)

      # if internal branch
      else
        ξ2 = Ξ[i2]
        # make between decoupled trees node update
        llc, dλ, ssλ, ssμ, λf =
          update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), lμ(lξi), lμ(ξ1), lμ(ξ2),
            e(lξi), e(ξ1), e(ξ2), fdt(lξi), fdt(ξ1), fdt(ξ2),
            α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt)

        # set fixed `λ(t)` in branch
        setλt!(bi, λf)
      end
    end

    # carry on updates in the daughters
    llc, dλ, ssλ, ssμ =
      _update_gbm!(ξ1, α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt,
        iszero(d1(idf[i1])))
    if i2 > 0
      llc, dλ, ssλ, ssμ =
        _update_gbm!(Ξ[i2], α, σλ, σμ, llc, dλ, ssλ, ssμ, δt, srδt, 
          iszero(d1(idf[i2])))
    end
  end

  return llc, dλ, ssλ, ssμ, mc
end





