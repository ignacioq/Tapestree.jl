#=

clads birth-death MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 28 07 2025
=#




"""
    insane_cladsfbd(tree    ::sTf_label;
                    λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                    μ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                    αλ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                    αμ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                    σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                    σμ_prior::NTuple{2,Float64}     = (3.0, 0.1),
                    ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    ψ_epoch ::Vector{Float64}       = Float64[],
                    f_epoch ::Vector{Int64}         = Int64[0],
                    niter   ::Int64                 = 1_000,
                    nthin   ::Int64                 = 10,
                    nburn   ::Int64                 = 200,
                    nflush  ::Int64                 = nthin,
                    ofile   ::String                = string(homedir(), "/cladsfbd"),
                    λi      ::Float64               = NaN,
                    μi      ::Float64               = NaN,
                    ϵi      ::Float64               = 0.2,
                    ψi      ::Float64               = NaN,
                    αλi     ::Float64               = 0.0,
                    αμi     ::Float64               = 0.0,
                    σλi     ::Float64               = 0.1,
                    σμi     ::Float64               = 0.1,
                    pupdp   ::NTuple{7,Float64}     = (1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 0.1, 0.2),
                    prints  ::Int64                 = 5,
                    stnλ    ::Float64               = 0.5,
                    stnμ    ::Float64               = 0.5,
                    survival::Bool                  = true,
                    mxthf   ::Float64               = 0.1,
                    tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for fossil clads.
"""
function insane_cladsfbd(tree    ::sTf_label;
                         λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                         μ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                         αλ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                         αμ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                         σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                         σμ_prior::NTuple{2,Float64}     = (3.0, 0.1),
                         ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                         ψ_epoch ::Vector{Float64}       = Float64[],
                         f_epoch ::Vector{Int64}         = Int64[0],
                         niter   ::Int64                 = 1_000,
                         nthin   ::Int64                 = 10,
                         nburn   ::Int64                 = 200,
                         nflush  ::Int64                 = nthin,
                         ofile   ::String                = string(homedir(), "/cladsfbd"),
                         λi      ::Float64               = NaN,
                         μi      ::Float64               = NaN,
                         ϵi      ::Float64               = 0.2,
                         ψi      ::Float64               = NaN,
                         αλi     ::Float64               = 0.0,
                         αμi     ::Float64               = 0.0,
                         σλi     ::Float64               = 0.1,
                         σμi     ::Float64               = 0.1,
                         pupdp   ::NTuple{7,Float64}     = (1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 0.1, 0.2),
                         prints  ::Int64                 = 5,
                         stnλ    ::Float64               = 0.5,
                         stnμ    ::Float64               = 0.5,
                         survival::Bool                  = true,
                         mxthf   ::Float64               = 0.1,
                         tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n  = ntips(tree)
  th = treeheight(tree)

  # turn to logarithmic terms
  λ0_prior = (log(λ0_prior[1]), 2.0*log(λ0_prior[2]))
  μ0_prior = (log(μ0_prior[1]), 2.0*log(μ0_prior[2]))

  # only include epochs where the tree occurs
  filter!(x -> x < th, ψ_epoch)
  sort!(ψ_epoch, rev = true)
  nep = lastindex(ψ_epoch) + 1

  # make initial fossils per epoch vector
  lep = lastindex(f_epoch)
  if lep !== nep
    if sum(f_epoch) > 0
      if lep > nep
        f_epoch = f_epoch[(end-nep+1):end]
      else 
        prepend!(f_epoch, zeros(Int64, nep-lep))
      end
    else
      f_epoch = zeros(Int64, nep)
    end
  end

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ, th * mxthf)

  # starting parameters
  λc, μc, ψc = λi, μi, ψi
  if isnan(λi) || isnan(μi) || isnan(ψi)
    # if only one tip
    if isone(n)
      λc = 1.0/th
      μc = 0.0
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
  end

  # make ψ vector
  ψc = fill(ψc, nep)

  # condition on first speciation event
  rmλ = iszero(e(tree)) && !isfossil(tree) ? 1.0 : 0.0

  surv = 0   # condition on survival of 0, 1, or 2 starting lineages
  if survival 
    if iszero(e(tree)) 
      if def1(tree)
        surv += Int64(anyalive(tree.d1))
        if def2(tree)
          surv += Int64(anyalive(tree.d2))
        end
      end
    else
      surv += Int64(anyalive(tree))
    end
  end

  # make a decoupled tree
  Ξ = make_Ξ(idf, λc, μc, cTfbd)

  # survival
  mc = m_surv_cladsfbd(th, log(λc), log(μc), αλi, αμi, σλi, σμi, 1_000, surv)

  # set end of fix branch speciation times and get vector of internal branches
  # and make epoch start vectors and indices for each `ξ`
  inodes = Int64[]
  eixi   = Int64[]
  eixf   = Int64[]
  bst    = Float64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    d1(bi) > 0 && push!(inodes, i)
    tib = ti(bi)
    ei  = findfirst(x -> x < tib, ψ_epoch)
    ei  = isnothing(ei) ? nep : ei
    ef  = findfirst(x -> x < tf(bi), ψ_epoch)
    ef  = isnothing(ef) ? nep : ef
    push!(bst, tib)
    push!(eixi, ei)
    push!(eixf, ef)
  end

  # parameter updates (1: αλ, 2: αμ, 3: σλ & σμ, 4: ψ, 5: scale, 6: internal, 7: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running cladogenetic fossilised birth-death"

  # burn-in phase
  Ξ, idf, llc, prc, αλc, αμc, σλc, σμc, mc, ns, ne, nf, 
  L, ddλ, ddμ, ssλ, ssμ, stnλ, stnμ =
    mcmc_burn_cladsfbd(Ξ, idf, 
      λ0_prior, μ0_prior, αλ_prior, αμ_prior, σλ_prior, σμ_prior, ψ_prior,
      ψ_epoch, f_epoch, nburn, αλi, αμi, σλi, σμi, ψc, mc, 
      th, rmλ, surv, stnλ, stnμ, bst, eixi, eixf, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cladsfbd(Ξ, idf, llc, prc, αλc, αμc, σλc, σμc, ψc, mc, th, rmλ, surv, 
      ns, ne, nf, L, ddλ, ddμ, ssλ, ssμ, stnλ, stnμ, λ0_prior, μ0_prior, 
      αλ_prior, αμ_prior, σλ_prior, σμ_prior, ψ_prior, ψ_epoch, f_epoch, 
      bst, eixi, eixf, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_cladsfbd(Ξ       ::Vector{cTfbd},
                       idf     ::Vector{iBffs},
                       λ0_prior::NTuple{2,Float64},
                       μ0_prior::NTuple{2,Float64},
                       αλ_prior::NTuple{2,Float64},
                       αμ_prior::NTuple{2,Float64},
                       σλ_prior::NTuple{2,Float64},
                       σμ_prior::NTuple{2,Float64},
                       ψ_prior ::NTuple{2,Float64},
                       ψ_epoch ::Vector{Float64},
                       f_epoch ::Vector{Int64},
                       nburn   ::Int64,
                       αλc     ::Float64,
                       αμc     ::Float64,
                       σλc     ::Float64,
                       σμc     ::Float64,
                       ψc      ::Vector{Float64},
                       mc      ::Float64,
                       th      ::Float64,
                       rmλ     ::Float64,
                       surv    ::Int64,
                       stnλ    ::Float64,
                       stnμ    ::Float64,
                       bst     ::Vector{Float64},
                       eixi    ::Vector{Int64},
                       eixf    ::Vector{Int64},
                       pup     ::Array{Int64,1},
                       prints  ::Int64)

MCMC burn-in chain for fossil clads.
"""
function mcmc_burn_cladsfbd(Ξ       ::Vector{cTfbd},
                            idf     ::Vector{iBffs},
                            λ0_prior::NTuple{2,Float64},
                            μ0_prior::NTuple{2,Float64},
                            αλ_prior::NTuple{2,Float64},
                            αμ_prior::NTuple{2,Float64},
                            σλ_prior::NTuple{2,Float64},
                            σμ_prior::NTuple{2,Float64},
                            ψ_prior ::NTuple{2,Float64},
                            ψ_epoch ::Vector{Float64},
                            f_epoch ::Vector{Int64},
                            nburn   ::Int64,
                            αλc     ::Float64,
                            αμc     ::Float64,
                            σλc     ::Float64,
                            σμc     ::Float64,
                            ψc      ::Vector{Float64},
                            mc      ::Float64,
                            th      ::Float64,
                            rmλ     ::Float64,
                            surv    ::Int64,
                            stnλ    ::Float64,
                            stnμ    ::Float64,
                            bst     ::Vector{Float64},
                            eixi    ::Vector{Int64},
                            eixf    ::Vector{Int64},
                            pup     ::Array{Int64,1},
                            prints  ::Int64)

  # starting likelihood and prior
  lλ0 = lλ(Ξ[1])
  llc = llik_clads(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi) - 
        rmλ*lλ0 + log(mc) + prob_ρ(idf)
  prc = logdnorm(lλ0,       λ0_prior[1], λ0_prior[2])   +
        logdnorm(lμ(Ξ[1]),  μ0_prior[1], μ0_prior[2])   +
        logdnorm(αλc,       αλ_prior[1], αλ_prior[2]^2) +
        logdnorm(αμc,       αμ_prior[1], αμ_prior[2]^2) +
        logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])   +
        logdinvgamma(σμc^2, σμ_prior[1], σμ_prior[2])

  L   = treelength(Ξ, ψ_epoch, bst, eixi)        # tree length
  el  = lastindex(idf)                           # number of branches
  ns  = sum(x -> Float64(d2(x) > 0), idf) - rmλ  # number of speciation events in likelihood
  ne  = Float64(ntipsextinct(Ξ))                 # number of extinction events in likelihood
  nf  = nfossils(idf, ψ_epoch, f_epoch)         # number of fossilization events per epoch
  nep = lastindex(ψc)                           # number of epochs
  λfs = Float64[]
  μfs = Float64[]

  # delta change, sum squares, path length and integrated rate
  ddλ, ddμ, ssλ, ssμ = _dd_ss(Ξ, idf, αλc, αμc)

  # for scale tuning
  ltn = zero(Int64)
  lup = lacλ = lacμ = zero(Float64)

  pbar = Progress(nburn, dt = prints, desc = "burning mcmc...", barlen = 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    ## parameter updates
    for pupi in pup

      # update speciation drift
      if pupi === 1

        llc, prc, αλc, mc = 
          update_αλ!(αλc, lλ(Ξ[1]), lμ(Ξ[1]), αμc, σλc, σμc, 
            2.0*(ns + rmλ), ddλ, llc, prc, mc, th, surv, αλ_prior)

        # update ssλ with new drift `αλ`
        ssλ = _ss(Ξ, idf, αλc, lλ, λt)

      # update extinction drift
      elseif pupi === 2

        llc, prc, αμc, mc = 
          update_αμ!(αμc, lλ(Ξ[1]), lμ(Ξ[1]), αλc, σλc, σμc, 
            2.0*(ns + rmλ), ddμ, llc, prc, mc, th, surv, αμ_prior)

       # update ssμ with new drift `αμ`
        ssμ = _ss(Ξ, idf, αμc, lμ, μt)

      # update speciation and extinction diffusion rate
      elseif pupi === 3

        llc, prc, σλc, σμc, mc =
          update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αλc, αμc, ssλ, ssμ, 
            2.0*(ns + rmλ), llc, prc, mc, th, surv, σλ_prior, σμ_prior)

      # update fossilization rate
      elseif pupi === 4

        llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # update scale
      elseif pupi === 5

        llc, prc, mc, accλ, accμ = 
          update_scale!(Ξ, idf, αλc, αμc, σλc, σμc, llc, prc, ns, ne, 
            stnλ, stnμ, mc, th, surv, λ0_prior, μ0_prior)

        lacλ += accλ
        lacμ += accμ
        lup  += 1.0

      # update internal
      elseif pupi === 6

        bix = fIrand(el) + 1

        llc, prc, ddλ, ddμ, ssλ, ssμ, mc =
          update_internal!(bix, Ξ, idf, αλc, αμc, σλc, σμc, llc, prc, 
            ddλ, ddμ, ssλ, ssμ, mc, th, λ0_prior, μ0_prior, surv)

      # forward simulation
      else

        bix = fIrand(el) + 1

        llc, ddλ, ddμ, ssλ, ssμ, ns, ne =
          update_fs!(bix, Ξ, idf, αλc, αμc, σλc, σμc, ψc, llc, L, 
            ddλ, ddμ, ssλ, ssμ, ns, ne, ψ_epoch, eixi, eixf, λfs, μfs)
      end
    end

    ltn += 1
    if ltn === 100
      stnλ = min(2.0, tune(stnλ, lacλ/lup))
      stnμ = min(2.0, tune(stnμ, lacμ/lup))
      ltn = zero(Int64)
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αλc, αμc, σλc, σμc, mc, ns, ne, nf,
         L, ddλ, ddμ, ssλ, ssμ, stnλ, stnμ
end




"""
    mcmc_cladsfbd(Ξ       ::Vector{cTfbd},
                  idf     ::Vector{iBffs},
                  llc     ::Float64,
                  prc     ::Float64,
                  αλc     ::Float64,
                  αμc     ::Float64,
                  σλc     ::Float64,
                  σμc     ::Float64,
                  ψc      ::Vector{Float64},
                  mc      ::Float64,
                  th      ::Float64,
                  rmλ     ::Float64,
                  surv    ::Int64,
                  ns      ::Float64,
                  ne      ::Float64,
                  nf      ::Vector{Float64},
                  L       ::Vector{Float64},
                  ddλ     ::Float64,
                  ddμ     ::Float64,
                  ssλ     ::Float64,
                  ssμ     ::Float64,
                  stnλ    ::Float64,
                  stnμ    ::Float64,
                  λ0_prior::NTuple{2,Float64},
                  μ0_prior::NTuple{2,Float64},
                  αλ_prior::NTuple{2,Float64},
                  αμ_prior::NTuple{2,Float64},
                  σλ_prior::NTuple{2,Float64},
                  σμ_prior::NTuple{2,Float64},
                  ψ_prior ::NTuple{2,Float64},
                  ψ_epoch ::Vector{Float64},
                  f_epoch ::Vector{Int64},
                  bst     ::Vector{Float64},
                  eixi    ::Vector{Int64},
                  eixf    ::Vector{Int64},
                  pup     ::Vector{Int64},
                  niter   ::Int64,
                  nthin   ::Int64,
                  nflush  ::Int64,
                  ofile   ::String,
                  prints  ::Int64)

MCMC chain for fossil clads.
"""
function mcmc_cladsfbd(Ξ       ::Vector{cTfbd},
                       idf     ::Vector{iBffs},
                       llc     ::Float64,
                       prc     ::Float64,
                       αλc     ::Float64,
                       αμc     ::Float64,
                       σλc     ::Float64,
                       σμc     ::Float64,
                       ψc      ::Vector{Float64},
                       mc      ::Float64,
                       th      ::Float64,
                       rmλ     ::Float64,
                       surv    ::Int64,
                       ns      ::Float64,
                       ne      ::Float64,
                       nf      ::Vector{Float64},
                       L       ::Vector{Float64},
                       ddλ     ::Float64,
                       ddμ     ::Float64,
                       ssλ     ::Float64,
                       ssμ     ::Float64,
                       stnλ    ::Float64,
                       stnμ    ::Float64,
                       λ0_prior::NTuple{2,Float64},
                       μ0_prior::NTuple{2,Float64},
                       αλ_prior::NTuple{2,Float64},
                       αμ_prior::NTuple{2,Float64},
                       σλ_prior::NTuple{2,Float64},
                       σμ_prior::NTuple{2,Float64},
                       ψ_prior ::NTuple{2,Float64},
                       ψ_epoch ::Vector{Float64},
                       f_epoch ::Vector{Int64},
                       bst     ::Vector{Float64},
                       eixi    ::Vector{Int64},
                       eixf    ::Vector{Int64},
                       pup     ::Vector{Int64},
                       niter   ::Int64,
                       nthin   ::Int64,
                       nflush  ::Int64,
                       ofile   ::String,
                       prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin = lit = sthin = zero(Int64)

  el    = lastindex(idf)   # number of branches
  nep   = lastindex(ψc)    # number of epochs
  λfs   = Float64[]
  μfs   = Float64[]

  # parameter results
  r     = Array{Float64,2}(undef, nlogs, 9 + nep)
  treev = cTfbd[]          # make tree log vector
  io    = IOBuffer()       # buffer 

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\tmu_root\talpha_lambda\talpha_mu\tsigma_lambda\tsigma_mu\t"*join(["psi"*(isone(nep) ? "" : string("_",i)) for i in 1:nep], '\t')*'\n')
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αλc = αλc, αμc = αμc, σλc = σλc, σμc = σμc, mc = mc, ns = ns, ne = ne, L = L, ddλ = ddλ, ddμ = ddμ, ssλ = ssλ, ssμ = ssμ, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update drift
            if pupi === 1

              llc, prc, αλc, mc = 
                update_αλ!(αλc, lλ(Ξ[1]), lμ(Ξ[1]), αμc, σλc, σμc, 
                  2.0*(ns + rmλ), ddλ, llc, prc, mc, th, surv, αλ_prior)

              # update ssλ with new drift `αλ`
              ssλ = _ss(Ξ, idf, αλc, lλ, λt)

              # ll0 = llik_clads(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #   @show ll0, llc, it, pupi
              #   return
              # end

            # update extinction drift
            elseif pupi === 2

              llc, prc, αμc, mc = 
                update_αμ!(αμc, lλ(Ξ[1]), lμ(Ξ[1]), αλc, σλc, σμc, 
                  2.0*(ns + rmλ), ddμ, llc, prc, mc, th, surv, αμ_prior)

             # update ssμ with new drift `αμ`
              ssμ = _ss(Ξ, idf, αμc, lμ, μt)

              # ll0 = llik_clads(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #   @show ll0, llc, it, pupi
              #   return
              # end

            # update speciation and extinction diffusion rate
            elseif pupi === 3

              llc, prc, σλc, σμc, mc =
                update_σ!(σλc, σμc, lλ(Ξ[1])[1], lμ(Ξ[1])[1], αλc, αμc, ssλ, ssμ, 
                  2.0*(ns + rmλ), llc, prc, mc, th, surv, σλ_prior, σμ_prior)

              # ll0 = llik_clads(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #   @show ll0, llc, it, pupi
              #   return
              # end

            # update fossilization rate
            elseif pupi === 4

              llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

              # ll0 = llik_clads(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #   @show ll0, llc, it, pupi
              #   return
              # end

            # update scale
            elseif pupi === 5

              llc, prc, mc, accλ, accμ = 
                update_scale!(Ξ, idf, αλc, αμc, σλc, σμc, llc, prc, ns, ne, 
                  stnλ, stnμ, mc, th, surv, λ0_prior, μ0_prior)

              # ll0 = llik_clads(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #   @show ll0, llc, it, pupi
              #   return
              # end

            # update internal λ & μ
            elseif pupi === 6

              bix = fIrand(el) + 1

              llc, prc, ddλ, ddμ, ssλ, ssμ, mc =
                update_internal!(bix, Ξ, idf, αλc, αμc, σλc, σμc, llc, prc, 
                  ddλ, ddμ, ssλ, ssμ, mc, th, λ0_prior, μ0_prior, surv)

              # ll0 = llik_clads(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #   @show ll0, llc, it, pupi
              #   return
              # end

            # update by forward simulation
            else

              bix = fIrand(el) + 1

              llc, ddλ, ddμ, ssλ, ssμ, ns, ne =
                update_fs!(bix, Ξ, idf, αλc, αμc, σλc, σμc, ψc, llc, L, 
                  ddλ, ddμ, ssλ, ssμ, ns, ne, ψ_epoch, eixi, eixf, λfs, μfs)

              # ll0 = llik_clads(Ξ, idf, αλc, αμc, σλc, σμc, ψc, ψ_epoch, bst, eixi) - rmλ*lλ(Ξ[1]) + log(mc) + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #   @show ll0, llc, it, pupi
              #   return
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
              r[lit,5] = exp(lμ(Ξ[1]))
              r[lit,6] = αλc
              r[lit,7] = αμc
              r[lit,8] = σλc
              r[lit,9] = σμc
              @turbo for i in Base.OneTo(nep)
                r[lit, 9 + i] = ψc[i]
              end
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
                  exp(lλ(Ξ[1])),'\t', exp(lμ(Ξ[1])),'\t', αλc, '\t', αμc, '\t', 
                  σλc, '\t', σμc, '\t', join(ψc, '\t'), '\n')
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
    update_αλ!(αλc     ::Float64,
               λ0      ::Float64,
               μ0      ::Float64,
               αμ      ::Float64,
               σλ      ::Float64,
               σμ      ::Float64,
               L       ::Float64,
               ddλ     ::Float64,
               llc     ::Float64,
               prc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               surv    ::Int64,
               αλ_prior::NTuple{2,Float64})

Gibbs update for `αλ`.
"""
function update_αλ!(αλc     ::Float64,
                    λ0      ::Float64,
                    μ0      ::Float64,
                    αμ      ::Float64,
                    σλ      ::Float64,
                    σμ      ::Float64,
                    L       ::Float64,
                    ddλ     ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    surv    ::Int64,
                    αλ_prior::NTuple{2,Float64})
  ν   = αλ_prior[1]
  τ2  = αλ_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2
  αλp = rnorm((ddλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  mp  = m_surv_cladsfbd(th, λ0, μ0, αλp, αμ, σλ, σμ, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σλ2*(αλc^2 - αλp^2 + 2.0*ddλ*(αλp - αλc)/L) + llr
    prc += llrdnorm_x(αλp, αλc, ν, τ2)
    αλc  = αλp
    mc   = mp
  end

  return llc, prc, αλc, mc
end




"""
    update_αμ!(αμc     ::Float64,
               λ0      ::Float64,
               μ0      ::Float64,
               αλ      ::Float64,
               σλ      ::Float64,
               σμ      ::Float64,
               L       ::Float64,
               ddμ     ::Float64,
               llc     ::Float64,
               prc     ::Float64,
               mc      ::Float64,
               th      ::Float64,
               surv    ::Int64,
               αμ_prior::NTuple{2,Float64})

Gibbs update for `αμ`.
"""
function update_αμ!(αμc     ::Float64,
                    λ0      ::Float64,
                    μ0      ::Float64,
                    αλ      ::Float64,
                    σλ      ::Float64,
                    σμ      ::Float64,
                    L       ::Float64,
                    ddμ     ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    mc      ::Float64,
                    th      ::Float64,
                    surv    ::Int64,
                    αμ_prior::NTuple{2,Float64})

  ν   = αμ_prior[1]
  τ2  = αμ_prior[2]^2
  σμ2 = σμ^2
  rs  = σμ2/τ2
  αμp = rnorm((ddμ + rs*ν)/(rs + L), sqrt(σμ2/(rs + L)))

  mp  = m_surv_cladsfbd(th, λ0, μ0, αλ, αμp, σλ, σμ, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += 0.5*L/σμ2*(αμc^2 - αμp^2 + 2.0*ddμ*(αμp - αμc)/L) + llr
    prc += llrdnorm_x(αμp, αμc, ν, τ2)
    αμc  = αμp
    mc   = mp
  end

  return llc, prc, αμc, mc
end




"""
    update_σ!(σλc     ::Float64,
              σμc     ::Float64,
              λ0      ::Float64,
              μ0      ::Float64,
              αλ      ::Float64,
              αμ      ::Float64,
              ssλ     ::Float64,
              ssμ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              mc      ::Float64,
              th      ::Float64,
              surv    ::Bool,
              σλ_prior::NTuple{2,Float64},
              σμ_prior::NTuple{2,Float64})

Gibbs update for `σλ` and `σμ`.
"""
function update_σ!(σλc     ::Float64,
                   σμc     ::Float64,
                   λ0      ::Float64,
                   μ0      ::Float64,
                   αλ      ::Float64,
                   αμ      ::Float64,
                   ssλ     ::Float64,
                   ssμ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   surv    ::Int64,
                   σλ_prior::NTuple{2,Float64},
                   σμ_prior::NTuple{2,Float64})

  # Gibbs update for σ
  σλp2 = rand(InverseGamma(σλ_prior[1] + 0.5 * n, σλ_prior[2] + ssλ))
  σμp2 = rand(InverseGamma(σμ_prior[1] + 0.5 * n, σμ_prior[2] + ssμ))
  σλp  = sqrt(σλp2)
  σμp  = sqrt(σμp2)

  mp  = m_surv_cladsfbd(th, λ0, μ0, αλ, αμ, σλp, σμp, 1_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += ssλ*(1.0/σλc^2 - 1.0/σλp2) - n*(log(σλp/σλc)) +
           ssμ*(1.0/σμc^2 - 1.0/σμp2) - n*(log(σμp/σμc)) +
           llr
    prc += llrdinvgamma(σλp2, σλc^2, σλ_prior[1], σλ_prior[2]) +
           llrdinvgamma(σμp2, σμc^2, σμ_prior[1], σμ_prior[2])
    σλc  = σλp
    σμc  = σμp
    mc   = mp
  end

  return llc, prc, σλc, σμc, mc
end




"""
    update_scale!(Ξ       ::Vector{cTfbd},
                  idf     ::Vector{iBffs},
                  αλ      ::Float64,
                  αμ      ::Float64,
                  σλ      ::Float64,
                  σμ      ::Float64,
                  llc     ::Float64,
                  prc     ::Float64,
                  ns      ::Float64,
                  ne      ::Float64,
                  stnλ    ::Float64,
                  stnμ    ::Float64,
                  mc      ::Float64,
                  th      ::Float64,
                  surv    ::Int64,
                  λ0_prior::NTuple{2,Float64},
                  μ0_prior::NTuple{2,Float64})

Update scale for speciation and extinction.
"""
function update_scale!(Ξ       ::Vector{cTfbd},
                       idf     ::Vector{iBffs},
                       αλ      ::Float64,
                       αμ      ::Float64,
                       σλ      ::Float64,
                       σμ      ::Float64,
                       llc     ::Float64,
                       prc     ::Float64,
                       ns      ::Float64,
                       ne      ::Float64,
                       stnλ    ::Float64,
                       stnμ    ::Float64,
                       mc      ::Float64,
                       th      ::Float64,
                       surv    ::Int64,
                       λ0_prior::NTuple{2,Float64},
                       μ0_prior::NTuple{2,Float64})

  # extract integrated rates
  irλ, irμ = _irbd(Ξ)

  accλ = accμ = 0.0

  lλ0, lμ0 = lλ(Ξ[1]), lμ(Ξ[1])

  ## start with speciation
  # sample log(scaling factor)
  s = randn()*stnλ

  # likelihood ratio
  llr = ns * s + (1.0 - exp(s)) * irλ

  # prior ratio
  prr = llrdnorm_x(lλ0 + s, lλ0, λ0_prior[1], λ0_prior[2]) 

  lU = -randexp()

  if lU < llr + prr + log(1000.0/mc)

    # survival
    mp   = m_surv_cladsfbd(th, lλ0 + s, lμ0, αλ, αμ, σλ, σμ, 1_000, surv)
    llr += log(mp/mc)

    if -randexp() < llr + prr
      accλ += 1.0
      llc  += llr
      prc  += prr
      mc   = mp
      scale_rate!(Ξ,   addlλ!, s)
      scale_rate!(idf, addlλ!, s)
      lλ0  += s
    end
  end

  ## continue with extinction
  # sample log(scaling factor)
  s = randn()*stnμ

  # likelihood ratio
  llr = ne * s + (1.0 - exp(s)) * irμ

  # prior ratio
  prr = llrdnorm_x(lμ0 + s, lμ0, μ0_prior[1], μ0_prior[2]) 

  lU = -randexp()

  if lU < llr + prr + log(1000.0/mc)

    # add survival ratio
    mp   = m_surv_cladsfbd(th, lλ0, lμ0 + s, αλ, αμ, σλ, σμ, 1_000, surv)
    llr += log(mp/mc)

    if lU < llr + prr
      accμ += 1.0
      llc  += llr
      prc  += prr
      mc    = mp
      scale_rate!(Ξ,   addlμ!, s)
      scale_rate!(idf, addlμ!, s)
    end
  end

  return llc, prc, mc, accλ, accμ
end




"""
    update_internal!(bix     ::Int64,
                     Ξ       ::Vector{cTfbd},
                     idf     ::Vector{iBffs},
                     αλ      ::Float64,
                     αμ      ::Float64,
                     σλ      ::Float64,
                     σμ      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ddλ     ::Float64,
                     ssλ     ::Float64,
                     ssμ     ::Float64,
                     mc      ::Float64,
                     th      ::Float64,
                     λ0_prior::NTuple{2,Float64},
                     μ0_prior::NTuple{2,Float64},
                     surv    ::Int64)

Make an update for an internal branch and its descendants.
"""
function update_internal!(bix     ::Int64,
                          Ξ       ::Vector{cTfbd},
                          idf     ::Vector{iBffs},
                          αλ      ::Float64,
                          αμ      ::Float64,
                          σλ      ::Float64,
                          σμ      ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          ddλ     ::Float64,
                          ddμ     ::Float64,
                          ssλ     ::Float64,
                          ssμ     ::Float64,
                          mc      ::Float64,
                          th      ::Float64,
                          λ0_prior::NTuple{2,Float64},
                          μ0_prior::NTuple{2,Float64},
                          surv    ::Int64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  it   = iszero(i1) # is terminal
  i2   = d2(bi)
  ia   = pa(bi)
  root = iszero(ia)
  λa = μa = NaN  # ancestral speciation and extinction

  # if crown root
  if root && iszero(e(ξi)) && !isfossil(bi)
    llc, prc, ddλ, ddμ, ssλ, ssμ, mc =
      _crown_update!(ξi, Ξ[i1], Ξ[i2], αλ, αμ, σλ, σμ, llc, prc, ddλ, ddμ, 
        ssλ, ssμ, mc, th, λ0_prior, μ0_prior, surv)
    λa, μa = lλ(ξi), lμ(ξi)
    setλt!(bi, λa)
    setμt!(bi, μa)
  else
    # if stem
    if root
      # if cladogenetic branch
      if i2 > 0
        ξ1, ξ2 = Ξ[i1], Ξ[i2]
        eds, λ1, λ2, μ1, μ2 = 0.0, lλ(ξ1), lλ(ξ2), lμ(ξ1), lμ(ξ2)
      # if fossil or mid branch
      elseif i1 > 0 || isfossil(bi)
        eds, λ1, λ2, μ1, μ2 = 
          downstreamλμs(bix, Ξ, idf, 0.0, NaN, NaN, NaN, NaN)
      end

      llc, prc, ddλ, ddμ, ssλ, ssμ, mc, λi, μi = 
        _stem_update!(ξi, eds, λ1, λ2, μ1, μ2, αλ, αμ, σλ, σμ, llc, prc, 
          ddλ, ddμ, ssλ, ssμ, mc, th, λ0_prior, μ0_prior, surv)

      # set new λ & μ downstream, if necessary
      setdownstreamλμ!(λi, μi, bix, Ξ, idf)

      # if there are speciation events in stem branch
      if !istip(ξi)
        eds, λ1, λ2, μ1, μ2 = 
          downstreamλμs(bix, Ξ, idf, 0.0, NaN, NaN, NaN, NaN)

        # updates within the parent branch
        llc, ddλ, ddμ, ssλ, ssμ, λx, μx = 
          _update_internal!(ξi.d1, bi, eas, λa, μa, αλ, αμ, σλ, σμ, 
            eds, λ1, λ2, μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ, false)
        llc, ddλ, ddμ, ssλ, ssμ, λx, μx = 
          _update_internal!(ξi.d2, bi, eas, λa, μa, αλ, αμ, σλ, σμ, 
            eds, λ1, λ2, μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ, false)
      end

    # if *not* root
    else

      # find cladogenetic ancestor
      eas, λa, μa, il = upstreamλμ(ia, Ξ, idf, 0.0, λa, μa)

      # check if mid branch does not lead to the stem branch
      if pa(idf[il]) > 0

        eds, λ1, λ2, μ1, μ2 = 0.0, NaN, NaN, NaN, NaN

        # if cladogenetic branch
        if i2 > 0
          ξ1, ξ2 = Ξ[i1], Ξ[i2]
          eds, λ1, λ2, μ1, μ2 = 0.0, lλ(ξ1), lλ(ξ2), lμ(ξ1), lμ(ξ2)
        # if mid or fossil branch
        elseif i1 > 0
          eds, λ1, λ2, μ1, μ2 = 
            downstreamλμs(i1, Ξ, idf, 0.0, NaN, NaN, NaN, NaN)
        end

        ll0 = llc

        # updates within the parent branch
        llc, ddλ, ddμ, ssλ, ssμ, λx, μx = 
          _update_internal!(ξi, bi, eas, λa, μa, αλ, αμ, σλ, σμ, 
            eds, λ1, λ2, μ1, μ2, llc, ddλ, ddμ, ssλ, ssμ, it)

        # if update, update up- and down-stream
        if ll0 != llc
          # if fossil tip
          if it && isfossil(bi)
            lξi = fixtip(ξi)
            λi, μi = lλ(lξi.d1), lμ(lξi.d1)
            setupstreamλμ!(λi, μi, bix, Ξ, idf)
          end

          λi, μi = lλ(ξi), lμ(ξi)
          setupstreamλμ!(λi, μi, ia, Ξ, idf)

          if !it && iszero(i2)
            lξi = fixtip(ξi)
            setdownstreamλμ!(lλ(lξi), lμ(lξi), i1, Ξ, idf)
          end
        end
      end
    end
  end

  return llc, prc, ddλ, ddμ, ssλ, ssμ, mc
end




"""
    update_fs!(bix ::Int64,
               Ξ   ::Vector{cTfbd},
               idf ::Vector{iBffs},
               αλ  ::Float64,
               αμ  ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               ψ   ::Vector{Float64},
               llc ::Float64,
               L   ::Vector{Float64},
               ddλ ::Float64,
               ddμ ::Float64,
               ssλ ::Float64,
               ssμ ::Float64,
               ns  ::Float64,
               ne  ::Float64,
               ψts ::Vector{Float64},
               eixi::Vector{Int64},
               eixf::Vector{Int64},
               λfs ::Vector{Float64},
               μfs ::Vector{Float64})

Forward simulation proposal for clads.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{cTfbd},
                    idf ::Vector{iBffs},
                    αλ  ::Float64,
                    αμ  ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Vector{Float64},
                    llc ::Float64,
                    L   ::Vector{Float64},
                    ddλ ::Float64,
                    ddμ ::Float64,
                    ssλ ::Float64,
                    ssμ ::Float64,
                    ns  ::Float64,
                    ne  ::Float64,
                    ψts ::Vector{Float64},
                    eixi::Vector{Int64},
                    eixf::Vector{Int64},
                    λfs ::Vector{Float64},
                    μfs ::Vector{Float64})

  bi  = idf[bix]
  ξc  = Ξ[bix]
  ixi = eixi[bix]
  ia  = pa(bi)

  λa = μa = NaN
  # if following a speciation event
  if ia > 0 
    ba = idf[ia]
    if d2(ba) > 0
      λa, μa = λt(ba), μt(ba)
    end
  end

  ddλr = ddμr = ssλr = ssμr = zero(Float64)
  llr  = NaN

  # terminal branch
  if iszero(d1(bi))

    # if fossil terminal branch
    if isfossil(bi)
      ixf = eixf[bix]

      ξp, llr = fsbi_t(bi, ξc, λa, μa, αλ, αμ, σλ, σμ, ψ, ψts, ixi, ixf)

      # if terminal but not successful proposal, update extinct
      if !isfinite(llr)
         ξp, llr = 
          fsbi_et(cTfbd_wofe(ξc), bi, αλ, αμ, σλ, σμ, ψ, ψts, ixf)
      end

    # non-fossil terminal branch
    else
      ξp, llr = fsbi_t(bi, ξc, λa, μa, αλ, αμ, σλ, σμ, ψ, ψts, ixi)
    end

  # internal non-bifurcating branch
  elseif iszero(d2(bi))

    ξp, llr, ddλr, ddμr, ssλr, ssμr =
      fsbi_m(bi, idf, ξc, Ξ, λa, μa, αλ, αμ, σλ, σμ, ψ, ψts, 
        ixi, eixf[bix], λfs, μfs)

  # internal bifurcating branch
  else

    ξ1, ξ2 = Ξ[d1(bi)], Ξ[d2(bi)]
    ξp, llr, ddλr, ddμr, ssλr, ssμr =
      fsbi_i(bi, ξc, λa, lλ(ξ1), lλ(ξ2), μa, lμ(ξ1), lμ(ξ2), 
        αλ, αμ, σλ, σμ, ψ, ψts, ixi, eixf[bix], λfs, μfs)
  end

  # if accepted
  if isfinite(llr)
    tii = ti(bi)
    nep = lastindex(ψts) + 1

    llc, ddλ, ddμ, ssλ, ssμ, ns, ne = 
      llik_cladsfbd_track!(ξc, αλ, αμ, σλ, σμ, ψ, tii, ψts, ixi, nep, 
        llc, L, ddλ, ddμ, ssλ, ssμ, ns, ne, -)
    llc, ddλ, ddμ, ssλ, ssμ, ns, ne = 
      llik_cladsfbd_track!(ξp, αλ, αμ, σλ, σμ, ψ, tii, ψts, ixi, nep, 
        llc, L, ddλ, ddμ, ssλ, ssμ, ns, ne, +)

    # first change from ancestor
    if isfinite(λa)
      λp, λc, = lλ(ξp), lλ(ξc)
      μp, μc, = lμ(ξp), lμ(ξc)
      llc += llrdnorm_x(λp, λc, λa + αλ, σλ^2) + 
             llrdnorm_x(μp, μc, μa + αμ, σμ^2)
      ddλ += λp - λc
      ddμ += μp - μc
      ssλ += 0.5*((λp - λa - αλ)^2 - (λc - λa - αλ)^2)
      ssμ += 0.5*((μp - μa - αμ)^2 - (μc - μa - αμ)^2)
    end

    # update quantities
    ddλ += ddλr
    ddμ += ddμr
    ssλ += ssλr
    ssμ += ssμr
    llc += llr

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ddμ, ssλ, ssμ, ns, ne
end




"""
    fsbi_t(bi ::iBffs,
           ξi ::cTfbd,
           λa ::Float64,
           μa ::Float64,
           αλ ::Float64,
           αμ ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           ψ  ::Vector{Float64},
           ψts::Vector{Float64},
           ix ::Int64)

Forward simulation for terminal branch `bi`.
"""
function fsbi_t(bi ::iBffs,
                ξi ::cTfbd,
                λa ::Float64,
                μa ::Float64,
                αλ ::Float64,
                αμ ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                ψ  ::Vector{Float64},
                ψts::Vector{Float64},
                ix ::Int64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # if does **not** come from a cladogenetic event
  if isnan(λa)
    λi, μi = lλ(ξi), lμ(ξi)
  else
    λi = rnorm(λa + αλ, σλ)
    μi = rnorm(μa + αμ, σμ)
  end

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nn, llr =
    _sim_cladsfbd_t(e(bi), λi, μi, αλ, αμ, σλ, σμ, ψ, ψts, 
      ix, nep, lc, lU, iρi, 0, 1, 500)

  if na > 0 && isfinite(llr)
    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  end

  return t0, NaN
end




"""
     fsbi_t(bi  ::iBffs,
            ξi  ::cTfbd,
            λa ::Float64,
            μa ::Float64,
            αλ  ::Float64,
            αμ  ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            ψ   ::Vector{Float64},
            ψts ::Vector{Float64},
            ixi ::Int64,
            ixf ::Int64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_t(bi ::iBffs,
                ξi ::cTfbd,
                λa ::Float64,
                μa ::Float64,
                αλ ::Float64,
                αμ ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                ψ  ::Vector{Float64},
                ψts::Vector{Float64},
                ixi::Int64,
                ixf::Int64)

  # if does **not** come from a cladogenetic event
  if isnan(λa)
    λi, μi = lλ(ξi), lμ(ξi)
  else
    λi = rnorm(λa + αλ, σλ)
    μi = rnorm(μa + αμ, σμ)
  end

   # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, af, nn =
    _sim_cladsfbd_i(ti(bi), tf(bi), λi, μi, αλ, αμ, σλ, σμ, ψ, ψts, 
      ixi, nep, 0, false, 1, 500)

  if na < 1 || af || nn > 499
    return t0, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  if lU < acr

    _fixrtip!(t0, na) # fix random tip

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, 
          ixf, nep, acr, lU, iρi, na, nn)
    end

    if lU < acr

      # fossilize extant tip
      fossilizefixedtip!(t0)

      # if terminal fossil branch
      tx, na, nn, acr =
        fossiltip_sim!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, ixf, nep,
          acr, lU, iρi, na, nn)

      if lU < acr

        llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
        setnt!(bi, ntp)      # set new nt
        setni!(bi, na)       # set new ni

        return t0, llr
      end
    end
  end

  return t0, NaN
end




"""
    fsbi_et(t0  ::cTfbd,
            bi  ::iBffs,
            αλ  ::Float64,
            αμ  ::Float64,
            σλ  ::Float64,
            σμ  ::Float64,
            ψ   ::Vector{Float64},
            ψts ::Vector{Float64},
            ixf ::Int64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_et(t0  ::cTfbd,
                 bi  ::iBffs,
                 αλ  ::Float64,
                 αμ  ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 ψ   ::Vector{Float64},
                 ψts ::Vector{Float64},
                 ixf ::Int64)

  nep = lastindex(ψts) + 1
  lU  = -randexp()            # log-probability
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # if terminal fossil branch
  tx, na, nn, acr =
    fossiltip_sim!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, ixf, nep,
      acr, lU, iρi, 1, 1)

  if lU < acr

    llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
    setni!(bi, na)       # set new ni

    return t0, llr
  end

  return t0, NaN
end




"""
    fsbi_m(bi ::iBffs,
           idf::Vector{iBffs},
           ξi ::cTfbd,
           Ξ  ::Vector{cTfbd},
           λa ::Float64,
           μa ::Float64,
           αλ  ::Float64,
           αμ  ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           ψ  ::Vector{Float64},
           ψts::Vector{Float64},
           ixi::Int64,
           ixf::Int64,
           λfs::Vector{Float64},
           μfs::Vector{Float64})

Forward simulation for mid branch `bi`.
"""
function fsbi_m(bi ::iBffs,
                idf::Vector{iBffs},
                ξi ::cTfbd,
                Ξ  ::Vector{cTfbd},
                λa ::Float64,
                μa ::Float64,
                αλ  ::Float64,
                αμ  ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                ψ  ::Vector{Float64},
                ψts::Vector{Float64},
                ixi::Int64,
                ixf::Int64,
                λfs::Vector{Float64},
                μfs::Vector{Float64})

  # if does **not** come from a cladogenetic event
  if isnan(λa)
    λi, μi = lλ(ξi), lμ(ξi)
  else
    λi = rnorm(λa + αλ, σλ)
    μi = rnorm(μa + αμ, σμ)
  end

  # forward simulation during branch length
  empty!(λfs)
  empty!(μfs)
  nep = lastindex(ψts) + 1
  t0, na, af, nn = _sim_cladsfbd_i(ti(bi), tf(bi), λi, μi, αλ, αμ, σλ, σμ, ψ, 
    ψts, ixi, nep, 0, false, 1, 500, λfs, μfs)

  if na < 1 || af || nn > 499
    return t0, NaN, NaN, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr  = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # search for next lλ1 and lλ2 if the exist
  i1 = d1(bi)
  eds, λ1, λ2, μ1, μ2 = downstreamλμs(i1, Ξ, idf, 0.0, NaN, NaN, NaN, NaN)

  ## choose most likely lineage to fix
  ddλr = ddμr = ssλr = ssμr = 0.0

  # if downstream is alive tip
  if isnan(λ1)
    # if downstream is extinct tip
    if isfinite(μ1)
      wt, λp, μp, pp, λc, μc, pc, acr = 
        wfix_m(ξi, e(bi), isfinite(μ1), λfs, μfs, eds, acr)
    else
      wt, λp, μp, pp, λc, μc, pc, acr = 
        wfix_m(ξi, e(bi), λfs, μfs, eds, acr)
    end
  # if downstream is cladogenetic
  else
    wt, λp, μp, pp, λc, μc, pc, acr, ddλr, ddμr, ssλr, ssμr = 
      wfix_m(ξi, e(bi), λfs, μfs, eds, λ1, λ2, μ1, μ2, αλ, αμ, σλ, σμ, acr)
  end

  if lU < acr

    # fix the tip
    if wt <= div(na,2)
      fixtip1!(t0, wt, 0)
    else
      fixtip2!(t0, na - wt + 1, 0)
    end

    # simulated remaining tips until the present
    tx, na, nn, acr =
      tip_sims!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, 
        ixf, nep, acr, lU, iρi, na, nn)

    if lU < acr
      na -= 1

      # fossilize 
      isfossil(bi) && fossilizefixedtip!(t0)

      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc)
      if isfinite(λ1)
        llr += λp - λc
      end
      setni!(bi, na)                     # set new ni

      # downstream change
      setdownstreamλμ!(λp, μp, i1, Ξ, idf)

      return t0, llr, ddλr, ddμr, ssλr, ssμr
    end
  end

  return t0, NaN, NaN, NaN, NaN, NaN
end




"""
    wfix_m(ξi ::cTfbd,
           ei ::Float64,
           λfs::Vector{Float64},
           μfs::Vector{Float64},
           eds::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for middle `i` branches with downstream **tips**.
"""
function wfix_m(ξi ::cTfbd,
                ei ::Float64,
                λfs::Vector{Float64},
                μfs::Vector{Float64},
                eds::Float64,
                acr::Float64)

  # select best from proposal
  sp, wt, λp, μp, pp = 0.0, 0, NaN, NaN, -Inf
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    p   = exp(- eds * (exp(λfi) + exp(μfi)))
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      μp  = μfi
      wt  = i
    end
  end

  # extract current λs and μs at time `t` and estimate ratio
  empty!(λfs)
  empty!(μfs)
  λc, μc = _λμat!(ξi, ei, λfs, μfs, 0.0, NaN, NaN)

  sc, pc = 0.0, NaN
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    p   = exp(- eds * (exp(λfi) + exp(μfs[i])))
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp/sc)

  return wt, λp, μp, pp, λc, μc, pc, acr
end




"""
    wfix_m(ξi ::cTfbd,
           ei ::Float64,
           eμ ::Bool,
           λfs::Vector{Float64},
           μfs::Vector{Float64},
           eds::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for middle `i` branches with downstream **tips**.
"""
function wfix_m(ξi ::cTfbd,
                ei ::Float64,
                eμ ::Bool,
                λfs::Vector{Float64},
                μfs::Vector{Float64},
                eds::Float64,
                acr::Float64)

  # select best from proposal
  sp, wt, λp, μp, pp = 0.0, 0, NaN, NaN, -Inf
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    μi  = exp(μfi)
    p   = exp(- eds * (exp(λfi) + μi)) * μi
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      μp  = μfi
      wt  = i
    end
  end

  # extract current λs and μs at time `t` and estimate ratio
  empty!(λfs)
  empty!(μfs)
  λc, μc = _λμat!(ξi, ei, λfs, μfs, 0.0, NaN, NaN)

  sc, pc = 0.0, NaN
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μi  = exp(μfs[i])
    p   = exp(- eds * (exp(λfi) + μi)) * μi 
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp/sc)

  return wt, λp, μp, pp, λc, μc, pc, acr
end




"""
    wfix_m(ξi ::cTfbd,
           ei ::Float64,
           λfs::Vector{Float64},
           μfs::Vector{Float64},
           eds::Float64,
           λ1 ::Float64,
           λ2 ::Float64,
           μ1 ::Float64, 
           μ2 ::Float64,
           αλ ::Float64,
           αμ ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for middle `i` branches with downstream **cladogenetic** daughters.
"""
function wfix_m(ξi ::cTfbd,
                ei ::Float64,
                λfs::Vector{Float64},
                μfs::Vector{Float64},
                eds::Float64,
                λ1 ::Float64,
                λ2 ::Float64,
                μ1 ::Float64, 
                μ2 ::Float64,
                αλ ::Float64,
                αμ ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                acr::Float64)

  # select best from proposal
  sp, wt, λp, μp, pp = 0.0, 0, NaN, NaN, -Inf
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    p   = dnorm2(λ1, λ2, λfi + αλ, σλ) * dnorm2(μ1, μ2, μfi + αμ, σμ) *
          exp(- eds * (exp(λfi) + exp(μfi)))
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      μp  = μfi
      wt  = i
    end
  end

  # extract current λs and μs at time `t` and estimate ratio
  empty!(λfs)
  empty!(μfs)
  λc, μc = _λμat!(ξi, ei, λfs, μfs, 0.0, NaN, NaN)

  sc, pc = 0.0, NaN
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    p   = dnorm2(λ1, λ2, λfi + αλ, σλ) * dnorm2(μ1, μ2, μfi + αμ, σμ) *
          exp(- eds * (exp(λfi) + exp(μfi)))
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood and acceptance ratio
  acr += log(sp/sc) + λp - λc
  ddλr = 2.0*(λc - λp)
  ddμr = 2.0*(μc - μp)
  ssλr = 0.5*((λ1 - λp - αλ)^2 + (λ2 - λp - αλ)^2 - 
              (λ1 - λc - αλ)^2 - (λ2 - λc - αλ)^2)
  ssμr = 0.5*((μ1 - μp - αμ)^2 + (μ2 - μp - αμ)^2 - 
              (μ1 - μc - αμ)^2 - (μ2 - μc - αμ)^2)

  return wt, λp, μp, pp, λc, μc, pc, acr, ddλr, ddμr, ssλr, ssμr
end




"""
    fsbi_i(bi ::iBffs,
           ξi ::cTfbd,
           λa ::Float64,
           λ1 ::Float64,
           λ2 ::Float64,
           μa ::Float64,
           μ1 ::Float64,
           μ2 ::Float64,
           αλ ::Float64,
           αμ ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           ψ  ::Vector{Float64},
           ψts::Vector{Float64},
           ixi::Int64,
           ixf::Int64,
           λfs::Vector{Float64},
           μfs::Vector{Float64})

Forward simulation for internal branch `bi`
"""
function fsbi_i(bi ::iBffs,
                ξi ::cTfbd,
                λa ::Float64,
                λ1 ::Float64,
                λ2 ::Float64,
                μa ::Float64,
                μ1 ::Float64,
                μ2 ::Float64,
                αλ ::Float64,
                αμ ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                ψ  ::Vector{Float64},
                ψts::Vector{Float64},
                ixi::Int64,
                ixf::Int64,
                λfs::Vector{Float64},
                μfs::Vector{Float64})

  # if does **not** come from a cladogenetic event
  if isnan(λa)
    λi, μi = lλ(ξi), lμ(ξi)
  else
    λi = rnorm(λa + αλ, σλ)
    μi = rnorm(μa + αμ, σμ)
  end

  # forward simulation during branch length
  empty!(λfs)
  empty!(μfs)
  nep = lastindex(ψts) + 1
  t0, na, af, nn = _sim_cladsfbd_i(ti(bi), tf(bi), λi, μi, αλ, αμ, σλ, σμ, ψ,
      ψts, ixi, nep, 0, false, 1, 500, λfs, μfs)

  if na < 1 || af || nn > 499
    return t0, NaN, NaN, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr  = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # choose most likely lineage to fix
  wt, λp, μp, pp, λc, μc, pc, acr, ddλr, ddμr, ssλr, ssμr = 
    wfix_i(ξi, e(bi), λfs, μfs, λ1, λ2, μ1, μ2, αλ, αμ, σλ, σμ, acr)

  if lU < acr

    # fix the tip
    if wt <= div(na,2)
      fixtip1!(t0, wt, 0)
    else
      fixtip2!(t0, na - wt + 1, 0)
    end

    # simulated remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), αλ, αμ, σλ, σμ, ψ, ψts, 
          ixf, nep, acr, lU, iρi, na, nn)
    end

    if lU < acr
      na -= 1
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc) + λp - λc
      setni!(bi, na)                     # set new ni
      setλt!(bi, λp)                     # set new λt
      setμt!(bi, μp)                     # set new λt

      return t0, llr, ddλr, ddμr, ssλr, ssμr
    end
  end

  return t0, NaN, NaN, NaN, NaN, NaN
end




"""
    wfix_i(ξi ::cTfbd,
           ei ::Float64,
           λfs::Vector{Float64},
           μfs::Vector{Float64},
           λ1 ::Float64,
           λ2 ::Float64,
           μ1 ::Float64,
           μ2 ::Float64,
           αλ ::Float64,
           αμ ::Float64,
           σλ ::Float64,
           σμ ::Float64,
           acr::Float64)

Choose most likely simulated lineage to fix with respect to daughter
for bifurcating `i` branches.
"""
function wfix_i(ξi ::cTfbd,
                ei ::Float64,
                λfs::Vector{Float64},
                μfs::Vector{Float64},
                λ1 ::Float64,
                λ2 ::Float64,
                μ1 ::Float64,
                μ2 ::Float64,
                αλ ::Float64,
                αμ ::Float64,
                σλ ::Float64,
                σμ ::Float64,
                acr::Float64)

  # select best from proposal
  sp, wt, λp, μp, pp = 0.0, 0, NaN, NaN, -Inf
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    μfi = μfs[i]
    p   = dnorm2(λ1, λ2, λfi + αλ, σλ) * dnorm2(μ1, μ2, μfi + αμ, σμ)
    sp += p
    if p > pp
      pp  = p
      λp  = λfi
      μp  = μfi
      wt  = i
    end
  end

  # extract current λs and μs at time `t` and estimate ratio
  empty!(λfs)
  empty!(μfs)
  λc, μc = _λμat!(ξi, ei, λfs, μfs, 0.0, NaN, NaN)

  sc, pc = 0.0, NaN
  for i in Base.OneTo(lastindex(λfs))
    λfi = λfs[i]
    p   = dnorm2(λ1, λ2, λfi + αλ, σλ) * dnorm2(μ1, μ2, μfs[i] + αμ, σμ)
    sc += p
    if λc === λfi
      pc = p
    end
  end

  # likelihood and acceptance ratio
  acr += log(sp/sc) + λp - λc
  ddλr = 2.0*(λc - λp)
  ddμr = 2.0*(μc - μp)
  ssλr = 0.5*((λ1 - λp - αλ)^2 + (λ2 - λp - αλ)^2 - 
              (λ1 - λc - αλ)^2 - (λ2 - λc - αλ)^2)
  ssμr = 0.5*((μ1 - μp - αμ)^2 + (μ2 - μp - αμ)^2 - 
              (μ1 - μc - αμ)^2 - (μ2 - μc - αμ)^2)

  return wt, λp, μp, pp, λc, μc, pc, acr, ddλr, ddμr, ssλr, ssμr
end




"""
    tip_sims!(tree::cTfbd,
              t   ::Float64,
              αλ  ::Float64,
              αμ  ::Float64,
              σλ  ::Float64,
              σμ  ::Float64,
              ψ   ::Vector{Float64},
              ψts ::Vector{Float64},
              ix  ::Int64,
              nep ::Int64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::cTfbd,
                   t   ::Float64,
                   αλ  ::Float64,
                   αμ  ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   ψ   ::Vector{Float64},
                   ψts ::Vector{Float64},
                   ix  ::Int64,
                   nep ::Int64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

 if lU < lr && na < 500

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nn, lr =
          _sim_cladsfbd_it(t, lλ(tree), lμ(tree), αλ, αμ, σλ, σμ, ψ, ψts, 
            ix, nep, lr, lU, iρi, na-1, nn, 500)

        if isnan(lr) || nn > 499
          return tree, na, nn, NaN
        end

        setproperty!(tree, :iμ, isextinct(stree))
        sete!(tree, e(tree) + e(stree))

        if def1(stree)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr = 
        tip_sims!(tree.d1, t, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, 
          lr, lU, iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, 
          lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::cTfbd,
                   t   ::Float64,
                   αλ  ::Float64,
                   αμ  ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   ψ   ::Vector{Float64},
                   ψts ::Vector{Float64},
                   ix  ::Int64,
                   nep ::Int64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

Continue simulation until time `t` for the fixed tip in `tree`.
"""
function fossiltip_sim!(tree::cTfbd,
                        t   ::Float64,
                        αλ  ::Float64,
                        αμ  ::Float64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        ψ   ::Vector{Float64},
                        ψts ::Vector{Float64},
                        ix  ::Int64,
                        nep ::Int64,
                        lr  ::Float64,
                        lU  ::Float64,
                        iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64)

  if lU < lr && nn < 500
    if istip(tree)

      nep = lastindex(ψts) + 1
      stree, na, nn, lr =
        _sim_cladsfbd_it(t, lλ(tree), lμ(tree), αλ, αμ, σλ, σμ, ψ, ψts, 
          ix, nep, lr, lU, iρi, na-1, nn, 500)

      if !isfinite(lr) || nn > 499
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr =
        fossiltip_sim!(tree.d1, t, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep,
          lr, lU, iρi, na, nn)
    else
      tree.d2, na, nn, lr =
        fossiltip_sim!(tree.d2, t, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep,
          lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




