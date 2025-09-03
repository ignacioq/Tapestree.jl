#=

constant occurrence birth-death MCMC using forward simulation

Jérémy Andréoletti

v(^-^v)

Created 11 02 2022
=#




"""
    insane_cobd(tree    ::sTf_label,
                ωtimes  ::Vector{Float64};
                λ_prior ::NTuple{2,Float64}     = (1.5, 1.0),
                μ_prior ::NTuple{2,Float64}     = (1.5, 1.0),
                ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                ω_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                ψω_epoch::Vector{Float64}       = Float64[],
                f_epoch ::Vector{Int64}         = Int64[0],
                niter   ::Int64                 = 1_000,
                nthin   ::Int64                 = 10,
                nburn   ::Int64                 = 200,
                nflushθ ::Int64                 = Int64(ceil(niter/5_000)),
                nflushΞ ::Int64                 = Int64(ceil(niter/100)),
                ofile   ::String                = homedir(),
                ϵi      ::Float64               = 0.4,
                λi      ::Float64               = NaN,
                μi      ::Float64               = NaN,
                ψi      ::Float64               = NaN,
                ωi      ::Float64               = NaN,
                pupdp   ::NTuple{5,Float64}     = (0.01, 0.01, 0.01, 0.01, 0.1),
                survival::Bool                  = true,
                prints  ::Int64                 = 5,
                mxthf   ::Float64               = 0.1,
                tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant occurrence birth-death.
"""
function insane_cobd(tree    ::sTf_label,
                     ωtimes  ::Vector{Float64};
                     λ_prior ::NTuple{2,Float64}     = (1.5, 1.0),
                     μ_prior ::NTuple{2,Float64}     = (1.5, 1.0),
                     ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     ω_prior ::NTuple{2,Float64}     = (1.0, 0.2),
                     ψω_epoch::Vector{Float64}       = Float64[],
                     f_epoch ::Vector{Int64}         = Int64[0],
                     niter   ::Int64                 = 1_000,
                     nthin   ::Int64                 = 10,
                     nburn   ::Int64                 = 200,
                     nflushθ ::Int64                 = Int64(ceil(niter/5_000)),
                     nflushΞ ::Int64                 = Int64(ceil(niter/100)),
                     ofile   ::String                = homedir(),
                     ϵi      ::Float64               = 0.4,
                     λi      ::Float64               = NaN,
                     μi      ::Float64               = NaN,
                     ψi      ::Float64               = NaN,
                     ωi      ::Float64               = NaN,
                     pupdp   ::NTuple{5,Float64}     = (0.01, 0.01, 0.01, 0.01, 0.1),
                     survival::Bool                  = true,
                     prints  ::Int64                 = 5,
                     mxthf   ::Float64               = 0.1,
                     tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n   = ntips(tree)
  th  = treeheight(tree)
  LTT = ltt(tree)

  # only include epochs where the tree occurs
  tix = findfirst(x -> x < th, ψω_epoch)
  if !isnothing(tix)
    ψω_epoch = ψω_epoch[tix:end]
  end
  nep = lastindex(ψω_epoch) + 1

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

  # make fix tree directory
  idf = make_idf(tree, tρ, th * mxthf)

  # starting parameters
  if isnan(λi) || isnan(μi)
    # if only one tip
    if isone(n)
      λc = λ_prior[1]/λ_prior[2]
      μc = μ_prior[1]/μ_prior[2]
    else
      λc, μc = moments(Float64(n), th, ϵi)
    end
  else
    λc, μc = λi, μi
  end
  
  if isnan(ψi)
    nψ = nfossils(tree)
    # if no sampled fossil
    if iszero(nψ)
      ψc = ψ_prior[1]/ψ_prior[2]
    else
      ψc = Float64(nψ)/treelength(tree)
    end
  else
    ψc = ψi
  end
  
  if isnan(ωi)
    nω = lastindex(ωtimes)
    # if no fossil occurrences
    if iszero(nω)
      # fω = fω_prior[1]/sum(fω_prior)  # expected fraction of fossils occurrences (i.e. without morphological characters)
      # ωc = ψc/(1/fω-1)
      ωc = ω_prior[1]/ω_prior[2]
    else
      ωc = Float64(nω)/treelength(tree)
    end
  else
    ωc = ωi
  end

  # make ψ and ω vectors for each epoch
  ψc = fill(ψc, nep)
  ωc = fill(ωc, nep)

  # make ω vector for each occurrence
  # sort!(ωtimes, rev=true)
  # ωvec = [ω0[end-something(findlast(ωtime .>= ψω_epochs), 0)] for ωtime in ωtimes]

  # count the number of occurrences in each epoch
  sort!(ωtimes,   rev=true)
  sort!(ψω_epoch, rev=true)

  if !isempty(ψω_epoch)
    nω = zeros(Int, nep)
    ep = 1

    for t in ωtimes
      while ep < nep && t < ψω_epoch[ep]
        ep += 1
      end
      nω[ep] += 1
    end
  else
    nω = [nω]
  end

  # condition on first speciation event
  rmλ = iszero(e(tree)) && !isfossil(tree) ? 1.0 : 0.0

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
  mc = m_surv_cbd(th, λc, μc, 5_000, surv)

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, sTfbd)

  # make epoch start vectors and indices for each `ξ`
  eixi = Int64[]
  eixf = Int64[]
  bst  = Float64[]
  for bi in idf
    tib = ti(bi)
    ei  = findfirst(x -> x < tib, ψω_epoch)
    ei  = isnothing(ei) ? nep : ei
    ef  = findfirst(x -> x < tf(bi), ψω_epoch)
    ef  = isnothing(ef) ? nep : ef
    push!(bst, tib)
    push!(eixi, ei)
    push!(eixf, ef)
  end

  # parameter updates (1: λ, 2: μ, 3: ψ, 4: ω, 5: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running constant occurrence birth-death"

  # adaptive phase
  llc, prc, λc, μc, ψc, ωc, mc, ns, L, LTT =
     mcmc_burn_cobd(Ξ, idf, ωtimes, LTT, λ_prior, μ_prior, ψ_prior, ω_prior, ψω_epoch,
        f_epoch, nburn, λc, μc, ψc, ωc, mc, nω, th, rmλ, surv, bst, eixi, eixf, pup, prints)

  # mcmc
  r, treev =
    mcmc_cobd(Ξ, idf, ωtimes, LTT, llc, prc, λc, μc, ψc, ωc, mc, ns, nω, L, 
      λ_prior, μ_prior, ψ_prior, ω_prior, ψω_epoch, f_epoch, th, rmλ, surv, bst, 
      eixi, eixf, pup, niter, nthin, nflushθ, nflushΞ, ofile, prints)

  return r, treev
end




"""
    mcmc_burn_cobd(Ξ       ::Vector{sTfbd},
                   idf     ::Array{iBffs,1},
                   ωtimes  ::Vector{Float64},
                   LTT     ::Ltt,
                   λ_prior ::NTuple{2,Float64},
                   μ_prior ::NTuple{2,Float64},
                   ψ_prior ::NTuple{2,Float64},
                   ω_prior ::NTuple{2,Float64},
                   ψω_epoch::Vector{Float64},
                   f_epoch ::Vector{Int64},
                   nburn   ::Int64,
                   λc      ::Float64,
                   μc      ::Float64,
                   ψc      ::Vector{Float64},
                   ωc      ::Vector{Float64},
                   mc      ::Float64,
                   nω      ::Vector{Int64},
                   th      ::Float64,
                   rmλ     ::Float64,
                   surv    ::Int64,
                   bst     ::Vector{Float64},
                   eixi    ::Vector{Int64},
                   eixf    ::Vector{Int64},
                   pup     ::Array{Int64,1},
                   prints  ::Int64)

Adaptive MCMC phase for da chain for constant occurrence birth-death using
forward simulation.
"""
function mcmc_burn_cobd(Ξ       ::Vector{sTfbd},
                        idf     ::Array{iBffs,1},
                        ωtimes  ::Vector{Float64},
                        LTT     ::Ltt,
                        λ_prior ::NTuple{2,Float64},
                        μ_prior ::NTuple{2,Float64},
                        ψ_prior ::NTuple{2,Float64},
                        ω_prior ::NTuple{2,Float64},
                        ψω_epoch::Vector{Float64},
                        f_epoch ::Vector{Int64},
                        nburn   ::Int64,
                        λc      ::Float64,
                        μc      ::Float64,
                        ψc      ::Vector{Float64},
                        ωc      ::Vector{Float64},
                        mc      ::Float64,
                        nω      ::Vector{Int64},
                        th      ::Float64,
                        rmλ     ::Float64,
                        surv    ::Int64,
                        bst     ::Vector{Float64},
                        eixi    ::Vector{Int64},
                        eixf    ::Vector{Int64},
                        pup     ::Array{Int64,1},
                        prints  ::Int64)

  el  = lastindex(idf)                     # number of branches
  L   = treelength(Ξ, ψω_epoch, bst, eixi) # tree length
  nψ  = nfossils(idf, ψω_epoch, f_epoch)   # number of fossilization events (in the tree) per epoch
  ns  = nnodesbifurcation(idf)             # number of speciation events
  ne  = Float64(ntipsextinct(Ξ))           # number of extinction events

  # likelihood
  llc = llik_cobd(Ξ, ωtimes, LTT, λc, μc, ψc, ωc, ns, ψω_epoch, bst, eixi) - 
        rmλ * log(λc) + log(mc) + prob_ρ(idf)
  prc = logdgamma(λc,      λ_prior[1], λ_prior[2])  +
        logdgamma(μc,      μ_prior[1], μ_prior[2])  +
        sum(logdgamma.(ψc, ψ_prior[1], ψ_prior[2]))  +
        sum(logdgamma.(ωc, ω_prior[1], ω_prior[2]))
        # sum(logdbeta.(ωc./(ψc.+ωc), fω_prior[1], fω_prior[2]))

  pbar = Progress(nburn, dt = prints, desc = "burn-in mcmc...", barlen = 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, sum(L), μc, mc, th, rmλ, surv, λ_prior)

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, sum(L), λc, mc, th, surv, μ_prior)

      # ψ proposal
      elseif p === 3

        llc, prc = update_ψ!(llc, prc, ψc, nψ, L, ψ_prior)

      # ω proposal
      elseif p === 4

        llc, prc = update_ω!(llc, prc, ωc, nω, L, ω_prior)

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L, LTT =
          update_fs!(bix, Ξ, idf, ωtimes, LTT, llc, λc, μc, ψc, ωc, 
            ψω_epoch, ns, ne, L, eixi, eixf)

      end

    end

    next!(pbar)
  end

  return llc, prc, λc, μc, ψc, ωc, mc, ns, L, LTT
end




"""
    mcmc_cobd(Ξ       ::Vector{sTfbd},
              idf     ::Array{iBffs,1},
              ωtimes  ::Vector{Float64},
              LTT     ::Ltt,
              llc     ::Float64,
              prc     ::Float64,
              λc      ::Float64,
              μc      ::Float64,
              ψc      ::Vector{Float64},
              ωc      ::Vector{Float64},
              mc      ::Float64,
              ns      ::Float64,
              nω      ::Vector{Int64},
              L       ::Vector{Float64},
              λ_prior ::NTuple{2,Float64},
              μ_prior ::NTuple{2,Float64},
              ψ_prior ::NTuple{2,Float64},
              ω_prior ::NTuple{2,Float64},
              ψω_epoch::Vector{Float64},
              f_epoch ::Vector{Int64},
              th      ::Float64,
              rmλ     ::Float64,
              surv    ::Int64,
              bst     ::Vector{Float64},
              eixi    ::Vector{Int64},
              eixf    ::Vector{Int64},
              pup     ::Array{Int64,1},
              niter   ::Int64,
              nthin   ::Int64,
              nflushθ ::Int64,
              nflushΞ ::Int64,
              ofile   ::String,
              prints  ::Int64)

MCMC da chain for constant occurrence birth-death using forward simulation.
"""
function mcmc_cobd(Ξ       ::Vector{sTfbd},
                   idf     ::Array{iBffs,1},
                   ωtimes  ::Vector{Float64},
                   LTT     ::Ltt,
                   llc     ::Float64,
                   prc     ::Float64,
                   λc      ::Float64,
                   μc      ::Float64,
                   ψc      ::Vector{Float64},
                   ωc      ::Vector{Float64},
                   mc      ::Float64,
                   ns      ::Float64,
                   nω      ::Vector{Int64},
                   L       ::Vector{Float64},
                   λ_prior ::NTuple{2,Float64},
                   μ_prior ::NTuple{2,Float64},
                   ψ_prior ::NTuple{2,Float64},
                   ω_prior ::NTuple{2,Float64},
                   ψω_epoch::Vector{Float64},
                   f_epoch ::Vector{Int64},
                   th      ::Float64,
                   rmλ     ::Float64,
                   surv    ::Int64,
                   bst     ::Vector{Float64},
                   eixi    ::Vector{Int64},
                   eixf    ::Vector{Int64},
                   pup     ::Array{Int64,1},
                   niter   ::Int64,
                   nthin   ::Int64,
                   nflushθ ::Int64,
                   nflushΞ ::Int64,
                   ofile   ::String,
                   prints  ::Int64)

  el  = lastindex(idf)                     # number of branches
  nψ  = nfossils(idf, ψω_epoch, f_epoch)   # number of fossilization events (in the tree) per epoch
  ne  = Float64(ntipsextinct(Ξ))           # number of extinction events
  nep = lastindex(ψc)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  r = Array{Float64,2}(undef, nlogs, 6 + 2*nep)

  treev  = sTfbd[]    # make tree vector
  sthinθ = sthinΞ = 0 # flush to file
  io     = IOBuffer() # buffer 

  function check_pr(pupi::Int64, i::Int64)
    pr0 = logdgamma(λc,      λ_prior[1], λ_prior[2])  +
          logdgamma(μc,      μ_prior[1], μ_prior[2])  +
          sum(logdgamma.(ψc, ψ_prior[1], ψ_prior[2])) +
          sum(logdgamma.(ωc, ω_prior[1], ω_prior[2]))
          # sum(logdbeta.(ωc./(ψc.+ωc), fω_prior[1], fω_prior[2]))
    if !isapprox(pr0, prc, atol = 1e-5)
       error(string("Wrong prior computation during the ", ["λ","μ","ψ","ω","forward simulation"][pupi], 
                    " update, at iteration ", i, ": pr0=", pr0, " and prc-pr0=", prc-pr0))
    end
  end

  function check_ll(pupi::Int64, i::Int64)
    ll0 = llik_cobd(Ξ, ωtimes, LTT, λc, μc, ψc, ωc, nnodesbifurcation(idf), ψω_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
    # @show llik_cfbd(Ξ, λc, μc, ψc, nnodesbifurcation(idf), ψω_epoch, bst, eixi)
    # @show ω_llik(ωtimes, ωc, ψω_epoch, LTT)
    # @show - rmλ, log(λc), log(mc)
    # @show prob_ρ(idf)
    if !isapprox(ll0, llc, atol = 1e-5)
       error(string("Wrong likelihood computation during the ", ["λ","μ","ψ","ω","forward simulation"][pupi], 
                    " update, at iteration ", i, ": ll0=", ll0, " and llc-ll0=", llc-ll0))
    end
  end

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda\tmu\tdiv\tns\tne\t"*join(["psi"*(isone(nep) ? "" : string("_",i)) for i in 1:nep], "\t")*"\t"*join(["omega"*(isone(nep) ? "" : string("_",i)) for i in 1:nep], "\t")*"\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

      for it in Base.OneTo(niter)

        # @show it
        shuffle!(pup)

        for p in pup
          # @show ["λ","μ","ψ","ω","forward simulation"][p]

          # λ proposal
          if p === 1

            llc, prc, λc, mc =
              update_λ!(llc, prc, λc, ns, sum(L), μc, mc, th, rmλ, surv, λ_prior)

          # μ proposal
          elseif p === 2

            llc, prc, μc, mc =
              update_μ!(llc, prc, μc, ne, sum(L), λc, mc, th, surv, μ_prior)

          # ψ proposal
          elseif p === 3

            llc, prc = update_ψ!(llc, prc, ψc, nψ, L, ψ_prior)

          # ω proposal
          elseif p === 4

            llc, prc = update_ω!(llc, prc, ωc, nω, L, ω_prior)

          # forward simulation proposal proposal
          else

            bix = ceil(Int64,rand()*el)

            llc, ns, ne, L, LTT =
              update_fs!(bix, Ξ, idf, ωtimes, LTT, llc, λc, μc, ψc, ωc, 
                ψω_epoch, ns, ne, L, eixi, eixf)

          end

          # check_pr(p, it)
          # check_ll(p, it)
        end

        # log parameters
        lthin += 1
        if lthin == nthin

          lit += 1
          @inbounds begin
            r[lit,1] = Float64(it)
            r[lit,2] = llc
            r[lit,3] = prc
            r[lit,4] = λc
            r[lit,5] = μc
            @turbo for i in Base.OneTo(nep)
              r[lit,6 + i] = ψc[i]
            end
            @turbo for i in Base.OneTo(nep)
              r[lit,6 + nep + i] = ωc[i]
            end
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthinθ += 1
        if sthinθ === nflushθ
          print(of, Float64(it), "\t", llc, "\t", prc, "\t", λc,"\t", μc, "\t", λc-μc, 
                "\t", ns, "\t", ne, "\t", join(ψc, "\t"), "\t", join(ωc, "\t"), "\n")
          flush(of)
          sthinθ = 0
        end
        sthinΞ += 1
        if sthinΞ === nflushΞ
          ibuffer(io, couple(Ξ, idf, 1))
          write(io, '\n')
          write(tf, take!(io))
          flush(tf)
          sthinΞ = 0
        end

        next!(pbar)
      end
    end
  end

  return r, treev
end




"""
    update_fs!(bix   ::Int64,
               Ξ     ::Vector{sTfbd},
               idf   ::Vector{iBffs},
               ωtimes::Vector{Float64},
               LTT   ::Ltt,
               llc   ::Float64,
               λ     ::Float64,
               μ     ::Float64,
               ψ     ::Vector{Float64},
               ω     ::Vector{Float64},
               ψωts  ::Vector{Float64},
               ns    ::Float64,
               ne    ::Float64,
               L     ::Vector{Float64},
               eixi  ::Vector{Int64},
               eixf  ::Vector{Int64})

Forward simulation proposal function for constant occurrence birth-death.
"""
function update_fs!(bix   ::Int64,
                    Ξ     ::Vector{sTfbd},
                    idf   ::Vector{iBffs},
                    ωtimes::Vector{Float64},
                    LTT   ::Ltt,
                    llc   ::Float64,
                    λ     ::Float64,
                    μ     ::Float64,
                    ψ     ::Vector{Float64},
                    ω     ::Vector{Float64},
                    ψωts  ::Vector{Float64},
                    ns    ::Float64,
                    ne    ::Float64,
                    L     ::Vector{Float64},
                    eixi  ::Vector{Int64},
                    eixf  ::Vector{Int64})

  bi  = idf[bix]
  ξc  = Ξ[bix]
  ixi = eixi[bix]

  if isfossil(bi)
    ixf = eixf[bix]
    # @show "fsbi_f"
    ξp, LTTp, llr = fsbi_f(ξc, bi, λ, μ, ψ, ω, ψωts, ωtimes, LTT, ixi, ixf)

    # if terminal but not successful proposal, update extinct
    if iszero(d1(bi)) && !isfinite(llr)
      # @show "fsbi_et"
      ξp, LTTp, llr = fsbi_et(sTfbd_wofe(ξc), bi, λ, μ, ψ, ω, ψωts, ωtimes, LTT, ixi, ixf)
    end
  else
    if iszero(d1(bi))
      # @show "fsbi_t"
      ξp, LTTp, llr = fsbi_t(ξc, bi, λ, μ, ψ, ω, ψωts, ωtimes, LTT, ixi)
    else
      # @show "fsbi_i"
      ξp, LTTp, llr = fsbi_i(ξc, bi, λ, μ, ψ, ω, ψωts, ωtimes, LTT, ixi, eixf[bix])
    end
  end

  if isfinite(llr)
    # @show bix
    # @show llr
    # llik_cfbdc = llik_cfbd(Ξ, λ, μ, ψ, nnodesbifurcation(idf), ψωts, bst, eixi)
    # ω_llikc = ω_llik(ωtimes, ω, ψωts, LTT)

    tii = ti(bi)

    nep = lastindex(ψωts) + 1
    # update llc, ns, ne & L
    # @show llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT)
    # @show ω_llik(ωtimes, ω, ψωts, LTTp) - ω_llik(ωtimes, ω, ψωts, LTT)
    # @show ξc, ξp
    # @show e(ξc)
    # @show e(ξp)
    # if def1(ξp) && def2(ξp)
    #   @show e(ξp.d1), e(ξp.d2)
    #   if def1(ξp.d1) && def2(ξp.d1)
    #     @show e(ξp.d1.d1), e(ξp.d1.d2)
    #   end
    #   if def1(ξp.d2) && def2(ξp.d2)
    #     @show e(ξp.d2.d1), e(ξp.d2.d2)
    #   end
    # end

    llr_cfbd = llik_cfbd(ξp, λ, μ, ψ, tii, ψωts, ixi, nep) - 
               llik_cfbd(ξc, λ, μ, ψ, tii, ψωts, ixi, nep)
    llc += llr_cfbd + llr

    ns  += Float64(nnodesbifurcation(ξp) - nnodesbifurcation(ξc))
    ne  += Float64(ntipsextinct(ξp)      - ntipsextinct(ξc))

    # update tree lengths
    Lc = zeros(nep)
    _treelength!(ξc, tii, Lc, ψωts, ixi, nep)
    _treelength!(ξp, tii, L,  ψωts, ixi, nep)
    @turbo for i in Base.OneTo(nep)
      L[i] -= Lc[i]
    end

    # set new decoupled tree
    Ξ[bix] = ξp

    # set new LTT
    # @assert ltt(couple(Ξ, idf, 1)).n == LTTp.n "[[ltt(couple(Ξ, idf, 1)).n, ltt(couple(Ξ, idf, 1)).t], [LTTp.n, LTTp.t]] = ", [[ltt(couple(Ξ, idf, 1)).n, ltt(couple(Ξ, idf, 1)).t], [LTTp.n, LTTp.t]]
    LTT = LTTp
 
    # llik_cfbdp = llik_cfbd(Ξ, λ, μ, ψ, nnodesbifurcation(idf), ψωts, bst, eixi)
    # ω_llikp = ω_llik(ωtimes, ω, ψωts, LTT)
    # @show llik_cfbdp
    # @show ω_llikp
    # @show llr_cfbd, llik_cfbdp - llik_cfbdc
    # @show llr, ω_llikp - ω_llikc
 end

  return llc, ns, ne, L, LTT
end




"""
    fsbi_t(ξc    ::sTfbd,
           bi    ::iBffs,
           λ     ::Float64,
           μ     ::Float64,
           ψ     ::Vector{Float64},
           ω     ::Vector{Float64},
           ψωts  ::Vector{Float64},
           ωtimes::Vector{Float64},
           LTT   ::Ltt,
           ix    ::Int64)

Forward simulation for terminal branch.
"""
function fsbi_t(ξc    ::sTfbd,
                bi    ::iBffs,
                λ     ::Float64,
                μ     ::Float64,
                ψ     ::Vector{Float64},
                ω     ::Vector{Float64},
                ψωts  ::Vector{Float64},
                ωtimes::Vector{Float64},
                LTT   ::Ltt,
                ix    ::Int64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  nep = lastindex(ψωts) + 1
  ξp, na, nn, llr =
    _sim_cfbd_t(e(bi), λ, μ, ψ, ψωts, ix, nep, lc, lU, iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr) && (treelength(ξc)!=treelength(ξp))
        
    llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ix)
    llr += llrLTTp

    if lU < llr != 0.0

      _fixrtip!(ξp, na) # fix random tip
      setni!(bi, na)    # set new ni

      return ξp, LTTp, llr
    end    
  end

  return ξp, LTT, -Inf
end




"""
    fsbi_f(ξc    ::sTfbd,
           bi    ::iBffs,
           λ     ::Float64,
           μ     ::Float64,
           ψ     ::Vector{Float64},
           ω     ::Vector{Float64},
           ψωts  ::Vector{Float64},
           ωtimes::Vector{Float64},
           LTT   ::Ltt,
           ixi   ::Int64,
           ixf   ::Int64)

Forward simulation for fossil branch `bi`.
"""
function fsbi_f(ξc    ::sTfbd,
                bi    ::iBffs,
                λ     ::Float64,
                μ     ::Float64,
                ψ     ::Vector{Float64},
                ω     ::Vector{Float64},
                ψωts  ::Vector{Float64},
                ωtimes::Vector{Float64},
                LTT   ::Ltt,
                ixi   ::Int64,
                ixf   ::Int64)

  # forward simulation during branch length
  nep = lastindex(ψωts) + 1
  ξp, na, nf, nn = 
    _sim_cfbd_i(ti(bi), tf(bi), λ, μ, ψ, ψωts, ixi, nep, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return ξp, LTT, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  if lU < acr

    _fixrtip!(ξp, na)

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(ξp, tf(bi), λ, μ, ψ, ψωts, ixf, acr, lU, iρi, na, nn)
    end

    if lU < acr
      # fossilize extant tip
      fossilizefixedtip!(ξp)

      if iszero(d1(bi))
        tx, na, nn, acr =
          fossiltip_sim!(ξp, tf(bi), λ, μ, ψ, ψωts, ixf, acr, lU, iρi, na, nn)
      else
        na -= 1
      end

      if lU < acr
        
        llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ixi)
        
        if lU < acr + llrLTTp != 0.0

          # @show llrLTTp
          llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + llrLTTp
          setnt!(bi, ntp)                # set new nt
          setni!(bi, na)                 # set new ni

          return ξp, LTTp, llr
        end
      end
    end
  end

  return ξp, LTT, NaN
end




"""
    fsbi_et(ξc    ::sTfbd,
            bi    ::iBffs,
            λ     ::Float64,
            μ     ::Float64,
            ψ     ::Vector{Float64},
            ω     ::Vector{Float64},
            ψωts  ::Vector{Float64},
            ωtimes::Vector{Float64},
            LTT   ::Ltt,
            ixi   ::Int64,
            ixf   ::Int64)

Forward simulation for extinct tip in terminal fossil branch.
"""
function fsbi_et(ξc    ::sTfbd,
                 bi    ::iBffs,
                 λ     ::Float64,
                 μ     ::Float64,
                 ψ     ::Vector{Float64},
                 ω     ::Vector{Float64},
                 ψωts  ::Vector{Float64},
                 ωtimes::Vector{Float64},
                 LTT   ::Ltt,
                 ixi   ::Int64,
                 ixf   ::Int64)

  lU  = -randexp()            # log-probability
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  ξp, na, nn, acr =
    fossiltip_sim!(ξc, tf(bi), λ, μ, ψ, ψωts, ixf, acr, lU, iρi, 1, 1)

  if lU < acr && (treelength(ξc)!=treelength(ξp))
    
    llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ixi)

    if lU < acr + llrLTTp != 0.0
      llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + llrLTTp
      setni!(bi, na)                 # set new ni

      return ξp, LTTp, llr
    end
  end

  return ξp, LTT, NaN
end




"""
    fsbi_i(ξc    ::sTfbd,
           bi    ::iBffs,
           λ     ::Float64,
           μ     ::Float64,
           ψ     ::Vector{Float64},
           ω     ::Vector{Float64},
           ψωts  ::Vector{Float64},
           ωtimes::Vector{Float64},
           LTT   ::Ltt,
           ixi   ::Int64,
           ixf   ::Int64)

Forward simulation for internal branch `bi`.
"""
function fsbi_i(ξc    ::sTfbd,
                bi    ::iBffs,
                λ     ::Float64,
                μ     ::Float64,
                ψ     ::Vector{Float64},
                ω     ::Vector{Float64},
                ψωts  ::Vector{Float64},
                ωtimes::Vector{Float64},
                LTT   ::Ltt,
                ixi   ::Int64,
                ixf   ::Int64)

  # forward simulation during branch length
  nep = lastindex(ψωts) + 1
  ξp, na, nf, nn = 
    _sim_cfbd_i(ti(bi), tf(bi), λ, μ, ψ, ψωts, ixi, nep, 0, 0, 1, 1_000)

  if na < 1 || nf > 0 || nn > 999
    return ξp, LTT, NaN
  end

  ntp = na

  lU = -randexp() # log-probability

  # acceptance probability
  acr  = log(Float64(ntp)/Float64(nt(bi)))
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  if lU < acr

    _fixrtip!(ξp, na)

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(ξp, tf(bi), λ, μ, ψ, ψωts, ixf, acr, lU, iρi, na, nn)
    end

    if lU < acr
      
      llrLTTp, LTTp = llrLTT(ξc, ξp, bi, ωtimes, ω, ψωts, LTT, ixi)

      if lU < acr + llrLTTp != 0.0

        na -= 1
        llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + llrLTTp
        setnt!(bi, ntp)                # set new nt
        setni!(bi, na)                 # set new ni

        return ξp, LTTp, llr
      end
    end
  end

  return ξp, LTT, NaN
end




"""
    update_ω!(llc    ::Float64,
              prc    ::Float64,
              ωc     ::Vector{Float64},
              LTT    ::Ltt,
              nω     ::Vector{Int64},
              L      ::Vector{Float64},
              ω_prior::NTuple{2,Float64})

Gibbs sampling of `ω` for constant occurrence birth-death.
"""
function update_ω!(llc    ::Float64,
                   prc    ::Float64,
                   ωc     ::Vector{Float64},
                   nω     ::Vector{Int64},
                   L      ::Vector{Float64},
                   ω_prior::NTuple{2,Float64})

  # MH steps for each epoch
  for i in Base.OneTo(lastindex(ωc))

    ωp  = rand(Gamma(ω_prior[1]+nω[i], ω_prior[2]+L[i]))
    ωci = ωc[i]
    
    prc += llrdgamma(ωp, ωci, ω_prior[1], ω_prior[2])
    llc += nω[i] * log(ωp/ωci) + L[i] * (ωci - ωp)
    ωc[i] = ωp
    
  end

  return llc, prc
end


