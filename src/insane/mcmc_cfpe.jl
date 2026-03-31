#=

constant fossil punkeek MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2025
=#



"""
    insane_cfpe(tree    ::sTf_label,
                xa      ::Dict{String, Float64};
                xs      ::Dict{String, Float64} = Dict{String,Float64}(),
                λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                σa_prior::NTuple{2,Float64}     = (0.05, 0.05),
                σk_prior::NTuple{2,Float64}     = (0.05, 0.05),
                ψ_epoch ::Vector{Float64}       = Float64[],
                f_epoch ::Vector{Int64}         = Int64[0],
                niter   ::Int64                 = 1_000,
                nthin   ::Int64                 = 10,
                nburn   ::Int64                 = 200,
                nflush  ::Int64                 = nthin,
                ofile   ::String                = string(homedir(), "/cfpe"),
                ϵi      ::Float64               = 0.4,
                λi      ::Float64               = NaN,
                μi      ::Float64               = NaN,
                ψi      ::Float64               = NaN,
                αai     ::Float64               = 0.0,
                αki     ::Float64               = 0.0,
                pupdp   ::NTuple{7,Float64}     = (1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-1, 0.2),
                survival::Bool                  = true,
                mxthf   ::Float64               = 0.1,
                tρ      ::Dict{String, Float64} = Dict("" => 1.0),
                prints  ::Int64                 = 5)

Run insane for constant fossilised birth-death punctuated equilibrium.
"""
function insane_cfpe(tree    ::sTf_label,
                     xa      ::Dict{String, Float64};
                     xs      ::Dict{String, Float64} = Dict{String,Float64}(),
                     λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                     σa_prior::NTuple{2,Float64}     = (0.05, 0.05),
                     σk_prior::NTuple{2,Float64}     = (0.05, 0.05),
                     ψ_epoch ::Vector{Float64}       = Float64[],
                     f_epoch ::Vector{Int64}         = Int64[0],
                     niter   ::Int64                 = 1_000,
                     nthin   ::Int64                 = 10,
                     nburn   ::Int64                 = 200,
                     nflush  ::Int64                 = nthin,
                     ofile   ::String                = string(homedir(), "/cfpe"),
                     ϵi      ::Float64               = 0.4,
                     λi      ::Float64               = NaN,
                     μi      ::Float64               = NaN,
                     ψi      ::Float64               = NaN,
                     αi      ::Float64               = 0.0,
                     pupdp   ::NTuple{8,Float64}     = (1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-1, 0.2),
                     survival::Bool                  = true,
                     mxthf   ::Float64               = 0.1,
                     tρ      ::Dict{String, Float64} = Dict("" => 1.0),
                     prints  ::Int64                 = 5)

  n  = ntips(tree)
  th = treeheight(tree)

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
        for i in Base.OneTo(nep - lep)
          pushfirst!(f_epoch, 0)
        end
      end
    else
      f_epoch = fill(0, nep)
    end
  end

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxi = make_idf(tree, tρ, xa, xs, th * mxthf)

  # starting parameters
  λc, μc, ψc = λi, μi, ψi
  if any(isnan, (λi, μi, ψi))
    # if only one tip
    if isone(n)
      λc = prod(λ_prior)
      μc = prod(μ_prior)
    else
      λc, μc = moments(Float64(n), th, ϵi)
    end
    # if no sampled fossil
    nf = nfossils(tree) + sum(f_epoch)
    if iszero(nf)
      ψc = prod(ψ_prior)
    else
      ψc = Float64(nf)/treelength(tree)
    end
  end

  # make ψ vector
  ψc = fill(ψc, nep)

  σac = σkc = σxi
  αc  = αi

  # if condition on first speciation event
  rmλ = Float64(iszero(e(tree)) && !isfossil(tree))

  # condition on survival of 0, 1, or 2 starting lineages
  surv = 0
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

  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 1_000, surv)

  # make epoch start vectors and indices for each `ξ`
  eixi = Int64[]
  eixf = Int64[]
  bst  = Float64[]
  for bi in idf
    tib = ti(bi)
    ei  = findfirst(x -> x < tib, ψ_epoch)
    ei  = isnothing(ei) ? nep : ei
    ef  = findfirst(x -> x < tf(bi), ψ_epoch)
    ef  = isnothing(ef) ? nep : ef
    push!(bst, tib)
    push!(eixi, ei)
    push!(eixf, ef)
  end

  # get vector of internal edges
  inodes = [i for i in Base.OneTo(lastindex(idf)) if d1(idf[i]) > 0]

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(8)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, xr, σkc, sTfpe)

  @info "running constant fossilised punctuated equilibrium"

  # adaptive phase
  llc, prc, λc, μc, ψc, αc, σac, σkc, mc, ns, ne, nf, L, dα, sσa, sσk, nσs =
    mcmc_burn_cfpe(Ξ, idf, 
        λ_prior, μ_prior, ψ_prior, α_prior, σa_prior, σk_prior, 
        ψ_epoch, f_epoch, nburn, λc, μc, ψc, αc, σac, σkc, mc, th, rmλ, 
        inodes, surv, bst, eixi, eixf, pup, prints)

  # mcmc
  r, treev = 
    mcmc_cfpe(Ξ, idf, llc, prc, λc, μc, ψc, αc, σac, σkc, mc, ns, ne, nf, L, 
      dα, sσa, sσk, nσs, th, rmλ, inodes, surv, bst, eixi, eixf, λ_prior, 
      μ_prior, ψ_prior, α_prior, σa_prior, σk_prior, ψ_epoch, f_epoch, 
      pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end



"""
    mcmc_burn_cfpe(Ξ       ::Vector{sTfpe},
                   idf     ::Array{iBffs,1},
                   λ_prior ::NTuple{2,Float64},
                   μ_prior ::NTuple{2,Float64},
                   ψ_prior ::NTuple{2,Float64},
                   α_prior ::NTuple{2,Float64},
                   σa_prior::NTuple{2,Float64},
                   σk_prior::NTuple{2,Float64},
                   ψ_epoch ::Vector{Float64},
                   f_epoch ::Vector{Int64},
                   nburn   ::Int64,
                   λc      ::Float64,
                   μc      ::Float64,
                   ψc      ::Vector{Float64},
                   αc      ::Float64,
                   σac     ::Float64,
                   σkc     ::Float64,
                   mc      ::Float64,
                   th      ::Float64,
                   rmλ     ::Float64,
                   inodes  ::Vector{Int64},
                   surv    ::Int64,
                   bst     ::Vector{Float64},
                   eixi    ::Vector{Int64},
                   eixf    ::Vector{Int64},
                   pup     ::Array{Int64,1},
                   prints  ::Int64)

Burn-in for constant birth-death punctuated equilibrium.
"""
function mcmc_burn_cfpe(Ξ       ::Vector{sTfpe},
                        idf     ::Array{iBffs,1},
                        λ_prior ::NTuple{2,Float64},
                        μ_prior ::NTuple{2,Float64},
                        ψ_prior ::NTuple{2,Float64},
                        α_prior ::NTuple{2,Float64},
                        σa_prior::NTuple{2,Float64},
                        σk_prior::NTuple{2,Float64},
                        ψ_epoch ::Vector{Float64},
                        f_epoch ::Vector{Int64},
                        nburn   ::Int64,
                        λc      ::Float64,
                        μc      ::Float64,
                        ψc      ::Vector{Float64},
                        αc      ::Float64,
                        σac     ::Float64,
                        σkc     ::Float64,
                        mc      ::Float64,
                        th      ::Float64,
                        rmλ     ::Float64,
                        inodes  ::Vector{Int64},
                        surv    ::Int64,
                        bst     ::Vector{Float64},
                        eixi    ::Vector{Int64},
                        eixf    ::Vector{Int64},
                        pup     ::Array{Int64,1},
                        prints  ::Int64)

  el  = lastindex(idf)
  L   = treelength(Ξ, ψ_epoch, bst, eixi) # tree length
  nf  = nfossils(idf, ψ_epoch, f_epoch)   # number of fossilization events per epoch
  ns  = nnodesbifurcation(idf)            # number of speciation events
  nin = lastindex(inodes)                 # number of internal nodes
  ne  = Float64(ntipsextinct(Ξ))          # number of extinction events

  # likelihood
  llc = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, 
          ns, ψ_epoch, f_epoch, bst, eixi) - 
        rmλ * log(λc) + log(mc) + prob_ρ(idf)

  # prior
  prc = logdgamma(λc, λ_prior[1], λ_prior[2])              +
        logdgamma(μc, μ_prior[1], μ_prior[2])              +
        logdnorm(αc,         α_prior[1],  α_prior[2]^2)    +
        logdinvgamma(σac^2, σa_prior[1], σa_prior[2])      +
        logdinvgamma(σkc^2, σk_prior[1], σk_prior[2])      +
        sum(x -> logdgamma(x, ψ_prior[1], ψ_prior[2]), ψc)

  # tracked quantities
  dα, sσa, sσk = gibbs_quanta(Ξ, idf, αc)

  # n number to sum to ns for σa updates
  nσs = nedgesF(Ξ) - 2.0*ns - Float64(iszero(e(idf[1])))

 # empty vectors
  xis = Float64[]
  xfs = Float64[]
  es  = Float64[]
  pv  = Float64[]

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

        llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # α proposal
      elseif p === 4

        llc, prc, αc, sσa = 
          update_α!(αc, σac, sum(L), dα, llc, prc, sσa, α_prior)

      # σa (anagenetic) proposal
      elseif p === 5

        llc, prc, σac = 
          update_σ!(σac, sσa, 2.0*ns + nσs, llc, prc, σa_prior)

      # σk (cladogenetic) proposal
      elseif p === 6

        llc, prc, σkc = update_σ!(σkc, sσk, ns, llc, prc, σk_prior)

      # update inner nodes traits
      elseif p === 7

        bix = inodes[fIrand(nin) + 1]
        llc, dα, sσa, sσk = 
          update_x!(bix, Ξ, idf, αc, σac, σkc, llc, dα, sσa, sσk)

      # forward simulation proposal proposal
      else

        bix = fIrand(el) + 1
        llc, ns, ne, dα, sσa, sσk = 
          update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ψ_epoch, αc, σac, σkc, 
            ns, ne, L, eixi, eixf, dα, sσa, sσk, xis, xfs, es, pv)

      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, ψc, αc, σac, σkc, 
           mc, ns, ne, nf, L, dα, sσa, sσk, nσs
end




"""
    mcmc_cfpe(Ξ       ::Vector{sTfpe},
              idf     ::Array{iBffs,1},
              llc     ::Float64,
              prc     ::Float64,
              λc      ::Float64,
              μc      ::Float64,
              ψc      ::Vector{Float64},
              αc      ::Float64,
              σac     ::Float64,
              σkc     ::Float64,
              mc      ::Float64,
              ns      ::Float64,
              ne      ::Float64,
              nf      ::Vector{Float64},
              L       ::Vector{Float64},
              dα      ::Float64, 
              sσa     ::Float64, 
              sσk     ::Float64,
              nσs     ::Float64,
              th      ::Float64,
              rmλ     ::Float64,
              inodes  ::Vector{Int64},
              surv    ::Int64,
              bst     ::Vector{Float64},
              eixi    ::Vector{Int64},
              eixf    ::Vector{Int64},
              λ_prior ::NTuple{2,Float64},
              μ_prior ::NTuple{2,Float64},
              ψ_prior ::NTuple{2,Float64},
              α_prior ::NTuple{2,Float64},
              σa_prior::NTuple{2,Float64},
              σk_prior::NTuple{2,Float64},
              ψ_epoch ::Vector{Float64},
              f_epoch ::Vector{Int64},
              pup     ::Vector{Int64},
              niter   ::Int64,
              nthin   ::Int64,
              nflush  ::Int64,
              ofile   ::String,
              prints  ::Int64)

Sampling for constant birth-death punctuated equilibrium.
"""
function mcmc_cfpe(Ξ       ::Vector{sTfpe},
                   idf     ::Array{iBffs,1},
                   llc     ::Float64,
                   prc     ::Float64,
                   λc      ::Float64,
                   μc      ::Float64,
                   ψc      ::Vector{Float64},
                   αc      ::Float64,
                   σac     ::Float64,
                   σkc     ::Float64,
                   mc      ::Float64,
                   ns      ::Float64,
                   ne      ::Float64,
                   nf      ::Vector{Float64},
                   L       ::Vector{Float64},
                   dα      ::Float64, 
                   sσa     ::Float64, 
                   sσk     ::Float64,
                   nσs     ::Float64,
                   th      ::Float64,
                   rmλ     ::Float64,
                   inodes  ::Vector{Int64},
                   surv    ::Int64,
                   bst     ::Vector{Float64},
                   eixi    ::Vector{Int64},
                   eixf    ::Vector{Int64},
                   λ_prior ::NTuple{2,Float64},
                   μ_prior ::NTuple{2,Float64},
                   ψ_prior ::NTuple{2,Float64},
                   α_prior ::NTuple{2,Float64},
                   σa_prior::NTuple{2,Float64},
                   σk_prior::NTuple{2,Float64},
                   ψ_epoch ::Vector{Float64},
                   f_epoch ::Vector{Int64},
                   pup     ::Vector{Int64},
                   niter   ::Int64,
                   nthin   ::Int64,
                   nflush  ::Int64,
                   ofile   ::String,
                   prints  ::Int64)

  el  = lastindex(idf)
  nin = lastindex(inodes)
  nep = lastindex(ψc)

  # logging
  nlogs = fld(niter,nthin)
  lthin = lit = sthin = zero(Int64)

  # parameter results
  r = Array{Float64,2}(undef, nlogs, 9 + nep)

  # empty vector
  xis = Float64[]
  xfs = Float64[]
  es  = Float64[]
  pv  = Float64[]

  treev = sTfpe[]     # make tree vector
  io    = IOBuffer() # buffer 

  open(ofile*".log", "w") do of 
    write(of, "iteration\tlikelihood\tprior\tlambda\tmu\t"*join(["psi"*(isone(nep) ? "" : string("_",i)) for i in 1:nep], '\t')*"\tx0\talpha\tsigma_a\tsigma_k\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, λc = λc, μc = μc, αc = αc, σac = σac, σkc = σkc, mc = mc, ns = ns, ne = ne, dα = dα, sσa = sσa, sσk = sσk, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for p in pup

             # λ proposal
            if p === 1

              llc, prc, λc, mc =
                update_λ!(llc, prc, λc, ns, sum(L), μc, mc, th, rmλ, surv, λ_prior)

              # llci = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, nnodesbifurcation(idf), ψ_epoch, f_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # μ proposal
            elseif p === 2

              llc, prc, μc, mc =
                update_μ!(llc, prc, μc, ne, sum(L), λc, mc, th, surv, μ_prior)

              # llci = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, nnodesbifurcation(idf), ψ_epoch, f_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # ψ proposal
            elseif p === 3

              llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

              # llci = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, nnodesbifurcation(idf), ψ_epoch, f_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # α proposal
            elseif p === 4

              llc, prc, αc, sσa = 
                update_α!(αc, σac, sum(L), dα, llc, prc, sσa, α_prior)

              # llci = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, nnodesbifurcation(idf), ψ_epoch, f_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # σa (anagenetic) proposal
            elseif p === 5

              llc, prc, σac = 
                update_σ!(σac, sσa, 2.0*ns + nσs, llc, prc, σa_prior)

              # llci = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, nnodesbifurcation(idf), ψ_epoch, f_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # σk (cladogenetic) proposal
            elseif p === 6

              llc, prc, σkc = update_σ!(σkc, sσk, ns, llc, prc, σk_prior)

              # llci = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, nnodesbifurcation(idf), ψ_epoch, f_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # update inner nodes traits
            elseif p === 7

              bix = inodes[fIrand(nin) + 1]

              llc, dα, sσa, sσk = 
                update_x!(bix, Ξ, idf, αc, σac, σkc, llc, dα, sσa, sσk)

              # llci = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, nnodesbifurcation(idf), ψ_epoch, f_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            # forward simulation proposal proposal
            else

              bix = fIrand(el) + 1

              llc, ns, ne, dα, sσa, sσk = 
                update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ψ_epoch, αc, σac, σkc, 
                  ns, ne, L, eixi, eixf, dα, sσa, sσk, xis, xfs, es, pv)

              # llci = llik_cfpe(Ξ, idf, λc, μc, ψc, αc, σac, σkc, nnodesbifurcation(idf), ψ_epoch, f_epoch, bst, eixi) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
              # if !isapprox(llci, llc, atol = 1e-6)
              #   @show llci, llc, it, p
              #   return
              # end

            end
          end

          # log parameters
          lthin += one(Int64)
          if lthin === nthin

            lit += 1
            @inbounds begin
              r[lit,1] = Float64(it)
              r[lit,2] = llc
              r[lit,3] = prc
              r[lit,4] = λc
              r[lit,5] = μc
              @turbo for i in Base.OneTo(nep)
                r[lit,5 + i] = ψc[i]
              end
              r[lit, 6 + nep] = xi(Ξ[1])
              r[lit, 7 + nep] = αc
              r[lit, 8 + nep] = σac
              r[lit, 9 + nep] = σkc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += one(Int64)
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', λc,'\t', μc, '\t', join(ψc, '\t'), '\t', xi(Ξ[1]), '\t', αc, '\t', σac, '\t', σkc, '\n')
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
    update_x!(bix::Int64,
              Ξ  ::Vector{sTfpe},
              idf::Vector{iBffs},
              α  ::Float64,
              σa ::Float64,
              σk ::Float64,
              ll ::Float64,
              dα ::Float64,
              sσa::Float64,
              sσk::Float64)

Perform a punkeek trait update for an internal branch and its descendants.
"""
function update_x!(bix::Int64,
                   Ξ  ::Vector{sTfpe},
                   idf::Vector{iBffs},
                   α  ::Float64,
                   σa ::Float64,
                   σk ::Float64,
                   ll ::Float64,
                   dα ::Float64,
                   sσa::Float64,
                   sσk::Float64)

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1, i2  = d1(bi), d2(bi) 
  b1   = idf[i1]
  ξ1   = Ξ[i1]
  root = iszero(pa(bi))

  ## update parent
  # if mrca
  if root && iszero(e(ξi))
    #if stem fossil
    if isfossil(bi)
      if !ifx(bi)
        ll, dα, sσa = _fstem_update!(ξi, ξ1, α, σa, ll, dα, sσa)
      end
    else
      # if crown
      ξ2  = Ξ[i2]
      ll, dα, sσa, sσk = _crown_update!(ξi, ξ1, ξ2, α, σa, σk, ll, dα, sσa, sσk)
    end
  else
    # if stem
    if root
      ll, dα, sσa = _stem_update!(ξi, α, σa, ll, dα, sσa)
    end

    # updates within the parent branch
    ll, dα, sσa, sσk = 
      _update_node!(ξi, NaN, NaN, α, σa, σk, ll, dα, sσa, sσk, false)

    # get fixed tip
    lξi = fixtip(ξi)

    isd = iszero(i2)
    # if duo
    if isd
      xavi = xavg(bi)
      if isnan(xavi)
        ll, sσa = _update_duo!(lξi, ξ1, xavi, xstd(bi), α, σa, ll, sσa)
      end
    # if triad
    else
      ll, dα, sσa, sσk = 
        _update_quartet!(lξi, ξ1, Ξ[i2], α, σa, σk, ll, dα, sσa, sσk)
    end

    ## update daughters
    b1 = idf[i1]
    ll, dα, sσa, sσk = 
      _update_node!(ξ1, xavg(b1), xstd(b1), 
        α, σa, σk, ll, dα, sσa, sσk, iszero(d1(b1)))

    # if triad
    if !isd
      b2 = idf[i2]
      ll, dα, sσa, sσk = 
        _update_node!(Ξ[i2], xavg(b2), xstd(b2), 
          α, σa, σk, ll, dα, sσa, sσk, iszero(d1(b2)))
    end
  end

  return ll, dα, sσa, sσk
end




"""
    update_fs!(bix ::Int64,
               Ξ   ::Vector{sTfpe},
               idf ::Vector{iBffs},
               llc ::Float64,
               λ   ::Float64,
               μ   ::Float64,
               ψ   ::Vector{Float64},
               ψts ::Vector{Float64},
               α   ::Float64,
               σa  ::Float64,
               σk  ::Float64,
               ns  ::Float64,
               ne  ::Float64,
               L   ::Vector{Float64},
               eixi::Vector{Int64},
               eixf::Vector{Int64},
               dα  ::Float64, 
               sσa ::Float64, 
               sσk ::Float64,
               xis ::Vector{Float64},
               xfs ::Vector{Float64},
               es  ::Vector{Float64},
               pv  ::Vector{Float64})

Forward simulation proposal function for constant punkeek.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{sTfpe},
                    idf ::Vector{iBffs},
                    llc ::Float64,
                    λ   ::Float64,
                    μ   ::Float64,
                    ψ   ::Vector{Float64},
                    ψts ::Vector{Float64},
                    α   ::Float64,
                    σa  ::Float64,
                    σk  ::Float64,
                    ns  ::Float64,
                    ne  ::Float64,
                    L   ::Vector{Float64},
                    eixi::Vector{Int64},
                    eixf::Vector{Int64},
                    dα  ::Float64, 
                    sσa ::Float64, 
                    sσk ::Float64,
                    xis ::Vector{Float64},
                    xfs ::Vector{Float64},
                    es  ::Vector{Float64},
                    pv  ::Vector{Float64})

  bi  = idf[bix]
  ξc  = Ξ[bix]
  ixi = eixi[bix]
  ixf = eixf[bix]

  llr = NaN
  dαr = sσar = sσkr = 0.0

   # if non-bifurcating branch
  if iszero(d2(bi))

    xav = xsd = NaN
    if !isnothing(xavg(bi))
      xav, xsd = xavg(bi), xstd(bi)
    end

    # if terminal branch
    if iszero(d1(bi))

      # if fossil terminal branch
      if isfossil(bi)
        ξp, llr = 
          fsbi_f(bi, xav, xsd, ξc, λ, μ, ψ, α, σa, σk, ψts, 
            ixi, ixf, xis, xfs, es, pv)

        # if not successful proposal, update extinct daughter
        if !isfinite(llr)
          ξp, llr = fsbi_et(sTfpe_wofe(ξc), bi, λ, μ, ψ, α, σa, σk, ψts, ixf)
        end

      # if terminal non-fossil branch
      else
        ξp, llr = fsbi_t(bi, xav, xsd, ξc, λ, μ, ψ, α, σa, σk, ψts, ixi, 
                    xis, xfs, es, pv)
      end

    # if mid (fossil or not) branch
    else
      ξp, llr, dαr, sσar = 
        fsbi_m(bi, xav, xsd, ξc, Ξ[d1(bi)], λ, μ, ψ, α, σa, σk, ψts, ixi, ixf, 
          xis, xfs, es, pv)
    end

  # if bifurcating branch
  elseif e(bi) > 0.0
    ξp, llr, dαr, sσar, sσkr = 
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], λ, μ, ψ, α, σa, σk, ψts, 
        ixi, ixf, xfs, pv)
  end

  if isfinite(llr)

    nep = lastindex(ψts) + 1
    σa2, σk2 = σa^2, σk^2
    llc, ns, ne, dα, sσa, sσk = 
      llik_cfpe_track!(ξc, λ, μ, ψ, α, σa2, σk2, llc, ns, ne, L, dα, sσa, sσk, 
        ti(bi), ψts, ixi, nep, -)
    llc, ns, ne, dα, sσa, sσk = 
      llik_cfpe_track!(ξp, λ, μ, ψ, α, σa2, σk2, llc, ns, ne, L, dα, sσa, sσk, 
        ti(bi), ψts, ixi, nep, +)

    llc += llr
    dα  += dαr
    sσa += sσar
    sσk += sσkr

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, ns, ne, dα, sσa, sσk
end





"""
    fsbi_f(bi ::iBffs,
           xav::Float64,
           xst::Float64,
           ξi ::sTfpe,
           λ  ::Float64,
           μ  ::Float64,
           ψ  ::Vector{Float64},
           α  ::Float64,
           σa ::Float64,
           σk ::Float64,
           ψts::Vector{Float64},
           ixi::Int64,
           ixf::Int64,
           xis::Vector{Float64},
           xfs::Vector{Float64},
           es ::Vector{Float64})

Forward simulation for **fossil** terminal branch.
"""
function fsbi_f(bi ::iBffs,
                xav::Float64,
                xst::Float64,
                ξi ::sTfpe,
                λ  ::Float64,
                μ  ::Float64,
                ψ  ::Vector{Float64},
                α  ::Float64,
                σa ::Float64,
                σk ::Float64,
                ψts::Vector{Float64},
                ixi::Int64,
                ixf::Int64,
                xis::Vector{Float64},
                xfs::Vector{Float64},
                es ::Vector{Float64}, 
                pv ::Vector{Float64})

  # forward simulation during branch length
  empty!(xis)
  empty!(xfs)
  empty!(es)
  nep = lastindex(ψts) + 1

  t0, na, nf, nn = 
    _sim_cfpe_i(ti(bi), tf(bi), λ, μ, ψ, xi(ξi), α, σa, σk, ψts, ixi, nep, 
      0, 0, 1, 500, xis, xfs, es)

  if na < 1 || nf > 0 || nn > 499
    return t0, NaN
  end

  lU = -randexp() # log-probability

  # add sampling fraction
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # inverse branch sampling fraction
  acr = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # if fixed node
  if ifx(bi)

    # propose trait value (if no uncertainty, then xp = xav)
    xp = xav
    if xst > 0.0
      xp = rnorm(xav, xst)
    end
    wt, acr, xp  = wfix_t(ξi, e(bi), xp, acr, xis, es, α, σa, na, pv)

    if lU < acr
      if wt <= div(na,2)
        fixtip1!(t0, wt, 0, xp)
      else
        fixtip2!(t0, na - wt + 1, 0, xp)
      end
    end
  # if unfixed node
  else
    _fixrtip!(t0, na)
  end

  if lU < acr
    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, α, σa, σk, ψts, ixf, nep, acr, lU, 
          iρi, na, nn)
    end

    if lU < acr

      # fossilize extant tip
      fossilizefixedtip!(t0)

      # forward simulate fixed tip daughter
      tx, na, nn, acr =
        fossiltip_sim!(t0, tf(bi), λ, μ, ψ, α, σa, σk, ψts, ixf, 
          nep, acr, lU, iρi, na, nn)

      if lU < acr

        llr = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
        setni!(bi, na)       # set new ni

        return t0, llr
      end
    end
  end

  return t0, NaN
end




"""
    fsbi_t(bi ::iBffs,
           xav::Float64,
           xst::Float64,
           ξi ::sTfpe,
           λ  ::Float64,
           μ  ::Float64,
           ψ  ::Vector{Float64},
           α  ::Float64,
           σa ::Float64,
           σk ::Float64,
           ψts::Vector{Float64},
           ix ::Int64,
           xis::Vector{Float64},
           xfs::Vector{Float64},
           es ::Vector{Float64})

Forward simulation for **non-fossil** terminal branch.
"""
function fsbi_t(bi ::iBffs,
                xav::Float64,
                xsd::Float64,
                ξi ::sTfpe,
                λ  ::Float64,
                μ  ::Float64,
                ψ  ::Vector{Float64},
                α  ::Float64,
                σa ::Float64,
                σk ::Float64,
                ψts::Vector{Float64},
                ix ::Int64,
                xis::Vector{Float64},
                xfs::Vector{Float64},
                es ::Vector{Float64}, 
                pv ::Vector{Float64})

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inverse branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  empty!(xis)
  empty!(xfs)
  empty!(es)
  nep = lastindex(ψts) + 1

  t0, na, nn, llr =
    _sim_cfpe_t(e(bi), λ, μ, ψ, xi(ξi), α, σa, σk, ψts, ix, nep, lc, iρi, 
      0, 1, 500, xis, xfs, es)

  if na < 1 || isnan(llr)
    return t0, NaN
  end

  # if fix node
  if ifx(bi)

    # propose trait value (if no uncertainty, then xp = xav)
    xp = xav
    if xsd > 0.0
      xp = rnorm(xav, xsd)
    end

    wt, acr, xp = wfix_pt(ξi, e(bi), xp, 0.0, xis, es, α, σa, na, nac, pv)

    if lU < acr + llr

      if wt <= div(na,2)
        fixtip1!(t0, wt, 0, xp)
      else
        fixtip2!(t0, na - wt + 1, 0, xp)
      end

      setni!(bi, na)    # set new ni
      return t0, llr
    end

  else
    if lU < llr
      _fixrtip!(t0, na)

      setni!(bi, na)    # set new ni
      return t0, llr
    end
  end

  return t0, NaN
end




"""
    wfix_t(ξi ::sTfpe,
           ei ::Float64,
           xav::Float64,
           acr::Float64,
           xis::Vector{Float64},
           es ::Vector{Float64},
           α  ::Float64,
           σa ::Float64,
           na ::Int64, 
           pv ::Vector{Float64})

Choose most likely simulated lineage to fix with respect to the
trait value **without uncertainty** of terminal branches.
"""
function wfix_t(ξi ::sTfpe,
                ei ::Float64,
                xav::Float64,
                acr::Float64,
                xis::Vector{Float64},
                es ::Vector{Float64},
                α  ::Float64,
                σa ::Float64,
                na ::Int64, 
                pv ::Vector{Float64})

  # sample from proposal
  empty!(pv)
  sp = 0.0
  for i in Base.OneTo(na)
    esi = es[i]
    p   = dnorm(xav, xis[i] + α*esi, sqrt(esi)*σa)
    push!(pv, p)
    sp += p
  end

  if iszero(sp)
    return 0, NaN, NaN
  end

  wt = _samplefast(pv, sp, na)

  # extract current `xis` and estimate ratio
  empty!(xis)
  empty!(es)
  nac, xic = _xatt!(ξi, ei, xis, es, 0.0, 0, NaN)

  sc = 0.0
  for i in Base.OneTo(nac)
    esi = es[i]
    p   = dnorm(xav, xis[i] + α*esi, sqrt(esi)*σa)
    sc += p
  end

  # likelihood ratio and acceptance
  acr += log(sp/sc)

  return wt, acr, xav
end




"""
    wfix_pt(ξi ::sTfpe,
            ei ::Float64,
            xav::Float64,
            acr::Float64,
            xis::Vector{Float64},
            es ::Vector{Float64},
            α  ::Float64,
            σa ::Float64,
            na ::Int64,
            nac::Int64,
            pv ::Vector{Float64})

Choose most likely simulated lineage to fix with respect to the
trait value **without uncertainty** of present terminal branches.
"""
function wfix_pt(ξi ::sTfpe,
                 ei ::Float64,
                 xav::Float64,
                 acr::Float64,
                 xis::Vector{Float64},
                 es ::Vector{Float64},
                 α  ::Float64,
                 σa ::Float64,
                 na ::Int64,
                 nac::Int64,
                 pv ::Vector{Float64})

  # sample from proposal
  wt, sp, pp = 0, 0.0, NaN
  if isone(na)
    empty!(pv)
    for i in Base.OneTo(na)
      esi = es[i]
      p   = dnorm(xav, xis[i] + α*esi, sqrt(esi)*σa)
      push!(pv, p)
      sp += p
    end

    if iszero(sp)
      return 0, NaN, NaN
    end

    wt = _samplefast(pv, sp, na)
    pp = pv[wt]
  else
    esi = es[1]
    pp = sp = dnorm(xav, xis[1] + α*esi, sqrt(esi)*σa)
    wt = 1
  end

  # extract current `xis` and estimate ratio
  sc, pc = 0.0, NaN
  if isone(nac)
    empty!(xis)
    empty!(es)
    nac, xic = _xatt!(ξi, ei, xis, es, 0.0, 0, NaN)

    for i in Base.OneTo(nac)
      p   = dnorm(xav, xis[i] + α*esi, sqrt(esi)*σa)
      sc += p
      if xic === xis[i]
        pc = p
      end
    end
  else
    pc = sc = 1.0
  end

  acr += log(pp) - log(pp/sp) + log(pc/sc)

  return wt, acr, xav
end




"""
    fsbi_et(t0  ::sTfpe,
            bi  ::iBffs,
            λ   ::Float64,
            μ   ::Float64,
            ψ   ::Vector{Float64},
            α   ::Float64,
            σa  ::Float64,
            σk  ::Float64,
            ψts ::Vector{Float64},
            ixf ::Int64)

Forward simulation for fossil terminal branch `bi`.
"""
function fsbi_et(t0  ::sTfpe,
                 bi  ::iBffs,
                 λ   ::Float64,
                 μ   ::Float64,
                 ψ   ::Vector{Float64},
                 α   ::Float64,
                 σa  ::Float64,
                 σk  ::Float64,
                 ψts ::Vector{Float64},
                 ixf ::Int64)

  nep = lastindex(ψts) + 1
  lU  = -randexp()            # log-probability
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # if terminal fossil branch
  tx, na, nn, acr =
    fossiltip_sim!(t0, tf(bi), λ, μ, ψ, α, σa, σk, ψts, ixf, nep, 
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
           xav::Float64,
           xst::Float64,
           ξi ::sTfpe,
           ξ1 ::sTfpe,
           λ  ::Float64,
           μ  ::Float64,
           ψ  ::Vector{Float64},
           α  ::Float64,
           σa ::Float64,
           σk ::Float64,
           ψts::Vector{Float64},
           ixi::Int64,
           ixf::Int64,
           xis::Vector{Float64},
           xfs::Vector{Float64},
           es ::Vector{Float64})

Forward simulation for internal branch.
"""
function fsbi_m(bi ::iBffs,
                xav::Float64,
                xst::Float64,
                ξi ::sTfpe,
                ξ1 ::sTfpe,
                λ  ::Float64,
                μ  ::Float64,
                ψ  ::Vector{Float64},
                α  ::Float64,
                σa ::Float64,
                σk ::Float64,
                ψts::Vector{Float64},
                ixi::Int64,
                ixf::Int64,
                xis::Vector{Float64},
                xfs::Vector{Float64},
                es ::Vector{Float64},
                pv ::Vector{Float64})

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  empty!(xis)
  empty!(xfs)
  empty!(es)

  t0, na, nf, nn = 
    _sim_cfpe_i(ti(bi), tf(bi), λ, μ, ψ, xi(ξi), α, σa, σk, ψts, ixi, nep, 
      0, 0, 1, 500, xis, xfs, es)

  if na < 1 || nf > 0 || nn > 499
    return t0, NaN, NaN, NaN
  end

  lU = -randexp() # log-probability

  # add sampling fraction
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # inverse branch sampling fraction
  acr = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  ## choose most likely lineage to fix
  pp = pc = 1.0
  # if fix node
  if ifx(bi)
    # if no uncertainty around trait value
    if iszero(xst)
      wt, acr, xp = wfix_t(ξi, e(bi), xav, acr, xis, es, α, σa, na, pv)
    # if uncertainty around trait value
    else
       xp, wt, pp, pc, acr = 
         wfix_m(ξi, ξ1, e(bi), xav, xst, acr, xfs, α, σa, na, pv)
    end
  # if non-fixed node
  else
    xp, wt, pp, pc, acr = wfix_m(ξi, ξ1, e(bi), acr, xfs, α, σa, na, pv)
  end

  if lU < acr

    # fix the tip
    if wt <= div(na, 2)
      fixtip1!(t0, wt, 0, xp)
    else
      fixtip2!(t0, na - wt + 1, 0, xp)
    end

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, α, σa, σk, ψts, ixf, nep, acr, lU, 
          iρi, na, nn)
    end

    if lU < acr

      # fossilize extant tip
      isfossil(bi) && fossilizefixedtip!(t0)

      # likelihood ratio
      na -= 1
      llr  = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + log(pp/pc)
      setni!(bi, na)   # set new ni
      xi1, xf1, e1   = xi(ξ1), xf(ξ1), e(ξ1)
      dαr  = xi1 - xp
      sσar = 0.5*((xf1 - xp - α*e1)^2 - (xf1 - xi1 - α*e1)^2)/e1
      setxi!(ξ1, xp)   # set new xp for initial x

      return t0, llr, dαr, sσar
    end
  end

  return t0, NaN, NaN, NaN
end




"""
    wfix_m(ξi ::sTfpe,
           ξ1 ::sTfpe,
           ei ::Float64,
           acr::Float64,
           xfs::Vector{Float64},
           α  ::Float64,
           σa ::Float64, 
           na ::Int64,
           pv ::Vector{Float64})

Choose most likely simulated lineage to fix with respect to daughter
for `mid` branches.
"""
function wfix_m(ξi ::sTfpe,
                ξ1 ::sTfpe,
                ei ::Float64,
                acr::Float64,
                xfs::Vector{Float64},
                α  ::Float64,
                σa ::Float64, 
                na ::Int64,
                pv ::Vector{Float64})

  e1 = e(ξ1)
  xf1, sre1σa = xf(ξ1), sqrt(e1)*σa

  empty!(pv)
  sp = 0.0
  for xfi in xfs
    p   = dnorm(xf1, xfi + α*e1, sre1σa)
    push!(pv, p)
    sp += p
  end

  if iszero(sp)
    return NaN, 0, NaN, NaN, NaN
  end

  wt = _samplefast(pv, sp, na)
  pp = pv[wt]
  xp = xfs[wt]

  # extract current xfs and estimate ratio
  empty!(xfs)
  xc, shc = _xatt!(ξi, σa^2, ei, xfs, 0.0, NaN, false)

  sc, pc = 0.0, NaN
  for xfi in xfs
    p   = dnorm(xf1, xfi + α*e1, sre1σa)
    sc += p
    if xc === xfi
      pc = p
    end
  end

  # likelihood ratio and acceptance
  acr += log(sp/sc)

  return xp, wt, pp, pc, acr
end




"""
    wfix_m(ξi ::sTfpe,
           ξ1 ::sTfpe,
           ei ::Float64,
           xav::Float64,
           xst::Float64,
           acr::Float64,
           xfs::Vector{Float64},
           α  ::Float64,
           σa ::Float64, 
           na ::Int64,
           pv ::Vector{Float64})

Choose most likely simulated lineage to fix with respect to the
trait value **with uncertainty** of mid branches.
"""
function wfix_m(ξi ::sTfpe,
                ξ1 ::sTfpe,
                ei ::Float64,
                xav::Float64,
                xst::Float64,
                acr::Float64,
                xfs::Vector{Float64},
                α  ::Float64,
                σa ::Float64, 
                na ::Int64,
                pv ::Vector{Float64})

  # select best from proposal
  e1 = e(ξ1)
  xf1, e1σ2a = xf(ξ1), e1*σa^2

  empty!(pv)
  sp = 0.0
  for xfi in xfs
    p   = duodnorm(xfi, xav, xf1 - α*e1, xst^2, e1σ2a)
    push!(pv, p)
    sp += p
  end

  if iszero(sp)
    return NaN, 0, NaN, NaN, NaN
  end

  wt = _samplefast(pv, sp, na)
  pp = pv[wt]
  xp = xfs[wt]

  # extract current xfs and estimate ratio
  empty!(xfs)
  xc, shc = _xatt!(ξi, σa^2, ei, xfs, 0.0, NaN, false)

  sc, pc = 0.0, NaN
  for xfi in xfs
    p   = duodnorm(xfi, xav, xf1 - α*e1, xst^2, e1σ2a)
    sc += p
    if xc === xfi
      pc = p
    end
  end

  # likelihoods ratio and acceptance
  acr += log(sp/sc)
  pp   = logdnorm(xf1, xp + α*e1, e1σ2a)
  pc   = logdnorm(xf1, xc + α*e1, e1σ2a)

  return xp, wt, pp, pc, acr
end




"""
    fsbi_i(bi ::iBffs,
           ξi ::sTfpe,
           ξ1 ::sTfpe,
           ξ2 ::sTfpe,
           λ  ::Float64,
           μ  ::Float64,
           ψ  ::Vector{Float64},
           α  ::Float64,
           σa ::Float64,
           σk ::Float64,
           ψts::Vector{Float64},
           ixi::Int64,
           ixf::Int64,
           xfs::Vector{Float64},
           xis::Vector{Float64})

Forward simulation for internal branch.
"""
function fsbi_i(bi ::iBffs,
                ξi ::sTfpe,
                ξ1 ::sTfpe,
                ξ2 ::sTfpe,
                λ  ::Float64,
                μ  ::Float64,
                ψ  ::Vector{Float64},
                α  ::Float64,
                σa ::Float64,
                σk ::Float64,
                ψts::Vector{Float64},
                ixi::Int64,
                ixf::Int64,
                xfs::Vector{Float64}, 
                pv ::Vector{Float64})

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  empty!(xfs)

  t0, na, nf, nn = 
    _sim_cfpe_i(ti(bi), tf(bi), λ, μ, ψ, xi(ξi), α, σa, σk, ψts, 
      ixi, nep, 0, 0, 1, 500, xfs)

  if na < 1 || nf > 0 || nn > 499
    return t0, NaN, NaN, NaN, NaN
  end

  lU = -randexp() #log-probability

  # add sampling fraction
  nac = ni(bi)                # current ni
  iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = - Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

  # choose most likely lineage to fix
  wt, xp, xkp, shp, pp, xc, shc, pc, acr = 
    wfix_i(ξi, ξ1, ξ2, e(bi), acr, xfs, α, σa^2, σk^2, na, pv)

  if lU < acr

    # fix the tip
    if wt <= div(na,2)
      fixtip1!(t0, wt, 0, shp)
    else
      fixtip2!(t0, na - wt + 1, 0, shp)
    end

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, α, σa, σk, ψts, ixf, nep, acr, lU, 
          iρi, na, nn)
    end

    if lU < acr
      na  -= 1
      llr  = (na - nac)*(iszero(iρi) ? 0.0 : log(iρi)) + pp - pc
      setni!(bi,  na)   # set new ni

      ξac, ξkc = if shc ξ2, ξ1 else ξ1, ξ2 end
      ξap, ξkp = if shp ξ2, ξ1 else ξ1, ξ2 end

      xikp, xfap, xfkp, eap, ekp = xi(ξkp), xf(ξap), xf(ξkp), e(ξap), e(ξkp)
      xikc, xfac, xfkc, eac, ekc = xi(ξkc), xf(ξac), xf(ξkc), e(ξac), e(ξkc)
      dαr  = (xfap -  xp)  - (xfac -   xc) + 
             (xfkp - xkp) - (xfkc - xikc)

      sσar = 0.5*((xfap - xp   - α*eap)^2/eap - 
                  (xfac - xc   - α*eac)^2/eac +
                  (xfkp - xkp  - α*ekp)^2/ekp - 
                  (xfkc - xikc - α*ekc)^2/ekc)
      sσkr = 0.5*((xp - xkp)^2 - (xc - xikc)^2)
      setxi!(ξap, xp)   # set new xp for initial anagenetic daughter
      setxi!(ξkp, xkp)   # set new xp for initial cladogenetic daughter

      return t0, llr, dαr, sσar, sσkr
    end
  end

  return t0, NaN, NaN, NaN, NaN
end



# """
#     wfix_i(ξi ::sTpe,
#            ξ1 ::sTpe,
#            ξ2 ::sTpe,
#            ei ::Float64,
#            acr::Float64,
#            xfs::Vector{Float64},
#            σa2::Float64,
#            σk2::Float64)

# Choose most likely simulated lineage to fix with respect to daughter
# for bifurcating `i` branches.
# """
# function wfix_i(ξi ::sTfpe,
#                 ξ1 ::sTfpe,
#                 ξ2 ::sTfpe,
#                 ei ::Float64,
#                 acr::Float64,
#                 xfs::Vector{Float64},
#                 α  ::Float64, 
#                 σa2::Float64,
#                 σk2::Float64)

#   xi1, xi2, xf1, xf2, e1, e2 = xi(ξ1), xi(ξ2), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2)

#   ntp = lastindex(xfs)
#   # select one tip at random
#   wt  = fIrand(ntp) + 1
#   xfp = xfs[wt]
#   pk1 = llik_cpe_dyad(xfp, xf2 - α*e2, xf1 - α*e1, e2, e1, σa2, σk2)
#   pk2 = llik_cpe_dyad(xfp, xf1 - α*e1, xf2 - α*e2, e1, e2, σa2, σk2)
#   o12 = exp(pk1 - pk2)  # odds
#   p1  = o12/(1.0 + o12) # probability of d1 cladogenetic
#   shp = rand() < p1
#   lpp = shp ? pk1 - log(p1) : pk2 - log(1.0 - p1)

#   # proposal cladogenetic and likelihood
#   xkp, ll3p = NaN, NaN
#   if shp
#     xkp  = duoprop(xfp, xf1, σk2, e1*σa2)
#     ll3p = llik_cpe_trio(xfp, xkp, xf2 - α*e2, xf1 - α*e1, e2, e1, σa2, σk2)
#   else
#     xkp  = duoprop(xfp, xf2, σk2, e2*σa2)
#     ll3p = llik_cpe_trio(xfp, xkp, xf1 - α*e1, xf2 - α*e2, e1, e2, σa2, σk2)
#   end

#   # extract current xis and estimate ratio
#   empty!(xfs)
#   xfc, shc = _xatt!(ξi, σa2, ei, xfs, 0.0, NaN, false)
#   ntc = Float64(lastindex(xfs))
#   pk1 = llik_cpe_dyad(xfc, xf2 - α*e2, xf1 - α*e1, e2, e1, σa2, σk2)
#   pk2 = llik_cpe_dyad(xfc, xf1 - α*e1, xf2 - α*e2, e1, e2, σa2, σk2)
#   o12 = exp(pk1 - pk2)  # odds
#   p1  = o12/(1.0 + o12) # probability of d1 cladogenetic

#   lpc, ll3c = NaN, NaN
#   if shc
#     lpc  = pk1 - log(p1)
#     ll3c = llik_cpe_trio(xfc, xi1, xf2 - α*e2, xf1 - α*e1, e2, e1, σa2, σk2)
#   else
#     lpc  = pk2 - log(1.0 - p1)
#     ll3c = llik_cpe_trio(xfc, xi2, xf1 - α*e1, xf2 - α*e2, e1, e2, σa2, σk2)
#   end

#   # add to acceptance ratio
#   acr += lpp - lpc + log(Float64(ntp)/ntc)

#   return wt, xfp, xkp, shp, ll3p, xfc, shc, ll3c, acr
# end



"""
    wfix_i(ξi ::sTfpe,
           ξ1 ::sTfpe,
           ξ2 ::sTfpe,
           ei ::Float64,
           acr::Float64,
           xfs::Vector{Float64},
           α  ::Float64, 
           σa2::Float64,
           σk2::Float64,
           na ::Int64,
           pv ::Vector{Float64})

Choose most likely simulated lineage to fix with respect to daughter
for bifurcating `i` branches.
"""
function wfix_i(ξi ::sTfpe,
                ξ1 ::sTfpe,
                ξ2 ::sTfpe,
                ei ::Float64,
                acr::Float64,
                xfs::Vector{Float64},
                α  ::Float64, 
                σa2::Float64,
                σk2::Float64,
                na ::Int64,
                pv ::Vector{Float64})

  xi1, xi2, xf1, xf2, e1, e2 = xi(ξ1), xi(ξ2), xf(ξ1), xf(ξ2), e(ξ1), e(ξ2)

  # select best from proposal
  empty!(pv)
  sp = 0.0
  for xfi in xfs
    pk1 = exp(llik_cpe_dyad(xfi, xf2 - α*e2, xf1 - α*e1, e2, e1, σa2, σk2))
    pk2 = exp(llik_cpe_dyad(xfi, xf1 - α*e1, xf2 - α*e2, e1, e2, σa2, σk2))
    push!(pv, pk1, pk2)
    sp += pk1 + pk2
  end

  if iszero(sp)
    return 0, NaN, NaN, false, NaN, NaN, false, NaN, NaN
  end

  wi  = _samplefast(pv, sp, 2*na)
  pkp = pv[wi]
  pp  = pkp/sp
  shp = isodd(wi)
  wt  = div(wi,2) + (shp ? 1 : 0)
  xp  = xfs[wt]

  # proposal cladogenetic and likelihood
  xkp, ll3p = NaN, NaN
  if shp
    xkp = duoprop(xp, xf1, σk2, e1*σa2)
    ll3p  = llik_cpe_trio(xp, xkp, xf2 - α*e2, xf1 - α*e1, e2, e1, σa2, σk2)
  else
    xkp = duoprop(xp, xf2, σk2, e2*σa2)
    ll3p  = llik_cpe_trio(xp, xkp, xf1 - α*e1, xf2 - α*e2, e1, e2, σa2, σk2)
  end

  # extract current xis and estimate ratio
  empty!(xfs)
  xc, shc = _xatt!(ξi, σa2, ei, xfs, 0.0, NaN, false)

  sc, pc, ll3c = 0.0, NaN, NaN
  for xfi in xfs
    pk1 = exp(llik_cpe_dyad(xfi, xf2 - α*e2, xf1 - α*e1, e2, e1, σa2, σk2))
    pk2 = exp(llik_cpe_dyad(xfi, xf1 - α*e1, xf2 - α*e2, e1, e2, σa2, σk2))
    sc += pk1 + pk2
    if xc === xfi
      if shc
        pc  = pk1
        ll3c = llik_cpe_trio(xfi, xi1, xf2 - α*e2, xf1 - α*e1, e2, e1, σa2, σk2)
      else
        pc  = pk2
        ll3c = llik_cpe_trio(xfi, xi2, xf1 - α*e1, xf2 - α*e2, e1, e2, σa2, σk2)
      end
    end
  end

  pkc = pc
  pc /= sc

  # likelihood ratio and acceptance
  acr += log(pkp/pp) - log(pkc/pc)

  return wt, xp, xkp, shp, ll3p, xc, shc, ll3c, acr
end




"""
    tip_sims!(tree::sTfpe,
              t   ::Float64,
              λ   ::Float64,
              μ   ::Float64,
              ψ   ::Vector{Float64},
              α   ::Float64,
              σa  ::Float64,
              σk  ::Float64,
              ψts ::Vector{Float64},
              ix  ::Int64,
              nep ::Int64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64,
              nn  ::Int64,
              nlim::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTfpe,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Vector{Float64},
                   α   ::Float64,
                   σa  ::Float64,
                   σk  ::Float64,
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
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nn, lr = 
          _sim_cfpe_it(t, λ, μ, ψ, xf(tree), α, σa, σk, ψts, ix, nep, 
            lr, lU, iρi, na-1, nn, 500)

        if isnan(lr) || nn > 499
          return tree, na, nn, NaN
        end

        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        setxf!(tree, xf(stree))
        setsh!(tree, sh(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr = 
        tip_sims!(tree.d1, t, λ, μ, ψ, α, σa, σk, ψts, ix, nep, 
          lr, lU, iρi, na, nn)
      tree.d2, na, nn, lr = 
        tip_sims!(tree.d2, t, λ, μ, ψ, α, σa, σk, ψts, ix, nep, 
          lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::sTfpe,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Vector{Float64},
                   α   ::Float64,
                   σa  ::Float64,
                   σk  ::Float64,
                   ψts ::Vector{Float64},
                   ix  ::Int64,
                   nep ::Int64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64)

Continue simulation until time `t` for the fixed fossil tip in `tree`.
"""
function fossiltip_sim!(tree::sTfpe,
                        t   ::Float64,
                        λ   ::Float64,
                        μ   ::Float64,
                        ψ   ::Vector{Float64},
                        α   ::Float64,
                        σa  ::Float64,
                        σk  ::Float64,
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

      stree, na, nn, lr = 
        _sim_cfpe_it(t, λ, μ, ψ, xf(tree), α, σa, σk, ψts, ix, nep, 
          lr, lU, iρi, na-1, nn, 500)

      if !isfinite(lr) || nn > 499
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr =
        fossiltip_sim!(tree.d1, t, λ, μ, ψ, α, σa, σk, ψts, ix, nep, 
          lr, lU, iρi, na, nn)
    else
      tree.d2, na, nn, lr =
        fossiltip_sim!(tree.d2, t, λ, μ, ψ, α, σa, σk, ψts, ix, nep, 
          lr, lU, iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end



