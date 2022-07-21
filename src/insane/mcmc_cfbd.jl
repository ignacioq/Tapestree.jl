#=

constant fossilized birth-death MCMC using forward simulation

Jérémy Andréoletti
Adapted from birth-death MCMC by Ignacio Quintero Mächler

v(°-°v)

Created 07 10 2021
=#




"""
    insane_cfbd(tree    ::sTf_label,
                out_file::String;
                λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                niter   ::Int64                 = 1_000,
                nthin   ::Int64                 = 10,
                nburn   ::Int64                 = 200,
                ϵi      ::Float64               = 0.4,
                λi      ::Float64               = NaN,
                μi      ::Float64               = NaN,
                ψi      ::Float64               = NaN,
                pupdp   ::NTuple{4,Float64}     = (0.2,0.2,0.2,0.2),
                prints  ::Int64                 = 5,
                tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant fossilized birth-death.
"""
function insane_cfbd(tree    ::sTf_label,
                     out_file::String;
                     λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     ψ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                     ψ_epoch ::Vector{Float64}       = Float64[],
                     niter   ::Int64                 = 1_000,
                     nthin   ::Int64                 = 10,
                     nburn   ::Int64                 = 200,
                     ϵi      ::Float64               = 0.4,
                     λi      ::Float64               = NaN,
                     μi      ::Float64               = NaN,
                     ψi      ::Float64               = NaN,
                     pupdp   ::NTuple{4,Float64}     = (0.01, 0.01, 0.01, 0.1),
                     survival::Bool                  = true,
                     prints  ::Int64                 = 5,
                     tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n   = ntips(tree)
  th  = treeheight(tree)
  nep = lastindex(ψ_epoch) + 1

  # set tips sampling fraction
  if isone(length(tρ))
    tl  = tiplabels(tree)
    tρu = tρ[""]
    tρ  = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ)

  # starting parameters
  if isnan(λi) || isnan(μi) || isnan(ψi)
    # if only one tip
    if isone(n)
      λc = prod(λ_prior)
      μc = prod(μ_prior)
    else
      λc, μc = moments(Float64(n), th, ϵi)
    end
    # if no sampled fossil
    if iszero(nfossils(tree))
      ψc = prod(ψ_prior)
    else
      ψc = Float64(nfossils(tree))/treelength(tree)
    end
  else
    λc, μc, ψc = λi, μi, ψi
  end

  # make ψ vector
  ψc = fill(ψc, nep)

  # survival conditioning
  if survival && ntipsalive(tree) > 0
    # if crown conditioning
    if def1(tree) && def2(tree) &&
       ntipsalive(tree.d1) > 0 && ntipsalive(tree.d2) > 0
      crown = 1
    # if crown conditioning
    else
      crown = 0
    end
  else
    crown = 2
  end

  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 5_000, crown)

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, sTfbd)

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

  # parameter updates (1: λ, 2: μ, 3: ψ, 4: forward simulation)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running constant fossilized birth-death"

  # adaptive phase
  llc, prc, λc, μc, ψc, mc =
     mcmc_burn_cfbd(Ξ, idf, λ_prior, μ_prior, ψ_prior, ψ_epoch, nburn,
        λc, μc, ψc, mc, th, crown, bst, eixi, eixf, pup, prints)

  # mcmc
  r, treev, λc, μc, ψc =
    mcmc_cfbd(Ξ, idf, llc, prc, λc, μc, ψc, mc, λ_prior, μ_prior, ψ_prior, 
      ψ_epoch, th, crown, bst, eixi, eixf, niter, nthin,  pup, prints)

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2))
  merge!(pardic, 
    Dict("psi"*(isone(nep) ? "" : string("_",i)) => 2+i 
           for i in Base.OneTo(nep)))

  write_ssr(r, pardic, out_file)

  return r, treev
end




"""
    mcmc_burn_cfbd(Ξ      ::Vector{sTfbd},
                   idf    ::Array{iBffs,1},
                   λ_prior::NTuple{2,Float64},
                   μ_prior::NTuple{2,Float64},
                   ψ_prior::NTuple{2,Float64},
                   nburn  ::Int64,
                   λc     ::Float64,
                   μc     ::Float64,
                   ψc     ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   crown   ::Bool,
                   pup    ::Array{Int64,1},
                   prints ::Int64)

Adaptive MCMC phase for da chain for constant fossilized birth-death using
forward simulation.
"""
function mcmc_burn_cfbd(Ξ      ::Vector{sTfbd},
                        idf    ::Array{iBffs,1},
                        λ_prior::NTuple{2,Float64},
                        μ_prior::NTuple{2,Float64},
                        ψ_prior::NTuple{2,Float64},
                        ψ_epoch::Vector{Float64},
                        nburn  ::Int64,
                        λc     ::Float64,
                        μc     ::Float64,
                        ψc     ::Vector{Float64},
                        mc     ::Float64,
                        th     ::Float64,
                        crown  ::Int64,
                        bst    ::Vector{Float64},
                        eixi   ::Vector{Int64},
                        eixf   ::Vector{Int64},
                        pup    ::Array{Int64,1},
                        prints ::Int64)

  el = lastindex(idf)                        # number of branches
  L  = treelength(Ξ, ψ_epoch, bst, eixi)     # tree length
  nf = nfossils(idf, ψ_epoch)                # number of fossilization events per epoch
  ns = nnodesbifurcation(Ξ) + Float64(crown) # number of speciation events
  ne = Float64(ntipsextinct(Ξ))              # number of extinction events

  # likelihood
  llc = llik_cfbd(Ξ, λc, μc, ψc, ψ_epoch, bst, eixi) + log(mc) + prob_ρ(idf)

  prc = logdgamma(λc, λ_prior[1], λ_prior[2])       +
        logdgamma(μc, μ_prior[1], μ_prior[2])       +
        sum(logdgamma.(ψc, ψ_prior[1], ψ_prior[2]))

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, sum(L), μc, mc, th, crown, λ_prior)

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, sum(L), λc, mc, th, crown, μ_prior)

      # ψ proposal
      elseif p === 3

        llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L =
          update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ψ_epoch, ns, ne, L, 
            eixi, eixf)

      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, ψc, mc
end




"""
    mcmc_cfbd(Ξ      ::Vector{sTfbd},
              idf    ::Array{iBffs,1},
              llc    ::Float64,
              prc    ::Float64,
              λc     ::Float64,
              μc     ::Float64,
              ψc     ::Float64,
              λ_prior::NTuple{2,Float64},
              μ_prior::NTuple{2,Float64},
              ψ_prior::NTuple{2,Float64},
              mc     ::Float64,
              th     ::Float64,
              crown   ::Bool,
              niter  ::Int64,
              nthin  ::Int64,
              pup    ::Array{Int64,1},
              prints ::Int64)

MCMC da chain for constant fossilized birth-death using forward simulation.
"""
function mcmc_cfbd(Ξ      ::Vector{sTfbd},
                   idf    ::Array{iBffs,1},
                   llc    ::Float64,
                   prc    ::Float64,
                   λc     ::Float64,
                   μc     ::Float64,
                   ψc     ::Vector{Float64},
                   mc     ::Float64,
                   λ_prior::NTuple{2,Float64},
                   μ_prior::NTuple{2,Float64},
                   ψ_prior::NTuple{2,Float64},
                   ψ_epoch::Vector{Float64},
                   th     ::Float64,
                   crown  ::Int64,
                   bst    ::Vector{Float64},
                   eixi   ::Vector{Int64},
                   eixf   ::Vector{Int64},
                   niter  ::Int64,
                   nthin  ::Int64,
                   pup    ::Array{Int64,1},
                   prints ::Int64)

  el  = lastindex(idf)                        # number of branches
  L   = treelength(Ξ, ψ_epoch, bst, eixi)     # tree length
  nf  = nfossils(idf, ψ_epoch)                # number of fossilization events per epoch
  ns  = nnodesbifurcation(Ξ) + Float64(crown) # number of speciation events
  ne  = Float64(ntipsextinct(Ξ))              # number of extinction events
  nep = lastindex(ψc)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  lastindex(ψc)
  R = Array{Float64,2}(undef, nlogs, 5 + nep)

  # make tree vector
  treev  = sTfbd[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, sum(L), μc, mc, th, crown, λ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, ψ_epoch, bst, eixi) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, sum(L), λc, mc, th, crown, μ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, ψ_epoch, bst, eixi) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # ψ proposal
      elseif p === 3

        llc, prc = update_ψ!(llc, prc, ψc, nf, L, ψ_prior)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, ψ_epoch, bst, eixi) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L =
          update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ψ_epoch, ns, ne, L, 
            eixi, eixf)

        # llci = llik_cfbd(Ξ, λc, μc, ψc, ψ_epoch, bst, eixi) + log(mc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end
      end
    end

    # log parameters
    lthin += 1
    if lthin == nthin

      lit += 1
      @inbounds begin
        R[lit,1] = Float64(lit)
        R[lit,2] = llc
        R[lit,3] = prc
        R[lit,4] = λc
        R[lit,5] = μc
        @avx for i in Base.OneTo(nep)
          R[lit,5 + i] = ψc[i]
        end
        push!(treev, couple(copy_Ξ(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev, λc, μc, ψc
end




"""
    update_fs!(bix::Int64,
               Ξ  ::Vector{sTfbd},
               idf::Vector{iBffs},
               llc::Float64,
               λ  ::Float64,
               μ  ::Float64,
               ψ  ::Float64,
               ns ::Float64,
               ne ::Float64,
               L  ::Float64)

Forward simulation proposal function for constant fossilized birth-death.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{sTfbd},
                    idf ::Vector{iBffs},
                    llc ::Float64,
                    λ   ::Float64,
                    μ   ::Float64,
                    ψ   ::Vector{Float64},
                    ψts ::Vector{Float64},
                    ns  ::Float64,
                    ne  ::Float64,
                    L   ::Vector{Float64},
                    eixi::Vector{Int64},
                    eixf::Vector{Int64})

  bi  = idf[bix]
  ixi = eixi[bix]

  if isfossil(bi)
    ixf = eixf[bix]
    ξp, llr = fsbi_f(bi, λ, μ, ψ, ψts, ixi, ixf)

    # if terminal but not successful proposal, update extinct
    if it(bi) && !isfinite(llr)
      ξp, llr = fsbi_et(sTfbd_wofe(Ξ[bix]), bi, λ, μ, ψ, ψts, ixf)
    end
  else
    if it(bi)
      ξp, llr = fsbi_t(bi, λ, μ, ψ, ψts, ixi)
    else
      ξp, llr = fsbi_i(bi, λ, μ, ψ, ψts, ixi, eixf[bix])
    end
  end

  if isfinite(llr)
    ξc  = Ξ[bix]
    tii = ti(bi)

    nep = lastindex(ψts) + 1
    # update llc, ns, ne & L
    llc += llik_cfbd(ξp, λ, μ, ψ, tii, ψts, ixi, nep) - 
           llik_cfbd(ξc, λ, μ, ψ, tii, ψts, ixi, nep) + llr

    ns  += Float64(nnodesbifurcation(ξp) - nnodesbifurcation(ξc))
    ne  += Float64(ntipsextinct(ξp)      - ntipsextinct(ξc))

    # update tree lengths
    Lc = zeros(nep)
    _treelength!(ξc, tii, Lc, ψts, ixi, nep)
    _treelength!(ξp, tii, L,  ψts, ixi, nep)
    @avx for i in Base.OneTo(nep)
      L[i] -= Lc[i]
    end

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, ns, ne, L
end




"""
    fsbi_t(bi::iBffs,
           λ::Float64,
           μ::Float64,
           ψ   ::Vector{Float64},
           ψts ::Vector{Float64},
           eixi::Vector{Float64})

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs,
                λ::Float64,
                μ::Float64,
                ψ   ::Vector{Float64},
                ψts ::Vector{Float64},
                ix  ::Int64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nn, llr =
    _sim_cfbd_t(e(bi), λ, μ, ψ, ψts, ix, nep, lc, lU, Iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr)

    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, -Inf
  end
end




"""
    fsbi_f(bi ::iBffs,
           λ  ::Float64,
           μ  ::Float64,
           ψ  ::Vector{Float64},
           ψts::Vector{Float64},
           ixi::Int64,
           ixf::Int64)

Forward simulation for fossil branch `bi`.
"""
function fsbi_f(bi ::iBffs,
                λ  ::Float64,
                μ  ::Float64,
                ψ  ::Vector{Float64},
                ψts::Vector{Float64},
                ixi::Int64,
                ixf::Int64)

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nf, nn = 
    _sim_cfbd_i(ti(bi), tf(bi), λ, μ, ψ, ψts, ixi, nep, 0, 0, 1, 1_000)

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

    _fixrtip!(t0, na)

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, ψts, ixf, acr, lU, Iρi, na, nn)
    end

    if lU < acr
      # fossilize extant tip
      fossilizefixedtip!(t0)

      if it(bi)
        tx, na, nn, acr =
          fossiltip_sim!(t0, tf(bi), λ, μ, ψ, ψts, ixf, acr, lU, Iρi, na, nn)
      else
        na -= 1
      end

      if lU < acr

        llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
        setnt!(bi, ntp)                # set new nt
        setni!(bi, na)                 # set new ni

        return t0, llr
      end
    end
  end

  return t0, NaN
end




"""
    fsbi_et(t0 ::sTfbd,
            bi ::iBffs,
            λ  ::Float64,
            μ  ::Float64,
            ψ  ::Vector{Float64},
            ψts::Vector{Float64},
            ixf::Int64)

Forward simulation for extinct tip in terminal fossil branch.
"""
function fsbi_et(t0 ::sTfbd,
                 bi ::iBffs,
                 λ  ::Float64,
                 μ  ::Float64,
                 ψ  ::Vector{Float64},
                 ψts::Vector{Float64},
                 ixf::Int64)

  lU  = -randexp()            # log-probability
  nac = ni(bi)                # current ni
  Iρi = (1.0 - ρi(bi))        # branch sampling fraction
  acr = Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  tx, na, nn, acr =
    fossiltip_sim!(t0, tf(bi), λ, μ, ψ, ψts, ixf, acr, lU, Iρi, 1, 1)

  if lU < acr

    llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
    setni!(bi, na)                 # set new ni

    return t0, llr
  end

  return t0, NaN
end




"""
    fsbi_i(bi::iBffs,
           λ::Float64,
           μ::Float64,
           ψ  ::Vector{Float64},
           ψts::Vector{Float64},
           ixi::Int64,
           ixf::Int64)

Forward simulation for internal branch `bi`.
"""
function fsbi_i(bi::iBffs,
                λ::Float64,
                μ::Float64,
                ψ  ::Vector{Float64},
                ψts::Vector{Float64},
                ixi::Int64,
                ixf::Int64)

  # forward simulation during branch length
  nep = lastindex(ψts) + 1
  t0, na, nf, nn = 
    _sim_cfbd_i(ti(bi), tf(bi), λ, μ, ψ, ψts, ixi, nep, 0, 0, 1, 1_000)

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

    _fixrtip!(t0, na)

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr =
        tip_sims!(t0, tf(bi), λ, μ, ψ, ψts, ixf, acr, lU, Iρi, na, nn)
    end

    if lU < acr

      na -= 1
      llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setnt!(bi, ntp)                # set new nt
      setni!(bi, na)                 # set new ni

      return t0, llr
    end
  end

  return t0, NaN
end




"""
    tip_sims!(tree::sTfbd,
              t   ::Float64,
              λ   ::Float64,
              μ   ::Float64,
              ψ   ::Vector{Float64},
              ψts ::Vector{Float64},
              ix  ::Int64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64,
              nn  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTfbd,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Vector{Float64},
                   ψts ::Vector{Float64},
                   ix  ::Int64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

  if lU < lr && nn < 1_000

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        nep = lastindex(ψts) + 1
        stree, na, nn, lr =
          _sim_cfbd_it(t, λ, μ, ψ, ψts, ix, nep, lr, lU, Iρi, na-1, nn, 1_000)

        if isnan(lr) || nn > 999
          return tree, na, nn, NaN
        end

        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        if def1(stree)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr =
        tip_sims!(tree.d1, t, λ, μ, ψ, ψts, ix, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr =
        tip_sims!(tree.d2, t, λ, μ, ψ, ψts, ix, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    fossiltip_sim!(tree::sTfbd,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64)

Continue simulation until time `t` for the fixed tip in `tree`.
"""
function fossiltip_sim!(tree::sTfbd,
                        t   ::Float64,
                        λ   ::Float64,
                        μ   ::Float64,
                        ψ   ::Vector{Float64},
                        ψts ::Vector{Float64},
                        ix  ::Int64,
                        lr  ::Float64,
                        lU  ::Float64,
                        Iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64)

  if isfinite(lr) && nn < 1_000

    if istip(tree)

      nep = lastindex(ψts) + 1
      stree, na, nn, lr =
        _sim_cfbd_it(t, λ, μ, ψ, ψts, ix, nep, lr, lU, Iρi, na-1, nn, 1_000)

      if isnan(lr) || nn > 999
        return tree, na, nn, NaN
      end

      # merge to current tip
      tree.d1 = stree
    elseif isfix(tree.d1)
      tree.d1, na, nn, lr =
        fossiltip_sim!(tree.d1, t, λ, μ, ψ, ψts, ix, lr, lU, Iρi, na, nn)
    else
      tree.d2, na, nn, lr =
        fossiltip_sim!(tree.d2, t, λ, μ, ψ, ψts, ix, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_ψ!(llc    ::Float64,
              prc    ::Float64,
              ψc     ::Float64,
              nf     ::Float64,
              L      ::Float64,
              ψ_prior::NTuple{2,Float64})

Gibbs sampling of `ψ` for constant fossilized birth-death.
"""
function update_ψ!(llc    ::Float64,
                   prc    ::Float64,
                   ψc     ::Vector{Float64},
                   nf     ::Vector{Float64},
                   L      ::Vector{Float64},
                   ψ_prior::NTuple{2,Float64})

  for i in Base.OneTo(lastindex(ψc))
    ψci   = ψc[i]
    ψp    = randgamma(ψ_prior[1] + nf[i], ψ_prior[2] + L[i])
    llc  += nf[i] * log(ψp/ψci) + L[i] * (ψci - ψp)
    prc  += llrdgamma(ψp, ψci, ψ_prior[1], ψ_prior[2])
    ψc[i] = ψp
  end

  return llc, prc
end




"""
    update_ψ!(llc    ::Float64,
              prc    ::Float64,
              ψc     ::Float64,
              nf     ::Float64,
              L      ::Float64,
              ψ_prior::NTuple{2,Float64})

Gibbs sampling of `ψ` for constant fossilized birth-death.
"""
function update_ψ!(llc    ::Float64,
                   prc    ::Float64,
                   ψc     ::Float64,
                   nf     ::Float64,
                   L      ::Float64,
                   ψ_prior::NTuple{2,Float64})

  ψp  = randgamma(ψ_prior[1] + nf, ψ_prior[2] + L)

  llc += nf * log(ψp/ψc) + L * (ψc - ψp)
  prc += llrdgamma(ψp, ψc, ψ_prior[1], ψ_prior[2])

  return llc, prc, ψp
end



