#=

constant birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cbd(tree    ::sT_label, 
               out_file::String;
               λprior  ::NTuple{2,Float64}     = (1.0, 1.0),
               μprior  ::NTuple{2,Float64}     = (1.0, 1.0),
               niter   ::Int64                 = 1_000,
               nthin   ::Int64                 = 10,
               nburn   ::Int64                 = 200,
               logZ    ::Bool                  = false,
               nitPP   ::Int64                 = 100, 
               nthPP   ::Int64                 = 10,
               K       ::Int64                 = 10,
               ϵi      ::Float64               = 0.4,
               λi      ::Float64               = NaN,
               μi      ::Float64               = NaN,
               pupdp   ::NTuple{3,Float64}     = (0.2,0.2,0.2),
               prints  ::Int64                 = 5,
               tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant pure-birth.
"""
function insane_cbd(tree    ::sT_label, 
                    out_file::String;
                    λprior  ::NTuple{2,Float64}     = (1.0, 1.0),
                    μprior  ::NTuple{2,Float64}     = (1.0, 1.0),
                    niter   ::Int64                 = 1_000,
                    nthin   ::Int64                 = 10,
                    nburn   ::Int64                 = 200,
                    logZ    ::Bool                  = false,
                    nitPP   ::Int64                 = 100, 
                    nthPP   ::Int64                 = 10,
                    K       ::Int64                 = 10,
                    ϵi      ::Float64               = 0.4,
                    λi      ::Float64               = NaN,
                    μi      ::Float64               = NaN,
                    pupdp   ::NTuple{3,Float64}     = (0.2,0.2,0.2),
                    prints  ::Int64                 = 5,
                    tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n  = ntips(tree)

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ)

  # starting parameters
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), ti(idf[1]), ϵi)
  else
    λc, μc = λi, μi
  end

  # make a decoupled tree and fix it
  Ψ = sTbd[]
  sTbd!(Ψ, tree)

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # conditioning functions
  sns = (BitVector(), BitVector(), BitVector())
  snodes! = make_snodes(idf, !iszero(e(tree)), sTbd)
  snodes!(Ψ, sns)
  scond, scond0 = make_scond(idf, !iszero(e(tree)), sTbd)

  @info "Running constant birth-death with forward simulation"

  # adaptive phase
  llc, prc, λc, μc = 
      mcmc_burn_cbd(Ψ, idf, λprior, μprior, nburn, λc, μc,
        pup, prints, sns, snodes!, scond, scond0)

  # mcmc
  R, treev = mcmc_cbd(Ψ, idf, llc, prc, λc, μc, λprior, μprior, niter, nthin, 
    pup, prints, sns, snodes!, scond, scond0)

  if logZ
    # powers
    βs = Vector{Float64}(undef, K-1)
    for i in Base.OneTo(K-1)
      βs[i] = (Float64(i)/Float64(K-1))^3.333333333333333
    end
    reverse!(βs)
    push!(βs, 0.0)

    # marginal likelihood
    PP = power_posterior(Ψ, idf, llc, prc, λc, μc, λprior, μprior, 
      nitPP, nthPP, βs, pup, sns, snodes!, scond, scond0)

    PP[1] = R[:,2]
    reverse!(PP)
    reverse!(βs)

    ps_ml = path_sampling(PP, βs)
    ss_ml = stepping_stone(PP, βs)
  else
    ps_ml = NaN
    ss_ml = NaN
  end

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2))

  write_ssr(R, pardic, out_file)

  return R, treev, ps_ml, ss_ml
end




"""
    mcmc_burn_cbd(Ψ       ::Vector{sTbd},
                  idf     ::Array{iBffs,1},
                  λprior  ::NTuple{2,Float64},
                  μprior  ::NTuple{2,Float64},
                  nburn   ::Int64,
                  λc      ::Float64,
                  μc      ::Float64,
                  pup     ::Array{Int64,1}, 
                  prints  ::Int64,
                  sns     ::NTuple{3, BitVector},
                  snodes! ::Function,
                  scond   ::Function,
                  scond0  ::Function)

Adaptive MCMC phase for da chain for constant birth-death using forward
simulation.
"""
function mcmc_burn_cbd(Ψ       ::Vector{sTbd},
                       idf     ::Array{iBffs,1},
                       λprior  ::NTuple{2,Float64},
                       μprior  ::NTuple{2,Float64},
                       nburn   ::Int64,
                       λc      ::Float64,
                       μc      ::Float64,
                       pup     ::Array{Int64,1}, 
                       prints  ::Int64,
                       sns     ::NTuple{3, BitVector},
                       snodes! ::Function,
                       scond   ::Function,
                       scond0  ::Function)

  el = lastindex(idf)
  L  = treelength(Ψ)     # tree length
  ns = Float64(el-1)/2.0 # number of speciation events
  ne = 0.0               # number of extinction events

  # likelihood
  llc = llik_cbd(Ψ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf)
  prc = logdgamma(λc, λprior[1], λprior[2]) + 
        logdgamma(μc, μprior[1], μprior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc = 
          update_λ!(llc, prc, λc, ns, L, μc, sns, λprior, scond, 1.0)

      # μ proposal
      elseif p === 2

        llc, prc, μc = 
          update_μ!(llc, prc, μc, ne, L, λc, sns, μprior, scond, 1.0)

      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L = update_fs!(bix, Ψ, idf, llc, λc, μc, ns, ne, L, sns,
                           snodes!, scond0, 1.0)
      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc
end




"""
    mcmc_cbd(Ψ      ::Vector{sTbd},
             idf    ::Array{iBffs,1},
             llc    ::Float64,
             prc    ::Float64,
             λc     ::Float64,
             μc     ::Float64,
             λprior ::NTuple{2,Float64},
             μprior ::NTuple{2,Float64},
             niter  ::Int64,
             nthin  ::Int64,
             pup    ::Array{Int64,1}, 
             prints ::Int64,
             sns    ::NTuple{3, BitVector},
             snodes!::Function,
             scond  ::Function,
             scond0 ::Function)

MCMC da chain for constant birth-death using forward simulation.
"""
function mcmc_cbd(Ψ      ::Vector{sTbd},
                  idf    ::Array{iBffs,1},
                  llc    ::Float64,
                  prc    ::Float64,
                  λc     ::Float64,
                  μc     ::Float64,
                  λprior ::NTuple{2,Float64},
                  μprior ::NTuple{2,Float64},
                  niter  ::Int64,
                  nthin  ::Int64,
                  pup    ::Array{Int64,1}, 
                  prints ::Int64,
                  sns    ::NTuple{3, BitVector},
                  snodes!::Function,
                  scond  ::Function,
                  scond0 ::Function)

  el = lastindex(idf)
  ns = Float64(nnodesinternal(Ψ))
  ne = Float64(ntipsextinct(Ψ))
  L  = treelength(Ψ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 5)

  # make tree vector
  treev  = sTbd[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc = 
          update_λ!(llc, prc, λc, ns, L, μc, sns, λprior, scond, 1.0)

        # llci = llik_cbd(Ψ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # μ proposal
      elseif p === 2

        llc, prc, μc = 
          update_μ!(llc, prc, μc, ne, L, λc, sns, μprior, scond, 1.0)

        # llci = llik_cbd(Ψ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)
        llc, ns, ne, L = update_fs!(bix, Ψ, idf, llc, λc, μc, ns, ne, L, sns,
                           snodes!, scond0, 1.0)

        # llci = llik_cbd(Ψ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf)
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
        push!(treev, couple(deepcopy(Ψ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev
end




"""
    power_posterior(Ψ      ::Vector{sTbd},
                    idf    ::Array{iBffs,1},
                    llc    ::Float64,
                    prc    ::Float64,
                    λc     ::Float64,
                    μc     ::Float64,
                    λprior ::Float64,
                    μprior ::Float64,
                    nitPP  ::Int64,
                    nthPP  ::Int64,
                    βs     ::Vector{Float64},
                    λtn    ::Float64,
                    μtn    ::Float64, 
                    pup    ::Array{Int64,1}, 
                    sns    ::NTuple{3, BitVector},
                    snodes!::Function,
                    scond  ::Function,
                    scond0 ::Function)

MCMC da chain for constant birth-death using forward simulation.
"""
function power_posterior(Ψ      ::Vector{sTbd},
                         idf    ::Array{iBffs,1},
                         llc    ::Float64,
                         prc    ::Float64,
                         λc     ::Float64,
                         μc     ::Float64,
                         λprior ::NTuple{2,Float64},
                         μprior ::NTuple{2,Float64},
                         nitPP  ::Int64,
                         nthPP  ::Int64,
                         βs     ::Vector{Float64},
                         pup    ::Array{Int64,1}, 
                         sns    ::NTuple{3, BitVector},
                         snodes!::Function,
                         scond  ::Function,
                         scond0 ::Function)

  K = lastindex(βs)

  # make log-likelihood table per power
  nlg = fld(nitPP, nthPP)
  PP  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

  el = lastindex(idf)
  ns = Float64(nnodesinternal(Ψ))
  ne = Float64(ntipsextinct(Ψ))
  L  = treelength(Ψ)

  for k in 2:K

    βi  = βs[k]
    llc = βi * (llik_cbd(Ψ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf))

    # logging
    lth, lit = 0, 0

    for it in Base.OneTo(nitPP)

      shuffle!(pup)

      for p in pup

        # λ proposal
        if p === 1

          llc, prc, λc = 
            update_λ!(llc, prc, λc, ns, L, μc, sns, λprior, scond, βi)

        # μ proposal
        elseif p === 2

          llc, prc, μc = 
            update_μ!(llc, prc, μc, ne, L, λc, sns, μprior, scond, βi)

        # forward simulation proposal proposal
        else 

          if k < K
            bix = ceil(Int64,rand()*el)
            llc, ns, ne, L = update_fs!(bix, Ψ, idf, llc, λc, μc, ns, ne, L, sns,
                               snodes!, scond0, βi)
          end
        end
      end

      # log log-likelihood
      lth += 1
      if lth === nthPP
        lit += 1
        PP[k][lit] = llik_cbd(Ψ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf)
        lth = 0
      end
    end

    @info string(βi," power done")
  end

  return PP
end




"""
    update_fs!(bix  ::Int64,
               Ψ    ::Vector{sTbd},
               idf  ::Vector{iBffs},
               llc  ::Float64,
               λ    ::Float64, 
               μ    ::Float64,
               ns   ::Float64,
               ne   ::Float64,
               L    ::Float64,
               scond::Function,
               pow  ::Float64)

Forward simulation proposal function for constant birth-death.
"""
function update_fs!(bix    ::Int64,
                    Ψ      ::Vector{sTbd},
                    idf    ::Vector{iBffs},
                    llc    ::Float64,
                    λ      ::Float64, 
                    μ      ::Float64,
                    ns     ::Float64,
                    ne     ::Float64,
                    L      ::Float64,
                    sns    ::NTuple{3, BitVector},
                    snodes!::Function,
                    scond0 ::Function,
                    pow    ::Float64)

  bi = idf[bix]

  # forward simulate an internal branch
  ψp, np, ntp = fsbi(bi, λ, μ, 1_000)

  itb = it(bi) # is it terminal
  ρbi = ρi(bi) # get branch sampling fraction
  nc  = ni(bi) # current ni
  ntc = nt(bi) # current nt

  if ntp > 0

    # current tree
    ψc  = Ψ[bix]

    # if terminal branch
    if itb
      llr = pow * (log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc)))
      acr = llr
    else
      np  -= 1
      llr = pow * log((1.0 - ρbi)^(np - nc))
      acr = llr + log(Float64(ntp)/Float64(ntc))
    end

    # MH ratio
    if -randexp() < acr

      # if survival conditioned
      scn = ((iszero(pa(bi)) && e(bi) > 0.0)) || 
             (isone(pa(bi)) && iszero(e(Ψ[1])))
      if scn
          llr += pow * (scond0(ψp, λ, μ, itb) - scond0(ψc, λ, μ, itb))
      end

      # update ns, ne & L
      ns += Float64(nnodesinternal(ψp) - nnodesinternal(ψc))
      ne += Float64(ntipsextinct(ψp)   - ntipsextinct(ψc))
      L  += treelength(ψp)             - treelength(ψc)

      # likelihood ratio
      llr += pow * (llik_cbd(ψp, λ, μ) - llik_cbd(ψc, λ, μ))

      Ψ[bix] = ψp     # set new decoupled tree
      llc += llr      # set new likelihood
      if scn
        snodes!(Ψ, sns) # set new sns
      end
      setni!(bi, np)  # set new ni
      setnt!(bi, ntp) # set new nt
    end
  end

  return llc, ns, ne, L
end





"""
    fsbi(bi::iBffs, λ::Float64, μ::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi(bi::iBffs, λ::Float64, μ::Float64, ntry::Int64)

  # times
  tfb = tf(bi)

  ext = 0
  # condition on non-extinction (helps in mixing)
  while ext < ntry 
    ext += 1

    # forward simulation during branch length
    t0, na = sim_cbd(e(bi), λ, μ, 0)

    nat = na

    if isone(na)
      fixalive!(t0)

      return t0, na, nat
    elseif na > 1
      # fix random tip
      fixrtip!(t0)

      if !it(bi)
        # add tips until the present
        tx, na = tip_sims!(t0, tfb, λ, μ, na)
      end

      return t0, na, nat
    end
  end

  return sTbd(), 0, 0
end




"""
    tip_sims!(tree::sTbd, t::Float64, λ::Float64, μ::Float64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::sTbd, t::Float64, λ::Float64, μ::Float64, na::Int64)

  if istip(tree) 
    if !isfix(tree) && isalive(tree)

      # simulate
      stree, na = sim_cbd(t, λ, μ, na-1)

      # merge to current tip
      sete!(tree, e(tree) + e(stree))
      setproperty!(tree, :iμ, isextinct(stree))
      if isdefined(stree, :d1)
        tree.d1 = stree.d1
        tree.d2 = stree.d2
      end
    end
  else
    tree.d1, na = tip_sims!(tree.d1, t, λ, μ, na)
    tree.d2, na = tip_sims!(tree.d2, t, λ, μ, na)
  end

  return tree, na
end





"""
    update_λ!(psi   ::Vector{sTpb},
              llc   ::Float64,
              prc   ::Float64,
              λc    ::Float64,
              lac   ::Float64,
              λtn   ::Float64,
              ns    ::Float64,
              L     ::Float64,
              λprior::Float64,
              pow   ::Float64)

`λ` proposal function for constant pure-birth in adaptive phase.
"""
function update_λ!(llc   ::Float64,
                   prc   ::Float64,
                   λc    ::Float64,
                   ns    ::Float64,
                   L     ::Float64,
                   μc    ::Float64,
                   sns   ::NTuple{3,BitVector},
                   λprior::NTuple{2,Float64},
                   scond ::Function,
                   pow   ::Float64)

  λp  = randgamma(λprior[1] + pow * ns, λprior[2] + pow * L)
  llr = pow * (scond(λp, μc, sns) - scond(λc, μc, sns))

  if -randexp() < llr
    llc    += pow * (ns*log(λp/λc) + L*(λc - λp)) + llr
    prc    += llrdgamma(λp, λc, λprior[1], λprior[2])
    λc      = λp
  end

  return llc, prc, λc
end




"""
    update_μ!(
            llc   ::Float64,
              prc   ::Float64,
              μc    ::Float64,
              lac   ::Array{Float64,1},
              μtn   ::Float64,
              λc    ::Float64,
              ne    ::Float64,
              L     ::Float64,
              sns   ::NTuple{3,BitVector},
              μprior::Float64,
              scond ::Function,
              pow   ::Float64)

`μ` proposal function for constant birth-death in adaptive phase.
"""
function update_μ!(llc   ::Float64,
                   prc   ::Float64,
                   μc    ::Float64,
                   ne    ::Float64,
                   L     ::Float64,
                   λc    ::Float64,
                   sns   ::NTuple{3,BitVector},
                   μprior::NTuple{2,Float64},
                   scond ::Function,
                   pow   ::Float64)

  μp  = randgamma(μprior[1] + pow * ne, μprior[2] + pow * L)
  llr = pow * (scond(λc, μp, sns) - scond(λc, μc, sns))

  if -randexp() < llr
    llc    += pow * (ne*log(μp/μc) + L*(μc - μp)) + llr
    prc    += llrdgamma(μp, μc, μprior[1], μprior[2])
    μc      = μp
  end

  return llc, prc, μc 
end



