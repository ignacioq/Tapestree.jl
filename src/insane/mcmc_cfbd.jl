#=

constant fossilized birth-death MCMC using graft and prune

Jérémy Andréoletti
Adapted from birth-death MCMC by Ignacio Quintero Mächler

v(°-°v)

Created 07 10 2021
=#




"""
    insane_cfbd(tree     ::sTf_label, 
                out_file ::String;
                λ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                μ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                λmμ_prior::NTuple{2,Float64}     = (1.0, 1.0),
                ψ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                niter    ::Int64                 = 1_000,
                nthin    ::Int64                 = 10,
                nburn    ::Int64                 = 200,
                logZ     ::Bool                  = false,
                nitpp    ::Int64                 = 100, 
                nthpp    ::Int64                 = 10,
                K        ::Int64                 = 10,
                ϵi       ::Float64               = 0.4,
                λi       ::Float64               = NaN,
                μi       ::Float64               = NaN,
                ψi       ::Float64               = NaN,
                pupdp    ::NTuple{4,Float64}     = (0.2,0.2,0.2,0.2),
                prints   ::Int64                 = 5,
                tρ       ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant fossilized birth-death.
"""
function insane_cfbd(tree     ::sTf_label, 
                     out_file ::String;
                     λ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                     μ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                     λmμ_prior::NTuple{2,Float64}     = (NaN, NaN),
                     ψ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
                     niter    ::Int64                 = 1_000,
                     nthin    ::Int64                 = 10,
                     nburn    ::Int64                 = 200,
                     marginal ::Bool                  = false,
                     nitpp    ::Int64                 = 100, 
                     nthpp    ::Int64                 = 10,
                     K        ::Int64                 = 11,
                     ϵi       ::Float64               = 0.4,
                     λi       ::Float64               = NaN,
                     μi       ::Float64               = NaN,
                     ψi       ::Float64               = NaN,
                     pupdp    ::NTuple{4,Float64}     = (0.2,0.2,0.2,0.2),
                     prints   ::Int64                 = 5,
                     tρ       ::Dict{String, Float64} = Dict("" => 1.0))

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
  if isnan(λi) && isnan(μi) && isnan(ψi)
    # if only one tip
    if isone(n)
      λc = λ_prior
      μc = isnan(μ_prior[1]) ? λ_prior-λmμ_prior : μ_prior
    else
      λc, μc = moments(Float64(n), ti(idf[1]), ϵi)
    end
    # if no sampled fossil
    if iszero(nfossils(tree))
      ψc = ψ_prior
    else
      ψc = nfossils(tree)/treelength(tree)
    end
  else
    λc, μc, ψc = λi, μi, ψi
  end

  # make a decoupled tree and fix it
  Ξ = sTfbd[]
  sTfbd!(Ξ, tree)

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # conditioning functions
  sns = (BitVector(), BitVector(), BitVector())
  snodes! = make_snodes(idf, !iszero(e(tree)), sTfbd)
  snodes!(Ξ, sns)
  scond, scond0 = make_scond(idf, !iszero(e(tree)), sTfbd)

  @info "Running constant fossilized birth-death with forward simulation"

  # adaptive phase
  @info "MCMC - Adaptive phase"
  llc, prc, λc, μc = mcmc_burn_cfbd(Ξ, idf, λ_prior, μ_prior, λmμ_prior, 
                                    ψ_prior, nburn, λc, μc, ψc, pup, prints, 
                                    sns, snodes!, scond, scond0)

  # mcmc
  @info "MCMC - Sampling iterations"
  r, treev, λc, μc = mcmc_cfbd(Ξ, idf, llc, prc, λc, μc, ψc, λ_prior, μ_prior, 
                               λmμ_prior, ψ_prior, niter, nthin, pup, prints, 
                               sns, snodes!, scond, scond0)

  pardic = Dict(("lambda"      => 1),
                ("mu"          => 2),
                ("psi"         => 3))

  write_ssr(r, pardic, out_file)

  if marginal

    #= # reference distribution
    βs = [range(0.0, 1.0, K)...]
    reverse!(βs)

    # make reference posterior for `λ`
    @views p = r[:,4]
    m     = mean(p)
    v     = var(p)
    λ_rdist = (m^2/v, m/v)

    # make reference posterior for `μ`
    @views p = r[:,5]
    m  = mean(p)
    sd = std(p)

    if sum(x -> x < 0.2, p) > sum(x -> 0.2 < x < 0.4, p)
      μ0 = 0.0
    else
      μ0 = m
    end

    σ0 = max(0.5, sd)

    x1 = run_newton(μ0, σ0, m, sd)

    μ_rdist = (x1[1], x1[2])

    # marginal likelihood
    pp = ref_posterior(Ξ, idf, λc, μc, v, λ_prior, μ_prior, λ_rdist, μ_rdist,
      nitpp, nthpp, βs, pup, sns, snodes!, scond, scond0)

    # process with reference distribution the posterior
    p1 = Vector{Float64}(undef, size(r,1))
    for i in Base.OneTo(size(r,1))
      p1[i] = r[i,2] + r[i,3] - 
              logdgamma(r[i,4], λ_rdist[1], λ_rdist[2]) -
              logdtnorm(r[i,5], μ_rdist[1], μ_rdist[2])
    end
    pp[1] = p1

    reverse!(pp)
    reverse!(βs)

    ml = gss(pp, βs)=#
  else
    ml = NaN
  end

  return r, treev, ml
end




"""
    mcmc_burn_cfbd(Ξ        ::Vector{sTfbd},
                   idf      ::Array{iBfffs,1},
                   λ_prior  ::NTuple{2,Float64},
                   μ_prior  ::NTuple{2,Float64},
                   λmμ_prior::NTuple{2,Float64},
                   ψ_prior  ::NTuple{2,Float64},
                   nburn    ::Int64,
                   λc       ::Float64,
                   μc       ::Float64,
                   ψc       ::Float64,
                   pup      ::Array{Int64,1}, 
                   prints   ::Int64,
                   sns      ::NTuple{3, BitVector},
                   snodes!  ::Function,
                   scond    ::Function,
                   scond0   ::Function)

Adaptive MCMC phase for da chain for constant fossilized birth-death using 
forward simulation.
"""
function mcmc_burn_cfbd(Ξ        ::Vector{sTfbd},
                        idf      ::Array{iBfffs,1},
                        λ_prior  ::NTuple{2,Float64},
                        μ_prior  ::NTuple{2,Float64},
                        λmμ_prior::NTuple{2,Float64},
                        ψ_prior  ::NTuple{2,Float64},
                        nburn    ::Int64,
                        λc       ::Float64,
                        μc       ::Float64,
                        ψc       ::Float64,
                        pup      ::Array{Int64,1}, 
                        prints   ::Int64,
                        sns      ::NTuple{3, BitVector},
                        snodes!  ::Function,
                        scond    ::Function,
                        scond0   ::Function)

  el = lastindex(idf)
  L  = treelength(Ξ)          # tree length
  nfos = Float64(nfossils(Ξ)) # number of fossilization events
  ns = (el-nfos-1)/2.0        # number of speciation events
  ne = 0.0                    # number of extinction events

  # add simulated subtrees to all fossil tips
  #=for bi in filter(x -> it(x) && ifos(x), idf)
    dri = dr(bi)
    ldr = lastindex(dri)
    t0 = sTfbd()
    ret = false
    while !ret
      t0, ret = fsψtip(bi, λc, μc, ψc)
    end
    swapfossil!(tree, t0, dri, ldr, 0)
  end=#

  # likelihood
  llc = llik_cfbd(Ξ, λc, μc, ψc) + scond(λc, μc, sns) + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2]) + 
        logdgamma(μc, μ_prior[1], μ_prior[2])

  if isnan(λmμ_prior[1])
    prc = logdgamma(λc, λ_prior[1], λ_prior[2]) + 
          logdgamma(μc, μ_prior[1], μ_prior[2]) + 
          logdgamma(ψc, ψ_prior[1], ψ_prior[2])
  else
    prc = logdgamma(λc,    λ_prior[1],   λ_prior[2]) + 
          logdgamma(λc-μc, λmμ_prior[1], λmμ_prior[2]) + 
          logdgamma(ψc,    ψ_prior[1],   ψ_prior[2])
  end

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      if p === 1
        if isnan(λmμ_prior[1])
          # λ proposal
          llc, prc, λc = update_λ!(llc, prc, λc, ns, L, μc, sns, 
                                   λ_prior, scond)
        else
          # parallel λ and μ proposal
          llc, prc, λc, μc = update_λμ!(llc, prc, λc, μc, ns, L, μc, ψc, sns, 
                                        λ_prior, scond)
        end

      elseif p === 2
        if isnan(λmμ_prior[1])
          # μ proposal
          llc, prc, μc = update_μ!(llc, prc, μc, ne, L, λc, sns, 
                                   μ_prior, scond)
        else
          # λ-μ proposal (μ proposal with constant λ)
          llc, prc, μc = update_λmμ!(llc, prc, μc, ne, L, λc, sns, 
                                     λmμ_prior, scond)
        end

      # ψ proposal
      elseif p === 3
        llc, prc, ψc = update_ψ!(llc, prc, ψc, nfos, L, ψ_prior)
      
      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)
        llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ns, ne, nfos, 
                                    L, sns, snodes!, scond0)
      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, ψc
end




"""
    mcmc_cfbd(Ξ      ::Vector{sTfbd},
              idf     ::Array{iBfffs,1},
              llc     ::Float64,
              prc     ::Float64,
              λc      ::Float64,
              μc      ::Float64,
              ψc      ::Float64,
              λ_prior ::NTuple{2,Float64},
              μ_prior ::NTuple{2,Float64},
              λmμ_prior::NTuple{2,Float64},
              ψ_prior  ::NTuple{2,Float64},
              niter   ::Int64,
              nthin   ::Int64,
              pup     ::Array{Int64,1}, 
              prints  ::Int64,
              sns     ::NTuple{3, BitVector},
              snodes! ::Function,
              scond   ::Function,
              scond0  ::Function)

MCMC da chain for constant fossilized birth-death using forward simulation.
"""
function mcmc_cfbd(Ξ      ::Vector{sTfbd},
                   idf     ::Array{iBfffs,1},
                   llc     ::Float64,
                   prc     ::Float64,
                   λc      ::Float64,
                   μc      ::Float64,
                   ψc      ::Float64,
                   λ_prior ::NTuple{2,Float64},
                   μ_prior ::NTuple{2,Float64},
                   λmμ_prior::NTuple{2,Float64},
                   ψ_prior  ::NTuple{2,Float64},
                   niter   ::Int64,
                   nthin   ::Int64,
                   pup     ::Array{Int64,1}, 
                   prints  ::Int64,
                   sns     ::NTuple{3, BitVector},
                   snodes! ::Function,
                   scond   ::Function,
                   scond0  ::Function)

  el = lastindex(idf)
  ns = Float64(nnodesbifurcation(Ξ))
  ne = Float64(ntipsextinct(Ξ))
  nfos = Float64(nfossils(Ξ))
  L  = treelength(Ξ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  R = Array{Float64,2}(undef, nlogs, 6)

  # make tree vector
  treev  = sTfbd[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      if p === 1
        if isnan(λmμ_prior[1])
          # λ proposal
          llc, prc, λc = update_λ!(llc, prc, λc, ns, L, μc, sns, 
                                   λ_prior, scond)
        else
          # parallel λ and μ proposal : TODO
          llc, prc, λc, μc = update_λμ!(llc, prc, λc, μc, ns, L, μc, ψc, sns, 
                                        λ_prior, scond)
        end

        llc, prc, λc = 
          update_λ!(llc, prc, λc, ns, L, μc, sns, λ_prior, scond)

        # llci = llik_cfbd(Ξ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      elseif p === 2
        if isnan(λmμ_prior[1])
          # μ proposal
          llc, prc, μc = update_μ!(llc, prc, μc, ne, L, λc, sns, 
                                   μ_prior, scond)
        else
          # λ-μ proposal (μ proposal with constant λ) : TODO
          llc, prc, μc = update_λmμ!(llc, prc, μc, ne, L, λc, sns, 
                                     λmμ_prior, scond)
        end

        # llci = llik_cfbd(Ξ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return 
        # end

      # ψ proposal
      elseif p === 3
        llc, prc, ψc = update_ψ!(llc, prc, ψc, nfos, L, ψ_prior)
      
      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)
        llc, ns, ne, nfos, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ψc, ns, ne, 
                                          nfos, L, sns, snodes!, scond0)

        # llci = llik_cfbd(Ξ, λc, μc) + scond(λc, μc, sns) + prob_ρ(idf)
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
        R[lit,6] = ψc
        push!(treev, couple(deepcopy(Ξ), idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, treev, λc, μc, ψc
end




#="""
    ref_posterior(Ξ        ::Vector{sTfbd},
                  idf      ::Array{iBfffs,1},
                  λc       ::Float64,
                  μc       ::Float64,
                  μtn      ::Float64,
                  ψtn      ::Float64,
                  λ_prior  ::NTuple{2,Float64},
                  μ_prior  ::NTuple{2,Float64},
                  λmμ_prior::NTuple{2,Float64},
                  ψ_prior  ::NTuple{2,Float64},
                  λ_rdist  ::NTuple{2,Float64},
                  μ_rdist  ::NTuple{2,Float64},
                  λmμ_rdist::NTuple{2,Float64},
                  ψ_rdist  ::NTuple{2,Float64},
                  nitpp    ::Int64,
                  nthpp    ::Int64,
                  βs       ::Vector{Float64},
                  pup      ::Array{Int64,1},
                  sns      ::NTuple{3, BitVector},
                  snodes!  ::Function,
                  scond    ::Function,
                  scond0   ::Function)

MCMC da chain for constant fossilized birth-death using forward simulation.
"""
function ref_posterior(Ξ        ::Vector{sTfbd},
                       idf      ::Array{iBfffs,1},
                       λc       ::Float64,
                       μc       ::Float64,
                       μtn      ::Float64,
                       ψtn      ::Float64,
                       λ_prior  ::NTuple{2,Float64},
                       μ_prior  ::NTuple{2,Float64},
                       λmμ_prior::NTuple{2,Float64},
                       ψ_prior  ::NTuple{2,Float64},
                       λ_rdist  ::NTuple{2,Float64},
                       μ_rdist  ::NTuple{2,Float64},
                       λmμ_rdist::NTuple{2,Float64},
                       ψ_rdist  ::NTuple{2,Float64},
                       nitpp    ::Int64,
                       nthpp    ::Int64,
                       βs       ::Vector{Float64},
                       pup      ::Array{Int64,1},
                       sns      ::NTuple{3, BitVector},
                       snodes!  ::Function,
                       scond    ::Function,
                       scond0   ::Function)

  K = lastindex(βs)

  # make log-likelihood table per power
  nlg = fld(nitpp, nthpp)
  pp  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

  el = lastindex(idf)
  ns = Float64(nnodesinternal(Ξ))
  ne = Float64(ntipsextinct(Ξ))
  L  = treelength(Ξ)

  llc = llik_cfbd(Ξ, λc, μc, ψc) + scond(λc, μc, sns) + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2]) + 
        logdgamma(μc, μ_prior[1], μ_prior[2]) + 
        logdgamma(ψc, ψ_prior[1], ψ_prior[2])

  for k in 2:K

    βi  = βs[k]
    rdc = logdgamma(λc, λ_rdist[1], λ_rdist[2]) + 
          logdtnorm(μc, μ_rdist[1], μ_rdist[2]) + 
          logdgamma(ψc, ψ_rdist[1], ψ_rdist[2])

    # logging
    lth, lit = 0, 0

    for it in Base.OneTo(nitpp)

      shuffle!(pup)

      for p in pup

        # λ proposal
        if p === 1

          llc, prc, rdc, λc = 
            update_λ!(llc, prc, rdc, λc, ns, L, μc, sns, λ_prior, λ_rdist, 
              scond, βi)

        # forward simulation proposal proposal
        elseif p === 2 

          llc, prc, rdc, μc = 
            update_μ!(llc, prc, rdc, μc, ne, L, μtn, λc, sns, μ_prior, μ_rdist, 
              scond, βi)

        else

          bix = ceil(Int64,rand()*el)
          llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L, sns,
                             snodes!, scond0)

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
end=#




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{sTfbd},
               idf  ::Vector{iBfffs},
               llc  ::Float64,
               λ    ::Float64, 
               μ    ::Float64,
               ψ    ::Float64,
               ns   ::Float64,
               ne   ::Float64,
               nfos ::Float64,
               L    ::Float64,
               scond::Function,
               pow  ::Float64)

Forward simulation proposal function for constant fossilized birth-death.
"""
function update_fs!(bix    ::Int64,
                    Ξ      ::Vector{sTfbd},
                    idf    ::Vector{iBfffs},
                    llc    ::Float64,
                    λ      ::Float64, 
                    μ      ::Float64,
                    ψ      ::Float64,
                    ns     ::Float64,
                    ne     ::Float64,
                    nfos   ::Float64,
                    L      ::Float64,
                    sns    ::NTuple{3, BitVector},
                    snodes!::Function,
                    scond0 ::Function)

  bi = idf[bix]

  # forward simulate an internal branch
  ξp, np, ntp = fsbi(bi, λ, μ, ψ, 100)
  
  # retained conditional on survival
  if ntp > 0

    itb = it(bi)   # is it terminal
    iψb = ifos(bi) # is it a fossil
    ρbi = ρi(bi)   # get branch sampling fraction
    nc  = ni(bi)   # current ni
    ntc = nt(bi)   # current nt

    # current tree
    ξc  = Ξ[bix]

    # if terminal branch
    if itb
      llr = log(Float64(np)/Float64(nc) * (1.0 - ρbi)^(np - nc))
      acr = 0.0
    else
      np  -= 1
      llr = log((1.0 - ρbi)^(np - nc))
      acr = log(Float64(ntp)/Float64(ntc))
    end

    # MH ratio
    if -randexp() < llr + acr

      # if survival conditioned
      scn = (iszero(pa(bi)) && e(bi)>0.0) || (isone(pa(bi)) && iszero(e(Ξ[1])))
      if scn
          llr += scond0(ξp, λ, μ, itb) - scond0(ξc, λ, μ, itb)
      end

      # update ns, ne, nfos & L
      ns +=   Float64(nnodesbifurcation(ξp) - nnodesbifurcation(ξc))
      ne +=   Float64(ntipsextinct(ξp)      - ntipsextinct(ξc))
      nfos += Float64(nfossils(ξp)          - nfossils(ξc))
      L  +=   treelength(ξp)                - treelength(ξc)

      # likelihood ratio
      llr += llik_cfbd(ξp, λ, μ, ψ) - llik_cfbd(ξc, λ, μ, ψ)

      Ξ[bix] = ξp     # set new decoupled tree
      llc += llr      # set new likelihood
      
      if scn
        snodes!(Ξ, sns) # set new sns
      end
      
      setni!(bi, np)  # set new ni
      setnt!(bi, ntp) # set new nt
    end
  end

  return llc, ns, ne, nfos, L
end





"""
    fsbi(bi::iBfffs, λ::Float64, μ::Float64, ψ::Float64, ntry::Int64)

Forward simulation for branch `bi`
"""
function fsbi(bi::iBfffs, λ::Float64, μ::Float64, ψ::Float64, ntry::Int64)

  # times
  tfb = tf(bi)

  ext = 0
  # condition on non-extinction (helps in mixing)
  while ext < ntry 

    # forward simulation during branch length
    t0, na, nfos = sim_cfbd(e(bi), λ, μ, ψ, 0, 0)

    if iszero(nfos) # Exclude if any fossil is sampled
      nat = na
      
      if isone(na)
        fixalive!(t0)
        if ifos(bi) maketipsfossil!(t0) end
        return t0, na, nat
      
      elseif na > 1
        # fix random tip
        fixrtip!(t0)

        if !it(bi)
          # add tips until the present
          tx, na, nfos = tip_sims!(t0, tfb, λ, μ, ψ, na, nfos)
          if !iszero(nfos)
            ext += 1
            continue
          end
        end

        if ifos(bi) maketipsfossil!(t0) end

        return t0, na, nat
      end
    end

    ext += 1
  end

  return sTfbd(), 0, 0
end




"""
    tip_sims!(tree::sTfbd, t::Float64, λ::Float64, μ::Float64, 
              ψ::Float64, na::Int64, nfos::Int64)

Continue simulation until time `t` for unfixed tips in `tree`. 
"""
function tip_sims!(tree::sTfbd, t::Float64, λ::Float64, μ::Float64, 
                   ψ::Float64, na::Int64, nfos::Int64)
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if !defd1 && !defd2
    # tips
    if !isfix(tree) && isalive(tree)

      # simulate
      stree, na, nfos = sim_cfbd(t, λ, μ, ψ, na-1, 0)

      if iszero(nfos)
        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    end
  else
    # bifurcations and sampled fossil ancestors
    if defd1 tree.d1, na, nfos = tip_sims!(tree.d1, t, λ, μ, ψ, na, nfos) end
    if !iszero(nfos) return tree, na, nfos end
    if defd2 tree.d2, na, nfos = tip_sims!(tree.d2, t, λ, μ, ψ, na, nfos) end
  end

  return tree, na, nfos
end




"""
    update_λmμ!(llc      ::Float64,
                prc      ::Float64,
                λc       ::Float64,
                μc       ::Float64,
                ns       ::Float64,
                ne       ::Float64,
                L        ::Float64,
                μc       ::Float64,
                sns      ::NTuple{3,BitVector},
                λmμ_prior::NTuple{2,Float64},
                scond    ::Function)

Gibbs sampling of `λ-μ` for constant fossilized birth-death.
"""
#=function update_λmμ!(llc      ::Float64,
                     prc      ::Float64,
                     λc       ::Float64,
                     μc       ::Float64,
                     ns       ::Float64,
                     ne       ::Float64,
                     L        ::Float64,
                     μc       ::Float64,
                     sns      ::NTuple{3,BitVector},
                     λmμ_prior::NTuple{2,Float64},
                     scond    ::Function)

  λmμc = λc-μc
  λmμp  = randgamma(λmμ_prior[1] + ns - ne, λmμ_prior[2] + L)
  μp = λc-λmμp
  llr = scond(λc, μp, sns) - scond(λc, μc, sns)

  if -randexp() < llr
    llc += (ns - ne) * log(λmμp/λmμc) + L * (λmμc - λmμp) + llr
    prc += llrdgamma(λmμp, λmμc, λmμ_prior[1], λmμ_prior[2])
    μc   = μp
  end

  return llc, prc, λc, μc
end

=> Probably doesn't work + change prior distribution

=#




"""
    update_ψ!(llc    ::Float64,
              prc    ::Float64,
              ψc     ::Float64,
              nfos   ::Float64,
              L      ::Float64,
              ψ_prior::NTuple{2,Float64})

Gibbs sampling of `ψ` for constant fossilized birth-death.
"""
function update_ψ!(llc    ::Float64,
                   prc    ::Float64,
                   ψc     ::Float64,
                   nfos   ::Float64,
                   L      ::Float64,
                   ψ_prior::NTuple{2,Float64})

  ψp  = randgamma(ψ_prior[1] + nfos, ψ_prior[2] + L)

  llc += nfos * log(ψp/ψc) + L * (ψc - ψp)
  prc += llrdgamma(ψp, ψc, ψ_prior[1], ψ_prior[2])
  ψc   = ψp

  return llc, prc, ψc 
end


