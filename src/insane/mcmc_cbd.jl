#=

constant birth-death MCMC using forward simulation

Ignacio Quintero Mächler

t(-_-t)

Created 25 08 2020
=#




"""
    insane_cbd(tree    ::sT_label;
               λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
               niter   ::Int64                 = 1_000,
               nthin   ::Int64                 = 10,
               nburn   ::Int64                 = 200,
               nflush  ::Int64                 = nthin,
               ofile   ::String                = string(homedir(), "/cbd"),
               ϵi      ::Float64               = 0.4,
               λi      ::Float64               = NaN,
               μi      ::Float64               = NaN,
               pupdp   ::NTuple{3,Float64}     = (0.2,0.2,0.2),
               prints  ::Int64                 = 5,
               survival::Bool                  = true,
               mxthf   ::Float64               = Inf,
               tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant birth-death.
"""
function insane_cbd(tree    ::sT_label;
                    λ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    μ_prior ::NTuple{2,Float64}     = (1.0, 1.0),
                    niter   ::Int64                 = 1_000,
                    nthin   ::Int64                 = 10,
                    nburn   ::Int64                 = 200,
                    nflush  ::Int64                 = nthin,
                    ofile   ::String                = string(homedir(), "/cbd"),
                    ϵi      ::Float64               = 0.4,
                    λi      ::Float64               = NaN,
                    μi      ::Float64               = NaN,
                    pupdp   ::NTuple{3,Float64}     = (0.2,0.2,0.2),
                    prints  ::Int64                 = 5,
                    survival::Bool                  = true,
                    mxthf   ::Float64               = Inf,
                    tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n     = ntips(tree)
  th    = treeheight(tree)

  surv = 0   # condition on survival of 0, 1, or 2 starting lineages
  rmλ  = 0.0 # condition on first speciation event
  if survival 
    if iszero(e(tree)) 
      surv += 2
      rmλ  += 1.0
    else
      surv += 1
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
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end
  # M attempts of survival
  mc = m_surv_cbd(th, λc, μc, 5_000, surv)

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, sTbd)

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(3)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant birth-death"

  # adaptive phase
  llc, prc, λc, μc, mc, ns, L =
      mcmc_burn_cbd(Ξ, idf, λ_prior, μ_prior, nburn, λc, μc, mc, th, rmλ, surv,
        pup, prints)

  # mcmc
  r, treev, λc, μc, mc = mcmc_cbd(Ξ, idf, llc, prc, λc, μc, mc, ns, L, 
    th, rmλ, surv, λ_prior, μ_prior, pup, niter, nthin, nflush, ofile, prints)

  # if marginal

  #    # reference distribution
  #   βs = [range(0.0, 1.0, K)...]
  #   reverse!(βs)

  #   # make reference posterior for `λ`
  #   @views p = r[:,4]
  #   m     = mean(p)
  #   v     = var(p)
  #   λ_rdist = (m^2/v, m/v)

  #   # make reference posterior for `μ`
  #   @views p = r[:,5]
  #   m  = mean(p)
  #   sd = std(p)

  #   if sum(x -> x < 0.2, p) > sum(x -> 0.2 < x < 0.4, p)
  #     μ0 = 0.0
  #   else
  #     μ0 = m
  #   end

  #   σ0 = max(0.5, sd)

  #   x1 = run_newton(μ0, σ0, m, sd)

  #   μ_rdist = (x1[1], x1[2])

  #   # marginal likelihood
  #   pp = ref_posterior(Ξ, idf, λc, μc, v, mc, th, crown,
  #     λ_prior, μ_prior, λ_rdist, μ_rdist, nitpp, nthpp, βs, pup)

  #   # process with reference distribution the posterior
  #   p1 = Vector{Float64}(undef, size(r,1))
  #   for i in Base.OneTo(size(r,1))
  #     p1[i] = r[i,2] + r[i,3] -
  #             logdgamma(r[i,4], λ_rdist[1], λ_rdist[2]) -
  #             logdtnorm(r[i,5], μ_rdist[1], μ_rdist[2])
  #   end
  #   pp[1] = p1

  #   reverse!(pp)
  #   reverse!(βs)

  #   ml = gss(pp, βs)
  # else
  #   ml = NaN
  # end

  return r, treev
end




"""
    mcmc_burn_cbd(Ξ      ::Vector{sTbd},
                  idf    ::Array{iBffs,1},
                  λ_prior::NTuple{2,Float64},
                  μ_prior::NTuple{2,Float64},
                  nburn  ::Int64,
                  λc     ::Float64,
                  μc     ::Float64,
                  mc     ::Float64,
                  th     ::Float64,
                  rmλ    ::Float64,
                  surv   ::Int64,
                  pup    ::Array{Int64,1},
                  prints ::Int64)

Adaptive MCMC phase for da chain for constant birth-death using forward
simulation.
"""
function mcmc_burn_cbd(Ξ      ::Vector{sTbd},
                       idf    ::Array{iBffs,1},
                       λ_prior::NTuple{2,Float64},
                       μ_prior::NTuple{2,Float64},
                       nburn  ::Int64,
                       λc     ::Float64,
                       μc     ::Float64,
                       mc     ::Float64,
                       th     ::Float64,
                       rmλ    ::Float64,
                       surv   ::Int64,
                       pup    ::Array{Int64,1},
                       prints ::Int64)

  el  = lastindex(idf)
  L   = treelength(Ξ)          # tree length
  ns  = nnodesbifurcation(idf) # number of speciation events
  ne  = 0.0                    # number of extinction events

  # likelihood
  llc = llik_cbd(Ξ, λc, μc, ns) - rmλ * log(λc) + 
        log(mc) + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2]) +
        logdgamma(μc, μ_prior[1], μ_prior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc, mc =
          update_λ!(llc, prc, λc, ns, L, μc, mc, th, rmλ, surv, λ_prior)

      # μ proposal
      elseif p === 2

        llc, prc, μc, mc =
          update_μ!(llc, prc, μc, ne, L, λc, mc, th, surv, μ_prior)

      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)

        llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L)
      end
    end

    next!(pbar)
  end

  return llc, prc, λc, μc, mc, ns, L
end




"""
    mcmc_cbd(Ξ      ::Vector{sTbd},
             idf    ::Array{iBffs,1},
             llc    ::Float64,
             prc    ::Float64,
             λc     ::Float64,
             μc     ::Float64,
             mc     ::Float64,
             ns     ::Float64,
             L      ::Float64,
             th     ::Float64,
             rmλ    ::Float64,
             surv   ::Int64,
             λ_prior::NTuple{2,Float64},
             μ_prior::NTuple{2,Float64},
             pup    ::Array{Int64,1},
             niter  ::Int64,
             nthin  ::Int64,
             nflush ::Int64,
             ofile  ::String,
             prints ::Int64)

MCMC da chain for constant birth-death using forward simulation.
"""
function mcmc_cbd(Ξ      ::Vector{sTbd},
                  idf    ::Array{iBffs,1},
                  llc    ::Float64,
                  prc    ::Float64,
                  λc     ::Float64,
                  μc     ::Float64,
                  mc     ::Float64,
                  ns     ::Float64,
                  L      ::Float64,
                  th     ::Float64,
                  rmλ    ::Float64,
                  surv   ::Int64,
                  λ_prior::NTuple{2,Float64},
                  μ_prior::NTuple{2,Float64},
                  pup    ::Array{Int64,1},
                  niter  ::Int64,
                  nthin  ::Int64,
                  nflush ::Int64,
                  ofile  ::String,
                  prints ::Int64)

  el = lastindex(idf)
  ne = Float64(ntipsextinct(Ξ))

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  # parameter results
  r = Array{Float64,2}(undef, nlogs, 5)

  # make tree vector
  treev  = sTbd[]

  # flush to file
  sthin = 0

  open(ofile*".log", "w") do of
    write(of, "iteration\tlikelihood\tprior\tlambda\tmu\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      pbar = Progress(niter, prints, "running mcmc...", 20)

      for it in Base.OneTo(niter)

        shuffle!(pup)

        for p in pup

          # λ proposal
          if p === 1

            llc, prc, λc, mc =
              update_λ!(llc, prc, λc, ns, L, μc, mc, th, rmλ, surv, λ_prior)

            # llci = llik_cbd(Ξ, λc, μc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
            # if !isapprox(llci, llc, atol = 1e-6)
            #    @show llci, llc, it, p
            #    return
            # end

          # μ proposal
          elseif p === 2

            llc, prc, μc, mc =
              update_μ!(llc, prc, μc, ne, L, λc, mc, th, surv, μ_prior)

            # llci = llik_cbd(Ξ, λc, μc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
            # if !isapprox(llci, llc, atol = 1e-6)
            #    @show llci, llc, it, p
            #    return
            # end

          # forward simulation proposal proposal
          else

            bix = ceil(Int64,rand()*el)
            llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L)

            # llci = llik_cbd(Ξ, λc, μc, nnodesbifurcation(idf)) - rmλ * log(λc) + log(mc) + prob_ρ(idf)
            # if !isapprox(llci, llc, atol = 1e-6)
            #    @show llci, llc, it, p
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
            r[lit,4] = λc
            r[lit,5] = μc
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthin += 1
        if sthin === nflush
          write(of, 
            string(Float64(it), "\t", llc, "\t", prc, "\t", λc,"\t", μc, "\n"))
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

  return r, treev, λc, μc, mc
end




# """
#     ref_posterior(Ξ      ::Vector{sTbd},
#                   idf    ::Array{iBffs,1},
#                   λc     ::Float64,
#                   μc     ::Float64,
#                   μtn    ::Float64,
#                   mc     ::Float64,
#                   th     ::Float64,
#                   crown   ::Bool,
#                   λ_prior::NTuple{2,Float64},
#                   μ_prior::NTuple{2,Float64},
#                   λ_rdist::NTuple{2,Float64},
#                   μ_rdist::NTuple{2,Float64},
#                   nitpp  ::Int64,
#                   nthpp  ::Int64,
#                   βs     ::Vector{Float64},
#                   pup    ::Array{Int64,1})

# MCMC da chain for constant birth-death using forward simulation.
# """
# function ref_posterior(Ξ      ::Vector{sTbd},
#                        idf    ::Array{iBffs,1},
#                        λc     ::Float64,
#                        μc     ::Float64,
#                        μtn    ::Float64,
#                        mc     ::Float64,
#                        th     ::Float64,
#                        crown   ::Bool,
#                        λ_prior::NTuple{2,Float64},
#                        μ_prior::NTuple{2,Float64},
#                        λ_rdist::NTuple{2,Float64},
#                        μ_rdist::NTuple{2,Float64},
#                        nitpp  ::Int64,
#                        nthpp  ::Int64,
#                        βs     ::Vector{Float64},
#                        pup    ::Array{Int64,1})

#   K = lastindex(βs)

#   # make log-likelihood table per power
#   nlg = fld(nitpp, nthpp)
#   pp  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

#   el = lastindex(idf)
#   ns = Float64(nnodesinternal(Ξ))
#   ne = Float64(ntipsextinct(Ξ))
#   L  = treelength(Ξ)

#   nsi = crown ? 0.0 : log(λc)

#   llc = llik_cbd(Ξ, λc, μc, ns) - nsi + log(mc) + prob_ρ(idf)
#   prc = logdgamma(λc, λ_prior[1], λ_prior[2]) +
#         logdgamma(μc, μ_prior[1], μ_prior[2])

#   for k in 2:K

#     βi  = βs[k]
#     rdc = logdgamma(λc, λ_rdist[1], λ_rdist[2]) +
#           logdtnorm(μc, μ_rdist[1], μ_rdist[2])

#     # logging
#     lth, lit = 0, 0

#     for it in Base.OneTo(nitpp)

#       shuffle!(pup)

#       for p in pup

#         # λ proposal
#         if p === 1

#           llc, prc, rdc, λc, mc =
#             update_λ!(llc, prc, rdc, λc, ns, L, μc, mc, th, crown,
#               λ_prior, λ_rdist, βi)

#         # forward simulation proposal proposal
#         elseif p === 2

#           llc, prc, rdc, μc, mc =
#             update_μ!(llc, prc, rdc, μc, ne, L, μtn, λc, mc, th, crown,
#               μ_prior, μ_rdist, βi)

#         else

#           bix = ceil(Int64,rand()*el)
#           llc, ns, ne, L = update_fs!(bix, Ξ, idf, llc, λc, μc, ns, ne, L)

#         end
#       end

#       # log log-likelihood
#       lth += 1
#       if lth === nthpp
#         lit += 1
#         pp[k][lit] = llc + prc - rdc
#         lth = 0
#       end
#     end

#     @info string(βi," power done")
#   end

#   return pp
# end




"""
    update_fs!(bix::Int64,
               Ξ  ::Vector{sTbd},
               idf::Vector{iBffs},
               llc::Float64,
               λ  ::Float64,
               μ  ::Float64,
               ns ::Float64,
               ne ::Float64,
               L  ::Float64)

Forward simulation proposal function for constant birth-death.
"""
function update_fs!(bix::Int64,
                    Ξ  ::Vector{sTbd},
                    idf::Vector{iBffs},
                    llc::Float64,
                    λ  ::Float64,
                    μ  ::Float64,
                    ns ::Float64,
                    ne ::Float64,
                    L  ::Float64)

  bi = idf[bix]

   # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, λ, μ)

  # if mid or internal
  else
    ξp, llr = fsbi_i(bi, λ, μ)
  end

  if isfinite(llr)
    ξc  = Ξ[bix]

    # update llc, ns, ne & L
    llc += llik_cbd(ξp, λ, μ) - llik_cbd(ξc, λ, μ) + llr
    ns  += Float64(nnodesinternal(ξp) - nnodesinternal(ξc))
    ne  += Float64(ntipsextinct(ξp)   - ntipsextinct(ξc))
    L   += treelength(ξp)             - treelength(ξc)

    # set new decoupled tree
    Ξ[bix] = ξp
  end

  return llc, ns, ne, L
end




"""
    fsbi(bi::iBffs, λ::Float64, μ::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs, λ::Float64, μ::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, na, nn, llr =
    _sim_cbd_t(e(bi), λ, μ, lc, lU, Iρi, 0, 1, 1_000)

  if na > 0 && isfinite(llr)
    _fixrtip!(t0, na) # fix random tip
    setni!(bi, na)    # set new ni

    return t0, llr
  else
    return t0, -Inf
  end
end




"""
    fsbi_i(bi::iBffs, λ::Float64, μ::Float64)

Forward simulation for internal branch.
"""
function fsbi_i(bi::iBffs, λ::Float64, μ::Float64)

  t0, na, nn = _sim_cbd_i(e(bi), λ, μ, 0, 1, 1_000)

  if na < 1 || nn > 999
    return t0, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

  if lU < acr

    _fixrtip!(t0, na) # fix random tip

    # simulate remaining tips until the present
    if na > 1
      tx, na, nn, acr = tip_sims!(t0, tf(bi), λ, μ, acr, lU, Iρi, na, nn)
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
    tip_sims!(tree::sTbd,
              t   ::Float64,
              λ   ::Float64,
              μ   ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTbd,
                   t   ::Float64,
                   λ   ::Float64,
                   μ   ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64,
                   nn ::Int64)

  if lU < lr && nn < 1_000

    if istip(tree)
      if !isfix(tree) && isalive(tree)

        # simulate
        stree, na, nn, lr = 
          _sim_cbd_it(t, λ, μ, lr, lU, Iρi, na-1, nn, 1_000)

        if isnan(lr) || nn > 999
          return tree, na, nn, NaN
        end

        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        setproperty!(tree, :iμ, isextinct(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, nn, lr = tip_sims!(tree.d1, t, λ, μ, lr, lU, Iρi, na, nn)
      tree.d2, na, nn, lr = tip_sims!(tree.d2, t, λ, μ, lr, lU, Iρi, na, nn)
    end

    return tree, na, nn, lr
  end

  return tree, na, nn, NaN
end




"""
    update_λ!(llc    ::Float64,
              prc    ::Float64,
              λc     ::Float64,
              ns     ::Float64,
              L      ::Float64,
              μc     ::Float64,
              mc     ::Float64,
              th     ::Float64,
              rmλ    ::Float64,
              surv   ::Int64,
              λ_prior::NTuple{2,Float64})

Mixed HM-Gibbs sampling of `λ` for constant birth-death.
"""
function update_λ!(llc    ::Float64,
                   prc    ::Float64,
                   λc     ::Float64,
                   ns     ::Float64,
                   L      ::Float64,
                   μc     ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   rmλ    ::Float64,
                   surv   ::Int64,
                   λ_prior::NTuple{2,Float64})

  λp  = randgamma(λ_prior[1] + ns - rmλ, λ_prior[2] + L)

  mp  = m_surv_cbd(th, λp, μc, 5_000, surv)
  llr = log(mp/mc)

  if -randexp() < llr
    llc += (ns - rmλ) * log(λp/λc) + L * (λc - λp) + llr
    prc += llrdgamma(λp, λc, λ_prior[1], λ_prior[2])
    λc   = λp
    mc   = mp
  end

  return llc, prc, λc, mc
end




# """
#     update_λ!(llc    ::Float64,
#               prc    ::Float64,
#               rdc    ::Float64,
#               λc     ::Float64,
#               ns     ::Float64,
#               L      ::Float64,
#               μc     ::Float64,
#               mc     ::Float64,
#               th     ::Float64,
#               crown   ::Bool,
#               λ_prior::NTuple{2,Float64},
#               λ_rdist::NTuple{2,Float64},
#               pow    ::Float64)

# Mixed HM-Gibbs of `λ` for constant birth-death with reference distribution.
# """
# function update_λ!(llc    ::Float64,
#                    prc    ::Float64,
#                    rdc    ::Float64,
#                    λc     ::Float64,
#                    ns     ::Float64,
#                    L      ::Float64,
#                    μc     ::Float64,
#                    mc     ::Float64,
#                    th     ::Float64,
#                    crown  ::Int64,
#                    λ_prior::NTuple{2,Float64},
#                    λ_rdist::NTuple{2,Float64},
#                    pow    ::Float64)

#   λp  = randgamma((λ_prior[1] + ns - rmλ) * pow + λ_rdist[1] * (1.0 - pow),
#                   (λ_prior[2] + L) * pow         + λ_rdist[2] * (1.0 - pow))
#   mp  = m_surv_cbd(th, λp, μc, 5_000, crown)
#   llr = log(mp/mc)

#   if -randexp() < (pow * llr)
#     llc += (ns - rmλ) * log(λp/λc) + L * (λc - λp) + llr
#     prc += llrdgamma(λp, λc, λ_prior[1], λ_prior[2])
#     rdc += llrdgamma(λp, λc, λ_rdist[1], λ_rdist[2])
#     λc   = λp
#     mc   = mp
#   end

#   return llc, prc, rdc, λc, mc
# end




"""
    update_μ!(llc    ::Float64,
              prc    ::Float64,
              μc     ::Float64,
              ne     ::Float64,
              L      ::Float64,
              λc     ::Float64,
              mc     ::Float64,
              th     ::Float64,
              crown  ::Int64,
              μ_prior::NTuple{2,Float64})

Mixed HM-Gibbs of `μ` for constant birth-death.
"""
function update_μ!(llc    ::Float64,
                   prc    ::Float64,
                   μc     ::Float64,
                   ne     ::Float64,
                   L      ::Float64,
                   λc     ::Float64,
                   mc     ::Float64,
                   th     ::Float64,
                   surv  ::Int64,
                   μ_prior::NTuple{2,Float64})

  μp  = randgamma(μ_prior[1] + ne, μ_prior[2] + L)

  mp   = m_surv_cbd(th, λc, μp, 5_000, surv)
  llr  = log(mp/mc)

  if -randexp() < llr
    llc += ne * log(μp/μc) + L * (μc - μp) + llr
    prc += llrdgamma(μp, μc, μ_prior[1], μ_prior[2])
    μc   = μp
    mc   = mp
  end

  return llc, prc, μc, mc
end




# """
#     update_μ!(llc    ::Float64,
#               prc    ::Float64,
#               rdc    ::Float64,
#               μc     ::Float64,
#               ne     ::Float64,
#               L      ::Float64,
#               μtn    ::Float64,
#               λc     ::Float64,
#               mc     ::Float64,
#               th     ::Float64,
#               crown  ::Int64,
#               μ_prior::NTuple{2,Float64},
#               μ_rdist::NTuple{2,Float64},
#               pow    ::Float64)

# Mixed HM-Gibbs of `μ` for constant birth-death with reference distribution.
# """
# function update_μ!(llc    ::Float64,
#                    prc    ::Float64,
#                    rdc    ::Float64,
#                    μc     ::Float64,
#                    ne     ::Float64,
#                    L      ::Float64,
#                    μtn    ::Float64,
#                    λc     ::Float64,
#                    mc     ::Float64,
#                    th     ::Float64,
#                    crown  ::Int64,
#                    μ_prior::NTuple{2,Float64},
#                    μ_rdist::NTuple{2,Float64},
#                    pow    ::Float64)

#   μp  = mulupt(μc, μtn)::Float64
#   mp  = m_surv_cbd(th, λc, μp, 5_000, crown)

#   μr  = log(μp/μc)
#   llr = ne * μr + L * (μc - μp) + log(mp/mc)
#   prr = llrdgamma(μp, μc, μ_prior[1], μ_prior[2])
#   rdr = llrdtnorm(μp, μc, μ_rdist[1], μ_rdist[2])

#   if -randexp() < (pow * (llr + prr) + (1.0 - pow) * rdr + μr)
#     llc += llr
#     prc += prr
#     rdc += rdr
#     μc   = μp
#     mc   = mp
#   end

#   return llc, prc, rdc, μc, mc
# end


