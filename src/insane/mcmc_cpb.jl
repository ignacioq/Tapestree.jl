#=

pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    insane_cpb(tree    ::sT_label,
               out_file::String;
               λ_prior  ::NTuple{2,Float64}     = (1.0, 1.0),
               niter   ::Int64                 = 1_000,
               nthin   ::Int64                 = 10,
               nburn   ::Int64                 = 200,
               tune_int::Int64                 = 100,
               marginal    ::Bool                  = false,
               nitpp   ::Int64                 = 100,
               nthpp   ::Int64                 = 10,
               K       ::Int64                 = 10,
               λi      ::Float64               = NaN,
               pupdp   ::NTuple{2,Float64}     = (0.2, 0.2),
               prints  ::Int64                 = 5,
               tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for constant pure-birth.
"""
function insane_cpb(tree    ::sT_label,
                    out_file::String;
                    λ_prior  ::NTuple{2,Float64}    = (1.0, 1.0),
                    niter   ::Int64                 = 1_000,
                    nthin   ::Int64                 = 10,
                    nburn   ::Int64                 = 200,
                    tune_int::Int64                 = 100,
                    # marginal ::Bool                 = false,
                    # nitpp   ::Int64                 = 100,
                    # nthpp   ::Int64                 = 10,
                    # K       ::Int64                 = 10,
                    λi      ::Float64               = NaN,
                    pupdp   ::NTuple{2,Float64}     = (0.2, 0.2),
                    prints  ::Int64                 = 5,
                    tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  stem = !iszero(e(tree))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ, Inf)

  # starting parameters
  if isnan(λi)
    λc = Float64(n-2)/treelength(tree)
  else
    λc = λi
  end

  # make a decoupled tree and fix it
  Ξ = make_Ξ(idf, sTpb)

  # make parameter updates scaling function for tuning
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(2)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "Running constant pure-birth"

  # adaptive phase
  llc, prc, λc =
    mcmc_burn_cpb(Ξ, idf, λ_prior, nburn, λc, pup, prints, stem)

  # mcmc
  r, treev, λc =
    mcmc_cpb(Ξ, idf, llc, prc, λc, λ_prior, niter, nthin, pup, prints, stem)

  pardic = Dict(("lambda" => 1))

  write_ssr(r, pardic, out_file)

  # if marginal
  #   # reference distribution
  #   βs = [range(0.0, 1.0, K)...]::Vector{Float64}
  #   reverse!(βs)

  #   @views p = r[:,4]

  #   # make reference posterior
  #   m     = mean(p)
  #   v     = var(p)
  #   λ_refd = (m^2/v, m/v)

  #   # marginal likelihood
  #   pp = ref_posterior(Ξ, idf, λc, λ_prior, λ_refd, nitpp, nthpp, βs, pup, stem)

  #   # process with reference distribution the posterior
  #   p1 = Vector{Float64}(undef, size(r,1))
  #   for i in Base.OneTo(size(r,1))
  #     p1[i] = r[i,2] + r[i,3] - logdgamma(r[i,4], λ_refd[1], λ_refd[2])
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
    mcmc_burn_cpb(Ξ      ::Vector{sTpb},
                  idf    ::Array{iBffs,1},
                  λ_prior::NTuple{2,Float64},
                  nburn  ::Int64,
                  λc     ::Float64,
                  pup    ::Array{Int64,1},
                  prints ::Int64,
                  stem   ::Bool)

MCMC chain for constant pure-birth.
"""
function mcmc_burn_cpb(Ξ      ::Vector{sTpb},
                       idf    ::Array{iBffs,1},
                       λ_prior::NTuple{2,Float64},
                       nburn  ::Int64,
                       λc     ::Float64,
                       pup    ::Array{Int64,1},
                       prints ::Int64,
                       stem   ::Bool)

  el  = lastindex(idf)
  L   = treelength(Ξ)     # tree length
  ns  = Float64(el-1)*0.5 # number of speciation events
  nsi = stem ? 0.0 : log(λc)

  #likelihood
  llc = llik_cpb(Ξ, λc) - nsi + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2])

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc = update_λ!(llc, prc, λc, ns, L, stem, λ_prior)

      # forward simulation proposal proposal
      else
        bix = ceil(Int64,rand()*el)

        llc, ns, L = update_fs!(bix, Ξ, idf, llc, λc, ns, L)
      end
    end

    next!(pbar)
  end

  return llc, prc, λc
end




"""
    mcmc_cpb(Ξ      ::Vector{sTpb},
             idf    ::Array{iBffs,1},
             llc    ::Float64,
             prc    ::Float64,
             λc     ::Float64,
             λ_prior::NTuple{2,Float64},
             niter  ::Int64,
             nthin  ::Int64,
             pup    ::Array{Int64,1},
             prints ::Int64,
             stem   ::Bool)

MCMC chain for constant pure-birth.
"""
function mcmc_cpb(Ξ      ::Vector{sTpb},
                  idf    ::Array{iBffs,1},
                  llc    ::Float64,
                  prc    ::Float64,
                  λc     ::Float64,
                  λ_prior::NTuple{2,Float64},
                  niter  ::Int64,
                  nthin  ::Int64,
                  pup    ::Array{Int64,1},
                  prints ::Int64,
                  stem   ::Bool)

  el = lastindex(idf)
  ns = Float64(nnodesinternal(Ξ))
  L  = treelength(Ξ)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  r = Array{Float64,2}(undef, nlogs, 4)

  # make tree vector
  treev  = sTpb[]

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for p in pup

      # λ proposal
      if p === 1

        llc, prc, λc = update_λ!(llc, prc, λc, ns, L, stem, λ_prior)

        # llci = llik_cpb(Ξ, λc) - !stem*log(λc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end

      # forward simulation proposal proposal
      else

        bix = ceil(Int64,rand()*el)
        llc, ns, L = update_fs!(bix, Ξ, idf, llc, λc, ns, L)

        # llci = llik_cpb(Ξ, λc) - !stem*log(λc) + prob_ρ(idf)
        # if !isapprox(llci, llc, atol = 1e-6)
        #    @show llci, llc, it, p
        #    return
        # end
      end
    end

    lthin += 1
    if lthin == nthin
      lit += 1
      @inbounds begin
        r[lit,1] = Float64(lit)
        r[lit,2] = llc
        r[lit,3] = prc
        r[lit,4] = λc
        push!(treev, couple(Ξ, idf, 1))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return r, treev, λc
end




"""
    ref_posterior(Ξ      ::Vector{sTpb},
                  idf    ::Array{iBffs,1},
                  λc     ::Float64,
                  λ_prior ::NTuple{2,Float64},
                  λ_refd  ::NTuple{2,Float64},
                  nitpp  ::Int64,
                  nthpp  ::Int64,
                  βs     ::Vector{Float64},
                  pup    ::Array{Int64,1},
                  stem   ::Bool)

MCMC da chain for constant birth-death using forward simulation.
"""
function ref_posterior(Ξ      ::Vector{sTpb},
                       idf    ::Array{iBffs,1},
                       λc     ::Float64,
                       λ_prior::NTuple{2,Float64},
                       λ_refd ::NTuple{2,Float64},
                       nitpp  ::Int64,
                       nthpp  ::Int64,
                       βs     ::Vector{Float64},
                       pup    ::Array{Int64,1},
                       stem   ::Bool)

  K = lastindex(βs)

  # make log-likelihood table per power
  nlg = fld(nitpp, nthpp)
  pp  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

  el = lastindex(idf)
  ns = Float64(nnodesinternal(Ξ))
  L  = treelength(Ξ)

  nsi = stem ? 0.0 : log(λc)

  llc = llik_cpb(Ξ, λc) - nsi + prob_ρ(idf)
  prc = logdgamma(λc, λ_prior[1], λ_prior[2])

  for k in 2:K

    βi  = βs[k]
    rdc = logdgamma(λc, λ_refd[1], λ_refd[2])

    # logging
    lth, lit = 0, 0

    for it in Base.OneTo(nitpp)

      shuffle!(pup)

      for p in pup

        # λ proposal
        if p === 1

          llc, prc, rdc, λc =
            update_λ!(llc, prc, rdc, λc, ns, L, stem, λ_prior, λ_refd, βi)

        # forward simulation proposal proposal
        else

          bix = ceil(Int64,rand()*el)
          llc, ns, L = update_fs!(bix, Ξ, idf, llc, λc, ns, L)

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
    update_fs!(bix    ::Int64,
               Ξ      ::Vector{sTpb},
               idf    ::Vector{iBffs},
               llc    ::Float64,
               λ      ::Float64,
               ns     ::Float64,
               L      ::Float64,
               pow    ::Float64)

Forward simulation proposal function for constant pure-birth.
"""
function update_fs!(bix    ::Int64,
                    Ξ      ::Vector{sTpb},
                    idf    ::Vector{iBffs},
                    llc    ::Float64,
                    λ      ::Float64,
                    ns     ::Float64,
                    L      ::Float64)

  bi = idf[bix]

  if iszero(d1(bi)) # is it terminal
    ξp, llr = fsbi_t(bi, λ)
  else
    ξp, llr = fsbi_i(bi, λ)
  end

  if isfinite(llr)
    ξc  = Ξ[bix]

    # update llc, ns & L
    llc += llik_cpb(ξp, λ) - llik_cpb(ξc, λ) + llr
    ns  += Float64(nnodesinternal(ξp) - nnodesinternal(ξc))
    L   += treelength(ξp)             - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ns, L
end




"""
    fsbi(bi::iBffs, λ::Float64)

Forward simulation for terminal branch.
"""
function fsbi_t(bi::iBffs, λ::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - 
       Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, nap, nn, llr =
    _sim_cpb_t(e(bi), λ, lc, lU, Iρi, 0, 1, 500)

  if isnan(llr) || nn >= 500
    return t0, -Inf
  else

    _fixrtip!(t0, nap) # fix random tip
    setni!(bi, nap)    # set new ni

    return t0, llr
  end
end




"""
    fsbi(bi::iBffs, λ::Float64, ntry::Int64)

Forward simulation for internal branch.
"""
function fsbi_i(bi::iBffs, λ::Float64)

  t0, na = _sim_cpb_i(e(bi), λ, 1, 500)

  if na >= 500
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

    # simulated remaining tips until the present
    if na > 1
      t0, na, acr = tip_sims!(t0, tf(bi), λ, acr, lU, Iρi, na)
    end

    if lU < acr
      na -= 1
      llr = (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      setnt!(bi, ntp)                    # set new nt
      setni!(bi, na)                 # set new ni

      return t0, llr
    end
  end

  return t0, NaN
end




"""
    tip_sims!(tree::sTpb,
              t   ::Float64,
              λ   ::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::sTpb,
                   t   ::Float64,
                   λ   ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64)

  if na < 500 && lU < lr

    if istip(tree)
      if !isfix(tree)

        # simulate
        stree, na, lr = _sim_cpb_it(t, λ, lr, lU, Iρi, na, 500)

        if isnan(lr) || na >= 500
          return tree, na, NaN
        end

        # merge to current tip
        sete!(tree, e(tree) + e(stree))
        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, lr = tip_sims!(tree.d1, t, λ, lr, lU, Iρi, na)
      tree.d2, na, lr = tip_sims!(tree.d2, t, λ, lr, lU, Iρi, na)
    end

    return tree, na, lr
  end

  return tree, na, NaN
end




"""
     update_λ!(llc    ::Float64,
               prc    ::Float64,
               λc     ::Float64,
               ns     ::Float64,
               L      ::Float64,
               stem   ::Bool
               λ_prior::NTuple{2,Float64})

Gibbs sampling of `λ` for constant pure-birth.
"""
function update_λ!(llc    ::Float64,
                   prc    ::Float64,
                   λc     ::Float64,
                   ns     ::Float64,
                   L      ::Float64,
                   stem   ::Bool,
                   λ_prior::NTuple{2,Float64})

  nsi = stem ? 0.0 : 1.0

  λp  = randgamma(λ_prior[1] + ns - nsi, λ_prior[2] + L)

  llc += (ns - nsi) * log(λp/λc) + L * (λc - λp)
  prc += llrdgamma(λp, λc, λ_prior[1], λ_prior[2])

  return llc, prc, λp
end




"""
    update_λ!(llc   ::Float64,
              prc   ::Float64,
              rdc   ::Float64,
              λc    ::Float64,
              ns    ::Float64,
              L     ::Float64,
              λ_prior::NTuple{2,Float64},
              λ_refd ::NTuple{2,Float64},
              pow   ::Float64)

Gibbs sampling of `λ` for constant pure-birth with reference distribution.
"""
function update_λ!(llc    ::Float64,
                   prc    ::Float64,
                   rdc    ::Float64,
                   λc     ::Float64,
                   ns     ::Float64,
                   L      ::Float64,
                   stem   ::Bool,
                   λ_prior::NTuple{2,Float64},
                   λ_refd ::NTuple{2,Float64},
                   pow    ::Float64)

  nsi = stem ? 0.0 : 1.0

  λp = randgamma((λ_prior[1] + ns - nsi)*pow + λ_refd[1] * (1.0 - pow),
                 (λ_prior[2] + L)*pow        + λ_refd[2] * (1.0 - pow))

  llc += (ns - nsi)*log(λp/λc) + L*(λc - λp)
  prc += llrdgamma(λp, λc, λ_prior[1], λ_prior[2])
  rdc += llrdgamma(λp, λc, λ_refd[1],  λ_refd[2])

  return llc, prc, rdc, λp
end




"""
  write_ssr(R       ::Array{Float64,2},
            pardic  ::Dict{String,Int64},
            out_file::String)

Write the samples from an MC sampler data frame
given a Dictionary of parameters.
"""
function write_ssr(R       ::Array{Float64,2},
                   pardic  ::Dict{String,Int64},
                   out_file::String)

  # column names
  col_nam = ["Iteration", "Likelihood", "Prior"]

  for (k,v) in sort!(collect(pardic), by = x -> x[2])
    push!(col_nam, k)
  end

  R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

  writedlm(out_file*".log", R)
end



