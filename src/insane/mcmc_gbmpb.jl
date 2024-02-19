#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#




"""
    insane_gbmpb(tree    ::sT_label;
                 λ0_prior::NTuple{2,Float64}     = (1.0, 0.0),
                 α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64}     = (3.0, 0.5),
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 nflush  ::Int64                 = nthin,
                 tune_int::Int64                 = 100,
                 ofile   ::String                = string(homedir(), "/ipb"),
                 αi      ::Float64               = 0.0,
                 σλi     ::Float64               = 0.1,
                 λtni    ::Float64               = 0.1,
                 obj_ar  ::Float64               = 0.234,
                 pupdp   ::NTuple{5,Float64}     = (0.01, 0.01, 0.1, 0.2),
                 δt      ::Float64               = 1e-3,
                 prints  ::Int64                 = 5,
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for `gbm-pb`.
"""
function insane_gbmpb(tree    ::sT_label;
                      λ0_prior::NTuple{2,Float64}     = (1.0, 0.0),
                      α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64}     = (3.0, 0.5),
                      niter   ::Int64                 = 1_000,
                      nthin   ::Int64                 = 10,
                      nburn   ::Int64                 = 200,
                      nflush  ::Int64                 = nthin,
                      tune_int::Int64                 = 100,
                      ofile   ::String                = string(homedir(), "/ipb"),
                      αi      ::Float64               = 0.0,
                      σλi     ::Float64               = 0.1,
                      λtni    ::Float64               = 0.1,
                      obj_ar  ::Float64               = 0.234,
                      pupdp   ::NTuple{5,Float64}     = (0.01, 0.01, 0.1, 0.2),
                      δt      ::Float64               = 1e-3,
                      prints  ::Int64                 = 5,
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf = make_idf(tree, tρ, Inf)

  # make a decoupled tree
  Ξ = make_Ξ(idf, λmle_cpb(tree), αi, σλi, δt, srδt, iTpb)

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf))  if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σ, 3: λ, 4: gbm, 5: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(5)
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  # make objecting scaling function for tuning
  scalef = makescalef(obj_ar)

  @info "running pure-birth gbm"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, λtn, ns =
    mcmc_burn_gbmpb(Ξ, idf, λ0_prior, α_prior, σλ_prior, nburn, tune_int, αi, σλi, λtni,
      δt, srδt, inodes, pup, prints, scalef)

  # mcmc
  r, treev = mcmc_gbmpb(Ξ, idf, llc, prc, αc, σλc, ns, λ0_prior, α_prior, σλ_prior, λtn,
      δt, srδt, inodes, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end



"""
    mcmc_burn_gbmpb(Ξ       ::Vector{iTpb},
                    idf     ::Vector{iBffs},
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    tune_int::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    λtn     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    prints  ::Int64,
                    scalef  ::Function)


MCMC burn-in chain for GBM pure-birth.
"""
function mcmc_burn_gbmpb(Ξ       ::Vector{iTpb},
                         idf     ::Vector{iBffs},
                         λ0_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         tune_int::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         λtn     ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         inodes  ::Array{Int64,1},
                         pup     ::Array{Int64,1},
                         prints  ::Int64,
                         scalef  ::Function)

  ltn = 0
  λlup = λlac = 0.0

  nsi = Float64(iszero(e(Ξ[1])))

  # starting likelihood and prior
  lλ0 = lλ(Ξ[1])[1]
  llc = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) - nsi*lλ0 + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2]) +
        logdgamma(exp(lλ0), λ0_prior[1], λ0_prior[2])  +
        logdnorm(αc,        α_prior[1], σλc^2)

  L       = treelength(Ξ)          # tree length
  dλ      = deltaλ(Ξ)              # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)         # sum squares in λ
  ns      = nnodesbifurcation(idf) # number of speciation events
  nin     = lastindex(inodes)      # number of internal nodes
  el      = lastindex(idf)         # number of branches

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      
      # # update drift
      # if pupi === 1

      #   #llc, prc, αc = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

      #   # update ssλ with new drift `α`
      #   ssλ, nλ = sss_gbm(Ξ, αc)

      # # update diffusion
      # elseif pupi === 2

      #   llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior, α_prior)

      # update drift and diffusion
      if pupi === 1 || pupi === 2

        #llc, prc, αc = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)
        llc, prc, αc, σλc = update_α_σ!(αc, σλc, L, dλ, ssλ, nλ, llc, prc, α_prior, σλ_prior)

        # update ssλ with new drift `α`
        ssλ, nλ = sss_gbm(Ξ, αc)

      # update all speciation rates through time simultaneously
      elseif pupi === 3

        llc, prc, Ξ, λlac = update_lλ!(Ξ, llc, prc, ns, nsi, λtn, λlac, δt, λ0_prior)

        λlup += 1.0 

      # update gbm
      elseif pupi === 4

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, dλ, ssλ =
          update_gbm!(bix, Ξ, idf, αc, σλc, llc, prc, dλ, ssλ, δt, srδt, λ0_prior)

      # forward simulation
      else
        bix = ceil(Int64,rand()*el)

        llc, dλ, ssλ, nλ, ns, L =
          update_fs!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, nλ, ns, L, δt, srδt)
      end
    end

    # log tuning parameters
    ltn += 1
    if ltn === tune_int
      λtn = scalef(λtn, λlac/λlup)
      ltn = 0
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, λtn, ns
end




"""
    mcmc_gbmpb(Ξ       ::Vector{iTpb},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               ns      ::Float64,
               λ0_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               λtn     ::Float64,
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Vector{Int64},
               niter   ::Int64,
               nthin   ::Int64,
               nflush  ::Int64,
               ofile   ::String,
               prints  ::Int64)

MCMC chain for GBM pure-birth.
"""
function mcmc_gbmpb(Ξ       ::Vector{iTpb},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    ns      ::Float64,
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    λtn     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Vector{Int64},
                    niter   ::Int64,
                    nthin   ::Int64,
                    nflush  ::Int64,
                    ofile   ::String,
                    prints  ::Int64)

  nsi = Float64(iszero(e(Ξ[1])))

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  r = Array{Float64,2}(undef, nlogs, 6)

  # make Ξ vector
  treev = iTpb[]

  L       = treelength(Ξ)      # tree length
  dλ      = deltaλ(Ξ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ξ, αc)     # sum squares in λ
  nin     = lastindex(inodes)  # number of internal nodes
  el      = lastindex(idf)     # number of branches

  # flush to file
  sthin = 0

  function check_pr(pupi::Int64, it::Int64)
   pr0 = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2]) +
         logdgamma(exp(lλ(Ξ[1])[1]), λ0_prior[1], λ0_prior[2])  +
         logdnorm(αc,        α_prior[1], σλc^2)
   if !isapprox(pr0, prc, atol = 1e-5)
      error(string("Wrong prior computation during the ", ["α","σλ","λ","gbm","forward simulation"][pupi], 
                   " update, at iteration ", it, ": pr0=", pr0, " and prc=", prc))
   end
  end

  function check_ll(pupi::Int64, it::Int64)
   ll0 = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) - nsi*lλ(Ξ[1])[1] + prob_ρ(idf)
   if !isapprox(ll0, llc, atol = 1e-5)
      error(string("Wrong likelihood computation during the ", ["α","σλ","λ","gbm","forward simulation"][pupi], 
                   " update, at iteration ", it, ": ll0=", ll0, " and llc=", llc))
   end
  end

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\talpha\tsigma_lambda\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      pbar = Progress(niter, prints, "running mcmc...", 20)

      for it in Base.OneTo(niter)

        shuffle!(pup)

        for pupi in pup

          ## parameter updates
          # update drift
          #if pupi === 1
          # update drift and diffusion
          if pupi === 1 || pupi === 2

            #llc, prc, αc = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)
            llc, prc, αc, σλc = update_α_σ!(αc, σλc, L, dλ, ssλ, nλ, llc, prc, α_prior, σλ_prior)

            # update ssλ with new drift `α`
            ssλ, nλ = sss_gbm(Ξ, αc)

          # # update diffusion rate
          # elseif pupi === 2
          #   llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior, α_prior)

          # update all speciation rates through time simultaneously
          elseif pupi === 3
            llc, prc, Ξ = update_lλ!(Ξ, llc, prc, ns, nsi, λtn, δt, λ0_prior)

          # update gbm
          elseif pupi === 4
            nix = ceil(Int64,rand()*nin)
            bix = inodes[nix]

            llc, prc, dλ, ssλ =
              update_gbm!(bix, Ξ, idf, αc, σλc, llc, prc, dλ, ssλ, δt, srδt, λ0_prior)

          # update by forward simulation
          else
            bix = ceil(Int64,rand()*el)

            llc, dλ, ssλ, nλ, ns, L =
              update_fs!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, nλ, ns, L, δt, srδt)

          end
          #check_pr(pupi, it)
          #check_ll(pupi, it)

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
            r[lit,5] = αc
            r[lit,6] = σλc
            push!(treev, couple(Ξ, idf, 1))
          end
          lthin = 0
        end

        # flush parameters
        sthin += 1
        if sthin === nflush
          write(of, 
            string(Float64(it), "\t", llc, "\t", prc, "\t", 
              exp(lλ(Ξ[1])[1]),"\t", αc, "\t", σλc, "\n"))
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





# """
#     ref_posterior(Ξ       ::Vector{iTpb},
#                   idf     ::Vector{iBffs},
#                   llc     ::Float64,
#                   prc     ::Float64,
#                   αc      ::Float64,
#                   σλc     ::Float64,
#                   λ0_prior::NTuple{2,Float64},
#                   α_prior ::NTuple{2,Float64},
#                   σλ_prior::NTuple{2,Float64},
#                   λ0rdist::NTuple{2,Float64},
#                   α_rdist ::NTuple{2,Float64},
#                   σλ_rdist::NTuple{2,Float64},
#                   nitpp   ::Int64,
#                   nthpp   ::Int64,
#                   βs      ::Vector{Float64},
#                   δt      ::Float64,
#                   srδt    ::Float64,
#                   inodes  ::Array{Int64,1},
#                   pup     ::Array{Int64,1},
#                   prints  ::Int64)

# MCMC chain for GBM pure-birth.
# """
# function ref_posterior(Ξ       ::Vector{iTpb},
#                        idf     ::Vector{iBffs},
#                        llc     ::Float64,
#                        prc     ::Float64,
#                        αc      ::Float64,
#                        σλc     ::Float64,
#                        α_prior ::NTuple{2,Float64},
#                        σλ_prior::NTuple{2,Float64},
#                        α_rdist ::NTuple{2,Float64},
#                        σλ_rdist::NTuple{2,Float64},
#                        nitpp   ::Int64,
#                        nthpp   ::Int64,
#                        βs      ::Vector{Float64},
#                        δt      ::Float64,
#                        srδt    ::Float64,
#                        inodes  ::Array{Int64,1},
#                        pup     ::Array{Int64,1})

#   # starting likelihood and prior
#   llc = llik_gbm(Ξ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
#   prc = logdnorm(αc,         α_prior[1], α_prior[2]^2) +
#         logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])

#   K = lastindex(βs)

#   # make log-likelihood table per power
#   nlg = fld(nitpp, nthpp)
#   pp  = [Vector{Float64}(undef,nlg) for i in Base.OneTo(K)]

#   L       = treelength(Ξ)      # tree length
#   dλ      = deltaλ(Ξ)          # delta change in λ
#   ssλ, nλ = sss_gbm(Ξ, αc)     # sum squares in λ
#   nin     = lastindex(inodes)  # number of internal nodes
#   el      = lastindex(idf)     # number of branches

#   for k in 2:K

#     βi  = βs[k]
#     rdc = logdnorm(       αc,  α_rdist[1], α_rdist[2]^2) +
#           logdinvgamma(σλc^2, σλ_rdist[1], σλ_rdist[2])

#     # logging
#     lth, lit = 0, 0

#     for it in Base.OneTo(nitpp)

#       shuffle!(pup)

#       for pupi in pup

#         ## parameter updates
#         # update drift
#         if pupi === 1

#           llc, prc, rdc, αc = update_α!(αc, σλc, L, dλ, llc, prc, rdc,
#             α_prior, α_rdist, βi)

#           # update ssλ with new drift `α`
#           ssλ, nλ = sss_gbm(Ξ, αc)

#         # update diffusion rate
#         elseif pupi === 2

#           llc, prc, rdc, σλc = update_σ!(σλc, ssλ, nλ, llc, prc, rdc,
#             σλ_prior, σλ_rdist, βi)

#         # update gbm
#         elseif pupi === 3

#           nix = ceil(Int64,rand()*nin)
#           bix = inodes[nix]

#           llc, dλ, ssλ =
#             update_gbm!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, δt, srδt)

#         # update by forward simulation
#         else
#           bix = ceil(Int64,rand()*el)

#           llc, dλ, ssλ, nλ, L =
#             update_fs!(bix, Ξ, idf, αc, σλc, llc, dλ, ssλ, nλ, L, δt, srδt)

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
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{iTpb},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               llc  ::Float64,
               dλ   ::Float64,
               ssλ  ::Float64,
               nλ   ::Float64,
               ns   ::Float64,
               L    ::Float64,
               δt   ::Float64,
               srδt ::Float64)

Forward simulation proposal function for constant birth-death.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{iTpb},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    llc  ::Float64,
                    dλ   ::Float64,
                    ssλ  ::Float64,
                    nλ   ::Float64,
                    ns   ::Float64,
                    L    ::Float64,
                    δt   ::Float64,
                    srδt ::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]

  # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, ξc, α, σλ, δt, srδt)
    drλ = ssrλ = 0.0
  # if internal
  else
    ξp, llr, drλ, ssrλ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, δt, srδt)
  end

  # if accepted
  if isfinite(llr)
    ll1, dλ1, ssλ1, nλ1 = llik_gbm_ssλ(ξp, α, σλ, δt, srδt)
    ll0, dλ0, ssλ0, nλ0 = llik_gbm_ssλ(ξc, α, σλ, δt, srδt)

    # update llr, ssλ, nλ, ns, L
    llc += ll1  - ll0 + llr
    dλ  += dλ1  - dλ0  + drλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    nλ  += nλ1  - nλ0
    ns  += nnodesinternal(ξp) - nnodesinternal(ξc)
    L   += treelength(ξp) - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, dλ, ssλ, nλ, ns, L
end




"""
    fsbi_t(bi  ::iBffs,
           ξc  ::iTpb,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTpb,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  Iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(Iρi) ? 0.0 : log(Iρi))

  # forward simulation during branch length
  t0, nap, nn, llr =
    _sim_gbmpb_t(e(bi), lλ(ξc)[1], α, σλ, δt, srδt, lc, lU, Iρi, 0, 1, 500)

  if isfinite(llr)
    _fixrtip!(t0, nap) # fix random tip
    setni!(bi, nap)    # set new ni

    return t0, llr
  else
    return t0, NaN
  end
end




"""
    fsbi_i(bi  ::iBffs,
           ξ1  ::iTpb,
           ξ2  ::iTpb,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTpb,
                ξ1  ::iTpb,
                ξ2  ::iTpb,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  t0, na = _sim_gbmpb(e(bi), lλ(ξc)[1], α, σλ, δt, srδt, 1, 1_000)

  if na >= 1_000
    return t0, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  Iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(Iρi) ? 0.0 : log(Iρi))

 # fix random tip
  λf = fixrtip!(t0, na, NaN)

  llrd, acrd, drλ, ssrλ, λ1p, λ2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulated remaining tips until the present
    t0, na, acr =
      tip_sims!(t0, tf(bi), α, σλ, δt, srδt, acr, lU, Iρi, na)

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(Iρi) ? 0.0 : log(Iρi))
      l1  = lastindex(λ1p)
      l2  = lastindex(λ2p)
      setnt!(bi, ntp)                    # set new nt
      setni!(bi, na)                     # set new ni
      setλt!(bi, λf)                     # set new λt
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 2 λ vector

      return t0, llr, drλ, ssrλ
    else
      return t0, NaN, NaN, NaN
    end
  end

  return t0, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTpb,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              Iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTpb,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   Iρi ::Float64,
                   na  ::Int64)

 if lU < lr && na < 1_000

    if istip(tree)
      if !isfix(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)

        # simulate
        stree, na, lr =
          _sim_gbmpb_it(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, δt, srδt,
            lr, lU, Iρi, na, 1_000)

        if isnan(lr) || na >= 1_000
          return tree, na, NaN
        end

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
      tree.d1, na, lr = tip_sims!(tree.d1, t, α, σλ, δt, srδt, lr, lU, Iρi, na)
      tree.d2, na, lr = tip_sims!(tree.d2, t, α, σλ, δt, srδt, lr, lU, Iρi, na)
    end

    return tree, na, lr
  end

  return tree, na, NaN
end




"""
    update_gbm!(bix     ::Int64,
                Ξ       ::Vector{iTpb},
                idf     ::Vector{iBffs},
                α       ::Float64,
                σλ      ::Float64,
                llc     ::Float64,
                prc     ::Float64,
                dλ      ::Float64,
                ssλ     ::Float64,
                δt      ::Float64,
                srδt    ::Float64,
                λ0_prior::NTuple{2,Float64})

Make a `gbm` update for an interna branch and its descendants.
"""
function update_gbm!(bix     ::Int64,
                     Ξ       ::Vector{iTpb},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     dλ      ::Float64,
                     ssλ     ::Float64,
                     δt      ::Float64,
                     srδt    ::Float64,
                     λ0_prior::NTuple{2,Float64})

  ξi   = Ξ[bix]
  bi   = idf[bix]
  i1   = d1(bi)
  i2   = d2(bi)
  ξ1   = Ξ[i1]
  ξ2   = Ξ[i2]
  root = iszero(pa(bi))

  # if crown root
  if root && iszero(e(ξi))
    llc, prc, dλ, ssλ =
      _crown_update!(ξi, ξ1, ξ2, α, σλ, llc, prc, dλ, ssλ, δt, srδt, λ0_prior)
    setλt!(bi, lλ(ξi)[1])
  else
    # if stem
    if root
      llc, prc, dλ, ssλ = _stem_update!(ξi, α, σλ, llc, prc, dλ, ssλ, δt, srδt, λ0_prior)
    end

    # updates within the parent branch
    llc, dλ, ssλ = _update_gbm!(ξi, α, σλ, llc, dλ, ssλ, δt, srδt, false)

    # get fixed tip
    lξi = fixtip(ξi)

    # make between decoupled trees node update
    llc, dλ, ssλ = update_triad!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2),
      fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, llc, dλ, ssλ, δt, srδt)

    # set fixed `λ(t)` in branch
    setλt!(bi, lλ(lξi)[end])
  end

  # # carry on updates in the daughters
  llc, dλ, ssλ = 
    _update_gbm!(ξ1, α, σλ, llc, dλ, ssλ, δt, srδt, iszero(d1(idf[i1])))
  llc, dλ, ssλ = 
    _update_gbm!(ξ2, α, σλ, llc, dλ, ssλ, δt, srδt, iszero(d1(idf[i2])))

  return llc, prc, dλ, ssλ
end





"""
    update_α!(αc     ::Float64,
              σλ     ::Float64,
              L      ::Float64,
              dλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              α_prior::NTuple{2,Float64})

Gibbs update for `α`.
"""
function update_α!(αc     ::Float64,
                   σλ     ::Float64,
                   L      ::Float64,
                   dλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   α_prior::NTuple{2,Float64})

  # ratio
  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  # τ2  = σλ^2
  σλ2 = σλ^2
  rs  = σλ2/τ2

  # gibbs update for σ
  αp = rnorm((dλ + rs*ν)/(rs + L), sqrt(σλ2/(rs + L)))

  # update prior
  prc += llrdnorm_x(αp, αc, ν, τ2)

  # update likelihood
  llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/L)

  return llc, prc, αp
end




"""
    update_α!(αc     ::Float64,
              σλ     ::Float64,
              L      ::Float64,
              dλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              rdc    ::Float64,
              α_prior::NTuple{2,Float64},
              α_rdist ::NTuple{2,Float64},
              pow    ::Float64)


Gibbs update for `α` given reference distribution.
"""
function update_α!(αc     ::Float64,
                   σλ     ::Float64,
                   L      ::Float64,
                   dλ     ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   rdc    ::Float64,
                   α_prior::NTuple{2,Float64},
                   α_rdist ::NTuple{2,Float64},
                   pow    ::Float64)

  # ratio
  ν   = α_prior[1]
  #τ2  = α_prior[2]^2
  τ2  = σλ^2
  σλ2 = σλ^2
  rs  = σλ2/τ2

  cpow = (1.0 - pow)

  # gibbs update for α
  m   = (dλ + rs*ν)/(rs + L)
  s2  = σλ2/(rs + L)
  m0  = α_rdist[1]
  s02 = α_rdist[2]^2
  αp  = rnorm((m0 * s2 * cpow + m * s02 * pow) / (pow * s02 + s2 * cpow),
              sqrt( s2 * s02 /  (pow * s02 + s2 * cpow)) )

  # update likelihood, prior and reference
  llc += 0.5*L/σλ2*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/L)
  prc += llrdnorm_x(αp, αc, ν, τ2)
  rdc += llrdnorm_x(αp, αc, m0, s02)

  return llc, prc, rdc, αp
end




"""
    update_σ!(σλc     ::Float64,
              αc       ::Float64,
              ssλ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64},
              α_prior ::NTuple{2,Float64})

Gibbs update for `σλ`.
"""
function update_σ!(σλc     ::Float64,
                   αc       ::Float64,
                   ssλ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   σλ_prior::NTuple{2,Float64},
                   α_prior ::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)

  # update prior
  prc += llrdinvgamma(σλp2, σλc^2, σλ_p1, σλ_p2)
  prc += llrdnorm_σ²(αc, α_prior[1], σλp2, σλc^2)

  σλp = sqrt(σλp2)

  # update likelihood
  llc += ssλ*(1.0/σλc^2 - 1.0/σλp2) - n*(log(σλp/σλc))

  return llc, prc, σλp
end




# """
#     update_σ!(σc     ::Float64,
#               ss     ::Float64,
#               n       ::Float64,
#               llc     ::Float64,
#               prc     ::Float64,
#               σλ_prior::NTuple{2,Float64})

# Gibbs update for `σλ` given reference distribution.
# """
# function update_σ!(σc     ::Float64,
#                    ss     ::Float64,
#                    n       ::Float64,
#                    llc     ::Float64,
#                    prc     ::Float64,
#                    rdc     ::Float64,
#                    σλ_prior::NTuple{2,Float64},
#                    σλ_rdist ::NTuple{2,Float64},
#                    pow     ::Float64)

#   σλ_p1 = σλ_prior[1]
#   σλ_p2 = σλ_prior[2]

#   # Gibbs update for σ
#   σλp2 = randinvgamma((σλ_p1 + 0.5 * n) * pow + σλ_rdist[1] * (1.0 - pow),
#                       (σλ_p2 + ss) * pow     + σλ_rdist[2] * (1.0 - pow))

#   # update likelihood, prior and reference
#   σλp = sqrt(σλp2)
#   llc += ss*(1.0/σc^2 - 1.0/σλp2) - n*(log(σλp/σc))
#   prc += llrdinvgamma(σλp2, σλc^2, σλ_p1, σλ_p2)
#   rdc += llrdinvgamma(σλp2, σλc^2, σλ_rdist[1], σλ_rdist[2])

#   return llc, prc, rdc, σλp
# end




"""
    update_α_σ!(αc      ::Float64,
                σλc     ::Float64,
                L       ::Float64,
                dλ      ::Float64,
                ssλ     ::Float64,
                n       ::Float64,
                llc     ::Float64,
                prc     ::Float64,
                α_prior ::NTuple{2,Float64},
                σλ_prior::NTuple{2,Float64})

Gibbs update for `α`.
"""
function update_α_σ!(αc      ::Float64,
                     σλc     ::Float64,
                     L       ::Float64,
                     dλ      ::Float64,
                     ssλ     ::Float64,
                     n       ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     α_prior ::NTuple{2,Float64},
                     σλ_prior::NTuple{2,Float64})

  # ratio
  ν   = α_prior[1]
  σλ2c = σλc^2
  rs  = 1.0

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]
  ν     = α_prior[1]

  # gibbs update for α and σ
  αp, σλ2p = randnorminvgamma((dλ + rs*ν)/(rs + L), 
                              rs + L,
                              σλ_p1 + 0.5 * n,
                              σλ_p2 + ssλ + rs*L/(rs + L)*(ν-dλ/n)^2)

  # update prior for α
  prc += llrdnorm_σ²(αc, α_prior[1], σλ2p, σλ2c) + llrdnorm_x(αp, αc, ν, σλ2p)
  # update prior for σ
  prc += llrdinvgamma(σλ2p, σλ2c, σλ_p1, σλ_p2)

  # update likelihood for α
  llc += 0.5*L/σλ2p*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/L)
  # update likelihood for σ
  σλp = sqrt(σλ2p)
  llc += ssλ*(1.0/σλ2c - 1.0/σλ2p) - n*(log(σλp/σλc))

  return llc, prc, αp, σλp
end




"""
     update_lλ!(Ξ       ::Vector{iTpb},
                llc     ::Float64,
                prc     ::Float64,
                ns      ::Float64,
                nsi     ::Float64,
                λtn     ::Float64,
                δt      ::Float64,
                λ0_prior::NTuple{2,Float64})

Mixed HM-Gibbs sampling, shifting all log-`λ` GBM rates simultaneously.
"""
function update_lλ!(Ξc      ::Vector{iTpb},
                    llc     ::Float64,
                    prc     ::Float64,
                    ns      ::Float64,
                    nsi     ::Float64,
                    λtn     ::Float64,
                    δt      ::Float64,
                    λ0_prior::NTuple{2,Float64})
  
  lλ0c = lλ(Ξc[1])[1]
  λ0c  = exp(lλ0c)
  lλ0p = rnorm(lλ0c, λtn)
  λ0p  = exp(lλ0p)
    
  lλshift = lλ0p-lλ0c

  llr = llr_gbm_lλshift(Ξc, δt, lλshift) + (ns-nsi)*lλshift
  prr = llrdgamma(λ0p, λ0c, λ0_prior[1], λ0_prior[2])
    
  if -randexp() < llr + prr
    llc += llr
    prc += prr
    for i in Base.OneTo(lastindex(Ξc))
      propagate_lλshift!(Ξc[i], lλshift)
    end
  end

  return llc, prc, Ξc
end




"""
     update_lλ!(Ξ       ::Vector{iTpb},
                llc     ::Float64,
                prc     ::Float64,
                ns      ::Float64,
                nsi     ::Float64,
                λtn     ::Float64,
                lac     ::Float64,
                δt      ::Float64,
                λ0_prior::NTuple{2,Float64})

Mixed HM-Gibbs sampling, shifting all log-`λ` GBM rates simultaneously.
"""
function update_lλ!(Ξc      ::Vector{iTpb},
                    llc     ::Float64,
                    prc     ::Float64,
                    ns      ::Float64,
                    nsi     ::Float64,
                    λtn     ::Float64,
                    lac     ::Float64,
                    δt      ::Float64,
                    λ0_prior::NTuple{2,Float64})
  
  lλ0c = lλ(Ξc[1])[1]
  λ0c  = exp(lλ0c)
  lλ0p = rnorm(lλ0c, λtn)
  λ0p  = exp(lλ0p)
    
  lλshift = lλ0p-lλ0c

  llr = llr_gbm_lλshift(Ξc, δt, lλshift) + (ns-nsi)*lλshift
  prr = llrdgamma(λ0p, λ0c, λ0_prior[1], λ0_prior[2])
    
  if -randexp() < llr + prr
    llc += llr
    prc += prr
    for i in Base.OneTo(lastindex(Ξc))
      propagate_lλshift!(Ξc[i], lλshift)
    end
    lac += 1.0
  end

  return llc, prc, Ξc, lac
end




"""
    propagate_lλshift!(tree   ::iTpb,
                       lλshift::Float64)

Propagate a shift in log-speciation across all GBM rates.
"""
function propagate_lλshift!(tree   ::T,
                            lλshift::Float64) where {T <: iT}
  setlλ!(tree, lλ(tree) .+ lλshift)
  if !istip(tree)
    propagate_lλshift!(tree.d1, lλshift)
    propagate_lλshift!(tree.d2, lλshift)
  end
end


