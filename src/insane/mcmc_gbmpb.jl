#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#





"""
    insane_gbmpb(tree    ::sT_label, 
                 out_file::String;
                 δt      ::Float64               = 1e-2,
                 niter   ::Int64                 = 1_000,
                 nthin   ::Int64                 = 10,
                 nburn   ::Int64                 = 200,
                 tune_int::Int64                 = 100,
                 σλi     ::Float64               = 0.1,
                 αi      ::Float64               = 0.0,
                 prints  ::Int64                 = 5,
                 pupdp   ::NTuple{4,Float64}     = (0.2, 0.2, 0.3, 0.3),
                 λa_prior::NTuple{2,Float64}     = (0.0, 100.0),
                 α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                 σλ_prior::NTuple{2,Float64}     = (0.05, 0.5),
                 tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for GBM pure-birth.
"""
function insane_gbmpb(tree    ::sT_label, 
                      out_file::String;
                      δt      ::Float64               = 1e-2,
                      niter   ::Int64                 = 1_000,
                      nthin   ::Int64                 = 10,
                      nburn   ::Int64                 = 200,
                      tune_int::Int64                 = 100,
                      σλi     ::Float64               = 0.1,
                      αi      ::Float64               = 0.0,
                      prints  ::Int64                 = 5,
                      pupdp   ::NTuple{4,Float64}     = (0.2, 0.2, 0.3, 0.3),
                      λa_prior::NTuple{2,Float64}     = (0.0, 100.0),
                      α_prior ::NTuple{2,Float64}     = (0.0, 10.0),
                      σλ_prior::NTuple{2,Float64}     = (0.05, 0.5),
                      tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  δt  *= treeheight(tree)
  srδt = sqrt(δt)
  n    = ntips(tree)

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # lλ root node
  lλa = log(λmle_cpb(tree))

  # make fix tree directory
  idf = make_idf(tree, tρ)

  # make an edges tree
  Ψ = iTgbmpb[]
  iTgbmpb!(Ψ, tree, δt, srδt, lλa, αi, σλi)

  # set end of fix branch speciation times and
  # get vector of internal branches
  inodes = Int64[]
  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    setλt!(bi, lλ(Ψ[i])[end])
    if !iszero(d1(bi))
      push!(inodes, i)
    end
  end

  # parameter updates (1: α, 2: σ, 3: gbm, 4: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(4) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running pure-birth gbm"

  # burn-in phase
  llc, prc, αc, σλc =
    mcmc_burn_gbmpb(Ψ, λa_prior, α_prior, σλ_prior, nburn, αi, σλi, 
      δt, srδt, idv, inodes, terminus, pup, prints)

  # mcmc
  R, Ψv = mcmc_gbmpb(Ψp, Ψc, llc, prc, αc, σλc, λa_prior, α_prior, σλ_prior, 
        niter, nthin, δt, srδt, idv, inodes, terminus, pup, prints)

  pardic = Dict(("lambda_root"  => 1,
                 "alpha" => 2,
                 "sigma_lambda" => 3))

  write_ssr(R, pardic, out_file)

  return R, Ψv
end




"""
    mcmc_burn_gbmpb(Ψp      ::iTgbmpb,
                    Ψc      ::iTgbmpb,
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    nburn   ::Int64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idv     ::Array{iBfb,1},
                    inodes  ::Array{Int64,1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

MCMC burn-in chain for GBM pure-birth.
"""
function mcmc_burn_gbmpb(Ψ       ::Vector{iTgbmpb},
                         λa_prior::NTuple{2,Float64},
                         α_prior ::NTuple{2,Float64},
                         σλ_prior::NTuple{2,Float64},
                         nburn   ::Int64,
                         αc      ::Float64,
                         σλc     ::Float64,
                         δt      ::Float64,
                         srδt    ::Float64,
                         idf     ::Array{iBffs,1},
                         inodes  ::Array{Int64,1},
                         terminus::Array{BitArray{1}},
                         pup     ::Array{Int64,1},
                         prints  ::Int64)

  # starting likelihood and prior
  llc = llik_gbm(Ψ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])      + 
        logdnorm(αc, α_prior[1], α_prior[2]^2)             +
        logdunif(exp(lλ(Ψ[1])[1]), λa_prior[1], λa_prior[2])

  # maximum bound in log space for uniform
  lλmxpr = log(λa_prior[2])

  L       = treelength(Ψ)      # tree length
  dλ      = deltaλ(Ψ)          # delta change in λ
  ssλ, nλ = sss_gbm(Ψ, αc)     # sum squares in λ
  nin     = lastindex(inodes)  # number of internal nodes

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update drift
      if pupi === 1

        llc, prc, αc = update_α!(αc, σλc, L, dλ, llc, prc, α_prior)

      # update diffusion
      elseif pupi === 2

        llc, prc, σλc = update_σ!(σλc, αc, ssλ, nλ, llc, prc, σλ_prior)

      # update gbm
      elseif pupi === 3

        """
        here
        """

        bix = ceil(Int64,rand()*nin)

        bix = ceil(Int64,rand()*nin)




        llc = update_gbm!(Ψp, Ψc, llc, αc, σλc, δt, srδt, lλmxpr, 
                dri, ldr, ter, 0)

      # forward simulation
      else



      end

    end

    next!(pbar)
  end

  return llc, prc, αc, σλc
end




"""
    mcmc_gbmpb(Ψp      ::iTgbmpb,
               Ψc      ::iTgbmpb,
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               λa_prior::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               niter   ::Int64,
               nthin   ::Int64,
               δt      ::Float64,
               srδt    ::Float64,
               idv     ::Array{iBfb,1},
               inodes  ::Array{Int64,1},
               terminus::Array{BitArray{1}},
               pup     ::Array{Int64,1},
               prints  ::Int64)

MCMC chain for GBM pure-birth.
"""
function mcmc_gbmpb(Ψp      ::iTgbmpb,
                    Ψc      ::iTgbmpb,
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    λa_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    niter   ::Int64,
                    nthin   ::Int64,
                    δt      ::Float64,
                    srδt    ::Float64,
                    idv     ::Array{iBfb,1},
                    inodes  ::Array{Int64,1},
                    terminus::Array{BitArray{1}},
                    pup     ::Array{Int64,1},
                    prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin, lit = 0, 0

  R = Array{Float64,2}(undef, nlogs, 6)

  # make Ψ vector
  Ψv = iTgbmpb[]

  ntr = lastindex(inodes)

  # prior
  lλmxpr = log(λa_prior[2])

  pbar = Progress(niter, prints, "running mcmc...", 20)

  for it in Base.OneTo(niter)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      if pupi === 1

        """
        here
        """

        # update drift
        llc, prc, αc = update_α!(αc, σλc, Ψc, llc, prc, α_prior)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, pupi, it
        #    return 
        # end

      elseif pupi === 2

        # update diffusion
        llc, prc, σλc = update_σ!(σλc, αc, Ψc, llc, prc, σλ_prior)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, pupi, it
        #    return 
        # end

      else 

        bix = ceil(Int64,rand()*ntr)

        ter  = terminus[tix]
        inod = inodes[tix]

        dri = dr(idv[inod])
        ldr = lastindex(dri)

        llc = update_gbm!(Ψp, Ψc, llc, αc, σλc, δt, srδt, lλmxpr, 
                dri, ldr, ter, 0)

        # ll0 = llik_gbm(Ψ, idf, αc, σλc, δt, srδt) + prob_ρ(idf)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, pupi, it
        #    return 
        # end

      end
    end

    # log parameters
    lthin += 1
    if lthin === nthin
      lit += 1
      @inbounds begin
        R[lit,1] = Float64(lit)
        R[lit,2] = llc
        R[lit,3] = prc
        R[lit,4] = exp(lλ(Ψc)[1])
        R[lit,5] = αc
        R[lit,6] = σλc
        push!(Ψv, deepcopy(Ψc))
      end
      lthin = 0
    end

    next!(pbar)
  end

  return R, Ψv
end




"""
    update_gbm!(bix ::Int64,
                psi ::Vector{iTgbmpb},
                idf ::Vector{iBffs},
                α   ::Float64,
                σλ  ::Float64,
                llc ::Float64,
                ssλ ::Float64,
                δt  ::Float64,
                srδt::Float64)

Make a `gbm` update for an interna branch and its descendants.
"""
function update_gbm!(bix ::Int64,
                     psi ::Vector{iTgbmpb},
                     idf ::Vector{iBffs},
                     α   ::Float64,
                     σλ  ::Float64,
                     llc ::Float64,
                     ssλ ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  ψi = psi[bix]
  bi = idf[bix]

  # if root
  if iszero(pa(bi))
    # updates within the parent branch
    llc, ssλ, lψi = _update_gbm!(ψi, α, σλ, llc, ssλ, δt, srδt, false)
  else
    # updates within the parent branch
    llc, ssλ, lψi = _update_gbm!(ψi, α, σλ, llc, ssλ, δt, srδt, false)
  end

  # make between decoupled trees node update
  ψ1   = psi[d1(bi)]
  ψ2   = psi[d2(bi)]
  ter1 = iszero(d1(idf[d1(bi)])) 
  ter2 = iszero(d1(idf[d2(bi)]))

  llc, ssλ = update_triad!(lλ(lψi), lλ(ψ1), lλ(ψ2) e(lψi), e(ψ1), e(ψ2), 
    fdt(lψi), fdt(ψ1), fdt(ψ2), α, σλ, llc, ssλ, δt, srδt)

  # carry on updates in the daughters
  llc, ssλ, tree = _update_gbm!(ψ1, α, σλ, llc, ssλ, δt, srδt, ter1)
  llc, ssλ, tree = _update_gbm!(ψ2, α, σλ, llc, ssλ, δt, srδt, ter2)

  return llc, ssλ
end





"""
    _update_gbm!(tree::iTgbmpb,
                 α   ::Float64,
                 σλ  ::Float64,
                 llc ::Float64,
                 ssλ ::Float64,
                 δt  ::Float64,
                 srδt::Float64)

Do gbm updates on a decoupled tree recursively.
"""
function _update_gbm!(tree::iTgbmpb,
                      α   ::Float64,
                      σλ  ::Float64,
                      llc ::Float64,
                      ssλ ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ter ::Bool)

  if isdefined(tree, :d1)
    llc, ssλ = update_triad!(tree, α, σλ, llc, ssλ, δt, srδt)

    llc, ssλ, tree.d1 = _update_gbm!(tree.d1, α, σλ, llc, ssλ, δt, srδt, ter)
    llc, ssλ, tree.d2 = _update_gbm!(tree.d2, α, σλ, llc, ssλ, δt, srδt, ter)
  else
    if !isfix(tree) || ter
      llc, ssλ = update_tip!(tree, α, σλ, μ, llc, ssλ, δt, srδt)
    end
  end

  return llc, ssλ, tree
end




"""
    update_tip!(tree::iTgbmpb
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ssλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` tip proposal.
"""
function update_tip!(tree::iTgbmpb
                     α   ::Float64,
                     σλ  ::Float64,
                     llc ::Float64,
                     ssλ ::Float64,
                     δt  ::Float64,
                     srδt::Float64)

  @inbounds begin

    λc   = lλ(treec)
    l    = lastindex(λc)
    fdtp = fdt(treec)
    λp   = Vector{Float64}(undef, l)

    bm!(λp, λc[1], α, σλ, δt, fdtp, srδt)

    llrbm, llrbd, ssrλ = llr_gbm_b_sep(λp, λc, α, σλ, δt, fdtp, srδt, false)

    acr = llrbd
    llr = acr + llrbm

    if -randexp() < acr
      llc += llr
      ssλ += ssrλ
      unsafe_copyto!(λc, 1, λp, l)
    end 
  end

  return llc, ssλ
end




"""
    update_triad!(λpc ::Vector{Float64},
                  λ1c ::Vector{Float64},
                  λ2c ::Vector{Float64},
                  ep  ::Float64,
                  e1  ::Float64,
                  e2  ::Float64,
                  fdtp::Float64,
                  fdt1::Float64,
                  fdt2::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ssλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` trio proposal.
"""
function update_triad!(λpc ::Vector{Float64},
                       λ1c ::Vector{Float64},
                       λ2c ::Vector{Float64},
                       ep  ::Float64,
                       e1  ::Float64,
                       e2  ::Float64,
                       fdtp::Float64,
                       fdt1::Float64,
                       fdt2::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       ssλ ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  @inbounds begin

    lp   = lastindex(λpc)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λpp  = Vector{Float64}(undef,lp)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    λp   = λpc[1]
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]

    # node proposal
    λn = trioprop(λpr + α*ep, λd1 - α*e1, λd2 - α*e2, ep, e1, e2, σλ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, σλ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, σλ, δt, fdt2, srδt)

    llr, acr, ssrλ = llr_propr(λpp, λ1p, λ2p, λpc, λ1c, λ2c, 
      α, σλ, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < acr
      llc += llr
      ssλ += ssrλ
      unsafe_copyto!(λpc, 1, λpp, lp)
      unsafe_copyto!(λ1c, 1, λ1p, l1)
      unsafe_copyto!(λ2c, 1, λ2p, l2)
    end
  end

  return llc, ssλ
end




"""
    update_triad!(tree::iTgbmpb,
                  α   ::Float64,
                  σλ  ::Float64,
                  llc ::Float64,
                  ssλ ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Make a `gbm` trio proposal.
"""
function update_triad!(tree::iTgbmpb,
                       α   ::Float64,
                       σλ  ::Float64,
                       llc ::Float64,
                       ssλ ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  @inbounds begin

    λpc  = lλ(tree)
    λ1c  = lλ(tree.d1)
    λ2c  = lλ(tree.d2)
    lp   = lastindex(λpc)
    l1   = lastindex(λ1c)
    l2   = lastindex(λ2c)
    λpp  = Vector{Float64}(undef,lp)
    λ1p  = Vector{Float64}(undef,l1)
    λ2p  = Vector{Float64}(undef,l2)
    λp   = λpc[1]
    λ1   = λ1c[l1]
    λ2   = λ2c[l2]
    ep   = e(tree)
    e1   = e(tree.d1)
    e1   = e(tree.d2)
    fdtp = fdt(tree)
    fdt1 = fdt(tree.d1)
    fdt2 = fdt(tree.d2)

    # node proposal
    λn = trioprop(λpr + α*ep, λd1 - α*e1, λd2 - α*e2, ep, e1, e2, σλ)

    # simulate fix tree vector
    bb!(λpp, λp, λn, σλ, δt, fdtp, srδt)
    bb!(λ1p, λn, λ1, σλ, δt, fdt1, srδt)
    bb!(λ2p, λn, λ2, σλ, δt, fdt2, srδt)

    llr, acr, ssrλ = llr_propr(λpp, λ1p, λ2p, λpc, λ1c, λ2c, 
      α, σλ, δt, fdtp, fdt1, fdt2, srδt)

    if -randexp() < acr
      llc += llr
      ssλ += ssrλ
      unsafe_copyto!(λpc, 1, λpp, lp)
      unsafe_copyto!(λ1c, 1, λ1p, l1)
      unsafe_copyto!(λ2c, 1, λ2p, l2)
    end
  end

  return llc, ssλ
end






"""
    triad_lλupdate_root!(treep ::iTgbmpb, 
                         treec ::iTgbmpb,
                         llc   ::Float64,
                         α     ::Float64, 
                         σλ    ::Float64, 
                         δt    ::Float64,
                         srδt  ::Float64,
                         lλmxpr::Float64)

Make a trio of Brownian motion MCMC updates when the root is involved.
"""
function triad_lλupdate_root!(treep ::iTgbmpb, 
                              treec ::iTgbmpb,
                              llc   ::Float64,
                              α     ::Float64, 
                              σλ    ::Float64, 
                              δt    ::Float64,
                              srδt  ::Float64,
                              lλmxpr::Float64)

  # get reference vectors
  λprv_p = lλ(treep)
  λd1v_p = lλ(treep.d1)
  λd2v_p = lλ(treep.d2)
  λprv_c = lλ(treec)
  λd1v_c = lλ(treec.d1)
  λd2v_c = lλ(treec.d2)
  λd1 = λd1v_c[end]
  λd2 = λd2v_c[end]

  # pendant edges
  pepr = e(treec)
  ped1 = e(treec.d1)
  ped2 = e(treec.d2)

  # final dt
  fdtpr = fdt(treec)
  fdtd1 = fdt(treec.d1)
  fdtd2 = fdt(treec.d2)

  # proposal given daughters
  lλp = duoprop(λd1 - α*ped1, λd2 - α*ped2, ped1, ped2, σλ)

  # propose for root
  lλrp = rnorm(lλp - α*pepr, sqrt(pepr)*σλ)

  # make Brownian bridge proposals
  bb!(λprv_p, lλrp, lλp, σλ, δt, fdtpr, srδt)
  bb!(λd1v_p, lλp,  λd1, σλ, δt, fdtd1, srδt)
  bb!(λd2v_p, lλp,  λd2, σλ, δt, fdtd2, srδt)

  ## make acceptance ratio 
  llr, acr, ssrλ = llr_propr(λprv_p, λd1v_p, λd2v_p, λprv_c, λd1v_c, λd2v_c, 
    α, σλ, δt, fdtpr, fdtd1, fdtd2, srδt, !istip(treec.d1), !istip(treec.d2))

  # prior ratio
  if lλrp > lλmxpr
    acr += -Inf
  end

  if -randexp() < acr 
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
end




"""
    llr_propr(λprv_p::Array{Float64,1},
              λd1v_p::Array{Float64,1},
              λd2v_p::Array{Float64,1},
              λprv_c::Array{Float64,1},
              λd1v_c::Array{Float64,1},
              λd2v_c::Array{Float64,1},
              α     ::Float64,
              σλ    ::Float64,
              δt    ::Float64,
              fdtpr ::Float64,
              fdtd1 ::Float64,
              fdtd2 ::Float64,
              srδt  ::Float64,
              itd1  ::Bool,
              itd2  ::Bool)

Return the likelihood and proposal ratio for pure-birth gbm.
"""
function llr_propr(λpp  ::Array{Float64,1},
                   λ1p  ::Array{Float64,1},
                   λ2p  ::Array{Float64,1},
                   λpc  ::Array{Float64,1},
                   λ1c  ::Array{Float64,1},
                   λ2c  ::Array{Float64,1},
                   α    ::Float64,
                   σλ   ::Float64,
                   δt   ::Float64,
                   fdtpr::Float64,
                   fdtd1::Float64,
                   fdtd2::Float64,
                   srδt ::Float64)

  # log likelihood ratios
  llrbmp, llrpbp, ssrλp = llr_gbm_b_sep(λpp, λpc, α, σλ, δt, fdtpr, srδt, true)
  llrbm1, llrpb1, ssrλ1 = llr_gbm_b_sep(λ1p, λ1c, α, σλ, δt, fdtd1, srδt, false)
  llrbm2, llrpb2, ssrλ2 = llr_gbm_b_sep(λ2p, λ2c, α, σλ, δt, fdtd2, srδt, false)

  acr  = llrpbp + llrpb1 + llrpb2
  llr  = llrbmp + llrbm1 + llrbm2 + acr
  ssrλ = ssrλp + ssrλ1 + ssrλ2

  return llr, acr, ssrλ
end




"""
    update_α!(αc     ::Float64,
              σλ     ::Float64,
              L      ::Float64,
              dλ     ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              α_prior::NTuple{2,Float64}) where {T <: iTgbm}

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
    update_σ!(σλc     ::Float64,
              α       ::Float64,
              ssλ     ::Float64,
              n       ::Float64,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64}) 

Gibbs update for `σλ`.
"""
function update_σ!(σλc     ::Float64,
                   α       ::Float64,
                   ssλ     ::Float64,
                   n       ::Float64,
                   llc     ::Float64,
                   prc     ::Float64,
                   σλ_prior::NTuple{2,Float64})

  σλ_p1 = σλ_prior[1]
  σλ_p2 = σλ_prior[2]

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_p1 + 0.5 * n, σλ_p2 + ssλ)

  # update prior
  prc += llrdinvgamma(σλp2, σλc^2, σλ_p1, σλ_p2)

  σλp = sqrt(σλp2)

  # update likelihood
  llc += ssλ*(1.0/σλc^2 - 1.0/σλp2) - n*(log(σλp/σλc))

  return llc, prc, σλp
end



