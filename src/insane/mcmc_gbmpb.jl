#=

Anagenetic GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#





"""
    insane_gbmpb(tree    ::sTpb, 
                 out_file::String;
                 δt      ::Float64  = 1e-2,
                 niter   ::Int64    = 1_000,
                 nthin   ::Int64    = 10,
                 nburn   ::Int64    = 200,
                 tune_int::Int64    = 100,
                 σλi     ::Float64  = 0.5,
                 prints  ::Int64    = 5,
                 pupdp   ::NTuple{2,Float64} = (0.9, 0.1),
                 λa_prior::NTuple{2,Float64} = (0.0,100.0),
                 σλ_prior::NTuple{2,Float64} = (0.05, 0.5))

Run insane for GBM pure-birth.
"""
function insane_gbmpb(tree    ::sTpb, 
                      out_file::String;
                      δt      ::Float64  = 1e-2,
                      niter   ::Int64    = 1_000,
                      nthin   ::Int64    = 10,
                      nburn   ::Int64    = 200,
                      tune_int::Int64    = 100,
                      σλi     ::Float64  = 0.1,
                      αi      ::Float64  = 0.0,
                      prints  ::Int64    = 5,
                      pupdp   ::NTuple{3,Float64} = (0.1,0.1, 0.8),
                      λa_prior::NTuple{2,Float64} = (0.0,100.0),
                      α_prior ::NTuple{2,Float64} = (0.0,10.0),
                      σλ_prior::NTuple{2,Float64} = (0.05, 0.5))

  δt  *= treeheight(tree)
  srδt = sqrt(δt)
  n    = sntn(tree, 0)

  # lλ root node
  lλa = log(λmle_cpb(tree))

  # make Ψ current and proposal parameters
  Ψc = iTgbmpb(tree, δt, srδt, lλa, αi, σλi)
  Ψp = deepcopy(Ψc)

  # make fix Ψ directory
  idv = iBfb[]
  bit = BitArray{1}()
  makeiBf!(Ψc, idv, bit)

  # make parent node directory to `iBfb`
  inodes, terminus = make_inodes(idv)

  # parameter updates (1: α, 2: σ, 3: gbm)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(3) 
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running pure-birth gbm"

  # burn-in phase
  llc, prc, αc, σλc =
    mcmc_burn_gbmpb(Ψp, Ψc, λa_prior, α_prior, σλ_prior, nburn, αi, σλi, 
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
function mcmc_burn_gbmpb(Ψp      ::iTgbmpb,
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

  # starting likelihood and prior
  llc = llik_gbm(Ψc, αc, σλc, δt, srδt)
  prc = logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])      + 
        logdnorm(αc, α_prior[1], α_prior[2]^2)             +
        logdunif(exp(lλ(Ψc)[1]), λa_prior[1], λa_prior[2])

  lλmxpr = log(λa_prior[2])

  ntr = lastindex(inodes)

  pbar = Progress(nburn, prints, "burning mcmc...", 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      if pupi === 1

        # update drift
        llc, prc, αc = update_α!(αc, σλc, Ψc, llc, prc, α_prior)

      elseif pupi === 2

        # update diffusion
        llc, prc, σλc = update_σ!(σλc, αc, Ψc, llc, prc, σλ_prior)

      else 

        tix = ceil(Int64,rand()*ntr)

        ter  = terminus[tix]
        inod = inodes[tix]

        dri = dr(idv[inod])
        ldr = lastindex(dri)

        llc = lλupdate!(Ψp, Ψc, llc, αc, σλc, δt, srδt, lλmxpr, 
                dri, ldr, ter, 0)

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

        # update drift
        llc, prc, αc = update_α!(αc, σλc, Ψc, llc, prc, α_prior)

        # ll0 = llik_gbm(Ψc, αc, σλc, δt, srδt)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, 1, it
        #    return 
        # end

      elseif pupi === 2

        # update diffusion
        llc, prc, σλc = update_σ!(σλc, αc, Ψc, llc, prc, σλ_prior)

        # ll0 = llik_gbm(Ψc, αc, σλc, δt, srδt)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, 1, it
        #    return 
        # end

      else 

        tix = ceil(Int64,rand()*ntr)

        ter  = terminus[tix]
        inod = inodes[tix]

        dri = dr(idv[inod])
        ldr = lastindex(dri)

        llc = lλupdate!(Ψp, Ψc, llc, αc, σλc, δt, srδt, lλmxpr, 
                dri, ldr, ter, 0)

        # ll0 = llik_gbm(Ψc, αc, σλc, δt, srδt)
        #  if !isapprox(ll0, llc, atol = 1e-4)
        #    @show ll0, llc, 1, it
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
    lλupdate!(Ψp      ::iTgbmpb,
              Ψc      ::iTgbmpb,
              llc     ::Float64, 
              prc     ::Float64,
              σλ      ::Float64, 
              δt      ::Float64, 
              srδt    ::Float64, 
              λa_prior::Tuple{Float64,Float64},
              dri     ::BitArray{1},
              ldr     ::Int64,
              ter     ::BitArray{1},
              ix      ::Int64)

Make a λgbm update for a triad.
"""
function lλupdate!(Ψp    ::iTgbmpb,
                   Ψc    ::iTgbmpb,
                   llc   ::Float64, 
                   α     ::Float64, 
                   σλ    ::Float64, 
                   δt    ::Float64, 
                   srδt  ::Float64, 
                   lλmxpr::Float64,
                   dri   ::BitArray{1},
                   ldr   ::Int64,
                   ter   ::BitArray{1},
                   ix    ::Int64)

  if ix === ldr 
    # if root
    if iszero(ldr)
      llc = triad_lλupdate_root!(Ψp, Ψc, llc, α, σλ, δt, srδt, lλmxpr)
    else
      llc = triad_lλupdate_trio!(Ψp, Ψc, llc, α, σλ, δt, srδt, ter)
    end
  elseif ix < ldr
    ix += 1
    if dri[ix]
      llc = 
        lλupdate!(Ψp.d1::iTgbmpb, Ψc.d1::iTgbmpb, 
          llc, α, σλ, δt, srδt, lλmxpr, dri, ldr, ter, ix)
    else
      llc = 
        lλupdate!(Ψp.d2::iTgbmpb, Ψc.d2::iTgbmpb, 
          llc, α, σλ, δt, srδt, lλmxpr, dri, ldr, ter, ix)
    end
  end

  return llc
end




"""
    triad_lλupdate_trio!(treep::iTgbmpb, 
                         treec::iTgbmpb,
                         llc  ::Float64,
                         σλ   ::Float64,
                         δt   ::Float64, 
                         srδt ::Float64,
                         ter  ::BitArray{1})

Make a trio of Brownian motion MCMC updates for an internal node.
"""
function triad_lλupdate_trio!(treep::iTgbmpb, 
                              treec::iTgbmpb,
                              llc  ::Float64,
                              α    ::Float64,
                              σλ   ::Float64,
                              δt   ::Float64, 
                              srδt ::Float64,
                              ter  ::BitArray{1})

  # speciation vectors
  λprv_p = lλ(treep)
  λd1v_p = lλ(treep.d1)
  λd2v_p = lλ(treep.d2)
  λprv_c = lλ(treec)
  λd1v_c = lλ(treec.d1)
  λd2v_c = lλ(treec.d2)
  λpr = λprv_c[1]
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

  if ter[1]
    if ter[2]
      # if both are terminal
      bm!(λprv_p, λpr, α, σλ, δt, fdtpr, srδt)
      lλp = λprv_p[end]
      bm!(λd1v_p, lλp, α, σλ, δt, fdtd1, srδt)
      bm!(λd2v_p, lλp, α, σλ, δt, fdtd2, srδt)

    else
      # if d1 is terminal
      lλp = duoprop(λpr + α*pepr, λd2 - α*ped2, pepr, ped2, σλ)

      # simulate fix tree vector
      bb!(λprv_p, λpr, lλp, fdtpr, σλ, δt, srδt)
      bm!(λd1v_p, lλp, α, σλ, δt, fdtd1, srδt)
      bb!(λd2v_p, lλp, λd2, fdtd2, σλ, δt, srδt)
    end
  elseif ter[2]
    # if d2 is terminal
    # node proposal
    lλp = duoprop(λpr + α*pepr, λd1 - α*ped1, pepr, ped1, σλ)

    # simulate fix tree vector
    bb!(λprv_p, λpr, lλp, fdtpr, σλ, δt, srδt)
    bb!(λd1v_p, lλp, λd1, fdtd1, σλ, δt, srδt)
    bm!(λd2v_p, lλp, α, σλ, δt, fdtd2, srδt)

  else
    # if no terminal branches involved
    # node proposal
    lλp  = trioprop(λpr + α*pepr, λd1 - α*ped1, λd2 - α*ped2, 
             pepr, ped1, ped2, σλ)

    # simulate fix tree vector
    bb!(λprv_p, λpr, lλp, fdtpr, σλ, δt, srδt)
    bb!(λd1v_p, lλp, λd1, fdtd1, σλ, δt, srδt)
    bb!(λd2v_p, lλp, λd2, fdtd2, σλ, δt, srδt)
  end

  # acceptance ratio
  llr, acr = llr_propr(λprv_p, λd1v_p, λd2v_p, λprv_c, λd1v_c, λd2v_c, 
    α, σλ, δt, fdtpr, fdtd1, fdtd2, srδt, !ter[1], !ter[2])

  if -randexp() < acr 
    llc += llr
    copyto!(λprv_c, λprv_p)
    copyto!(λd1v_c, λd1v_p)
    copyto!(λd2v_c, λd2v_p)
  end

  return llc
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
  lλrp = rnorm(lλp + α*pepr, sqrt(pepr)*σλ)

  # make Brownian bridge proposals
  bb!(λprv_p, lλrp, lλp, fdtpr, σλ, δt, srδt)
  bb!(λd1v_p, lλp,  λd1, fdtd1, σλ, δt, srδt)
  bb!(λd2v_p, lλp,  λd2, fdtd2, σλ, δt, srδt)

  ## make acceptance ratio 
  llr, acr = llr_propr(λprv_p, λd1v_p, λd2v_p, λprv_c, λd1v_c, λd2v_c, 
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
function llr_propr(λprv_p::Array{Float64,1},
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

  # log likelihood ratio functions
  llrbm_pr, llrpb_pr = llr_gbm_b_sep(λprv_p, λprv_c, α, σλ, 
                          δt, fdtpr, srδt, true)
  llrbm_d1, llrpb_d1 = llr_gbm_b_sep(λd1v_p, λd1v_c, α, σλ, 
                          δt, fdtd1, srδt, itd1)
  llrbm_d2, llrpb_d2 = llr_gbm_b_sep(λd2v_p, λd2v_c, α, σλ, 
                          δt, fdtd2, srδt, itd2)

  acr = llrpb_pr + llrpb_d1 + llrpb_d2
  llr = llrbm_pr + llrbm_d1 + llrbm_d2 + acr

  return llr, acr
end




"""
    update_σ!(σλc     ::Float64,
              α       ::Float64,
              Ψ       ::iTgbmpb,
              llc     ::Float64,
              prc     ::Float64,
              σλ_prior::NTuple{2,Float64}) where {T <: iTgbm}

Gibbs update for `σλ`.
"""
function update_σ!(σλc     ::Float64,
                   α       ::Float64,
                   Ψ       ::T,
                   llc     ::Float64,
                   prc     ::Float64,
                   σλ_prior::NTuple{2,Float64}) where {T <: iTgbm}

  # standardized sum of squares
  ssλ, n = sss_gbm(Ψ, α, 0.0, 0.0)

  # Gibbs update for σ
  σλp2 = randinvgamma(σλ_prior[1] + 0.5 * n, σλ_prior[2] + ssλ)

  # update prior
  prc += llrdinvgamma(σλp2, σλc^2, σλ_prior[1], σλ_prior[2])

  σλp = sqrt(σλp2)

  # update likelihood
  llc += ssλ*(1.0/σλc^2 - 1.0/σλp^2) - n*(log(σλp/σλc))

  return llc, prc, σλp
end




"""
    update_α!(αc     ::Float64,
              σλ     ::Float64,
              Ψ      ::T,
              llc    ::Float64,
              prc    ::Float64,
              α_prior::NTuple{2,Float64}) where {T <: iTgbm}

Gibbs update for `α`.
"""
function update_α!(αc     ::Float64,
                   σλ     ::Float64,
                   Ψ      ::T,
                   llc    ::Float64,
                   prc    ::Float64,
                   α_prior::NTuple{2,Float64}) where {T <: iTgbm}

  # difference sum and tree length
  dλ, l = treelength_dλ(Ψ, 0.0, 0.0)

  # ratio 
  ν   = α_prior[1]
  τ2  = α_prior[2]^2
  σλ2 = σλ^2
  rs  = σλ2/τ2

  # gibbs update for σ
  αp = rnorm((dλ + rs*ν)/(rs + l), σλ/(rs + l))

  # update prior
  prc += llrdnorm_x(αp, αc, ν, τ2)

  # update likelihood
  llc += 0.5*l/σλ2*(αc^2 - αp^2 + 2.0*dλ*(αp - αc)/l)

  return llc, prc, αp
end



