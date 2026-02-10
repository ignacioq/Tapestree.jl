#=

GBM pure-birth MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 14 09 2020
=#




"""
    insane_tb(tree    ::sT_label;
                λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                α_prior ::NTuple{2,Float64}     = (0.0, 1.0),
                σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                niter   ::Int64                 = 1_000,
                nthin   ::Int64                 = 10,
                nburn   ::Int64                 = 200,
                nflush  ::Int64                 = nthin,
                ofile   ::String                = string(homedir(), "/ib"),
                αi      ::Float64               = 0.0,
                σλi     ::Float64               = 0.1,
                pupdp   ::NTuple{5,Float64}     = (1e-3, 1e-3, 1e-3, 0.2, 0.2),
                δt      ::Float64               = 1e-3,
                prints  ::Int64                 = 5,
                stn     ::Float64               = 0.5,
                tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for pure `b`.
"""
function insane_tb(tree    ::sT_label,
                   xa      ::Dict{String, Float64};
                   xs      ::Dict{String, Float64} = Dict{String,Float64}(),
                   ασ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                   σσ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                   λ0_prior::NTuple{2,Float64}     = (0.05, 148.41),
                   αλ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                   βλ_prior::NTuple{2,Float64}     = (0.0, 1.0),
                   σλ_prior::NTuple{2,Float64}     = (0.05, 0.05),
                   niter   ::Int64                 = 1_000,
                   nthin   ::Int64                 = 10,
                   nburn   ::Int64                 = 200,
                   nflush  ::Int64                 = nthin,
                   ofile   ::String                = string(homedir(), "/tb"),
                   αxi     ::Float64               = 0.0,
                   ασi     ::Float64               = 0.0,
                   σσi     ::Float64               = 0.01,
                   αλi     ::Float64               = 0.0,
                   βλi     ::Float64               = 0.0,
                   σλi     ::Float64               = 0.01,
                   pupdp   ::NTuple{5,Float64}     = (1e-3, 1e-3, 1e-3, 0.2, 0.2),
                   δt      ::Float64               = 1e-3,
                   prints  ::Int64                 = 5,
                   stn     ::Float64               = 0.1,
                   mxthf    ::Float64              = Inf,
                   tρ      ::Dict{String, Float64} = Dict("" => 1.0))

  n    = ntips(tree)
  th   = treeheight(tree)
  δt  *= max(0.1, round(th, RoundDown, digits = 2))
  srδt = sqrt(δt)

  # turn to logarithmic terms
  λ0_prior = (log(λ0_prior[1]), 2.0*log(λ0_prior[2]))

  # set tips sampling fraction
  if isone(length(tρ))
    tl = tiplabels(tree)
    tρu = tρ[""]
    tρ = Dict(tl[i] => tρu for i in 1:n)
  end

  # make fix tree directory
  idf, xr, σxi = make_idf(tree, tρ, xa, xs, Inf)

  # make a decoupled tree
  Ξ = make_Ξ(idf, xr, σxi, σσi, λmle_cb(tree), σλi, δt, srδt, iTxb)

  # get vector of internal branches
  inodes = [i for i in Base.OneTo(lastindex(idf))  if d1(idf[i]) > 0]

  # parameter updates (1: α, 2: σ, 3: scale, 4: gbm, 5: fs)
  spup = sum(pupdp)
  pup  = Int64[]
  for i in Base.OneTo(lastindex(pupdp))
    append!(pup, fill(i, ceil(Int64, Float64(2*n - 1) * pupdp[i]/spup)))
  end

  @info "running trait driven pure-birth diffusion"

  # burn-in phase
  Ξ, idf, llc, prc, αc, σλc, ns, stn =
    mcmc_burn_tb(Ξ, idf, λ0_prior, α_prior, σλ_prior, nburn, αi, σλi, stn,
      δt, srδt, inodes, pup, prints)

  # mcmc
  r, treev = 
    mcmc_tb(Ξ, idf, llc, prc, αc, σλc, ns, stn, λ0_prior, α_prior, σλ_prior,
      δt, srδt, inodes, pup, niter, nthin, nflush, ofile, prints)

  return r, treev
end



"""
    mcmc_burn_tb(Ξ       ::Vector{iTxb},
                 idf     ::Vector{iBffs},
                 ασ_prior::NTuple{2,Float64},
                 σσ_prior::NTuple{2,Float64},
                 λ0_prior::NTuple{2,Float64},
                 αλ_prior::NTuple{2,Float64},
                 βλ_prior::NTuple{2,Float64},
                 σλ_prior::NTuple{2,Float64},
                 nburn   ::Int64,
                 ασc     ::Float64,
                 σσc     ::Float64,
                 αλc     ::Float64,
                 βλc     ::Float64,
                 σλc     ::Float64,
                 stn     ::Float64,
                 δt      ::Float64,
                 srδt    ::Float64,
                 inodes  ::Array{Int64,1},
                 pup     ::Array{Int64,1},
                 prints  ::Int64)

MCMC burn-in chain for trait driven pure birth.
"""
function mcmc_burn_tb(Ξ       ::Vector{iTxb},
                      idf     ::Vector{iBffs},
                      ασ_prior::NTuple{2,Float64},
                      σσ_prior::NTuple{2,Float64},
                      λ0_prior::NTuple{2,Float64},
                      αλ_prior::NTuple{2,Float64},
                      βλ_prior::NTuple{2,Float64},
                      σλ_prior::NTuple{2,Float64},
                      nburn   ::Int64,
                      ασc     ::Float64,
                      σσc     ::Float64,
                      αλc     ::Float64,
                      βλc     ::Float64,
                      σλc     ::Float64,
                      stn     ::Float64,
                      δt      ::Float64,
                      srδt    ::Float64,
                      inodes  ::Array{Int64,1},
                      pup     ::Array{Int64,1},
                      prints  ::Int64)

  nsi = Float64(iszero(e(Ξ[1])))

  # starting likelihood and prior
  lλ0 = lλ(Ξ[1])[1]
  llc = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - nsi*lλ0 + prob_ρ(idf)
  prc = logdnorm(ασc,       ασ_prior[1], ασ_prior[2]^2) + 
        logdinvgamma(σσc^2, σσ_prior[1], σσ_prior[2])   + 
        logdnorm(lλ0,       λ0_prior[1], λ0_prior[2])   +
        logdnorm(αλc,       αλ_prior[1], αλ_prior[2]^2) +
        logdnorm(βλc,       βλ_prior[1], βλ_prior[2]^2) +
        logdinvgamma(σλc^2, σλ_prior[1], σλ_prior[2])

  L   = treelength(Ξ)                           # tree length
  nin = lastindex(inodes)                       # number of internal nodes
  el  = lastindex(idf)                          # number of branches
  ns  = sum(x -> Float64(d2(x) > 0), idf) - nsi # number of speciation events in likelihood

  # delta change, sum squares, path length and integrated rate
  dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ = 
    _gibbs_quanta(Ξ, ασc, αλc, βλc)

  # for scale tuning
  ltn = 0
  lup = lac = 0.0

  pbar = Progress(nburn, dt = prints, desc = "burn-in mcmc...", barlen = 20)

  for it in Base.OneTo(nburn)

    shuffle!(pup)

    for pupi in pup

      ## parameter updates
      # update `ασ` evolutionary rates drift
      if pupi === 1

        llc, prc, ασc, ssσ = 
          update_α!(ασc, σσc, L, ddσ, llc, prc, ssσ, ασ_prior)

      # update `σσ` evolutionary rates rate
      elseif pupi === 2

        llc, prc, σσc = update_σ!(σσc, ssσ, nλ, llc, prc, σσ_prior)

      # update `αλ` speciation rates drift
      elseif pupi === 3

        llc, prc, αλc, ssλ = 
          update_α!(αλc, σλc, L, ddλ - βλc*ddx, llc, prc, ssλ, αλ_prior)

      # update `βλ` speciation rates trait effect
      elseif pupi === 4

        llc, prc, βλc, ssλ = 
          update_α!(βλc, σλc, dxs, dxl - αλc*ddx, llc, prc, ssλ, βλ_prior)

      # update `σλ` speciation rates trait effect
      elseif pupi === 5

        llc, prc, σλc = update_σ!(σλc, ssλ, nλ, llc, prc, σλ_prior)

      # update scale
      elseif pupi === 6

        llc, prc, irλ, acc = 
          update_scale!(Ξ, idf, llc, prc, irλ, ns, stn, λ0_prior)

        lac += acc
        lup += 1.0

      # update gbm
      elseif pupi === 4

"""
here
"""

        nix = ceil(Int64,rand()*nin)
        bix = inodes[nix]

        llc, prc, ddλ, ssλ, irλ =
          update_gbm!(bix, Ξ, idf, αc, σλc, llc, prc, ddλ, ssλ, irλ, 
            δt, srδt, λ0_prior)

      # forward simulation
      else

        bix = ceil(Int64,rand()*el)

        llc, ddλ, ssλ, nλ, irλ, ns, L =
          update_fs!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, nλ, irλ, ns, L, 
            δt, srδt)

      end
    end

    ltn += 1
    if ltn === 100
      stn = tune(stn, lac/lup)
      ltn = 0
    end

    next!(pbar)
  end

  return Ξ, idf, llc, prc, αc, σλc, ns, stn
end




"""
    mcmc_tb( Ξ       ::Vector{iTxb},
               idf     ::Vector{iBffs},
               llc     ::Float64,
               prc     ::Float64,
               αc      ::Float64,
               σλc     ::Float64,
               ns      ::Float64,
               λ0_prior::NTuple{2,Float64},
               α_prior ::NTuple{2,Float64},
               σλ_prior::NTuple{2,Float64},
               δt      ::Float64,
               srδt    ::Float64,
               inodes  ::Array{Int64,1},
               pup     ::Vector{Int64},
               niter   ::Int64,
               nthin   ::Int64,
               nflush  ::Int64,
               ofile   ::String,
               prints  ::Int64)

MCMC chain for for trait driven pure birth.
"""
function mcmc_tb( Ξ       ::Vector{iTxb},
                    idf     ::Vector{iBffs},
                    llc     ::Float64,
                    prc     ::Float64,
                    αc      ::Float64,
                    σλc     ::Float64,
                    ns      ::Float64,
                    stn     ::Float64,
                    λ0_prior::NTuple{2,Float64},
                    α_prior ::NTuple{2,Float64},
                    σλ_prior::NTuple{2,Float64},
                    δt      ::Float64,
                    srδt    ::Float64,
                    inodes  ::Array{Int64,1},
                    pup     ::Vector{Int64},
                    niter   ::Int64,
                    nthin   ::Int64,
                    nflush  ::Int64,
                    ofile   ::String,
                    prints  ::Int64)

  # logging
  nlogs = fld(niter,nthin)
  lthin = lit = sthin = zero(Int64)

  r = Array{Float64,2}(undef, nlogs, 6)

  L   = treelength(Ξ)      # tree length
  nin = lastindex(inodes)  # number of internal nodes
  el  = lastindex(idf)     # number of branches

  # delta change, sum squares, path length and integrated rate
  ddλ, ssλ, nλ, irλ = 
    _ss_ir_dd(Ξ, lλ, αc)

  treev = iTxb[]  # make Ξ vector
  io = IOBuffer() # buffer 

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tlambda_root\talpha\tsigma_lambda\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, αc = αc, σλc = σλc, ns = ns, nλ = nλ, ssλ = ssλ, ddλ = ddλ, irλ = irλ, L = L, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update drift
            if pupi === 1

              llc, prc, αc = update_α!(αc, σλc, L, ddλ, llc, prc, α_prior)

              # update ssλ with new drift `α`
              ssλ = _ss(Ξ, lλ, αc)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update diffusion rate
            elseif pupi === 2

              llc, prc, σλc = update_σ!(σλc, ssλ, nλ, llc, prc, σλ_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update scale
            elseif pupi === 3

              llc, prc, irλ, acc = 
                update_scale!(Ξ, idf, llc, prc, irλ, ns, stn, λ0_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update gbm
            elseif pupi === 4

              nix = ceil(Int64,rand()*nin)
              bix = inodes[nix]

              llc, prc, ddλ, ssλ, irλ =
                update_gbm!(bix, Ξ, idf, αc, σλc, llc, prc, ddλ, ssλ, irλ, 
                  δt, srδt, λ0_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update by forward simulation
            else

              bix = ceil(Int64,rand()*el)

              llc, ddλ, ssλ, nλ, irλ, ns, L =
                update_fs!(bix, Ξ, idf, αc, σλc, llc, ddλ, ssλ, nλ, irλ, ns, L, 
                  δt, srδt)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
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
              r[lit,4] = exp(lλ(Ξ[1])[1])
              r[lit,5] = αc
              r[lit,6] = σλc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
                  exp(lλ(Ξ[1])[1]),'\t', αc, '\t', σλc, '\n')
            flush(of)
            ibuffer(io, couple(Ξ, idf, 1))
            write(io, '\n')
            write(tf, take!(io))
            flush(tf)
            sthin = zero(Int64)
          end

          next!(pbar)
        end

        return r, treev
      end
    end
  end
end




"""
    update_α!(αc     ::Float64,
              σ      ::Float64,
              L      ::Float64,
              dd     ::Float64,
              ll     ::Float64,
              pr     ::Float64,
              ss     ::Float64,
              α_prior::NTuple{2,Float64})

Gibbs update for Normal conjugacy `α`.
"""
function update_α!(αc     ::Float64,
                   σ      ::Float64,
                   L      ::Float64,
                   dd     ::Float64,
                   ll     ::Float64,
                   pr     ::Float64,
                   ss     ::Float64,
                   α_prior::NTuple{2,Float64})

  # ratio
  ν  = α_prior[1]
  τ2 = α_prior[2]^2
  σ2 = σ^2
  rs = σ2/τ2

  # gibbs update for σ
  αp = rnorm((dd + rs*ν)/(rs + L), sqrt(σ2/(rs + L)))

  # update prior
  pr += llrdnorm_x(αp, αc, ν, τ2)

  # update likelihood
  ll += 0.5*L/σ2*(αc^2 - αp^2 + 2.0*dd*(αp - αc)/L)

  # update residual ss
  ss += 0.5*L*(αp^2 - αc^2) - (αp - αc)*dd

  return ll, pr, αp, ss
end





"""
    update_σ!(σc     ::Float64,
              ss     ::Float64,
              n      ::Float64,
              llc    ::Float64,
              prc    ::Float64,
              σ_prior::NTuple{2,Float64})

Gibbs update for variance `σ`.
"""
function update_σ!(σc     ::Float64,
                   ss     ::Float64,
                   n      ::Float64,
                   llc    ::Float64,
                   prc    ::Float64,
                   σ_prior::NTuple{2,Float64})

  σ_p1, σ_p2 = σ_prior

  # Gibbs update for σ
  σp2 = rand(InverseGamma(σ_p1 + 0.5 * n, σ_p2 + ss))

  # update prior
  prc += llrdinvgamma(σp2, σc^2, σ_p1, σ_p2)

  σp = sqrt(σp2)

  # update likelihood
  llc += ss*(1.0/σc^2 - 1.0/σp2) - n*(log(σp/σc))

  return llc, prc, σp
end




# """
#     update_scale!(Ξ       ::Vector{T},
#                   idf     ::Vector{iBffs},
#                   llc     ::Float64,
#                   prc     ::Float64,
#                   ir      ::Float64,
#                   ns      ::Float64,
#                   stn     ::Float64,
#                   λ0_prior::NTuple{2,Float64}) where {T <: iTree}

# Update scale for speciation.
# """
# function update_scale!(Ξ       ::Vector{T},
#                        idf     ::Vector{iBffs},
#                        llc     ::Float64,
#                        prc     ::Float64,
#                        ir      ::Float64,
#                        ns      ::Float64,
#                        stn     ::Float64,
#                        λ0_prior::NTuple{2,Float64}) where {T <: iTree}

#   # sample log(scaling factor)
#   s = randn()*stn

#   # likelihood ratio
#   iri = (1.0 - exp(s)) * ir
#   llr = ns * s + iri

#   lλ0 = lλ(Ξ[1])[1]

#   # prior ratio
#   prr = llrdnorm_x(lλ0 + s, lλ0, λ0_prior[1], λ0_prior[2]) 

#   acc = 0.0

#   if -randexp() < llr + prr
#     acc += 1.0
#     llc += llr
#     prc += prr
#     ir  -= iri
#     scale_rate!(Ξ, lλ, s)
#     scale_rate!(idf, s)
#   end

#   return llc, prc, ir, acc
# end




"""
    update_gbm!(bix     ::Int64,
                Ξ       ::Vector{iTxb},
                idf     ::Vector{iBffs},
                α       ::Float64,
                σλ      ::Float64,
                llc     ::Float64,
                prc     ::Float64,
                ddλ     ::Float64,
                ssλ     ::Float64,
                irλ     ::Float64,
                δt      ::Float64,
                srδt    ::Float64,
                λ0_prior::NTuple{2,Float64})

Make a `gbm` update for an internal branch and its descendants.
"""
function update_gbm!(bix     ::Int64,
                     Ξ       ::Vector{iTxb},
                     idf     ::Vector{iBffs},
                     α       ::Float64,
                     σλ      ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     ddλ     ::Float64,
                     ssλ     ::Float64,
                     irλ     ::Float64,
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
    llc, prc, ddλ, ssλ, irλ =
      _crown_update!(ξi, ξ1, ξ2, α, σλ, llc, prc, ddλ, ssλ, irλ, δt, srδt, λ0_prior)
    setλt!(bi, lλ(ξi)[1])
  else
    # if stem
    if root
      llc, prc, ddλ, ssλ, irλ = 
        _stem_update!(ξi, α, σλ, llc, prc, ddλ, ssλ, irλ, δt, srδt, λ0_prior)
    end

    # updates within the parent branch
    llc, ddλ, ssλ, irλ = 
      _update_gbm!(ξi, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt, false)

    # get fixed tip
    lξi = fixtip(ξi)

    # make between decoupled trees node update
    llc, ddλ, ssλ, irλ = 
      update_triad_b!(lλ(lξi), lλ(ξ1), lλ(ξ2), e(lξi), e(ξ1), e(ξ2),
        fdt(lξi), fdt(ξ1), fdt(ξ2), α, σλ, llc, ddλ, ssλ, irλ, δt, srδt)

    # set fixed `λ(t)` in branch
    setλt!(bi, lλ(lξi)[end])
  end

  # # carry on updates in the daughters
  llc, ddλ, ssλ, irλ = 
    _update_gbm!(ξ1, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt, iszero(d1(idf[i1])))
  llc, ddλ, ssλ, irλ = 
    _update_gbm!(ξ2, α, σλ, llc, ddλ, ssλ, irλ, δt, srδt, iszero(d1(idf[i2])))

  return llc, prc, ddλ, ssλ, irλ
end




"""
    update_fs!(bix  ::Int64,
               Ξ    ::Vector{iTxb},
               idf  ::Vector{iBffs},
               α    ::Float64,
               σλ   ::Float64,
               llc  ::Float64,
               ddλ  ::Float64,
               ssλ  ::Float64,
               nλ   ::Float64,
               irλ  ::Float64,
               ns   ::Float64,
               L    ::Float64,
               δt   ::Float64,
               srδt ::Float64)

Forward simulation proposal function for pure birth diffusion.
"""
function update_fs!(bix  ::Int64,
                    Ξ    ::Vector{iTxb},
                    idf  ::Vector{iBffs},
                    α    ::Float64,
                    σλ   ::Float64,
                    llc  ::Float64,
                    ddλ  ::Float64,
                    ssλ  ::Float64,
                    nλ   ::Float64,
                    irλ  ::Float64,
                    ns   ::Float64,
                    L    ::Float64,
                    δt   ::Float64,
                    srδt ::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]

  # if terminal
  if iszero(d1(bi))
    ξp, llr = fsbi_t(bi, ξc, α, σλ, δt, srδt)
    ddrλ = ssrλ = irrλ = 0.0
  # if internal
  else
    ξp, llr, ddrλ, ssrλ, irrλ =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], α, σλ, δt, srδt)
  end

  # if accepted
  if isfinite(llr)
    ll1, ddλ1, ssλ1, nλ1, irλ1, ns1 = llik_gbm_ssλ(ξp, α, σλ, δt, srδt, 0.0)
    ll0, ddλ0, ssλ0, nλ0, irλ0, ns0 = llik_gbm_ssλ(ξc, α, σλ, δt, srδt, 0.0)

    # update quantities
    llc += ll1  - ll0 + llr
    ddλ += ddλ1 - ddλ0 + ddrλ
    ssλ += ssλ1 - ssλ0 + ssrλ
    nλ  += nλ1  - nλ0
    irλ += irλ1 - irλ0 + irrλ
    ns  += ns1  - ns0
    L   += treelength(ξp) - treelength(ξc)

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, ddλ, ssλ, nλ, irλ, ns, L
end




"""
    fsbi_t(bi  ::iBffs,
           ξc  ::iTxb,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                ξc  ::iTxb,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  t0, nap, nn, llr =
    _sim_tb_t(e(bi), lλ(ξc)[1], α, σλ, δt, srδt, lc, lU, iρi, 0, 1, 500)

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
           ξ1  ::iTxb,
           ξ2  ::iTxb,
           λ0  ::Float64,
           α   ::Float64,
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`
"""
function fsbi_i(bi  ::iBffs,
                ξc  ::iTxb,
                ξ1  ::iTxb,
                ξ2  ::iTxb,
                α   ::Float64,
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  t0, na = _sim_tb(e(bi), lλ(ξc)[1], α, σλ, δt, srδt, 1, 1_000)

  if na >= 1_000
    return t0, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(ntp/nt(bi))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

 # fix random tip
  λf = fixrtip!(t0, na, NaN)

  llrd, acrd, drλ, ssrλ, irrλ, λ1p, λ2p =
    _daughters_update!(ξ1, ξ2, λf, α, σλ, δt, srδt)

  acr += acrd

  if lU < acr

    # simulated remaining tips until the present
    t0, na, acr =
      tip_sims!(t0, tf(bi), α, σλ, δt, srδt, acr, lU, iρi, na)

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
      l1  = lastindex(λ1p)
      l2  = lastindex(λ2p)
      setnt!(bi, ntp)                    # set new nt
      setni!(bi, na)                     # set new ni
      setλt!(bi, λf)                     # set new λt
      unsafe_copyto!(lλ(ξ1), 1, λ1p, 1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2), 1, λ2p, 1, l2) # set new daughter 2 λ vector

      return t0, llr, drλ, ssrλ, irrλ
    else
      return t0, NaN, NaN, NaN, NaN
    end
  end

  return t0, NaN, NaN, NaN, NaN
end




"""
    tip_sims!(tree::iTxb,
              t   ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64)

Continue simulation until time `t` for unfixed tips in `tree`.
"""
function tip_sims!(tree::iTxb,
                   t   ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64)

 if lU < lr && na < 1_000

    if istip(tree)
      if !isfix(tree)

        fdti = fdt(tree)
        lλ0  = lλ(tree)

        # simulate
        stree, na, lr =
          _sim_tb_it(max(δt-fdti, 0.0), t, lλ0[end], α, σλ, δt, srδt,
            lr, lU, iρi, na, 1_000)

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
      tree.d1, na, lr = tip_sims!(tree.d1, t, α, σλ, δt, srδt, lr, lU, iρi, na)
      tree.d2, na, lr = tip_sims!(tree.d2, t, α, σλ, δt, srδt, lr, lU, iρi, na)
    end

    return tree, na, lr
  end

  return tree, na, NaN
end





