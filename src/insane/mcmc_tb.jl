#=

Trait driven pure-birth speciation diffusion

Ignacio Quintero Mächler

t(-_-t)

Created 09 02 2026
=#




"""
    insane_tb(tree    ::sT_label,
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
              ασi     ::Float64               = 0.0,
              σσi     ::Float64               = 0.01,
              αλi     ::Float64               = 0.0,
              βλi     ::Float64               = 0.0,
              σλi     ::Float64               = 0.01,
              pupdp   ::NTuple{8,Float64}     = (1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0.1, 0.2),
              δt      ::Float64               = 1e-3,
              prints  ::Int64                 = 5,
              stn     ::Float64               = 0.1,
              tρ      ::Dict{String, Float64} = Dict("" => 1.0))

Run insane for trait driven pure-birth.
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
                   ασi     ::Float64               = 0.0,
                   σσi     ::Float64               = 0.01,
                   αλi     ::Float64               = 0.0,
                   βλi     ::Float64               = 0.0,
                   σλi     ::Float64               = 0.01,
                   pupdp   ::NTuple{8,Float64}     = (1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0.1, 0.2),
                   δt      ::Float64               = 1e-3,
                   prints  ::Int64                 = 5,
                   stn     ::Float64               = 0.1,
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
  llc, prc, ασc, σσc, αλc, βλc, σλc, stn,
    dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns =
      mcmc_burn_tb(Ξ, idf, 
        ασ_prior, σσ_prior, λ0_prior, αλ_prior, βλ_prior, σλ_prior, 
        nburn, ασi, σσi, αλi, βλi, σλi, stn, δt, srδt, inodes, pup, prints)

  # mcmc
  r, treev = 
   mcmc_tb(Ξ, idf, llc, prc, ασc, σσc, αλc, βλc, σλc, stn, 
      dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, 
      ασ_prior, σσ_prior, λ0_prior, αλ_prior, βλ_prior, σλ_prior, 
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

MCMC burn-in chain for trait driven pure-birth.
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
      elseif pupi === 7

        bix = inodes[fIrand(nin) + 1]

        llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ =
          update_internal!(bix, Ξ, idf, ασc, σσc, αλc, βλc, σλc, llc, prc, 
            dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, λ0_prior)

      # forward simulation
      else

        bix = fIrand(el) + 1

        llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L =
          update_fs!(bix, Ξ, idf, ασc, σσc, αλc, βλc, σλc, llc, 
            dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L, δt, srδt)

      end
    end

    ltn += 1
    if ltn === 100
      stn = tune(stn, lac/lup)
      ltn = 0
    end

    next!(pbar)
  end

  return llc, prc, ασc, σσc, αλc, βλc, σλc, stn,
           dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns
end




"""
    mcmc_tb(Ξ       ::Vector{iTxb},
            idf     ::Vector{iBffs},
            llc     ::Float64,
            prc     ::Float64,
            ασc     ::Float64,
            σσc     ::Float64,
            αλc     ::Float64,
            βλc     ::Float64,
            σλc     ::Float64,
            stn     ::Float64,
            dxs     ::Float64,
            dxl     ::Float64,
            ddx     ::Float64,
            ddσ     ::Float64,
            ssσ     ::Float64,
            ddλ     ::Float64,
            ssλ     ::Float64,
            nλ      ::Float64,
            irλ     ::Float64,
            ns      ::Float64,
            ασ_prior::NTuple{2,Float64},
            σσ_prior::NTuple{2,Float64},
            λ0_prior::NTuple{2,Float64},
            αλ_prior::NTuple{2,Float64},
            βλ_prior::NTuple{2,Float64},
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

MCMC chain for for trait driven pure-birth.
"""
function mcmc_tb(Ξ       ::Vector{iTxb},
                 idf     ::Vector{iBffs},
                 llc     ::Float64,
                 prc     ::Float64,
                 ασc     ::Float64,
                 σσc     ::Float64,
                 αλc     ::Float64,
                 βλc     ::Float64,
                 σλc     ::Float64,
                 stn     ::Float64,
                 dxs     ::Float64,
                 dxl     ::Float64,
                 ddx     ::Float64,
                 ddσ     ::Float64,
                 ssσ     ::Float64,
                 ddλ     ::Float64,
                 ssλ     ::Float64,
                 nλ      ::Float64,
                 irλ     ::Float64,
                 ns      ::Float64,
                 ασ_prior::NTuple{2,Float64},
                 σσ_prior::NTuple{2,Float64},
                 λ0_prior::NTuple{2,Float64},
                 αλ_prior::NTuple{2,Float64},
                 βλ_prior::NTuple{2,Float64},
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

  L   = treelength(Ξ)                           # tree length
  nin = lastindex(inodes)                       # number of internal nodes
  el  = lastindex(idf)                          # number of branches

  r = Array{Float64,2}(undef, nlogs, 11)

  treev = iTxb[]  # make Ξ vector
  io = IOBuffer() # buffer 

  open(ofile*".log", "w") do of

    write(of, "iteration\tlikelihood\tprior\tx_root\tsigma2_root\tlambda_root\talpha_sigma\tsigma_sigma\talpha_lambda\tbeta_lambda\tsigma_lambda\n")
    flush(of)

    open(ofile*".txt", "w") do tf

      let llc = llc, prc = prc, ασc = ασc, σσc = σσc, αλc = αλc, βλc = βλc, σλc = σλc, stn = stn, dxs = dxs, dxl = dxl, ddx = ddx, ddσ = ddσ, ssσ = ssσ, ddλ = ddλ, ssλ = ssλ, nλ = nλ, irλ = irλ, ns = ns, L = L, lthin = lthin, lit = lit, sthin = sthin

        pbar = Progress(niter, dt = prints, desc = "running mcmc...", barlen = 20)

        for it in Base.OneTo(niter)

          shuffle!(pup)

          for pupi in pup

            ## parameter updates
            # update drift
            if pupi === 1

              llc, prc, ασc, ssσ = 
                update_α!(ασc, σσc, L, ddσ, llc, prc, ssσ, ασ_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update `σσ` evolutionary rates rate
            elseif pupi === 2

              llc, prc, σσc = update_σ!(σσc, ssσ, nλ, llc, prc, σσ_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update `αλ` speciation rates drift
            elseif pupi === 3

              llc, prc, αλc, ssλ = 
                update_α!(αλc, σλc, L, ddλ - βλc*ddx, llc, prc, ssλ, αλ_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update `βλ` speciation rates trait effect
            elseif pupi === 4

              llc, prc, βλc, ssλ = 
                update_α!(βλc, σλc, dxs, dxl - αλc*ddx, llc, prc, ssλ, βλ_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update `σλ` speciation rates trait effect
            elseif pupi === 5

              llc, prc, σλc = update_σ!(σλc, ssλ, nλ, llc, prc, σλ_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

              # update scale
            elseif pupi === 6

              llc, prc, irλ, acc = 
                update_scale!(Ξ, idf, llc, prc, irλ, ns, stn, λ0_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # update gbm
            elseif pupi === 7

              bix = inodes[fIrand(nin) + 1]

              llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ =
                update_internal!(bix, Ξ, idf, ασc, σσc, αλc, βλc, σλc, llc, prc, 
                  dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, λ0_prior)

              # ll0 = llik_xb(Ξ, idf, ασc, σσc, αλc, βλc, σλc, δt) - Float64(iszero(e(Ξ[1])))*lλ(Ξ[1])[1] + prob_ρ(idf)
              # if !isapprox(ll0, llc, atol = 1e-4)
              #    @show ll0, llc, it, pupi
              #    return
              # end

            # forward simulation
            else

              bix = fIrand(el) + 1

              llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L =
                update_fs!(bix, Ξ, idf, ασc, σσc, αλc, βλc, σλc, llc, 
                  dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L, δt, srδt)

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
              r[lit,1]  = Float64(it)
              r[lit,2]  = llc
              r[lit,3]  = prc
              ξ1 = Ξ[1]
              r[lit,4]  = xv(ξ1)[1]
              r[lit,5]  = exp(lλ(ξ1)[1])
              r[lit,6]  = exp(lσ2(ξ1)[1])
              r[lit,7]  = ασc
              r[lit,8]  = σσc
              r[lit,9]  = αλc
              r[lit,10] = βλc
              r[lit,11] = σλc
              push!(treev, couple(Ξ, idf, 1))
            end
            lthin = zero(Int64)
          end

          # flush parameters
          sthin += 1
          if sthin === nflush
            ξ1 = Ξ[1]
            print(of, Float64(it), '\t', llc, '\t', prc, '\t', 
              xv(ξ1)[1], '\t', exp(lλ(ξ1)[1]), '\t', exp(lσ2(ξ1)[1]), '\t',
              ασc, '\t', σσc, '\t', αλc, '\t', βλc, '\t', σλc, '\n')
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
    update_internal!(bix     ::Int64,
                     Ξ       ::Vector{iTxb},
                     idf     ::Vector{iBffs},
                     ασc     ::Float64, 
                     σσc     ::Float64, 
                     αλc     ::Float64, 
                     βλc     ::Float64, 
                     σλc     ::Float64,
                     llc     ::Float64,
                     prc     ::Float64,
                     dxs    ::Float64,
                     dxl    ::Float64,
                     ddx    ::Float64,
                     ddσ    ::Float64,
                     ssσ    ::Float64,
                     ddλ    ::Float64,
                     ssλ    ::Float64,
                     irλ    ::Float64,
                     δt      ::Float64,
                     srδt    ::Float64,
                     λ0_prior::NTuple{2,Float64})

Make a `gbm` update for an internal branch and its descendants.
"""
function update_internal!(bix     ::Int64,
                          Ξ       ::Vector{iTxb},
                          idf     ::Vector{iBffs},
                          ασ      ::Float64, 
                          σσ      ::Float64, 
                          αλ      ::Float64, 
                          βλ      ::Float64, 
                          σλ      ::Float64,
                          llc     ::Float64,
                          prc     ::Float64,
                          dxs     ::Float64,
                          dxl     ::Float64,
                          ddx     ::Float64,
                          ddσ     ::Float64,
                          ssσ     ::Float64,
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
    llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ =
      _update_crown!(ξi, ξ1, ξ2, ασ, σσ, αλ, βλ, σλ, 
        llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, λ0_prior)
    setλt!(bi, lλ(ξi)[1])
  else
    # if stem
    if root
      llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ = 
        _update_stem!(ξi, ασ, σσ, αλ, βλ, σλ, llc, prc, 
          dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, λ0_prior)
    end

    # updates within the parent branch
    llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ = 
      _update_node!(ξi, xavg(bi), xstd(bi), ασ, σσ, αλ, βλ, σλ, llc, 
        dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, false)

    # get fixed tip
    lξi = fixtip(ξi)

    # make between decoupled trees node update
    llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ = 
      update_triad!(lξi, ξ1, ξ2, ασ, σσ, αλ, βλ, σλ, llc, 
        dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt)

    # set fixed `λ(t)` in branch
    setλt!(bi, lλ(ξ1)[1])
  end

  # # carry on updates in the daughters
  b1 = idf[i1]
  b2 = idf[i2]
  llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ = 
    _update_node!(ξ1, xavg(b1), xstd(b1), ασ, σσ, αλ, βλ, σλ, llc, 
        dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, iszero(d1(b1)))
  llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ = 
    _update_node!(ξ2, xavg(b2), xstd(b2), ασ, σσ, αλ, βλ, σλ, llc, 
        dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ, δt, srδt, iszero(d1(b2)))

  return llc, prc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, irλ
end




"""
    update_fs!(bix ::Int64,
               Ξ   ::Vector{iTxb},
               idf ::Vector{iBffs},
               ασ  ::Float64, 
               σσ  ::Float64, 
               αλ  ::Float64, 
               βλ  ::Float64, 
               σλ  ::Float64,
               llc ::Float64,
               dxs ::Float64,
               dxl ::Float64,
               ddx ::Float64,
               ddσ ::Float64,
               ssσ ::Float64,
               ddλ ::Float64,
               ssλ ::Float64,
               nλ  ::Float64,
               irλ ::Float64,
               ns  ::Float64,
               L   ::Float64,
               δt  ::Float64,
               srδt::Float64)

Forward simulation proposal function for trait driven pure-birth diffusion.
"""
function update_fs!(bix ::Int64,
                    Ξ   ::Vector{iTxb},
                    idf ::Vector{iBffs},
                    ασ  ::Float64, 
                    σσ  ::Float64, 
                    αλ  ::Float64, 
                    βλ  ::Float64, 
                    σλ  ::Float64,
                    llc ::Float64,
                    dxs ::Float64,
                    dxl ::Float64,
                    ddx ::Float64,
                    ddσ ::Float64,
                    ssσ ::Float64,
                    ddλ ::Float64,
                    ssλ ::Float64,
                    nλ  ::Float64,
                    irλ ::Float64,
                    ns  ::Float64,
                    L   ::Float64,
                    δt  ::Float64,
                    srδt::Float64)

  bi  = idf[bix]
  ξc  = Ξ[bix]

  dxsr = dxlr = ddxr = ddσr = ssσr = ddλr = ssλr = irλr = 0.0

  # if terminal node
  if iszero(d1(bi))
    xav = xsd = NaN
    if !isnothing(xavg(bi))
      xav, xsd = xavg(bi), xstd(bi)
    end

    ξp, llr = fsbi_t(bi, xav, xsd, ξc, ασ, σσ, αλ, βλ, σλ, δt, srδt)
  # if internal node
  else
    ξp, llr, dxsr, dxlr, ddxr, ddσr, ssσr, ddλr, ssλr, irλr =
      fsbi_i(bi, ξc, Ξ[d1(bi)], Ξ[d2(bi)], ασ, σσ, αλ, βλ, σλ, δt, srδt)
  end

  # if accepted
  if isfinite(llr)

    ll1, dxs1, dxl1, ddx1, ddσ1, ssσ1, ddλ1, ssλ1, nλ1, irλ1, ns1, L1 = 
      ll_gibbs_xb!(ξp, ασ, σσ, αλ, βλ, σλ, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 

    ll0, dxs0, dxl0, ddx0, ddσ0, ssσ0, ddλ0, ssλ0, nλ0, irλ0, ns0, L0 = 
      ll_gibbs_xb!(ξc, ασ, σσ, αλ, βλ, σλ, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 

    # update quantities
    llc += ll1  - ll0  + llr
    dxs += dxs1 - dxs0 + dxsr
    dxl += dxl1 - dxl0 + dxlr
    ddx += ddx1 - ddx0 + ddxr
    ddσ += ddσ1 - ddσ0 + ddσr
    ssσ += ssσ1 - ssσ0 + ssσr
    ddλ += ddλ1 - ddλ0 + ddλr
    ssλ += ssλ1 - ssλ0 + ssλr
    irλ += irλ1 - irλ0 + irλr
    nλ  += nλ1  - nλ0
    ns  += ns1  - ns0
    L   += L1   - L0

    # set new tree
    Ξ[bix] = ξp
  end

  return llc, dxs, dxl, ddx, ddσ, ssσ, ddλ, ssλ, nλ, irλ, ns, L
end




"""
    fsbi_t(bi  ::iBffs,
           xav ::Float64,
           xsd ::Float64,
           ξc  ::iTxb,
           ασ  ::Float64, 
           σσ  ::Float64, 
           αλ  ::Float64, 
           βλ  ::Float64, 
           σλ  ::Float64,
           δt  ::Float64,
           srδt::Float64)

Forward simulation for branch `bi`.
"""
function fsbi_t(bi  ::iBffs,
                xav ::Float64,
                xsd ::Float64,
                ξc  ::iTxb,
                ασ  ::Float64, 
                σσ  ::Float64, 
                αλ  ::Float64, 
                βλ  ::Float64, 
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  nac = ni(bi)         # current ni
  iρi = (1.0 - ρi(bi)) # inv branch sampling fraction
  lU  = -randexp()     # log-probability

  # current ll
  lc = - log(Float64(nac)) - Float64(nac - 1) * (iszero(iρi) ? 0.0 : log(iρi))

  # forward simulation during branch length
  ξp, nap, nn, llr =
    _sim_tb_t(e(bi), xv(ξc)[1], lσ2(ξc)[1], ασ, σσ,
      lλ(ξc)[1], αλ, βλ, σλ, δt, srδt, lc, lU, iρi, 0, 1, 1_000)

  """
  here: add fixed terminal
  """

  if isfinite(llr)
   # if fix node
    if ifx(bi)

      # propose trait value (if no uncertainty, then xpf = xav)
      xpf = xav
      if xsd > 0.0
        xpf = rnorm(xav, xsd)
      end

      # fix a random tip
      _fixrtip!(ξp, nap)
      lξp  = fixtip(ξp)
      xvp  = xv(lξp)
      xpi  = xvp[1]
      ep   = e(lξp)
      lλp  = lλ(lξp)
      lλpi = lλv[1]
      lξc  = fixtip(ξc)
      xvc  = xv(lξc)

      """
      here, do this well!
      """


      # log-likelihood ratio
      acr = logdnorm(xvpi,   xpf,       intσ2(lσ2(lξp), δt, fdt(lξp))) -
            logdnorm(xvc[1], xvc[end], intσ2(lσ2(lξc), δt, fdt(lξc)))

      lλfp = rnorm(lλpi + αλ*ep + βλ*(xpf - xvpi), sqrt(ep)*σλ)

      cbb!(xvp, xpi, xpf, lλv, lλpi, lλfp,  βλ, σλ, δt, fdti, srδt)



      if lU < acr + llr

        setni!(bi, na)    # set new ni
        return ξp, llr
      end

    else
      if lU < llr
        _fixrtip!(ξp, na)

        setni!(bi, na)    # set new ni
        return ξp, llr
      end
    end
  else
    return ξp, NaN
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
                ασ  ::Float64, 
                σσ  ::Float64, 
                αλ  ::Float64, 
                βλ  ::Float64, 
                σλ  ::Float64,
                δt  ::Float64,
                srδt::Float64)

  # forward simulation during branch length
  t0, na = _sim_tb(e(bi), xv(ξc)[1], lσ2(ξc)[1], ασ, σσ, 
             lλ(ξc)[1], αλ, βλ, σλ, δt, srδt, 1, 1_000)

  if na > 999
    return t0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN
  end

  ntp = na

  lU = -randexp() #log-probability

  # continue simulation only if acr on sum of tip rates is accepted
  acr  = log(Float64(ntp)/Float64(nt(bi)))

  # add sampling fraction
  nac  = ni(bi)                # current ni
  iρi  = (1.0 - ρi(bi))        # branch sampling fraction
  acr -= Float64(nac) * (iszero(iρi) ? 0.0 : log(iρi))

 # fix random tip
  #=
  Look more efficiently selecting which tip based on trait, rate and speciation 
  =#

  xf, lσ2f, lλf = fixrtip!(t0, na, NaN, NaN, NaN)

  llrd, acrd, dxsr, dxlr, ddxr, ddσr, ssσr, ddλr, ssλr, irλr, 
  x1p, x2p, lσ21p, lσ22p, lλ1p, lλ2p =
    _daughters_update!(ξ1, ξ2, xf, lσ2f, lλf, ασ, σσ, αλ, βλ, σλ, δt, srδt)

  acr += acrd

  if lU < acr
    # simulated remaining tips until the present
    t0, na, acr =
      tip_sims!(t0, tf(bi), ασ, σσ, αλ, βλ, σλ, δt, srδt, acr, lU, iρi, na)

    if lU < acr
      na -= 1

      llr = llrd + (na - nac)*(iszero(iρi) ? 0.0 : log(iρi))
      l1  = lastindex(x1p)
      l2  = lastindex(x2p)
      setnt!(bi, ntp)                          # set new nt
      setni!(bi, na)                           # set new ni
      setλt!(bi, lλf)                          # set new lλt
      unsafe_copyto!(xv(ξ1),  1, x1p,   1, l1) # set new daughter 1 x vector
      unsafe_copyto!(xv(ξ2),  1, x2p,   1, l2) # set new daughter 2 x vector
      unsafe_copyto!(lσ2(ξ1), 1, lσ21p, 1, l1) # set new daughter 1 σ vector
      unsafe_copyto!(lσ2(ξ2), 1, lσ22p, 1, l2) # set new daughter 2 σ vector
      unsafe_copyto!(lλ(ξ1),  1, lλ1p,  1, l1) # set new daughter 1 λ vector
      unsafe_copyto!(lλ(ξ2),  1, lλ2p,  1, l2) # set new daughter 2 λ vector

      return t0, llr, dxsr, dxlr, ddxr, ddσr, ssσr, ddλr, ssλr, irλr
    else
      return t0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN
    end
  end

  return t0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN
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
                   ασ  ::Float64, 
                   σσ  ::Float64, 
                   αλ  ::Float64, 
                   βλ  ::Float64, 
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
        x0   = xv(tree)
        lσ20 = lσ2(tree)
        lλ0  = lλ(tree)
        l0   = lastindex(x0)

        # simulate
        stree, na, lr =
          _sim_tb_it(max(δt-fdti, 0.0), t, x0[l0], lσ20[l0], ασ, σσ, 
            lλ0[l0], αλ, βλ, σλ, δt, srδt, lr, lU, iρi, na, 1_000)

        if isnan(lr) || na > 999
          return tree, na, NaN
        end

        sete!(tree, e(tree) + e(stree))

        xs   = xv(stree)
        lσ2s = lσ2(stree)
        lλs  = lλ(stree)

        if lastindex(xs) === 2
          setfdt!(tree, fdt(tree) + fdt(stree))
        else
          setfdt!(tree, fdt(stree))
        end

        pop!(x0)
        popfirst!(xs)
        append!(x0, xs)

        pop!(lσ20)
        popfirst!(lσ2s)
        append!(lσ20, lσ2s)

        pop!(lλ0)
        popfirst!(lλs)
        append!(lλ0, lλs)

        if isdefined(stree, :d1)
          tree.d1 = stree.d1
          tree.d2 = stree.d2
        end
      end
    else
      tree.d1, na, lr = 
        tip_sims!(tree.d1, t, ασ, σσ, αλ, βλ, σλ, δt, srδt, lr, lU, iρi, na)
      tree.d2, na, lr = 
        tip_sims!(tree.d2, t, ασ, σσ, αλ, βλ, σλ, δt, srδt, lr, lU, iρi, na)
    end

    return tree, na, lr
  end

  return tree, na, NaN
end





