#=

Burning phase for tribe.

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#




"""
    burn_compete(...)

Burning & adaptive phase for MCMC.
"""
function burn_tribe(total_llf     ::Function,
                    ω1upd_llr     ::Function,
                    ω0upd_llr     ::Function,
                    λ1upd_llr     ::Function,
                    λ0upd_llr     ::Function,
                    Xupd_llr      ::Function,
                    Rupd_llr      ::Function,
                    ωxupd_llr     ::Function,
                    σ²upd_llr     ::Function,
                    bgiid         ::Function,
                    bgiid_br      ::Function,
                    mhr_upd_Xbr   ::Function,
                    mhr_upd_Xtrio ::Function,
                    mhr_upd_Ybr   ::Function,
                    mhr_upd_Ytrio ::Function,
                    mhr_upd_Ystem ::Function,
                    mhr_upd_XYbr  ::Function,
                    mhr_upd_XYtrio::Function,
                    Xc      ::Array{Float64,2},
                    Yc      ::Array{Int64,3},
                    δXc     ::Array{Float64,3},
                    δYc     ::Array{Float64,3},
                    LApc    ::Array{Float64,2},
                    LAnc    ::Array{Float64,2},
                    LDc     ::Array{Float64,3},
                    σ²c     ::Float64,
                    ωxc     ::Float64,
                    ω1c     ::Float64,
                    ω0c     ::Float64,
                    λ1c     ::Float64,
                    λ0c     ::Float64,
                    Xnc1    ::Array{Int64,1},
                    Xnc2    ::Array{Int64,1},
                    brl     ::Array{Float64,1},
                    wcol    ::Array{Array{Int64,1},1},
                    bridx_a ::Array{Array{UnitRange{Int64},1},1},
                    brs     ::Array{Int64,3},
                    stemevc ::Array{Array{Float64,1},1},
                    trios   ::Array{Array{Int64,1},1},
                    wXp     ::Array{Int64,1},
                    λprior  ::Float64,
                    ωxprior ::NTuple{2,Float64},
                    ω1prior ::NTuple{2,Float64},
                    ω0prior ::NTuple{2,Float64},
                    σ²prior ::Float64,
                    np      ::Int64,
                    parvec  ::Array{Int64,1},
                    nburn   ::Int64,
                    obj_ar  ::Float64 = 0.234,
                    tune_int::Int64   = 100)

  m, ntip, narea = size(Yc)

  nedge = size(brs, 1) 

  # likelihood and prior
  if ωxc >= 0.0
    llc = total_llf(Xc, Yc, LApc, LDc, ωxc, ω1c, ω0c, λ1c, λ0c,
                    stemevc, brs, σ²c)
  else 
    llc = total_llf(Xc, Yc, LAnc, LDc, ωxc, ω1c, ω0c, λ1c, λ0c,
                    stemevc, brs, σ²c)
  end

  prc = allλpr(   λ1c, λ0c, λprior)           +
        logdexp(  σ²c, σ²prior)               +
        logdnorm(ωxc, ωxprior[1], ωxprior[2]) +
        logdnorm(ω1c, ω1prior[1], ω1prior[2]) +
        logdnorm(ω0c, ω0prior[1], ω0prior[2])

  # make scaling function
  scalef = makescalef(obj_ar)

  # rest of tuning parameters
  ptn = fill(0.1, np) 

  # initialize acceptance log
  ltn = zeros(Int64, np)
  lup = zeros(Int64, np)
  lac = zeros(Int64, np)

  # row i proposals for X
  xpi  = fill(NaN, ntip)             # proposal x slice
  δxi  = fill(NaN, ntip, ntip)       # lineage pairwise differences
  lapi = fill(NaN, ntip)             # lineage average for ωx > 0
  lani = fill(NaN, ntip)             # lineage average for ωx < 0
  ldi  = fill(NaN, ntip, narea)      # lineage difference

  rj = wcol[1][2]

  # progress bar
  p = Progress(nburn, dt=5, desc="burn...", barlen=20, color=:green)

  # print number of parameters
  printstyled(
    "\nωx updates per iter = ", lastindex(filter(x -> x == 2,parvec)),
    "\nω1 updates per iter = ", lastindex(filter(x -> x == 3,parvec)),
    "\nω0 updates per iter = ", lastindex(filter(x -> x == 4,parvec)),
    "\nσ² updates per iter = ", lastindex(filter(x -> x == 1,parvec)),
    "\nλ1 updates per iter = ", lastindex(filter(x -> x == 5,parvec)),
    "\nλ0 updates per iter = ", lastindex(filter(x -> x == 6,parvec)),
    "\nθ  updates per iter = ", length(parvec), "\n", color = :green)

  Xcidx = CartesianIndices(Xc)

  #start burnin
  for it = Base.OneTo(nburn)

    # Update vector
    shuffle!(parvec)

    for up = parvec

      # update X
      if up > 6

        # X updates
        @inbounds begin

          upx = wXp[up - 6]

          # if root
          if upx == 1

            # allocate
            @simd for i = Base.OneTo(ntip)
              xpi[i] =  Xc[1,i]
            end

            # update xi
            addupt!(xpi, ptn, 1, up)

            xpi[rj] = xpi[1]::Float64

            llr = Rupd_llr(xpi, Xc, σ²c)::Float64

            if -Random.randexp() < llr
              llc     += llr::Float64
              Xc[1,:]  = xpi::Array{Float64,1}
              lac[up] += 1   # log acceptance
            end

          else

            xi, xj = Xcidx[upx].I

            # allocate
            for j = Base.OneTo(ntip)
              xpi[j]  = Xc[xi,j]
              lapi[j] = LApc[xi,j]
              lani[j] = LAnc[xi,j]
              @simd for k = Base.OneTo(narea)
                ldi[j,k] = LDc[xi,j,k]
              end
              @simd for i = Base.OneTo(ntip)
                δxi[i,j] = δXc[i,j,xi]
              end
            end

            # update xi
            addupt!(xpi, ptn, xj, up)

            if in(upx, Xnc1)        # if an internal node
              xpi[Xcidx[Xnc2[findfirst(isequal(upx),Xnc1)]][2]] = xpi[xj]
            end

            # calculate new averages
            Xupd_linavg!(δxi, lapi, lani, ldi, wcol, 
                         xpi, xi, xj, Yc, δYc, narea)

            if ωxc >= 0.0
              llr = Xupd_llr(xi, xpi, Xc, lapi, ldi, LApc, LDc, Yc, 
                             ωxc, ω1c, ω0c, λ1c, λ0c, σ²c)::Float64
            else
              llr = Xupd_llr(xi, xpi, Xc, lani, ldi, LAnc, LDc, Yc, 
                             ωxc, ω1c, ω0c, λ1c, λ0c, σ²c)::Float64
            end

            if -Random.randexp() < llr
              llc        += llr::Float64
              Xc[xi,:]    = xpi::Array{Float64,1}
              δXc[:,:,xi] = δxi::Array{Float64,2}
              LApc[xi,:]  = lapi::Array{Float64,1}
              LAnc[xi,:]  = lani::Array{Float64,1}
              LDc[xi,:,:] = ldi::Array{Float64,2}
              lac[up]    += 1   # log acceptance
            end
          end
        end

      # update λ1
      elseif up == 5

        # update λ
        λ1p = mulupt(λ1c, ptn[5])::Float64

        llr = λ1upd_llr(Yc, λ1c, λ1p, λ0c, ω1c, LDc, stemevc, brs)

        prr = llrdexp_x(λ1p, λ1c, λprior)

        if -Random.randexp() < (llr + prr + log(λ1p/λ1c))
          llc    += llr::Float64
          prc    += prr::Float64
          λ1c     = λ1p::Float64
          lac[5] += 1
        end

        # update internal node
        if rand() < 0.4
          λϕ1, λϕ0 = λϕprop()
          llc = mhr_upd_Ytrio(rand(trios), Xc, Yc, 
                      λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                      LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
        end

      # update λ0
      elseif up == 6

        # update λ
        λ0p = mulupt(λ0c, ptn[6])::Float64

        llr = λ0upd_llr(Yc, λ1c, λ0c, λ0p, ω0c, LDc, stemevc, brs)::Float64

        prr = llrdexp_x(λ0p, λ0c, λprior)

        if -Random.randexp() < (llr + prr + log(λ0p/λ0c))
          llc     += llr::Float64
          prc     += prr::Float64
          λ0c      = λ0p::Float64
          lac[6] += 1
        end

        # update internal node
        if rand() < 0.4
          λϕ1, λϕ0 = λϕprop()
          llc = mhr_upd_Ytrio(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
        end

      # if σ² is updated
      elseif up == 1

        σ²p = mulupt(σ²c, ptn[1])::Float64

        #likelihood ratio
        if ωxc >= 0.0
          llr = σ²upd_llr(Xc, LApc, ωxc, σ²c, σ²p)::Float64
        else
          llr = σ²upd_llr(Xc, LAnc, ωxc, σ²c, σ²p)::Float64
        end

        # prior ratio
        prr = llrdexp_x(σ²p, σ²c, σ²prior)

        if -Random.randexp() < (llr + prr + log(σ²p/σ²c))
          llc += llr::Float64
          prc += prr::Float64
          σ²c  = σ²p::Float64
          lac[1] += 1
        end

      #update ωx
      elseif up == 2

        ωxp = addupt(ωxc, ptn[2])::Float64

        #likelihood ratio
        llr = ωxupd_llr(Xc, LApc, LAnc, ωxc, ωxp, σ²c)

        # prior ratio
        prr = llrdnorm_x(ωxp, ωxc, ωxprior[1], ωxprior[2])

        if -Random.randexp() < (llr + prr)
          llc += llr::Float64
          prc += prr::Float64
          ωxc  = ωxp::Float64
          lac[2] += 1
        end

      #update ω1
      elseif up == 3

        ω1p = addupt(ω1c, ptn[3])::Float64

        # proposal likelihood and prior
        llr = ω1upd_llr(Yc, λ1c, ω1c, ω1p, LDc)

        # prior ratio
        prr = llrdnorm_x(ω1p, ω1c, ω1prior[1], ω1prior[2])

        if -Random.randexp() < (llr + prr)
          llc += llr::Float64
          prc += prr::Float64
          ω1c  = ω1p::Float64
          lac[3] += 1
        end

        # update internal node
        if rand() < 0.4
          λϕ1, λϕ0 = λϕprop()
          llc = mhr_upd_Ytrio(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
        end

      # update ω0
      else

        ω0p = addupt(ω0c, ptn[4])

        # proposal likelihood and prior
        llr = ω0upd_llr(Yc, λ0c, ω0c, ω0p, LDc)::Float64

        # prior ratio
        prr = llrdnorm_x(ω0p, ω0c, ω0prior[1], ω0prior[2])

        if -Random.randexp() < (llr + prr)
          llc += llr
          prc += prr
          ω0c  = ω0p
          lac[4] += 1
        end

        # update internal node
        if rand() < 0.4
          λϕ1, λϕ0 = λϕprop()
          llc = mhr_upd_Ytrio(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
        end

      end

      ## make DA updates with some Pr
      # make Y branch update
      if rand() < 2e-3
        λϕ1, λϕ0 = λϕprop()
        llc = mhr_upd_Ybr(rand(Base.OneTo(nedge)), 
                        Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0,
                        llc, prc, LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
      end

      # make X branch update
      if rand() < 2e-3
        σ²ϕ = σ²ϕprop()
        llc = mhr_upd_Xbr(rand(Base.OneTo(nedge-1)),
                          Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, σ²ϕ, llc, 
                          LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
      end

      # make X trio update
      if rand() < 2e-3
        σ²ϕ = σ²ϕprop()
        llc = mhr_upd_Xtrio(rand(trios),
                            Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, σ²ϕ, llc, 
                            LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
      end

      # make joint X & Y branch update
      if rand() < 2e-3
        λϕ1, λϕ0 = λϕprop()
        σ²ϕ = σ²ϕprop()
        llc = mhr_upd_XYbr(rand(Base.OneTo(nedge-1)), 
                           Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, σ²ϕ, λϕ1, λϕ0,
                           llc, prc, LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
      end

      # make joint X & Y trio update
      if rand() < 2e-3
        λϕ1, λϕ0 = λϕprop()
        σ²ϕ = σ²ϕprop()
        llc      = mhr_upd_XYtrio(rand(trios), 
                           Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, σ²ϕ, λϕ1, λϕ0,
                           llc, prc, LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
      end

      # update stem branch
      if rand() < 2e-3
        λϕ1, λϕ0 = λϕprop()
        llc = mhr_upd_Ystem(λ1c, λ0c, λϕ1, λϕ0, llc, stemevc, brs)
      end

      # log number of updates
      ltn[up] += 1
      lup[up] += 1

      if (in(tune_int,ltn))
        wts = findall(isequal(tune_int),ltn)      # which to scale
        for j = wts
          ar     = lac[j]/lup[j]
          ptn[j] = scalef(ptn[j],ar)
          ltn[j] = 0
        end
      end
    end

    next!(p)
  end

  return llc, prc, Xc, Yc, LApc, LAnc, LDc, δXc, δYc, stemevc, 
         brs, σ²c, ωxc, ω1c, ω0c, λ1c, λ0c, ptn
end


