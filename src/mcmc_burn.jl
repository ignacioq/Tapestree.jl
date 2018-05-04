#=

Burning phase of the
Biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#




"""
    burn_compete(...)

Burning & adaptive phase for MCMC.
"""
function burn_compete(total_llf    ::Function,
                      λupd_llr     ::Function,
                      ω10upd_llr   ::Function,
                      Xupd_llr     ::Function,
                      Rupd_llr     ::Function,
                      σ²ωxupd_llr  ::Function,
                      bgiid        ::Function,
                      bgiid_br     ::Function,
                      mhr_upd_Xbr  ::Function,
                      mhr_upd_Y    ::Function,
                      mhr_upd_Ybr  ::Function,
                      mhr_upd_Ystem::Function,
                      mhr_upd_XYbr ::Function,
                      Xc       ::Array{Float64,2},
                      Yc       ::Array{Int64,3},
                      areavg   ::Array{Float64,2},
                      areaoc   ::Array{Int64,2},
                      linavg   ::Array{Float64,2},
                      lindiff  ::Array{Float64,3},
                      σ²c      ::Float64,
                      ωxc      ::Float64,
                      ω1c      ::Float64,
                      ω0c      ::Float64,
                      λ1c      ::Float64,
                      λ0c      ::Float64,
                      Xnc1     ::Array{Int64,1},
                      Xnc2     ::Array{Int64,1},
                      brl      ::Array{Float64,1},
                      wcol     ::Array{Array{Int64,1},1},
                      bridx_a  ::Array{Array{UnitRange{Int64},1},1},
                      brδt     ::Array{Array{Float64,1},1},
                      brs      ::Array{Int64,3},
                      stemevc  ::Array{Array{Float64,1},1},
                      trios    ::Array{Array{Int64,1},1},
                      wXp      ::Array{Int64,1},
                      λprior   ::Float64,
                      λϕprior  ::Float64,
                      ωxprior  ::NTuple{2,Float64},
                      ω1prior  ::NTuple{2,Float64},
                      ω0prior  ::NTuple{2,Float64},
                      σ²prior  ::Float64,
                      np       ::Int64,
                      parvec   ::Array{Int64,1},
                      nburn    ::Int64,
                      obj_ar   ::Float64 = 0.234,
                      tune_int ::Int64   = 100)

  const m, ntip, narea  = size(Yc)

  const nedge = size(brs,1) 

  # likelihood and prior
  llc = total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λ1c, λ0c,
                  stemevc, brs[nedge,1,:], σ²c)
  prc = allλpr(λ1c, λ0c, λprior)              +
        logdexp(σ²c, σ²prior)                 +
        logdnorm(ωxc, ωxprior[1], ωxprior[2]) +
        logdnorm(ω1c, ω1prior[1], ω1prior[2]) +
        logdnorm(ω0c, ω0prior[1], ω0prior[2])

  # make scaling function
  scalef = makescalef(obj_ar)

  # rest of tuning parameters
  const ptn = fill(.1, np) 

  # initialize acceptance log
  const ltn = zeros(Int64, np)
  const lup = zeros(Int64, np)
  const lac = zeros(Int64, np)

  # row i proposals for X
  const aai = zeros(Float64, narea)       # area average
  const lai = fill(NaN, ntip)             # lineage average
  const ldi = fill(NaN, ntip, narea)      # lineage difference

  # progress bar
  p = Progress(nburn, dt=5, desc="burn...", barlen=20, color=:green)

  # print number of parameters
  print_with_color(:green,
    "\nσ² updates per iter = ", endof(filter(x -> x == 1,parvec)),
    "\nωx updates per iter = ", endof(filter(x -> x == 2,parvec)),
    "\nω1 updates per iter = ", endof(filter(x -> x == 3,parvec)),
    "\nω0 updates per iter = ", endof(filter(x -> x == 4,parvec)),
    "\nλ1 updates per iter = ", endof(filter(x -> x == 5,parvec)),
    "\nλ0 updates per iter = ", endof(filter(x -> x == 6,parvec)),
    "\nθ  updates per iter = ", length(parvec), "\n")

  #start burnin
  for it = Base.OneTo(nburn)

    # Update vector
    shuffle!(parvec)

    for up = parvec

      # update X
      if up > 6

        upx = wXp[up - 6]::Int64                 # X indexing

        xi, xj = ind2sub(Xc, upx)

        const xpi = Xc[xi,:]::Array{Float64,1}

        xpi[xj] = addupt(xpi[xj], ptn[up])::Float64      # update X

        if in(upx, Xnc1)        # if an internal node
          xpi[ind2sub(Xc, Xnc2[findfirst(Xnc1, upx)])[2]] = xpi[xj]::Float64
        end

        # calculate new averages
        Xupd_linavg!(aai, lai, ldi, areaoc, xi, wcol[xi], xpi, Yc, narea)

        if upx == 1  # if root
          @inbounds llr = Rupd_llr(wcol[1], 
                                   xpi[wcol[1]], 
                                   Xc[1,wcol[1]], Xc[2,wcol[1]], 
                                   lai[wcol[1]], ldi[wcol[1],:], 
                                   linavg[1,wcol[1]], lindiff[1,wcol[1],:],
                                   Yc, 
                                   ωxc, ω1c, ω0c, λ1c, λ0c, σ²c)::Float64
        else
          @inbounds llr = Xupd_llr(xi, wcol[xi], wcol[xi-1], 
                                   xpi, 
                                   Xc[xi,:], Xc[xi-1,:], Xc[xi+1,:], 
                                   lai, ldi, 
                                   linavg[xi,:], linavg[xi-1,:], 
                                   lindiff[xi,:,:],
                                   Yc, 
                                   ωxc, ω1c, ω0c, λ1c, λ0c, σ²c)::Float64
        end

        if -randexp() < llr
          @inbounds begin
            Xc[xi,:]        = xpi::Array{Float64,1}
            llc            += llr::Float64
            areavg[xi,:]    = aai::Array{Float64,1}
            linavg[xi,:]    = lai::Array{Float64,1}
            lindiff[xi,:,:] = ldi::Array{Float64,2}
            lac[up]        += 1   # log acceptance
          end
        end

      # update λ1
      elseif up == 5

        # update λ
        λ1p = mulupt(λ1c, ptn[5])::Float64

        llr = λupd_llr(Yc, λ1c, λ0c, λ1p, λ0c, ω1c, ω0c, 
                       lindiff, stemevc, brs[nedge,1,:])::Float64

        prr = logdexp(λ1p, λprior) - logdexp(λ1c, λprior)::Float64

        if -randexp() < (llr + prr + log(λ1p)  - log(λ1c))
          llc     += llr::Float64
          prc     += prr::Float64
          λ1c      = λ1p::Float64
          lac[5] += 1
        end

        # which internal node to update
        if rand() < 0.4
          λϕ1, λϕ0 = λϕprop()

          llc = mhr_upd_Y(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    areavg, areaoc, linavg, lindiff, brs, stemevc)
        end

      # update λ0
      elseif up == 6

        # update λ
        λ0p = mulupt(λ0c, ptn[6])::Float64

        llr = λupd_llr(Yc, λ1c, λ0c, λ1c, λ0p, ω1c, ω0c, 
                       lindiff, stemevc, brs[nedge,1,:])::Float64

        prr = logdexp(λ0p, λprior) - logdexp(λ0c, λprior)::Float64

        if -randexp() < (llr + prr + log(λ0p)  - log(λ0c))
          llc     += llr::Float64
          prc     += prr::Float64
          λ0c      = λ0p::Float64
          lac[6] += 1
        end

        # which internal node to update
        if rand() < 0.4
          λϕ1, λϕ0 = λϕprop()

          llc = mhr_upd_Y(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    areavg, areaoc, linavg, lindiff, brs, stemevc)
        end

      # if σ² is updated
      elseif up == 1         

        σ²p = mulupt(σ²c, ptn[1])::Float64

        #likelihood ratio
        llr = σ²ωxupd_llr(Xc, linavg, ωxc, ωxc, σ²c, σ²p)::Float64

        # prior ratio
        prr = logdexp(σ²p, σ²prior) - logdexp(σ²c, σ²prior)::Float64

        if -randexp() < (llr + prr + 
                         log(σ²p) - log(σ²c))
          llc += llr::Float64
          prc += prr::Float64
          σ²c  = σ²p::Float64
          lac[1] += 1
        end

      #update ωx
      elseif up == 2
        ωxp = addupt(ωxc, ptn[2])::Float64

        #likelihood ratio
        llr = σ²ωxupd_llr(Xc, linavg, ωxc, ωxp, σ²c, σ²c)

        # prior ratio
        prr = logdnorm(ωxp, ωxprior[1], ωxprior[2]) -
              logdnorm(ωxc, ωxprior[1], ωxprior[2])::Float64

        if -randexp() < (llr + prr)
          llc += llr::Float64
          prc += prr::Float64
          ωxc  = ωxp::Float64
          lac[2] += 1
        end

      #update ω1
      elseif up == 3

        ω1p = addupt(ω1c, ptn[3])::Float64

        # proposal likelihood and prior
        llr = ω10upd_llr(Yc, λ1c, λ0c, ω1c, ω0c, ω1p, ω0c, lindiff)::Float64

        # prior ratio
        prr = logdnorm(ω1p, ω1prior[1], ω1prior[2]) -
              logdnorm(ω1c, ω1prior[1], ω1prior[2])::Float64

        if -randexp() < (llr + prr)
          llc += llr::Float64
          prc += prr::Float64
          ω1c  = ω1p::Float64
          lac[3] += 1
        end

        # which internal node to update
        if rand() < 0.4
          λϕ1, λϕ0 = λϕprop()

          llc = mhr_upd_Y(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    areavg, areaoc, linavg, lindiff, brs, stemevc)
        end

      # update ω0
      else
        ω0p = addupt(ω0c, ptn[4])

        # proposal likelihood ratio
        llr = ω10upd_llr(Yc, λ1c, λ0c, ω1c, ω0c, ω1c, ω0p, lindiff)::Float64

        # prior ratio
        prr = logdnorm(ω0p, ω0prior[1], ω0prior[2])::Float64 -
              logdnorm(ω0c, ω0prior[1], ω0prior[2])::Float64

        if -randexp() < (llr + prr)
          llc += llr
          prc += prr
          ω0c  = ω0p
          lac[4] += 1
        end

        # which internal node to update
        if rand() < 0.4
          λϕ1, λϕ0 = λϕprop()

          llc = mhr_upd_Y(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    areavg, areaoc, linavg, lindiff, brs, stemevc)
        end

      end

      ## make a branch updates with some Pr
      # make Y branch update
      if rand() < 2e-3
        λϕ1, λϕ0 = λϕprop()

        llc = mhr_upd_Ybr(rand(Base.OneTo(nedge)), 
                          Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, 
                          λϕ1, λϕ0,
                          llc, prc, areavg, areaoc, linavg, lindiff, 
                          brs, stemevc)
      end

      if rand() < 2e-3
        llc = mhr_upd_Xbr(rand(Base.OneTo(nedge-1)),
                          Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, llc, 
                          areavg, linavg, lindiff, areaoc, brs, stemevc)
      end

      # make joint X & Y branch update
      if rand() < 2e-3
        λϕ1, λϕ0 = λϕprop()

        llc = mhr_upd_XYbr(rand(Base.OneTo(nedge-1)), 
                           Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, 
                           λϕ1, λϕ0,
                           llc, prc, areavg, areaoc, linavg, lindiff, 
                           brs, stemevc)
      end

      # update stem node
      if rand() < 2e-3
        λϕ1, λϕ0 = λϕprop()

        llc = mhr_upd_Ystem(λ1c, λ0c, λϕ1, λϕ0, llc, stemevc, brs)
      end

      # log number of updates
      ltn[up] += 1
      lup[up] += 1

      if (in(tune_int,ltn))
        wts = find(map(x -> x == tune_int,ltn))      # which to scale
        for j = wts
          ar     = lac[j]/lup[j]
          ptn[j] = scalef(ptn[j],ar)
          ltn[j] = 0
        end
      end
    end

    next!(p)
  end

  return llc, prc, Xc, Yc, areavg, areaoc, linavg, lindiff, stemevc, 
         brs, σ²c, ωxc, ω1c, ω0c, λ1c, λ0c, ptn
end


