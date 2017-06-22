
"""

Burning phase of the
Biogeographic competition model

Ignacio Quintero

t(-_-t)

April 27 2017

"""




# burning & adaptive phase
function burn_compete(total_llf,
                      λupd_llf,
                      ωλμupd_llf,
                      Xupd_llf,
                      Rupd_llf,
                      σ²ωxupd_llf,
                      biogeo_upd_iid,
                      Xc       ::Array{Float64,2},
                      Yc       ::Array{Int64,3},
                      areavg   ::Array{Float64,2},
                      linavg   ::Array{Float64,2},
                      lindiff  ::Array{Float64,3},
                      λc       ::Array{Float64,2},
                      ωxc      ::Float64,
                      ωλc      ::Float64,
                      ωμc      ::Float64,
                      σ²c      ::Float64,
                      Ync1     ::Array{Int64,1},
                      Ync2     ::Array{Int64,1},
                      Xnc1     ::Array{Int64,1},
                      Xnc2     ::Array{Int64,1},
                      brl      ::Array{Float64,1},
                      wcol     ::Array{Array{Int64,1},1},
                      bridx_a  ::Array{Array{Array{Int64,1},1},1},
                      brδt     ::Array{Array{Float64,1},1},
                      brs      ::Array{Int64,3},
                      stemevc  ::Array{Array{Float64,1},1},
                      trios    ::Array{Array{Int64,1},1},
                      wXp      ::Array{Int64,1},
                      λprior   ::Float64,
                      ωxprior  ::Tuple{Float64,Float64},
                      ωλprior  ::Tuple{Float64,Float64},
                      ωμprior  ::Tuple{Float64,Float64},
                      σ²prior  ::Float64,
                      fix_ωλ_ωμ::Bool,
                      np       ::Int64,
                      nburn    ::Int64,
                      obj_ar   ::Float64 = 0.234,
                      tune_int ::Int64   = 100)

  const m, ntip, narea  = size(Yc)

  const nedge = size(brs,1) 

  # likelihood and prior
  llc = total_llf(Xc, Yc, linavg, lindiff, ωxc, ωλc, ωμc, λc, 
                  stemevc, brs[nedge,1,:], σ²c)
  prc = allλpr(λc, λprior)    + 
        logdexp(σ²c, σ²prior) + 
        logdnorm(ωxc, ωxprior[1], ωxprior[2]) +
        logdnorm(ωλc, ωλprior[1], ωλprior[2]) +
        logdnorm(ωμc, ωμprior[1], ωμprior[2])

  # make scaling function
  scalef = makescalef(obj_ar)

  # rest of tuning parameters
  ptn = fill(.1,np) 

  # initialize acceptance log
  ltn = zeros(Int64, np)
  lup = zeros(Int64, np)
  lac = zeros(Int64, np)

  # number of internal nodes
  const nin = length(trios) + 1

  # parameter location for λ
  const λlessthan = 2*narea + 4

  # progress bar
  p = Progress(nburn + 1, 5, "burning...", 20)

  if fix_ωλ_ωμ
    const pv  = append!(collect(1:np),fill(1, floor(Int64,np*.15)))
    const parvec = setdiff(pv,(3:4))
    const lparvec = length(parvec)
  else
    const parvec  = append!(collect(1:np),fill(1, floor(Int64,np*.15)))
    const lparvec = length(parvec)
  end

  #start brunin
  for it = Base.OneTo(nburn)

    # Update vector
    upvector = rand(parvec,lparvec)

    for up = upvector

      # update X
      if up > λlessthan

        Xp = copy(Xc)

        upx     = wXp[up - λlessthan]           # X indexing
        Xp[upx] = addupt(Xc[upx], ptn[up])      # update X

        k     = rowind(upx, m)
        wck   = wcol[k]

        if in(upx, Xnc1)        # if an internal node
          Xp[Xnc2] = Xp[Xnc1]
        end

        # calculate new averages
        aak, lak, ldk = Xupd_linavg(k, wck, Xp, Yc, narea)

        if upx == 1  # if root
          llr = Rupd_llf(k, wck, Xp, Yc, lak, ldk, ωxc, ωλc, ωμc, λc, σ²c) - 
                Rupd_llf(k, wck, Xc, Yc, linavg[k,wck], lindiff[k,wck,:], 
                          ωxc, ωλc, ωμc, λc, σ²c)
        else
          wckm1 = wcol[k-1]
          lakm1 = linavg[(k-1),wckm1]

          llr = Xupd_llf(k, wck, wckm1, Xp, Yc, lak, lakm1, 
                         ldk, ωxc, ωλc, ωμc, λc, σ²c) - 
                Xupd_llf(k, wck, wckm1, Xc, Yc, linavg[k,wck], lakm1,
                         lindiff[k,wck,:], ωxc, ωλc, ωμc, λc, σ²c)
        end
        
        if log(rand()) < llr
          Xc   = Xp
          llc += llr
          @inbounds begin
            areavg[k,:]      = aak
            linavg[k,wck]    = lak
            lindiff[k,wck,:] = ldk
            lac[up] += 1   # log acceptance
          end
        end

      #randomly select λ to update and branch histories
      elseif up > 4 && up <= λlessthan

        upλ = up - 4

        λp = copy(λc)

        # update λ
        λp[upλ] = logupt(λc[upλ], ptn[up])

        # proposal likelihood and prior
        llr = λupd_llf(Yc, λp, ωλc, ωμc, lindiff, stemevc, brs[nedge,1,:]) -
              λupd_llf(Yc, λc, ωλc, ωμc, lindiff, stemevc, brs[nedge,1,:])

        prr = logdexp(λp[upλ], λprior) - logdexp(λc[upλ], λprior)

        if log(rand()) < (llr + prr + 
                          log(λp[upλ])  - log(λc[upλ]))
          llc += llr
          prc += prr
          λc   = λp
          lac[up] += 1
        end

        # which internal node to update
        bup = rand(Base.OneTo(nin))
        if bup < nin

          Yp = copy(Yc)

          aa = copy(areavg)
          la = copy(linavg)
          ld = copy(lindiff)

          upnode!(λc, trios[bup], Yp, bridx_a, brδt, brl, brs, narea, nedge)
          
          Yp[Ync2] = Yp[Ync1]

          area_lineage_means!(aa, la, Xc, Yp, wcol, m)
          linarea_diff!(ld, Xc, aa, narea, ntip, m)

          llr = total_llf(Xc, Yp, la, ld, ωxc, ωλc, ωμc, λc, 
                          stemevc, brs[nedge,1,:], σ²c) - 
                total_llf(Xc, Yc, linavg, lindiff, ωxc, ωλc, ωμc, λc, 
                                stemevc, brs[nedge,1,:], σ²c)

          propr_iid = biogeo_upd_iid(Yc, λc, trios[bup]) - 
                      biogeo_upd_iid(Yp, λc, trios[bup])

          if log(rand()) < (llr + propr_iid)
            llc    += llr
            Yc      = Yp
            areavg  = aa
            linavg  = la
            lindiff = ld
          end

        else

          # update stem
            llr = 0.0
            for j=Base.OneTo(narea)
              @inbounds llr -= brll(stemevc[j], λc[j,1], λc[j,2], brs[nedge,1,j])
            end

            stemevc = upstem(λc, nedge, brs, brl, narea)

            for j=Base.OneTo(narea)
              @inbounds llr += brll(stemevc[j], λc[j,1], λc[j,2], brs[nedge,1,j])
            end

            llc += llr
        end


      elseif up == 1         # if σ² is updated

        σ²p = logupt(σ²c, ptn[1])

        #likelihood ratio
        llr = σ²ωxupd_llf(Xc, linavg, ωxc, σ²p) - 
              σ²ωxupd_llf(Xc, linavg, ωxc, σ²c)

        # prior ratio
        prr = logdexp(σ²p, σ²prior) - logdexp(σ²c, σ²prior)

        if log(rand()) < (llr + prr + 
                          log(σ²p) - log(σ²c))
          llc += llr
          prc += prr
          σ²c  = σ²p
          lac[1] += 1
        end

      #update ωx
      elseif up == 2
        ωxp = -abs(addupt(ωxc, ptn[2]))

        #likelihood ratio
        llr = σ²ωxupd_llf(Xc, linavg, ωxp, σ²c) - 
              σ²ωxupd_llf(Xc, linavg, ωxc, σ²c)

        # prior ratio
        prr = logdnorm(ωxp, ωxprior[1], ωxprior[2]) -
              logdnorm(ωxc, ωxprior[1], ωxprior[2])

        if log(rand()) < (llr + prr)
          llc += llr
          prc += prr
          ωxc  = ωxp
          lac[2] += 1
        end

      #update ωλ
      elseif up == 3

        ωλp = addupt(ωλc, ptn[3])

            # proposal likelihood and prior
        llr = ωλμupd_llf(Yc, λc, ωλp, ωμc, lindiff) - 
              ωλμupd_llf(Yc, λc, ωλc, ωμc, lindiff)

        # prior ratio
        prr = logdnorm(ωλp, ωλprior[1], ωλprior[2]) -
              logdnorm(ωλc, ωλprior[1], ωλprior[2])

        if log(rand()) < (llr + prr)
          llc += llr
          prc += prr
          ωλc  = ωλp
          lac[3] += 1
        end

      # update ωμ
      else
        ωμp = addupt(ωμc, ptn[4])

        llr = ωλμupd_llf(Yc, λc, ωλc, ωμp, lindiff) - 
              ωλμupd_llf(Yc, λc, ωλc, ωμc, lindiff)

        # prior ratio
        prr = logdnorm(ωμp, ωμprior[1], ωμprior[2]) -
              logdnorm(ωμc, ωμprior[1], ωμprior[2])

        if log(rand()) < (llr + prr)
          llc += llr
          prc += prr
          ωμc  = ωμp
          lac[4] += 1
        end
      end

      # log number of updates
      ltn[up] += 1
      lup[up] += 1

      if (in(tune_int,ltn))
        wts = find(ltn .== tune_int)      # which to scale
        for j = wts
          ar     = lac[j]/lup[j]
          ptn[j] = scalef(ptn[j],ar)
          ltn[j] = 0
        end
      end

    end
    
    next!(p)
  end

  return llc, prc, Xc, Yc, areavg, linavg, lindiff, stemevc, brs, λc, ωxc, ωλc, ωμc, σ²c, ptn

end



