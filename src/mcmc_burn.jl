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
                      areaoc   ::Array{Int64,2},
                      linavg   ::Array{Float64,2},
                      lindiff  ::Array{Float64,3},
                      avg_Δx   ::Array{Float64,2},
                      λc       ::Array{Float64,1},
                      ωxc      ::Float64,
                      ω1c      ::Float64,
                      ω0c      ::Float64,
                      σ²c      ::Float64,
                      Ync1     ::Array{Int64,1},
                      Ync2     ::Array{Int64,1},
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
                      weight   ::NTuple{3,Float64},
                      λprior   ::Float64,
                      ωxprior  ::NTuple{2,Float64},
                      ω1prior  ::NTuple{2,Float64},
                      ω0prior  ::NTuple{2,Float64},
                      σ²prior  ::Float64,
                      fix_ωλ_ωμ::Bool,
                      np       ::Int64,
                      nburn    ::Int64,
                      obj_ar   ::Float64 = 0.234,
                      tune_int ::Int64   = 100)

  const m, ntip, narea  = size(Yc)

  const nedge = size(brs,1) 

  # likelihood and prior
  llc = total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λc, 
                  stemevc, brs[nedge,1,:], σ²c)
  prc = allλpr(λc, λprior)    + 
        logdexp(σ²c, σ²prior) + 
        logdnorm(ωxc, ωxprior[1], ωxprior[2]) +
        logdnorm(ω1c, ω1prior[1], ω1prior[2]) +
        logdnorm(ω0c, ω0prior[1], ω0prior[2])

  # make scaling function
  scalef = makescalef(obj_ar)

  # rest of tuning parameters
  const ptn = fill(.1,np) 

  # initialize acceptance log
  const ltn = zeros(Int64, np)
  const lup = zeros(Int64, np)
  const lac = zeros(Int64, np)

  # number of internal nodes
  const nin = length(trios) + 1

  # parameter location for λ
  const λlessthan = 6

  # progress bar
  p = Progress(nburn + 1, 5, "burning...", 20)

  if fix_ωλ_ωμ
    const pv      = append!(collect(1:np),
                            repeat(1:2, inner = ceil(Int64,np*weight[1])))
    append!(pv, repeat(5:6, inner = ceil(Int64,np*weight[3])))
    const parvec  = setdiff(pv,(3:4))
    const lparvec = length(parvec)

    print_with_color(:green,
      "\n ωx & σ² updates per iter = ", ceil(Int64,np*weight[1]) + 1,
      "\n λ1 & λ0 updates per iter = ", ceil(Int64,np*weight[3]) + 1, "\n")
  else
    const parvec  = append!(collect(1:np),
                            repeat(1:2, inner = ceil(Int64,np*weight[1])))
    append!(parvec, repeat(3:4, inner = ceil(Int64,np*weight[2])))
    append!(parvec, repeat(5:6, inner = ceil(Int64,np*weight[3])))
    const lparvec = length(parvec)

    print_with_color(:green,
      "\n ωx & σ² updates per iter = ", ceil(Int64,np*weight[1]) + 1,
      "\n ω1 & ω0 updates per iter = ", ceil(Int64,np*weight[2]) + 1,
      "\n λ1 & λ0 updates per iter = ", ceil(Int64,np*weight[3]) + 1, "\n")
  end


  #start brunin
  for it = Base.OneTo(nburn)

    # Update vector
    upvector = rand(parvec,lparvec)

    for up = upvector

      # update X
      if up > λlessthan

        Xp = copy(Xc)::Array{Float64,2}

        upx     = wXp[up - λlessthan]::Int64             # X indexing
        Xp[upx] = addupt(Xc[upx], ptn[up])::Float64      # update X

        k     = rowind(upx, m)::Int64
        wck   = wcol[k]::Array{Int64,1}

        if in(upx, Xnc1)        # if an internal node
          Xp[Xnc2] = Xp[Xnc1]::Array{Float64,1}
        end

        # calculate new averages
        aak, lak, ldk = Xupd_linavg(k, wck, Xp, Yc, narea)

        if upx == 1  # if root
          llr = Rupd_llf(k, wck, Xp, Yc, lak, ldk, ωxc, ω1c, ω0c, λc, σ²c) - 
                Rupd_llf(k, wck, Xc, Yc, linavg[k,wck], lindiff[k,wck,:], 
                          ωxc, ω1c, ω0c, λc, σ²c)::Float64
        else
          wckm1 = wcol[k-1]::Array{Int64,1}
          lakm1 = linavg[(k-1),wckm1]::Array{Float64,1}

          llr = Xupd_llf(k, wck, wckm1, Xp, Yc, lak, lakm1, 
                         ldk, ωxc, ω1c, ω0c, λc, σ²c) - 
                Xupd_llf(k, wck, wckm1, Xc, Yc, linavg[k,wck], lakm1,
                         lindiff[k,wck,:], ωxc, ω1c, ω0c, λc, σ²c)::Float64
        end
        
        if log(rand()) < llr
          Xc   = Xp::Array{Float64,2}
          llc += llr::Float64
          @inbounds begin
            areavg[k,:]      = aak
            linavg[k,wck]    = lak
            lindiff[k,wck,:] = ldk
            lac[up] += 1   # log acceptance
          end
        end

      #randomly select λ to update and branch histories
      elseif up > 4 && up <= λlessthan

        upλ = (up - 4)::Int64

        λp = copy(λc)::Array{Float64,1}

        # update λ
        λp[upλ] = logupt(λc[upλ], ptn[up])::Float64

        # proposal likelihood and prior
        llr = λupd_llf(Yc, λp, ω1c, ω0c, lindiff, stemevc, brs[nedge,1,:]) -
              λupd_llf(Yc, λc, ω1c, ω0c, lindiff, stemevc, brs[nedge,1,:])::Float64

        prr = logdexp(λp[upλ], λprior) - logdexp(λc[upλ], λprior)::Float64

        if log(rand()) < (llr + prr + 
                          log(λp[upλ])  - log(λc[upλ]))
          llc += llr::Float64
          prc += prr::Float64
          λc   = λp ::Array{Float64,1}
          lac[up] += 1
        end

        # which internal node to update
        if rand() < 0.4

          bup = rand(Base.OneTo(nin))::Int64
          if bup < nin

            Yp = copy(Yc)::Array{Int64,3}

            aa = copy(areavg)::Array{Float64,2}
            ao = copy(areaoc)::Array{Int64,2}
            la = copy(linavg)::Array{Float64,2}
            ld = copy(lindiff)::Array{Float64,3}

            linarea_branch_avg!(avg_Δx, lindiff, bridx_a, narea, nedge)

            upnode!(λc, ω1c, ω0c, avg_Δx, trios[bup], 
                    Yp, bridx_a, brδt, brl, brs, narea, nedge)

            Yp[Ync2] = Yp[Ync1]::Array{Int64,1}

            area_lineage_means!(aa, la, ao, Xc, Yp, wcol, m)
            linarea_diff!(ld, Xc, aa, ao, narea, ntip, m)

            llr = total_llf(Xc, Yp, la, ld, ωxc, ω1c, ω0c, λc, 
                            stemevc, brs[nedge,1,:], σ²c) - 
                  total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λc, 
                            stemevc, brs[nedge,1,:], σ²c)::Float64

            propr_iid = biogeo_upd_iid(Yc, λc, ω1c, ω0c, avg_Δx, trios[bup]) - 
                        biogeo_upd_iid(Yp, λc, ω1c, ω0c, avg_Δx, trios[bup])::Float64

            if log(rand()) < (llr + propr_iid)
              llc    += llr::Float64
              Yc      = Yp::Array{Int64,3}
              areavg  = aa::Array{Float64,2}
              areaoc  = ao::Array{Int64,2}
              linavg  = la::Array{Float64,2}
              lindiff = ld::Array{Float64,3}
            end

          else

            # update stem
            llr = 0.0
            for j=Base.OneTo(narea)
              @inbounds llr -= brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])::Float64
            end

            stemevc = upstem(λc, nedge, brs, brl, narea)

            for j=Base.OneTo(narea)
              @inbounds llr += brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])::Float64
            end

            llc += llr::Float64
          end
        end

      elseif up == 1         # if σ² is updated

        σ²p = logupt(σ²c, ptn[1])::Float64

        #likelihood ratio
        llr = σ²ωxupd_llf(Xc, linavg, ωxc, σ²p) - 
              σ²ωxupd_llf(Xc, linavg, ωxc, σ²c)::Float64

        # prior ratio
        prr = logdexp(σ²p, σ²prior) - logdexp(σ²c, σ²prior)::Float64

        if log(rand()) < (llr + prr + 
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
        llr = σ²ωxupd_llf(Xc, linavg, ωxp, σ²c) - 
              σ²ωxupd_llf(Xc, linavg, ωxc, σ²c)::Float64

        # prior ratio
        prr = logdnorm(ωxp, ωxprior[1], ωxprior[2]) -
              logdnorm(ωxc, ωxprior[1], ωxprior[2])::Float64

        if log(rand()) < (llr + prr)
          llc += llr::Float64
          prc += prr::Float64
          ωxc  = ωxp::Float64
          lac[2] += 1
        end

      #update ωλ
      elseif up == 3

        ωλp = addupt(ω1c, ptn[3])::Float64

            # proposal likelihood and prior
        llr = ωλμupd_llf(Yc, λc, ωλp, ω0c, lindiff) - 
              ωλμupd_llf(Yc, λc, ω1c, ω0c, lindiff)::Float64

        # prior ratio
        prr = logdnorm(ωλp, ω1prior[1], ω1prior[2]) -
              logdnorm(ω1c, ω1prior[1], ω1prior[2])::Float64

        if log(rand()) < (llr + prr)
          llc += llr::Float64
          prc += prr::Float64
          ω1c  = ωλp::Float64
          lac[3] += 1
        end

        # which internal node to update
        if rand() < 0.4

          bup = rand(Base.OneTo(nin))
          if bup < nin

            Yp = copy(Yc)::Array{Int64,3}

            aa = copy(areavg)::Array{Float64,2}
            ao = copy(areaoc)::Array{Int64,2}
            la = copy(linavg)::Array{Float64,2}
            ld = copy(lindiff)::Array{Float64,3}

            linarea_branch_avg!(avg_Δx, lindiff, bridx_a, narea, nedge)

            upnode!(λc, ω1c, ω0c, avg_Δx, trios[bup], 
                    Yp, bridx_a, brδt, brl, brs, narea, nedge)

            Yp[Ync2] = Yp[Ync1]::Array{Int64,1}

            area_lineage_means!(aa, la, ao, Xc, Yp, wcol, m)
            linarea_diff!(ld, Xc, aa, ao, narea, ntip, m)

            llr = total_llf(Xc, Yp, la, ld, ωxc, ω1c, ω0c, λc, 
                            stemevc, brs[nedge,1,:], σ²c) - 
                  total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λc, 
                            stemevc, brs[nedge,1,:], σ²c)::Float64

            propr_iid = biogeo_upd_iid(Yc, λc, ω1c, ω0c, avg_Δx, trios[bup]) - 
                        biogeo_upd_iid(Yp, λc, ω1c, ω0c, avg_Δx, trios[bup])::Float64

            if log(rand()) < (llr + propr_iid)
              llc    += llr::Float64
              Yc      = Yp::Array{Int64,3}
              areavg  = aa::Array{Float64,2}
              areaoc  = ao::Array{Int64,2}
              linavg  = la::Array{Float64,2}
              lindiff = ld::Array{Float64,3}
            end

          else

            # update stem
            llr = 0.0
            for j=Base.OneTo(narea)
              @inbounds llr -= brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])::Float64
            end

            stemevc = upstem(λc, nedge, brs, brl, narea)

            for j=Base.OneTo(narea)
              @inbounds llr += brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])::Float64
            end

            llc += llr::Float64
          end
        end

      # update ωμ
      else
        ωμp = addupt(ω0c, ptn[4])

        llr = ωλμupd_llf(Yc, λc, ω1c, ωμp, lindiff) - 
              ωλμupd_llf(Yc, λc, ω1c, ω0c, lindiff)

        # prior ratio
        prr = logdnorm(ωμp, ω0prior[1], ω0prior[2]) -
              logdnorm(ω0c, ω0prior[1], ω0prior[2])

        if log(rand()) < (llr + prr)
          llc += llr
          prc += prr
          ω0c  = ωμp
          lac[4] += 1
        end

                # which internal node to update
        if rand() < 0.4

          bup = rand(Base.OneTo(nin))
          if bup < nin

            Yp = copy(Yc)::Array{Int64,3}

            aa = copy(areavg)::Array{Float64,2}
            ao = copy(areaoc)::Array{Int64,2}
            la = copy(linavg)::Array{Float64,2}
            ld = copy(lindiff)::Array{Float64,3}

            linarea_branch_avg!(avg_Δx, lindiff, bridx_a, narea, nedge)

            upnode!(λc, ω1c, ω0c, avg_Δx, trios[bup], 
                    Yp, bridx_a, brδt, brl, brs, narea, nedge)

            Yp[Ync2] = Yp[Ync1]::Array{Int64,1}

            area_lineage_means!(aa, la, ao, Xc, Yp, wcol, m)
            linarea_diff!(ld, Xc, aa, ao, narea, ntip, m)

            llr = total_llf(Xc, Yp, la, ld, ωxc, ω1c, ω0c, λc, 
                            stemevc, brs[nedge,1,:], σ²c) - 
                  total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λc, 
                            stemevc, brs[nedge,1,:], σ²c)::Float64

            propr_iid = biogeo_upd_iid(Yc, λc, ω1c, ω0c, avg_Δx, trios[bup]) - 
                        biogeo_upd_iid(Yp, λc, ω1c, ω0c, avg_Δx, trios[bup])::Float64

            if log(rand()) < (llr + propr_iid)
              llc    += llr::Float64
              Yc      = Yp::Array{Int64,3}
              areavg  = aa::Array{Float64,2}
              areaoc  = ao::Array{Int64,2}
              linavg  = la::Array{Float64,2}
              lindiff = ld::Array{Float64,3}
            end

          else
            
            # update stem
            llr = 0.0
            for j=Base.OneTo(narea)
              @inbounds llr -= brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])::Float64
            end

            stemevc = upstem(λc, nedge, brs, brl, narea)

            for j=Base.OneTo(narea)
              @inbounds llr += brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])::Float64
            end

            llc += llr::Float64
          end
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

  return llc, prc, Xc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx, stemevc, brs, λc, ωxc, ω1c, ω0c, σ²c, ptn

end



