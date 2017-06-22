
"""

Updates for MCMC

Ignacio Quintero

t(-_-t)

June 20 2017

"""




# update λ
function make_mhr_upd_λ(nedge ::Int64, 
                        λprior::Float64,
                        ptn   ::Array{Float64},
                        λupd_llf)

  function f(up     ::Int64,
             Yc     ::Array{Int64,3},
             λc     ::Array{Float64,2},
             llc    ::Float64,
             prc    ::Float64,
             ωλc    ::Float64,
             ωμc    ::Float64,
             lindiff::Array{Float64,3},
             stemevc::Array{Array{Float64,1},1},
             stemss ::Array{Int64,1})
    upλ = up - 4
    
    λp = copy(λc)

    # update λ
    λp[upλ] = logupt(λc[upλ], rand() < 0.5 ? ptn[up] : 10*ptn[up])

    # proposal likelihood and prior
    llr = λupd_llf(Yc, λp, ωλc, ωμc, lindiff, stemevc, stemss) -
          λupd_llf(Yc, λc, ωλc, ωμc, lindiff, stemevc, stemss)

    prr = logdexp(λp[upλ], λprior) - logdexp(λc[upλ], λprior)

    if log(rand()) < (llr + prr + 
                      log(λp[upλ])  - log(λc[upλ]))
      llc += llr
      prc += prr
      λc   = λp
    end
  end

end




# update trio in Y
function make_mhr_upd_Y(narea  ::Int64,
                        nedge  ::Int64,
                        m      ::Int64,
                        ntip   ::Int64,
                        bridx_a::Vector{Vector{Vector{Int64}}},
                        brδt   ::Array{Array{Float64,1},1},
                        brl    ::Array{Float64,1},
                        wcol   ::Array{Array{Int64,1},1},
                        Ync1   ::Array{Int64,1},
                        Ync2   ::Array{Int64,1},
                        total_llf,
                        biogeo_upd_iid)

  function f(triad  ::Array{Int64,1},
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λc     ::Array{Float64,2},
             ωxc    ::Float64, 
             ωλc    ::Float64, 
             ωμc    ::Float64,
             σ²c    ::Float64,
             llc    ::Float64,
             prc    ::Float64,
             areavg ::Array{Float64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    Yp = copy(Yc)
    aa = copy(areavg)
    la = copy(linavg)
    ld = copy(lindiff)

    upnode!(λc, triad, Yp, bridx_a, brδt, brl, brs, narea, nedge)

    Yp[Ync2] = Yp[Ync1]

    area_lineage_means!(aa, la, Xc, Yp, wcol, m)
    linarea_diff!(ld, Xc, aa, narea, ntip, m)

    @views begin
      llr = total_llf(Xc, Yp, la, ld, ωxc, ωλc, ωμc, λc, 
                      stemevc, brs[nedge,1,:], σ²c) - 
            total_llf(Xc, Yc, linavg, lindiff, ωxc, ωλc, ωμc, λc, 
                      stemevc, brs[nedge,1,:], σ²c)
    end

    propr_iid = biogeo_upd_iid(Yc, λc, triad) - 
                biogeo_upd_iid(Yp, λc, triad)

    if log(rand()) < (llr + propr_iid)
      llc    += llr
      Yc      = Yp
      areavg  = aa
      linavg  = la
      lindiff = ld
    end

  end

end




# Update X
function make_mhr_upd_X(Xnc1     ::Array{Int64,1},
                        Xnc2     ::Array{Int64,1},
                        wcol     ::Array{Array{Int64,1},1},
                        m        ::Int64,
                        ptn      ::Array{Float64,1},
                        wXp      ::Array{Int64,1},
                        λlessthan::Int64,
                        narea    ::Int64,
                        Xupd_llf,
                        Rupd_llf)

  function f(up     ::Int64,
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λc     ::Array{Float64,2},
             ωxc    ::Float64, 
             ωλc    ::Float64, 
             ωμc    ::Float64,
             σ²c    ::Float64,
             llc    ::Float64,
             areavg ::Array{Float64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3})

    upx = wXp[up - λlessthan] 

    Xp      = copy(Xc)
    Xp[upx] = addupt(Xc[upx], ptn[up])      # update X

    k     = rowind(upx, m)
    wck   = wcol[k]

    if ∈(upx, Xnc1)        # if an internal node
      Xp[Xnc2] = Xp[Xnc1]
    end

    # calculate new averages
    aak, lak, ldk = Xupd_linavg(k, wck, Xp, Yc, narea)

    if upx == 1  # if root

      @views begin
        llr = Rupd_llf(k, wck, Xp, Yc, lak, ldk, ωxc, ωλc, ωμc, λc, σ²c) - 
              Rupd_llf(k, wck, Xc, Yc, linavg[k,wck], lindiff[k,wck,:], 
                       ωxc, ωλc, ωμc, λc, σ²c)
      end

    else

      wckm1 = wcol[k-1]
      lakm1 = linavg[(k-1),wckm1]

      @views begin
        llr = Xupd_llf(k, wck, wckm1, Xp, Yc, lak, lakm1, 
                       ldk, ωxc, ωλc, ωμc, λc, σ²c) - 
              Xupd_llf(k, wck, wckm1, Xc, Yc, linavg[k,wck], lakm1,
                       lindiff[k,wck,:], ωxc, ωλc, ωμc, λc, σ²c)
      end

    end

    if log(rand()) < llr
      Xc   = Xp
      llc += llr
      @inbounds begin
        areavg[k,:]      = aak
        linavg[k,wck]    = lak
        lindiff[k,wck,:] = ldk
      end
    end

  end

end





# update σ²
function mhr_upd_σ²!(σ²c    ::Float64,
                     Xc     ::Array{Float64,2},
                     ωxc    ::Float64,
                     llc    ::Float64,
                     prc    ::Float64,
                     σ²tn   ::Float64,
                     linavg ::Array{Float64,2},
                     σ²prior::Float64,
                     σ²ωxupd_llf)

  σ²p = logupt(σ²c, rand() < 0.5 ? σ²tn : 10*σ²tn)

  #likelihood ratio
  llr = σ²ωxupd_llf(Xc, linavg, ωxc, σ²p) - 
        σ²ωxupd_llf(Xc, linavg, ωxc, σ²c)

  # prior ratio
  prr = logdexp(σ²p, σ²prior) - logdexp(σ²c, σ²prior)

  if log(rand()) < (llr + 
                    prr + 
                    log(σ²p) - log(σ²c))
    llc += llr
    prc += prr
    σ²c  = σ²p
  end

end





# update ωx
function mhr_upd_ωx!(ωxc    ::Float64,
                     Xc     ::Array{Float64,2},
                     σ²c    ::Float64,
                     llc    ::Float64,
                     prc    ::Float64,
                     ωxtn   ::Float64,
                     linavg ::Array{Float64,2},
                     ωxprior::Tuple{Float64,Float64},
                     σ²ωxupd_llf)

  ωxp = addupt(ωxc, rand() < 0.5 ? ωxtn : 10*ωxtn)

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
  end

end





#update ωλ
function mhr_upd_ωλ!(ωλc    ::Float64,
                     λc     ::Array{Float64,2},
                     ωμc    ::Float64,
                     Yc     ::Array{Int64,3},
                     llc    ::Float64,
                     prc    ::Float64,
                     ωλtn   ::Float64,
                     linavg ::Array{Float64,2},
                     lindiff::Array{Float64,3},
                     ωλprior::Tuple{Float64,Float64},
                     ωλμupd_llf)

  ωλp = addupt(ωλc, rand() < 0.5 ? ωλtn : 10*ωλtn)

  # likelihood ratio
  llr = ωλμupd_llf(Yc, λc, ωλp, ωμc, lindiff) - 
        ωλμupd_llf(Yc, λc, ωλc, ωμc, lindiff)

  # prior ratio
  prr = logdnorm(ωλp, ωλprior[1], ωλprior[2]) -
        logdnorm(ωλc, ωλprior[1], ωλprior[2])

  if log(rand()) < (llr + prr)
    llc += llr
    prc += prr
    ωλc  = ωλp
  end

end





# update ωμ
function mhr_upd_ωμ!(ωμc    ::Float64,
                     λc     ::Array{Float64,2},
                     ωλc    ::Float64,
                     Yc     ::Array{Int64,3},
                     llc    ::Float64,
                     prc    ::Float64,
                     ωμtn   ::Float64,
                     linavg ::Array{Float64,2},
                     lindiff::Array{Float64,3},
                     ωμprior::Tuple{Float64,Float64},
                     ωλμupd_llf)

  ωμp = addupt(ωμc, rand() < 0.5 ? ωμtn : 10*ωμtn)

  # likelihood ratio
  llr = ωλμupd_llf(Yc, λc, ωλc, ωμp, lindiff) - 
        ωλμupd_llf(Yc, λc, ωλc, ωμc, lindiff)

  # prior ratio
  prr = logdnorm(ωμp, ωμprior[1], ωμprior[2]) -
        logdnorm(ωμc, ωμprior[1], ωμprior[2])

  if log(rand()) < (llr + prr)
    llc += llr
    prc += prr
    ωμc  = ωμp
  end

end


