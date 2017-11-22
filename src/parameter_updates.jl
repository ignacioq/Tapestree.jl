#=

Updates for MCMC

Ignacio Quintero Mächler

t(-_-t)

June 20 2017

=#



"""
    make_mhr_upd_λ(nedge::Int64, λprior::Float64, ptn::Array{Float64}, λupd_llf)

Make function to update λ.
"""
function make_mhr_upd_λ(nedge   ::Int64, 
                        λprior  ::Float64,
                        ptn     ::Array{Float64},
                        λupd_llf::Function)

  function f(up     ::Int64,
             Yc     ::Array{Int64,3},
             λc     ::Array{Float64,1},
             llc    ::Float64,
             prc    ::Float64,
             ω1c    ::Float64,
             ω0c    ::Float64,
             lindiff::Array{Float64,3},
             stemevc::Array{Array{Float64,1},1},
             stemss ::Array{Int64,1})
    
    upλ = (up - 4)::Int64

    λp = copy(λc)::Array{Float64,1}

    # update λ
    λp[upλ] = logupt(λc[upλ], rand() < 0.5 ? ptn[up] : 4*ptn[up])::Float64

    # proposal likelihood and prior
    llr = λupd_llf(Yc, λp, ω1c, ω0c, lindiff, stemevc, stemss) -
          λupd_llf(Yc, λc, ω1c, ω0c, lindiff, stemevc, stemss)::Float64

    prr = logdexp(λp[upλ], λprior) - logdexp(λc[upλ], λprior)::Float64

    if log(rand()) < (llr + prr + 
                      log(λp[upλ])  - log(λc[upλ]))
      llc += llr::Float64
      prc += prr::Float64
      λc   = λp ::Array{Float64,1}
    end
  
    return (llc, prc, λc)::Tuple{Float64,Float64,Array{Float64,1}}
  end

end





"""
    make_mhr_upd_Y(narea::Int64, nedge::Int64, m::Int64, ntip::Int64, bridx_a::Vector{Vector{Vector{Int64}}}, brδt::Array{Array{Float64,1},1}, brl::Array{Float64,1}, wcol::Array{Array{Int64,1},1}, Ync1::Array{Int64,1}, Ync2::Array{Int64,1}, total_llf, biogeo_upd_iid)

Make function to update trio in Y.
"""
function make_mhr_upd_Y(narea          ::Int64,
                        nedge          ::Int64,
                        m              ::Int64,
                        ntip           ::Int64,
                        bridx_a        ::Array{Array{UnitRange{Int64},1},1},
                        brδt           ::Array{Array{Float64,1},1},
                        brl            ::Array{Float64,1},
                        wcol           ::Array{Array{Int64,1},1},
                        Ync1           ::Array{Int64,1},
                        Ync2           ::Array{Int64,1},
                        total_llf      ::Function,
                        biogeo_upd_iid ::Function)

  function f(triad  ::Array{Int64,1},
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λc     ::Array{Float64,1},
             ωxc    ::Float64, 
             ω1c    ::Float64, 
             ω0c    ::Float64,
             σ²c    ::Float64,
             llc    ::Float64,
             prc    ::Float64,
             areavg ::Array{Float64,2},
             areaoc ::Array{Int64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             avg_Δx ::Array{Float64,2},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})


    Yp = copy(Yc)::Array{Int64,3}
    aa = copy(areavg)::Array{Float64,2}
    ao = copy(areaoc)::Array{Int64,2}
    la = copy(linavg)::Array{Float64,2}
    ld = copy(lindiff)::Array{Float64,3}

    linarea_branch_avg!(avg_Δx, lindiff, bridx_a, narea, nedge)

    upnode!(λc, ω1c, ω0c, avg_Δx, triad, 
            Yp, bridx_a, brδt, brl, brs, narea, nedge)

    Yp[Ync2] = Yp[Ync1]::Array{Int64,1}

    area_lineage_means!(aa, la, ao, Xc, Yp, wcol, m)
    linarea_diff!(ld, Xc, aa, ao, narea, ntip, m)

    llr = (total_llf(Xc, Yp, la, ld, ωxc, ω1c, ω0c, λc, 
                     stemevc, brs[nedge,1,:], σ²c) - 
           total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λc, 
                     stemevc, brs[nedge,1,:], σ²c))::Float64

    propr_iid = (biogeo_upd_iid(Yc, λc, ω1c, ω0c, avg_Δx, triad) - 
                 biogeo_upd_iid(Yp, λc, ω1c, ω0c, avg_Δx, triad))::Float64

    if log(rand()) < (llr + propr_iid)::Float64
      llc    += llr::Float64
      Yc      = Yp ::Array{Int64,3}
      areavg  = aa ::Array{Float64,2}
      areaoc  = ao ::Array{Int64,2}
      linavg  = la ::Array{Float64,2}
      lindiff = ld ::Array{Float64,3}
    end

    return llc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx
  end

end




"""
    make_mhr_upd_X(Xnc1::Array{Int64,1}, Xnc2::Array{Int64,1}, wcol::Array{Array{Int64,1},1}, m::Int64, ptn::Array{Float64,1}, wXp::Array{Int64,1}, λlessthan::Int64, narea::Int64, Xupd_llf, Rupd_llf)

Make DA update X.
"""
function make_mhr_upd_X(Xnc1     ::Array{Int64,1},
                        Xnc2     ::Array{Int64,1},
                        wcol     ::Array{Array{Int64,1},1},
                        m        ::Int64,
                        ptn      ::Array{Float64,1},
                        wXp      ::Array{Int64,1},
                        λlessthan::Int64,
                        narea    ::Int64,
                        Xupd_llf ::Function,
                        Rupd_llf ::Function)

  function f(up     ::Int64,
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λc     ::Array{Float64,1},
             ωxc    ::Float64, 
             ω1c    ::Float64, 
             ω0c    ::Float64,
             σ²c    ::Float64,
             llc    ::Float64,
             areavg ::Array{Float64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3})

    upx = wXp[up - λlessthan]::Int64

    Xp      = copy(Xc)::Array{Float64,2}
    Xp[upx] = addupt(Xc[upx], ptn[up])::Float64      # update X

    k     = rowind(upx, m)::Int64
    wck   = wcol[k]::Array{Int64,1}

    if ∈(upx, Xnc1)                         # if an internal node
      Xp[Xnc2] = Xp[Xnc1]::Array{Float64,1} 
    end

    # calculate new averages
    aak, lak, ldk = Xupd_linavg(k, wck, Xp, Yc, narea)

    if upx == 1  # if root
      llr = (Rupd_llf(k, wck, Xp, Yc, lak, ldk, ωxc, ω1c, ω0c, λc, σ²c) - 
             Rupd_llf(k, wck, Xc, Yc, linavg[k,wck], lindiff[k,wck,:], 
                      ωxc, ω1c, ω0c, λc, σ²c))::Float64 
    else

      wckm1 = wcol[k-1]::Array{Int64,1}
      lakm1 = linavg[(k-1),wckm1]::Array{Float64,1}

      llr = (Xupd_llf(k, wck, wckm1, Xp, Yc, lak, lakm1, 
                      ldk, ωxc, ω1c, ω0c, λc, σ²c) - 
             Xupd_llf(k, wck, wckm1, Xc, Yc, linavg[k,wck], lakm1,
                      lindiff[k,wck,:], ωxc, ω1c, ω0c, λc, σ²c))::Float64 

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

    return Xc, llc, areavg, linavg, lindiff
  end

end





"""
    mhr_upd_σ²(σ²c::Float64, Xc::Array{Float64,2}, ωxc::Float64, llc::Float64, prc::Float64, σ²tn::Float64, linavg::Array{Float64,2}, σ²prior::Float64, σ²ωxupd_llf)


MHR update for σ².
"""
function mhr_upd_σ²(σ²c        ::Float64,
                    Xc         ::Array{Float64,2},
                    ωxc        ::Float64,
                    llc        ::Float64,
                    prc        ::Float64,
                    σ²tn       ::Float64,
                    linavg     ::Array{Float64,2},
                    σ²prior    ::Float64,
                    σ²ωxupd_llf::Function)

  σ²p = logupt(σ²c, rand() < 0.5 ? σ²tn : 4*σ²tn)::Float64

  #likelihood ratio
  llr = (σ²ωxupd_llf(Xc, linavg, ωxc, σ²p) - 
         σ²ωxupd_llf(Xc, linavg, ωxc, σ²c))::Float64

  # prior ratio
  prr = (logdexp(σ²p, σ²prior) - logdexp(σ²c, σ²prior))::Float64

  if log(rand()) < (llr + 
                    prr + 
                    log(σ²p) - log(σ²c))
    llc += llr::Float64
    prc += prr::Float64
    σ²c  = σ²p::Float64
  end

  return (llc, prc, σ²c)::Tuple{Float64,Float64,Float64}
end





"""
    mhr_upd_ωx(ωxc::Float64, Xc::Array{Float64,2}, σ²c::Float64, llc::Float64, prc::Float64, ωxtn::Float64, linavg::Array{Float64,2}, ωxprior::Tuple{Float64,Float64}, σ²ωxupd_llf)

MHR update for ωx.
"""
function mhr_upd_ωx(ωxc         ::Float64,
                    Xc         ::Array{Float64,2},
                    σ²c        ::Float64,
                    llc        ::Float64,
                    prc        ::Float64,
                    ωxtn       ::Float64,
                    linavg     ::Array{Float64,2},
                    ωxprior    ::Tuple{Float64,Float64},
                    σ²ωxupd_llf::Function)

  ωxp = addupt(ωxc, rand() < 0.5 ? ωxtn : 4*ωxtn)::Float64

  #likelihood ratio
  llr = (σ²ωxupd_llf(Xc, linavg, ωxp, σ²c) - 
         σ²ωxupd_llf(Xc, linavg, ωxc, σ²c))::Float64

  # prior ratio
  prr = (logdnorm(ωxp, ωxprior[1], ωxprior[2]) -
         logdnorm(ωxc, ωxprior[1], ωxprior[2]))::Float64

  if log(rand()) < (llr + prr)
    llc += llr::Float64
    prc += prr::Float64
    ωxc  = ωxp::Float64
  end

  return (llc, prc, ωxc)::Tuple{Float64,Float64,Float64}
end





"""
    mhr_upd_ωλ(ω1c::Float64, λc::Array{Float64,2}, ω0c::Float64, Yc::Array{Int64,3}, llc::Float64, prc::Float64, ωλtn::Float64, linavg::Array{Float64,2}, lindiff::Array{Float64,3}, ωλprior::Tuple{Float64,Float64}, ωλμupd_llf)

MHR update for ωλ.
"""
function mhr_upd_ωλ(ω1c       ::Float64,
                    λc        ::Array{Float64,1},
                    ω0c       ::Float64,
                    Yc        ::Array{Int64,3},
                    llc       ::Float64,
                    prc       ::Float64,
                    ωλtn      ::Float64,
                    linavg    ::Array{Float64,2},
                    lindiff   ::Array{Float64,3},
                    ωλprior   ::Tuple{Float64,Float64},
                    ωλμupd_llf::Function)

  ω1p = addupt(ω1c, rand() < 0.5 ? ωλtn : 4*ωλtn)::Float64

  # likelihood ratio
  llr = (ωλμupd_llf(Yc, λc, ω1p, ω0c, lindiff) - 
         ωλμupd_llf(Yc, λc, ω1c, ω0c, lindiff))::Float64

  # prior ratio
  prr = (logdnorm(ω1p, ωλprior[1], ωλprior[2]) -
         logdnorm(ω1c, ωλprior[1], ωλprior[2]))::Float64

  if log(rand()) < (llr + prr)
    llc += llr::Float64
    prc += prr::Float64
    ω1c  = ω1p::Float64
  end

  return (llc, prc, ω1c)::Tuple{Float64,Float64,Float64}
end





"""
    mhr_upd_ωμ(ω0c::Float64, λc::Array{Float64,2}, ω1c::Float64, Yc::Array{Int64,3}, llc::Float64, prc::Float64, ωμtn::Float64, linavg::Array{Float64,2}, lindiff::Array{Float64,3}, ωμprior::Tuple{Float64,Float64}, ωλμupd_llf)

MHR update for ωμ.
"""
function mhr_upd_ωμ(ω0c       ::Float64,
                    λc        ::Array{Float64,1},
                    ω1c       ::Float64,
                    Yc        ::Array{Int64,3},
                    llc       ::Float64,
                    prc       ::Float64,
                    ωμtn      ::Float64,
                    linavg    ::Array{Float64,2},
                    lindiff   ::Array{Float64,3},
                    ωμprior   ::Tuple{Float64,Float64},
                    ωλμupd_llf::Function)

  ω0p = addupt(ω0c, rand() < 0.5 ? ωμtn : 4*ωμtn)::Float64

  # likelihood ratio
  llr = (ωλμupd_llf(Yc, λc, ω1c, ω0p, lindiff) - 
         ωλμupd_llf(Yc, λc, ω1c, ω0c, lindiff))::Float64

  # prior ratio
  prr = (logdnorm(ω0p, ωμprior[1], ωμprior[2]) -
         logdnorm(ω0c, ωμprior[1], ωμprior[2]))::Float64

  if log(rand()) < (llr + prr)
    llc += llr::Float64
    prc += prr::Float64
    ω0c  = ω0p::Float64
  end

  return (llc, prc, ω0c)::Tuple{Float64,Float64,Float64}
end


