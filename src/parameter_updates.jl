#=

Updates for MCMC

Ignacio Quintero Mächler

t(-_-t)

June 20 2017

=#





"""
    make_mhr_upd_X(Xnc1::Array{Int64,1}, Xnc2::Array{Int64,1}, wcol::Array{Array{Int64,1},1}, m::Int64, ptn::Array{Float64,1}, wXp::Array{Int64,1}, λlessthan::Int64, narea::Int64, Xupd_llr, Rupd_llr)

Make DA update X.
"""
function make_mhr_upd_X(Xnc1     ::Array{Int64,1},
                        Xnc2     ::Array{Int64,1},
                        wcol     ::Array{Array{Int64,1},1},
                        m        ::Int64,
                        ptn      ::Array{Float64,1},
                        wXp      ::Array{Int64,1},
                        narea    ::Int64,
                        ntip     ::Int64,
                        Xupd_llr ::Function,
                        Rupd_llr ::Function)

  const aai = fill(NaN, narea)
  const lai = fill(NaN, ntip)
  const ldi = fill(NaN, ntip, narea)

  function f(up     ::Int64,
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λ1c    ::Float64,
             λ0c    ::Float64,
             ωxc    ::Float64, 
             ω1c    ::Float64, 
             ω0c    ::Float64,
             σ²c    ::Float64,
             llc    ::Float64,
             areavg ::Array{Float64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             areaoc ::Array{Int64,2})

    @inbounds begin
      upx = wXp[up - 6]::Int64                 # X indexing

      xi, xj = ind2sub(Xc, upx)

      xpi = Xc[xi,:]::Array{Float64,1}

      xpi[xj] = addupt(xpi[xj], ptn[up])::Float64      # update X

      if in(upx, Xnc1)        # if an internal node
        xpi[ind2sub(Xc, Xnc2[findfirst(Xnc1, upx)])[2]] = xpi[xj]::Float64
      end

      # calculate new averages
      Xupd_linavg!(aai, lai, ldi, areaoc, xi, wcol[xi], xpi, Yc, narea)

      if upx == 1  # if root
        llr = Rupd_llr(wcol[1], 
                       xpi[wcol[1]], 
                       Xc[1,wcol[1]], Xc[2,wcol[1]], 
                       lai[wcol[1]], ldi[wcol[1],:], 
                       linavg[1,wcol[1]], lindiff[1,wcol[1],:],
                       Yc, 
                       ωxc, ω1c, ω0c, λ1c, λ0c, σ²c)::Float64
      else
        llr = Xupd_llr(xi, wcol[xi], wcol[xi-1], 
                       xpi, 
                       Xc[xi,:], Xc[xi-1,:], Xc[xi+1,:], 
                       lai, ldi, 
                       linavg[xi,:], linavg[xi-1,:], 
                       lindiff[xi,:,:],
                       Yc, 
                       ωxc, ω1c, ω0c, λ1c, λ0c, σ²c)::Float64
      end

      if -randexp() < llr
        llc            += llr::Float64
        Xc[xi,:]        = xpi::Array{Float64,1}
        areavg[xi,:]    = aai::Array{Float64,1}
        linavg[xi,:]    = lai::Array{Float64,1}
        lindiff[xi,:,:] = ldi::Array{Float64,2}
      end
    end

    return llc::Float64
  end
end





"""
    make_mhr_upd_Xbr(Xnc1::Array{Int64,1}, Xnc2::Array{Int64,1}, wcol::Array{Array{Int64,1},1}, m::Int64, ptn::Array{Float64,1}, wXp::Array{Int64,1}, λlessthan::Int64, narea::Int64, Xupd_llr, Rupd_llr)

Make X branch DA update for a single branch using BB proposals.
"""
function make_mhr_upd_Xbr(wcol               ::Array{Array{Int64,1},1},
                          m                  ::Int64,
                          narea              ::Int64,
                          ntip               ::Int64,
                          nedge              ::Int64,
                          bridx              ::Array{UnitRange{Int64},1},
                          brδt               ::Array{Array{Float64,1},1},
                          total_llf          ::Function,
                          σ²prior            ::Float64)

  const Xp = zeros(m, ntip)
  const aa = zeros(m, narea)
  const la = zeros(m, ntip)
  const ld = zeros(m, ntip, narea)

  function f(br     ::Int64,
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λ1c    ::Float64,
             λ0c    ::Float64,
             ωxc    ::Float64, 
             ω1c    ::Float64, 
             ω0c    ::Float64,
             σ²c    ::Float64,
             llc    ::Float64,
             areavg ::Array{Float64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             areaoc ::Array{Int64,2},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Xp, Xc)
    copy!(aa, areavg)
    copy!(la, linavg)
    copy!(ld, lindiff)

    upbranchX!(br, Xp, bridx, brδt, σ²c)

    area_lineage_means!(aa, la, areaoc, Xp, Yc, wcol, m, narea)
    linarea_diff!(ld, Xp, aa, areaoc, narea, ntip, m)

    llr = (total_llf(Xp, Yc, la, ld, ωxc, ω1c, ω0c, λ1c, λ0c,
                     stemevc, brs[nedge,1,:], σ²c) - 
           total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λ1c, λ0c,
                     stemevc, brs[nedge,1,:], σ²c))::Float64

    if -randexp() < (llr + llr_bm(Xc, Xp, bridx[br], brδt[br], σ²c))::Float64
      llc += llr::Float64
      copy!(Xc,      Xp)
      copy!(areavg,  aa)
      copy!(linavg,  la)
      copy!(lindiff, ld)
    end

    return llc::Float64
  end

end





"""
    make_mhr_upd_Y(narea::Int64, nedge::Int64, m::Int64, ntip::Int64, bridx_a::Vector{Vector{Vector{Int64}}}, brδt::Array{Array{Float64,1},1}, brl::Array{Float64,1}, wcol::Array{Array{Int64,1},1}, Ync1::Array{Int64,1}, Ync2::Array{Int64,1}, total_llf, biogeo_upd_iid)

Make function to update trio in Y.
"""
function make_mhr_upd_Y(narea              ::Int64,
                        nedge              ::Int64,
                        m                  ::Int64,
                        ntip               ::Int64,
                        bridx_a            ::Array{Array{UnitRange{Int64},1},1},
                        brδt               ::Array{Array{Float64,1},1},
                        brl                ::Array{Float64,1},
                        wcol               ::Array{Array{Int64,1},1},
                        total_llf          ::Function,
                        bgiid              ::Function)

  const Yp      = zeros(Int64, m, ntip, narea)
  const aa      = zeros(m, narea)
  const ao      = zeros(Int64, m, narea)
  const la      = zeros(m, ntip)
  const ld      = zeros(m, ntip, narea)
  const brsp    = zeros(Int64, nedge, 2, narea)
  const stemevp = [[rand()] for i in 1:narea]

  function f(triad  ::Array{Int64,1},
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λ1c    ::Float64,
             λ0c    ::Float64,
             ωxc    ::Float64, 
             ω1c    ::Float64, 
             ω0c    ::Float64,
             σ²c    ::Float64,
             λϕ1   ::Float64,
             λϕ0   ::Float64,
             llc    ::Float64,
             prc    ::Float64,
             areavg ::Array{Float64,2},
             areaoc ::Array{Int64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Yp,      Yc)
    copy!(aa,      areavg)
    copy!(ao,      areaoc)
    copy!(la,      linavg)
    copy!(ld,      lindiff)
    copy!(brsp,    brs)

    # allocate stemevp
    for k in Base.OneTo(narea) 
      stemevp[k] = copy(stemevc[k])
    end

    # if an efficient sample
    if upnode!(λϕ1, λϕ0, triad, Yp, stemevp,
               bridx_a, brδt, brl, brsp, narea, nedge) === nothing

      area_lineage_means!(aa, la, ao, Xc, Yp, wcol, m, narea)
      linarea_diff!(ld, Xc, aa, ao, narea, ntip, m)

      llr = (total_llf(Xc, Yp, la, ld, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevp, brsp[nedge,1,:], σ²c) - 
             total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevc, brs[nedge,1,:], σ²c))::Float64

      if -randexp() < (llr + 
                       bgiid(Yc, stemevc, brs[nedge,1,:],  triad, λϕ1, λϕ0) - 
                       bgiid(Yp, stemevp, brsp[nedge,1,:], triad, λϕ1, λϕ0))

        llc += llr
        copy!(Yc,           Yp)
        copy!(areavg,       aa)
        copy!(areaoc,       ao)
        copy!(linavg,       la)
        copy!(lindiff,      ld)
        copy!(brs,        brsp)
        # allocate stemevc
        for k in Base.OneTo(narea) 
          stemevc[k] = copy(stemevp[k])
        end
      end

      return llc
    else
      return llc
    end
  end

end





"""
    make_mhr_upd_Ybr(narea::Int64, nedge::Int64, m::Int64, ntip::Int64, bridx_a::Vector{Vector{Vector{Int64}}}, brδt::Array{Array{Float64,1},1}, brl::Array{Float64,1}, wcol::Array{Array{Int64,1},1}, Ync1::Array{Int64,1}, Ync2::Array{Int64,1}, total_llf, biogeo_upd_iid)

Make function to update a single branch in Y.
"""
function make_mhr_upd_Ybr(narea              ::Int64,
                          nedge              ::Int64,
                          m                  ::Int64,
                          ntip               ::Int64,
                          bridx_a            ::Array{Array{UnitRange{Int64},1},1},
                          brδt               ::Array{Array{Float64,1},1},
                          brl                ::Array{Float64,1},
                          wcol               ::Array{Array{Int64,1},1},
                          total_llf          ::Function,
                          bgiid_br           ::Function)

  const Yp      = zeros(Int64, m, ntip, narea)
  const aa      = zeros(m, narea)
  const ao      = zeros(Int64,m, narea)
  const la      = zeros(m, ntip)
  const ld      = zeros(m, ntip, narea)
  const stemevp = [[rand()] for i in 1:narea]

  function f(br     ::Int64,
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λ1c    ::Float64,
             λ0c    ::Float64,
             ωxc    ::Float64, 
             ω1c    ::Float64, 
             ω0c    ::Float64,
             σ²c    ::Float64,
             λϕ1   ::Float64,
             λϕ0   ::Float64,
             llc    ::Float64,
             prc    ::Float64,
             areavg ::Array{Float64,2},
             areaoc ::Array{Int64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Yp,           Yc)
    copy!(aa,       areavg)
    copy!(ao,       areaoc)
    copy!(la,       linavg)
    copy!(ld,      lindiff)
    # allocate stemevc
    for k in Base.OneTo(narea) 
      stemevp[k] = copy(stemevc[k])
    end

    # if an efficient sample
    if upbranchY!(λϕ1, λϕ0, br, Yp, stemevp, 
                  bridx_a, brδt, brl[nedge], brs, narea, nedge) === nothing

      area_lineage_means!(aa, la, ao, Xc, Yp, wcol, m, narea)
      linarea_diff!(ld, Xc, aa, ao, narea, ntip, m)

      llr = (total_llf(Xc, Yp, la, ld, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevp, brs[nedge,1,:], σ²c) - 
             total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevc, brs[nedge,1,:], σ²c))::Float64

      if -randexp() < (llr + 
                       bgiid_br(Yc, stemevc, brs[nedge,1,:], br, λϕ1, λϕ0) - 
                       bgiid_br(Yp, stemevp, brs[nedge,1,:], br, λϕ1, λϕ0))
        llc  += llr::Float64
        copy!(Yc,      Yp)
        copy!(areavg,  aa)
        copy!(areaoc,  ao)
        copy!(linavg,  la)
        copy!(lindiff, ld)
        # allocate stemevc
        for k in Base.OneTo(narea) 
          stemevc[k] = copy(stemevp[k])
        end      
      end
      return llc
    else
      return llc
    end
  end

end





"""
    make_mhr_upd_Ybr(narea::Int64, nedge::Int64, m::Int64, ntip::Int64, bridx_a::Vector{Vector{Vector{Int64}}}, brδt::Array{Array{Float64,1},1}, brl::Array{Float64,1}, wcol::Array{Array{Int64,1},1}, Ync1::Array{Int64,1}, Ync2::Array{Int64,1}, total_llf, biogeo_upd_iid)

Make function to simultaneously update a single branch in `X` & `Y`.
"""
function make_mhr_upd_XYbr(narea              ::Int64,
                           nedge              ::Int64,
                           m                  ::Int64,
                           ntip               ::Int64,
                           bridx              ::Array{UnitRange{Int64},1},
                           bridx_a            ::Array{Array{UnitRange{Int64},1},1},
                           brδt               ::Array{Array{Float64,1},1},
                           brl                ::Array{Float64,1},
                           wcol               ::Array{Array{Int64,1},1},
                           total_llf          ::Function,
                           bgiid_br           ::Function,
                           σ²prior            ::Float64)

  const Xp = zeros(m, ntip)
  const Yp = zeros(Int64, m, ntip, narea)
  const aa = zeros(m, narea)
  const ao = zeros(Int64, m, narea)
  const la = zeros(m, ntip)
  const ld = zeros(m, ntip, narea)

  function f(br     ::Int64,
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λ1c    ::Float64,
             λ0c    ::Float64,
             ωxc    ::Float64,
             ω1c    ::Float64,
             ω0c    ::Float64,
             σ²c    ::Float64,
             λϕ1   ::Float64, 
             λϕ0   ::Float64,
             llc    ::Float64,
             prc    ::Float64,
             areavg ::Array{Float64,2},
             areaoc ::Array{Int64,2},
             linavg ::Array{Float64,2},
             lindiff::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Xp, Xc)
    copy!(Yp, Yc)
    copy!(aa, areavg)
    copy!(ao, areaoc)
    copy!(la, linavg)
    copy!(ld, lindiff)

  # if an efficient sample
    if upbranchY!(λϕ1, λϕ0, br, Yp, stemevc, 
                  bridx_a, brδt, brl[nedge], brs, narea, nedge) === nothing

      upbranchX!(br, Xp, bridx, brδt, σ²c)

      area_lineage_means!(aa, la, ao, Xp, Yp, wcol, m, narea)
      linarea_diff!(ld, Xp, aa, ao, narea, ntip, m)

      llr = total_llf(Xp, Yp, la, ld, ωxc, ω1c, ω0c, λ1c, λ0c,
                      stemevc, brs[nedge,1,:], σ²c) - 
            total_llf(Xc, Yc, linavg, lindiff, ωxc, ω1c, ω0c, λ1c, λ0c,
                      stemevc, brs[nedge,1,:], σ²c)

      if -randexp() < (llr + 
                       bgiid_br(Yc, stemevc, brs[nedge,1,:], br, λϕ1, λϕ0) - 
                       bgiid_br(Yp, stemevc, brs[nedge,1,:], br, λϕ1, λϕ0) +
                       llr_bm(Xc, Xp, bridx[br], brδt[br], σ²c))::Float64
        llc += llr::Float64
        copy!(Xc,      Xp)
        copy!(Yc,      Yp)
        copy!(areavg,  aa)
        copy!(areaoc,  ao)
        copy!(linavg,  la)
        copy!(lindiff, ld)
      end

      return llc::Float64
    else
      return llc::Float64
    end
  end

end





"""
    make_mhr_upd_Ystem(stbrl::Float64,
                       narea::Int64,
                       nedge::Int64)

Perform MH update for stem node and branch.
"""
function make_mhr_upd_Ystem(stbrl::Float64,
                            narea::Int64,
                            nedge::Int64)

  const brsp    = zeros(Int64, nedge, 2, narea)
  const stemevp = [[rand()] for i in 1:narea]

  function f(λ1c    ::Float64,
             λ0c    ::Float64,
             λϕ1   ::Float64, 
             λϕ0   ::Float64,
             llc    ::Float64,
             stemevc::Array{Array{Float64,1},1},
             brs    ::Array{Int64,3})

    copy!(brsp, brs)
    # allocate stemevc
    for k in Base.OneTo(narea) 
      stemevp[k] = copy(stemevc[k])
    end

    # update stem node and branch
    if upstemnode!(λϕ1, λϕ0, nedge, stemevp, brsp, stbrl, narea) === nothing

      llr = stem_llr(λ1c, λ0c, 
                     brs[nedge,1,:], brsp[nedge,1,:], 
                     stemevc, stemevp, narea)

      # likelihood ratio
      if -randexp() < (llr + 
                       stemiid_propr(λϕ1, λϕ0, brs[nedge,1,:], brsp[nedge,1,:], 
                                     stemevc, stemevp, narea))
        
        llc += llr::Float64
        copy!(brs, brsp)
        for k in Base.OneTo(narea) 
          stemevc[k] = copy(stemevp[k])
        end
      end

      return llc::Float64
    else
      return llc::Float64
    end
  end
  
  return f::Function
end





#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Single governing parameter updates
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-





"""
    mhr_upd_σ²(σ²c::Float64, Xc::Array{Float64,2}, ωxc::Float64, llc::Float64, prc::Float64, σ²tn::Float64, linavg::Array{Float64,2}, σ²prior::Float64, σ²ωxupd_llr)


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
                    σ²ωxupd_llr::Function)

  const σ²p = mulupt(σ²c, rand() < .3 ? σ²tn : 4.*σ²tn)::Float64

  #likelihood ratio
  const llr = σ²ωxupd_llr(Xc, linavg, ωxc, ωxc, σ²c, σ²p)::Float64

  # prior ratio
  const prr = (logdexp(σ²p, σ²prior) - 
               logdexp(σ²c, σ²prior))::Float64

  if -randexp() < (llr + prr + log(σ²p) - log(σ²c))
    llc += llr::Float64
    prc += prr::Float64
    σ²c  = σ²p::Float64
  end

  return (llc, prc, σ²c)::Tuple{Float64,Float64,Float64}
end





"""
    mhr_upd_ωx(ωxc::Float64, Xc::Array{Float64,2}, σ²c::Float64, llc::Float64, prc::Float64, ωxtn::Float64, linavg::Array{Float64,2}, ωxprior::Tuple{Float64,Float64}, σ²ωxupd_llr)

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
                    σ²ωxupd_llr::Function)

  const ωxp = addupt(ωxc, rand() < .3 ? ωxtn : 4.*ωxtn)::Float64

  #likelihood ratio
  const llr = σ²ωxupd_llr(Xc, linavg, ωxc, ωxp, σ²c, σ²c)::Float64

  # prior ratio
  const prr = (logdnorm_tc(ωxp, ωxprior[1], ωxprior[2]) -
               logdnorm_tc(ωxc, ωxprior[1], ωxprior[2]))::Float64

  if -randexp() < (llr + prr)
    llc += llr::Float64
    prc += prr::Float64
    ωxc  = ωxp::Float64
  end

  return (llc, prc, ωxc)::Tuple{Float64,Float64,Float64}
end





"""
    mhr_upd_ω1(ω1c::Float64, λc::Array{Float64,2}, ω0c::Float64, Yc::Array{Int64,3}, llc::Float64, prc::Float64, ω1tn::Float64, linavg::Array{Float64,2}, lindiff::Array{Float64,3}, ω1prior::Tuple{Float64,Float64}, ω10upd_llr)

MHR update for ω1.
"""
function mhr_upd_ω1(ω1c       ::Float64,
                    λ1c       ::Float64,
                    λ0c       ::Float64,
                    ω0c       ::Float64,
                    Yc        ::Array{Int64,3},
                    llc       ::Float64,
                    prc       ::Float64,
                    ω1tn      ::Float64,
                    linavg    ::Array{Float64,2},
                    lindiff   ::Array{Float64,3},
                    ω1prior   ::Tuple{Float64,Float64},
                    ω10upd_llr::Function)

  const ω1p = addupt(ω1c, rand() < .3 ? ω1tn : 4.*ω1tn)::Float64

  # likelihood ratio
  const llr = ω10upd_llr(Yc, λ1c, λ0c, ω1c, ω0c, ω1p, ω0c, lindiff)::Float64

  # prior ratio
  const prr = (logdnorm_tc(ω1p, ω1prior[1], ω1prior[2]) -
               logdnorm_tc(ω1c, ω1prior[1], ω1prior[2]))::Float64

  if -randexp() < (llr + prr)
    llc += llr::Float64
    prc += prr::Float64
    ω1c  = ω1p::Float64
  end

  return (llc, prc, ω1c)::Tuple{Float64,Float64,Float64}
end





"""
    mhr_upd_ω0(ω0c::Float64, λc::Array{Float64,2}, ω1c::Float64, Yc::Array{Int64,3}, llc::Float64, prc::Float64, ω0tn::Float64, linavg::Array{Float64,2}, lindiff::Array{Float64,3}, ω0prior::Tuple{Float64,Float64}, ω10upd_llr)

MHR update for ω0.
"""
function mhr_upd_ω0(ω0c       ::Float64,
                    λ1c       ::Float64,
                    λ0c       ::Float64,
                    ω1c       ::Float64,
                    Yc        ::Array{Int64,3},
                    llc       ::Float64,
                    prc       ::Float64,
                    ω0tn      ::Float64,
                    linavg    ::Array{Float64,2},
                    lindiff   ::Array{Float64,3},
                    ω0prior   ::Tuple{Float64,Float64},
                    ω10upd_llr::Function)

  const ω0p = addupt(ω0c, rand() < .3 ? ω0tn : 4.*ω0tn)::Float64

  # likelihood ratio
  const llr = ω10upd_llr(Yc, λ1c, λ0c, ω1c, ω0c, ω1c, ω0p, lindiff)::Float64

  # prior ratio
  const prr = (logdnorm_tc(ω0p, ω0prior[1], ω0prior[2]) -
               logdnorm_tc(ω0c, ω0prior[1], ω0prior[2]))::Float64

  if -randexp() < (llr + prr)
    llc += llr::Float64
    prc += prr::Float64
    ω0c  = ω0p::Float64
  end

  return (llc, prc, ω0c)::Tuple{Float64,Float64,Float64}
end





"""
    mhr_upd_λ1(...)

Update λ1.
"""
function mhr_upd_λ1(λ1c     ::Float64,
                    Yc      ::Array{Int64,3},
                    λ0c     ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    ω1c     ::Float64,
                    ω0c     ::Float64,
                    lindiff ::Array{Float64,3},
                    stemevc ::Array{Array{Float64,1},1},
                    stemss  ::Array{Int64,1},
                    λprior  ::Float64,
                    λ1tn    ::Float64,
                    λupd_llr::Function)

  # update λ
  const λ1p = mulupt(λ1c, rand() < .3 ? λ1tn : 4.*λ1tn)::Float64

  # proposal likelihood and prior
  const llr = λupd_llr(Yc, λ1c, λ0c, λ1p, λ0c, ω1c, ω0c, 
                       lindiff, stemevc, stemss)::Float64

  const prr = (logdexp(λ1p, λprior) - 
               logdexp(λ1c, λprior))::Float64

  if -randexp() < (llr + prr + log(λ1p) - log(λ1c))
    llc += llr::Float64
    prc += prr::Float64
    λ1c  = λ1p::Float64
  end

  return (llc, prc, λ1c)::Tuple{Float64,Float64,Float64}
end





"""
    mhr_upd_λ0(...)

Update λ0.
"""
function mhr_upd_λ0(λ0c     ::Float64,
                    Yc      ::Array{Int64,3},
                    λ1c     ::Float64,
                    llc     ::Float64,
                    prc     ::Float64,
                    ω1c     ::Float64,
                    ω0c     ::Float64,
                    lindiff ::Array{Float64,3},
                    stemevc ::Array{Array{Float64,1},1},
                    stemss  ::Array{Int64,1},
                    λprior  ::Float64,
                    λ0tn    ::Float64,
                    λupd_llr::Function)

  # update λ
  const λ0p = mulupt(λ0c, rand() < .3 ? λ0tn : 4.*λ0tn)::Float64

  # proposal likelihood and prior
  const llr = λupd_llr(Yc, λ1c, λ0c, λ1c, λ0p, ω1c, ω0c, 
                       lindiff, stemevc, stemss)::Float64

  const prr = (logdexp(λ0p, λprior) - 
               logdexp(λ0c, λprior))::Float64

  if -randexp() < (llr + prr + log(λ0p) - log(λ0c))
    llc += llr::Float64
    prc += prr::Float64
    λ0c  = λ0p::Float64
  end

  return (llc, prc, λ0c)::Tuple{Float64,Float64,Float64}
end




