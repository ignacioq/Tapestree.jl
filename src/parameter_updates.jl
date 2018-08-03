#=

Updates for MCMC

Ignacio Quintero Mächler

t(-_-t)

June 20 2017

=#





"""
    make_mhr_upd_X(Xnc1     ::Array{Int64,1},
                   Xnc2     ::Array{Int64,1},
                   wcol     ::Array{Array{Int64,1},1},
                   ptn      ::Array{Float64,1},
                   wXp      ::Array{Int64,1},
                   narea    ::Int64,
                   ntip     ::Int64,
                   Xupd_llr ::Function,
                   Rupd_llr ::Function)

Make DA update X.
"""
function make_mhr_upd_X(Xnc1     ::Array{Int64,1},
                        Xnc2     ::Array{Int64,1},
                        wcol     ::Array{Array{Int64,1},1},
                        ptn      ::Array{Float64,1},
                        wXp      ::Array{Int64,1},
                        m        ::Int64,
                        narea    ::Int64,
                        ntip     ::Int64,
                        Xupd_llr ::Function,
                        Rupd_llr ::Function)

  const rj = ind2sub((m,ntip), Xnc2[findfirst(Xnc1, 1)])[2]

  const xpi  = fill(NaN, ntip)
  const δxi  = fill(NaN, ntip, ntip)
  const lapi = fill(NaN, ntip)
  const lani = fill(NaN, ntip)
  const ldi  = fill(NaN, ntip, narea)

  function f(up ::Int64,
             Xc ::Array{Float64,2},
             Yc ::Array{Int64,3},
             δXc::Array{Float64,3},
             δYc::Array{Float64,3},
             λ1c::Float64,
             λ0c::Float64,
             ωxc::Float64, 
             ω1c::Float64, 
             ω0c::Float64,
             σ²c::Float64,
             llc::Float64,
             LAc::Array{Float64,2},
             LDc::Array{Float64,3})

    @inbounds begin

      upx = wXp[up - 6]::Int64                 # X indexing

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

        if -randexp() < llr
          llc    += llr::Float64
          Xc[1,:] = xpi::Array{Float64,1}
        end

      else

        xi, xj = ind2sub(Xc, upx)

        # allocate
        for j = Base.OneTo(ntip)
          xpi[j] =  Xc[xi,j]
          lai[j] = LAc[xi,j]
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
          xpi[ind2sub(Xc, Xnc2[findfirst(Xnc1, upx)])[2]] = xpi[xj]::Float64
        end

        # calculate new averages
        Xupd_linavg!(δxi, lai, ldi, wcol, xpi, xi, xj, Yc, δYc, narea)

        llr = Xupd_llr(xi, xpi, Xc, lai, ldi, LAc, LDc, Yc, 
                       ωxc, ω1c, ω0c, λ1c, λ0c, σ²c)::Float64

        if -randexp() < llr
          llc        += llr::Float64
          Xc[xi,:]    = xpi::Array{Float64,1}
          δXc[:,:,xi] = δxi::Array{Float64,2}
          LAc[xi,:]   = lai::Array{Float64,1}
          LDc[xi,:,:] = ldi::Array{Float64,2}
        end
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
                          total_llf          ::Function)

  const Xp  = zeros(m, ntip)
  const δXp = fill( NaN, ntip, ntip, m)
  const LAp = zeros(m, ntip)
  const LDp = zeros(m, ntip, narea)

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
             LAc    ::Array{Float64,2},
             LDc    ::Array{Float64,3},
             δXc    ::Array{Float64,3},
             δYc    ::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Xp, Xc)

    upbranchX!(br, Xp, bridx, brδt, σ²c)

    deltaX!(δXp, Xp, wcol, m, ntip, narea)
    sde!(LAp, δXp, δYc, wcol, m, ntip)
    lindiff!(LDp, δXp, Yc, wcol, m, ntip, narea)

    llr = (total_llf(Xp, Yc, LAp, LDp, ωxc, ω1c, ω0c, λ1c, λ0c,
                     stemevc, brs, σ²c) - 
           total_llf(Xc, Yc, LAc, LDc, ωxc, ω1c, ω0c, λ1c, λ0c,
                     stemevc, brs, σ²c))::Float64

    if -randexp() < (llr + llr_bm(Xc, Xp, bridx[br], brδt[br], σ²c))::Float64
      llc += llr::Float64
      copy!(Xc,   Xp)
      copy!(δXc, δXp)
      copy!(LAc, LAp)
      copy!(LDc, LDp)
    end

    return llc::Float64
  end

end





"""
    make_mhr_upd_Xtrio(Xnc1::Array{Int64,1}, Xnc2::Array{Int64,1}, wcol::Array{Array{Int64,1},1}, m::Int64, ptn::Array{Float64,1}, wXp::Array{Int64,1}, λlessthan::Int64, narea::Int64, Xupd_llr, Rupd_llr)

Make X trio DA update for a node and adjacent branches using BB proposals.
"""
function make_mhr_upd_Xtrio(wcol               ::Array{Array{Int64,1},1},
                            m                  ::Int64,
                            narea              ::Int64,
                            ntip               ::Int64,
                            nedge              ::Int64,
                            bridx              ::Array{UnitRange{Int64},1},
                            brδt               ::Array{Array{Float64,1},1},
                            total_llf          ::Function)

  const Xp  = zeros(m, ntip)
  const δXp = fill( NaN, ntip, ntip, m)
  const LAp = zeros(m, ntip)
  const LDp = zeros(m, ntip, narea)

  function f(trio   ::Array{Int64,1},
             Xc     ::Array{Float64,2},
             Yc     ::Array{Int64,3},
             λ1c    ::Float64,
             λ0c    ::Float64,
             ωxc    ::Float64, 
             ω1c    ::Float64, 
             ω0c    ::Float64,
             σ²c    ::Float64,
             llc    ::Float64,
             LAc    ::Array{Float64,2},
             LDc    ::Array{Float64,3},
             δXc    ::Array{Float64,3},
             δYc    ::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Xp, Xc)

    pr, d1, d2 = trio

    uptrioX!(pr, d1, d2, Xp, bridx, brδt, σ²c, nedge)

    deltaX!(δXp, Xp, wcol, m, ntip, narea)
    sde!(LAp, δXp, δYc, wcol,  m, ntip)
    lindiff!(LDp, δXp, Yc, wcol, m, ntip, narea)

    llr = (total_llf(Xp, Yc, LAp, LDp, ωxc, ω1c, ω0c, λ1c, λ0c,
                     stemevc, brs, σ²c) - 
           total_llf(Xc, Yc, LAc, LDc, ωxc, ω1c, ω0c, λ1c, λ0c,
                     stemevc, brs, σ²c))::Float64

    if -randexp() < (llr + 
                     (pr != nedge ? 
                      llr_bm(Xc, Xp, bridx[pr], brδt[pr], σ²c):0.0) +
                      llr_bm(Xc, Xp, bridx[d1], brδt[d1], σ²c) +
                      llr_bm(Xc, Xp, bridx[d2], brδt[d2], σ²c))::Float64
      llc += llr::Float64
      copy!(Xc,   Xp)
      copy!(δXc, δXp)
      copy!(LAc, LAp)
      copy!(LDc, LDp)
    end

    return llc::Float64
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
  const δYp     = fill( NaN, ntip, ntip, m)
  const LAp     = zeros(m, ntip)
  const LDp     = zeros(m, ntip, narea)
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
             LAc    ::Array{Float64,2},
             LDc    ::Array{Float64,3},
             δXc    ::Array{Float64,3},
             δYc    ::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Yp, Yc)

    # if a successful sample
    if upbranchY!(λϕ1, λϕ0, br, Yp, stemevp, 
                  bridx_a, brδt, brl[nedge], brs, narea, nedge)

      deltaY!(δYp, Yp, wcol, m, ntip, narea)
      sde!(LAp, δXc, δYp, wcol, m, ntip)
      lindiff!(LDp, δXc, Yp, wcol, m, ntip, narea)

      llr = (total_llf(Xc, Yp, LAp, LDp, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevp, brs, σ²c) - 
             total_llf(Xc, Yc, LAc, LDc, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevc, brs, σ²c))::Float64

      if -randexp() < (llr + 
                       bgiid_br(Yc, stemevc, brs, br, λϕ1, λϕ0) - 
                       bgiid_br(Yp, stemevp, brs, br, λϕ1, λϕ0))
        llc  += llr::Float64
        copy!(Yc,  Yp)
        copy!(δYc, δYp)
        copy!(LAc, LAp)
        copy!(LDc, LDp)
        # allocate stemevc
        @simd for k in Base.OneTo(narea) 
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
    make_mhr_upd_Ytrio(narea    ::Int64,
                       nedge    ::Int64,
                       m        ::Int64,
                       ntip     ::Int64,
                       bridx_a  ::Array{Array{UnitRange{Int64},1},1},
                       brδt     ::Array{Array{Float64,1},1},
                       brl      ::Array{Float64,1},
                       wcol     ::Array{Array{Int64,1},1},
                       total_llf::Function,
                       bgiid    ::Function)

Make function to update trio in Y.
"""
function make_mhr_upd_Ytrio(narea    ::Int64,
                            nedge    ::Int64,
                            m        ::Int64,
                            ntip     ::Int64,
                            bridx_a  ::Array{Array{UnitRange{Int64},1},1},
                            brδt     ::Array{Array{Float64,1},1},
                            brl      ::Array{Float64,1},
                            wcol     ::Array{Array{Int64,1},1},
                            total_llf::Function,
                            bgiid    ::Function)

  const Yp      = zeros(Int64, m, ntip, narea)
  const δYp     = fill( NaN, ntip, ntip, m)
  const LAp     = zeros(m, ntip)
  const LDp     = zeros(m, ntip, narea)
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
             LAc    ::Array{Float64,2},
             LDc    ::Array{Float64,3},
             δXc    ::Array{Float64,3},
             δYc    ::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Yp,   Yc)
    copy!(brsp, brs)

    # if an efficient sample
    if upnode!(λϕ1, λϕ0, triad, Yp, stemevp,
               bridx_a, brδt, brl, brsp, narea, nedge)

      deltaY!(δYp, Yp, wcol, m, ntip, narea)
      sde!(LAp, δXc, δYp, wcol, m, ntip)
      lindiff!(LDp, δXc, Yp, wcol, m, ntip, narea)

      llr = (total_llf(Xc, Yp, LAp, LDp, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevp, brsp, σ²c) - 
             total_llf(Xc, Yc, LAc, LDc, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevc, brs, σ²c))::Float64

      if -randexp() < (llr + 
                       bgiid(Yc, stemevc, brs,  triad, λϕ1, λϕ0) - 
                       bgiid(Yp, stemevp, brsp, triad, λϕ1, λϕ0))
        llc += llr
        copy!(Yc,  Yp)
        copy!(δYc, δYp)
        copy!(LAc, LAp)
        copy!(LDc, LDp)
        copy!(brs, brsp)
        # allocate stemevc
        @simd for k in Base.OneTo(narea) 
          stemevc[k] = copy(stemevp[k])
        end
      end

      return llc::Float64
    else
      return llc::Float64
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
                           bgiid_br           ::Function)

  const Xp  = zeros(m, ntip)
  const Yp  = zeros(Int64, m, ntip, narea)
  const δXp = fill( NaN, ntip, ntip, m)
  const δYp = fill( NaN, ntip, ntip, m)
  const LAp = zeros(m, ntip)
  const LDp = zeros(m, ntip, narea)

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
             LAc    ::Array{Float64,2},
             LDc    ::Array{Float64,3},
             δXc    ::Array{Float64,3},
             δYc    ::Array{Float64,3},
             brs    ::Array{Int64,3},
             stemevc::Array{Array{Float64,1},1})

    copy!(Xp, Xc)
    copy!(Yp, Yc)

  # if an successful sample
  if upbranchY!(λϕ1, λϕ0, br, Yp, stemevc, 
                  bridx_a, brδt, brl[nedge], brs, narea, nedge)

      upbranchX!(br, Xp, bridx, brδt, σ²c)

      deltaXY!(δXp, δYp, Xp, Yp, wcol, m, ntip, narea)
      sde!(LAp, δXp, δYp, wcol, m, ntip)
      lindiff!(LDp, δXp, Yp, wcol, m, ntip, narea)

      llr = (total_llf(Xp, Yp, LAp, LDp, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevc, brs, σ²c) - 
             total_llf(Xc, Yc, LAc, LDc, ωxc, ω1c, ω0c, λ1c, λ0c,
                       stemevc, brs, σ²c))::Float64

      if -randexp() < (llr + 
                       bgiid_br(Yc, stemevc, brs, br, λϕ1, λϕ0) - 
                       bgiid_br(Yp, stemevc, brs, br, λϕ1, λϕ0) +
                       llr_bm(Xc, Xp, bridx[br], brδt[br], σ²c))::Float64
        llc += llr::Float64
        copy!(Xc,  Xp)
        copy!(Yc,  Yp)
        copy!(δXc, δXp)
        copy!(δYc, δYp)
        copy!(LAc, LAp)
        copy!(LDc, LDp)
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

    # update stem node and branch
    if upstemnode!(λϕ1, λϕ0, nedge, stemevp, brsp, stbrl, narea)

      llr = stem_llr(λ1c, λ0c, brs, brsp, stemevc, stemevp, narea, nedge)

      # likelihood ratio
      if -randexp() < (llr + stemiid_propr(λϕ1, λϕ0, brs, brsp, 
                                           stemevc, stemevp, narea, nedge))
        llc += llr::Float64
        copy!(brs, brsp)
        @simd for k in Base.OneTo(narea) 
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
                    LAc        ::Array{Float64,2},
                    σ²prior    ::Float64,
                    σ²ωxupd_llr::Function)

  const σ²p = mulupt(σ²c, rand() < .3 ? σ²tn : 4.*σ²tn)::Float64

  #likelihood ratio
  const llr = σ²ωxupd_llr(Xc, LAc, ωxc, ωxc, σ²c, σ²p)::Float64

  # prior ratio
  const prr = (logdexp(σ²p, σ²prior) - logdexp(σ²c, σ²prior))::Float64

  if -randexp() < (llr + prr + Base.Math.JuliaLibm.log(σ²p) - 
                               Base.Math.JuliaLibm.log(σ²c))
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
function mhr_upd_ωx(ωxc        ::Float64,
                    Xc         ::Array{Float64,2},
                    σ²c        ::Float64,
                    llc        ::Float64,
                    prc        ::Float64,
                    ωxtn       ::Float64,
                    LAc        ::Array{Float64,2},
                    ωxprior    ::Tuple{Float64,Float64},
                    σ²ωxupd_llr::Function)

  const ωxp = addupt(ωxc, rand() < .3 ? ωxtn : 4.*ωxtn)::Float64

  #likelihood ratio
  const llr = σ²ωxupd_llr(Xc, LAc, ωxc, ωxp, σ²c, σ²c)::Float64

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
                    LDc       ::Array{Float64,3},
                    ω1prior   ::Tuple{Float64,Float64},
                    ω10upd_llr::Function)

  const ω1p = addupt(ω1c, rand() < .3 ? ω1tn : 4.*ω1tn)::Float64

  # likelihood ratio
  const llr = ω10upd_llr(Yc, λ1c, λ0c, ω1c, ω0c, ω1p, ω0c, LDc)::Float64

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
                    LDc       ::Array{Float64,3},
                    ω0prior   ::Tuple{Float64,Float64},
                    ω10upd_llr::Function)

  const ω0p = addupt(ω0c, rand() < .3 ? ω0tn : 4.*ω0tn)::Float64

  # likelihood ratio
  const llr = ω10upd_llr(Yc, λ1c, λ0c, ω1c, ω0c, ω1c, ω0p, LDc)::Float64

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
                    LDc     ::Array{Float64,3},
                    stemevc ::Array{Array{Float64,1},1},
                    brs     ::Array{Int64,3},
                    λprior  ::Float64,
                    λ1tn    ::Float64,
                    λupd_llr::Function)

  # update λ
  const λ1p = mulupt(λ1c, rand() < .3 ? λ1tn : 4.*λ1tn)::Float64

  # proposal likelihood and prior
  const llr = λupd_llr(Yc, λ1c, λ0c, λ1p, λ0c, ω1c, ω0c, LDc, stemevc, brs)

  const prr = (logdexp(λ1p, λprior) - 
               logdexp(λ1c, λprior))::Float64

  if -randexp() < (llr + prr + Base.Math.JuliaLibm.log(λ1p) - 
                               Base.Math.JuliaLibm.log(λ1c))
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
                    LDc     ::Array{Float64,3},
                    stemevc ::Array{Array{Float64,1},1},
                    brs     ::Array{Int64,3},
                    λprior  ::Float64,
                    λ0tn    ::Float64,
                    λupd_llr::Function)

  # update λ
  const λ0p = mulupt(λ0c, rand() < .3 ? λ0tn : 4.*λ0tn)::Float64

  # proposal likelihood and prior
  const llr = λupd_llr(Yc, λ1c, λ0c, λ1c, λ0p, ω1c, ω0c, LDc, stemevc, brs)

  const prr = (logdexp(λ0p, λprior) - logdexp(λ0c, λprior))::Float64

  if -randexp() < (llr + prr + Base.Math.JuliaLibm.log(λ0p) - 
                               Base.Math.JuliaLibm.log(λ0c))
    llc += llr::Float64
    prc += prr::Float64
    λ0c  = λ0p::Float64
  end

  return (llc, prc, λ0c)::Tuple{Float64,Float64,Float64}
end




