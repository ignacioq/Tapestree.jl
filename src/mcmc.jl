#=

MCMC for biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#







"""
    tribe_mcmc(...)

Run MCMC for trait and range interspecitific biogeographic evolution (tribe).
"""
function tribe_mcmc(Xc      ::Array{Float64,2},
                    Yc      ::Array{Int64,3},
                    ncoup   ::Array{Int64,2},
                    δt      ::Array{Float64,1},
                    edges   ::Array{Int64,2},
                    brl     ::Array{Float64,1},
                    B       ::Array{Float64,2};
                    niter   ::Int64             = 500_000,
                    nthin   ::Int64             = 1_000,
                    nburn   ::Int64             = 500_000,
                    saveXY  ::Tuple{Bool,Int64} = (false, 1_000),
                    saveDM  ::Tuple{Bool,Int64} = (false, 1_000),
                    ωxprior ::NTuple{2,Float64} = (0.,10.),
                    ω1prior ::NTuple{2,Float64} = (0.,10.),
                    ω0prior ::NTuple{2,Float64} = (0.,10.),
                    σ²prior ::Float64           = 1e-1,
                    λprior  ::Float64           = 1e-1,
                    out_file::String            = "tribe_results",
                    weight  ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),
                    σ²i     ::Float64           = 1.,
                    ωxi     ::Float64           = 0.,
                    ω1i     ::Float64           = 0.,
                    ω0i     ::Float64           = 0.,
                    λ1i     ::Float64           = 1.0,
                    λ0i     ::Float64           = 0.2,
                    stbrl   ::Float64           = 1.,
                    fix_ωx  ::Bool              = false,
                    fix_ω1  ::Bool              = false,
                    fix_ω0  ::Bool              = false)

  printstyled("Data successfully processed", bold = true,color=:green)

  # dims
  m, ntip, narea  = size(Yc)

  # coupled nodes for X
  Xnc1 = ncoup[:,1]
  Xnc2 = ncoup[:,2]

  # which nodes are not NaN in Xc
  wXp = setdiff(LinearIndices(Xc)[findall(!isnan, Xc)], m:m:length(Xc))

  # tie trait coupled nodes
  Xc[Xnc2] = Xc[Xnc1]

  #create object with column indices
  wcol = create_wcol(Xc)

  # add a branch as long as the tree
  stbrl = isone(stbrl) ? sum(δt) : stbrl

  # add long stem branch
  edges = cat(edges, [2*ntip ntip + 1], dims = 1)

  push!(brl, stbrl)

  # number of edges
  nedge = size(edges,1) 

  # make edge triads
  trios = maketriads(edges)

  # make ragged array with index for each edge in Yc
  bridx = make_edgeind(edges[:,2], B, ntip)

  # expand bridx for each area
  bridx_a = Array{UnitRange{Int64},1}[]
  push!(bridx_a, bridx)

  for j = 2:narea
    bridxj = copy(bridx)
    for i in Base.OneTo(nedge)
      bridxj[i] = (bridxj[i]) .+ (j-1)*(m*ntip)
    end
    push!(bridx_a, bridxj)
  end

  # make ragged array of cumulative delta times for each branch
  brct = make_edgeδt(bridx, δt, m)

  # array for states at start and end of branches
  brs = ones(Int64, nedge, 2, narea)

  # filter tips in edges
  wtips = map(x -> x <= ntip, edges[:,2])

  ## assign tip area states to brs according Yc
  for j in Base.OneTo(narea)
    brs[wtips,2,j] = Yc[end,:,j]
  end

  # create stem events
  stemevc = [[rand()] for i in 1:narea]

  # Sample all internal node values according to Pr transitions
  for triad in trios
    λϕ1, λϕ0 = λϕprop()
    while !(upnode!(λϕ1, λϕ0, triad, Yc, stemevc, bridx_a, 
                   brct, brl, brs, narea, nedge))
      λϕ1, λϕ0 = λϕprop()
    end
  end

  ## allocate averages for X and Y
  # X and Y distance matrix
  δXc = fill(NaN, ntip, ntip, m)
  δYc = fill(NaN, ntip, ntip, m)
  deltaXY!(δXc, δYc, Xc, Yc, wcol, m, ntip, narea)

  # lineage averages
  LApc = fill(NaN, m, ntip)
  LAnc = fill(NaN, m, ntip)
  sde!(LApc, LAnc, δXc, δYc, wcol, m, ntip)

  # area lineage distances
  LDc = fill(NaN, m, ntip, narea)
  lindiff!(LDc, δXc, Yc, wcol, m, ntip, narea)

  ## make likelihood and prior functions
  total_llf, total_llr = makellf(δt, Yc, ntip, narea, m, nedge)
  ω1upd_llr, ω0upd_llr, λ1upd_llr, λ0upd_llr = 
    makellr_biogeo(Yc, δt, narea, ntip, m, nedge)
  Xupd_llr,  Rupd_llr  = makellr_XRupd(δt, narea, wcol)
  ωxupd_llr, σ²upd_llr = makellr_ωxσupd(δt, Yc, ntip)
  bgiid, bgiid_br      = makellf_bgiid(bridx_a, δt, narea, nedge, m)

  # number of xnodes + ωx + ω1 + ω0 + λ1 + λ0 + σ² 
  np = length(wXp) + 6

  # parameter update vector
  parvec = collect(1:np)

  # add to parameter update vector according to weights
  append!(parvec, fill(1,ceil(Int64,np*weight[1])))
  append!(parvec, fill(2,ceil(Int64,np*weight[2])))
  append!(parvec, repeat(3:4,  inner = ceil(Int64,np*weight[3])))
  append!(parvec, repeat(5:6,  inner = ceil(Int64,np*weight[4])))
  append!(parvec, repeat(7:np, inner = ceil(Int64,np*weight[5])))

  # if fixed ωx, ω1 or ω0, then remove
  fix_ωx && filter!(x -> x ≠ 2, parvec)
  fix_ω1 && filter!(x -> x ≠ 3, parvec)
  fix_ω0 && filter!(x -> x ≠ 4, parvec)

  # create update functions for Xbr, Y, Ybr and XYbr
  mhr_upd_Xbr     = make_mhr_upd_Xbr(wcol, m, narea, ntip, nedge, 
                                     bridx, brct, total_llr)
  mhr_upd_Xtrio   = make_mhr_upd_Xtrio(wcol, m, narea, ntip, nedge,
                                       brl, bridx, brct, total_llr)
  mhr_upd_Ybr     = make_mhr_upd_Ybr(narea, nedge, m, ntip, bridx_a, 
                                     brct, brl, wcol,total_llr, bgiid_br)
  mhr_upd_Ytrio   = make_mhr_upd_Ytrio(narea, nedge, m, ntip, bridx_a,  brct, brl, 
                                        wcol, total_llr, bgiid)
  mhr_upd_Ystem   = make_mhr_upd_Ystem(stbrl, narea, nedge)
  mhr_upd_XYbr    = make_mhr_upd_XYbr(narea, nedge, m, ntip, 
                                      bridx, bridx_a, brct, brl, wcol, 
                                      total_llr, bgiid_br)
  mhr_upd_XYtrio  = make_mhr_upd_XYtrio(narea, nedge, m, ntip, 
                                        bridx, bridx_a, brct, brl, 
                                        wcol, total_llr, bgiid)

  # burning phase
  llc, prc, Xc, Yc, LApc, LAnc, LDc, δXc, δYc, 
  stemevc, brs, σ²c, ωxc, ω1c, ω0c, λ1c, λ0c, ptn = 
    burn_tribe(total_llf, ω1upd_llr, ω0upd_llr, λ1upd_llr, λ0upd_llr, 
      Xupd_llr, Rupd_llr, ωxupd_llr, σ²upd_llr, 
      bgiid, bgiid_br, mhr_upd_Xbr, mhr_upd_Xtrio, 
      mhr_upd_Ybr, mhr_upd_Ytrio, mhr_upd_Ystem, mhr_upd_XYbr, mhr_upd_XYtrio,
      Xc, Yc, δXc, δYc, LApc, LAnc, LDc, σ²i, ωxi, ω1i, ω0i, λ1i, λ0i,
      Xnc1, Xnc2, brl, wcol, bridx_a, brs, stemevc, 
      trios, wXp,
      λprior, ωxprior, ω1prior, ω0prior, σ²prior, np, parvec, nburn)

  # log probability of collision
  max_δt = maximum(δt)::Float64

  # log for nthin
  lthin = 0

  # variables to save X and Y 
  if saveXY[1]
    XYsav  = 0
    XYlit  = 0
    xylogs = fld(niter,saveXY[2])
    Xlog   = zeros(Float64, m, ntip, xylogs)
    Ylog   = zeros(Int64,   m, ntip, narea, xylogs)

    if saveDM[1]
      DMsav  = 0
      DMlit  = 0
      dmlogs = fld(niter,saveDM[2])
      LAlog = zeros(Float64, m, ntip, dmlogs)
      LDlog  = zeros(Float64, m, ntip, narea, dmlogs)
    end
  end

  # progress bar
  p = Progress(niter, dt=5, desc="mcmc...", barlen=20, color=:green)

  # create X parameter update function
  mhr_upd_X = make_mhr_upd_X(Xc, Xnc1, Xnc2, wcol, ptn, wXp, m,
                             narea, ntip, Xupd_llr, Rupd_llr)

  # write to file
  open(out_file*".log", "w")

  open(out_file*".log", "a") do f

    print(f, "iteration\tlikelihood\tprior\tomega_x\tomega_1\tomega_0\tsigma2\tlambda_1\tlambda_0\tcollision_probability\n")

    #=
    start MCMC
    =#

    for it = Base.OneTo(niter)

      # update vector order
      shuffle!(parvec)

      for up = parvec

        # update X[i]
        if up > 6
          llc = mhr_upd_X(up, Xc, Yc, δXc, δYc, 
                          λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, llc, LApc, LAnc, LDc)

        # update λ1 
        elseif up == 5
          llc, prc, λ1c = mhr_upd_λ1(λ1c, Yc, λ0c, llc, prc, ω1c, 
                                     LDc, stemevc, brs, λprior,
                                     ptn[5], λ1upd_llr)

          # update internal node
          if rand() < 0.4
            λϕ1, λϕ0 = λϕprop()
            llc = mhr_upd_Ytrio(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
          end

        # if λ0 is updated
        elseif up == 6
          llc, prc, λ0c = mhr_upd_λ0(λ0c, Yc, λ1c, llc, prc, ω0c, 
                                     LDc, stemevc, brs, λprior,
                                     ptn[6], λ0upd_llr)

          # update internal node
          if rand() < 0.4
            λϕ1, λϕ0 = λϕprop()
            llc = mhr_upd_Ytrio(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
          end

        # if σ² is updated
        elseif up == 1
          if ωxc >= 0.0
            llc, prc, σ²c = mhr_upd_σ²(σ²c, Xc, ωxc, llc, prc, ptn[1], 
                                       LApc, σ²prior, σ²upd_llr)
          else
            llc, prc, σ²c = mhr_upd_σ²(σ²c, Xc, ωxc, llc, prc, ptn[1], 
                                       LAnc, σ²prior, σ²upd_llr)
          end

        # update ωx
        elseif up == 2
          llc, prc, ωxc = mhr_upd_ωx(ωxc, Xc, σ²c, llc, prc, ptn[2], 
                                     LApc, LAnc, ωxprior, ωxupd_llr)

        #update ω1
        elseif up == 3
          llc, prc, ω1c = mhr_upd_ω1(ω1c, λ1c, Yc, llc, prc, ptn[3], 
                                     LDc, ω1prior, ω1upd_llr)

          # update internal node
          if rand() < 0.4
            λϕ1, λϕ0 = λϕprop()
            llc = mhr_upd_Ytrio(rand(trios), Xc, Yc, 
                    λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, λϕ1, λϕ0, llc, prc,
                    LApc, LAnc, LDc, δXc, δYc, brs, stemevc)
          end

        # update ω0
        else
          llc, prc, ω0c = mhr_upd_ω0(ω0c, λ0c, Yc, llc, prc, ptn[4],
                                     LDc, ω0prior, ω0upd_llr)

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

      end # end parameter loop

      # log parameters
      lthin += 1
      if lthin == nthin
        pci = Pc(f_λ(λ1c,ω1c,1.0), f_λ(λ0c,ω0c,1.0), max_δt)
        print(f, it,"\t", llc, "\t", prc,"\t",ωxc,"\t",ω1c,"\t",ω0c,"\t",
             σ²c,"\t",λ1c,"\t",λ0c,"\t",pci,"\n")
        lthin = 0
      end

      # log X & Y
      if saveXY[1]
        XYsav += 1
        if XYsav == saveXY[2]
          @inbounds begin
            XYlit += 1
            Xlog[:,:,  XYlit] = Xc
            Ylog[:,:,:,XYlit] = Yc
          end
          XYsav = 0
        end
      end

      # log deterministic matrices
      if saveDM[1]
        DMsav += 1
        if DMsav == saveDM[2]
          @inbounds begin
            DMlit += 1
            LAlog[:,:,  DMlit] = ωxc > 0 ? LApc * ωxc : LAnc * ωxc
            LDlog[:,:,:,DMlit] = LDc 
          end
          DMsav = 0
        end
      end

      next!(p)
    end # end MCMC
  end # end print loop

  if saveXY[1]
    # save X and Y as R objects
    @rput Xlog
    @rput Ylog
    @rput δt
    @rput B

    if saveDM[1]
      @rput LAlog
      @rput LDlog
      reval("""
        delta.t <- δt
        save(Xlog, Ylog, B, delta.t, LAlog, LDlog, 
            file = '$out_file.rda')
      """)
    else
      reval("""
        delta.t <- δt
        save(Xlog, Ylog, B, delta.t, file = '$out_file.rda')
      """)
    end

  end

  return llc, prc, ωxc, ω1c, ω0c, σ²c, λ1c, λ0c, Xc, Yc

end



