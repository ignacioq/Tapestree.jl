#=

MCMC for biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#




"""
    compete_mcmc(...)

Run MCMC for join inference of trait
and biogeographic evolution and competition.
"""
function compete_mcmc(Xc      ::Array{Float64,2},
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
                      ωxprior ::NTuple{2,Float64} = (0.,10.),
                      ω1prior ::NTuple{2,Float64} = (0.,10.),
                      ω0prior ::NTuple{2,Float64} = (0.,10.),
                      σ²prior ::Float64           = 1e-1,
                      λprior  ::Float64           = 1e-1,
                      out_file::String            = "compete_results",
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

  print_with_color(:green, "Data successfully processed", bold = true)

  # dims
  const m, ntip, narea  = size(Yc)

  # coupled nodes for Y
  const Ync1 = ncoup[:,1]
  const Ync2 = ncoup[:,2]

  # expand Y coupled nodes to multidimensional matrix 
  for i = 2:narea
    append!(Ync1, ncoup[:,1] + (i-1)*(m*ntip))
    append!(Ync2, ncoup[:,2] + (i-1)*(m*ntip))
  end

  # tie biogeographic coupled nodes
  Yc[Ync2] = Yc[Ync1]

  # coupled nodes for X
  const Xnc1 = ncoup[:,1]
  const Xnc2 = ncoup[:,2]

  # which nodes are not NaN in Xc
  const wXp = setdiff(find(map(x -> !isnan(x), Xc)), m:m:length(Xc))

  # tie trait coupled nodes
  Xc[Xnc2] = Xc[Xnc1]

  #create object with column indices
  const wcol = create_wcol(Xc)

  # add long stem branch
  edges = cat(1, edges, [2*ntip ntip + 1])
  push!(brl, stbrl)

  # number of edges
  const nedge = size(edges,1) 

  # make edge triads
  const trios = maketriads(edges)

  # make ragged array with index for each edge in Yc
  const bridx = make_edgeind(edges[:,2], B, ntip)

  # expand bridx for each area
  const bridx_a = Array{UnitRange{Int64},1}[]
  push!(bridx_a, bridx)

  for j = 2:narea
    bridxj = deepcopy(bridx)
    for i in Base.OneTo(nedge)
      bridxj[i] += (j-1)*(m*ntip)
    end
    push!(bridx_a, bridxj)    
  end

  # make ragged array of cumulative delta times for each branch
  const brδt = make_edgeδt(bridx, δt, m)

  # array for states at start and end of branches
  brs = ones(Int64, nedge, 2, narea)

  # filter tips in edges
  wtips = map(x -> x <= ntip, edges[:,2])

  ## assign tip area states to brs according Yc
  for j in Base.OneTo(narea)
    brs[wtips,2,j] = Yc[end,:,j]
  end

  # Sample all internal node values according to Pr transitions
  for triad in trios
    upnode!(λ1i, λ0i, triad, Yc, bridx_a, brδt, brl, brs, narea, nedge)
  end

  # tie biogeographic nodes equal
  Yc[Ync2] = Yc[Ync1]

  # assign same value as mrca
  brs[nedge,1,:] = brs[nedge,2,:]

  # create stem events
  stemevc = br_samp(brs[nedge,1,:], brs[nedge,2,:], λ1i, λ0i, brl[nedge], narea)

  # estimate current area & lineage means and area occupancy
  areavg = zeros(m, narea)
  areaoc = zeros(Int64, m, narea)
  linavg = fill(NaN, m, ntip)

  area_lineage_means!(areavg, linavg, areaoc, Xc, Yc, wcol, m, narea)

  # estimate current lineage specific means
  lindiff = fill(NaN, m, ntip, narea)
  linarea_diff!(lindiff, Xc, areavg, areaoc, narea, ntip, m)

  # estimate average branch lineage specific means
  avg_Δx = zeros((nedge-1), narea)
  linarea_branch_avg! = make_la_branch_avg(bridx_a, length(Yc), m, narea, nedge)
  linarea_branch_avg!(avg_Δx, lindiff)


  # make likelihood and prior functions
  total_llf   = makellf(δt, Yc, ntip, narea, m)
  λupd_llr    = makellr_λ_upd(Yc, δt, narea, ntip, m)
  ω10upd_llr  = makellr_ω10_upd(Yc, δt, narea, ntip, m)
  Xupd_llr    = makellr_Xupd(δt, narea)
  Rupd_llr    = makellr_Rupd(δt[1], narea)
  σ²ωxupd_llr = makellr_σ²ωxupd(δt, Yc, ntip)  
  bgiid       = makellf_bgiid(bridx_a, δt, narea, nedge, m)
  bgiid_br    = makellf_bgiid_br(bridx_a, δt, narea, nedge, m)

  # number of xnodes + ωx + ω1 + ω0 + λ1 + λ0 + σ² 
  np = length(wXp) + 6

  # parameter update vector
  const parvec = collect(1:np)

  # add to parameter update vector according to weights
  append!(parvec, fill(1,ceil(Int64,np*weight[1])))
  append!(parvec, fill(2,ceil(Int64,np*weight[2])))
  append!(parvec, repeat(3:4,  inner = ceil(Int64,np*weight[3])))
  append!(parvec, repeat(5:6,  inner = ceil(Int64,np*weight[4])))
  append!(parvec, repeat(7:np, inner = ceil(Int64,np*weight[5])))

  # if fixed ωx, ω1 or ω0 remove
  fix_ωx && filter!(x -> x ≠ 2, parvec)
  fix_ω1 && filter!(x -> x ≠ 3, parvec)
  fix_ω0 && filter!(x -> x ≠ 4, parvec)

  # create update functions for Xbr, Y, Ybr and XYbr
  mhr_upd_Xbr  = make_mhr_upd_Xbr(wcol, m, narea, ntip, nedge, total_llf)
  mhr_upd_Y    = make_mhr_upd_Y(narea, nedge, m, ntip, bridx_a, 
                                brδt, brl, wcol, Ync1, Ync2, 
                                total_llf, bgiid, linarea_branch_avg!)
  mhr_upd_Ybr  = make_mhr_upd_Ybr(narea, nedge, m, ntip, bridx_a, 
                                  brδt, brl, wcol, Ync1, Ync2, 
                                  total_llf, bgiid_br, linarea_branch_avg!)
  mhr_upd_XYbr = make_mhr_upd_XYbr(narea, nedge, m, ntip, 
                                   bridx_a, brδt, brl, wcol, 
                                   total_llf, bgiid_br, linarea_branch_avg!)

  # burning phase
  llc, prc, Xc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx,
  stemevc, brs, σ²c, ωxc, ω1c, ω0c, λ1c, λ0c, ptn = burn_compete(
    total_llf, λupd_llr, ω10upd_llr, Xupd_llr, Rupd_llr, σ²ωxupd_llr, 
    bgiid, bgiid_br, linarea_branch_avg!, 
    mhr_upd_Xbr, mhr_upd_Y, mhr_upd_Ybr, mhr_upd_XYbr,
    Xc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx,
    σ²i, ωxi, ω1i, ω0i, λ1i, λ0i,
    Ync1, Ync2, Xnc1, Xnc2, brl, wcol, bridx_a, brδt, brs, stemevc, 
    trios, wXp,
    λprior, ωxprior, ω1prior, ω0prior, σ²prior, np, parvec, nburn)

  # log likelihood and prior
  h[1] = llc::Float64
  o[1] = prc::Float64

  # log probability of collision
  const max_δt = maximum(δt)::Float64
  pc[1] = Pc(λ1c, λ0c, max_δt)::Float64

  # log for nthin
  lthin = 0

  # variables to save X and Y 
  if saveXY[1]
    const XYsav = 0
    const XYlit = 0
    xylogs = fld(niter,saveXY[2])
    Xlog   = zeros(Float64, m, ntip, xylogs)
    Ylog   = zeros(Int64,   m, ntip, narea, xylogs)
  end

  # number of internal nodes (to perform updates on)
  const nin = length(trios) + 1

  # progess bar
  p = Progress(niter, dt=5, desc="running MCMC...", barlen=20, color=:green)

  # create parameter update functions
  mhr_upd_X    = make_mhr_upd_X(Xnc1, Xnc2, wcol, m, ptn, wXp, 
                                narea, ntip, Xupd_llr, Rupd_llr)

  #=
  start MCMC
  =#


  open(out_file*".log", "w")

  open(out_file*".log", "a") do f

    print(f, "iteration\tlikelihood\tprior\tomega_x\tomega_1\tomega_0\tsigma2\tlambda_1\tlambda_0\tcollision_probability\n")

    for it = Base.OneTo(niter)

      # update vector order
      shuffle!(parvec)

      for up = parvec

        # update X[i]
        if up > 6

          llc = mhr_upd_X(up, Xc, Yc, λ1c, λ0c, 
                          ωxc, ω1c, ω0c, σ²c, llc, 
                          areavg, linavg, lindiff, areaoc)

        # update λ1 
        elseif up == 5

          llc, prc, λ1c = mhr_upd_λ1(λ1c, Yc, λ0c, llc, prc, ω1c, ω0c, 
                                     lindiff, stemevc, brs[nedge,1,:], λprior,
                                     ptn[5], λupd_llr)

          # which internal node to update
          if rand() < 0.4
            bup = rand(Base.OneTo(nin))
            # update a random internal node, including the mrca
            if bup < nin
              llc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx = 
                mhr_upd_Y(trios[bup], 
                          Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, llc, prc, 
                          areavg, areaoc, linavg, lindiff, avg_Δx, brs, stemevc)
            else
              # update stem
              llr = 0.0
              for j=Base.OneTo(narea)
                @inbounds llr -= brll(stemevc[j], λ1c, λ0c, brs[nedge,1,j])
              end

              stemevc = upstem(λ1c, λ0c, nedge, brs, brl, narea)

              for j=Base.OneTo(narea)
                @inbounds llr += brll(stemevc[j], λ1c, λ0c, brs[nedge,1,j])
              end

              llc += llr
            end
          end

        # if λ0 is updated
        elseif up == 6

          llc, prc, λ0c = mhr_upd_λ0(λ0c, Yc, λ1c, llc, prc, ω1c, ω0c, 
                                     lindiff, stemevc, brs[nedge,1,:], λprior,
                                     ptn[6], λupd_llr)

          # which internal node to update
          if rand() < 0.4
            bup = rand(Base.OneTo(nin))
            # update a random internal node, including the mrca
            if bup < nin
              llc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx = 
                mhr_upd_Y(trios[bup], 
                          Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, llc, prc, 
                          areavg, areaoc, linavg, lindiff, avg_Δx, brs, stemevc)
            else
              # update stem
              llr = 0.0
              for j=Base.OneTo(narea)
                @inbounds llr -= brll(stemevc[j], λ1c, λ0c, brs[nedge,1,j])
              end

              stemevc = upstem(λ1c, λ0c, nedge, brs, brl, narea)

              for j=Base.OneTo(narea)
                @inbounds llr += brll(stemevc[j], λ1c, λ0c, brs[nedge,1,j])
              end

              llc += llr
            end
          end

        # if σ² is updated
        elseif up == 1
          llc, prc, σ²c = mhr_upd_σ²(σ²c, Xc, ωxc, llc, prc, ptn[1], 
                                     linavg, σ²prior, σ²ωxupd_llr)

        # update ωx
        elseif up == 2
          llc, prc, ωxc = mhr_upd_ωx(ωxc, Xc, σ²c, llc, prc, ptn[2], 
                                     linavg, ωxprior, σ²ωxupd_llr)

        #update ω1
        elseif up == 3
          llc, prc, ω1c = mhr_upd_ω1(ω1c, λ1c, λ0c, ω0c, Yc, llc, prc, ptn[3], 
                                     linavg, lindiff, ω1prior, ω10upd_llr)

          # which internal node to update
          if rand() < 0.4
            bup = rand(Base.OneTo(nin))
            # update a random internal node, including the mrca
            if bup < nin
              llc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx = 
                mhr_upd_Y(trios[bup], 
                          Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, llc, prc, 
                          areavg, areaoc, linavg, lindiff, avg_Δx, brs, stemevc)
            else
              # update stem
              llr = 0.0
              for j=Base.OneTo(narea)
                @inbounds llr -= brll(stemevc[j], λ1c, λ0c, brs[nedge,1,j])
              end

              stemevc = upstem(λ1c, λ0c, nedge, brs, brl, narea)

              for j=Base.OneTo(narea)
                @inbounds llr += brll(stemevc[j], λ1c, λ0c, brs[nedge,1,j])
              end

              llc += llr
            end
          end

        # update ω0      
        else
          llc, prc, ω0c = mhr_upd_ω0(ω0c, λ1c, λ0c, ω1c, Yc, llc, prc, ptn[4],
                                     linavg, lindiff, ω0prior, ω10upd_llr)

          # which internal node to update
          if rand() < 0.4
            bup = rand(Base.OneTo(nin))
            # update a random internal node, including the mrca
            if bup < nin
              llc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx = 
                mhr_upd_Y(trios[bup], 
                          Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, llc, prc, 
                          areavg, areaoc, linavg, lindiff, avg_Δx, brs, stemevc)
            else
              # update stem
              llr = 0.0
              for j=Base.OneTo(narea)
                @inbounds llr -= brll(stemevc[j], λ1c, λ0c, brs[nedge,1,j])
              end

              stemevc = upstem(λ1c, λ0c, nedge, brs, brl, narea)

              for j=Base.OneTo(narea)
                @inbounds llr += brll(stemevc[j], λ1c, λ0c, brs[nedge,1,j])
              end

              llc += llr
            end
          end
        end

        ## make a branch updates with Pr = 0.01
        # make X branch update
        if 0.01 < rand()
          llc, Xc, areavg, areaoc, linavg, lindiff = 
            mhr_upd_Xbr(Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, llc, 
                        areavg, linavg, lindiff, areaoc)
        end

        # make Y branch update
        if 0.01 < rand()
          llc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx = 
            mhr_upd_Ybr(rand(Base.OneTo(nedge-1)), Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, 
                        llc, prc, areavg, areaoc, linavg, lindiff, avg_Δx, 
                        brs, stemevc)
        end

        # make joint X & Y branch update
        if 0.01 < rand()
          llc, Xc, Yc, areavg, areaoc, linavg, lindiff, avg_Δx =
            mhr_upd_XYbr(rand(Base.OneTo(nedge-1)), Xc, Yc, λ1c, λ0c, ωxc, ω1c, ω0c, σ²c, 
                         llc, prc, areavg, areaoc, linavg, lindiff, avg_Δx, 
                         brs, stemevc)
        end
      end

      # log parameters
      lthin += 1
      if lthin == nthin
        pci = Pc(f_λ(λ1c,ω1c,1.0), f_λ(λ0c,ω0c,1.0), max_δt)
        print(f, "$it\t$llc\t$prc\t$ωxc\t$ω1c\t$ω0c\t$σ²c\t$λ1c\t$λ0c\t$pci\n")
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

      next!(p)
    end #end MCMC
  end # end print loop

  if saveXY[1]
    # save X and Y as R objects
    @rput Xlog
    @rput Ylog
    @rput δt
    @rput B
    reval("""
      delta.t <- δt
      save(Xlog, Ylog, B, delta.t, file = '$out_file.rda')
    """)
  end

  return llc, prc, ωxc, ω1c, ω0c, σ²c, λ1c, λ0c, Xc, Yc

end



