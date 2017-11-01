#=

MCMC for biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#




"""
    compete_mcmc(...)

Run MCMC for join inference of trait
and biogeograhic evolution and competition.
"""
function compete_mcmc(Xc       ::Array{Float64,2},
                      Yc       ::Array{Int64,3},
                      ncoup    ::Array{Int64,2},
                      δt       ::Array{Float64,1},
                      edges    ::Array{Int64,2},
                      brl      ::Array{Float64,1},
                      B        ::Array{Float64,2};
                      niter    ::Int64             = 500_000,
                      nthin    ::Int64             = 1_000,
                      nburn    ::Int64             = 500_000,
                      ωxprior  ::NTuple{2,Float64} = (0.,10.),
                      ωλprior  ::NTuple{2,Float64} = (0.,10.),
                      ωμprior  ::NTuple{2,Float64} = (0.,10.),
                      σ²prior  ::Float64           = 1e-1,
                      λprior   ::Float64           = 1e-1,
                      out_file ::String            = "compete_results",
                      weight   ::NTuple{3,Float64} = (0.1,0.02,0.02),
                      λi       ::Array{Float64,1}  = [1.0, 0.01],
                      ωxi      ::Float64           = 0.,
                      ω1i      ::Float64           = 0.,
                      ω0i      ::Float64           = 0.,
                      σ²i      ::Float64           = 1.,
                      stbrl    ::Float64           = 1.,
                      fix_ωλ_ωμ::Bool              = true)

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
  const wNaN   = find(.!isnan.(Xc))
  const wspp   = m:m:length(Xc)
  const wXp    = setdiff(wNaN,wspp)

  # tie trait coupled nodes
  Xc[Xnc2] = Xc[Xnc1]

  #create object with column indices
  const wcol = create_wcol(Xc)

  # initialize result arrays
  const nlogs = fld(niter,nthin)         # number of logged iterations
  const iter  = zeros(Float64, nlogs)             # iterations
  const ωx    = zeros(Float64, nlogs)             # trait competition parameter
  const ωλ    = zeros(Float64, nlogs)             # colonization competition parameter
  const ωμ    = zeros(Float64, nlogs)             # extinction competition parameter
  const σ²    = zeros(Float64, nlogs)             # drift parameter
  const λs    = zeros(Float64, nlogs, 2)          # rate parameters
  const h     = zeros(Float64, nlogs)             # likelihood
  const o     = zeros(Float64, nlogs)             # prior
  const pc    = zeros(Float64, nlogs)             # collision probability

  # add long stem branch
  edges = cat(1, edges, [2*ntip ntip + 1])
  push!(brl, stbrl)

  const nedge = size(edges,1) 

  # make edge triads
  const trios = maketriads(edges)

  # make ragged array with index for each edge in Yc
  const bridx = make_edgeind(edges[:,2], B, ntip)

  # expand bridx for each area
  const bridx_a = Vector{Vector{Int64}}[]
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
  brs = zeros(Int64, nedge, 2, narea)

  # assign states to brs according Yc
  for i = Base.OneTo(nedge), j in Base.OneTo(narea)
    wh1 = find(edges[:,1] .== edges[i,2])
    brs[i,2,j] = brs[wh1,1,j] = Yc[bridx_a[j][i]][end]
  end

  # Sample all internal node values according to Pr transitions
  for triad in trios
    upnode!(λi, triad, Yc, bridx_a, brδt, brl, brs, narea, nedge)
  end

  # assign same value as mrca
  brs[nedge,1,:] = brs[nedge,2,:]

  # create stem events
  stemevc = br_samp(brs[nedge,1, :], brs[nedge,2,:], λi, brl[nedge], narea)

  # make values at nodes equal
  Yc[Ync2] = Yc[Ync1]

  # estimate current area & lineage means
  const areavg = zeros(m,narea)
  const linavg = zeros(m, ntip)

  area_lineage_means!(areavg, linavg, Xc, Yc, wcol, m)

  # estimate current lineage specific means
  const lindiff = zeros(m, ntip, narea)
  linarea_diff!(lindiff, Xc, areavg, narea, ntip, m)

  # estimate average branch lineage specific means
  const avg_Δx = zeros(nedge, narea)
  linarea_branch_avg!(avg_Δx, lindiff, bridx_a, narea, nedge)

  # make likelihood and prior functions
  total_llf      = makellf(δt, Yc, ntip, wcol, narea)
  λupd_llf       = makellf_λ_upd(Yc, δt, narea)
  ωλμupd_llf     = makellf_ωλμ_upd(Yc, δt, narea)
  Xupd_llf       = makellf_Xupd(δt, narea)
  Rupd_llf       = makellf_Rupd(δt, narea)
  σ²ωxupd_llf    = makellf_σ²ωxupd(δt, Yc, ntip)  
  biogeo_upd_iid = makellf_biogeo_upd_iid(bridx_a, δt, narea, nedge, m)

  # number of free parameters
  # number of xnodes + λ1 + λ0 + σ² + ωx + ωλ + ωμ
  const np = length(wXp) + 6

  # burning phase
  llc, prc, Xc, Yc, areavg, linavg, lindiff, avg_Δx,
  stemevc, brs, λc, ωxc, ω1c, ω0c, σ²c, ptn = burn_compete(total_llf, 
      λupd_llf, ωλμupd_llf, Xupd_llf, Rupd_llf, σ²ωxupd_llf, biogeo_upd_iid, 
      Xc, Yc, areavg, linavg, lindiff, avg_Δx,
      λi, ωxi, ω1i, ω0i, σ²i, 
      Ync1, Ync2, Xnc1, Xnc2, brl, wcol, bridx_a, brδt, brs, stemevc, trios, wXp, 
      weight, λprior, ωxprior, ωλprior, ωμprior, σ²prior, fix_ωλ_ωμ, np, nburn)

  # log likelihood and prior
  h[1]  = llc
  o[1]  = prc

  # log probability of collision
  const max_δt = maximum(δt)
  pc[1] = Pc(λc[1], λc[2], max_δt)

  # log for nthin
  const lit   = 0
  const lthin = 0

  # number of internal nodes (to perform updates on)
  const nin = length(trios) + 1

  # parameter location for λ
  const λlessthan = 6

  # progess bar
  p  = Progress(niter, 5, "running MCMC...", 20)

  if fix_ωλ_ωμ
    const pv      = append!(collect(1:np),
                            repeat(1:2, inner = ceil(Int64,np*weight[1])))
    append!(pv, repeat(5:6, inner = ceil(Int64,np*weight[3])))
    const parvec  = setdiff(pv,(3:4))
    const lparvec = length(parvec)
  else
    const parvec  = append!(collect(1:np),
                            repeat(1:2, inner = ceil(Int64,np*weight[1])))
    append!(parvec, repeat(3:4, inner = ceil(Int64,np*weight[2])))
    append!(parvec, repeat(5:6, inner = ceil(Int64,np*weight[3])))
    const lparvec = length(parvec)
  end

  # create parameter update functions

  mhr_upd_λ = make_mhr_upd_λ(nedge, λprior, ptn, λupd_llf)
  mhr_upd_Y = make_mhr_upd_Y(narea, nedge, m, ntip, bridx_a, 
                             brδt, brl, wcol, Ync1, Ync2, 
                             total_llf, biogeo_upd_iid)
  mhr_upd_X = make_mhr_upd_X(Xnc1, Xnc2, wcol, m, ptn, wXp, 
                             λlessthan, narea, Xupd_llf, Rupd_llf)

  #start MCMC
  for it = Base.OneTo(niter)

    # Update vector
    upvector = rand(parvec,lparvec)

    for up = upvector

      # update X[i]
      if up > λlessthan
        Xc, llc, areavg, linavg, lindiff = mhr_upd_X(up, Xc, Yc, λc, 
                                            ωxc, ω1c, ω0c, σ²c, llc, 
                                            areavg, linavg, lindiff)

      #randomly select λ to update and branch histories
      elseif up > 4 && up <= λlessthan

        llc, prc, λc = mhr_upd_λ(up, Yc, λc, llc, prc, ω1c, ω0c, 
                                 lindiff, stemevc, brs[nedge,1,:])

        # which internal node to update
        if rand() < 0.2
          bup = rand(Base.OneTo(nin))
          # update a random internal node, including the mrca
          if bup < nin
            llc, Yc, areavg, linavg, lindiff, avg_Δx = mhr_upd_Y(trios[bup], 
                       Xc, Yc, λc, ωxc, ω1c, ω0c, σ²c, llc, prc, 
                       areavg, linavg, lindiff, avg_Δx, brs, stemevc)
          else
            # update stem
            llr = 0.0
            for j=Base.OneTo(narea)
              @inbounds llr -= brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])
            end

            stemevc = upstem(λc, nedge, brs, brl, narea)

            for j=Base.OneTo(narea)
              @inbounds llr += brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])
            end

            llc += llr
          end
        end

      # if σ² is updated
      elseif up == 1
        llc, prc, σ²c = mhr_upd_σ²(σ²c, Xc, ωxc, llc, prc, ptn[1], 
                                   linavg, σ²prior, σ²ωxupd_llf)

      # update ωx
      elseif up == 2
        llc, prc, ωxc = mhr_upd_ωx(ωxc, Xc, σ²c, llc, prc, ptn[2], 
                                   linavg, ωxprior, σ²ωxupd_llf)

      #update ωλ
      elseif up == 3
        llc, prc, ω1c = mhr_upd_ωλ(ω1c, λc, ω0c, Yc, llc, prc, ptn[3], 
                                   linavg, lindiff, ωλprior, ωλμupd_llf)

        # which internal node to update
        if rand() < 0.2
          bup = rand(Base.OneTo(nin))
          # update a random internal node, including the mrca
          if bup < nin
            llc, Yc, areavg, linavg, lindiff, avg_Δx = mhr_upd_Y(trios[bup], 
                       Xc, Yc, λc, ωxc, ω1c, ω0c, σ²c, llc, prc, 
                       areavg, linavg, lindiff, avg_Δx, brs, stemevc)
          else
            # update stem
            llr = 0.0
            for j=Base.OneTo(narea)
              @inbounds llr -= brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])
            end

            stemevc = upstem(λc, nedge, brs, brl, narea)

            for j=Base.OneTo(narea)
              @inbounds llr += brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])
            end

            llc += llr
          end
        end

      # update ωμ      
      else
        llc, prc, ω0c = mhr_upd_ωμ(ω0c, λc, ω1c, Yc, llc, prc, ptn[4],
                                    linavg, lindiff, ωμprior, ωλμupd_llf)
     
        # which internal node to update
        if rand() < 0.2
          bup = rand(Base.OneTo(nin))
          # update a random internal node, including the mrca
          if bup < nin
            llc, Yc, areavg, linavg, lindiff, avg_Δx = mhr_upd_Y(trios[bup], 
                       Xc, Yc, λc, ωxc, ω1c, ω0c, σ²c, llc, prc, 
                       areavg, linavg, lindiff, avg_Δx, brs, stemevc)
          else
            # update stem
            llr = 0.0
            for j=Base.OneTo(narea)
              @inbounds llr -= brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])
            end

            stemevc = upstem(λc, nedge, brs, brl, narea)

            for j=Base.OneTo(narea)
              @inbounds llr += brll(stemevc[j], λc[1], λc[2], brs[nedge,1,j])
            end

            llc += llr
          end
        end
      end

    end

    lthin += 1
    if lthin == nthin
      @inbounds begin
        lit += 1
        setindex!(iter,  it, lit)
        setindex!(h,    llc, lit)
        setindex!(o,    prc, lit)
        setindex!(ωx,   ωxc, lit)
        setindex!(ωλ,   ω1c, lit)
        setindex!(ωμ,   ω0c, lit)
        setindex!(σ²,   σ²c, lit)
        setindex!(pc, Pc(λc[1], λc[2], max_δt), lit)
        for j = eachindex(λc)
          setindex!(λs, λc[j], lit, j)
        end
      end
      lthin = 0
    end

    next!(p)
  end

  R = hcat(iter, h, o, ωx, ωλ, ωμ, σ², λs, pc)

  # add column names
  col_nam = ["Iteration", "Likelihood", "Prior", "Trait_competition", 
             "Colonization_competition", "Extinction_competition", 
             "Sigma2"]
  for i = 1:-1:0
    xn = *("Lambda_", string(i))
    push!(col_nam, xn)
  end
  
  push!(col_nam, "Collision_probability")

  R = vcat(reshape(col_nam, 1, endof(col_nam)), R)

  writedlm(out_file*".log", R)

  return R
end
