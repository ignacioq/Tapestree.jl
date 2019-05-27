#=

Utility functions for simulations for tribe

Ignacio Quintero Mächler

August 29 2017

t(-_-t)

=#





"""
    simulate_tribe(X_initial::Float64,
                   nareas   ::Int64,
                   tree_file::String;
                   ωx       = 0.0,
                   σ²       = 0.5,
                   λ1       = 0.5,
                   λ0       = 0.2,
                   ω1       = 0.0,
                   ω0       = 0.0,
                   const_δt = 1e-4)

Simulate tribe model.
"""
function simulate_tribe(X_initial::Float64,
                        nareas   ::Int64,
                        tree_file::String;
                        ωx       = 0.0,
                        σ²       = 0.5,
                        λ1       = 0.5,
                        λ0       = 0.2,
                        ω1       = 0.0,
                        ω0       = 0.0,
                        const_δt = 1e-4)

  # bounds checks for parameters
  0.0 >= σ² && error("σ² has to be > 0.0")
  0.0 >= λ1 && error("λ1 has to be > 0.0")
  0.0 >= λ0 && error("λ0 has to be > 0.0")

  Y_initial = 
    [rand() < λ1/(λ1 + λ0) ? 1 : 0 for i in Base.OneTo(nareas)]::Array{Int64,1}

  while iszero(sum(Y_initial))
    Y_initial = 
      [rand() < λ1/(λ1 + λ0) ? 1 : 0 for i in Base.OneTo(nareas)]::Array{Int64,1}
  end

  tree, bts = read_tree(tree_file)

  br = branching_times(tree)

  # sort according to branching times
  brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)

  # sort branching times
  sort!(bts, rev = true)

  # add present
  push!(bts, 0.0)

  # number of speciation events
  nbt = lastindex(bts) - 1

  # calculate speciation waiting times
  swt = Array{Float64,1}(undef,nbt)
  for i in Base.OneTo(nbt)
    swt[i] = bts[i] - bts[i+1]
  end

  # initial values
  Xt = fill(X_initial, 2)

  Y_initial = reshape(vec(Y_initial), 1, :)
  Yt = vcat(Y_initial, Y_initial)

  # start of alive
  alive = sortslices(br, dims = 1, by = x -> x[1])[1:2,2]

  nalive = lastindex(alive)

  # is ωx positive? (for lineage averages)
  isωxP = ωx >= 0.0

  rate = sqrt(const_δt*σ²)

  # loop through waiting times
  for j in Base.OneTo(nbt)

    nreps = reps_per_period(swt[j], const_δt)

    # simulate during the speciation waiting time
    nconst_sim!(Xt, Yt, nreps, const_δt, ωx, ω1, ω0, rate, λ1, λ0, isωxP)

    if j == nbt
      break
    end

    # which lineage speciates
    wsp = brs[j,2]

    # where to insert
    wti = findfirst(isequal(wsp), alive)

    # index to insert
    idx = sort(push!(collect(1:nalive), wti))

    # insert in Xt & Yt
    Xt = Xt[idx]
    Yt = Yt[idx,:]

    # update alive
    chs = brs[wsp .== brs[:,1],2]

    alive[wti] = chs[1]
    insert!(alive, wti+1, chs[2])

    nalive = length(alive)

  end

  tip_traits = 
    Dict(convert(Int64, alive[i]) => Xt[i]   for i = Base.OneTo(nalive))
  tip_areas  = 
    Dict(convert(Int64, alive[i]) => Yt[i,:] for i = Base.OneTo(nalive))

  pop!(bts)

  return tip_traits, tip_areas, tree, bts

end





"""
    reps_per_period(br_length::Float64, const_δt::Float64)

Estimate number of repetitions for each speciation waiting time.
"""
reps_per_period(br_length::Float64, const_δt::Float64) = 
  round(Int64,cld(br_length,const_δt))





"""
    nconst_sim!(Xt   ::Array{Float64,1}, 
                Yt   ::Array{Int64,2},
                nreps::Int64,
                δt   ::Float64,
                rate ::Float64,
                ωx   ::Float64, 
                λ1   ::Float64, 
                λ0   ::Float64, 
                ω1   ::Float64, 
                ω0   ::Float64,
                isωxP::Bool)

Simulate biogeographic and trait evolution along a 
speciation waiting time.
"""
function nconst_sim!(Xt   ::Array{Float64,1}, 
                     Yt   ::Array{Int64,2},
                     nreps::Int64,
                     δt   ::Float64,
                     ωx   ::Float64, 
                     ω1   ::Float64, 
                     ω0   ::Float64,
                     rate ::Float64,
                     λ1   ::Float64, 
                     λ0   ::Float64, 
                     isωxP::Bool)

  # n species and narea areas
  n, narea = size(Yt)
  Ytp      = similar(Yt)
  nch      = zeros(Int64, n)

  # allocate memory for lineage averages and differences
  δX = zeros(n, n)      # X pairwise differences
  δY = zeros(n, n)      # Y pairwise differences
  la = zeros(n)         # lineage averages
  ld = zeros(n, narea)  # area specific lineage rates

  for i in Base.OneTo(nreps)

    # estimate area and lineage averages and area occupancy
    δXY_la_ld!(δX, δY, la, ld, Xt, Yt, n, narea, isωxP)

    # trait step
    traitsam_1step!(Xt, la, δt, rate, ωx, n)

    # biogeographic step
    copyto!(Ytp, Yt)
    biogeosam_1step!(ω1, ω0, λ1, λ0, ld, Ytp, nch, δt, n, narea)

    while check_sam(Ytp, nch, n, narea)
      copyto!(Ytp, Yt)
      biogeosam_1step!(ω1, ω0, λ1, λ0, ld, Ytp, nch, δt, n, narea)
    end

    copyto!(Yt,Ytp)
  end

  return nothing
end





"""
    δXY_la_ld!(δX   ::Array{Float64,2},
               δY   ::Array{Float64,2},
               la   ::Array{Float64,1},
               ld   ::Array{Float64,2},
               Xt   ::Array{Float64,1}, 
               Yt   ::Array{Int64,2}, 
               n    ::Int64,
               narea::Int64)

Estimate area and lineage specific averages given sympatry configuration.
"""
function δXY_la_ld!(δX   ::Array{Float64,2},
                    δY   ::Array{Float64,2},
                    la   ::Array{Float64,1},
                    ld   ::Array{Float64,2},
                    Xt   ::Array{Float64,1}, 
                    Yt   ::Array{Int64,2}, 
                    n    ::Int64,
                    narea::Int64,
                    isωxP::Bool)
  @inbounds begin

    # estimate pairwise distances
    for j = Base.OneTo(n), i = Base.OneTo(n)
      i == j && continue
      # δX differences
      δX[i,j] = Xt[i] - Xt[j]
      # δY differences
      sk      = 0.0
      δY[i,j] = 0.0
      @simd for k = Base.OneTo(narea)
        if Yt[j,k] == 1
            sk += 1.0
            δY[i,j] += Float64(Yt[i,k])
        end
      end
      δY[i,j] /= sk
    end

    # estimate lineage averages
    la[:] .= 0.0
    if isωxP
      for j = Base.OneTo(n), i = Base.OneTo(n)
          i == j && continue
          la[i] += δX[j,i] * δY[j,i]
      end
    else
      for j = Base.OneTo(n), i = Base.OneTo(n)
          i == j && continue
          la[i] += sign(δX[j,i]) * δY[j,i] * exp(-abs(δX[j,i]))
      end
    end

    # estimate lineage sum of distances
    ld[:] .= NaN
    for k = Base.OneTo(narea), i = Base.OneTo(n)
      xmin = 1.0e20
      for j = Base.OneTo(n)
        j == i && continue
        if Yt[j,k] == 1
          x = abs(δX[j,i])
          iszero(x) && continue
          xmin = x < xmin ? x : xmin
        end
      end
      ld[i,k] = xmin == 1.0e20 ? 0.0 : xmin
    end
  end

  return nothing
end





"""
    traitsam_1step!(Xt  ::Array{Float64,1}, 
                    la  ::Array{Float64,1}, 
                    δt  ::Float64, 
                    rate::Float64, 
                    ωx  ::Float64, 
                    n   ::Int64)

Sample one step for trait evolution history: `X(t + δt)`.
"""
function traitsam_1step!(Xt  ::Array{Float64,1}, 
                         la  ::Array{Float64,1}, 
                         δt  ::Float64, 
                         rate::Float64, 
                         ωx  ::Float64, 
                         n   ::Int64)

  @inbounds begin

    for i in Base.OneTo(n)
      Xt[i] += Eδx(la[i], ωx, δt) + randn()*rate
    end

  end
end





"""
    biogeosam_1step!(ω1   ::Float64,
                     ω0   ::Float64,
                     λ1   ::Float64,
                     λ0   ::Float64,
                     ld   ::Array{Float64,2},
                     Ytp  ::Array{Int64,2},
                     nch  ::Array{Int64,1},
                     δt   ::Float64,
                     n    ::Int64,
                     narea::Int64)

Sample one step for biogeographic history: `Y(t + δt)`.
"""
function biogeosam_1step!(ω1   ::Float64,
                          ω0   ::Float64,
                          λ1   ::Float64,
                          λ0   ::Float64,
                          ld   ::Array{Float64,2},
                          Ytp  ::Array{Int64,2},
                          nch  ::Array{Int64,1},
                          δt   ::Float64,
                          n    ::Int64,
                          narea::Int64)

  @inbounds begin

    nch[:] .= 0
    for k = Base.OneTo(narea), i = Base.OneTo(n)
      if iszero(Ytp[i,k])
        if rand() < f_λ(λ1,ω1,ld[i,k])*δt
          setindex!(Ytp,1,i,k)
          nch[i] += 1
        end
      else 
        if rand() < f_λ(λ0,ω0,ld[i,k])*δt
          setindex!(Ytp,0,i,k)
          nch[i] += 1 
        end
      end
    end

  end

  return nothing
end





"""
    check_sam(Ytc::Array{Int64,2}, nch::Array{Int64,1}, n::Int64, k::Int64)

Returns false if biogeographic step consists of *only one* change 
or if the species does *not* go globally extinct. 
"""
function check_sam(Ytp   ::Array{Int64,2}, 
                   nch   ::Array{Int64,1}, 
                   n     ::Int64, 
                   nareas::Int64)

  @inbounds begin

    # check lineage did not go extinct
    for i in Base.OneTo(n)
      s = 0
      @simd for k in Base.OneTo(nareas)
        s += Ytp[i,k]
      end
      if iszero(s)
        return true
      end
    end

    # check that, at most, only one event happened per lineage
    for i in nch
      if i > 1
        return true
      end
    end

    return false
  end
end









