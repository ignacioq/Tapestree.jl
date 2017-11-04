#=

Utility functions for simulations in Compete

Ignacio Quintero Mächler

August 29 2017

t(-_-t)

=#



"""
    simulate_compete(X_initial::Float64, Y_initial::Array{Int64,1}, tree_file::String;  ωx = 0.0, σ = 0.5, λ1 = 0.5, λ0 = 0.2, ω1 = 0.0, ω0 = 0.0, const_δt = 1e-4)

Simulate biogeographic and trait evolution according to Compete model.
"""
function simulate_compete(X_initial::Float64,
                          nareas   ::Int64,
                          tree_file::String;
                          ωx       = 0.0,
                          σ        = 0.5,
                          λ1       = 0.5,
                          λ0       = 0.2,
                          ω1       = 0.0,
                          ω0       = 0.0,
                          const_δt = 1e-4)

  Y_initial::Array{Int64,1} = 
    [rand() < λ1/(λ1 + λ0) ? 1 : 0 for i in Base.OneTo(nareas)]

  while sum(Y_initial) == 0
    Y_initial = 
      [rand() < λ1/(λ1 + λ0) ? 1 : 0 for i in Base.OneTo(nareas)]
  end

  tree, bts = read_tree(tree_file)

  const br = branching_times(tree)

  # sort according to branching times
  const brs = sortrows(br, by = x->(x[5]), rev = true)

  # sort branching times
  sort!(bts, rev = true)

  # add present
  push!(bts, 0.0)

  # number of speciation events
  nbt = endof(bts)-1

  # calculate speciation waiting times
  const swt = Array{Float64,1}(nbt)
  for i in Base.OneTo(nbt)
    swt[i] = bts[i] - bts[i+1]
  end

  # initial values
  Xt = fill(X_initial, 2)

  Y_initial = reshape(Y_initial, 1, :)
  Yt = vcat(Y_initial, Y_initial)

  # start of alive
  alive  = sortrows(br, by = x->(x[1]))[1:2,2]
  nalive = length(alive)

  # loop through waiting times
  for j in Base.OneTo(nbt)

    nreps = reps_per_period(swt[j], const_δt)

    # simulate durin the waiting time
    Yt = branch_sim!(Xt, Yt, nreps, const_δt, ωx, σ, λ1, λ0, ω1, ω0)

    if j == nbt
      break
    end

    # which lineage speciates
    wsp = brs[j,2]

    # where to insert
    wti = find(wsp .== alive)[1]

    # index to insert
    idx = sort(push!(collect(1:nalive), wti))

    # insert in Xt & Yt
    Xt = Xt[idx]
    Yt = Yt[idx,:]

    # update alive
    chs = brs[find(wsp .== brs[:,1]),2]

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
    branch_sim!(Xt::Array{Float64,1}, Yt::Array{Int64,2}, br_δts::Array{Float64,1})

Simulate biogeographic and trait evolution in a branch.
"""
function branch_sim!(Xt    ::Array{Float64,1}, 
                     Yt    ::Array{Int64,2},
                     nreps ::Int64,
                     δt    ::Float64,
                     ωx    ::Float64, 
                     σ     ::Float64, 
                     λ1    ::Float64, 
                     λ0    ::Float64, 
                     ω1    ::Float64, 
                     ω0    ::Float64)

  # n species and k areas
  const n, k = size(Yt)

  # allocate memory for area and lineage averages 
  arav = zeros(k)        # area averages
  liav = zeros(n)        # lineage averages
  aroc = zeros(Int64,k)  # area occupancy
  alλ1 = zeros(n,k)      # area specific lineage rates
  alλ0 = zeros(n,k)      # area specific lineage rates

  for i in Base.OneTo(nreps)

    # estimate area and lineage averages and area occupancy
    arav, liav, aroc = ar_lin_avg(Xt, Yt, arav, liav, aroc, n, k)

    # estimate area colonization/loss rates
    alλ1, alλ0 = lin_rates(Xt, arav, aroc, alλ1, alλ0, λ1, λ0, ω1, ω0, n, k)

    # trait step
    traitsam_1step!(Xt, liav, δt, ωx, σ, n)

    # biogeographic step
    Ytn = copy(Yt)
    nch = biogeosam_1step(alλ1, alλ0, δt, Ytn, n, k)

    while check_sam(Ytn, nch, n, k)
      Ytn = copy(Yt)
      nch = biogeosam_1step(alλ1, alλ0, δt, Ytn, n, k)
    end

    Yt = Ytn
  
  end

  return Yt
end





"""
    f_λ(λ::Float64, ω::Float64, Δx::Float64)

Estimate rates for area colonization/loss based 
on the difference between lineage traits and area averages.
"""
f_λ(λ::Float64, ω::Float64, Δx::Float64) = @fastmath λ * exp(ω*Δx)





"""
    lin_rates(Xt  ::Array{Float64,1}, arav::Array{Float64,1}, alλ1::Array{Float64,2}, alλ0::Array{Float64,2}, λ1::Float64, λ0::Float64, ω1::Float64, ω0::Float64, n::Int64, k::Int64)

Estimate lineage specific area colonization/loss rates based on differences between each each lineages trait and area averages.
"""
function lin_rates(Xt  ::Array{Float64,1}, 
                   arav::Array{Float64,1}, 
                   aroc::Array{Int64,1},
                   alλ1::Array{Float64,2},
                   alλ0::Array{Float64,2},
                   λ1  ::Float64,
                   λ0  ::Float64,
                   ω1  ::Float64,
                   ω0  ::Float64,
                   n   ::Int64,
                   k   ::Int64)
  @inbounds begin

    for j in Base.OneTo(k)
      
      if aroc[j] < 1 

        for i in Base.OneTo(n)
          alλ1[i,j] = λ1
          alλ0[i,j] = λ0
        end

      else

        for i in Base.OneTo(n)
          Δx = abs(Xt[i] - arav[j])
          alλ1[i,j] = f_λ(λ1, ω1, Δx)
          alλ0[i,j] = f_λ(λ0, ω0, Δx)
        end

      end
    end

    return alλ1, alλ0
  end
end





"""
    ar_lin_avg(Xt::Array{Float64,1}, Yt::Array{Int64,2}, arav::Array{Float64,1}, 
              liav::Array{Float64,1}, n::Int64, k::Int64)

Estimate area and lineage specific averages given sympatry configuration.
"""
function ar_lin_avg(Xt  ::Array{Float64,1}, 
                    Yt  ::Array{Int64,2}, 
                    arav::Array{Float64,1},
                    liav::Array{Float64,1},
                    aroc::Array{Int64,1},
                    n   ::Int64,
                    k   ::Int64)

  @inbounds begin
    # estimate area averages
    for j in Base.OneTo(k)
      arav[j] = 0.0
      na      = 0.0
      aroc[j] = 0
      for i in Base.OneTo(n)
        if Yt[i,j] == 1
          arav[j] += Xt[i]
          na      += 1.0
          aroc[j]  = 1
        end
      end

      arav[j] /= (na == 0.0 ? 1.0 : na) 
    end

    # estimate lineage averages
    for i in Base.OneTo(n)
      liav[i] = 0.0
      na      = 0.0 
      for j in Base.OneTo(k)
        if Yt[i,j] == 1
          liav[i] += arav[j]
          na      += 1.0
        end
      end
      
      liav[i] /= na
    end

  end

  return arav, liav, aroc
end





"""
    traitsam_1step((Xt::Array{Float64,1}, μ ::Array{Float64,1}, δt::Float64, ωx::Float64, σ::Float64, n::Int64)

Sample one step for trait evolution history: X(t + δt).
"""
function traitsam_1step!(Xt::Array{Float64,1}, 
                         μ ::Array{Float64,1}, 
                         δt::Float64, 
                         ωx::Float64, 
                         σ ::Float64,
                         n ::Int64)

  @inbounds @fastmath begin
    
    for i in Base.OneTo(n)
      Xt[i] += E_sde(Xt[i], μ[i], ωx, δt) + randn()*σ*sqrt(δt)
    end
  
  end
end





"""
    biogeosam_1step(λ1::Float64, λ0::Float64, δt::Float64, v1::Array{Int64,1})

Sample one step for biogeographic history Y(t + δt).
"""
function biogeosam_1step(alλ1 ::Array{Float64,2}, 
                         alλ0 ::Array{Float64,2}, 
                         δt   ::Float64, 
                         Ytn  ::Array{Int64,2},
                         n    ::Int64,
                         k    ::Int64)

  @inbounds begin

    nch = zeros(Int64, n)
    for j in Base.OneTo(k), i in Base.OneTo(n)
      if Ytn[i,j] == 0
        if rand() < alλ1[i,j]*δt
          setindex!(Ytn,1,i,j)
          nch[i] += 1
        end
      else 
        if rand() < alλ0[i,j]*δt
          setindex!(Ytn,0,i,j)
          nch[i] += 1 
        end
      end
    end

  end

  return nch::Array{Int64,1}
end





"""
    check_sam(Ytc::Array{Int64,2}, nch::Array{Int64,1}, n::Int64, k::Int64)

Returns false if biogeographic step consists of *only one* change 
or if the species does *not* go globally extinct. 
"""
function check_sam(Ytc::Array{Int64,2}, nch::Array{Int64,1}, n::Int64, k::Int64)

  @fastmath @inbounds begin

    # check lineage did not go extinct
    for i in Base.OneTo(n)
      s = 0
      for j in Base.OneTo(k)
        s += Ytc[i,j]
      end
      if s == 0
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









