#=

Proposal functions for joint
Biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

May 16 2017

=#





#=
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
`Y` IID proposal functions
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#





"""
    upnode!(λ1     ::Float64,
            λ0     ::Float64,
            triad  ::Array{Int64,1},
            Y      ::Array{Int64,3},
            stemevc::Array{Array{Float64,1},1},
            bridx_a::Array{Array{UnitRange{Int64},1},1},
            brδt   ::Vector{Vector{Float64}},
            brl    ::Vector{Float64},
            brs    ::Array{Int64,3},
            narea  ::Int64,
            nedge  ::Int64)

Update node and incident branches using discrete 
Data Augmentation for all areas using a non-competitive 
mutual-independence markov model.
"""
function upnode!(λ1     ::Float64,
                 λ0     ::Float64,
                 triad  ::Array{Int64,1},
                 Y      ::Array{Int64,3},
                 stemevc::Array{Array{Float64,1},1},
                 bridx_a::Array{Array{UnitRange{Int64},1},1},
                 brδt   ::Vector{Vector{Float64}},
                 brl    ::Vector{Float64},
                 brs    ::Array{Int64,3},
                 narea  ::Int64,
                 nedge  ::Int64)

  @inbounds begin
   
    # define branch triad
    pr, d1, d2 = triad

    # sample
    samplenode!(λ1, λ0, pr, d1, d2, brs, brl, narea)

    # save extinct
    while sum(brs[pr,2,:]) < 1
      samplenode!(λ1, λ0, pr, d1, d2, brs, brl, narea)
    end

    # set new node in Y
    for k in Base.OneTo(narea)
      Y[bridx_a[k][pr][end]] = 
      Y[bridx_a[k][d1][1]] = 
      Y[bridx_a[k][d2][1]] = brs[pr,2,k]
    end

    # sample a consistent history
    createhists!(λ1, λ0, Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge,
                 stemevc, brl[nedge])

    while ifextY(Y, triad, narea, bridx_a)
      createhists!(λ1, λ0, Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge,
                   stemevc, brl[nedge])
    end
  end

  return nothing
end




"""
    createhists!(λ::Array{Float64,1}, Y::Array{Int64,3}, pr::Int64, d1::Int64, d2::Int64, brs::Array{Int64,3}, brδt::Array{Array{Float64,1},1}, bridx_a::Array{Array{Array{Int64,1},1},1}, narea::Int64)

Create bit histories for all areas for the branch trio.
"""
function createhists!(λ1     ::Float64,
                      λ0     ::Float64,
                      Y      ::Array{Int64,3},
                      pr     ::Int64,
                      d1     ::Int64,
                      d2     ::Int64,
                      brs    ::Array{Int64,3},
                      brδt   ::Array{Array{Float64,1},1},
                      bridx_a::Array{Array{UnitRange{Int64},1},1},
                      narea  ::Int64,
                      nedge  ::Int64,
                      stemevc::Array{Array{Float64,1},1},
                      stbrl  ::Float64)

  @inbounds begin

    # if stem branch do continuous DA
    if pr == nedge
      br_samp!(stemevc, brs[pr,1,:], brs[pr,2,:], λ1, λ0, stbrl, narea)

      for j = Base.OneTo(narea), idx = (d1,d2)
        bit_rejsam!(Y, bridx_a[j][idx], brs[idx,2,j], 
                    λ1, λ0, brδt[idx])
      end
    else

      for j = Base.OneTo(narea), idx = (pr,d1,d2)
        bit_rejsam!(Y, bridx_a[j][idx], brs[idx,2,j], 
                    λ1, λ0, brδt[idx])
      end
    end
  end

  return nothing
end





"""
    samplenode!(λ::Array{Float64,1}, pr::Int64, d1::Int64, d2::Int64, brs::Array{Int64,3}, brl::Array{Float64,1}, narea::Int64)

Sample one internal node according to 
mutual-independence model transition probabilities.
"""
function samplenode!(λ1   ::Float64, 
                     λ0   ::Float64,
                     pr   ::Int64,
                     d1   ::Int64,
                     d2   ::Int64,
                     brs  ::Array{Int64,3},
                     brl  ::Array{Float64,1},
                     narea::Int64)
  @inbounds begin

    for j = Base.OneTo(narea)

      # transition probabilities for the trio
      ppr_1, ppr_2 = 
        Ptrfast_start(λ1, λ0, brl[pr], brs[pr,1,j])
      pd1_1, pd1_2 = 
        Ptrfast_end(  λ1, λ0, brl[d1], brs[d1,2,j])
      pd2_1, pd2_2 = 
        Ptrfast_end(  λ1, λ0, brl[d2], brs[d2,2,j])

      # normalize probability
      tp = normlize(*(ppr_1, pd1_1, pd2_1),
                    *(ppr_2, pd1_2, pd2_2))::Float64

      # sample the node's character
      brs[pr,2,j] = brs[d1,1,j] = brs[d2,1,j] = coinsamp(tp)::Int64
    end
  end

  return nothing
end





#=
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
`Y` with `ΔX` proposal functions
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#





"""
    upnode!(λ1     ::Float64,
            λ0     ::Float64,
            ω1     ::Float64,
            ω0     ::Float64,
            avg_Δx ::Array{Float64,2},
            triad  ::Array{Int64,1},
            Y      ::Array{Int64,3},
            bridx_a::Array{Array{UnitRange{Int64},1},1},
            brδt   ::Vector{Vector{Float64}},
            brl    ::Vector{Float64},
            brs    ::Array{Int64,3},
            narea  ::Int64,
            nedge  ::Int64)

Update node and incident branches using discrete 
Data Augmentation for all areas with mutual-independence 
proposals taking into account `Δx` and `ω1` & `ω0`.
"""
function upnode!(λ1     ::Float64,
                 λ0     ::Float64,
                 ω1     ::Float64,
                 ω0     ::Float64,
                 avg_Δx ::Array{Float64,2},
                 triad  ::Array{Int64,1},
                 Y      ::Array{Int64,3},
                 stemevc::Array{Array{Float64,1},1},
                 bridx_a::Array{Array{UnitRange{Int64},1},1},
                 brδt   ::Vector{Vector{Float64}},
                 brl    ::Vector{Float64},
                 brs    ::Array{Int64,3},
                 narea  ::Int64,
                 nedge  ::Int64)

  @inbounds begin
   
    # define branch triad
    pr, d1, d2 = triad

    # sample
    samplenode!(λ1, λ0, ω1, ω0, avg_Δx, pr, d1, d2, brs, brl, narea)

    # save extinct
    while sum(brs[pr,2,:]) < 1
       samplenode!(λ1, λ0, ω1, ω0, avg_Δx, pr, d1, d2, brs, brl, narea)
    end

    # set new node in Y
    for k = Base.OneTo(narea)
      Y[bridx_a[k][pr][end]] = 
      Y[bridx_a[k][d1][1]]   = 
      Y[bridx_a[k][d2][1]]   = brs[pr,2,k]
    end

    # sample a consistent history
    createhists!(λ1, λ0, ω1, ω0, avg_Δx, 
                 Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge,
                 stemevc, brl[nedge])

    ntries = 1
    # save extinct
    while ifextY(Y, triad, narea, bridx_a)
      createhists!(λ1, λ0, ω1, ω0, avg_Δx, 
                   Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge,
                   stemevc, brl[nedge])

      ntries += 1
      if ntries == 50_000
        return true
      end

    end
  end

  return nothing
end





"""
    ifextY(Y::Array{Int64,3}, triad::Array{Int64,1}, narea::Int64, bridx_a::Array{Array{Array{Int64,1},1},1})

Return `true` if at some point the species
goes extinct and/or more than one change is 
observed after some **δt**, otherwise returns `false`.
"""
function ifextY(Y      ::Array{Int64,3},
                triad  ::Array{Int64,1},
                narea  ::Int64,
                bridx_a::Array{Array{UnitRange{Int64},1},1})

  @inbounds begin

    for k in triad

      for i = Base.OneTo((length(bridx_a[1][k]::UnitRange{Int64})-1)::Int64)
        s_e = 0::Int64            # count current areas
        s_c = 0::Int64            # count area changes

        for j = Base.OneTo(narea)
          s_e += Y[bridx_a[j][k][i]]::Int64
          if Y[bridx_a[j][k][i]]::Int64 != Y[bridx_a[j][k][i+1]]::Int64
            s_c += 1::Int64
          end
        end

        if s_e::Int64 == 0::Int64 || s_c::Int64 > 1::Int64
          return true::Bool
        end
      end
    end

  end

  return false::Bool
end





"""
    createhists!(λ1     ::Float64,
                 λ0     ::Float64,
                 ω1     ::Float64,  
                 ω0     ::Float64, 
                 avg_Δx ::Array{Float64,2},
                 Y      ::Array{Int64,3},
                 pr     ::Int64,
                 d1     ::Int64,
                 d2     ::Int64,
                 brs    ::Array{Int64,3},
                 brδt   ::Array{Array{Float64,1},1},
                 bridx_a::Array{Array{UnitRange{Int64},1},1},
                 narea  ::Int64,
                 nedge  ::Int64)

Create bit histories for all areas for the branch 
trio taking into account `Δx` and `ω1` & `ω0`.
"""
function createhists!(λ1     ::Float64,
                      λ0     ::Float64,
                      ω1     ::Float64,  
                      ω0     ::Float64, 
                      avg_Δx ::Array{Float64,2},
                      Y      ::Array{Int64,3},
                      pr     ::Int64,
                      d1     ::Int64,
                      d2     ::Int64,
                      brs    ::Array{Int64,3},
                      brδt   ::Array{Array{Float64,1},1},
                      bridx_a::Array{Array{UnitRange{Int64},1},1},
                      narea  ::Int64,
                      nedge  ::Int64, 
                      stemevc::Array{Array{Float64,1},1},
                      stbrl  ::Float64)

  @inbounds begin
    if pr == nedge
      # if stem branch do continuous DA
      br_samp!(stemevc, brs[pr,1,:], brs[pr,2,:], λ1, λ0, stbrl, narea)

      for j = Base.OneTo(narea), idx = (d1,d2)
        bit_rejsam!(Y, bridx_a[j][idx], brs[idx,2,j], 
                    λ1, λ0, ω1, ω0, avg_Δx[idx,j], brδt[idx])
      end
    else

      for j = Base.OneTo(narea), idx = (pr,d1,d2)
        bit_rejsam!(Y, bridx_a[j][idx], brs[idx,2,j], 
                    λ1, λ0, ω1, ω0, avg_Δx[idx,j], brδt[idx])
      end
    end
  end

  return nothing
end




"""
    samplenode!(λ::Array{Float64,1}, pr::Int64, d1::Int64, d2::Int64, brs::Array{Int64,3}, brl::Array{Float64,1}, narea::Int64)

Sample one internal node according to 
mutual-independence model transition probabilities
taking into account `Δx` and `ω1` & `ω0`.
"""
function samplenode!(λ1    ::Float64, 
                     λ0    ::Float64,
                     ω1    ::Float64,
                     ω0    ::Float64,
                     avg_Δx::Array{Float64,2},
                     pr    ::Int64,
                     d1    ::Int64,
                     d2    ::Int64,
                     brs   ::Array{Int64,3},
                     brl   ::Array{Float64,1},
                     narea ::Int64)
  @inbounds begin

    for j = Base.OneTo(narea)

      # transition probabilities for the trio
      ppr_1, ppr_2 = 
        Ptrfast_start(λ1, λ0, ω1, ω0, avg_Δx[pr,j], brl[pr], brs[pr,1,j])
      pd1_1, pd1_2 = 
        Ptrfast_end(  λ1, λ0, ω1, ω0, avg_Δx[d1,j], brl[d1], brs[d1,2,j])
      pd2_1, pd2_2 = 
        Ptrfast_end(  λ1, λ0, ω1, ω0, avg_Δx[d2,j], brl[d2], brs[d2,2,j])

      # normalize probabilitye
      tp = normlize(*(ppr_1, pd1_1, pd2_1),
                    *(ppr_2, pd1_2, pd2_2))::Float64

      # sample the node's character
      brs[pr,2,j] = brs[d1,1,j] = brs[d2,1,j] = coinsamp(tp)::Int64
    end
  end

  return nothing
end





#=
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Y branch proposal functions
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#





"""
upbranchY!(λ1     ::Float64,
           λ0     ::Float64,
           ω1     ::Float64,
           ω0     ::Float64,
           avg_Δx ::Array{Float64,2},
           br     ::Int64,
           Y      ::Array{Int64,3},
           stemevc::Array{Array{Float64,1},1},
           wareas ::Array{Int64,1},
           bridx_a::Array{Array{UnitRange{Int64},1},1},
           brδt   ::Vector{Vector{Float64}},
           brl    ::Vector{Float64},
           brs    ::Array{Int64,3},
           narea  ::Int64,
           nedge  ::Int64)

Update one branch using discrete Data Augmentation 
for all areas with independent 
proposals taking into account `Δx` and `ω1` & `ω0`.
"""
function upbranchY!(λ1     ::Float64,
                    λ0     ::Float64,
                    ω1     ::Float64,
                    ω0     ::Float64,
                    avg_Δx ::Array{Float64,2},
                    br     ::Int64,
                    Y      ::Array{Int64,3},
                    stemevc::Array{Array{Float64,1},1},
                    bridx_a::Array{Array{UnitRange{Int64},1},1},
                    brδt   ::Vector{Vector{Float64}},
                    brl    ::Vector{Float64},
                    brs    ::Array{Int64,3},
                    narea  ::Int64,
                    nedge  ::Int64)

  # if stem branch
  if br == nedge
    br_samp!(stemevc, brs[nedge,1,:], brs[nedge,2,:], λ1, λ0, 
             brl[nedge], narea)
  else 
    # sample a consistent history
    createhists!(λ1, λ0, ω1, ω0, avg_Δx, 
                 Y, br, brs, brδt, bridx_a, narea)

    ntries = 1
    # check if extinct
    while ifextY(Y, br, narea, bridx_a)
      createhists!(λ1, λ0, ω1, ω0, avg_Δx, 
                   Y, br, brs, brδt, bridx_a, narea)
      ntries += 1
      if ntries == 50_000
        return true
      end
    end
  end

  return nothing
end





"""
    createhists!(λ1     ::Float64,
                 λ0     ::Float64,
                 ω1     ::Float64,  
                 ω0     ::Float64, 
                 avg_Δx ::Array{Float64,2},
                 Y      ::Array{Int64,3},
                 br     ::Int64,
                 brs    ::Array{Int64,3},
                 brδt   ::Array{Array{Float64,1},1},
                 bridx_a::Array{Array{UnitRange{Int64},1},1},
                 narea  ::Int64)

Create bit histories for all areas for one single branch 
taking into account `Δx` and `ω1` & `ω0` for all areas.
"""
function createhists!(λ1     ::Float64,
                      λ0     ::Float64,
                      ω1     ::Float64,  
                      ω0     ::Float64, 
                      avg_Δx ::Array{Float64,2},
                      Y      ::Array{Int64,3},
                      br     ::Int64,
                      brs    ::Array{Int64,3},
                      brδt   ::Array{Array{Float64,1},1},
                      bridx_a::Array{Array{UnitRange{Int64},1},1},
                      narea  ::Int64)

  @inbounds begin
    for j = Base.OneTo(narea)
      bit_rejsam!(Y, bridx_a[j][br], brs[br,2,j], 
                  λ1, λ0, ω1, ω0, avg_Δx[br,j], brδt[br])
    end
  end

  return nothing
end





"""
    ifextY(Y      ::Array{Int64,3},
           br     ::Int64,
           narea  ::Int64,
           bridx_a::Array{Array{UnitRange{Int64},1},1})

Return `true` if at some point the species
goes extinct and/or more than one change is 
observed after some **δt**, otherwise returns `false`. 
This specific method is for single branch updates.
"""
function ifextY(Y      ::Array{Int64,3},
                br     ::Int64,
                narea  ::Int64,
                bridx_a::Array{Array{UnitRange{Int64},1},1})

  @inbounds begin

    for i = Base.OneTo((length(bridx_a[1][br]::UnitRange{Int64})-1)::Int64)
      s_e = 0::Int64            # count current areas
      s_c = 0::Int64            # count area changes

      for j = Base.OneTo(narea)
        s_e += Y[bridx_a[j][br][i]]::Int64
        if Y[bridx_a[j][br][i]]::Int64 != Y[bridx_a[j][br][i+1]]::Int64
          s_c += 1::Int64
        end
      end

      if s_e::Int64 == 0::Int64 || s_c::Int64 > 1::Int64
        return true::Bool
      end
    end

  end

  return false::Bool
end





#=
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
X proposal functions
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#





"""
    upbranchX!(j    ::Int64, 
               X    ::Array{Float64,2}, 
               bridx::Array{UnitRange{Int64},1},
               brδt ::Array{Array{Float64,1},1}, 
               σ²c  ::Float64)

Update a branch j in X using a Brownian bridge.
"""
function upbranchX!(j    ::Int64, 
                    X    ::Array{Float64,2}, 
                    bridx::Array{UnitRange{Int64},1},
                    brδt ::Array{Array{Float64,1},1},
                    s2   ::Float64)

  @inbounds bbX!(X, bridx[j], brδt[j], s2)

  return nothing
end





"""
    bbX!(x::Array{Float64,1}, t::::Array{Float64,1}, σ::Float64)

Brownian bridge simulation function for updating a branch in X in place.
"""
function bbX!(X  ::Array{Float64,2}, 
              idx::UnitRange,
              t  ::Array{Float64,1},
              s2 ::Float64)

  @inbounds begin

    xf::Float64 = X[idx[end]]

    for i = Base.OneTo(endof(t)-1)
      X[idx[i+1]] = (X[idx[i]] + randn()*sqrt((t[i+1] - t[i])*s2))::Float64
    end

    for i = Base.OneTo(endof(t))
      X[idx[i]] = (X[idx[i]] - t[i]/t[end] * (X[idx[end]] - xf))::Float64
    end
  end

  return nothing
end




