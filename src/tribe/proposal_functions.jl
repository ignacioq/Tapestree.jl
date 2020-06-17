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
            stemevs::Array{Array{Float64,1},1},
            bridx_a::Array{Array{UnitRange{Int64},1},1},
            brδt   ::Vector{Vector{Float64}},
            brl    ::Vector{Float64},
            brs    ::Array{Int64,3},
            narea  ::Int64,
            nedge  ::Int64)

Update node and incident branches using discrete 
Data Augmentation for all areas using a non-competitive 
mutual-independence Markov model.
"""
function upnode!(λ1     ::Float64,
                 λ0     ::Float64,
                 triad  ::Array{Int64,1},
                 Y      ::Array{Int64,3},
                 stemevs::Array{Array{Float64,1},1},
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
    ntries = 1
    while iszero(sum(view(brs,pr,2,:)))
      samplenode!(λ1, λ0, pr, d1, d2, brs, brl, narea)

      ntries += 1
      if ntries == 500
        return false
      end
    end

    # set new node in Y
    @simd for k in Base.OneTo(narea)
      Y[bridx_a[k][d1][1]] = Y[bridx_a[k][d2][1]] = brs[pr,2,k]
    end

    # sample a consistent history
    createhists!(λ1, λ0, Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge,
                 stemevs, brl[nedge])

    ntries = 1
    while ifextY(Y, stemevs, triad, brs, brl[nedge], narea, bridx_a, nedge)
      createhists!(λ1, λ0, Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge,
                   stemevs, brl[nedge])

      ntries += 1
      if ntries == 500
        return false
      end
    end

  end

  return true
end






"""
    samplenode!(λ1   ::Float64, 
                λ0   ::Float64,
                pr   ::Int64,
                d1   ::Int64,
                d2   ::Int64,
                brs  ::Array{Int64,3},
                brl  ::Array{Float64,1},
                narea::Int64)

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

    # estimate transition probabilities
    pr0_1, pr0_2 = Ptrfast_start(λ1, λ0, brl[pr], Val{0})
    pr1_1, pr1_2 = Ptrfast_start(λ1, λ0, brl[pr], Val{1})
    d10_1, d10_2 = Ptrfast_end(  λ1, λ0, brl[d1], Val{0})
    d11_1, d11_2 = Ptrfast_end(  λ1, λ0, brl[d1], Val{1})
    d20_1, d20_2 = Ptrfast_end(  λ1, λ0, brl[d2], Val{0})
    d21_1, d21_2 = Ptrfast_end(  λ1, λ0, brl[d2], Val{1})

    for k = Base.OneTo(narea)

      if iszero(brs[pr,1,k])
        ppr_1, ppr_2 = pr0_1, pr0_2
      else 
        ppr_1, ppr_2 = pr1_1, pr1_2
      end

      if iszero(brs[d1,2,k])
        pd1_1, pd1_2 = d10_1, d10_2
      else 
        pd1_1, pd1_2 = d11_1, d11_2
      end

      if iszero(brs[d2,2,k])
        pd2_1, pd2_2 = d20_1, d20_2
      else 
        pd2_1, pd2_2 = d21_1, d21_2
      end

      tp = normlize(*(ppr_1, pd1_1, pd2_1),
                    *(ppr_2, pd1_2, pd2_2))::Float64

      # sample the node's character
      brs[pr,2,k] = brs[d1,1,k] = brs[d2,1,k] = coinsamp(tp)::Int64
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
                      stemevs::Array{Array{Float64,1},1},
                      stbrl  ::Float64)

  @inbounds begin

    if pr == nedge
      # if stem branch do continuous DA
      mult_rejsam!(stemevs, brs, λ1, λ0, stbrl, narea, nedge)
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
    mult_rejsam!(evs  ::Array{Array{Float64,1},1},
                 brs  ::Array{Int64,3}, 
                 λ1   ::Float64,
                 λ0   ::Float64,
                 t    ::Float64,
                 narea::Int64,
                 nedge::Int64)

  Multi-area branch rejection independent model sampling.
"""
function mult_rejsam!(evs  ::Array{Array{Float64,1},1},
                      brs  ::Array{Int64,3}, 
                      λ1   ::Float64,
                      λ0   ::Float64,
                      t    ::Float64,
                      narea::Int64,
                      nedge::Int64)

  @simd for k = Base.OneTo(narea)
    rejsam!(evs[k], brs[nedge,1,k], brs[nedge,2,k], λ1, λ0, t)
  end

  return nothing
end





"""
    ifextY(Y      ::Array{Int64,3},
           triad  ::Array{Int64,1},
           narea  ::Int64,
           bridx_a::Array{Array{UnitRange{Int64},1},1})

Return `true` if at some point the species
goes extinct and/or more than one change is 
observed after some **δt**, otherwise returns `false`.
"""
function ifextY(Y      ::Array{Int64,3},
                stemevs::Array{Array{Float64,1},1},
                triad  ::Array{Int64,1},
                brs    ::Array{Int64,3},
                stbrl  ::Float64,
                narea  ::Int64,
                bridx_a::Array{Array{UnitRange{Int64},1},1},
                nedge  ::Int64)

  @inbounds begin

    if triad[1] == nedge
      ifext_cont(stemevs, brs, stbrl, narea, nedge) && return true::Bool

      for k in (triad[2],triad[3])
        ifext_disc(Y, k, narea, bridx_a) && return true::Bool
      end
    else 

      for k in triad
        ifext_disc(Y, k, narea, bridx_a) && return true::Bool
      end
    end
  end

  return false::Bool
end





"""
    ifext_cont(t_hist::Array{Array{Float64,1},1},
               brs   ::Array{Int64,3}, 
               t     ::Float64,
               narea ::Int64,
               nedge ::Int64)

Return true if lineage goes extinct.
"""
function ifext_cont(t_hist::Array{Array{Float64,1},1},
                    brs   ::Array{Int64,3}, 
                    t     ::Float64,
                    narea ::Int64,
                    nedge ::Int64)

  @inbounds begin

    ioc = 0
    # initial occupancy time
    for k in Base.OneTo(narea)
      if brs[nedge,1,k] == 1
        ioc = k
        break
      end
    end

    ioct = t_hist[ioc][1]

    ntries = 0
    while !isapprox(ioct, t, atol = 1.0e-12)

      if ioc == narea
        ioc = 1
      else 
        ioc += 1
      end

      tc = 0.0
      cs = brs[nedge,1,ioc]
      for ts in t_hist[ioc]
        tc += ts
        if ioct < tc 
          if cs == 1
            ioct  = tc 
            ntries = 0
            break
          else
            ntries += 1
            if ntries > narea
              return true
            end
            break
          end
        end
        cs = 1 - cs
      end

    end
  end

  return false::Bool
end






"""
    ifext_disc(Y      ::Array{Int64,3},
               br     ::Int64,
               narea  ::Int64,
               bridx_a::Array{Array{UnitRange{Int64},1},1})

Return `true` if at some point the species
goes extinct and/or more than one change is 
observed after some **δt**, otherwise returns `false`. 
This specific method is for single branch updates.
"""
function ifext_disc(Y      ::Array{Int64,3},
                    br     ::Int64,
                    narea  ::Int64,
                    bridx_a::Array{Array{UnitRange{Int64},1},1})

  @inbounds begin

    for i = Base.OneTo(length(bridx_a[1][br]::UnitRange{Int64})-1)
      s_e::Int64 = 0            # count current areas
      s_c::Int64 = 0            # count area changes

      for k = Base.OneTo(narea)
        s_e += Y[bridx_a[k][br][i]]::Int64
        if Y[bridx_a[k][br][i]]::Int64 != Y[bridx_a[k][br][i+1]]::Int64
          s_c += 1
        end
      end

      if s_e == 0 || s_c > 1
        return true::Bool
      end
    end

  end

  return false::Bool
end





#=
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Y stem node proposal function (continuous DA)
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#

"""
    upstemnode!(λ1   ::Float64, 
                λ0   ::Float64,
                nedge::Int64,
                stemevs::Array{Array{Float64,1},1},
                brs  ::Array{Int64,3},
                brl  ::Array{Float64,1},
                narea::Int64)

Update stem node.
"""
function upstemnode!(λ1     ::Float64, 
                     λ0     ::Float64,
                     nedge  ::Int64,
                     stemevs::Array{Array{Float64,1},1},
                     brs    ::Array{Int64,3},
                     stbrl  ::Float64,
                     narea  ::Int64)

  @inbounds begin

    # sample
    samplestem!(λ1, λ0, nedge, brs, stbrl, narea)

    # save extinct
    ntries = 1
    while sum(view(brs,nedge,1,:)) < 1
      samplestem!(λ1, λ0, nedge, brs, stbrl, narea)
      ntries += 1
      if ntries == 500
        return false::Bool
      end
    end

    # sample a congruent history
    mult_rejsam!(stemevs, brs, λ1, λ0, stbrl, narea, nedge)

    ntries = 1
    # check if extinct

    while ifext_cont(stemevs, brs, stbrl, narea, nedge)
      mult_rejsam!(stemevs, brs, λ1, λ0, stbrl, narea, nedge)

      ntries += 1
      if ntries == 500
        return false::Bool
      end
    end
  end

  return true::Bool
end





"""
    samplestem!(λ1   ::Float64, 
                λ0   ::Float64,
                nedge::Int64,
                brs  ::Array{Int64,3},
                brl  ::Array{Float64,1},
                narea::Int64)

Sample stem node.
"""
function samplestem!(λ1   ::Float64, 
                     λ0   ::Float64,
                     nedge::Int64,
                     brs  ::Array{Int64,3},
                     stbrl::Float64,
                     narea::Int64)

 @inbounds begin

    # estimate transition probabilities
    p0 = normlize(Ptrfast_end(λ1, λ0, stbrl, Val{0}))::Float64
    p1 = normlize(Ptrfast_end(λ1, λ0, stbrl, Val{1}))::Float64

    for k = Base.OneTo(narea)
      # sample the node's character
      if iszero(brs[nedge,2,k])
        brs[nedge,1,k] = coinsamp(p0)::Int64
      else 
        brs[nedge,1,k] = coinsamp(p1)::Int64
      end
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
                    br     ::Int64,
                    Y      ::Array{Int64,3},
                    stemevs::Array{Array{Float64,1},1},
                    bridx_a::Array{Array{UnitRange{Int64},1},1},
                    brδt   ::Vector{Vector{Float64}},
                    stbrl  ::Float64,
                    brs    ::Array{Int64,3},
                    narea  ::Int64,
                    nedge  ::Int64)

  ntries = 1

  # if stem branch
  if br == nedge
    mult_rejsam!(stemevs, brs, λ1, λ0, stbrl, narea, nedge)

    # check if extinct
    while ifext_cont(stemevs, brs, stbrl, narea, nedge)
      mult_rejsam!(stemevs, brs, λ1, λ0, stbrl, narea, nedge)

      ntries += 1
      if ntries == 500
        return false
      end
    end

  else
    createhists!(λ1, λ0, Y, br, brs, brδt, bridx_a, narea)

    # check if extinct
    while ifext_disc(Y, br, narea, bridx_a)
      createhists!(λ1, λ0, Y, br, brs, brδt, bridx_a, narea)

      ntries += 1
      if ntries == 500
        return false
      end
    end
  end

  return true
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
                      Y      ::Array{Int64,3},
                      br     ::Int64,
                      brs    ::Array{Int64,3},
                      brδt   ::Array{Array{Float64,1},1},
                      bridx_a::Array{Array{UnitRange{Int64},1},1},
                      narea  ::Int64)
  @inbounds begin
    for j = Base.OneTo(narea)
      bit_rejsam!(Y, bridx_a[j][br], brs[br,2,j], λ1, λ0, brδt[br])
    end
  end

  return nothing
end




#=
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
X proposal functions
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    uptrioX!(pr   ::Int64, 
             d1   ::Int64,
             d2   ::Int64, 
             X    ::Array{Float64,2}, 
             bridx::Array{UnitRange{Int64},1},
             brδt ::Array{Array{Float64,1},1}, 
             σ²c  ::Float64)

Update the node and adjoining branhces of `trio` using Brownian bridges.
"""
function uptrioX!(pr   ::Int64, 
                  d1   ::Int64,
                  d2   ::Int64, 
                  X    ::Array{Float64,2}, 
                  bridx::Array{UnitRange{Int64},1},
                  brδt ::Array{Array{Float64,1},1},
                  brl  ::Array{Float64,1},
                  σ²ϕ  ::Float64, 
                  nedge::Int64)
  @inbounds begin

    # if not root
    if pr != nedge

      ipr = bridx[pr]
      id1 = bridx[d1]
      id2 = bridx[d2]

      # update node
      X[id1[1]] = 
      X[id2[1]] = trioupd(X[ipr[1]], 
                          X[id1[end]], 
                          X[id2[end]],
                          brl[pr], brl[d1], brl[d2], σ²ϕ)

      #update branches
      bbX!(X, ipr, brδt[pr], σ²ϕ)
      bbX!(X, id1, brδt[d1], σ²ϕ)
      bbX!(X, id2, brδt[d2], σ²ϕ)

    else
      id1 = bridx[d1]
      id2 = bridx[d2]

      # update node
      X[id1[1]] = 
      X[id2[1]] = duoupd(X[id1[end]],
                         X[id2[end]], 
                         brl[d1], brl[d2], σ²ϕ)

      # update branches
      bbX!(X, id1, brδt[d1], σ²ϕ)
      bbX!(X, id2, brδt[d2], σ²ϕ)
    end

  end

  return nothing
end





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
                    σ²ϕ  ::Float64)

  @inbounds bbX!(X, bridx[j], brδt[j], σ²ϕ)

  return nothing
end





"""
    bbX!(X::Array{Float64,2}, idx::UnitRange, t::Array{Float64,1}, σ::Float64)

Brownian bridge simulation function for updating a branch in X in place.
"""
function bbX!(X  ::Array{Float64,2}, 
              idx::UnitRange,
              t  ::Array{Float64,1},
              σ²ϕ::Float64)

  @inbounds begin

    xf::Float64 = X[idx[end]]

    for i = Base.OneTo(lastindex(t)-1)
      X[idx[i+1]] = (X[idx[i]] + randn()*sqrt((t[i+1] - t[i])*σ²ϕ))::Float64
    end

    invte::Float64 = 1.0/t[end]
    xdif ::Float64 = (X[idx[end]] - xf)

    @simd for i = Base.OneTo(lastindex(t))
      X[idx[i]] = (X[idx[i]] - t[i] * invte * xdif)::Float64
    end
  end

  return nothing
end



