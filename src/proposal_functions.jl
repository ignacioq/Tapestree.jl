
#=

Proposal functions for joint
Biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

May 16 2017

=#




"""
  upnode!(λ::Array{Float64,2},
          triad::Vector{Int64},
          Y::Array{Int64,3},
          bridx_a::Vector{Vector{Vector{Int64}}},
          brδt::Vector{Vector{Float64}},
          brl::Vector{Float64},
          brs::Array{Int64,3},
          narea::Int64,
          nedge::Int64)

  update node and incident branches
"""
function upnode!(λ      ::Array{Float64,2},
                 triad  ::Vector{Int64},
                 Y      ::Array{Int64,3},
                 bridx_a::Vector{Vector{Vector{Int64}}},
                 brδt   ::Vector{Vector{Float64}},
                 brl    ::Vector{Float64},
                 brs    ::Array{Int64,3},
                 narea  ::Int64,
                 nedge  ::Int64)

  @inbounds begin
   
    # define branch triad
    pr, d1, d2 = triad

    # sample
    samplenode!(λ, pr, d1, d2, brs, brl, narea)

    # save extinct
    while sum(brs[pr,2,:]) == 0
       samplenode!(λ, pr, d1, d2, brs, brl, narea)
    end


    # sample a consistent history
    createhists!(λ, Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge)
  
    # save extinct
    while ifextY(Y,  triad, narea, bridx_a)
      createhists!(λ, Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge)
    end

  end
end


#=
  *************************************
  allow only one change at a time for δt
  *************************************
=#


"""
    ifextY(Y::Array{Int64,3}, triad::Array{Int64,1}, narea::Int64, bridx_a::Array{Array{Array{Int64,1},1},1})

Return `true` if at some point the species
goes extinct and/or more than one change is 
observed after some **δt**, otherwise returns `false`.
"""
function ifextY(Y      ::Array{Int64,3},
                triad  ::Array{Int64,1},
                narea  ::Int64,
                bridx_a::Array{Array{Array{Int64,1},1},1})

  @inbounds begin

    for k=triad, i=2:(length(Y[bridx_a[1][k]])-1),
      s::Int64 = 0
      for j=Base.OneTo(narea)
        s += Y[bridx_a[j][k]][i]::Int64
      end
      if s == 0
        return true
      end
    end
    

  end

  return false
end





"""
    createhists!(λ::Array{Float64,2}, Y::Array{Int64,3}, pr::Int64, d1::Int64, d2::Int64, brs::Array{Int64,3}, brδt::Array{Array{Float64,1},1}, bridx_a::Array{Array{Array{Int64,1},1},1}, narea::Int64)

Create bit histories for all areas for the branch trio.
"""
function createhists!(λ      ::Array{Float64,2}, 
                      Y      ::Array{Int64,3},
                      pr     ::Int64,
                      d1     ::Int64,
                      d2     ::Int64,
                      brs    ::Array{Int64,3},
                      brδt   ::Array{Array{Float64,1},1},
                      bridx_a::Array{Array{Array{Int64,1},1},1},
                      narea  ::Int64,
                      nedge  ::Int64)

  @inbounds begin

    for j = Base.OneTo(narea)

      # set new node in Y
      setindex!(Y, brs[pr,2,j], bridx_a[j][pr][end])
      setindex!(Y, brs[pr,2,j], bridx_a[j][d1][1])
      setindex!(Y, brs[pr,2,j], bridx_a[j][d2][1])

      λj1 = λ[j,1]::Float64
      λj2 = λ[j,2]::Float64

      if pr < nedge
        # for parent branch
        bit_rejsam!(Y, bridx_a[j][pr], brs[pr,2,j], λj1, λj2, brδt[pr])
      end

      # for daughter branch 1
      bit_rejsam!(Y, bridx_a[j][d1], brs[d1,2,j], λj1, λj2, brδt[d1])

      # for daughter branch 2
      bit_rejsam!(Y, bridx_a[j][d2], brs[d2,2,j], λj1, λj2, brδt[d2])

    end

  end
end





"""
  samplenode!(λ::Array{Float64,2},
              pr::Int64,
              d1::Int64,
              d2::Int64,
              brs::Array{Int64,3},
              brl::Array{Float64,1},
              narea::Int64)

  sample one internal node according to 
  mutual-independence model
  transition probabilities
"""
function samplenode!(λ    ::Array{Float64,2},
                     pr   ::Int64,
                     d1   ::Int64,
                     d2   ::Int64,
                     brs  ::Array{Int64,3},
                     brl  ::Array{Float64,1},
                     narea::Int64)
  @inbounds begin
    
    brl_pr = brl[pr]::Float64
    brl_d1 = brl[d1]::Float64 
    brl_d2 = brl[d2]::Float64 

    for j=Base.OneTo(narea)

      # transition probabilities for the trio
      ppr_1, ppr_2 = 
        Ptrfast_start(λ[j,1], λ[j,2], brl_pr, brs[pr,1,j])
      pd1_1, pd1_2 = 
        Ptrfast_end(  λ[j,1], λ[j,2], brl_d1, brs[d1,2,j])
      pd2_1, pd2_2 = 
        Ptrfast_end(  λ[j,1], λ[j,2], brl_d2, brs[d2,2,j])

      # normalize probability
      tp = normlize(*(ppr_1, pd1_1, pd2_1),
                             *(ppr_2, pd1_2, pd2_2))::Float64

      # sample the node's character
      brs[pr,2,j] = brs[d1,1,j] = brs[d2,1,j] = coinsamp(tp)::Int64
    end
  
  end
end





"""
  update stem branch

  *********
  move indexing outside
  *********
"""
function upstem(λ    ::Array{Float64,2}, 
                idx  ::Int64,
                brs  ::Array{Int64,3}, 
                brl  ::Vector{Float64},
                narea::Int64)

  @inbounds begin
    
    for j in Base.OneTo(narea)
      # transition probabilities
      p1::Tuple{Float64,Float64} = 
        Ptrfast_end(λ[j,1], λ[j,2], brl[idx], brs[idx,2,j])

      # sample the stem node character
      brs[idx,1,j] = coinsamp(normlize(p1[1],p1[2]))
    end

    while sum(brs[idx,1,:]) < 1
      for j in Base.OneTo(narea)
        # transition probabilities
        p1::Tuple{Float64,Float64} = 
          Ptrfast_end(λ[j,1], λ[j,2], brl[idx], brs[idx,2,j])

        # sample the stem node character
        brs[idx,1,j] = coinsamp(normlize(p1[1],p1[2]))
      end
    end

  end

  # sample new history for stem branch
  br_samp(brs[idx,1,:], brs[idx,2,:], λ, brl[idx], narea)

end




