"""
Proposal Functions for joint
Biogeographic competition model


Ignacio Quintero

t(-_-t)

May 16 2017
"""




# update node and incident branches
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

    createhists!(λ, Y, pr, d1, d2, brs, brl, brδt, bridx_a, narea, nedge)

    # save extinct
    while ifextY(Y,  triad, narea, bridx_a)
      createhists!(λ, Y, pr, d1, d2, brs, brl, brδt, bridx_a, narea, nedge)
    end

  end
end



#=

Remove multiple events per unit time

=#


# returns true if at some point the species
# goes extinct
function ifextY(Y      ::Array{Int64,3},
                triad  ::Array{Int64,1},
                narea  ::Int64,
                bridx_a::Array{Array{Array{Int64,1},1},1})

  @inbounds begin

    for k=triad
      for i=eachindex(Y[bridx_a[1][k]]),
        s::Int64 = 0
        for j=Base.OneTo(narea)
          s += Y[bridx_a[j][k]][i]::Int64
        end
        if s == 0
          return true
        end
      end
    end

  end

  return false
end





# create and assign to Yc the discrete DA histories 
function createhists!(λ      ::Array{Float64,2}, 
                      Y      ::Array{Int64,3},
                      pr     ::Int64,
                      d1     ::Int64,
                      d2     ::Int64,
                      brs    ::Array{Int64,3},
                      brl    ::Vector{Float64},
                      brδt   ::Vector{Vector{Float64}},
                      bridx_a::Vector{Vector{Vector{Int64}}},
                      narea  ::Int64,
                      nedge  ::Int64)

  @inbounds begin

    for j = Base.OneTo(narea)

      prs::Int64 = brs[pr,2,j]
      setindex!(Y, prs, bridx_a[j][pr][end])

      λj1::Float64 = λ[j,1]
      λj2::Float64 = λ[j,2]

      # sample trio branches event times
      cspr::Array{Float64,1} = 
        rejsam_cumsum(brs[pr,1,j], prs, λj1, λj2, brl[pr]) 
      csd1::Array{Float64,1} = 
        rejsam_cumsum(prs, brs[d1,2,j], λj1, λj2, brl[d1]) 
      csd2::Array{Float64,1} = 
        rejsam_cumsum(prs, brs[d2,2,j], λj1, λj2, brl[d2])

      # discretize such events into current δt
      if pr < nedge 
        assigndisceve!(brs[pr,1,j], Y, cspr, bridx_a[j][pr], brδt[pr])
      end
      assigndisceve!(prs, Y, csd2, bridx_a[j][d1], brδt[d1])
      assigndisceve!(prs, Y, csd2, bridx_a[j][d2], brδt[d2])
    end

  end
end




# assigns discrete values according to 
# the continuous sampling to Yc
function assigndisceve!(si     ::Int64, 
                        Y      ::Array{Int64,3}, 
                        contsam::Array{Float64,1}, 
                        bridx  ::Array{Int64,1}, 
                        δtvec  ::Array{Float64,1})
  s    ::Int64 = 1
  cur_s::Int64 = si
 
  @inbounds begin

    lbr = endof(bridx)
    
    for i=eachindex(contsam)
      f = indmindif_sorted(δtvec, contsam[i])
      setindex!(Y, cur_s, bridx[s:f]) 
      cur_s = 1 - cur_s
      s     = f == lbr ? f : (f + 1)
    end

  end
end




# sample one internal node according to 
# mutual-independent independence model
# transition probabilities
function samplenode!(λ    ::Array{Float64,2},
                     pr   ::Int64,
                     d1   ::Int64,
                     d2   ::Int64,
                     brs  ::Array{Int64,3},
                     brl  ::Array{Float64,1},
                     narea::Int64)
  @inbounds begin
    
    blr_pr::Float64 = brl[pr]
    blr_d1::Float64 = brl[d1]
    blr_d2::Float64 = brl[d2]

    for j=Base.OneTo(narea)
      # transition probabilities for the trio
      ppr_1, ppr_2 = 
        Ptrfast_start(λ[j,1], λ[j,2], blr_pr, brs[pr,1,j])
      pd1_1, pd1_2 = 
        Ptrfast_end(  λ[j,1], λ[j,2], blr_d1, brs[d1,2,j])
      pd2_1, pd2_2 = 
        Ptrfast_end(  λ[j,1], λ[j,2], blr_d2, brs[d2,2,j])

      # normalize probability
      tp::Float64 = normlize(*(ppr_1, pd1_1, pd2_1),
                             *(ppr_2, pd1_2, pd2_2))

      # sample the node's character
      brs[pr,2,j] = brs[d1,1,j] = brs[d2,1,j] = coinsamp(tp)::Int64
    end
  
  end
end




# update stem branch
function upstem(λ   ::Array{Float64,2}, 
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




