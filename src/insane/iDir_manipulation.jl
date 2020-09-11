#=

iDir manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 07 2020
=#


"""
    addda!(id::iDir)

Add `1` to number of data augmented branches. 
"""
addda!(id::iDir) = id.da[] += 1




"""
    rmda!(id::iDir)

Substract `1` to number of data augmented branches. 
"""
rmda!(id::iDir) = id.da[] -= 1




"""
    randbranch(th   ::Float64,
               sth  ::Float64,
               idv  ::Array{iDir,1},
               wbr  ::BitArray{1})

Randomly select branch to graft a subtree with height `sth` onto tree 
with height `th`.
"""
function randbranch(th   ::Float64,
                    sth  ::Float64,
                    idv  ::Array{iDir,1},
                    wbr  ::BitArray{1})

  # uniformly set the height at which to insert `stree`
  h = sth + rand()*(th - sth)

  # get which branches are cut by `h` and sample one at random
  nb  = branchescut!(wbr, h, idv)
  rn  = rand(Base.OneTo(nb))
  br, i  = getbranch(wbr, rn, idv)

  return h, br, i, nb
end




"""
    getbranch(wbr::BitArray{1}, rn::Int64, idv::Array{iDir,1})

Sample one branch among those in `wbr` that are `true` according to 
randomly picked `rn`.
"""
function getbranch(wbr::BitArray{1}, 
                   rn ::Int64, 
                   idv::Array{iDir,1})::Tuple{iDir,Int64}

  i::Int64 = 0
  n::Int64 = 0
  for bit in wbr
    i += 1
    if bit
      n += 1
      if n == rn
        return (idv[i], i)::Tuple{iDir,Int64}
      end
    end
  end
end




"""
    branchescut!(wbr::BitArray{1}, 
                 h  ::Float64, 
                 idv::Array{iDir,1}; 
                 n  ::Int64 = 0)::Int64

Fill the bit array with branches that are cut by  `h`.
"""
function branchescut!(wbr::BitArray{1}, 
                      h  ::Float64, 
                      idv::Array{iDir,1})::Int64
  n::Int64 = 0
  i::Int64 = 0
  for b in idv
    i += 1
    if ti(b) > h > tf(b)
      wbr[i] = true
      n += 1
    else
      wbr[i] = false
    end
  end

  return n
end




"""
    make_triads(idv::Array{iDir, 1})

Return parent and two daughter indices for a given `iDir` vector and a vector
that is true for which ever daughter is a tip.
"""
function make_triads(idv::Array{iDir, 1})

  triads   = Array{Int64,1}[]
  terminus = BitArray{1}[]

  for i in 1:length(idv)
    pr  = i
    drpr = dr(idv[i])
    drd1 = push!(copy(drpr), true)
    drd2 = push!(copy(drpr), false)

    # d1
    d1 = findfirst(x -> dr(x) == drd1, idv)
    # d2
    d2 = findfirst(x -> dr(x) == drd2, idv) 

    if !isnothing(d1) && !isnothing(d2)
      push!(triads, [pr, d1, d2])

      bit = BitArray{1}([false, false])

      # check if either of the tips are terminal
      drd11 = push!(copy(drd1), true)
      drd12 = push!(copy(drd1), false)
      # d1
      d11 = findfirst(x -> dr(x) == drd11, idv)
      # d2
      d12 = findfirst(x -> dr(x) == drd12, idv) 

      if isnothing(d11) && isnothing(d12) 
        bit[1] = true
      end

      drd21 = push!(copy(drd2), true)
      drd22 = push!(copy(drd2), false)
      # d1
      d21 = findfirst(x -> dr(x) == drd21, idv)
      # d2
      d22 = findfirst(x -> dr(x) == drd22, idv) 

      if isnothing(d21) && isnothing(d22) 
        bit[2] = true
      end
      push!(terminus, bit)
    end
  end

  return triads, terminus
end


