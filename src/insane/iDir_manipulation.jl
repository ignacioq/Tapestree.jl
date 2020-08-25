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

  return h, br, i
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
