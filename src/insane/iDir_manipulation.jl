#=

iDir manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 07 2020
=#



"""
    getbranch(wbr::BitArray{1}, rn::Int64, idv::Array{iDir,1})

Sample one branch among those in `wbr` that are `true` according to 
randomly picked `rn`.
"""
function getbranch(wbr::BitArray{1}, 
                   rn ::Int64, 
                   idv::Array{iDir,1})::iDir
  i::Int64 = 0
  n::Int64 = 0
  for bit in wbr
    i += 1
    if bit
      n += 1
      if n == rn
        return idv[i]::iDir
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
