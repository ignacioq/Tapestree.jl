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
function getbranch(wbr::BitArray{1}, rn::Int64, idv::Array{iDir,1})
  n = 0
  for (i, bit) in enumerate(wbr)
    if bit
      n += 1
      if n == rn
        return idv[i]
      end
    end
  end
end




"""
    branchescut!(wbr::BitArray{1}, h::Float64, idv::Array{iDir,1})

Fill the bit array with branches that are cut by  `h`.
"""
function branchescut!(wbr::BitArray{1}, h::Float64, idv::Array{iDir,1})
  n = 0
  for (i, b) in enumerate(idv)
    if ti(b) > h > tf(b)
      wbr[i] = true
      n += 1
    else
      wbr[i] = false
    end
  end

  return n
end
