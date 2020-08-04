#=

insane tree address structure `iDir`

Ignacio Quintero Mächler

t(-_-t)

Created 03 07 2020
=#



"""
    iDir

A Composite type representing node address for a branch in `iTree`. This is 
intended for the fixed tree.

  `dr`: BitArray address where `true` = iTree.d1 and `false` = iTree.d2.
  `ti`: initial absolute time.
  `tf`: final absolute time.

    iDir()

Constructs an empty `iDir` object.
"""
struct iDir
  dr::BitArray{1}
  ti::Float64
  tf::Float64

  # constructors
  iDir() = new(BitArray{1}(), 0.0, 0.0)
  iDir(dr::BitArray{1}, ti::Float64, tf::Float64) = new(dr, ti, tf)
end




"""
    dr(id::iDir)

Return bit directory.
"""
dr(id::iDir) = getproperty(id, :dr)




"""
    ti(id::iDir)

Return initial absolute time.
"""
ti(id::iDir) = getproperty(id, :ti)




"""
    tf(id::iDir)

Return final absolute time.
"""
tf(id::iDir) = getproperty(id, :tf)




"""
    makeiDir(tree::iTree, idv ::Array{iDir,1}, bit ::BitArray{1})

Make `iDir` vector for an `iTree`.
"""
function makeiDir(tree::iTree, idv ::Array{iDir,1}, bit ::BitArray{1})

  push!(idv, iDir(bit, treeheight(tree), treeheight(tree) - pe(tree)))

  bit1 = copy(bit)
  bit2 = copy(bit)
  
  push!(bit1, true)
  makeiDir(tree.d1, idv, bit1)

  push!(bit2, false)
  makeiDir(tree.d2, idv, bit2)

  return nothing
end

"""
    makeiDir(::Nothing, idv ::Array{iDir,1}, bit ::BitArray{1})

Make `iDir` vector for an `iTree`.
"""
makeiDir(::Nothing, idv::Array{iDir,1}, bit::BitArray{1}) = nothing







