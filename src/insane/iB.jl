#=

insane tree address structure `iB`

Ignacio Quintero Mächler

t(-_-t)

Created 03 07 2020
=#





"""
    iB

An abstract type for all branches in `iTree`.
"""
abstract type iB end




"""
    iBf

A Composite type representing node address for a **fixed** branch in `iTree`:

  `dr`: BitArray address where `true` = iTree.d1 and `false` = iTree.d2.
  `da`: mutable scalar denoting the number of grafted data augmented branches.
  `ti`: initial absolute time.
  `tf`: final absolute time.

    iBf()

Constructs an empty `iBf` object.
"""
struct iBf <: iB
  dr::BitArray{1}
  da::Base.RefValue{Int64}
  ti::Float64
  tf::Float64

  # constructors
  iBf() = new(BitArray{1}(), Ref(0), 0.0, 0.0)
  iBf(dr::BitArray{1}, da::Int64, ti::Float64, tf::Float64) = 
    new(dr, Ref(da), ti, tf)
end


# pretty-printing
Base.show(io::IO, id::iBf) = 
  print(io, "fixed ibranch (", ti(id), ", ", tf(id), "), ", dr(id), 
    " with ", da(id), " graft", isone(da(id)) ? "" : "s")




"""
    makeiBf!(tree::iTree, idv ::Array{iBf,1}, bit ::BitArray{1})

Make `iBf` vector for an `iTree`.
"""
function makeiBf!(tree::T, 
                  idv ::Array{iBf,1}, 
                  bit ::BitArray{1}) where {T <: iTree} 

  push!(idv, iBf(bit, 0, treeheight(tree), treeheight(tree) - pe(tree)))

  bit1 = copy(bit)
  bit2 = copy(bit)

  push!(bit1, true)
  makeiBf!(tree.d1, idv, bit1)

  push!(bit2, false)
  makeiBf!(tree.d2, idv, bit2)

  return nothing
end




"""
    makeiB(::Nothing, idv ::Array{iB,1}, bit ::BitArray{1})

Make `iBf` vector for an `iTree`.
"""
makeiBf!(::Nothing, idv::Array{iBf,1}, bit::BitArray{1}) = nothing




"""
    iBa

A Composite type representing node address for an **augmented** branch in `iTree`:

  `dr`: BitArray address where `true` = iTree.d1 and `false` = iTree.d2.
  `fB`: Link to `iBf` array specifying which fixed branch it attaches to.
  `ti`: initial absolute time.
  `tf`: final absolute time.

    iBa()

Constructs an empty `iBa` object.
"""
struct iBa <: iB
  dr::BitArray{1}
  fB::Int64
  ti::Float64
  tf::Float64

  # constructors
  iBa() = new(BitArray{1}(), 0, 0.0, 0.0)
  iBa(dr::BitArray{1}, fB::Int64, ti::Float64, tf::Float64) = 
    new(dr, fB, ti, tf)
end


# pretty-printing
Base.show(io::IO, id::iBa) = 
  print(io, "augmented ibranch (", ti(id), ", ", tf(id), "), ", dr(id), 
    " attached to ", fB(id))




"""
    dr(id::iB)

Return bit directory.
"""
dr(id::iB) = getproperty(id, :dr)




"""
    da(id::iB)

Return number of data augmented grafts.
"""
da(id::iBf) = getproperty(id, :da)[]




"""
    da(id::iB)

Return index of fixed branch it is attached to.
"""
fB(id::iBa) = getproperty(id, :fB)[]




"""
    ti(id::iB)

Return initial absolute time.
"""
ti(id::iB) = getproperty(id, :ti)




"""
    tf(id::iB)

Return final absolute time.
"""
tf(id::iB) = getproperty(id, :tf)








