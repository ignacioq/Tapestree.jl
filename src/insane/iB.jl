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
    iB

An abstract type for all fixed branches in `iTree`.
"""
abstract type iBf <: iB end




"""
    iBfb

A Composite type representing node address for a **fixed** branch in `iTree`:

  `dr`: BitArray address where `true` = iTree.d1 and `false` = iTree.d2.
  `ti`: initial absolute time.
  `tf`: final absolute time.
  `it`: `true` if a terminal branch.

    iBfb()

Constructs an empty `iBfb` object.
"""
struct iBfb <: iBf
  dr::BitArray{1}
  da::Base.RefValue{Int64}
  ti::Float64
  tf::Float64
  it::Bool

  # constructors
  iBfb() = new(BitArray{1}(), Ref(0), 0.0, 0.0, false)
  iBfb(dr::BitArray{1}, da::Int64, ti::Float64, tf::Float64, it::Bool) = 
    new(dr, Ref(da), ti, tf, it)
end


# pretty-printing
Base.show(io::IO, id::iBfb) = 
  print(io, "fixed", it(id) ? " terminal" : ""," ibranch (", ti(id), ", ", tf(id), "), ", dr(id))




"""
    makeiBf!(tree::iTree, idv ::Array{iBfb,1}, bit ::BitArray{1})

Make `iBfb` vector for an `iTree`.
"""
function makeiBf!(tree::T, 
                  idv ::Array{iBfb,1}, 
                  bit ::BitArray{1}) where {T <: iTree} 

  itb = istip(tree)

  push!(idv, iBfb(bit, 0, treeheight(tree), treeheight(tree) - e(tree), itb))

  bit1 = copy(bit)
  bit2 = copy(bit)


  if isdefined(tree, :d1)
    push!(bit1, true)
    makeiBf!(tree.d1, idv, bit1)
    push!(bit2, false)
    makeiBf!(tree.d2, idv, bit2)
  end
  return nothing
end




"""
    makeiBf!(tree::sTfbd, idv ::Array{iBfb,1}, bit ::BitArray{1})

Make `iBfb` vector for an `iTree`.
"""
function makeiBf!(tree::sTfbd, idv ::Array{iBfb,1}, bit ::BitArray{1})

  itb = istip(tree)

  push!(idv, iBfb(bit, 0, treeheight(tree), treeheight(tree) - e(tree), itb))

  bit1 = copy(bit)
  bit2 = copy(bit)


  if isdefined(tree, :d1)
    push!(bit1, true)
    makeiBf!(tree.d1, idv, bit1)
  end
  if isdefined(tree, :d2)
    push!(bit2, false)
    makeiBf!(tree.d2, idv, bit2)
  end
  return nothing
end




"""
    iBfgp

A Composite type representing node address for a **fixed** branch in `iTree`:

  `dr`: BitArray address where `true` = iTree.d1 and `false` = iTree.d2.
  `da`: mutable scalar denoting the number of grafted data augmented branches.
  `ti`: initial absolute time.
  `tf`: final absolute time.
  `it`: `true` if a terminal branch.
  `ie`: `true` if a terminal branch and extinct.

    iBfgp()

Constructs an empty `iBfgp` object.
"""
struct iBfgp <: iBf
  dr::BitArray{1}
  da::Base.RefValue{Int64}
  ti::Float64
  tf::Float64
  it::Bool
  ie::Bool

  # constructors
  iBfgp() = new(BitArray{1}(), Ref(0), 0.0, 0.0, false, false)
  iBfgp(dr::BitArray{1}, da::Int64, ti::Float64, tf::Float64, it::Bool, ie::Bool) = 
    new(dr, Ref(da), ti, tf, it, ie)
end


# pretty-printing
Base.show(io::IO, id::iBfgp) = 
  print(io, "fixed", it(id) ? " terminal" : ""," ibranch (", ti(id), ", ", tf(id), "), ", dr(id), 
    " with ", da(id), " graft", isone(da(id)) ? "" : "s")




"""
    makeiBf!(tree::iTree, idv ::Array{iBfgp,1}, bit ::BitArray{1})

Make `iBfgp` vector for an `iTree`.
"""
function makeiBf!(tree::T, 
                  idv ::Array{iBfgp,1}, 
                  bit ::BitArray{1}) where {T <: iTree} 

  itb = istip(tree)
  ieb = isextinct(tree)

  push!(idv, 
    iBfgp(bit, 0, treeheight(tree), treeheight(tree) - e(tree), itb, ieb))

  bit1 = copy(bit)
  bit2 = copy(bit)

  if isdefined(tree, :d1)
    push!(bit1, true)
    makeiBf!(tree.d1, idv, bit1)
    push!(bit2, false)
    makeiBf!(tree.d2, idv, bit2)
  end

  return nothing
end




"""
    makeiBf!(tree::sTfbd, idv ::Array{iBfgp,1}, bit ::BitArray{1})

Make `iBfgp` vector for an `iTree`.
"""
function makeiBf!(tree::sTfbd, 
                  idv ::Array{iBfgp,1}, 
                  bit ::BitArray{1}) where {T <: iTree} 

  itb = istip(tree)
  ieb = isextinct(tree)

  push!(idv, 
    iBfgp(bit, 0, treeheight(tree), treeheight(tree) - e(tree), itb, ieb))

  bit1 = copy(bit)
  bit2 = copy(bit)

  if isdefined(tree, :d1)
    push!(bit1, true)
    makeiBf!(tree.d1, idv, bit1)
  end
  if isdefined(tree, :d2)
    push!(bit2, false)
    makeiBf!(tree.d2, idv, bit2)
  end

  return nothing
end




"""
    iBffs

A Composite type representing node address for a **fixed** branch in `iTree`:

  `dr`: BitArray address where `true` = iTree.d1 and `false` = iTree.d2.
  `ti`: initial absolute time.
  `tf`: final absolute time.
  `it`: `true` if a terminal branch.
  `ie`: `true` if an extinct branch.
  `sc`: is `0` if stem branch, `1` if either of the crown branches and `23` if 
        another plebeian branch.

    iBffs()

Constructs an empty `iBf` object.
"""
struct iBffs <: iBf
  dr::BitArray{1}
  ti::Float64
  tf::Float64
  it::Bool
  ie::Bool
  sc::Int64

  # constructors
  iBffs() = new(BitArray{1}(), 0.0, 0.0, false, false, 23)
  iBffs(dr::BitArray{1}, ti::Float64, tf::Float64, it::Bool, ie::Bool, 
    sc::Int64) = 
    new(dr, ti, tf, it, ie, sc)
end


# pretty-printing
Base.show(io::IO, id::iBffs) = 
  print(io, "fixed", 
    it(id)     ? " terminal" : "", 
    iszero(sc(id)) ? " stem" : "", 
    isone(sc(id))  ? " crown" : "", 
    " ibranch (", ti(id), ", ", tf(id), "), ", dr(id))




"""
    makeiBf!(tree::iTree, idv ::Array{iBf,1}, bit ::BitArray{1})

Make `iBf` vector for an `iTree`.
"""
function makeiBf!(tree::T, 
                  idv ::Array{iBffs,1}, 
                  bit ::BitArray{1}) where {T <: iTree} 

  itb = istip(tree)
  ieb = isextinct(tree)

  lb = lastindex(bit)

  sc = 23
  if iszero(lb)
    sc = 0
  elseif isone(lb)
    sc = 1
  end

  push!(idv, 
    iBffs(bit, treeheight(tree), treeheight(tree) - e(tree), itb, ieb, sc))

  bit1 = copy(bit)
  bit2 = copy(bit)

  if isdefined(tree, :d1)
    push!(bit1, true)
    makeiBf!(tree.d1, idv, bit1)
    push!(bit2, false)
    makeiBf!(tree.d2, idv, bit2)
  end

  return nothing
end




"""
    iBfffs

A Composite type representing node address for a **fixed** branch in `iTree`:

  `dr`: BitArray address where `true` = iTree.d1 and `false` = iTree.d2.
  `ti`: initial absolute time.
  `tf`: final absolute time.
  `it`: `true` if a terminal branch.
  `ie`: `true` if an extinct branch.
  `iψ`: `true` if a fossil branch.
  `sc`: is `0` if stem branch, `1` if either of the crown branches and `23` if 
        another plebeian branch.

    iBfffs()

Constructs an empty `iBf` object.
"""
struct iBfffs <: iBf
  dr::BitArray{1}
  ti::Float64
  tf::Float64
  it::Bool
  ie::Bool
  iψ::Bool
  sc::Int64

  # constructors
  iBfffs() = new(BitArray{1}(), 0.0, 0.0, false, false, false, 23)
  iBfffs(dr::BitArray{1}, ti::Float64, tf::Float64, it::Bool, ie::Bool, iψ::Bool, 
    sc::Int64) = 
    new(dr, ti, tf, it, ie, iψ, sc)
end


# pretty-printing
Base.show(io::IO, id::iBfffs) = 
  print(io, "fixed", 
    it(id)     ? " terminal" : "",
    ifos(id)     ? " fossil" : "", 
    iszero(sc(id)) ? " stem" : "", 
    isone(sc(id))  ? " crown" : "", 
    " ibranch (", ti(id), ", ", tf(id), "), ", dr(id))




"""
    makeiBf!(tree::sTfbd, idv ::Array{iBf,1}, bit ::BitArray{1})

Make `iBf` vector for an `iTree`.
"""
function makeiBf!(tree::sTfbd, 
                  idv::Array{iBfffs,1}, 
                  bit::BitArray{1}) where {T <: iTree}
  
  return _makeiBf!(tree, idv, bit, treeheight(tree), 0)
end




"""
    _makeiBf!(tree::sTfbd, 
              idv ::Array{iBf,1}, 
              bit ::BitArray{1}, 
              ti::Float64, 
              sc::Int64)

Make `iBf` vector for an `iTree`, initialized at time `ti` and status `sc`.
"""
function _makeiBf!(tree::sTfbd, 
                   idv ::Array{iBfffs,1}, 
                   bit ::BitArray{1},
                   ti::Float64,
                   sc::Int64) where {T <: iTree}

  itb = istip(tree)
  ieb = isextinct(tree)
  iψb = isfossil(tree)

  lb = lastindex(bit)


  tf = ti-e(tree)

  push!(idv, iBfffs(bit, ti, tf, itb, ieb, iψb, sc))

  bit1 = copy(bit)
  bit2 = copy(bit)

  survd1 = isdefined(tree, :d1) && survives(tree.d1)
  survd2 = isdefined(tree, :d2) && survives(tree.d2)
  
  if isdefined(tree, :d1)
    push!(bit1, true)
    _makeiBf!(tree.d1, idv, bit1, tf, sc + (survd1&&survd2))
  end
  if isdefined(tree, :d2)
    push!(bit2, false)
    _makeiBf!(tree.d2, idv, bit2, tf, sc + (survd1&&survd2))
  end

  return nothing
end




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
    findsubtree(tree::sTfbd, dri::BitArray{1})

Return the subtree in `tree` localized by the directory `dri`.
"""
function findsubtree(tree::sTfbd, dri::BitArray{1})
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  fixd1 = defd1 && isfix(tree.d1)
  fixd2 = defd2 && isfix(tree.d2)
  if (fixd1 && fixd2) || (fixd1 && !defd2) || (!defd1 && fixd2)
    if isone(lastindex(dri))
      return dri[1] ? tree.d1 : tree.d2
    else
      return findsubtree(dri[1] ? tree.d1 : tree.d2, dri[2:end])
    end
  else
    return findsubtree(fixd1 ? tree.d1 : tree.d2, dri::BitArray{1})
  end
end




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




"""
    it(id::iBf)

Return if is terminal.
"""
it(id::iBf) = getproperty(id, :it)




"""
    sc(id::iBffs)
    sc(id::iBfffs)

Return `0` if stem branch, `1` if either of the crown branches and `23` if 
another plebeian branch.
"""
sc(id::iBffs) = getproperty(id, :sc)
sc(id::iBfffs) = getproperty(id, :sc)




"""
    ifos(id::iBfffs)

Return if is a fossil.
"""
ifos(id::iBfffs) = getproperty(id, :iψ)





