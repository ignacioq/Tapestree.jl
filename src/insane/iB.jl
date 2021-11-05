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
   `ρ`: specific sampling fraction.
  `ti`: initial absolute time.
  `tf`: final absolute time.
  `it`: `true` if a terminal branch.

    iBfb()

Constructs an empty `iBfb` object.
"""
struct iBfb <: iBf
  dr::BitArray{1}
  ti::Float64
  tf::Float64
  ρi::Float64
  it::Bool

  # constructors
  iBfb() = new(BitArray{1}(), 0.0, 0.0, 1.0, false)
  iBfb(dr::BitArray{1}, ti::Float64, tf::Float64, ρi::Float64, it::Bool) = 
    new(dr, ρi, ti, tf, it)
end


# pretty-printing
Base.show(io::IO, id::iBfb) = 
  print(io, "fixed", it(id) ? " terminal" : ""," ibranch (", 
    ti(id), ", ", tf(id), "), ", dr(id))




"""
    tree::sT_label, 
    idv ::Array{iBfb,1}, 
    bitv::BitArray{1},
    tρ  ::Dict{String, Float64}

Make `iBfb` vector for an `iTree` taking into account 
species-specific sampling fraction `ρ`.
"""
function makeiBf!(tree::sT_label, 
                  idv ::Array{iBfb,1}, 
                  bitv::BitArray{1},
                  tρ  ::Dict{String, Float64})

  if istip(tree)
    lab = l(tree)
    ρi  = tρ[lab]
    push!(idv, 
      iBfb(bitv, treeheight(tree), treeheight(tree) - e(tree), ρi, true))
    return ρi, 1, bitv
  end

  bitv1 = copy(bitv)
  bitv2 = copy(bitv)

  if isdefined(tree, :d1)
    push!(bitv1, true)
    ρ1, n1, bitv1 = makeiBf!(tree.d1, idv, bitv1, tρ)
  end

  if isdefined(tree, :d2)
    push!(bitv2, false)
    ρ2, n2, bitv2 = makeiBf!(tree.d2, idv, bitv2, tρ)
  end

  bitv = copy(bitv1)
  pop!(bitv)

  n  = n1 + n2 
  ρi = n / (n1/ρ1 + n2/ρ2)
  push!(idv, 
    iBfb(bitv, treeheight(tree), treeheight(tree) - e(tree), ρi, false))

  return ρi, n, bitv
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

  `t` : edge length.
  `pa`: parent node
  `d1`: daughter 1 node
  `d2`: daughter 2 node
  `ti`: initial absolute time.
  `tf`: final absolute time.
  `ie`: `true` if an extinct branch.
  `ni`: Current direct alive descendants.

    iBffs()

Constructs an empty `iBf` object.
"""
struct iBffs <: iBf
  t ::Float64
  pa::Base.RefValue{Int64}
  d1::Base.RefValue{Int64}
  d2::Base.RefValue{Int64}
  ti::Float64
  tf::Float64
  ρi::Float64
  ie::Bool
  ni::Base.RefValue{Int64}

  # constructors
  iBffs() = new(0., Ref(0), Ref(0), Ref(0), 
                0., 0., 1., false, Ref(0))
  iBffs(t::Float64, pa::Int64, d1::Int64, d2::Int64, 
    ti::Float64, tf::Float64, ρi::Float64, it::Bool, ie::Bool, ni::Int64) = 
    new(t, Ref(pa), Ref(d1), Ref(d2), ti, tf, ρi, ie, Ref(ni))
end


# pretty-printing
Base.show(io::IO, id::iBffs) = 
  print(io, "fixed", 
    iszero(d1(id))   ? " terminal" : "", 
    # iszero(sc(id)) ? " stem" : "", 
    # isone(sc(id))  ? " crown" : "", 
    " ibranch (", ti(id), ", ", tf(id), "), p:", 
    pa(id), ", d1:", d1(id), ", d2:", d2(id))




"""
    makeiBf!(tree::sT_label, 
             idv ::Array{iBffs,1}, 
             n2v ::Array{Int64,1}, 
             tρ  ::Dict{String, Float64})

Make `iBf` vector for an `iTree`.
"""
function makeiBf!(tree::sT_label, 
                  idv ::Array{iBffs,1}, 
                  n2v ::Array{Int64,1}, 
                  tρ  ::Dict{String, Float64})

  if istip(tree)
    lab = l(tree)
    ρi  = tρ[lab]
    push!(idv, 
      iBffs(e(tree), 0, 0, 0, treeheight(tree), 
            treeheight(tree) - e(tree), ρi, true, false, 1))
    push!(n2v, 0)
    return ρi, 1
  end

  if isdefined(tree, :d1)
    ρ1, n1 = makeiBf!(tree.d1, idv, n2v, tρ)
  end

  if isdefined(tree, :d2)
    ρ2, n2 = makeiBf!(tree.d2, idv, n2v, tρ)
  end

  n  = n1 + n2 
  ρi = n / (n1/ρ1 + n2/ρ2)

  push!(idv, 
    iBffs(e(tree), 0, 1, 1, treeheight(tree), 
      treeheight(tree) - e(tree), ρi, false, false, 0))
  push!(n2v, n2)

  return ρi, n
end




"""
    make_idf(tree::sT_label, tρ::Dict{String, Float64})

Make the edge dictionary.
"""
function make_idf(tree::sT_label, tρ::Dict{String, Float64})

  idf = iBffs[]
  n2v = Int64[]
  makeiBf!(tree, idf, n2v, tρ)

  reverse!(idf)
  reverse!(n2v)

  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    n2 = n2v[i]

    if n2 > 0
      setd1!(bi, n2*2 + i)
      setd2!(bi, i + 1)
      setpa!(idf[d1(bi)], i)
      setpa!(idf[d2(bi)], i)
    end
  end

  return idf
end




"""
    prob_ρ(idv::Array{iBffs,1})

Estimate initial sampling fraction probability without augmented data.
"""
function prob_ρ(idv::Array{iBffs,1})
  ll = 0.0
  for bi in idv
    nbi = ni(bi)
    if iszero(d1(bi))
      ll += log(Float64(nbi) * (1.0 - ρi(bi))^(nbi - 1))
    else
      ll += log((1.0 - ρi(bi))^(nbi))
    end
  end
  return ll
end





"""
    makeiBf!(tree::sTfbd, idv ::Array{iBf,1}, bit ::BitArray{1})

Make `iBf` vector for an `iTree`.
"""
function makeiBf!(tree::sTfbd, 
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
  end
  if isdefined(tree, :d2)
    push!(bit2, false)
    makeiBf!(tree.d2, idv, bit2)
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
    pa(id::iB)

Return parent edge.
"""
pa(id::iBffs) = getproperty(id, :pa)[]



"""
    d1(id::iB)

Return daughter edge.
"""
d1(id::iBffs) = getproperty(id, :d1)[]



"""
    d2(id::iB)

Return daughter edge.
"""
d2(id::iBffs) = getproperty(id, :d2)[]




"""
    e(id::iB)

Return initial absolute time.
"""
e(id::iB) = getproperty(id, :t)




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

Return `0` if stem branch, `1` if either of the crown branches and `23` if 
another plebeian branch.
"""
sc(id::iBffs) = getproperty(id, :sc)




"""
    ρi(id::iBffs)

Return the branch-specific sampling fraction. 
"""
ρi(id::iBffs) = getproperty(id, :ρi)




"""
    ni(id::iBffs)

Return the current number of direct descendants alive at the present.
"""
ni(id::iBffs) = getproperty(id, :ni)[]


