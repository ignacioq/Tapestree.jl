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
  `it`: `true` if a terminal branch.
  `ρi`: branch specific sampling fraction.
  `ie`: `true` if an extinct branch.
  `ni`: current direct alive descendants.
  `nt`: current alive descendants at time `t`.
  `λt`: final speciation rate for fixed at time `t`.
  `μt`: final extinction rate for fixed at time `t`.

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
  it::Bool
  ρi::Float64
  ie::Bool
  ni::Base.RefValue{Int64}
  nt::Base.RefValue{Int64}
  λt::Base.RefValue{Float64}
  μt::Base.RefValue{Float64}

  # constructors
  iBffs() = new(0., Ref(0), Ref(0), Ref(0), 
                0., 0., false, 1., false, Ref(0), Ref(0), Ref(0.0), Ref(0.0))
  iBffs(t::Float64, pa::Int64, d1::Int64, d2::Int64, ti::Float64, tf::Float64, 
    it::Bool, ρi::Float64, ie::Bool, ni::Int64, nt::Int64, 
    λt::Float64, μt::Float64) = 
    new(t, Ref(pa), Ref(d1), Ref(d2), ti, tf, it, ρi, ie, 
      Ref(ni), Ref(nt), Ref(λt), Ref(μt))
end


# pretty-printing
Base.show(io::IO, id::iBffs) = 
  print(io, "fixed", 
    it(id)         ? " terminal" : "", 
    iszero(pa(id)) ? " stem" : "", 
    isone(pa(id))  ? " crown" : "", 
    " ibranch (", ti(id), ", ", tf(id), 
    "), p:", pa(id), ", d1:", d1(id), ", d2:", d2(id))




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

  th = treeheight(tree)
  el = e(tree)

  if istip(tree)
    lab = l(tree)
    ρi  = tρ[lab]
    push!(idv, 
      iBffs(el, 0, 0, 0, th, th - el, true, ρi, false, 1, 1, 0.0, 0.0))
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
    iBffs(el, 0, 1, 1, th, th - el, false, ρi, false, 0, 1, 0.0, 0.0))
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
    if it(bi)
      ll += log(Float64(nbi) * ρi(bi) * (1.0 - ρi(bi))^(nbi - 1))
    else
      ll += log((1.0 - ρi(bi))^(nbi))
    end
  end
  return ll
end




"""
    iBfffs

A Composite type representing node address for a **fixed** branch in `iTree` with fossils:

  `t` : edge length.
  `pa`: parent node
  `d1`: daughter 1 node
  `d2`: daughter 2 node
  `ti`: initial absolute time.
  `tf`: final absolute time.
  `it`: `true` if a terminal branch.
  `ρi`: branch specific sampling fraction.
  `ie`: `true` if an extinct branch.
  `iψ`: `true` if a fossil branch.
  `ni`: current direct alive descendants.
  `nt`: current alive descendants at time `t`.
  `λt`: final speciation rate for fixed at time `t`.
  `μt`: final extinction rate for fixed at time `t`.
  `ψt`: final fossilization rate for fixed at time `t`.

    iBfffs()

Constructs an empty `iBf` object.
"""
struct iBfffs <: iBf
  t ::Float64
  pa::Base.RefValue{Int64}
  d1::Base.RefValue{Int64}
  d2::Base.RefValue{Int64}
  ti::Float64
  tf::Float64
  it::Bool
  ρi::Float64
  ie::Bool
  iψ::Bool
  ni::Base.RefValue{Int64}
  nt::Base.RefValue{Int64}
  λt::Base.RefValue{Float64}
  μt::Base.RefValue{Float64}
  ψt::Base.RefValue{Float64}

  # constructors
  iBfffs() = new(0., Ref(0), Ref(0), Ref(0), 0., 0., false, 1., false, false, 
                 Ref(0), Ref(0), Ref(0.0), Ref(0.0), Ref(0.0))
  iBfffs(t::Float64, pa::Int64, d1::Int64, d2::Int64, ti::Float64, tf::Float64, 
         it::Bool, ρi::Float64, ie::Bool, iψ::Bool, ni::Int64, nt::Int64, 
         λt::Float64, μt::Float64, ψt::Float64) = 
         new(t, Ref(pa), Ref(d1), Ref(d2), ti, tf, it, ρi, ie, iψ, 
             Ref(ni), Ref(nt), Ref(λt), Ref(μt), Ref(ψt))
end


# pretty-printing
Base.show(io::IO, id::iBfffs) = 
  print(io, "fixed", 
    it(id)         ? " terminal" : "", 
    ifos(id)       ? " fossil" : "", 
    iszero(pa(id)) ? " stem" : "", 
    isone(pa(id))  ? " crown" : "", 
    " ibranch (", ti(id), ", ", tf(id), 
    "), p:", pa(id), ", d1:", d1(id), ", d2:", d2(id))




#="""
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
end=#




"""
    makeiBf!(tree::sTf_label, 
             idv ::Array{iBfffs,1}, 
             ti::Float64,
             n1v ::Array{Int64,1}, 
             n2v ::Array{Int64,1}, 
             sa2v ::Array{Int64,1}, 
             tρ  ::Dict{String, Float64})

Make `iBf` vector for an `iTree` with fossils.
"""
function makeiBf!(tree::sTf_label, 
                  idv ::Array{iBfffs,1}, 
                  ti::Float64,
                  n1v ::Array{Int64,1}, 
                  n2v ::Array{Int64,1}, 
                  sa2v ::Array{Int64,1}, 
                  tρ  ::Dict{String, Float64})

  #th = treeheight(tree)
  el = e(tree)
  tf = ti - el

  if istip(tree)
    lab = l(tree)
    ρi  = tρ[lab]
    push!(idv, iBfffs(el, 0, 0, 0, ti, tf, true, ρi, isextinct(tree), 
                      isfossil(tree), 1, 1, 0.0, 0.0, 0.0))
    push!(n1v, 0)
    push!(n2v, 0)
    push!(sa2v, 0)
    return ρi, 1, 0
  end

  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  ρ1,n1,sa1 = defd1 ? makeiBf!(tree.d1, idv, tf, n1v, n2v, sa2v, tρ) : (1,0,0)
  ρ2,n2,sa2 = defd2 ? makeiBf!(tree.d2, idv, tf, n1v, n2v, sa2v, tρ) : (1,0,0)


  n  = n1 + n2
  ρi = n / (n1/ρ1 + n2/ρ2)
  sa = sa1 + sa2 + isfossil(tree)

  push!(idv, iBfffs(el, 0, 1, 1, ti, tf, false, ρi, false, 
                    isfossil(tree), 0, 1, 0.0, 0.0, 0.0))
  push!(n1v, n1)
  push!(n2v, n2)
  push!(sa2v, sa2)

  return ρi, n, sa
end




"""
    make_idf(tree::sTf_label, tρ::Dict{String, Float64})

Make the edge dictionary.
"""
function make_idf(tree::sTf_label, tρ::Dict{String, Float64})

  idf = iBfffs[]
  n1v = Int64[]
  n2v = Int64[]
  sa2v = Int64[]
  makeiBf!(tree, idf, treeheight(tree), n1v, n2v, sa2v, tρ)

  reverse!(idf)
  reverse!(n1v)
  reverse!(n2v)
  reverse!(sa2v)

  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    n1 = n1v[i]
    n2 = n2v[i]
    sa2 = sa2v[i]
    #@show i; @show bi; @show n1; @show n2; @show sa2

    if n1>0 && n2>0 # bifurcation
      setd1!(bi, n2*2 + sa2 + i)
      setd2!(bi, i + 1)
      setpa!(idf[d1(bi)], i)
      setpa!(idf[d2(bi)], i)
    elseif n1>0    # fossil sampled ancestor of d1
      setd1!(bi, i + 1)
      setd2!(bi, 0)
      setpa!(idf[d1(bi)], i)
    elseif n2>0    # fossil sampled ancestor of d2
      setd1!(bi, 0)
      setd2!(bi, i + 1)
      setpa!(idf[d2(bi)], i)
    end
  end

  return idf
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
    pa(id::iBffs)
    pa(id::iBfffs)

Return parent edge.
"""
pa(id::iBffs) = getproperty(id, :pa)[]
pa(id::iBfffs) = getproperty(id, :pa)[]




"""
    d1(id::iBffs)
    d1(id::iBfffs)

Return daughter edge.
"""
d1(id::iBffs) = getproperty(id, :d1)[]
d1(id::iBfffs) = getproperty(id, :d1)[]




"""
    d2(id::iBffs)
    d2(id::iBfffs)

Return daughter edge.
"""
d2(id::iBffs) = getproperty(id, :d2)[]
d2(id::iBfffs) = getproperty(id, :d2)[]




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
    it(id::iBf)

Return if is extinct.
"""
ie(id::iBf) = getproperty(id, :ie)




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




"""
    ρi(id::iBffs)
    ρi(id::iBfffs)

Return the branch-specific sampling fraction. 
"""
ρi(id::iBffs) = getproperty(id, :ρi)
ρi(id::iBfffs) = getproperty(id, :ρi)




"""
    ni(id::iBffs)
    ni(id::iBfffs)

Return the current number of direct descendants alive at the present.
"""
ni(id::iBffs) = getproperty(id, :ni)[]
ni(id::iBfffs) = getproperty(id, :ni)[]




"""
    nt(id::iBffs)
    nt(id::iBfffs)

Return the current number of direct descendants alive at time `t`.
"""
nt(id::iBffs) = getproperty(id, :nt)[]
nt(id::iBfffs) = getproperty(id, :nt)[]




"""
    λt(id::iBffs)
    λt(id::iBfffs)

Return final speciation rate for fixed at time `t
"""
λt(id::iBffs) = getproperty(id, :λt)[]
λt(id::iBfffs) = getproperty(id, :λt)[]




"""
    μt(id::iBffs)
    μt(id::iBfffs)

Return final extinction rate for fixed at time `t`.
"""
μt(id::iBffs) = getproperty(id, :μt)[]
μt(id::iBfffs) = getproperty(id, :μt)[]




"""
    ψt(id::iBfffs)

Return final fossil sampling rate for fixed at time `t`.
"""
ψt(id::iBfffs) = getproperty(id, :ψt)[]




