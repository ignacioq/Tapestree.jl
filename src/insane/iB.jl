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
    makeiBf!(tree::sTfbd, idv::Array{iBfb,1}, bit::BitArray{1})

Make `iBfb` vector for an `iTree`.
"""
function makeiBf!(tree::sTfbd, idv::Array{iBfb,1}, bit::BitArray{1})

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
    makeiBf!(tree::T, 
             idv ::Array{iBfgp,1}, 
             bit ::BitArray{1}) where {T <: iTree} 

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
  `ie`: `true` if in a lineage that goes extinct.
  `iψ`: `true` if a fossil branch.
  `ni`: current direct alive descendants (≠ fixed ones in daugnter branches).
  `nt`: current alive descendants at time `t`.
  `λt`: final speciation rate for fixed at time `t`.
  `μt`: final extinction rate for fixed at time `t`.
  `ψt`: final fossilization rate for fixed at time `t`.

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
  iψ::Bool
  ni::Base.RefValue{Int64}
  nt::Base.RefValue{Int64}
  λt::Base.RefValue{Float64}
  μt::Base.RefValue{Float64}
  ψt::Base.RefValue{Float64}

  # constructors
  iBffs() = new(0., Ref(0), Ref(0), Ref(0), 0., 0., false, 1., false, false, 
                Ref(0), Ref(0), Ref(0.0), Ref(0.0), Ref(0.0))
  iBffs(t::Float64, pa::Int64, d1::Int64, d2::Int64, ti::Float64, tf::Float64, 
        it::Bool, ρi::Float64, ie::Bool, ni::Int64, nt::Int64, 
        λt::Float64, μt::Float64,) = 
        new(t, Ref(pa), Ref(d1), Ref(d2), ti, tf, it, ρi, ie, false, 
            Ref(ni), Ref(nt), Ref(λt), Ref(μt), Ref(0.0))
  iBffs(t::Float64, pa::Int64, d1::Int64, d2::Int64, ti::Float64, tf::Float64, 
        it::Bool, ρi::Float64, ie::Bool, iψ::Bool, ni::Int64, nt::Int64, 
        λt::Float64, μt::Float64, ψt::Float64) = 
        new(t, Ref(pa), Ref(d1), Ref(d2), ti, tf, it, ρi, ie, iψ, 
            Ref(ni), Ref(nt), Ref(λt), Ref(μt), Ref(ψt))
end

# pretty-printing
Base.show(io::IO, id::iBffs) = 
  print(io, 
    ie(id)         ? "extinct " : "", 
    ifos(id)       ? "fossil " : "", 
    it(id)         ? "terminal " : "", 
    iszero(pa(id)) ? "stem " : "", 
    isone(pa(id))  ? "crown " : "", 
    "ibranch (", ti(id), ", ", tf(id), 
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
    push!(idv, iBffs(el, 0, 0, 0, th, th - el, true, ρi, false, 1, 1, 0.0, 0.0))
    push!(n2v, 0)
    return ρi, 1
  end

  ρ1, n1 = makeiBf!(tree.d1, idv, n2v, tρ)
  ρ2, n2 = makeiBf!(tree.d2, idv, n2v, tρ)

  n  = n1 + n2 
  ρi = n / (n1/ρ1 + n2/ρ2)

  push!(idv, iBffs(el, 0, 1, 1, th, th - el, false, ρi, false, 0, 1, 0.0, 0.0))
  push!(n2v, n2)

  return ρi, n
end





"""
    makeiBf!(tree::sT_label, 
                  idv ::Array{iBffs,1}, 
                  n2v ::Array{Int64,1}, 
                  tρ  ::Dict{String, Float64},
                  sc  ::Array{Float64,1},
                  xr  ::Array{Float64,1},
                  X   ::Dict{String, Float64})

Make `iBf` vector for an `iTreeX` and estimate phylogenetic independent 
contrasts `sc` and ancestors `xr` given `tree` and data `X`.
"""
function makeiBf!(tree::sT_label, 
                  idv ::Array{iBffs,1}, 
                  n2v ::Array{Int64,1}, 
                  tρ  ::Dict{String, Float64},
                  sc  ::Array{Float64,1},
                  xr  ::Array{Float64,1},
                  X   ::Dict{String, Float64})

  th = treeheight(tree)
  el = e(tree)

  if istip(tree)
    lab = l(tree)
    ρi  = tρ[lab]
    push!(idv, iBffs(el, 0, 0, 0, th, th - el, true, ρi, false, 1, 1, 0.0, 0.0))
    push!(n2v, 0)
    xi  = X[lab]
    push!(xr, xi)
    return ρi, 1, xi, el
  end

  ρ1, n1, x1, e1 = makeiBf!(tree.d1, idv, n2v, tρ, sc, xr, X)
  ρ2, n2, x2, e2 = makeiBf!(tree.d2, idv, n2v, tρ, sc, xr, X)

  # tree order
  n  = n1 + n2 
  ρi = n / (n1/ρ1 + n2/ρ2)

  push!(idv, iBffs(el, 0, 1, 1, th, th - el, false, ρi, false, 0, 1, 0.0, 0.0))
  push!(n2v, n2)

  # pic
  scn = (x2 - x1)/(e1 + e2)
  xn = (x1/e1 + x2/e2) / (1.0/e1 + 1.0/e2)
  en = el + e1*e2/(e1 + e2)

  push!(sc, scn)
  push!(xr,  xn)

  return ρi, n, xn, en
end




"""
    makeiBf!(tree::sTf_label, 
             idv ::Array{iBffs,1}, 
             ti  ::Float64,
             n1v ::Array{Int64,1}, 
             n2v ::Array{Int64,1}, 
             ft1v::Array{Int64,1}, 
             ft2v::Array{Int64,1}, 
             sa2v::Array{Int64,1}, 
             tρ  ::Dict{String, Float64})

Make `iBf` vector for an `iTree` with fossils.
"""
function makeiBf!(tree::sTf_label, 
                  idv ::Array{iBffs,1}, 
                  ti  ::Float64,
                  n1v ::Array{Int64,1}, 
                  n2v ::Array{Int64,1}, 
                  ft1v::Array{Int64,1}, 
                  ft2v::Array{Int64,1}, 
                  sa2v::Array{Int64,1}, 
                  tρ  ::Dict{String, Float64})

  #th = treeheight(tree)
  el = e(tree)
  tf = ti - el

  if istip(tree)
    lab = l(tree)
    ρi  = tρ[lab]
    ie = iψ = isfossil(tree)       # count fossil tips as extinct
    push!(idv, iBffs(el, 0, 0, 0, ti, tf, true, ρi, ie, iψ, 
                     Int64(!iψ), 1, 0.0, 0.0, 0.0))
    
    push!(n1v, 0); push!(n2v, 0); push!(ft1v, 0); push!(ft2v, 0); push!(sa2v, 0)
    
    return iψ ? (ρi, 0, 1, 0) : (ρi, 1, 0, 0)
  end

  if isdefined(tree, :d1)
    ρ1,n1,ft1,sa1 = makeiBf!(tree.d1, idv, tf, n1v, n2v, ft1v, ft2v, sa2v, tρ)
  else
    ρ1,n1,ft1,sa1 = (1,0,0,0)
  end

  if isdefined(tree, :d2)
    ρ2,n2,ft2,sa2 = makeiBf!(tree.d2, idv, tf, n1v, n2v, ft1v, ft2v, sa2v, tρ)
  else
    ρ2,n2,ft2,sa2 = (1,0,0,0)
  end

  n  = n1 + n2                     # number of alive descendants
  ft = ft1 + ft2                   # number of fossil tips
  if n>0
    ρi = n / (n1/ρ1 + n2/ρ2)       # branch specific sampling fraction
  else
    ρi = ft / (ft1/ρ1 + ft2/ρ2)
  end
  sa = sa1 + sa2 + isfossil(tree)  # number of sampled ancestors

  push!(idv, iBffs(el, 0, 1, 1, ti, tf, false, ρi, false, 
                   isfossil(tree), 0, 1, 0.0, 0.0, 0.0))
  push!(n1v, n1)
  push!(n2v, n2)
  push!(ft1v, ft1)
  push!(ft2v, ft2)
  push!(sa2v, sa2)

  return ρi, n, ft, sa
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
    make_idf(tree::sT_label,
             tρ  ::Dict{String, Float64},
             X   ::Dict{String, Float64})

Make the edge dictionary and estimate ancestors and evolutionary rates given X
using phylogenetic independent contrasts.
"""
function make_idf(tree::sT_label,
                  tρ  ::Dict{String, Float64},
                  X   ::Dict{String, Float64})

  idf = iBffs[]
  n2v = Int64[]
  sc  = Float64[]
  xr  = Float64[]
  makeiBf!(tree, idf, n2v, tρ, sc, xr, X)

  σxi = sum(abs2, sc) / Float64(lastindex(sc)-1)

  reverse!(idf)
  reverse!(n2v)
  reverse!(xr)

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

  return idf, xr, σxi
end




"""
    make_idf(tree::sTf_label, tρ::Dict{String, Float64})

Make the edge dictionary.
"""
function make_idf(tree::sTf_label, tρ::Dict{String, Float64})

  idf  = iBffs[]
  n1v  = Int64[]; n2v = Int64[]    # vector of nb of alive descendants (d1 & d2)
  ft1v = Int64[]; ft2v = Int64[]  # vector of nb of fossil tips (d1 & d2)
  sa2v = Int64[]                  # vector of nb of sampled ancestors (d2)
  makeiBf!(tree, idf, treeheight(tree), n1v, n2v, ft1v, ft2v, sa2v, tρ)

  reverse!(idf)
  reverse!(n1v)
  reverse!(n2v)
  reverse!(ft1v)
  reverse!(ft2v)
  reverse!(sa2v)

  for i in Base.OneTo(lastindex(idf))
    bi = idf[i]
    n1 = n1v[i]; n2 = n2v[i]
    ft1 = ft1v[i]; ft2 = ft2v[i]
    sa2 = sa2v[i]
    #@show i; @show bi; @show n1; @show n2; @show sa2

    if n1+ft1>0 && n2+ft2>0 # bifurcation
      setd1!(bi, (n2+ft2)*2 + sa2 + i)
      setd2!(bi, i + 1)
      setpa!(idf[d1(bi)], i)
      setpa!(idf[d2(bi)], i)
    elseif n1+ft1>0    # fossil sampled ancestor of d1
      setd1!(bi, i + 1)
      setd2!(bi, 0)
      setpa!(idf[d1(bi)], i)
    elseif n2+ft2>0    # fossil sampled ancestor of d2
      setd1!(bi, 0)
      setd2!(bi, i + 1)
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
    pa(id::iBffs)

Return parent edge.
"""
pa(id::iBffs) = getproperty(id, :pa)[]




"""
    d1(id::iBffs)

Return daughter edge.
"""
d1(id::iBffs) = getproperty(id, :d1)[]




"""
    d2(id::iBffs)

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
    it(id::iBf)

Return if is extinct.
"""
ie(id::iBf) = getproperty(id, :ie)




"""
    ifos(id::iBffs)

Return if is a fossil.
"""
ifos(id::iBffs) = getproperty(id, :iψ)




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




"""
    nt(id::iBffs)

Return the current number of direct descendants alive at time `t`.
"""
nt(id::iBffs) = getproperty(id, :nt)[]




"""
    λt(id::iBffs)

Return final speciation rate for fixed at time `t
"""
λt(id::iBffs) = getproperty(id, :λt)[]




"""
    μt(id::iBffs)

Return final extinction rate for fixed at time `t`.
"""
μt(id::iBffs) = getproperty(id, :μt)[]




"""
    ψt(id::iBffs)

Return final fossil sampling rate for fixed at time `t`.
"""
ψt(id::iBffs) = getproperty(id, :ψt)[]




