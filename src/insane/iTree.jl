#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    iTree

An abstract type for all composite recursive types representing
a binary phylogenetic tree for `insane` use
"""
abstract type iTree end




"""
    sT

An abstract type for all composite recursive types representing
a simple binary phylogenetic tree for `insane` use
"""
abstract type sT <: iTree end




"""
    sT_label

A composite recursive type of supertype `sT`
representing a labelled binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  l:  label

    sT_label()

Constructs an empty `sT_label` object.

    sT_label(e::Float64)

Constructs an empty `sT_label` object with edge `e`.

    sT_label(d1::sT_label, d2::sT_label, e::Float64)

Constructs an `sT_label` object with two `sT_label` daughters and edge `e`.
"""
mutable struct sT_label <: sT
  d1::sT_label
  d2::sT_label
  e ::Float64
  l ::String

  sT_label() = new()
  sT_label(e::Float64, l::AbstractString) =
    (x = new(); x.e = e; x.l = l; x)
  sT_label(d1::sT_label, d2::sT_label, e::Float64, l::AbstractString) =
    new(d1, d2, e, l)
end

# pretty-printing
function Base.show(io::IO, t::sT_label)
  nt = ntips(t)
  print(io, "insane simple labelled tree with ", nt, " tip",
    (isone(nt) ? "" : "s"))
end



"""
    sT_label(tree::sT_label)

Copies a tree to `sT_label`.
"""
function sT_label(tree::sT_label)
  if def1(tree)
    sT_label(sT_label(tree.d1), sT_label(tree.d2), e(tree), label(tree))
  else
    sT_label(e(tree), label(tree))
  end
end



"""
    sT_label(tree::T) where {T <: iTree}

Demotes a tree to `sT_label`.
"""
function sT_label(tree::T) where {T <: iTree}
  _sT_label(tree::T, 0)[1]
end



"""
    _sT_label(tree::T, i::Int64) where {T <: iTree}

Demotes a tree to `sT_label`.
"""
function _sT_label(tree::T, i::Int64) where {T <: iTree}
  if def1(tree)
    t1, i = _sT_label(tree.d1, i)
    t2, i = _sT_label(tree.d2, i)
    lab = isdefined(tree, :l) ? label(tree) : ""
    tree  = sT_label(t1, t2, e(tree), lab)
  else
    i += 1
    lab = isdefined(tree, :l) ? label(tree) : string("t",i)
    tree = sT_label(e(tree), lab)
  end
  return tree, i
end



"""
    sT_label(tree::T, reftree::sT_label)

Demotes a `tree` to `sT_label` but retains the labels from `reftree`.
"""
function sT_label(tree::T, reftree::sT_label) where {T <: iTree}
  _sT_label(tree::T, reftree::sT_label, 0)[1]
end


"""
    _sT_label(tree::T, i::Int64) where {T <: iTree}

Demotes a `tree` to `sT_label` but retains the labels from `reftree`.
"""
function _sT_label(tree::T, reftree::sT_label, i::Int64) where {T <: iTree}

  if def1(tree)
    if isfix(tree.d1) && isfix(tree.d2)
      t1, i = _sT_label(tree.d1, reftree.d1, i)
      t2, i = _sT_label(tree.d2, reftree.d2, i)
      tree  = sT_label(t1, t2, e(tree), label(reftree))
    else
      t1, i = _sT_label(tree.d1, reftree, i)
      t2, i = _sT_label(tree.d2, reftree, i)
      tree  = sT_label(t1, t2, e(tree), "")
    end
  else
    if isfix(tree)
      tree = sT_label(e(tree), label(reftree))
    else
      i += 1
      tree = sT_label(e(tree), string("t",i))
    end
  end
  return tree, i
end



"""
    _sT_label(tree::rtree)

Convert an `rtree` object to an `sT_label` tree
"""
function sT_label(tree::rtree)

  ed   = tree.ed
  ed0  = ed[:,1]
  ed1  = ed[:,2]
  el   = tree.el
  ntip = tree.nnod + 1
  tl   = tree.tlab

  i1 = findfirst(isequal(ntip + 1), ed0)
  i2 = findnext(isequal(ntip + 1), ed0, i1+1)
  # if stem
  if isnothing(i2)
    stree = _sT_label(i1, ed0, ed1, el, tl)
  #if crown
  else
    stree = sT_label(_sT_label(i1, ed0, ed1, el, tl),
                     _sT_label(i2, ed0, ed1, el, tl), 
              0.0, "")
  end

  return stree
end


"""
    _sT_label(ix ::Int64, 
              ed0::Vector{Int64}, 
              ed1::Vector{Int64}, 
              el ::Vector{Float64}, 
              tl ::Vector{String})

Recursive conversion of an `rtree` object to an `sT_label` tree.
"""
function _sT_label(ix ::Int64, 
                   ed0::Vector{Int64}, 
                   ed1::Vector{Int64}, 
                   el ::Vector{Float64}, 
                   tl ::Vector{String})

  # inner node
  ei = el[ix]
  if iszero(ei) 
    ei = 1e-16
  end

  i1 = findfirst(isequal(ed1[ix]), ed0)
  if !isnothing(i1)
    i2 = findnext(isequal(ed1[ix]), ed0, i1 + 1)
    return sT_label(_sT_label(i1, ed0, ed1, el, tl),
                    _sT_label(i2, ed0, ed1, el, tl), 
            ei, "")
  # if tip
  else
    return sT_label(ei, tl[ed1[ix]])
  end
end



"""
    rtree(tree::sT_label)

Create a r tree object from a `sT_label`.
"""
function rtree(tree::sT_label)

  e0   = Int64[]
  e1   = Int64[]
  el   = Float64[]
  tlab = String[]

  n = ntips(tree)

  _rtree!(tree, e0, e1, el, tlab, n+1, n+1, 0)

  return rtree(hcat(e0, e1), el, tlab, n-1)
end



"""
    _rtree!(tree::sT_label, 
            e0  ::Vector{Int64}, 
            e1  ::Vector{Int64}, 
            el  ::Vector{Float64}, 
            tlab::Vector{String},
            pa  ::Int64,
            i   ::Int64,
            j   ::Int64)

Create a `rtree` object from a `sT_label` recursively.
"""
function _rtree!(tree::sT_label, 
                 e0  ::Vector{Int64}, 
                 e1  ::Vector{Int64}, 
                 el  ::Vector{Float64}, 
                 tlab::Vector{String},
                 pa  ::Int64,
                 i   ::Int64,
                 j   ::Int64)
  ei = e(tree)

  if def1(tree)
    if def2(tree)

      if ei > 0.0
        i += 1
        push!(e0, pa)
        pa = i
        push!(e1, i)
        push!(el, ei)
      end

      i, j = _rtree!(tree.d1, e0, e1, el, tlab, pa, i, j)
      i, j = _rtree!(tree.d2, e0, e1, el, tlab, pa, i, j)

    end
  else
    j += 1

    push!(e0, pa)
    push!(e1, j)
    push!(el, ei)
    push!(tlab, label(tree))
  end

  return i, j
end




"""
    sTf_label

A composite recursive type of supertype `sT` representing a labelled binary
phylogenetic tree with fossils for `insane` use, with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: if it is extinct
  iψ: if it is a fossil
  l:  label

    sTf_label(e::Float64, iμ::Bool, iψ::Bool, l::String)

Constructs an empty `sTf_label` object with edge `e` and label `l`.

    sTf_label(d1::sTf_label e::Float64, l::String)

Constructs an `sTf_label` object with one `sTf_label` daughter and edge `e`.

    sTf_label(d1::sTf_label, d2::sTf_label, e::Float64, l ::String)

Constructs an `sTf_label` object with two `sTf_label` daughters and edge `e`.
"""
mutable struct sTf_label <: sT
  d1::sTf_label
  d2::sTf_label
  e ::Float64
  iμ::Bool
  iψ::Bool
  l ::String

  sTf_label() = new()
  sTf_label(e::Float64, l::AbstractString) =
    (x=new(); x.e=e; x.iμ=false; x.iψ=false; x.l=l; x)
  sTf_label(e::Float64, iμ::Bool, iψ::Bool) =
    (x=new(); x.e=e; x.iμ=iμ; x.iψ=iψ; x.l=""; x)
  sTf_label(e::Float64, iμ::Bool, iψ::Bool, l::AbstractString) =
    (x=new(); x.e=e; x.iμ=iμ; x.iψ=iψ; x.l=l; x)
  sTf_label(d1::sTf_label, e::Float64, l::AbstractString) =
    (x=new(); x.d1=d1; x.e=e; x.iμ=false; x.iψ=true; x.l=l; x)
  sTf_label(d1::sTf_label, d2::sTf_label, e::Float64, l::AbstractString) =
    new(d1, d2, e, false, false, l)
end

# pretty-printing
Base.show(io::IO, t::sTf_label) =
  print(io, "insane simple fossil labelled tree with ", ntips(t),
    " tips (", ntipsextinct(t)," extinct) and ", nfossils(t)," fossils")



"""
    _sTf_label(tree::T, i::Int64) where {T <: iTree}

Copies a tree to `sTf_label`
"""
function sTf_label(tree::sTf_label)
  if def1(tree)
    if def2(tree)
      sTf_label(sTf_label(tree.d1), sTf_label(tree.d2), e(tree), label(tree))
    else
      sTf_label(sTf_label(tree.d1), e(tree), label(tree))
    end
  else
    sTf_label(e(tree), isextinct(tree), isfossil(tree), label(tree))
  end
end



"""
    sTf_label(tree::T) where {T <: iTree}

Demotes a tree to `sTf_label`.
"""
function sTf_label(tree::T) where {T <: iTree}
  _sTf_label(tree::T, 0, 0)[1]
end



"""
    _sTf_label(tree::T, i::Int64, j::Int64) where {T <: iTree}

Demotes a tree to `sTf_label`, initialized with label i.
"""
function _sTf_label(tree::T, n::Int64, nf::Int64) where {T <: iTree}

  if def1(tree)
    t1, n, nf = _sTf_label(tree.d1, n, nf)
    if def2(tree)
      t2, n, nf = _sTf_label(tree.d2, n, nf)
      tree  = sTf_label(t1, t2, e(tree), "")
    else
      nf += 1
      tree = sTf_label(t1, e(tree), string("f",nf))
    end
  else
    n += 1
    tree = sTf_label(e(tree), isextinct(tree), isfossil(tree), string("t",n))
  end

  return tree, n, nf
end




"""
    sTf_label(tree::T) where {T <: iTree}

Demotes a tree to `sTf_label`.
"""
function sTf_label(tree::T, reftree::sTf_label) where {T <: iTree}
  _sTf_label(tree::T, reftree, 0, 0)[1]
end

"""
    _sTf_label(tree::T, reftree::sTf_label, n::Int64, nf::Int64) where {T <: iTree}

Demotes a `tree` to `sTf_label` but retains the labels from `reftree`.
"""
function _sTf_label(tree::T, reftree::sTf_label, n::Int64, nf::Int64) where {T <: iTree}

  if def1(tree)
    if def2(tree)
      if isfix(tree.d1) && isfix(tree.d2)
        t1, n, nf = _sTf_label(tree.d1, reftree.d1, n, nf)
        t2, n, nf = _sTf_label(tree.d2, reftree.d2, n, nf)
        tree  = sTf_label(t1, t2, e(tree), label(reftree))
      else
        t1, n, nf = _sTf_label(tree.d1, reftree, n, nf)
        t2, n, nf = _sTf_label(tree.d2, reftree, n, nf)
        tree  = sTf_label(t1, t2, e(tree), "")
      end
    else
      if isfix(tree.d1)
        t1, n, nf = _sTf_label(tree.d1, reftree.d1, n, nf)
        tree = sTf_label(t1, e(tree), label(reftree))
      else
        nf += 1
        t1, n, nf = _sTf_label(tree.d1, reftree, n, nf)
        tree = sTf_label(t1, e(tree), string("f", nf))
      end
    end
  else
    if isfix(tree)
      tree = sTf_label(e(tree), label(reftree))
    else
      n += 1
      tree = sTf_label(e(tree), string("t", n))
    end
  end
  return tree, n, nf
end




"""
    sTpb

The simplest composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  fx: if fix

    sTpb(e::Float64)

Constructs an empty `sTpb` object with edge `e`.

    sTpb(d1::sTpb, d2::sTpb, e::Float64)

Constructs an `sTpb` object with two `sTpb` daughters and edge `e`.
"""
mutable struct sTpb <: sT
  d1::sTpb
  d2::sTpb
  e ::Float64
  fx::Bool

  sTpb() = new()
  sTpb(e::Float64) = (x = new(); x.e = e; x.fx = false; x)
  sTpb(e::Float64, fx::Bool) = (x = new(); x.e = e; x.fx = fx; x)
  sTpb(d1::sTpb, d2::sTpb, e::Float64, fx::Bool) = new(d1, d2, e, fx)
end

# pretty-printing
Base.show(io::IO, t::sTpb) =
  print(io, "insane simple pure-birth tree with ", ntips(t), " tips")



"""
    sTpb(tree::sT_label)

Transform a tree of type `sT_label` to `sTpb`.
"""
function sTpb(tree::sT_label)
  if def1(tree)
    sTpb(sTpb(tree.d1), sTpb(tree.d2), e(tree), false)
  else
    sTpb(e(tree), isfix(tree))
  end
end



"""
    sTpb(tree::sTpb)

Creates a new `sTpb` copy.
"""
function sTpb(tree::sTpb)
  if def1(tree)
    sTpb(sTpb(tree.d1), sTpb(tree.d2), e(tree), isfix(tree))
  else
    sTpb(e(tree), isfix(tree))
  end
end



"""
    sTbd

The simplest composite recursive type of supertype `sT`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: is an extinction node
  fx: if it is fix


    sTbd(e::Float64)

Constructs an empty `sTbd` object with edge `e`.

    sTbd(d1::sTbd, d2::sTbd, e::Float64)

Constructs an `sTbd` object with two `sTbd` daughters and edge `e`.
"""
mutable struct sTbd <: sT
  d1::sTbd
  d2::sTbd
  e ::Float64
  iμ::Bool
  fx::Bool

  sTbd() = new()
  sTbd(e::Float64, iμ::Bool, fx::Bool) =
    (x = new(); x.e = e; x.iμ = iμ; x.fx = fx; x)
  sTbd(d1::sTbd, d2::sTbd, e::Float64) =
    new(d1, d2, e, false, false)
  sTbd(d1::sTbd, d2::sTbd, e::Float64, iμ::Bool, fx::Bool) =
    new(d1, d2, e, iμ, fx)
end

# pretty-printing
Base.show(io::IO, t::sTbd) =
  print(io, "insane simple birth-death tree with ", ntips(t), " tips (",
    ntipsextinct(t)," extinct)")




"""
    sTbd(tree::sT_label)

Transform a tree of type `sT_label` to `sTbd`.
"""
function sTbd(tree::sT_label)
  if def1(tree)
    sTbd(sTbd(tree.d1), sTbd(tree.d2), e(tree), false, false)
  else
    sTbd(e(tree), false)
  end
end




"""
    sTbd(tree::sTbd)

Produce a new copy of `sTbd`.
"""
function sTbd(tree::sTbd)
  if def1(tree)
    sTbd(sTbd(tree.d1), sTbd(tree.d2), e(tree), isextinct(tree), isfix(tree))
  else
    sTbd(e(tree), isextinct(tree), isfix(tree))
  end
end




"""
    sTfbd

The simplest composite recursive type of supertype `sTf`
representing a binary phylogenetic tree for `insane` use,
with the following fields:

  d1: daughter tree 1
  d2: daughter tree 2
  e:  edge
  iμ: is an extinction node
  iψ: is a fossil node
  fx: if it is fix

    sTfbd(e::Float64)

Constructs an empty `sTfbd` object with edge `e`.

    sTfbd(d1::sTfbd, d2::sTfbd, e::Float64)

Constructs an `sTfbd` object with two `sTfbd` daughters and edge `e`.

    sTfbd(d1::sTfbd, e::Float64)

Constructs an `sTfbd` object with one sampled ancestor, one `sTfbd` daughter and
edge `e`.
"""
mutable struct sTfbd <: sT
  d1::sTfbd
  d2::sTfbd
  e ::Float64
  iμ::Bool
  iψ::Bool
  fx::Bool

  sTfbd() = new()
  sTfbd(e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
    (x = new(); x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x)
  sTfbd(d1::sTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
    (x = new(); x.d1 = d1; x.e = e; x.iμ = iμ; x.iψ = iψ; x.fx = fx; x)
  sTfbd(d1::sTfbd, d2::sTfbd, e::Float64, iμ::Bool, iψ::Bool, fx::Bool) =
    new(d1, d2, e, iμ, iψ, fx)
end

# pretty-printing
function Base.show(io::IO, t::sTfbd)
  nt = ntips(t)
  nf = nfossils(t)

  print(io, "insane simple fossil tree with ",
    nt , " tip",  (isone(nt) ? "" : "s" ),
    ", (", ntipsextinct(t)," extinct) and ",
    nf," fossil", (isone(nf) ? "" : "s" ))
end




"""
    sTfbd(tree::sTfbd)

Creates a copy of a `sTfbd` tree.
"""
function sTfbd(tree::sTfbd)
  if def1(tree)
    if def2(tree)
      sTfbd(sTfbd(tree.d1), sTfbd(tree.d2), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree))
    else
      sTfbd(sTfbd(tree.d1), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree))
    end
  else
    sTfbd(e(tree), isextinct(tree), isfossil(tree), isfix(tree))
  end
end




"""
    sTfbd_wofe(tree::sTfbd)

Creates a copy of a `sTfbd` tree without fossils extinct tips.
"""
function sTfbd_wofe(tree::sTfbd)
  if def1(tree) && def2(tree)
      sTfbd(sTfbd_wofe(tree.d1), sTfbd_wofe(tree.d2), e(tree),
        isextinct(tree), isfossil(tree), isfix(tree))
  else
    sTfbd(e(tree), isextinct(tree), isfossil(tree), isfix(tree))
  end
end




