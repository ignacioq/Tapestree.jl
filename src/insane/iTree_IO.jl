#=

insane tree read and write

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#



"""
    read_newick(in_file::String; fossil = false)

Reads a newick tree into `sT` if fossil is false and `sTf` if fossil
is true from `in_file`.
"""
function read_newick(in_file::String; fossil = false)

  io = open(in_file)
  s = readlines(io)[1]
  close(io)

  s = s[2:(findfirst(isequal(';'), s)-2)]

  np = 0     # number of parenthesis
  nlp = 0    # number of left parenthesis yet to be closed
  nlin = 1   # number of distinct lineages at the origin
  nbif = 0   # number of bifurcations (= speciations)
  for (i,v) in enumerate(s)
    if v == '('
      nlp += 1
      np += 1
    elseif v == ')'
      nlp -= 1
      np += 1
    elseif v == ','
      nbif += 1
      if iszero(nlp) nlin += 1 end
    end
  end

  if fossil
    tree = from_string(s, np, nlp, nlin, nbif, sTf_label)
  else
    tree = from_string(s, np, nlp, nlin, nbif, sT_label)
  end

  return tree::Union{sT_label, sTf_label}
end




"""
    from_string(s      ::String,
                np     ::Int64,
                nlp    ::Int64,
                nlin   ::Int64,
                nbif   ::Int64,
                ::Type{sT_label})

Takes a string and turns it into a `sT_label` tree.
"""
function from_string(s      ::String,
                     np     ::Int64,
                     nlp    ::Int64,
                     nlin   ::Int64,
                     nbif   ::Int64,
                     ::Type{sT_label})

  # if root (starts with one stem lineage)
  if isone(nlin)
    tree = _from_string(s, sT_label)
  # if no root (starts with two crown lineage)
  else
    nrp = 0
    nlp = 0
    ci  = 0
    for (i,v) in enumerate(s)
      if v == '('
        nlp += 1
      elseif v == ')'
        nrp += 1
      elseif v == ',' && nlp == nrp
        ci = i
        break
      end
    end

    s1 = s[1:(ci-1)]
    s2 = s[(ci+1):end]

    tree = sT_label(_from_string(s1, sT_label),
                    _from_string(s2, sT_label),
                    0.0, "")
  end

  return tree
end




"""
    from_string(s   ::String,
                np  ::Int64,
                nlp ::Int64,
                nlin::Int64,
                nbif::Int64,
                ::Type{sTf_label})

Takes a string and turns it into a `sTf_label` tree.
"""
function from_string(s   ::String,
                     np  ::Int64,
                     nlp ::Int64,
                     nlin::Int64,
                     nbif::Int64,
                     ::Type{sTf_label})

  # if root (starts with one stem lineage)
  if isone(nlin)
    tree = _from_string(s, sTf_label)
  # if no root (starts with two crown lineage)
  else
    nrp = 0
    nlp = 0
    ci  = 0
    for (i,v) in enumerate(s)
      if v == '('
        nlp += 1
      elseif v == ')'
        nrp += 1
      elseif v == ',' && nlp == nrp
        ci = i
        break
      end
    end

    s1 = s[1:(ci-1)]
    s2 = s[(ci+1):end]

    tree = sTf_label(_from_string(s1, sTf_label),
                     _from_string(s2, sTf_label),
                     0.0, "")
  end

  # fossilize tips
  fossilizepasttips!(tree)

  return tree
end




"""
    _from_string(s::String, ::Type{T}) where {T <: sT}

Returns a tree of type `T` from newick string.
"""
function _from_string(s::String, ::Type{T}) where {T <: sT}

  # find pendant edge
  wd  = findlast(isequal(':'), s)
  pei = parse(Float64, s[(wd+1):end])
  lab = s[1:(wd-1)]

  # if tip
  if isone(count(isequal(':'), s))
      return T(pei, lab)
  else
    while s[wd-1] != ')'
      wd -= 1
    end

    s = s[2:(wd-2)]

    nrp = 0
    nlp = 0
    ci  = 0
    for (i,v) in enumerate(s)
      if v == '('
        nlp += 1
      elseif v == ')'
        nrp += 1
      elseif v == ',' && nlp == nrp
        ci = i
        break
      end
    end
  end

  s1 = s[1:(ci-1)]
  s2 = s[(ci+1):end]

  if isempty(s1)
    return T(_from_string(s2, T), pei, "")
  elseif isempty(s2)
    return T(_from_string(s1, T), pei, "")
  else
    return T(_from_string(s1, T), _from_string(s2, T), pei, "")
  end
end





"""
    write_newick(tree::T, out_file::String)

Writes `iTsimple` as a newick tree to `out_file`.
"""
function write_newick(tree::T, out_file::String) where {T <: iTree}

  s = to_string(tree)

  # if no stem branch
  if last(s, 4) == ":0.0"
    s = s[1:(end-4)]*";"
  else
    s = string("(", s, ");")
  end

  io = open(out_file*".tre", "w")
  write(io, s)
  close(io)

  return nothing
end




"""
    to_string(tree::T; n::Int64 = 0) where {T <: iTree})

Returns newick string.
"""
to_string(tree::T) where {T <: iTree} = _to_string(tree, 0)




"""
    _to_string(tree::T; n::Int64 = 0) where {T <: iTree})

Returns newick string.
"""
function _to_string(tree::T, n::Int64) where {T <: iTree}

  if istip(tree)
    return(string("t1:",e(tree)))
  end

  if istip(tree.d1)
    if istip(tree.d2)
      n += 1
      s1 = string("(t",n,":",e(tree.d1),",")
      n += 1
      s2 = string("t",n,":",e(tree.d2),"):",e(tree))
      return s1*s2
    else
      n += 1
      return string("(t",n,":",e(tree.d1), ",",
              _to_string(tree.d2, n),"):", e(tree))
    end
  elseif istip(tree.d2)
    n += 1
    return string("(", _to_string(tree.d1, n),
      ",t",n,":",e(tree.d2), "):", e(tree))
  else
    return string("(",_to_string(tree.d1, n),",",
               _to_string(tree.d2, ntips(tree.d1) + n),"):",e(tree))
  end
end




"""
    to_string(tree::T; n::Int64 = 0) where {T <: iTree})

Returns newick string.
"""
to_string(tree::T) where {T <: sTf} = _to_string(tree, 0, 0)




"""
    _to_string(tree::T, n::Int64, sa::Int64) where {T <: sTf}

Returns newick string.
"""
function _to_string(tree::T, n::Int64, sa::Int64) where {T <: sTf}

  if istip(tree)
    return(string("t1:",e(tree)))
  end

  if def1(tree) && def2(tree)
    if istip(tree.d1)
      if istip(tree.d2)
        n += 1
        s1 = string("(t",n,":",e(tree.d1),",")
        n += 1
        s2 = string("t",n,":",e(tree.d2),"):",e(tree))
        return s1*s2
      else
        n += 1
        return string("(t",n,":",e(tree.d1), ",",
                _to_string(tree.d2, n, sa),"):", e(tree))
      end
    elseif istip(tree.d2)
      n += 1
      return string("(", _to_string(tree.d1, n, sa),
        ",t",n,":",e(tree.d2), "):", e(tree))
    else
      return string("(",_to_string(tree.d1, n, sa),",",
                        _to_string(tree.d2, ntips(tree.d1) + n,
                                  nfossils(tree.d1) + sa),"):",
                        e(tree))
    end

  # sampled ancestors
  elseif def1(tree)
    sa += 1
    if istip(tree.d1)
      n += 1
      return string("(t",n,":",e(tree.d1),")sa",sa,":", e(tree))
    else
      return string("(",_to_string(tree.d1, n, sa),")sa",sa,":",e(tree))
    end
  else
    sa += 1
    if istip(tree.d2)
      n += 1
      return string("(t",n,":",e(tree.d2),")sa",sa,":", e(tree))
    else
      return string("(",_to_string(tree.d2, n, sa),")sa",sa,":",e(tree))
    end
  end
end



