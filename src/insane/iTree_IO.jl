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

  nop  = 0    # number of open parenthesis yet to be closed
  stem = true # if stem tree
  for (i,v) in enumerate(s)
    if v === '('
      nop += 1
    elseif v === ')'
      nop -= 1
    elseif v === ','
      if iszero(nop) 
        stem = false 
      end
    end
  end

  if fossil
    tree = from_string(s, stem, sTf_label)
  else
    tree = from_string(s, stem, sT_label)
  end

  return tree::Union{sT_label, sTf_label}
end




"""
    from_string(s::String, stem::Bool, ::Type{sT_label})

Takes a string and turns it into a `sT_label` tree.
"""
function from_string(s::String, stem::Bool, ::Type{sT_label})

  # if root (starts with one stem lineage)
  if stem
    tree = _from_string(s, sT_label)
  # if no root (starts with two crown lineage)
  else
    nop = 0
    ci  = 0
    for (i,v) in enumerate(s)
      if v == '('
        nop += 1
      elseif v == ')'
        nop -= 1
      elseif v == ',' && iszero(nop)
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
    from_string(s::String, stem::Bool, ::Type{sTf_label})

Takes a string and turns it into a `sTf_label` tree.
"""
function from_string(s::String, stem::Bool, ::Type{sTf_label})

  # if root (starts with one stem lineage)
  if stem
    tree = _from_string(s, sTf_label)
  # if no root (starts with two crown lineage)
  else
    nop = 0
    ci  = 0
    for (i,v) in enumerate(s)
      if v === '('
        nop += 1
      elseif v === ')'
        nop -= 1
      elseif v === ',' && iszero(nop)
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
  ei  = parse(Float64, s[(wd+1):end])
  lp  = findlast(isequal(')'), s)

  # if tip
  if isnothing(lp)
    return T(ei, s[1:(wd-1)])
  else

    lab = s[(lp+1):(wd-1)]
    s   = s[2:(lp-1)]

    nop = 0
    ci  = 0
    for (i,v) in enumerate(s)
      if v === '('
        nop += 1
      elseif v === ')'
        nop -= 1
      elseif v === ',' && iszero(nop)
        ci = i
        break
      end
    end
  end

  s1 = s[1:(ci-1)]
  s2 = s[(ci+1):end]

  if isempty(s1)
    return T(_from_string(s2, T), ei, lab)
  elseif isempty(s2)
    return T(_from_string(s1, T), ei, lab)
  else
    return T(_from_string(s1, T), _from_string(s2, T), ei, lab)
  end
end



"""
    onlyone(s::String, c::Char)

Returns true if there is only one of 'c' in string `s`.
"""
function onlyone(s::String, c::Char)
  n = 0
  for i in s
    if i === c
      n += 1
      if n > 1
        return false
      end
    end
  end

  return true
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
to_string(tree::T) where {T <: iTree} = _to_string(tree, 0)[1]




"""
    _to_string(tree::T; n::Int64 = 0) where {T <: iTree})

Returns newick string.
"""
function _to_string(tree::T, n::Int64) where {T <: iTree}

  if def1(tree)
    s1, n = _to_string(tree.d1, n)
    s2, n = _to_string(tree.d2, n)

    return string("(",s1,",", s2,"):",e(tree)), n
  else
    n += 1

    return string("t",n,":",e(tree)), n 
  end
end




"""
    to_string(tree::T; n::Int64 = 0) where {T <: iTree})

Returns newick string.
"""
to_string(tree::T) where {T <: sTf} = _to_string(tree, 0, 0)[1]




"""
    _to_string(tree::T, n::Int64, nf::Int64) where {T <: sTf}

Returns newick string.
"""
function _to_string(tree::T, n::Int64, nf::Int64) where {T <: sTf}

  if def1(tree)
    s1, n, nf = _to_string(tree.d1, n, nf)

    if def2(tree)
      s2, n, nf = _to_string(tree.d2, n, nf)
      s = string("(",s1,",", s2,"):",e(tree))
    else
      nf += 1
      s = string("(",s1,")f",nf,":", e(tree))
    end

    return s, n, nf
  else
    if isfossil(tree)
      nf += 1
      return string("f",nf,":",e(tree)), n, nf
    else
      n += 1
      return string("t",n,":",e(tree)), n, nf
    end
  end
end



