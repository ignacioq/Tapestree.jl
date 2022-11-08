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
function read_newick(in_file::String; fossil::Bool = false)

  io = open(in_file)
  s = readlines(io)[1]
  close(io)

  # if 1 tree
  if onlyone(s, ';')
    return _parse_newick(s, fossil)
  # if more than 1 tree
  else
    allsc = findall(';', s)

    t1 = _parse_newick(s[1:allsc[1]], fossil)
    tv = typeof(t1)[t1]

    for i in 2:lastindex(allsc)
      push!(tv, _parse_newick(s[(allsc[i-1] + 1):(allsc[i])], fossil))
    end

    return tv
  end
end




"""
    _parse_newick(in_file::String; fossil = false)

Reads a newick tree into `sT` if fossil is false and `sTf` if fossil
is true from `in_file`.
"""
function _parse_newick(s::String, fossil::Bool)

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
    return from_string(s, stem, sTf_label)
  else
    return from_string(s, stem, sT_label)
  end
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
    ci = find_ci(s)

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

  # if root
  if stem
    tree = _from_string(s, sTf_label)
  # if crown
  else
    ci = find_ci(s)

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
  ei  = Pparse(Float64, s[(wd+1):end])
  lp  = findlast(isequal(')'), s)

  # if tip
  if isnothing(lp)
    return T(ei, s[1:(wd-1)])
  else

    lab = s[(lp+1):(wd-1)]
    s   = s[2:(lp-1)]

    ci = find_ci(s)
  end

  s1 = s[1:(ci-1)]
  s2 = s[(ci+1):end]

  # if fossils are coded as 0 edge tip
  if last(s1, 4) == ":0.0" && onlyone(s1, ':')
    wd = findlast(isequal(':'), s1)
    return T(_from_string(s2, T), ei, s1[1:(wd-1)])
  elseif last(s2, 4) == ":0.0" && onlyone(s2, ':')
    wd = findlast(isequal(':'), s2)
    return T(_from_string(s1, T), ei, s2[1:(wd-1)])
  elseif isempty(s1)
    return T(_from_string(s2, T), ei, lab)
  elseif isempty(s2)
    return T(_from_string(s1, T), ei, lab)
  else
    return T(_from_string(s1, T), _from_string(s2, T), ei, lab)
  end
end




"""
    find_ci(s::String)

Find comma index in string within parentheses.
"""
function find_ci(s::String)
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

  return ci
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
to_string(tree::T) where {T <: iTf} = _to_string(tree, 0, 0)[1]




"""
    _to_string(tree::T, n::Int64, nf::Int64) where {T <: iTf}

Returns newick string.
"""
function _to_string(tree::T, n::Int64, nf::Int64) where {T <: iTf}

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




"""
    nsignif(x::String)

Return the number of significant digits in `x`, a string representing a number.
"""
function nsignif(x::String)

  pix = findfirst('.', x)

  if isnothing(pix)
    return lastindex(x)
  else
    bp  = Pparse(Float64,x[1:(pix-1)])
    # if less than 1
    if iszero(bp)
      l0 = findfirst(x -> x !== '0', x[(pix+1):end])
      if isnothing(l0)
        return 1
      else
        return lastindex(x[(pix+l0):end])
      end
    else
      return lastindex(x) - 1
    end
  end
end




"""
    istring(tree::T)

`sTbd` to istring.
"""
function istring(tree::T) where {T <: iTree}
  return string(T, '-', _istring(tree))
end




"""
    _istring(tree::sTpb)

`sTpb` to istring.
"""
function _istring(tree::sTpb)
  if def1(tree)
    return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
        e(tree), ',', 
        short(isfix(tree)), ')')
  else
    return string('(', 
        e(tree), ',', 
        short(isfix(tree)), ')')
  end
end




"""
    _istring(tree::sTbd)

`sTbd` to istring.
"""
function _istring(tree::sTbd)
  if def1(tree)
    return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
        e(tree), ',', 
        short(isextinct(tree)), ',', 
        short(isfix(tree)), ')')
  else
    return string('(', 
        e(tree), ',', 
        short(isextinct(tree)), ',', 
        short(isfix(tree)), ')')
  end
end




"""
    _istring(tree::sTfbd)

`sTfbd` to istring.
"""
function _istring(tree::sTfbd)
  if def1(tree)
    if def2(tree)
      return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
          e(tree), ',', 
          short(isextinct(tree)), ',', 
          "0,", 
          short(isfix(tree)), ')')
    else
      return string('(', _istring(tree.d1), ',', 
          e(tree), ',', 
          short(isextinct(tree)), ',', 
          "1,", 
          short(isfix(tree)), ')')
    end
  else
    return string('(', 
        e(tree), ',', 
        short(isextinct(tree)), ',', 
        short(isfossil(tree)), ',', 
        short(isfix(tree)), ')')
  end
end




"""
    _istring(tree::iTpb)

`iTpb` to istring.
"""
function _istring(tree::iTpb)
  if def1(tree)
    return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
             e(tree), ',',
             dt(tree), ',',
             fdt(tree), ',',
             short(isfix(tree)), ',', 
             lλ(tree), ')')
  else
    return string('(', 
             e(tree), ',',
             dt(tree), ',',
             fdt(tree), ',',
             short(isfix(tree)), ',', 
             lλ(tree), ')')
  end
end




"""
    _istring(tree::iT)

`iT` to istring.
"""
function _istring(tree::iT)
  if def1(tree)
    return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
             e(tree), ',',
             dt(tree), ',',
             fdt(tree), ',',
             short(isextinct(tree)), ',', 
             short(isfix(tree)), ',', 
             lλ(tree), ')')
  else
    return string('(', 
             e(tree), ',',
             dt(tree), ',',
             fdt(tree), ',',
             short(isextinct(tree)), ',', 
             short(isfix(tree)), ',', 
             lλ(tree), ')')
  end
end




"""
    _istring(tree::iTbd)

`iTbd` to istring.
"""
function _istring(tree::iTbd)
  if def1(tree)
    return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
             e(tree), ',',
             dt(tree), ',',
             fdt(tree), ',',
             short(isextinct(tree)), ',', 
             short(isfix(tree)), ',', 
             lλ(tree), ',', 
             lμ(tree), ')')
  else
    return string('(', 
             e(tree), ',',
             dt(tree), ',',
             fdt(tree), ',',
             short(isextinct(tree)), ',', 
             short(isfix(tree)), ',', 
             lλ(tree), ',', 
             lμ(tree), ')')
  end
end




"""
    _istring(tree::iTfbd)

`iTfbd` to istring.
"""
function _istring(tree::iTfbd)
  if def1(tree)
    if def2(tree)
      return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
               e(tree), ',',
               dt(tree), ',',
               fdt(tree), ',',
               short(isextinct(tree)), ',', 
               "0,", 
               short(isfix(tree)), ',', 
               lλ(tree), ',', 
               lμ(tree), ')')
    else
      return string('(', _istring(tree.d1), ',', 
               e(tree), ',',
               dt(tree), ',',
               fdt(tree), ',',
               short(isextinct(tree)), ',', 
               "1,", 
               short(isfix(tree)), ',', 
               lλ(tree), ',', 
               lμ(tree), ')')
    end
  else
    return string('(', 
             e(tree), ',',
             dt(tree), ',',
             fdt(tree), ',',
             short(isextinct(tree)), ',', 
             short(isfossil(tree)), ',', 
             short(isfix(tree)), ',', 
             lλ(tree), ',', 
             lμ(tree), ')')
  end
end




"""
    iread(file::String; ix::OrdinalRange{Int64,Int64} = 0:0)

Read a tree file exported by insane in `file` and with optional OrdinalRange
specifying which trees to sample.
"""
function iread(file::String;
               ix::OrdinalRange{Int64,Int64} = 0:0)

  # read all
  if iszero(ix[1])
    s  = readlines(file)
  # read according to ix
  else
    iix = first(ix)
    lix = last(ix)
    six = step(ix)

    s = String[]

    ii = 0
    it = six
    for line in eachline(file)
      ii += 1

      if ii < iix
        continue
      end

      # read trees
      if it === six
        push!(s, line)
        it = 0
      end

      if ii >= lix
        break
      end
      it += 1
    end
  end

  ls = lastindex(s)
  t0 = iparse(s[1])
  tv::Vector{typeof(t0)} = typeof(t0)[t0]

  if ls > 1
    for i in 2:ls
      push!(tv, iparse(s[i])::typeof(t0))
    end
  end

  return tv
end




"""
    iparse(s::String)

from istring to `iTree`.
"""
function iparse(s::String)
  ls = lastindex(s)
  i  = findfirst('-', s)
  st = SubString(s,1:(i-1))
  T  = iTd[st]
  si = s[i+2:ls-1]

  t0, ix = _iparse(si, 1, ls - lastindex(st) - 3, T)

  return t0 
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{sTpb})

parse istring to `sTpb`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{sTpb})

 @inbounds begin

    inode = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, sTpb)
      inode = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, sTpb)
    end

    i1 = findnext(',', s, i + 1)

    if inode
      tree = sTpb(sd1, sd2, 
                  Pparse(Float64, s[i:i1-1]), 
                  long(s[i1+1]))
    else
      tree = sTpb(Pparse(Float64, s[i:i1-1]), 
                  long(s[i1+1]))
    end

    i = i1 + 2

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{sTbd})

parse istring to `sTbd`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{sTbd})

  @inbounds begin

    inode = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, sTbd)
      inode = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, sTbd)
    end

    i1 = findnext(',', s, i + 1)

    if inode
      tree = sTbd(sd1, sd2, 
                  Pparse(Float64, s[i:i1-1]), 
                  long(s[i1+1]), long(s[i1+3]))
    else
      tree = sTbd(Pparse(Float64, s[i:i1-1]), 
                  long(s[i1+1]), long(s[i1+3]))
    end

    i = i1 + 4

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{sTfbd})

parse istring to `sTfbd`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{sTfbd})

  @inbounds begin

    in1 = false
    in2 = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, sTfbd)
      in1 = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, sTfbd)
      in2 = true
    end

    i1 = findnext(',', s, i + 1)

    if in1
      if in2
        tree = sTfbd(sd1, sd2, 
                    Pparse(Float64, s[i:i1-1]), 
                    long(s[i1+1]), long(s[i1+3]), long(s[i1+5]))
      else
        tree = sTfbd(sd1, Pparse(Float64, s[i:i1-1]), 
                    long(s[i1+1]), long(s[i1+3]), long(s[i1+5]))
      end
    else
      tree = sTfbd(Pparse(Float64, s[i:i1-1]), 
                  long(s[i1+1]), long(s[i1+3]), long(s[i1+5]))
    end

    i = i1 + 6

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{iTpb})

parse istring to `iTpb`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{iTpb})

  @inbounds begin

    inode = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, iTpb)
      inode = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, iTpb)
    end

    i1 = findnext(',', s, i  + 1)
    i2 = findnext(',', s, i1 + 1)
    i3 = findnext(',', s, i2 + 1)
    i4 = findnext(']', s, i3 + 1)

    if inode
      tree = iTpb(sd1, sd2,
                  Pparse(Float64, s[i:i1-1]),
                  Pparse(Float64, s[i1+1:i2-1]),
                  Pparse(Float64, s[i2+1:i3-1]),
                  long(s[i3+1]), 
                  _iparse_v(s[i3+4:i4-1]))
    else
      tree = iTpb(Pparse(Float64, s[i:i1-1]),
                  Pparse(Float64, s[i1+1:i2-1]),
                  Pparse(Float64, s[i2+1:i3-1]),
                  long(s[i3+1]), 
                  _iparse_v(s[i3+4:i4-1]))
    end

    i = i4 + 1

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{T}) where {T <: iT}

parse istring to `iT`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{T}) where {T <: iT}

  @inbounds begin

    inode = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, T)
      inode = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, T)
    end

    i1 = findnext(',', s, i  + 1)
    i2 = findnext(',', s, i1 + 1)
    i3 = findnext(',', s, i2 + 1)
    i4 = findnext(']', s, i3 + 1)

    if inode
      tree = T(sd1, sd2,
               Pparse(Float64, s[i:i1-1]),
               Pparse(Float64, s[i1+1:i2-1]),
               Pparse(Float64, s[i2+1:i3-1]),
               long(s[i3+1]), 
               long(s[i3+3]), 
               _iparse_v(s[i3+6:i4-1]))
    else
      tree = T(Pparse(Float64, s[i:i1-1]),
               Pparse(Float64, s[i1+1:i2-1]),
               Pparse(Float64, s[i2+1:i3-1]),
               long(s[i3+1]), 
               long(s[i3+3]), 
               _iparse_v(s[i3+6:i4-1]))
    end

    i = i4 + 1

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{iTbd})

parse istring to `iT`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{iTbd})

  @inbounds begin

    inode = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, iTbd)
      inode = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, iTbd)
    end

    i1 = findnext(',', s, i  + 1)
    i2 = findnext(',', s, i1 + 1)
    i3 = findnext(',', s, i2 + 1)
    i4 = findnext(']', s, i3 + 1)
    i5 = findnext(']', s, i4 + 1)

    if inode
      tree = iTbd(sd1, sd2,
                  Pparse(Float64, s[i:i1-1]),
                  Pparse(Float64, s[i1+1:i2-1]),
                  Pparse(Float64, s[i2+1:i3-1]),
                  long(s[i3+1]), 
                  long(s[i3+3]), 
                  _iparse_v(s[i3+6:i4-1]),
                  _iparse_v(s[i4+3:i5-1]))
    else
      tree = iTbd(Pparse(Float64, s[i:i1-1]),
                  Pparse(Float64, s[i1+1:i2-1]),
                  Pparse(Float64, s[i2+1:i3-1]),
                  long(s[i3+1]), 
                  long(s[i3+3]), 
                  _iparse_v(s[i3+6:i4-1]),
                  _iparse_v(s[i4+3:i5-1]))
    end

    i = i5 + 1

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{iTfbd})

parse istring to `iT`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{iTfbd})

  @inbounds begin

    in1 = false
    in2 = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, iTfbd)
      in1 = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, iTfbd)
      in2 = true
    end

    i1 = findnext(',', s, i  + 1)
    i2 = findnext(',', s, i1 + 1)
    i3 = findnext(',', s, i2 + 1)
    i4 = findnext(']', s, i3 + 1)
    i5 = findnext(']', s, i4 + 1)

    if in1
      if in2
        tree = iTfbd(sd1, sd2,
                     Pparse(Float64, s[i:i1-1]),
                     Pparse(Float64, s[i1+1:i2-1]),
                     Pparse(Float64, s[i2+1:i3-1]),
                     long(s[i3+1]), 
                     long(s[i3+3]), 
                     long(s[i3+5]), 
                     _iparse_v(s[i3+8:i4-1]),
                     _iparse_v(s[i4+3:i5-1]))
      else
        tree = iTfbd(sd1,
                     Pparse(Float64, s[i:i1-1]),
                     Pparse(Float64, s[i1+1:i2-1]),
                     Pparse(Float64, s[i2+1:i3-1]),
                     long(s[i3+1]), 
                     long(s[i3+3]), 
                     long(s[i3+5]), 
                     _iparse_v(s[i3+8:i4-1]),
                     _iparse_v(s[i4+3:i5-1]))
      end
    else
      tree = iTfbd(Pparse(Float64, s[i:i1-1]),
                   Pparse(Float64, s[i1+1:i2-1]),
                   Pparse(Float64, s[i2+1:i3-1]),
                   long(s[i3+1]), 
                   long(s[i3+3]), 
                   long(s[i3+5]), 
                   _iparse_v(s[i3+8:i4-1]),
                   _iparse_v(s[i4+3:i5-1]))
    end

    i = i5 + 1

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse_v(s::String)

Parse a string into a `Float64` vector.
"""
_iparse_v(s::String) = Parsers.parse.(Float64, split(s, ','))




"""
    short(x::Bool)

Return 0 or 1 for false or true
"""
short(x::Bool) = x ? '1' : '0'




"""
    long(x::Char)

Return 0 or 1 for false or true
"""
long(x::Char) = x === '1' ? true : false



