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
  ei  = parse(Float64, s[(wd+1):end])
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
    bp  = parse(Float64,x[1:(pix-1)])
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
    iread(in_file::String)

Read a tree file exported by insane
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

    ii = 1
    it = six
    for line in eachline(file)
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

      ii += 1
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
  i = findfirst('-', s)
  T = iTd[s[1:(i-1)]]

  return _iparse(s[(i+2):end-1], T)
end




"""
    _istring(tree::sTpb)

parse istring to `sTpb`.
"""
function _iparse(s::String, ::Type{sTpb})

  lp = findlast(')', s)

  # if tip
  if isnothing(lp)
    ci = findfirst(',', s)
    return sTpb(parse(Float64, s[1:(ci-1)]), long(s[ci+1]))
  end

  si = s[(lp+2):end]
  s  = s[1:lp]

  ci = find_ci(s)

  s1 = s[2:(ci-2)]
  s2 = s[(ci+2):(end-1)]

  ci = findfirst(',', si)

  return sTpb(_iparse(s1, sTpb),
              _iparse(s2, sTpb),
              parse(Float64, si[1:(ci-1)]), long(si[ci+1]))
end




"""
    _istring(tree::sTbd)

parse istring to `sTbd`.
"""
function _iparse(s::String, ::Type{sTbd})

  lp = findlast(')', s)

  # if tip
  if isnothing(lp)
    ci = findfirst(',', s)
    return sTbd(parse(Float64, s[1:(ci-1)]), long(s[ci+1]), long(s[ci+3]))
  end

  si = s[(lp+2):end]
  s  = s[1:lp]

  ci = find_ci(s)

  s1 = s[2:(ci-2)]
  s2 = s[(ci+2):(end-1)]

  ci = findfirst(',', si)

  return sTbd(_iparse(s1, sTbd),
              _iparse(s2, sTbd),
              parse(Float64, si[1:(ci-1)]), long(si[ci+1]), long(si[ci+3]))
end




"""
    _istring(tree::sTfbd)

parse istring to `sTfbd`.
"""
function _iparse(s::String, ::Type{sTfbd})

  lp = findlast(')', s)

  # if tip
  if isnothing(lp)
    ci = findfirst(',', s)
    return sTfbd(parse(Float64, s[1:(ci-1)]), 
             long(s[ci+1]), long(s[ci+3]), long(s[ci+5]))
  end

  si = s[(lp+2):end]
  s  = s[1:lp]

  ci = find_ci(s)

  s1 = s[2:(ci-2)]
  s2 = s[(ci+2):(end-1)]

  ci = findfirst(',', si)

  if isempty(s1)
    return sTfbd(_iparse(s2, sTfbd),
                 parse(Float64, si[1:(ci-1)]), 
                 long(si[ci+1]), long(si[ci+3]), long(si[ci+5]))
  else
    return sTfbd(_iparse(s1, sTfbd),
                 _iparse(s2, sTfbd),
                 parse(Float64, si[1:(ci-1)]), 
                 long(si[ci+1]), long(si[ci+3]), long(si[ci+5]))
  end
end




"""
    _iparse(s::String, ::Type{iTpb})

parse istring to `iTpb`.
"""
function _iparse(s::String, ::Type{iTpb})

  lp = findlast(')', s)

  # if tip
  if isnothing(lp)
    c1 = findfirst(',', s)
    c2 = findnext(',', s, c1 + 1)
    c3 = findnext(',', s, c2 + 1)
    c4 = findnext(',', s, c3 + 1)

    return iTpb(parse(Float64, s[1:c1-1]), 
                parse(Float64, s[c1+1:c2-1]),
                parse(Float64, s[c2+1:c3-1]),
                long(s[c3+1]), 
                _iparse_v(s[c4+1:end]))
  end

  si = s[(lp+2):end]
  s  = s[1:lp]

  ci = find_ci(s)

  s1 = s[2:(ci-2)]
  s2 = s[(ci+2):(end-1)]

  c1 = findfirst(',', si)
  c2 = findnext(',', si, c1 + 1)
  c3 = findnext(',', si, c2 + 1)
  c4 = findnext(',', si, c3 + 1)

  return iTpb(_iparse(s1, iTpb),
              _iparse(s2, iTpb),
              parse(Float64, si[1:c1-1]), 
              parse(Float64, si[c1+1:c2-1]),
              parse(Float64, si[c2+1:c3-1]),
              long(si[c3+1]), 
              _iparse_v(si[c4+1:end]))
end




"""
    _iparse(s::String, ::Type{T}) where {T <: iT}

parse istring to `iT`.
"""
function _iparse(s::String, ::Type{T}) where {T <: iT}

  lp = findlast(')', s)

  # if tip
  if isnothing(lp)
    c1 = findfirst(',', s)
    c2 = findnext(',', s, c1 + 1)
    c3 = findnext(',', s, c2 + 1)
    c4 = findnext(',', s, c3 + 1)
    c5 = findnext(',', s, c4 + 1)

    return T(parse(Float64, s[1:c1-1]), 
             parse(Float64, s[c1+1:c2-1]),
             parse(Float64, s[c2+1:c3-1]),
             long(s[c3+1]), 
             long(s[c4+1]),
             _iparse_v(s[c5+1:end]))
  end

  si = s[(lp+2):end]
  s  = s[1:lp]

  ci = find_ci(s)

  s1 = s[2:(ci-2)]
  s2 = s[(ci+2):(end-1)]

  c1 = findfirst(',', si)
  c2 = findnext(',', si, c1 + 1)
  c3 = findnext(',', si, c2 + 1)
  c4 = findnext(',', si, c3 + 1)
  c5 = findnext(',', si, c4 + 1)

  return T(_iparse(s1, T),
           _iparse(s2, T),
           parse(Float64, si[1:c1-1]), 
           parse(Float64, si[c1+1:c2-1]),
           parse(Float64, si[c2+1:c3-1]),
           long(si[c3+1]), 
           long(si[c4+1]),
           _iparse_v(si[c5+1:end]))
end




"""
    _iparse(s::String, ::Type{iTbd})

parse istring to `iTbd`.
"""
function _iparse(s::String, ::Type{iTbd})

  lp = findlast(')', s)

  # if tip
  if isnothing(lp)
    c1 = findfirst(',', s)
    c2 = findnext(',', s, c1 + 1)
    c3 = findnext(',', s, c2 + 1)
    c4 = findnext(',', s, c3 + 1)
    c5 = findnext(',', s, c4 + 1)
    c6 = findnext(']', s, c5 + 1)

    return iTbd(parse(Float64, s[1:c1-1]), 
                parse(Float64, s[c1+1:c2-1]),
                parse(Float64, s[c2+1:c3-1]),
                long(s[c3+1]), 
                long(s[c4+1]),
                _iparse_v(s[c5+1:c6]),
                _iparse_v(s[c6+2:end]))
  end

  si = s[(lp+2):end]
  s  = s[1:lp]

  ci = find_ci(s)

  s1 = s[2:(ci-2)]
  s2 = s[(ci+2):(end-1)]

  c1 = findfirst(',', si)
  c2 = findnext(',', si, c1 + 1)
  c3 = findnext(',', si, c2 + 1)
  c4 = findnext(',', si, c3 + 1)
  c5 = findnext(',', si, c4 + 1)
  c6 = findnext(']', si, c5 + 1)

  return iTbd(_iparse(s1, iTbd),
              _iparse(s2, iTbd),
              parse(Float64, si[1:c1-1]), 
              parse(Float64, si[c1+1:c2-1]),
              parse(Float64, si[c2+1:c3-1]),
              long(si[c3+1]), 
              long(si[c4+1]),
              _iparse_v(si[c5+1:c6]),
              _iparse_v(si[c6+2:end]))
end




"""
    _iparse(s::String, ::Type{iTfbd})

parse istring to `iTfbd`.
"""
function _iparse(s::String, ::Type{iTfbd})

  lp = findlast(')', s)

  # if tip
  if isnothing(lp)
    c1 = findfirst(',', s)
    c2 = findnext(',', s, c1 + 1)
    c3 = findnext(',', s, c2 + 1)
    c4 = findnext(',', s, c3 + 1)
    c5 = findnext(',', s, c4 + 1)
    c6 = findnext(',', s, c5 + 1)
    c7 = findnext(']', s, c6 + 1)

    return iTfbd(parse(Float64, s[1:c1-1]), 
                 parse(Float64, s[c1+1:c2-1]),
                 parse(Float64, s[c2+1:c3-1]),
                 long(s[c3+1]), 
                 long(s[c4+1]),
                 long(s[c5+1]),
                 _iparse_v(s[c6+1:c7]),
                 _iparse_v(s[c7+2:end]))
  end

  si = s[(lp+2):end]
  s  = s[1:lp]

  ci = find_ci(s)

  s1 = s[2:(ci-2)]
  s2 = s[(ci+2):(end-1)]

  c1 = findfirst(',', si)
  c2 = findnext(',', si, c1 + 1)
  c3 = findnext(',', si, c2 + 1)
  c4 = findnext(',', si, c3 + 1)
  c5 = findnext(',', si, c4 + 1)
  c6 = findnext(',', si, c5 + 1)
  c7 = findnext(']', si, c6 + 1)

  if isempty(s1)
    return iTfbd(_iparse(s2, iTfbd),
                 parse(Float64, si[1:c1-1]), 
                 parse(Float64, si[c1+1:c2-1]),
                 parse(Float64, si[c2+1:c3-1]),
                 long(si[c3+1]), 
                 long(si[c4+1]),
                 long(si[c5+1]),
                 _iparse_v(si[c6+1:c7]),
                 _iparse_v(si[c7+2:end]))
  else
    return iTfbd(_iparse(s1, iTfbd),
                 _iparse(s2, iTfbd),
                 parse(Float64, si[1:c1-1]), 
                 parse(Float64, si[c1+1:c2-1]),
                 parse(Float64, si[c2+1:c3-1]),
                 long(si[c3+1]), 
                 long(si[c4+1]),
                 long(si[c5+1]),
                 _iparse_v(si[c6+1:c7]),
                 _iparse_v(si[c7+2:end]))
  end
end




"""
    _iparse_v(s::String)

Parse a string into a Float64 vector.
"""
function _iparse_v(s::String)
  v  = Float64[]
  ci = findfirst(',', s)
  i  = 1
  while !isnothing(ci)
    push!(v, 
      parse(Float64, s[(i+1):(ci-1)]))
    i  = ci + 1
    ci = findnext(',', s, ci + 1)
  end
  push!(v, 
    parse(Float64, s[(i+1):(end-1)]))
end




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



