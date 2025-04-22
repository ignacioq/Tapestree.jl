#=

insane tree read and write

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#




"""
    read_newick(in_file::String; ix::OrdinalRange{Int64,Int64} = 0:0)

Reads a newick tree into `sT` from `in_file` at lines `ix`.
"""
function read_newick(in_file::String; ix::OrdinalRange{Int64,Int64} = 0:0)

  io = open(in_file)

  iix = first(ix)
  lix = iszero(ix[1]) ? typemax(Int64) : last(ix)
  six = step(ix)

  tv = sT_label[]

  ii = 0
  it = six
  for line in eachline(io)

    iszero(lastindex(line)) && continue

    # read trees
    if it === six
      if onlyone(line, ';')
        ii += 1
        ii < iix && continue
        push!(tv, _parse_newick(line))
      else
        allsc = findall(';', line)
        pushfirst!(allsc, 0)

        for i in Base.OneTo(lastindex(allsc)-1)
          ii += 1
          ii < iix && continue
          if it === six
            push!(tv, _parse_newick(line[(allsc[i] + 1):(allsc[i+1])]))
            it = 0
          end
          ii >= lix && break
          it += 1
        end
      end
      it = 0
    end

    ii >= lix && break
    it += 1
  end

  close(io)

  if isone(lastindex(tv))
    return tv[1]
  else
    return tv
  end
end




"""
    _parse_newick(in_file::String)

Reads a newick tree into `sT` from `in_file`.
"""
function _parse_newick(s::String)

  s = s[2:(findlast(isequal(';'), s)-1)]

  # find break if crown tree
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

  tree = from_string(s, ci, sT_label)

  return tree
end




"""
    read_newick(in_file::String,
                fossil ::Bool;
                ix     ::OrdinalRange{Int64,Int64} = 0:0,
                ne     ::Float64                   = accer)

Reads a newick tree into `sTf` from `in_file` at lines `ix`.
"""
function read_newick(in_file::String,
                     fossil ::Bool;
                     ix     ::OrdinalRange{Int64,Int64} = 0:0,
                     ne     ::Float64                   = accerr)

  io = open(in_file)

  iix = first(ix)
  lix = iszero(iix) ? typemax(Int64) : last(ix)
  six = step(ix)

  tv = sTf_label[]

  ii = 0
  it = six
  for line in eachline(io)

    iszero(lastindex(line)) && continue

    # read trees
    if it === six
      if onlyone(line, ';')
        ii += 1
        ii < iix && continue
        push!(tv, _parse_newick(line, ne))
      else
        allsc = findall(';', line)
        pushfirst!(allsc, 0)

        for i in Base.OneTo(lastindex(allsc)-1)
          ii += 1
          ii < iix && continue
          if it === six
            push!(tv, _parse_newick(line[(allsc[i] + 1):(allsc[i+1])], ne))
            it = 0
          end
          ii >= lix && break
          it += 1
        end
      end
      it = 0
    end

    ii >= lix && break
    it += 1
  end

  close(io)

  if isone(lastindex(tv))
    return tv[1]
  else
    return tv
  end
end




"""
    _parse_newick(in_file::String; fossil = false)

Reads a newick tree into `sT` if fossil is false and `sTf` if fossil
is true from `in_file`.
"""
function _parse_newick(s::String, ne::Float64)

  s = s[2:(findlast(isequal(';'), s)-1)]

  # find break if crown tree
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

  tree = from_string(s, ci, sTf_label)
  fossilizepasttips!(tree, ne)

  return tree
end




"""
    from_string(s::String, stem::Bool, ::Type{sT_label})

Takes a string and turns it into a `sT_label` tree.
"""
function from_string(s::String, ci::Int64, ::Type{T}) where {T <: sT}

  # if root (starts with one stem lineage)
  if iszero(ci)
    tree, i = _from_string(s, 1, T)
  # if no root (starts with two crown lineage)
  else

    sd1, i = _from_string(s[1:ci],       1, T)
    sd2, i = _from_string(s[(ci+1):end], 1, T)

    if e(sd1) === 0.0
      tree = T(sd2, 0.0, label(sd1))
    elseif e(sd2) === 0.0
      tree = T(sd1, 0.0, label(sd2))
    else
      tree = T(sd1, sd2, 0.0, "")
    end
  end

  return tree
end



"""
    _from_string(s::String, ::Type{T}) where {T <: sT}

Returns a tree of type `T` from newick string.
"""
function _from_string(s::String, i::Int64, ::Type{T}) where {T <: sT}

  @inbounds begin

    in1 = false
    in2 = false

    if s[i] === '('
      sd1, i = _from_string(s, i + 1, T)
      in1 = true
    end

    if s[i] === ','
      sd2, i = _from_string(s, i + 1, T)
      in2 = true
    end

    i1 = findnext(':', s, i)
    i2 = find_cp(s, i1 + 1)

    if in1
      if in2
        if e(sd1) === 0.0
          tree = T(sd2, Pparse(Float64, s[i1+1:i2-1]), label(sd1))
        elseif e(sd2) === 0.0
          tree = T(sd1, Pparse(Float64, s[i1+1:i2-1]), label(sd2))
        else
          tree = T(sd1, sd2, Pparse(Float64, s[i1+1:i2-1]), s[i:i1-1])
        end
      else
        tree = T(sd1, Pparse(Float64, s[i1+1:i2-1]), s[i:i1-1])
      end
    else
      tree = T(Pparse(Float64, s[i1+1:i2-1]), s[i:i1-1])
    end

    i = i2

    while s[i] === ')'
      i += 1
    end
  end

  return tree, i
end




"""
    find_cp(s::String, i::Int64)

Find next ',' or ')' after index `i`.
"""
function find_cp(s::String, i::Int64)

  f1 = findnext(',', s, i)
  f2 = findnext(')', s, i)

  if isnothing(f1)
    if isnothing(f2)
      return lastindex(s)
    else 
      return f2
    end
  elseif isnothing(f2)
    return f1
  else
    return min(f1, f2)
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
    write_newick(tree::T, ofile::String)

Writes an `iTree` as a newick tree to `ofile`.
"""
function write_newick(tree::T, ofile::String) where {T <: iTree}

  io = IOBuffer()
  ic = iszero(e(tree))
  !ic && write(io, '(')
  nw_buffer(io, tree, ic)
  !ic && write(io, ')')
  write(io, ';')
  write(ofile*".tre", take!(io))

  return nothing
end




"""
    write_newick(treev::Vector{T}, ofile::String)

Writes an `iTree` as a newick tree to `ofile`.
"""
function write_newick(treev::Vector{T}, ofile::String) where {T <: iTree}

  to = open(ofile*".trees", "w")
  io = IOBuffer()

  for t in treev
    ic = iszero(e(t))
    !ic && write(io, '(')
    nw_buffer(io, t, ic)
    !ic && write(io, ')')
    write(io, ';', '\n')
    write(to, take!(io))
  end

  close(to)

  return nothing
end




"""
    nw_buffer(io::IOBuffer, tree::T) where {T <: iTree})

Writes an `iTree` to IOBuffer `io`.
"""
nw_buffer(io::IOBuffer, tree::T, ic::Bool) where {T <: iTree} = 
  _nw_buffer(io, tree, 0, ic)

"""
    _nw_buffer(io::IOBuffer, tree::T, n::Int64) where {T <: iTree})

Writes an `iTree` to IOBuffer `io`.
"""
function _nw_buffer(io::IOBuffer, tree::T, n::Int64, ic::Bool) where {T <: iTree}

  if def1(tree)
    write(io, '(')
    n = _nw_buffer(io, tree.d1, n, false)
    write(io, ',')
    n = _nw_buffer(io, tree.d2, n, false)
    write(io, ')')
    !ic && print(io, ':', e(tree))
  else
    n += 1
    print(io, 't', n, ':', e(tree))
  end

  return n
end




"""
    nw_buffer(io::IOBuffer, tree::T) where {T <: iTree})

Writes a fossil tree `uTf` to IOBuffer `io`.
"""
nw_buffer(io::IOBuffer, tree::T, ic::Bool) where {T <: uTf} = 
  _nw_buffer(io, tree, 0, 0, ic)

"""
    _nw_buffer(io. ::IOBuffer, 
               tree::T, 
               n   ::Int64, 
               nf  ::Int64, 
               ic  ::Bool) where {T <: uTf}

Writes a fossil tree `uTf` to IOBuffer `io`.
"""
function _nw_buffer(io  ::IOBuffer, 
                    tree::T, 
                    n   ::Int64, 
                    nf  ::Int64, 
                    ic  ::Bool) where {T <: uTf}

  if def1(tree)
    write(io, '(')
    n, nf = _nw_buffer(io, tree.d1, n, nf, false)

    if def2(tree)
      write(io, ',')
      n, nf = _nw_buffer(io, tree.d2, n, nf,  false)
      write(io, ')')
      !ic && print(io, ':', e(tree))
    else
      nf += 1
      print(io, ")f", nf, ':', e(tree))
    end
  else
    if isfossil(tree)
      nf += 1
      print(io, 'f', nf, ':', e(tree))
    else
      n += 1
      print(io, 't', n, ':', e(tree))
    end
  end

  return n, nf
end




"""
    nw_buffer(io::IOBuffer, tree::T)

Writes a labelled tree `Tlabel` to IOBuffer `io`.
"""
nw_buffer(io::IOBuffer, tree::T, ic::Bool) where {T <: Tlabel} = 
  _nw_buffer(io, tree, ic)

"""
    _nw_buffer(io::IOBuffer, tree::T) where {T <: Tlabel}

Writes a labelled tree `Tlabel` to IOBuffer `io`.
"""
function _nw_buffer(io::IOBuffer, tree::T, ic::Bool) where {T <: Tlabel}

  if def1(tree)
    write(io, '(')
    _nw_buffer(io, tree.d1, false)

    if def2(tree)
      write(io, ',')
      _nw_buffer(io, tree.d2, false)
      write(io, ')')
      !ic && print(io, ':', e(tree))
    else
      print(io, ")", label(tree), ':', e(tree))
    end
  else
    print(io, label(tree), ':', e(tree))
  end
end




"""
    write_nexus(tree::T, reftree::sT_label, ofile::String) where {T <: iTree}

Writes an `iTree` as a extensive nexus tree to `ofile`.
"""
function write_nexus(tree::T, reftree::Tl, ofile::String) where {T <: iTree, Tl <: Tlabel}

  io = IOBuffer()
  write(io, "#NEXUS\n\nBegin trees;\ntree 1 = ")
  ic = iszero(e(tree))
  !ic && write(io, '(')
  nx_buffer(io, tree, reftree, ic)
  !ic && write(io, ')')
  write(io, ";\nEnd;")

  write(ofile*".nex", take!(io))

  return nothing
end




"""
    write_nexus(treev::Vector{T}, reftree::sT_label, ofile::String) where {T <: iTree}

Writes an `iTree` as a extensive nexus tree to `ofile`.
"""
function write_nexus(treev::Vector{T}, 
                     reftree::Tl, 
                     ofile::String) where {T <: iTree, Tl <: Tlabel}

  to = open(ofile*".nex", "w")
  io = IOBuffer()
  write(io, "#NEXUS\n\nBegin trees;\n")

  for (i,t) in enumerate(treev)
    print(io, "tree ", i, " = ")
    ic = iszero(e(t))
    !ic && write(io, '(')
    nx_buffer(io, t, reftree, ic)
    !ic && write(io, ')')
    write(io, ';', '\n')
  end

  write(io, "End;")
  write(to, take!(io))
  close(to)

  return nothing
end




"""
    nx_buffer(io::IOBuffer, tree::T) where {T <: iTree})

Writes an `iTree` to IOBuffer `io`.
"""
nx_buffer(io::IOBuffer, tree::T, reftree::Tl, ic::Bool) where {T <: iTree, Tl <: Tlabel} = 
  _nx_buffer(io, tree, reftree, 0, ic)

"""
    _nx_buffer(io     ::IOBuffer, 
               tree   ::T, 
               reftree::sT_label, 
               n      ::Int64, 
               ic     ::Bool) where {T <: iT}

Writes an `iTree` to IOBuffer `io`.
"""
function _nx_buffer(io     ::IOBuffer, 
                    tree   ::T, 
                    reftree::sT_label, 
                    n      ::Int64, 
                    ic     ::Bool) where {T <: iT}

  if def1(tree)
    write(io, '(')
    if isfix(tree.d1) && isfix(tree.d2)
      n = _nx_buffer(io, tree.d1, reftree.d1, n, false)
      write(io, ',')
      n = _nx_buffer(io, tree.d2, reftree.d2, n, false)
    else
      n = _nx_buffer(io, tree.d1, reftree, n, false)
      write(io, ',')
      n = _nx_buffer(io, tree.d2, reftree, n, false)
    end
    write(io, ")[&sr=")
    nx_printv(io, lλ(tree))
    print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), ",da=", !isfix(tree), ']')
    !ic && print(io, ':', e(tree))
  else
    if isfix(tree)
      print(io, label(reftree))
    else
      n += 1
      print(io, 't', n)
    end
    write(io, "[&sr=")
    nx_printv(io, lλ(tree))
    print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), ",da=", !isfix(tree), ']', 
      ':', e(tree))
  end

  return n
end




"""
    _nx_buffer(io     ::IOBuffer, 
               tree   ::iTbd, 
               reftree::sT_label, 
               n      ::Int64, 
               ic     ::Bool)

Writes an `iTree` to IOBuffer `io`.
"""
function _nx_buffer(io     ::IOBuffer, 
                    tree   ::iTbd, 
                    reftree::sT_label, 
                    n      ::Int64, 
                    ic     ::Bool)

  if def1(tree)
    write(io, '(')
    if isfix(tree.d1) && isfix(tree.d2)
      n = _nx_buffer(io, tree.d1, reftree.d1, n, false)
      write(io, ',')
      n = _nx_buffer(io, tree.d2, reftree.d2, n, false)
    else
      n = _nx_buffer(io, tree.d1, reftree, n, false)
      write(io, ',')
      n = _nx_buffer(io, tree.d2, reftree, n, false)
    end
    write(io, ")[&sr=")
    nx_printv(io, lλ(tree))
    write(io, ",er=")
    nx_printv(io, lμ(tree))
    print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), ",da=", !isfix(tree), ']')
    !ic && print(io, ':', e(tree))
  else
    if isfix(tree)
      print(io, label(reftree))
    else
      n += 1
      print(io, 't', n)
    end
    write(io, "[&sr=")
    nx_printv(io, lλ(tree))
    write(io, ",er=")
    nx_printv(io, lμ(tree))
    print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), ",da=", !isfix(tree), ']', 
      ':', e(tree))
  end

  return n
end




"""
    nx_buffer(io::IOBuffer, tree::iTfbd, reftree::sTf_label, ic::Bool)

Writes an `iTree` to IOBuffer `io`.
"""
nx_buffer(io::IOBuffer, tree::iTfbd, reftree::sTf_label, ic::Bool) = 
  _nx_buffer(io, tree, reftree, 0, 0, ic)

"""
    _nx_buffer(io     ::IOBuffer, 
               tree   ::iTfbd, 
               reftree::sTf_label, 
               n      ::Int64, 
               ic     ::Bool)

Writes an `iTree` to IOBuffer `io`.
"""
function _nx_buffer(io     ::IOBuffer, 
                    tree   ::iTfbd, 
                    reftree::sTf_label, 
                    n      ::Int64, 
                    nf     ::Int64, 
                    ic     ::Bool)

  if def1(tree)
    write(io, '(')
    if def2(tree)
      if isfix(tree.d1) && isfix(tree.d2)
        n, nf = _nx_buffer(io, tree.d1, reftree.d1, n, nf, false)
        write(io, ',')
        n, nf = _nx_buffer(io, tree.d2, reftree.d2, n, nf, false)
      else
        n, nf = _nx_buffer(io, tree.d1, reftree, n, nf, false)
        write(io, ',')
        n, nf = _nx_buffer(io, tree.d2, reftree, n, nf, false)
      end
      write(io, ")[&sr=")
      nx_printv(io, lλ(tree))
      write(io, ",er=")
      nx_printv(io, lμ(tree))
      print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), ",da=", !isfix(tree), ']')
      !ic && print(io, ':', e(tree))
    else
      if isfix(tree.d1)
        n, nf = _nx_buffer(io, tree.d1, reftree.d1, n, nf, false)
      else
        n, nf = _nx_buffer(io, tree.d1, reftree, n, nf, false)
      end
      write(io, ')')
      if isfix(tree)
        print(io, label(reftree))
      else
        nf += 1
        print(io, 'f', nf)
      end
      write(io, "[&sr=")
      nx_printv(io, lλ(tree))
      write(io, ",er=")
      nx_printv(io, lμ(tree))
      print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), ",da=", !isfix(tree), 
        "]:", e(tree))
    end

  else
    if isfix(tree)
      print(io, label(reftree))
    else
      n += 1
      print(io, 't', n)
    end
    write(io, "[&sr=")
    nx_printv(io, lλ(tree))
    write(io, ",er=")
    nx_printv(io, lμ(tree))
    print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), ",da=", !isfix(tree), ']', 
      ':', e(tree))
  end

  return n, nf
end




"""
    nx_printv(io::IOBuffer, x::Vector{Float64}) 

Print vector for nexus format
"""
function nx_printv(io::IOBuffer, x::Vector{Float64}) 
  write(io, '{')
  for xi in x
    print(io, xi, ',')
  end
  write(io, '}')
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
    iwrite(tree::T, ofile::String) where {T <: iTree}

Write iTree to file.
"""
function iwrite(tree::T, ofile::String) where {T <: iTree}
  io = IOBuffer()
  ibuffer(io, tree)
  write(string(ofile, ".txt"), take!(io))
  return nothing
end




"""
    iwrite(tree::Vector{T}, ofile::String) where {T <: iTree}

Write a vector of trees.
"""
function iwrite(tree::Vector{T}, ofile::String) where {T <: iTree}
  open(ofile*".txt", "w") do to
    for tri in tree
      io = IOBuffer()
      ibuffer(io, tri)
      write(io, '\n')
      write(to, take!(io))
      flush(to)
    end
  end
end




"""
    ibuffer(tree::T, io::IOBuffer)

Write iTree to IOBuffer.
"""
function ibuffer(io::IOBuffer, tree::T) where {T <: iTree}
  write(io, string(T), '-') 
  _ibuffer(io, tree)
end




"""
    _ibuffer(io::IOBuffer, tree::sTpb)

Write `sTpb` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::sTpb)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1)
    write(io, ',')
    _ibuffer(io, tree.d2)
    print(io, ',', e(tree), ',', short(isfix(tree)), ')')
  else
    print(io, '(', e(tree), ',', short(isfix(tree)), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::sTbd)

Write `sTbd` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::sTbd)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1)
    write(io, ',')
    _ibuffer(io, tree.d2)
    print(io, ',', e(tree), ',', 
          short(isextinct(tree)), ',', 
          short(isfix(tree)), ')')
  else
    print(io, '(', e(tree), ',', 
        short(isextinct(tree)), ',', 
        short(isfix(tree)), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::sTfbd)

Write `sTfbd` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::sTfbd)
  if def1(tree)
    write(io, '(')
    if def2(tree)
      _ibuffer(io, tree.d1), 
      write(io, ',') 
      _ibuffer(io, tree.d2), 
      print(io, ',', e(tree), ',', 
          short(isextinct(tree)), ',', "0,", 
          short(isfix(tree)), ')')
    else
      _ibuffer(io, tree.d1), 
      print(io, ',', e(tree), ',', 
          short(isextinct(tree)), ',', "1,", 
          short(isfix(tree)), ')')
    end
  else
    print(io, '(', e(tree), ',', 
        short(isextinct(tree)), ',', 
        short(isfossil(tree)), ',',
        short(isfix(tree)), ')')
  end
end



"""
    _istring(tree::sTxs)

`sTxs` to istring.
"""
function _istring(tree::sTxs)
  if def1(tree)
    if def2(tree)
      return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
          e(tree), ',', 
          dt(tree), ',',
          fdt(tree), ',',
          xv(tree), ',',
          lσ2(tree), ')')
    else
      return string('(', _istring(tree.d1), ',', 
          e(tree), ',', 
          dt(tree), ',',
          fdt(tree), ',',
          xv(tree), ',',
          lσ2(tree), ')')
    end
  else
    return string('(', 
          e(tree), ',', 
          dt(tree), ',',
          fdt(tree), ',',
          xv(tree), ',',
          lσ2(tree), ')')
  end
end



"""
    _ibuffer(io::IOBuffer, tree::iTpb)

Write `iTpb` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::iTpb)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1), 
    write(io, ',')
    _ibuffer(io, tree.d2), 
    print(io, ',', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          short(isfix(tree)), ',', lλ(tree), ')')
  else
    print(io, '(', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          short(isfix(tree)), ',', lλ(tree), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::iT)

Write `iT` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::iT)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1), 
    write(io, ',')
    _ibuffer(io, tree.d2), 
    print(io, ',', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          short(isextinct(tree)), ',', short(isfix(tree)), ',', 
          lλ(tree), ')')
  else
    print(io, '(', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          short(isextinct(tree)), ',', short(isfix(tree)), ',', 
          lλ(tree), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::iTbd)

Write `iTbd` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::iTbd)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1), 
    write(io, ',')
    _ibuffer(io, tree.d2), 
    print(io, ',', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          short(isextinct(tree)), ',',  short(isfix(tree)), ',', 
          lλ(tree), ',', lμ(tree), ')')
  else
    print(io, '(', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          short(isextinct(tree)), ',',  short(isfix(tree)), ',', 
          lλ(tree), ',', lμ(tree), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::iTfbd)

Write `iTfbd` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::iTfbd)
  if def1(tree)
    write(io, '(')
    if def2(tree)
      _ibuffer(io, tree.d1), 
      write(io, ',')
      _ibuffer(io, tree.d2), 
      print(io, ',', e(tree), ',', dt(tree), ',', fdt(tree), ',',
            short(isextinct(tree)), ',', "0,", short(isfix(tree)), ',', 
            lλ(tree), ',', lμ(tree), ')')
    else
      _ibuffer(io, tree.d1), 
      print(io, ',', e(tree), ',', dt(tree), ',', fdt(tree), ',',
            short(isextinct(tree)), ',', "1,", short(isfix(tree)), ',', 
            lλ(tree), ',', lμ(tree), ')')
    end
  else
    print(io, '(', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          short(isextinct(tree)), ',', short(isfossil(tree)), ',', 
          short(isfix(tree)), ',', lλ(tree), ',', lμ(tree), ')')
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
  # getfield(Tapestree.INSANE, Symbol(st))
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
    _iparse(s::String, i::Int64, ls::Int64, ::Type{sTxs})

parse istring to `iT`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{sTxs})

  @inbounds begin

    in1 = false
    in2 = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, sTxs)
      in1 = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, sTxs)
      in2 = true
    end

    i1 = findnext(',', s, i  + 1)
    i2 = findnext(',', s, i1 + 1)
    i3 = findnext(',', s, i2 + 1)
    i4 = findnext(']', s, i3 + 1)
    i5 = findnext(']', s, i4 + 1)

    if in1
      if in2
        tree = sTxs(sd1, sd2,
                    Pparse(Float64, s[i:i1-1]),
                    Pparse(Float64, s[i1+1:i2-1]),
                    Pparse(Float64, s[i2+1:i3-1]),
                    _iparse_v(s[i3+2:i4-1]),
                    _iparse_v(s[i4+3:i5-1]))
      else
        tree = sTxs(sd1,
                    Pparse(Float64, s[i:i1-1]),
                    Pparse(Float64, s[i1+1:i2-1]),
                    Pparse(Float64, s[i2+1:i3-1]),
                    _iparse_v(s[i3+2:i4-1]),
                    _iparse_v(s[i4+3:i5-1]))
      end
    else
      tree = sTxs(Pparse(Float64, s[i:i1-1]),
                  Pparse(Float64, s[i1+1:i2-1]),
                  Pparse(Float64, s[i2+1:i3-1]),
                  _iparse_v(s[i3+2:i4-1]),
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
_iparse_v(s::String) = Pparse.(Float64, split(s, ','))




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



