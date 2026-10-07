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

  io = open(in_file, "r")

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
            push!(tv, 
              _parse_newick(SubString(line, (allsc[i] + 1),(allsc[i+1]))))
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
            push!(tv, 
              _parse_newick(SubString(line, (allsc[i] + 1), (allsc[i+1])), ne))
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
    _parse_newick(in_file::AbstractString)

Reads a newick tree into `sT` from `in_file`.
"""
function _parse_newick(s::AbstractString) 
  tree, i = _from_string(s, 1, sT_label)
  return tree
end




"""
    _parse_newick(s::AbstractString, ne::Float64)

Reads a newick tree into `sT` if fossil is false and `sTf` if fossil
is true from `in_file`.
"""
function _parse_newick(s::AbstractString, ne::Float64)
  tree, i = _from_string(s, 1, sTf_label)
  fossilizepasttips!(tree, ne)

  return tree
end




"""
    _from_string(s::AbstractString, ::Type{T}) where {T <: sT}

Returns a tree of type `T` from newick string.
"""
function _from_string(s::AbstractString, i::Int64, ::Type{T}) where {T <: sT}

  @inbounds begin

    in1 = in2 = false

    if s[i] === '('
      sd1, i = _from_string(s, i + 1, T)
      in1 = true
    end

    if s[i] === ','
      sd2, i = _from_string(s, i + 1, T)
      in2 = true
    end

    i1 = findnext(':', s, i)

    # if root
    if isnothing(i1)
      # if stem tree
      if !in2
        tree = sd1
      # if fossil 1
      elseif e(sd1) === 0.0
        tree = T(sd2, 0.0, label(sd1))
      # if fossil 2
      elseif e(sd2) === 0.0
        tree = T(sd1, 0.0, label(sd2))
      # if crown
      else
        tree = T(sd1, sd2, 0.0, "")
      end
      return tree, i
    end

    i2 = find_cp(s, i1 + 1)

    if in1
      if in2
        if e(sd1) === 0.0
          tree = T(sd2, parse(Float64, SubString(s, i1+1, i2-1)), label(sd1))
        elseif e(sd2) === 0.0
          tree = T(sd1, parse(Float64, SubString(s, i1+1, i2-1)), label(sd2))
        else
          tree = T(sd1, sd2, 
                   parse(Float64, SubString(s, i1+1, i2-1)), 
                   SubString(s, i+1, i1-1))
        end
      else
         tree = T(sd1, 
                  parse(Float64, SubString(s, i1+1, i2-1)), 
                  SubString(s, i+1, i1-1))
      end
    else
      tree = T(parse(Float64, SubString(s, i1+1, i2-1)), SubString(s, i, i1-1))
    end

    i = i2
  end

  return tree, i
end




"""
    find_cp(s::AbstractString, i::Int64)

Find next ',', ')' or ';' after index `i`.
"""
function find_cp(s::AbstractString, i::Int64)

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
    onlyone(s::String, c::Char)

Returns true if there is only one of 'c' in string `s`.
"""
function onlyone(s::AbstractString, c::Char)
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
    _to_string(tree::T, n::Int64, nf::Int64) where {T <: iTf}

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
function write_nexus(treev  ::Vector{T}, 
                     reftree::Tl, 
                     ofile  ::String) where {T <: iTree, Tl <: Tlabel}

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
    nx_buffer(io::IOBuffer, tree::T, reftree::Tl, ic::Bool) where {T <: Tpe, Tl <: Tlabel}

Writes an `iTree` to IOBuffer `io`.
"""
nx_buffer(io::IOBuffer, tree::T, reftree::Tl, ic::Bool) where {T <: Tpe, Tl <: Tlabel} = 
  _nx_buffer(io, tree, reftree, 0, 0, ic)

"""
    _nx_buffer(io     ::IOBuffer, 
               tree   ::T, 
               reftree::Tl, 
               n      ::Int64, 
               nf     ::Int64, 
               ic     ::Bool) where {T <: Tpe, Tl <: Tlabel}

Writes an `iTree` to IOBuffer `io`.
"""
function _nx_buffer(io     ::IOBuffer, 
                    tree   ::T, 
                    reftree::Tl, 
                    n      ::Int64, 
                    nf     ::Int64, 
                    ic     ::Bool) where {T <: Tpe, Tl <: Tlabel}

  if def1(tree)
    write(io, '(')
    if  def2(tree)
      if isfix(tree.d1) && isfix(tree.d2)
        n, nf = _nx_buffer(io, tree.d1, reftree.d1, n, nf, false)
        write(io, ',')
        n, nf = _nx_buffer(io, tree.d2, reftree.d2, n, nf, false)
      else
        n, nf = _nx_buffer(io, tree.d1, reftree, n, nf, false)
        write(io, ',')
        n, nf = _nx_buffer(io, tree.d2, reftree, n, nf, false)
      end
      print(io, ")[&xi=", xi(tree), ",xf=", xf(tree), ",sh1=", sh(tree), 
                ",isextinct=false,isfossil=false", ",da=", !isfix(tree), ']')
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
      print(io, "[&xi=", xi(tree), ",xf=", xf(tree), ",sh1=", sh(tree), 
                ",isextinct=false,isfossil=true", 
                ",da=", !isfix(tree), "]:", e(tree))
    end
  else
    if isfix(tree)
      print(io, label(reftree))
    else
      n += 1
      print(io, 't', n)
    end
    print(io, "[&xi=", xi(tree), ",xf=", xf(tree), 
              ",sh1=", sh(tree), ",isextinct=", isextinct(tree), 
              ",isfossil=", isfossil(tree), ",da=", !isfix(tree), "]:", e(tree))
  end

  return n, nf
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
               nf     ::Int64, 
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
    nx_buffer(io::IOBuffer, tree::iTfbd, reftree::sTf_label, ic::Bool)

Writes an `iTree` to IOBuffer `io`.
"""
nx_buffer(io::IOBuffer, tree::sTxs, reftree::T, ic::Bool) where {T <: Tlabel} = 
  _nx_buffer(io, tree, reftree, ic)

"""
    _nx_buffer(io     ::IOBuffer, 
               tree   ::sTxs, 
               reftree::T, 
               ic     ::Bool) where {T <: Tlabel}

Writes an `iTree` to IOBuffer `io`.
"""
function _nx_buffer(io     ::IOBuffer, 
                    tree   ::sTxs, 
                    reftree::T, 
                    ic     ::Bool) where {T <: Tlabel}

  if def1(tree)
    write(io, '(')
    if def2(tree)
      _nx_buffer(io, tree.d1, reftree.d1, false)
      write(io, ',')
      _nx_buffer(io, tree.d2, reftree.d2, false)
      write(io, ")[&x=")
      nx_printv(io, xv(tree))
      write(io, ",logs=")
      nx_printv(io, lσ2(tree))
      print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), ']')
      !ic && print(io, ':', e(tree))
    else
      _nx_buffer(io, tree.d1, reftree.d1, false)
      write(io, ')')
      print(io, label(reftree))
      write(io, "[&x=")
      nx_printv(io, xv(tree))
      write(io, ",logs=")
      nx_printv(io, lσ2(tree))
      print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), "]:", e(tree))
    end
  else
    print(io, label(reftree))
    write(io, "[&x=")
    nx_printv(io, xv(tree))
    write(io, ",logs=")
    nx_printv(io, lσ2(tree))
    print(io, ",dt=", dt(tree), ",fdt=", fdt(tree), "]:", e(tree))
  end

  return nothing
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
    bp  = parse(Float64, SubString(x, 1, (pix-1)))
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
    write(tree::T, ofile::String) where {T <: iTree}

Write iTree to file.
"""
function Base.write(tree::T, ofile::String) where {T <: iTree}
  open(ofile*".insane", "w") do io
    typestring = TtoS[T]
    write(io, lastindex(typestring))
    write(io, typestring)
    write(io, tree)
  end
  return nothing
end


"""
    write(trees::Vector{T}, ofile::String) where {T <: iTree}

Write a vector of trees.
"""
function Base.write(trees::Vector{T}, ofile::String) where {T <: iTree}
  open(ofile*".insane", "w") do io
    typestring = TtoS[T]
    write(io, lastindex(typestring))
    write(io, typestring)
    for tree in trees
      write(io, tree)
    end
  end
end


"""
    read(file::String; ix::OrdinalRange{Int64,Int64} = 0:0)

Read an `.insane` file with optional OrdinalRange specifying which 
data augmented trees to sample.
"""
function Base.read(file::String)

  open(file, "r") do io
    typelength = read(io, Int64)
    T          = StoT[String(read(io, typelength))]
    tree       = read(io, T)

    if eof(io)
      return tree
    else
      trees = T[tree]
      while !eof(io)
        try
          push!(trees, read(io, T))
        catch e
          if e isa EOFError
            break
          else
            rethrow(e)
          end
        end
      end
    end
    return trees
  end
end




"""
    write(io::IOBuffer, tree::sTb)

Write `sTb` to IO.
"""
function Base.write(io::IO, tree::sTb)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, isfix(tree))
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isfix(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{sTb})

Read `sTb` from IO.
"""
function Base.read(io::IO, ::Type{sTb})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    return sTb(ei, ifx)
  else
    d1  = read(io, sTb)
    d2  = read(io, sTb)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    return sTb(d1, d2, ei, ifx)
  end
end




"""
    write(io::IOBuffer, tree::sTbd)

Write `sTbd` to IO.
"""
function Base.write(io::IO, tree::sTbd)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, isfix(tree))
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, isfix(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{sTbd})

Read `sTbd` from IO.
"""
function Base.read(io::IO, ::Type{sTbd})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    ifx = read(io, Bool)
    return sTbd(ei, iμ, ifx)
  else
    d1  = read(io, sTbd)
    d2  = read(io, sTbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    return sTbd(d1, d2, ei, false, ifx)
  end
end




"""
    write(io::IOBuffer, tree::sTpe)

Write `sTpe` to IO.
"""
function Base.write(io::IO, tree::sTpe)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, xi(tree))
    write(io, xf(tree))
    write(io, sh(tree))
    write(io, isfix(tree))
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, xi(tree))
    write(io, xf(tree))
    write(io, isfix(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{sTpe})

Read `sTpe` from IO.
"""
function Base.read(io::IO, ::Type{sTpe})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    xi  = read(io, Bool)
    xf  = read(io, Bool)
    ifx = read(io, Bool)
    return sTpe(ei, iμ, xi, xf, false, ifx)
  else
    d1  = read(io, sTpe)
    d2  = read(io, sTpe)
    ei  = read(io, Float64)
    xi  = read(io, Bool)
    xf  = read(io, Bool)
    sh  = read(io, Bool)
    ifx = read(io, Bool)
    return sTpe(d1, d2, ei, false, xi, xf, sh, ifx)
  end
end




"""
    write(io::IOBuffer, tree::sTfbd)

Write `sTfbd` to IO.
"""
function Base.write(io::IO, tree::sTfbd)
  if def1(tree)
    if def2(tree)
      write(io, 2)
      write(io, tree.d1)
      write(io, tree.d2)
      write(io, e(tree))
      write(io, isfix(tree))
    else
      write(io, 1)
      write(io, tree.d1)
      write(io, e(tree))
      write(io, isfix(tree))
    end
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, isfossil(tree))
    write(io, isfix(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{sTfbd})

Read `sTfbd` from IO.
"""
function Base.read(io::IO, ::Type{sTfbd})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    iψ  = read(io, Bool)
    ifx = read(io, Bool)
    return sTfbd(ei, iμ, iψ, ifx)
  elseif isone(nds)
    d1  = read(io, sTfbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    return sTfbd(d1, ei, false, true, ifx)
  else
    d1  = read(io, sTfbd)
    d2  = read(io, sTfbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    return sTfbd(d1, d2, ei, false, false, ifx)
  end
end




"""
    write(io::IOBuffer, tree::sTfpe)

Write `sTfpe` to IO.
"""
function Base.write(io::IO, tree::sTfpe)
  if def1(tree)
    if def2(tree)
      write(io, 2)
      write(io, tree.d1)
      write(io, tree.d2)
      write(io, e(tree))
      write(io, xi(tree))
      write(io, xf(tree))
      write(io, sh(tree))
      write(io, isfix(tree))
    else
      write(io, 1)
      write(io, tree.d1)
      write(io, e(tree))
      write(io, xi(tree))
      write(io, xf(tree))
      write(io, isfix(tree))
    end
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, isfossil(tree))
    write(io, xi(tree))
    write(io, xf(tree))
    write(io, isfix(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{sTfpe})

Read `sTfpe` from IO.
"""
function Base.read(io::IO, ::Type{sTfpe})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    iψ  = read(io, Bool)
    xi  = read(io, Bool)
    xf  = read(io, Bool)
    ifx = read(io, Bool)
    return sTfpe(ei, iμ, iψ, xi, xf, false, ifx)
  elseif isone(nds)
    d1  = read(io, sTfpe)
    ei  = read(io, Float64)
    xi  = read(io, Bool)
    xf  = read(io, Bool)
    ifx = read(io, Bool)
    return sTfpe(d1, ei, false, true, xi, xf, false, ifx)
  else
    d1  = read(io, sTfpe)
    d2  = read(io, sTfpe)
    ei  = read(io, Float64)
    xi  = read(io, Bool)
    xf  = read(io, Bool)
    sh  = read(io, Bool)
    ifx = read(io, Bool)
    return sTfpe(d1, d2, ei, false, false, xi, xf, sh, ifx)
  end
end




"""
    write(io::IOBuffer, tree::cTb)

Write `cTb` to IO.
"""
function Base.write(io::IO, tree::cTb)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{cTb})

Read `cTb` from IO.
"""
function Base.read(io::IO, ::Type{cTb})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    return cTb(ei, ifx, lλ)
  else
    d1  = read(io, cTb)
    d2  = read(io, cTb)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    return cTb(d1, d2, ei, ifx, lλ)
  end
end




"""
    write(io::IOBuffer, tree::cTce)

Write `cTce` to IO.
"""
function Base.write(io::IO, tree::cTce)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{cTce})

Read `cTce` from IO.
"""
function Base.read(io::IO, ::Type{cTce})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    return cTce(ei, iμ, ifx, lλ)
  else
    d1  = read(io, cTce)
    d2  = read(io, cTce)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    return cTce(d1, d2, ei, false, ifx, lλ)
  end
end




"""
    write(io::IOBuffer, tree::cTct)

Write `cTct` to IO.
"""
function Base.write(io::IO, tree::cTct)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{cTct})

Read `cTct` from IO.
"""
function Base.read(io::IO, ::Type{cTct})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    return cTct(ei, iμ, ifx, lλ)
  else
    d1  = read(io, cTct)
    d2  = read(io, cTct)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    return cTct(d1, d2, ei, false, ifx, lλ)
  end
end




"""
    write(io::IOBuffer, tree::iTb)

Write `iTb` to IO.
"""
function Base.write(io::IO, tree::iTb)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    write(io, lastindex(lλv))
    write(io, lλv)
  else
    write(io, 0)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    write(io, lastindex(lλv))
    write(io, lλv)
  end
end

"""
    read(io::IOBuffer, ::Type{iTb})

Read `iTb` from IO.
"""
function Base.read(io::IO, ::Type{iTb})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    ifx = read(io, Bool)
    lλv = Vector{Float64}(undef, read(io, Int64))
    read!(io, lλv)
    return iTb(ei, dt, fdt, ifx, lλv)
  else
    d1  = read(io, iTb)
    d2  = read(io, iTb)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    ifx = read(io, Bool)
    lλv = Vector{Float64}(undef, read(io, Int64))
    read!(io, lλv)
    return iTb(d1, d2, ei, dt, fdt, ifx, lλv)
  end
end




"""
    write(io::IOBuffer, tree::iTce)

Write `iTce` to IO.
"""
function Base.write(io::IO, tree::iTce)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    write(io, lastindex(lλv))
    write(io, lλv)
  else
    write(io, 0)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isextinct(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    write(io, lastindex(lλv))
    write(io, lλv)
  end
end

"""
    read(io::IOBuffer, ::Type{iTce})

Read `iTce` from IO.
"""
function Base.read(io::IO, ::Type{iTce})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    iμ  = read(io, Bool)
    ifx = read(io, Bool)
    lλv = Vector{Float64}(undef, read(io, Int64))
    read!(io, lλv)
    return iTce(ei, dt, fdt, iμ, ifx, lλv)
  else
    d1  = read(io, iTce)
    d2  = read(io, iTce)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    ifx = read(io, Bool)
    lλv = Vector{Float64}(undef, read(io, Int64))
    read!(io, lλv)
    return iTce(d1, d2, ei, dt, fdt, false, ifx, lλv)
  end
end




"""
    write(io::IOBuffer, tree::iTct)

Write `iTct` to IO.
"""
function Base.write(io::IO, tree::iTct)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    write(io, lastindex(lλv))
    write(io, lλv)
  else
    write(io, 0)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isextinct(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    write(io, lastindex(lλv))
    write(io, lλv)
  end
end

"""
    read(io::IOBuffer, ::Type{iTct})

Read `iTct` from IO.
"""
function Base.read(io::IO, ::Type{iTct})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    iμ  = read(io, Bool)
    ifx = read(io, Bool)
    lλv = Vector{Float64}(undef, read(io, Int64))
    read!(io, lλv)
    return iTct(ei, dt, fdt, iμ, ifx, lλv)
  else
    d1  = read(io, iTct)
    d2  = read(io, iTct)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    ifx = read(io, Bool)
    lλv = Vector{Float64}(undef, read(io, Int64))
    read!(io, lλv)
    return iTct(d1, d2, ei, dt, fdt, false, ifx, lλv)
  end
end




"""
    write(io::IOBuffer, tree::cTbd)

Write `cTbd` to IO.
"""
function Base.write(io::IO, tree::cTbd)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
    write(io, lμ(tree))
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
    write(io, lμ(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{cTbd})

Read `cTbd` from IO.
"""
function Base.read(io::IO, ::Type{cTbd})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    lμ  = read(io, Float64)
    return cTbd(ei, iμ, ifx, lλ, lμ)
  else
    d1  = read(io, cTbd)
    d2  = read(io, cTbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    lμ  = read(io, Float64)
    return cTbd(d1, d2, ei, false, ifx, lλ, lμ)
  end
end




"""
    write(io::IOBuffer, tree::iTbd)

Write `iTbd` to IO.
"""
function Base.write(io::IO, tree::iTbd)
  if def1(tree)
    write(io, 2)
    write(io, tree.d1)
    write(io, tree.d2)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    n   = lastindex(lλv)
    lμv = lμ(tree)
    write(io, n)
    write(io, lλv)
    write(io, lμv)
  else
    write(io, 0)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isextinct(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    n   = lastindex(lλv)
    lμv = lμ(tree)
    write(io, n)
    write(io, lλv)
    write(io, lμv)
  end
end

"""
    read(io::IOBuffer, ::Type{iTbd})

Read `iTbd` from IOBuffer.
"""
function Base.read(io::IO, ::Type{iTbd})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    iμ  = read(io, Bool)
    ifx = read(io, Bool)
    n   = read(io, Int64)
    lλv = Vector{Float64}(undef, n)
    read!(io, lλv)
    lμv = Vector{Float64}(undef, n)
    read!(io, lμv)
    return iTbd(ei, dt, fdt, iμ, ifx, lλv, lμv)
  else
    d1  = read(io, iTbd)
    d2  = read(io, iTbd)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    ifx = read(io, Bool)
    n   = read(io, Int64)
    lλv = Vector{Float64}(undef, n)
    read!(io, lλv)
    lμv = Vector{Float64}(undef, n)
    read!(io, lμv)
    return iTbd(d1, d2, ei, dt, fdt, false, ifx, lλv, lμv)
  end
end




"""
    write(io::IOBuffer, tree::cTfbd)

Write `cTfbd` to IO.
"""
function Base.write(io::IO, tree::cTfbd)
  if def1(tree)
    if def2(tree)
      write(io, 2)
      write(io, tree.d1)
      write(io, tree.d2)
      write(io, e(tree))
      write(io, isfix(tree))
      write(io, lλ(tree))
      write(io, lμ(tree))
    else
      write(io, 1)
      write(io, tree.d1)
      write(io, e(tree))
      write(io, isfix(tree))
      write(io, lλ(tree))
      write(io, lμ(tree))
    end
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, isfossil(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
    write(io, lμ(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{cTfbd})

Read `cTfbd` from IOBuffer.
"""
function Base.read(io::IO, ::Type{cTfbd})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    iψ  = read(io, Bool)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    lμ  = read(io, Float64)
    return cTfbd(ei, iμ, iψ, ifx, lλ, lμ)
  elseif isone(nds)
    d1  = read(io, cTfbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    lμ  = read(io, Float64)
    return cTfbd(d1, ei, false, true, ifx, lλ, lμ)
  else
    d1  = read(io, cTfbd)
    d2  = read(io, cTfbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    lμ  = read(io, Float64)
    return cTfbd(d1, d2, ei, false, false, ifx, lλ, lμ)
  end
end




"""
    write(io::IOBuffer, tree::acTfbd)

Write `acTfbd` to IO.
"""
function Base.write(io::IO, tree::acTfbd)
  if def1(tree)
    if def2(tree)
      write(io, 2)
      write(io, tree.d1)
      write(io, tree.d2)
      write(io, e(tree))
      write(io, isfix(tree))
      write(io, sh(tree))
      write(io, lλ(tree))
      write(io, lμ(tree))
    else
      write(io, 1)
      write(io, tree.d1)
      write(io, e(tree))
      write(io, isfix(tree))
      write(io, lλ(tree))
      write(io, lμ(tree))
    end
  else
    write(io, 0)
    write(io, e(tree))
    write(io, isextinct(tree))
    write(io, isfossil(tree))
    write(io, isfix(tree))
    write(io, lλ(tree))
    write(io, lμ(tree))
  end
end

"""
    read(io::IOBuffer, ::Type{acTfbd})

Read `acTfbd` from IOBuffer.
"""
function Base.read(io::IO, ::Type{acTfbd})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    iμ  = read(io, Bool)
    iψ  = read(io, Bool)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    lμ  = read(io, Float64)
    return acTfbd(ei, iμ, iψ, false, ifx, lλ, lμ)
  elseif isone(nds)
    d1  = read(io, acTfbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    lμ  = read(io, Float64)
    return acTfbd(d1, ei, false, true, false, ifx, lλ, lμ)
  else
    d1  = read(io, acTfbd)
    d2  = read(io, acTfbd)
    ei  = read(io, Float64)
    sh  = read(io, Bool)
    ifx = read(io, Bool)
    lλ  = read(io, Float64)
    lμ  = read(io, Float64)
    return acTfbd(d1, d2, ei, false, false, sh, ifx, lλ, lμ)
  end
end





"""
    write(io::IOBuffer, tree::iTfbd)

Write `iTfbd` to IO.
"""
function Base.write(io::IO, tree::iTfbd)
  if def1(tree)
    if def2(tree)
      write(io, 2)
      write(io, tree.d1)
      write(io, tree.d2)
      write(io, e(tree))
      write(io, dt(tree))
      write(io, fdt(tree))
      write(io, isfix(tree))
      lλv = lλ(tree)
      n   = lastindex(lλv)
      lμv = lμ(tree)
      write(io, n)
      write(io, lλv)
      write(io, lμv)
    else
      write(io, 1)
      write(io, tree.d1)
      write(io, e(tree))
      write(io, dt(tree))
      write(io, fdt(tree))
      write(io, isfix(tree))
      lλv = lλ(tree)
      n   = lastindex(lλv)
      lμv = lμ(tree)
      write(io, n)
      write(io, lλv)
      write(io, lμv)
    end
  else
    write(io, 0)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    write(io, isextinct(tree))
    write(io, isfossil(tree))
    write(io, isfix(tree))
    lλv = lλ(tree)
    n   = lastindex(lλv)
    lμv = lμ(tree)
    write(io, n)
    write(io, lλv)
    write(io, lμv)
  end
end

"""
    read(io::IOBuffer, ::Type{iTfbd})

Read `iTfbd` from IOBuffer.
"""
function Base.read(io::IO, ::Type{iTfbd})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    iμ  = read(io, Bool)
    iψ  = read(io, Bool)
    ifx = read(io, Bool)
    n   = read(io, Int64)
    lλv = Vector{Float64}(undef, n)
    read!(io, lλv)
    lμv = Vector{Float64}(undef, n)
    read!(io, lμv)
    return iTfbd(ei, dt, fdt, iμ, iψ, ifx, lλv, lμv)
  elseif isone(nds)
    d1  = read(io, iTfbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    n   = read(io, Int64)
    lλv = Vector{Float64}(undef, n)
    read!(io, lλv)
    lμv = Vector{Float64}(undef, n)
    read!(io, lμv)
    return iTfbd(d1, ei, dt, fdt, false, true, ifx, lλv, lμv)
  else
    d1  = read(io, iTfbd)
    d2  = read(io, iTfbd)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    n   = read(io, Int64)
    lλv = Vector{Float64}(undef, n)
    read!(io, lλv)
    lμv = Vector{Float64}(undef, n)
    read!(io, lμv)
    return iTfbd(d1, d2, ei, dt, fdt, false, false, ifx, lλv, lμv)
  end
end




"""
    write(io::IOBuffer, tree::sTxs)

Write `sTxs` to IO.
"""
function Base.write(io::IO, tree::sTxs)
  if def1(tree)
    if def2(tree)
      write(io, 2)
      write(io, tree.d1)
      write(io, tree.d2)
      write(io, e(tree))
      write(io, dt(tree))
      write(io, fdt(tree))
      xv   = xv(tree)
      n    = lastindex(xv)
      lσ2v = lσ2(tree)
      write(io, n)
      write(io, xv)
      write(io, lσ2v)
    else
      write(io, 1)
      write(io, tree.d1)
      write(io, e(tree))
      write(io, dt(tree))
      write(io, fdt(tree))
      xv   = xv(tree)
      n    = lastindex(xv)
      lσ2v = lσ2(tree)
      write(io, n)
      write(io, xv)
      write(io, lσ2v)
    end
  else
    write(io, 0)
    write(io, e(tree))
    write(io, dt(tree))
    write(io, fdt(tree))
    xv = xv(tree)
    n   = lastindex(xv)
    lσ2v = lσ2(tree)
    write(io, n)
    write(io, xv)
    write(io, lσ2v)
  end
end

"""
    read(io::IOBuffer, ::Type{sTxs})

Read `sTxs` from IOBuffer.
"""
function Base.read(io::IO, ::Type{sTxs})
  nds = read(io, Int64)
  if iszero(nds)
    ei  = read(io, Float64)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    n   = read(io, Int64)
    xv = Vector{Float64}(undef, n)
    read!(io, xv)
    lσ2v = Vector{Float64}(undef, n)
    read!(io, lσ2v)
    return sTxs(ei, dt, fdt, iμ, iψ, ifx, xv, lσ2v)
  elseif isone(nds)
    d1  = read(io, sTxs)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    n   = read(io, Int64)
    xv = Vector{Float64}(undef, n)
    read!(io, xv)
    lσ2v = Vector{Float64}(undef, n)
    read!(io, lσ2v)
    return sTxs(d1, ei, dt, fdt, false, true, ifx, xv, lσ2v)
  else
    d1  = read(io, sTxs)
    d2  = read(io, sTxs)
    ei  = read(io, Float64)
    ifx = read(io, Bool)
    dt  = read(io, Float64)
    fdt = read(io, Float64)
    n   = read(io, Int64)
    xv = Vector{Float64}(undef, n)
    read!(io, xv)
    lσ2v = Vector{Float64}(undef, n)
    read!(io, lσ2v)
    return sTxs(d1, d2, ei, dt, fdt, false, false, ifx, xv, lσ2v)
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
    _ibuffer(io::IOBuffer, tree::sTb)

Write `sTb` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::sTb)
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
    _ibuffer(io::IOBuffer, tree::sTpe)

Write `sTpe` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::sTpe)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1)
    write(io, ',')
    _ibuffer(io, tree.d2)
    print(io, ',', e(tree), ',', 
          short(isextinct(tree)), ',', 
          xi(tree), ',',
          xf(tree), ',',
          short(sh(tree)), ',',
          short(isfix(tree)), ')')
  else
    print(io, '(', e(tree), ',', 
        short(isextinct(tree)), ',',
        xi(tree), ',',
        xf(tree), ',',
        short(sh(tree)), ',',
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
    _ibuffer(io::IOBuffer, tree::sTfpe)

Write `sTfpe` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::sTfpe)
  if def1(tree)
    write(io, '(')
    if def2(tree)
      _ibuffer(io, tree.d1), 
      write(io, ',') 
      _ibuffer(io, tree.d2), 
      print(io, ',', e(tree), ',', 
          short(isextinct(tree)), ',', "0,", 
          xi(tree), ',',
          xf(tree), ',',
          short(sh(tree)), ',',
          short(isfix(tree)), ')')
    else
      _ibuffer(io, tree.d1), 
      print(io, ',', e(tree), ',', 
          short(isextinct(tree)), ',', "1,", 
          xi(tree), ',',
          xf(tree), ',',
          short(sh(tree)), ',',
          short(isfix(tree)), ')')
    end
  else
    print(io, '(', e(tree), ',', 
        short(isextinct(tree)), ',', 
        short(isfossil(tree)), ',',
        xi(tree), ',',
        xf(tree), ',',
        short(sh(tree)), ',',
        short(isfix(tree)), ')')
  end
end



"""
    _ibuffer(io::IOBuffer, tree::sTxs)

`sTxs` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::sTxs)
  if def1(tree)
    write(io, '(')
    if def2(tree)
      _ibuffer(io, tree.d1)
      write(io, ',')
      _ibuffer(io, tree.d2)
      print(io, ',', e(tree), ',', dt(tree), ',', fdt(tree), ',', 
                xv(tree), ',', lσ2(tree), ')')
    else
      _ibuffer(io, tree.d1)
      print(io, ',', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          xv(tree), ',', lσ2(tree), ')')
    end
  else
    print(io, '(', e(tree), ',', dt(tree), ',', fdt(tree), ',', 
              xv(tree), ',', lσ2(tree), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::cTb)

Write `cTb` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::cTb)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1), 
    write(io, ',')
    _ibuffer(io, tree.d2), 
    print(io, ',', e(tree), ',', short(isfix(tree)), ',', lλ(tree), ')')
  else
    print(io, '(', e(tree), ',', short(isfix(tree)), ',', lλ(tree), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::T) where {T <: cT}

Write `cTce` or `cTct` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::T) where {T <: cT}
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1), 
    write(io, ',')
    _ibuffer(io, tree.d2), 
    print(io, ',', e(tree), ',', short(isextinct(tree)), ',',
      short(isfix(tree)), ',', lλ(tree), ')')
  else
    print(io, '(', e(tree), ',', short(isextinct(tree)), ',', 
      short(isfix(tree)), ',', lλ(tree), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::cTbd)

Write `cTbd` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::cTbd)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1), 
    write(io, ',')
    _ibuffer(io, tree.d2), 
    print(io, ',', e(tree), ',', short(isextinct(tree)), ',',
      short(isfix(tree)), ',', lλ(tree), ',', lμ(tree), ')')
  else
    print(io, '(', e(tree), ',', short(isextinct(tree)), ',', 
      short(isfix(tree)), ',', lλ(tree), ',', lμ(tree), ')')
  end
end




"""
    _ibuffer(io::IOBuffer, tree::cTfbd)

Write `cTfbd` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::cTfbd)
  if def1(tree)
    write(io, '(')
    if def2(tree)
      _ibuffer(io, tree.d1), 
      write(io, ',')
      _ibuffer(io, tree.d2), 
      print(io, ',', e(tree), ',',
            short(isextinct(tree)), ',', "0,", short(isfix(tree)), ',', 
            lλ(tree), ',', lμ(tree), ')')
    else
      _ibuffer(io, tree.d1), 
      print(io, ',', e(tree), ',',
            short(isextinct(tree)), ',', "1,", short(isfix(tree)), ',', 
            lλ(tree), ',', lμ(tree), ')')
    end
  else
    print(io, '(', e(tree), ',',
          short(isextinct(tree)), ',', short(isfossil(tree)), ',', 
          short(isfix(tree)), ',', lλ(tree), ',', lμ(tree), ')')
  end
end





"""
    _ibuffer(io::IOBuffer, tree::iTb)

Write `iTb` to IOBuffer.
"""
function _ibuffer(io::IOBuffer, tree::iTb)
  if def1(tree)
    write(io, '(')
    _ibuffer(io, tree.d1), 
    write(io, ',')
    _ibuffer(io, tree.d2), 
    print(io, ',', e(tree), ',', dt(tree), ',', fdt(tree), ',',
          short(isextinct(tree)), ',', short(isfix(tree)), ',', lλ(tree), ')')
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
    _istring(tree::iTpbd)
`iTpbd` to istring.
"""
function _istring(tree::iTpbd)
  if def1(tree)
    if def2(tree)
      return string('(', _istring(tree.d1), ',', _istring(tree.d2), ',', 
               e(tree), ',',
               dt(tree), ',',
               fdt(tree), ',',
               short(isextinct(tree)), ',', 
               short(iscomplete(tree)), ',', 
               short(isfix(tree)), ',', 
               lb(tree), ',', 
               lλ(tree), ',', 
               lμ(tree), ')')
    else
      return string('(', _istring(tree.d1), ',', 
               e(tree), ',',
               dt(tree), ',',
               fdt(tree), ',',
               short(isextinct(tree)), ',', 
               short(iscomplete(tree)), ',', 
               short(isfix(tree)), ',', 
               lb(tree), ',', 
               lλ(tree), ',', 
               lμ(tree), ')')
    end
  else
    return string('(', 
             e(tree), ',',
             dt(tree), ',',
             fdt(tree), ',',
             short(isextinct(tree)), ',', 
             short(iscomplete(tree)), ',', 
             short(isfix(tree)), ',', 
             lb(tree), ',', 
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
  st = SubString(s, 1, (i-1))
  T  = iTd[st]
  # getfield(Tapestree.INSANE, Symbol(st))
  si = s[i+2:ls-1]

  t0, ix = _iparse(si, 1, ls - lastindex(st) - 3, T)

  return t0 
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{sTb})
parse istring to `sTb`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{sTb})

 @inbounds begin

    inode = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, sTb)
      inode = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, sTb)
    end

    i1 = findnext(',', s, i + 1)

    if inode
      tree = sTb(sd1, sd2, 
                 parse(Float64, SubString(s, i, i1-1)), 
                 long(s[i1+1]))
    else
      tree = sTb(parse(Float64, SubString(s, i, i1-1)), 
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
                  parse(Float64, SubString(s, i, i1-1)), 
                  long(s[i1+1]), long(s[i1+3]))
    else
      tree = sTbd(parse(Float64, SubString(s, i, i1-1)), 
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
                     parse(Float64, SubString(s, i, i1-1)), 
                     long(s[i1+1]), long(s[i1+3]), long(s[i1+5]))
      else
        tree = sTfbd(sd1, parse(Float64, SubString(s, i, i1-1)), 
                     long(s[i1+1]), long(s[i1+3]), long(s[i1+5]))
      end
    else
      tree = sTfbd(parse(Float64, SubString(s, i, i1-1)), 
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
    _iparse(s::String, i::Int64, ls::Int64, ::Type{sTpe})

parse istring to `sTpe`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{sTpe})

  @inbounds begin

    inode = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, sTpe)
      inode = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, sTpe)
    end

    i1 = findnext(',', s, i + 1)
    i2 = findnext(',', s, i1 + 3)
    i3 = findnext(',', s, i2 + 1)

    if inode
      tree = sTpe(sd1, sd2, 
                  parse(Float64, SubString(s, i, i1-1)), long(s[i1+1]), 
                  parse(Float64, SubString(s, i1+3, i2-1)), 
                  parse(Float64, SubString(s, i2+1, i3-1)),
                  long(s[i3+1]), long(s[i3+3]))
    else
      tree = sTpe(parse(Float64, SubString(s, i, i1-1)), long(s[i1+1]), 
                  parse(Float64, SubString(s, i1+3, i2-1)), 
                  parse(Float64, SubString(s, i2+1, i3-1)),
                  long(s[i3+1]), long(s[i3+3]))
    end

    i = i3 + 4

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{sTfpe})

parse istring to `sTfpe`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{sTfpe})

  @inbounds begin

    in1 = false
    in2 = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, sTfpe)
      in1 = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, sTfpe)
      in2 = true
    end

    i1 = findnext(',', s, i + 1)
    i2 = findnext(',', s, i1 + 5)
    i3 = findnext(',', s, i2 + 1)

    if in1
      if in2
      tree = sTfpe(sd1, sd2, 
                   parse(Float64, SubString(s, i, i1-1)), 
                   long(s[i1+1]), long(s[i1+3]), 
                   parse(Float64, SubString(s, i1+5, i2-1)), 
                   parse(Float64, SubString(s, i2+1, i3-1)),
                   long(s[i3+1]), long(s[i3+3]))
      else
        tree = sTfpe(sd1,
                     parse(Float64, SubString(s, i, i1-1)), 
                     long(s[i1+1]), long(s[i1+3]), 
                     parse(Float64, SubString(s, i1+5, i2-1)), 
                     parse(Float64, SubString(s, i2+1, i3-1)),
                     long(s[i3+1]), long(s[i3+3]))
      end
    else
      tree = sTfpe(parse(Float64, SubString(s, i, i1-1)), 
                   long(s[i1+1]), long(s[i1+3]), 
                   parse(Float64, SubString(s, i1+5, i2-1)), 
                   parse(Float64, SubString(s, i2+1, i3-1)),
                   long(s[i3+1]), long(s[i3+3]))
    end

    i = i3 + 4

    if i < ls
      while s[i] === ')'
        i += 1
      end
    end
  end

  return tree, i + 1
end




"""
    _iparse(s::String, i::Int64, ls::Int64, ::Type{iTb})
parse istring to `iTb`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{iTb})

  @inbounds begin

    inode = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, iTb)
      inode = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, iTb)
    end

    i1 = findnext(',', s, i  + 1)
    i2 = findnext(',', s, i1 + 1)
    i3 = findnext(',', s, i2 + 1)
    i4 = findnext(']', s, i3 + 1)

    if inode
      tree = iTb(sd1, sd2,
                 parse(Float64, SubString(s, i, i1-1)),
                 parse(Float64, SubString(s, i1+1, i2-1)),
                 parse(Float64, SubString(s, i2+1, i3-1)),
                 long(s[i3+1]), 
                 _iparse_v(s, i3+4, i4-1))
    else
      tree = iTb(parse(Float64, SubString(s, i, i1-1)),
                 parse(Float64, SubString(s, i1+1, i2-1)),
                 parse(Float64, SubString(s, i2+1, i3-1)),
                 long(s[i3+1]), 
                 _iparse_v(s, i3+4, i4-1))
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
               parse(Float64, SubString(s, i, i1-1)),
               parse(Float64, SubString(s, i1+1, i2-1)),
               parse(Float64, SubString(s, i2+1, i3-1)),
               long(s[i3+1]), 
               long(s[i3+3]), 
               _iparse_v(s, i3+6, i4-1))
    else
      tree = T(parse(Float64, SubString(s, i, i1-1)),
               parse(Float64, SubString(s, i1+1, i2-1)),
               parse(Float64, SubString(s, i2+1, i3-1)),
               long(s[i3+1]), 
               long(s[i3+3]), 
               _iparse_v(s, i3+6, i4-1))
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
                  parse(Float64, SubString(s, i, i1-1)),
                  parse(Float64, SubString(s, i1+1, i2-1)),
                  parse(Float64, SubString(s, i2+1, i3-1)),
                  long(s[i3+1]), 
                  long(s[i3+3]), 
                  _iparse_v(s, i3+6, i4-1),
                  _iparse_v(s, i4+3, i5-1))
    else
      tree = iTbd(parse(Float64, SubString(s, i, i1-1)),
                  parse(Float64, SubString(s, i1+1, i2-1)),
                  parse(Float64, SubString(s, i2+1, i3-1)),
                  long(s[i3+1]), 
                  long(s[i3+3]), 
                  _iparse_v(s, i3+6, i4-1),
                  _iparse_v(s, i4+3, i5-1))
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
                     parse(Float64, SubString(s, i, i1-1)),
                     parse(Float64, SubString(s, i1+1, i2-1)),
                     parse(Float64, SubString(s, i2+1,i3-1)),
                     long(s[i3+1]), 
                     long(s[i3+3]), 
                     long(s[i3+5]), 
                     _iparse_v(s, i3+8, i4-1),
                     _iparse_v(s, i4+3, i5-1))
      else
        tree = iTfbd(sd1,
                     parse(Float64, SubString(s, i, i1-1)),
                     parse(Float64, SubString(s, i1+1, i2-1)),
                     parse(Float64, SubString(s, i2+1,i3-1)),
                     long(s[i3+1]), 
                     long(s[i3+3]), 
                     long(s[i3+5]), 
                     _iparse_v(s, i3+8, i4-1),
                     _iparse_v(s, i4+3, i5-1))
      end
    else
      tree = iTfbd(parse(Float64, SubString(s, i, i1-1)),
                   parse(Float64, SubString(s, i1+1, i2-1)),
                   parse(Float64, SubString(s, i2+1,i3-1)),
                   long(s[i3+1]), 
                   long(s[i3+3]), 
                   long(s[i3+5]), 
                   _iparse_v(s, i3+8, i4-1),
                   _iparse_v(s, i4+3, i5-1))
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
                    parse(Float64, SubString(s, i, i1-1)),
                    parse(Float64, SubString(s, i1+1,i2-1)),
                    parse(Float64, SubString(s, i2+1, i3-1)),
                    _iparse_v(s, i3+2, i4-1),
                    _iparse_v(s, i4+3, i5-1))
      else
        tree = sTxs(sd1,
                    parse(Float64, SubString(s, i, i1-1)),
                    parse(Float64, SubString(s, i1+1,i2-1)),
                    parse(Float64, SubString(s, i2+1, i3-1)),
                    _iparse_v(s, i3+2, i4-1),
                    _iparse_v(s, i4+3, i5-1))
      end
    else
      tree = sTxs(parse(Float64, SubString(s, i, i1-1)),
                  parse(Float64, SubString(s, i1+1,i2-1)),
                  parse(Float64, SubString(s, i2+1, i3-1)),
                  _iparse_v(s, i3+2, i4-1),
                  _iparse_v(s, i4+3, i5-1))
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
    _iparse(s::String, i::Int64, ls::Int64, ::Type{iTpbd})
parse istring to `iT`.
"""
function _iparse(s::String, i::Int64, ls::Int64, ::Type{iTpbd})

  @inbounds begin

    in1 = false
    in2 = false

    if s[i] === '('
      sd1, i = _iparse(s, i + 1, ls, iTpbd)
      in1 = true
    end

    if s[i] === '('
      sd2, i = _iparse(s, i + 1, ls, iTpbd)
      in2 = true
    end

    i1 = findnext(',', s, i  + 1)
    i2 = findnext(',', s, i1 + 1)
    i3 = findnext(',', s, i2 + 1)
    i4 = findnext(']', s, i3 + 1)
    i5 = findnext(']', s, i4 + 1)
    i6 = findnext(']', s, i5 + 1)

    if in1
      if in2
        tree = iTpbd(sd1, sd2,
                     Pparse(Float64, s[i:i1-1]),
                     Pparse(Float64, s[i1+1:i2-1]),
                     Pparse(Float64, s[i2+1:i3-1]),
                     long(s[i3+1]), 
                     long(s[i3+3]), 
                     long(s[i3+5]), 
                     _iparse_v(s[i3+8:i4-1]),
                     _iparse_v(s[i4+3:i5-1]),
                     _iparse_v(s[i5+3:i6-1]))
      else
        tree = iTpbd(sd1,
                     Pparse(Float64, s[i:i1-1]),
                     Pparse(Float64, s[i1+1:i2-1]),
                     Pparse(Float64, s[i2+1:i3-1]),
                     long(s[i3+1]), 
                     long(s[i3+3]), 
                     long(s[i3+5]), 
                     _iparse_v(s[i3+8:i4-1]),
                     _iparse_v(s[i4+3:i5-1]),
                     _iparse_v(s[i5+3:i6-1]))
      end
    else
      tree = iTpbd(Pparse(Float64, s[i:i1-1]),
                   Pparse(Float64, s[i1+1:i2-1]),
                   Pparse(Float64, s[i2+1:i3-1]),
                   long(s[i3+1]), 
                   long(s[i3+3]), 
                   long(s[i3+5]), 
                   _iparse_v(s[i3+8:i4-1]),
                   _iparse_v(s[i4+3:i5-1]),
                   _iparse_v(s[i5+3:i6-1]))
    end

    i = i6 + 1

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
function _iparse_v(s::String, from::Int64, to::Int64)
  v = Float64[]
  i = from
  f = findnext(',', s, i)
  while !isnothing(f) && f < to
    push!(v, parse(Float64, SubString(s, i, f-1)))
    i = f + 1
    f = findnext(',', s, i)
  end
  push!(v, parse(Float64, SubString(s, i, to)))
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



