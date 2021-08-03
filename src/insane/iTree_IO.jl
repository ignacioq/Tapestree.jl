#=

insane tree read and write

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#




"""
    read_newick(in_file::String, ::Type{T}) where {T <: iTree} 

Reads a newick tree into `iTsimple` from `in_file`.
"""
function read_newick(in_file::String, ::Type{T}) where {T <: iTree} 

  io = open(in_file)
  s = readlines(io)[1]
  close(io)

  s = s[2:(findfirst(isequal(';'), s)-2)]

  np = 1
  nnp = 0
  for (i,v) in enumerate(s)
    if v == '('
      np += 1
    elseif v == ')'
      np -= 1
    end
    if v == ',' && isone(np)
      nnp += 1
    end
  end

  # if no root
  if isone(nnp)
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

    tree = T(from_string(s1, T), from_string(s2, T), 0.0)
  # if root
  else
    tree = from_string(s, T)
  end

  return tree
end




"""
    from_string(s::String, ::Type{T}) where {T <: iTree} )

Returns `iTree` from newick string.
"""
function from_string(s::String, ::Type{T}) where {T <: iTree}

  # find pendant edge
  wd  = findlast(isequal(':'), s)
  pei = parse(Float64, s[(wd+1):end])
  s = s[2:(wd-2)]

  # if tip
  if !(occursin('(', s) || occursin(',', s))
      return T(pei)
  else
    # estimate number of parentheses (when np returns to 1)
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

  T(from_string(s1, T), from_string(s2, T), pei)
end




"""
    write_newick(tree::iTsimple, out_file::String)

Writes `iTsimple` as a newick tree to `out_file`.
"""
function write_newick(tree::T, out_file::String) where {T <: iTree}

  s = to_string(tree, n = 0)
  s = string("(", s, ");")

  io = open(out_file*".tre", "w")
  write(io, s)
  close(io)

  return nothing
end




"""
    to_string(tree::T; n::Int64 = 0) where {T <: iTree})

Returns newick string.
"""
function to_string(tree::T; n::Int64 = 0) where {T <: iTree}

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
              to_string(tree.d2, n = n),"):", e(tree))
    end
  elseif istip(tree.d2)
    n += 1
    return string("(", to_string(tree.d1, n = n), 
      ",t",n,":",e(tree.d2), "):", e(tree))
  else
    return string("(",to_string(tree.d1, n = n),",",
               to_string(tree.d2, n = sntn(tree.d1, 0) + n),"):",e(tree))
  end
end





