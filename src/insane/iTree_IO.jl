#=

insane tree read and write

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#



"""
    read_newick(in_file::String)

Reads a newick tree into `iTsimple` from `in_file`.
"""
function read_newick(in_file::String)

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

    tree = sT_label(from_string(s1), from_string(s2), 0.0, "")
  # if root
  else
    tree = from_string(s)
  end

  return tree
end




"""
    from_string(s::String)

Returns `iTree` from newick string.
"""
function from_string(s::String)

  # find pendant edge
  wd  = findlast(isequal(':'), s)
  pei = parse(Float64, s[(wd+1):end])
  lab = s[1:(wd-1)]
  s   = s[2:(wd-2)]

  # if tip
  if !(occursin('(', s) || occursin(',', s))
      return sT_label(pei, lab)
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

  sT_label(from_string(s1), from_string(s2), pei, "")
end




"""
    from_string(s::String, ::Type{T}) where {T <: sTfbd}

Returns `sTfbd` from newick string.
"""
function from_string(s::String, ::Type{T}) where {T <: sTfbd}

  # find pendant edge
  wd = findlast(isequal(':'), s)
  pei = parse(Float64, s[(wd+1):end])

  # if tip
  if isone(count(i->(i==':'), s))
      return sTfbd(pei)
  else
    while s[wd-1] != ')'  wd -= 1  end
    s = s[2:(wd-2)]

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

  if isempty(s1)
    return sTfbd(from_string(s2, sTfbd), pei, false, true, false)
  elseif isempty(s2)
    return sTfbd(from_string(s1, sTfbd), pei, false, true, false)
  else
    return sTfbd(from_string(s1, sTfbd), from_string(s2, sTfbd), pei)
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
function to_string(tree::T; n::Int64 = 0) where {T <: iTree}

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
              to_string(tree.d2, n = n),"):", e(tree))
    end
  elseif istip(tree.d2)
    n += 1
    return string("(", to_string(tree.d1, n = n), 
      ",t",n,":",e(tree.d2), "):", e(tree))
  else
    return string("(",to_string(tree.d1, n = n),",",
               to_string(tree.d2, n = ntips(tree.d1) + n),"):",e(tree))
  end
end




"""
    to_string(tree::sTfbd; n::Int64=0, sa::Int64=0)

Returns newick string.
"""
function to_string(tree::sTfbd; n::Int64=0, sa::Int64=0)

  if istip(tree)
    return(string("t1:",e(tree)))
  end

  if isdefined(tree, :d1) && isdefined(tree, :d2)
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
                to_string(tree.d2, n=n, sa=sa),"):", e(tree))
      end
    elseif istip(tree.d2)
      n += 1
      return string("(", to_string(tree.d1, n=n, sa=sa), 
        ",t",n,":",e(tree.d2), "):", e(tree))
    else
      return string("(",to_string(tree.d1, n=n, sa=sa),",",
                        to_string(tree.d2, n=ntips(tree.d1) + n, 
                                  sa=nsampledancestors(tree.d1) + sa),"):",
                        e(tree))
    end
  
  # sampled ancestors
  elseif isdefined(tree, :d1)
    sa += 1
    if istip(tree.d1)
      n += 1
      return string("(t",n,":",e(tree.d1),")sa",sa,":", e(tree))
    else
      return string("(",to_string(tree.d1, n=n, sa=sa),")sa",sa,":",e(tree))
    end
  else
    sa += 1
    if istip(tree.d2)
      n += 1
      return string("(t",n,":",e(tree.d2),")sa",sa,":", e(tree))
    else
      return string("(",to_string(tree.d2, n=n, sa=sa),")sa",sa,":",e(tree))
    end
  end
end





