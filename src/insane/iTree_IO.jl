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

  # change the tree type if it has fossils (as sampled ancestors)
  T = (2*(nbif-nlin+1) == np) ? sT_label : sTf_label

  # if root (starts with one stem lineage)
  if isone(nlin)
    tree = from_string(s, T)
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

    tree = T(from_string(s1, T), from_string(s2, T), 0.0, "")
  end

  return tree
end




"""
    from_string(s::String, ::Type{T}) where {T <: sT}

Returns `iTree` from newick string.
"""
function from_string(s::String, ::Type{T}) where {T <: sT}

  # find pendant edge
  wd  = findlast(isequal(':'), s)
  pei = parse(Float64, s[(wd+1):end])
  lab = s[1:(wd-1)]
  #s   = s[2:(wd-2)]

  # if tip
  if isone(count(i->(i==':'), s))
      return T(pei, lab)
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
    return T(from_string(s2, T), pei, "")
  elseif isempty(s2)
    return T(from_string(s1, T), pei, "")
  else
    return T(from_string(s1, T), from_string(s2, T), pei, "")
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





