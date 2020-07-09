#=

insane tree read and write

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#


"""
    read_newick(out_file::String)

Writes `iTree` as a newick tree to `out_file`.
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

    tree = iTree(from_string(s1), from_string(s2), 0.0)
  # if root
  else
    tree = from_string(s)
  end

  return tree
end

# io = open(homedir()*"/Desktop/tree.tre")
# sr = readlines(io)[1]
# close(io)
# sr = sr[2:(findfirst(isequal(';'), sr)-2)]

# io = open(homedir()*"/Desktop/t.tre")
# sn = readlines(io)[1]
# close(io)
# sn = sn[2:(findfirst(isequal(';'), sn)-2)]



# s = sn

function from_string(s::String)

  # find pendant edge
  wd  = findlast(isequal(':'), s)
  pei = parse(Float64, s[(wd+1):end])
  s = s[2:(wd-2)]

  # if tip
  if !(occursin('(', s) || occursin(',', s))
      return iTree(pei)
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

  iTree(from_string(s1), from_string(s2), pei)
end








"""
    write_newick(tree::iTree, out_file::String)

Writes `iTree` as a newick tree to `out_file`.
"""
function write_newick(tree::iTree, out_file::String)

  s = to_string(tree, n = 0)
  s = string("(", s, ");")

  io = open(out_file*".tre", "w")
  write(io, s)
  close(io)

  return nothing
end




"""
    to_string(tree::iTree; n::Int64 = 0)

Returns newick string.
"""
function to_string(tree::iTree; n::Int64 = 0)

  if istip(tree.d1)
    if istip(tree.d2)
      n += 1
      s1 = string("(t",n,":",pe(tree.d1),",")
      n += 1
      s2 = string("t",n,":",pe(tree.d2),"):",pe(tree))
      return s1*s2
    else 
      n += 1
      return string("(t",n,":",pe(tree.d1), ",",
              to_string(tree.d2, n = n),"):", pe(tree))
    end
  elseif istip(tree.d2)
    n += 1
    return string("(", to_string(tree.d1, n = n), 
      ",t",n,":",pe(tree.d2), "):", pe(tree))
  else
    return string("(",to_string(tree.d1, n = n),",",
               to_string(tree.d2, n = sntn(tree.d1) + n),"):",pe(tree))
  end
end





