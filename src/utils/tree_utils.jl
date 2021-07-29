#=

Tree utilities for ESSE

Ignacio Quintero Mächler

t(-_-t)

September 19 2017

=#




"""
Immutable type of an R tree `phylo` object type.
"""
struct rtree
  ed  ::Array{Int64,2}
  el  ::Array{Float64,1}
  tlab::Array{String,1}
  nnod::Int64
end





"""
    read_tree(tree_file::String; 
              order::String = "cladewise", 
              branching_times::Bool = true)

Function to read a tree using `RCall`
to call **ape** tree reading capabilities. 
"""
function read_tree(tree_file      ::String; 
                   order          ::String = "cladewise", 
                   branching_times::Bool = true)

  str = reval("""
                library(\"ape\")
                tree     <- read.tree('$tree_file') 
                tree     <- reorder(tree, order = '$order')
                edge     <- .subset2(tree,'edge')
                Nnode    <- .subset2(tree,'Nnode')
                tiplabel <- .subset2(tree,'tip.label')
                edlength <- .subset2(tree,'edge.length')
                list(edge,Nnode,tiplabel,edlength)
              """)

  edge     = rcopy(str[1])
  edge     = convert(Array{Int64},edge)
  Nnode    = rcopy(str[2])
  Nnode    = convert(Int64,Nnode)
  tiplabel = rcopy(str[3])
  edlength = rcopy(str[4])
  edlength = convert(Array{Float64},edlength)

  tree = rtree(edge, edlength, tiplabel, Nnode)

  if branching_times
    brtimes = reval("""
                      brtimes <- branching.times(tree)
                    """)
    brtimes = rcopy(brtimes)
    return tree, brtimes
  else
    return tree
  end
end






"""
    make_ape_tree(n::Int64, 
                  λ::Float64, 
                  μ::Float64; 
                  order::String = "cladewise", 
                  branching_times::Bool = true)

Make a phylogenetic tree using `phytools` in R. 
"""
function make_ape_tree(n              ::Int64, 
                       λ              ::Float64, 
                       μ              ::Float64; 
                       order          ::String = "cladewise", 
                       branching_times::Bool   = true)

  str = reval("""
                library(ape)
                library(phytools)
                tree     <- pbtree(n = $n, b = $λ, d = $μ, extant.only = TRUE)
                tree     <- reorder(tree, order = '$order')
                edge     <- .subset2(tree,'edge')
                Nnode    <- .subset2(tree,'Nnode')
                tiplabel <- .subset2(tree,'tip.label')
                edlength <- .subset2(tree,'edge.length')
                list(edge,Nnode,tiplabel,edlength)
              """)

  edge     = rcopy(str[1])
  edge     = convert(Array{Int64},edge)
  Nnode    = rcopy(str[2])
  Nnode    = convert(Int64,Nnode)
  tiplabel = rcopy(str[3])
  edlength = rcopy(str[4])
  edlength = convert(Array{Float64},edlength)

  tree = rtree(edge, edlength, tiplabel, Nnode)

  if branching_times
    brtimes = reval("""
                      brtimes <- branching.times(tree)
                    """)
    brtimes = rcopy(brtimes)
    return tree, brtimes
  else
    return tree
  end
end





"""
    maketriads(ed::Array{Int64,2})

Make edge triads given the tree. The first number is the parent, 
the second and third the children.
"""
function maketriads(ed::Array{Int64,2}; rev::Bool = false)

  # internal nodes
  ins = unique(ed[:,1])[1:(end-1)]::Array{Int64,1}

  rev && sort!(ins, rev = true)

  ed1 = ed[:,1]
  ed2 = ed[:,2]

  trios = Array{Int64,1}[]

  # for all internal nodes
  for i in ins
    daus = findall(isequal(i), ed1)
    pushfirst!(daus, findfirst(isequal(i), ed2))
    push!(trios, daus)
  end

  return trios::Array{Array{Int64,1},1}
end





"""
    abs_time_branches(el  ::Array{Float64,1}, 
                      ed  ::Array{Int64,2},
                      ntip::Int64)

Make array with absolute initial and end time for each 
branch. Time goes backwards with the present being `0.0`.
"""
function abs_time_branches(el  ::Array{Float64,1}, 
                           ed  ::Array{Int64,2},
                           ntip::Int64)

  @inbounds begin
    # make real time edges
    elrt = zeros(size(ed))

    # edges 1
    ed1 = ed[:,1]

    for i in axes(elrt,1)
      # is not tip
      if ed[i,2] > ntip
        elrt[i,2] = elrt[findfirst(isequal(ed[i,2]), ed1)::Int64,1]
        elrt[i,1] = elrt[i,2] + el[i]
      else
        elrt[i,1] = el[i]
      end
    end

  end

  return elrt
end





"""
    brts(el  ::Array{Float64,1}, 
                ed  ::Array{Int64,2},
                ntip::Int64)

Get branching times for a tree in "cladewise" order. 
Time goes backwards with the present being `0.0`.
"""
function brts(el  ::Array{Float64,1}, 
              ed  ::Array{Int64,2},
              ntip::Int64)

  @inbounds begin
    # make real time edges

    # edges 1
    ed1 = ed[:,1]
    ed2 = ed[:,2]

    ne   = lastindex(ed1)
    nn   = zeros(ntip-1)
    intn = findall(map(x -> x > ntip, ed2))

    for i in intn
      nn[ed2[i]-ntip] = nn[ed1[i] - ntip] + el[i]
    end

    # tree height
    trh = nn[ed1[ne] - ntip] + el[ne]
    for i in Base.OneTo(ntip-1)
      nn[i] = trh - nn[i]
    end

  end

  return nn
end





"""
    tree_height(el  ::Array{Float64,1}, 
                ed  ::Array{Int64,2},
                ntip::Int64)

Estimate tree height.
"""
function tree_height(el  ::Array{Float64,1}, 
                     ed  ::Array{Int64,2},
                     ntip::Int64)

  @inbounds begin
    # tree height
    th  = 0.0
    ed2 = ed[:,2]

    da::Int64 = findfirst(isequal(1), ed2)

    # if the first branch reaches the tree height
    ed[da,1] == (ntip+1) && return el[da]

    while true
      th += el[da]

      pr = findfirst(isequal(ed[da,1]), ed2)::Int64

      if ed[pr,1] == (ntip+1)
        th += el[pr]
        break
      end

      da = findfirst(isequal(ed[pr,2]), ed2)::Int64
    end
  end

  return th
end





"""
    tip_dictionary(tS::Array{Int64,1})

Create a dictionary. WARNING: ONLY FOR USE WITHOUT CARING ABOUT TOPOLOGY.
"""
function tip_dictionary(tS::Array{Int64,1})

  # make tip values Dictionary
  tv = Dict{Int64, Int64}()

  for i in Base.OneTo(lastindex(tv))
    push!(tv, i => tS[i])
  end

  return tv
end




"""
    numberedges(ed::Array{Int64,2}, tN::Array{Int64,1})

Change numbering scheme so that tips are `1:numerofspecies` followed
by node numbering. MRCA is `numerofspecies+1`.
"""
function numberedges(ed::Array{Int64,2}, tN::Array{Int64,1}, tS::Array{Int64,1})
  nt = 1

  ni = lastindex(tN) + 1

  edc = zeros(Int64,size(ed))
  edc[1,1] = edc[2,1] = ni
  ni += 1

  e1 = ed[:,1]

  # make tip values Dictionary
  tv = Dict{Int64, Int64}()

  for i in axes(ed,1)

    nn = ed[i,2]
    ww = findfirst(isequal(nn), tN)

    # if tip
    if !isnothing(ww)
      edc[i,2] = nt
      push!(tv, nt => tS[ww])
      nt += 1
    else
      ii = findfirst(isequal(nn), e1)
      edc[ii,1] = edc[ii+1,1] = edc[i,2] = ni
      ni  += 1 
    end

  end

  return edc, tv
end




"""
    postorderedges(ed  ::Array{Int64,2},
                   el  ::Array{Float64,1},
                   ntip::Int64)

Organize edges, edge lengths in postorder traversal fashion.
"""
function postorderedges(ed  ::Array{Int64,2},
                        el  ::Array{Float64,1},
                        ntip::Int64)

  # post-order transversal using 2 stacks
  s1 = [ed[1]]
  s2 = Int64[]

  while lastindex(s1) > 0
    nod = pop!(s1)
    push!(s2, nod)

    wn = findfirst(isequal(nod), ed[:,1])
    if isnothing(wn)
      continue
    else
      push!(s1, ed[wn,2],ed[(wn+1),2])
    end
  end

  # rearrange edges accordingly
  indx = deleteat!(indexin(reverse(s2), ed[:,2]), size(ed,1)+1)

  ed = ed[indx,:]
  el = el[indx]

  # advance nodes with only daughter tips
  tnd = Int64[]
  ndp = Int64[]
  for nd in unique(ed[:,1])
    fed = findall(ed[:,1 ] .== nd)
    if length(filter(x -> x <= ntip, ed[fed,2])) == 2
      push!(tnd, nd)
      push!(ndp, fed[1], fed[2])
    end
  end

  append!(ndp,setdiff(1:(2ntip-2), ndp))

  ed[ndp,:]

  return ed[ndp,:], el[ndp]
end





"""
    remove_extinct(ed::Array{Int64,2}, 
                   el::Array{Float64,1}, 
                   ee::Array{Int64,1})

Remove extinct nodes from tree given extinct edges `ee`.
"""
function remove_extinct(ed::Array{Int64,2}, 
                        el::Array{Float64,1}, 
                        ee::Array{Int64,1})

  # identify extinct nodes
  nse = ed[ee,2]

  for i in nse

    @views e1 = ed[:,1]
    @views e2 = ed[:,2]
    r = findfirst(x -> x === i, e2)

    # parent node
    pn = e1[r] 

    # other row
    if lastindex(e1) === r
      or = r-1
    else
      or = e1[r+1] === pn ? r+1 : r-1
    end

    # ancestral row
    ar = findfirst(x -> x === pn, e2) 

    if !isnothing(ar)
      # assign node
      e2[ar]  = e2[or]
      # add edge length
      el[ar] += el[or]
    end

    # remove
    news = setdiff(1:length(e2), r, or)
    ed   = ed[news,:]   # remove from edges
    el   = el[news]     # remove from edge lengths
  end

  return ed, el
end



