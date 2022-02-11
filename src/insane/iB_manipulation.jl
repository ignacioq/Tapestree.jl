#=

iB manipulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 07 2020
=#



"""
    setpa!(id::iBffs, pa::Int64)

Set parent branch.
"""
setpa!(id::iBffs, pa::Int64) = id.pa[] = pa




"""
    setd1!(id::iBffs, d1::Int64)

Set daughter 1 branch.
"""
setd1!(id::iBffs, d1::Int64) = id.d1[] = d1




"""
    setd2!(id::iBffs, d2::Int64)

Set daughter 2 branch
"""
setd2!(id::iBffs, d2::Int64) = id.d2[] = d2




"""
    setni!(id::iBffs, ni::Int64) 

Set number of alive lineages at present. 
"""
setni!(id::iBffs, ni::Int64) = id.ni[] = ni




"""
    setnt!(id::iBffs, nt::Int64)

Set number of alive lineages at time `t`. 
"""
setnt!(id::iBffs, nt::Int64) = id.nt[] = nt




"""
    setλt!(id::iBffs, λt::Float64) 

Set number of alive lineages at time `t`. 
"""
setλt!(id::iBffs, λt::Float64) = id.λt[] = λt




"""
    setμt!(id::iBffs, μt::Float64)

Set number of alive lineages at time `t`. 
"""
setμt!(id::iBffs, μt::Float64) = id.μt[] = μt




"""
    setψt!(id::iBffs, ψt::Float64)

Set number of alive lineages at time `t`. 
"""
setψt!(id::iBffs, ψt::Float64) = id.μt[] = ψt




"""
    addda!(id::iBf)

Add `1` to number of data augmented branches. 
"""
addda!(id::iBf) = id.da[] += 1




"""
    rmda!(id::iBf)

Substract `1` to number of data augmented branches. 
"""
rmda!(id::iBf) = id.da[] -= 1




"""
    randbranch(th   ::Float64,
               sth  ::Float64,
               idf  ::Array{iB,1},
               wbr  ::BitArray{1})

Randomly select branch to graft a subtree with height `sth` onto tree 
with height `th`.
"""
function randbranch(th   ::Float64,
                    sth  ::Float64,
                    idf  ::Array{iBf,1},
                    wbr  ::BitArray{1})

  # uniformly set the height at which to insert `stree`
  h = sth + rand()*(th - sth)

  # get which branches are cut by `h` and sample one at random
  nb  = branchescut!(h, wbr, idf)
  rn  = rand(Base.OneTo(nb))
  br, i  = getbranch(rn, wbr, idf)

  return h, br, i, nb
end





"""
    randbranch(th ::Float64,
               sth::Float64,
               wbf::BitArray{1},
               wba::BitArray{1},
               idf::Array{iBf,1},
               ida::Array{iBa,1})

Randomly select branch to graft a subtree with height `sth` onto tree 
with height `th`.
"""
function randbranch(th ::Float64,
                    sth::Float64,
                    wbf::BitArray{1},
                    wba::BitArray{1},
                    idf::Array{iBf,1},
                    ida::Array{iBa,1})

  # uniformly set the height at which to insert `stree`
  h = sth + rand()*(th - sth)

  # get which branches are cut by `h` and sample one at random
  nf, na = branchescut!(h, wbf, wba, idf, ida)
  # which branch
  rn  = rand(1:(nf + na))

  return h, nf, na, rn
end




"""
    getbranch(rn ::Int64, 
              wbf::BitArray{1}, 
              idf::Array{iBf,1})::Tuple{iBf,Int64}

Sample one branch among those in `wbf` that are `true` according to 
randomly picked `rn`.
"""
function getbranch(rn ::Int64, 
                   wbf::BitArray{1}, 
                   idf::Array{iBf,1})::Tuple{iBf,Int64}

  i::Int64 = 0
  n::Int64 = 0
  for bit in wbf
    i += 1
    if bit
      n += 1
      if n == rn
        return (idf[i], i)::Tuple{iBf,Int64}
      end
    end
  end
end




"""
    getbranch(rn ::Int64, 
              wba::BitArray{1}, 
              idf::Array{iBa,1})::Tuple{iBa,Int64}

Sample one branch among those in `wba` that are `true` according to 
randomly picked `rn`.
"""
function getbranch(rn ::Int64, 
                   wba::BitArray{1}, 
                   ida::Array{iBa,1})::Tuple{iBa,Int64}

  i::Int64 = 0
  n::Int64 = 0
  for bit in wba
    i += 1
    if bit
      n += 1
      if n == rn
        return (ida[i], i)::Tuple{iBa,Int64}
      end
    end
  end
end




"""
    branchescut!(h  ::Float64, 
                 wbf::BitArray{1},
                 idf::Array{iBf,1})::Int64

Fill the bit array with branches that are cut by  `h`.
"""
function branchescut!(h  ::Float64, 
                      wbf::BitArray{1},
                      idf::Array{iBf,1})::Int64

  @inbounds begin
    n::Int64 = 0
    i::Int64 = 0
    for b in idf
      i += 1
      if ti(b) > h > tf(b)
        wbf[i] = true
        n += 1
      else
        wbf[i] = false
      end
    end
  end

  return n
end





"""
    branchescut!(h  ::Float64,
                 wbf::BitArray{1}, 
                 wba::BitArray{1},
                 idf::Array{iBf,1},
                 ida::Array{iBa,1})::Int64

Fill the bit arrays with both fixed and augmented branches 
that are cut by  `h`.
"""
function branchescut!(h  ::Float64,
                      wbf::BitArray{1}, 
                      wba::BitArray{1},
                      idf::Array{iBf,1},
                      ida::Array{iBa,1})

  @inbounds begin

    i ::Int64 = 0
    nf::Int64 = 0
    for b in idf
      i += 1
      if ti(b) > h > tf(b)
        wbf[i] = true
        nf += 1
      else
        wbf[i] = false
      end
    end

    i = 0
    na::Int64 = 0
    for b in ida
      i += 1
      if ti(b) > h > tf(b)
        wba[i] = true
        na += 1
      else
        wba[i] = false
      end
    end

  end

  return nf, na
end




"""
    make_inodes(idf::Array{B, 1}) where B <: iBf

Return all the internal node indices for a given `iBf` vector and a vector
that is true for which ever daughter is a tip.
"""
function make_inodes(idf::Array{B, 1}) where B <: iBf

  inodes   = Int64[]
  terminus = BitArray{1}[]

  for i in 1:length(idf)
    pr  = i
    drpr = dr(idf[i])
    drd1 = push!(copy(drpr), true)
    drd2 = push!(copy(drpr), false)

    # d1
    d1 = findfirst(x -> dr(x) == drd1, idf)
    # d2
    d2 = findfirst(x -> dr(x) == drd2, idf) 

    if !isnothing(d1) && !isnothing(d2)
      push!(inodes, pr)

      bit = BitArray{1}([false, false])

      # check if either of the tips are terminal
      drd11 = push!(copy(drd1), true)
      drd12 = push!(copy(drd1), false)
      # d1
      d11 = findfirst(x -> dr(x) == drd11, idf)
      # d2
      d12 = findfirst(x -> dr(x) == drd12, idf) 

      if isnothing(d11) && isnothing(d12) 
        bit[1] = true
      end

      drd21 = push!(copy(drd2), true)
      drd22 = push!(copy(drd2), false)
      # d1
      d21 = findfirst(x -> dr(x) == drd21, idf)
      # d2
      d22 = findfirst(x -> dr(x) == drd22, idf) 

      if isnothing(d21) && isnothing(d22) 
        bit[2] = true
      end
      push!(terminus, bit)
    end
  end

  return inodes, terminus
end




"""
    make_triads(idf::Array{iBf, 1})

Return parent and two daughter indices for a given `iB` vector and a vector
that is true for which ever daughter is a tip.
"""
function make_triads(idf::Array{iBffs, 1})

  triads   = Array{Int64,1}[]
  terminus = BitArray{1}[]
  btotriad = Int64[]

  ii = 0

  for i in 1:length(idf)
    pr  = i
    drpr = dr(idf[i])
    drd1 = push!(copy(drpr), true)
    drd2 = push!(copy(drpr), false)

    d1 = findfirst(x -> dr(x) == drd1, idf)
    d2 = findfirst(x -> dr(x) == drd2, idf) 

    if !isnothing(d1) && !isnothing(d2)
      # push to triads
      push!(triads, [pr, d1, d2])

      # check terminus
      bit = BitArray{1}([it(idf[d1]), it(idf[d2])])
      push!(terminus, bit)

      # push to btotriad
      ii += 1
      push!(btotriad, ii)
    else
      # push to btotriad
      push!(btotriad, 0)
    end
  end

  return triads, terminus, btotriad
end


