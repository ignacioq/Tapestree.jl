#=

Survival nodes for conditioning

Ignacio Quintero Mächler

t(-_-t)

Created 16 11 2021
=#




"""
    sum_alone_stem!(tree::T, 
                     tna ::Float64, 
                     ll  ::Float64, 
                     μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem!(tree::T, tna::Float64, sn::BitVector) where {T <: iTree}

  if isdefined(tree, :d1)

    if tna < e(tree)
      push!(sn, true)
    else
      push!(sn, false)
    end
    tna -= e(tree)

    if isfix(tree.d1::T)
      tnx = treeheight(tree.d2::T)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem!(tree.d1::T, tna, sn)
    else
      tnx = treeheight(tree.d1::T)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem!(tree.d2::T, tna, sn)
    end
  end
end




"""
    sum_alone_stem_p!(tree::T, 
                     tna ::Float64, 
                     ll  ::Float64, 
                     μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem_p!(tree::T, 
                          tna ::Float64, 
                          sn  ::BitVector) where {T <: iTree}

  if tna < e(tree)
    push!(sn, true)
  else
    push!(sn, false)
  end
  tna -= e(tree)

  if isdefined(tree, :d1)
    if isfix(tree.d1::T)
      tnx = treeheight(tree.d2::T)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem_p!(tree.d1::T, tna, sn)
    else
      tnx = treeheight(tree.d1::T)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem_p!(tree.d2::T, tna, sn)
    end
  end
end



"""
    make_snodes(idf::Vector{iBffs}, stem::Bool, ::Type{T})

Make closure for conditioning function
"""
function make_snodes(idf::Vector{iBffs}, stem::Bool, ::Type{T}) where {T <: iTree}

  # conditioning
  if stem
    function f(psi::Vector{T}, sns::NTuple{3,BitVector})
      sum_alone_stem_p!(psi[1], 0.0, sns[1])
    end
  else
    b1  = idf[1]
    d1i = d1(b1)
    d2i = d2(b1)
    t1  = it(idf[d1i])
    t2  = it(idf[d2i])

    if t1 
      if t2
        f = let d1i = d1i, d2i = d2i
          function (psi::Vector{T}, sns::NTuple{3,BitVector})
            empty!(sns[2])
            empty!(sns[3])
            sum_alone_stem!(psi[d1i], 0.0, sns[2])
            sum_alone_stem!(psi[d2i], 0.0, sns[3])
          end
        end
      else
        f = let d1i = d1i, d2i = d2i
          function (psi::Vector{T}, sns::NTuple{3,BitVector})
            empty!(sns[2])
            empty!(sns[3])
            sum_alone_stem!(  psi[d1i], 0.0, sns[2])
            sum_alone_stem_p!(psi[d2i], 0.0, sns[3])
          end
        end 
      end
    elseif t2
      f = let d1i = d1i, d2i = d2i
        function (psi::Vector{T}, sns::NTuple{3,BitVector})
          empty!(sns[2])
          empty!(sns[3])
          sum_alone_stem_p!(psi[d1i], 0.0, sns[2]) 
          sum_alone_stem!(  psi[d2i], 0.0, sns[3])
        end
      end
    else
      f = let d1i = d1i, d2i = d2i
        function (psi::Vector{T}, sns::NTuple{3,BitVector})
          empty!(sns[2])
          empty!(sns[3])
          sum_alone_stem_p!(psi[d1i], 0.0, sns[2]) 
          sum_alone_stem_p!(psi[d2i], 0.0, sns[3])
        end
      end
    end
  end

  return f
end
