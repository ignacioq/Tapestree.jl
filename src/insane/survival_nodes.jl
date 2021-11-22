#=

Survival nodes for conditioning

Ignacio Quintero Mächler

t(-_-t)

Created 16 11 2021
=#




"""
    sum_alone_stem!(tree::T, 
                    tna ::Float64, 
                    exx ::Bool,
                    sn  ::BitVector) where {T <: iTree}

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem!(tree::T, 
                         tna ::Float64, 
                         exx ::Bool,
                         sn  ::BitVector) where {T <: iTree}

  if isdefined(tree, :d1)

    et = e(tree)
    if isapprox(tna, et, atol = 1e-6)
      if exx
        push!(sn, true)
      else
        push!(sn, false)
      end
    elseif tna < et
      push!(sn, true)
    else
      push!(sn, false)
    end
    tna -= et

    if isfix(tree.d1::T)
      tnx = treeheight(tree.d2::T)
      if isapprox(tna, tnx, atol = 1e-6)
        exx = iszero(ntipsalive(tree.d2)) && exx
      elseif  tnx > tna
        exx = iszero(ntipsalive(tree.d2))
        tna = tnx
      end
      sum_alone_stem!(tree.d1::T,tna, exx, sn)
    else
      tnx = treeheight(tree.d1::T)
      if isapprox(tna, tnx, atol = 1e-6)
        exx = iszero(ntipsalive(tree.d1)) && exx
      elseif tnx > tna
        exx = iszero(ntipsalive(tree.d1))
        tna = tnx
      end
      sum_alone_stem!(tree.d2::T,tna, exx, sn)
    end
  end
end




"""
    sum_alone_stem_p!(tree::T, 
                      tna ::Float64, 
                      exx ::Bool,
                      sn  ::BitVector) where {T <: iTree}

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem_p!(tree::T, 
                           tna ::Float64, 
                           exx ::Bool,
                           sn  ::BitVector) where {T <: iTree}

  et = e(tree)
  if isapprox(tna, et, atol = 1e-6)
    if exx
      push!(sn, true)
    else
      push!(sn, false)
    end
  elseif tna < et
    push!(sn, true)
  else
    push!(sn, false)
  end
  tna -= et

  if isdefined(tree, :d1)
    if isfix(tree.d1::T)
      tnx = treeheight(tree.d2::T)
      if isapprox(tna, tnx, atol = 1e-6)
        exx = iszero(ntipsalive(tree.d2)) && exx
      elseif  tnx > tna
        exx = iszero(ntipsalive(tree.d2))
        tna = tnx
      end
      sum_alone_stem_p!(tree.d1::T,tna, exx, sn)
    else
      tnx = treeheight(tree.d1::T)
      if isapprox(tna, tnx, atol = 1e-6)
        exx = iszero(ntipsalive(tree.d1)) && exx
      elseif tnx > tna
        exx = iszero(ntipsalive(tree.d1))
        tna = tnx
      end
      sum_alone_stem_p!(tree.d2::T,tna, exx, sn)
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
      empty!(sns[1])
      sum_alone_stem_p!(psi[1], 0.0, false, sns[1])
    end
  else
    b1  = idf[1]
    d1i = d1(b1)
    d2i = d2(b1)

    if it(idf[d1i]) 
      if it(idf[d2i])
        f = let d1i = d1i, d2i = d2i
          function (psi::Vector{T}, sns::NTuple{3,BitVector})
            empty!(sns[2])
            empty!(sns[3])
            sum_alone_stem!(psi[d1i], 0.0, false, sns[2])
            sum_alone_stem!(psi[d2i], 0.0, false, sns[3])
          end
        end
      else
        f = let d1i = d1i, d2i = d2i
          function (psi::Vector{T}, sns::NTuple{3,BitVector})
            empty!(sns[2])
            empty!(sns[3])
            sum_alone_stem!(  psi[d1i], 0.0, false, sns[2])
            sum_alone_stem_p!(psi[d2i], 0.0, false, sns[3])
          end
        end 
      end
    elseif it(idf[d2i])
      f = let d1i = d1i, d2i = d2i
        function (psi::Vector{T}, sns::NTuple{3,BitVector})
          empty!(sns[2])
          empty!(sns[3])
          sum_alone_stem_p!(psi[d1i], 0.0, false, sns[2]) 
          sum_alone_stem!(  psi[d2i], 0.0, false, sns[3])
        end
      end
    else
      f = let d1i = d1i, d2i = d2i
        function (psi::Vector{T}, sns::NTuple{3,BitVector})
          empty!(sns[2])
          empty!(sns[3])
          sum_alone_stem_p!(psi[d1i], 0.0, false, sns[2]) 
          sum_alone_stem_p!(psi[d2i], 0.0, false, sns[3])
        end
      end
    end
  end

  return f
end




"""
    cond_surv_crown(tree::sTbd, λ::Float64, μ::Float64)

Log-probability of at least two lineage surviving for 
birth-death process with `λ` and `μ` for crown age.
"""
function cond_surv_crown(tree::sTbd, λ::Float64, μ::Float64)
  n = sum_alone_stem(tree.d1::sTbd, 0.0, 0.0) +
      sum_alone_stem(tree.d2::sTbd, 0.0, 0.0)
  return n*log((λ + μ)/λ) - log(λ)
end




"""
    cond_surv_stem(tree::sTbd, λ::Float64, μ::Float64)

Log-probability of at least one lineage surviving for 
birth-death process with `λ` and `μ` for stem age.
"""
function cond_surv_stem(tree::sTbd, λ::Float64, μ::Float64)
  n = sum_alone_stem(tree, 0.0, 0.0)
  return n*log((λ + μ)/λ)
end




"""
    sum_alone_stem(tree::sTbd, tna::Float64, n::Float64)

Count nodes in stem lineage when a diversification event could have 
returned an overall extinction.
"""
function sum_alone_stem(tree::sTbd, tna::Float64, n::Float64)

  if istip(tree)
    return n
  end

  if tna < e(tree)
    n += 1.0
  end
  tna -= e(tree)

  if isfix(tree.d1::sTbd)
    tnx = treeheight(tree.d2::sTbd)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem(tree.d1::sTbd, tna, n)
  else
    tnx = treeheight(tree.d1::sTbd)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem(tree.d2::sTbd, tna, n)
  end

end




"""
    cond_surv_stem_p(tree::sTbd, λ::Float64, μ::Float64)

Log-probability of at least one lineage surviving after time `t` for 
birth-death process with `λ` and `μ` for stem age.
"""
function cond_surv_stem_p(tree::sTbd, λ::Float64, μ::Float64)
  n = sum_alone_stem_p(tree, 0.0, 0.0)
  return n*log((λ + μ)/λ)
end




"""
    sum_alone_stem_p(tree::sTbd, tna::Float64, n::Float64)

Count nodes in stem lineage when a diversification event could have 
returned an overall extinction.
"""
function sum_alone_stem_p(tree::sTbd, tna::Float64, n::Float64)

  if tna < e(tree)
    n += 1.0
  end
  tna -= e(tree)

  if istip(tree)
    return n
  end

  if isfix(tree.d1::sTbd)
    tnx = treeheight(tree.d2::sTbd)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem_p(tree.d1::sTbd, tna, n)
  else
    tnx = treeheight(tree.d1::sTbd)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem_p(tree.d2::sTbd, tna, n)
  end
end




"""
    sum_alone_stem(tree::iTgbmce, 
                   tna ::Float64,
                   exx ::Bool, 
                   ll  ::Float64,
                   μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem(tree::iTgbmce, 
                        tna ::Float64,
                        exx ::Bool, 
                        ll  ::Float64,
                        μ   ::Float64)

  if istip(tree)
    return ll
  end

  et = e(tree)
  if isapprox(tna, et, atol = 1e-6)
    if exx
      λi  = lλ(tree)[end]
      ll += log(exp(λi) + μ) - λi
    end
  elseif tna < et
    λi  = lλ(tree)[end]
    ll += log(exp(λi) + μ) - λi
  end
  tna -= et

  if isfix(tree.d1::iTgbmce)
    tnx = treeheight(tree.d2::iTgbmce)
    if isapprox(tna, tnx, atol = 1e-6)
      exx = iszero(ntipsalive(tree.d2)) && exx
    elseif  tnx > tna
      exx = iszero(ntipsalive(tree.d2))
      tna = tnx
    end
    sum_alone_stem(tree.d1::iTgbmce, tna, exx, ll, μ)
  else
    tnx = treeheight(tree.d1::iTgbmce)
    if isapprox(tna, tnx, atol = 1e-6)
      exx = iszero(ntipsalive(tree.d1)) && exx
    elseif tnx > tna
      exx = iszero(ntipsalive(tree.d1))
      tna = tnx
    end
    sum_alone_stem(tree.d2::iTgbmce, tna, exx, ll, μ)
  end
end




"""
    sum_alone_stem_p(tree::iTgbmce, 
                     tna ::Float64, 
                     exx ::Bool,
                     ll  ::Float64, 
                     μ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem_p(tree::iTgbmce, 
                          tna ::Float64, 
                          exx ::Bool,
                          ll  ::Float64, 
                          μ   ::Float64)

  et = e(tree)
  if isapprox(tna, et, atol = 1e-6)
    if exx
      λi  = lλ(tree)[end]
      ll += log(exp(λi) + μ) - λi
    end
  elseif tna < et
    λi  = lλ(tree)[end]
    ll += log(exp(λi) + μ) - λi
  end
  tna -= et

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmce)
    tnx = treeheight(tree.d2::iTgbmce)
    if isapprox(tna, tnx, atol = 1e-6)
      exx = iszero(ntipsalive(tree.d2)) && exx
    elseif  tnx > tna
      exx = iszero(ntipsalive(tree.d2))
      tna = tnx
    end
    sum_alone_stem_p(tree.d1::iTgbmce, tna, exx, ll, μ)
  else
    tnx = treeheight(tree.d1::iTgbmce)
    if isapprox(tna, tnx, atol = 1e-6)
      exx = iszero(ntipsalive(tree.d1)) && exx
    elseif tnx > tna
      exx = iszero(ntipsalive(tree.d1))
      tna = tnx
    end
    sum_alone_stem_p(tree.d2::iTgbmce, tna, exx, ll, μ)
  end
end



