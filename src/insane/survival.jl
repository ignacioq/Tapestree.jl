#=

Survival nodes for conditioning

Ignacio Quintero Mächler

t(-_-t)

Created 16 11 2021
=#

"""
    surv_m(t::Float64, λ::Float64, μ::Float64, ntry::Int64)

Sample the total number of `m` trials until both simulations survive.
"""
function surv_m(t::Float64, λ::Float64, μ::Float64, ntry::Int64)

  ntries = 0
  m      = 0.0
  while true
    m      += 1.0
    ntries += 1

    t1, s1, n1 = sim_cbd_surv(t, λ, μ, false, 1)
    t2, s2, n2 = sim_cbd_surv(t, λ, μ, false, 1)

    s1 && s2 && break

  
    ntries == ntry && break
  end

  return m
end



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
                         sn  ::BitVector) where {T <: iTree}

  if isdefined(tree, :d1)

    et = e(tree)
    if tna < et
      push!(sn, true)
    else
      push!(sn, false)
    end
    tna -= et

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
                      sn  ::BitVector) where {T <: iTree}

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem_p!(tree::T, 
                           tna ::Float64, 
                           sn  ::BitVector) where {T <: iTree}

  et = e(tree)
  if tna < et
    push!(sn, true)
  else
    push!(sn, false)
  end
  tna -= et

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
      empty!(sns[1])
      sum_alone_stem_p!(psi[1], 0.0, sns[1])
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
    elseif it(idf[d2i])
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
                        ll  ::Float64,
                        μ   ::Float64)

  if istip(tree)
    return ll
  end

  et = e(tree)
  if tna < et
    λi  = lλ(tree)[end]
    ll += log(exp(λi) + μ) - λi
  end
  tna -= et

  if isfix(tree.d1::iTgbmce)
    tnx = treeheight(tree.d2::iTgbmce)
    if  tnx > tna
      tna = tnx
    end
    sum_alone_stem(tree.d1::iTgbmce, tna, ll, μ)
  else
    tnx = treeheight(tree.d1::iTgbmce)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem(tree.d2::iTgbmce, tna, ll, μ)
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
                          ll  ::Float64, 
                          μ   ::Float64)

  et = e(tree)
  if tna < et
    λi  = lλ(tree)[end]
    ll += log(exp(λi) + μ) - λi
  end
  tna -= et

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmce)
    tnx = treeheight(tree.d2::iTgbmce)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem_p(tree.d1::iTgbmce, tna, ll, μ)
  else
    tnx = treeheight(tree.d1::iTgbmce)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem_p(tree.d2::iTgbmce, tna, ll, μ)
  end
end




"""
    sum_alone_stem(tree::iTgbmct, 
                   tna ::Float64,
                   ll  ::Float64,
                   ϵ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem(tree::iTgbmct, 
                        tna ::Float64,
                        ll  ::Float64,
                        ϵ   ::Float64)

  if istip(tree)
    return ll
  end

  et = e(tree)
  if tna < et
    ll += log(1.0 + ϵ)
  end
  tna -= et

  if isfix(tree.d1::iTgbmct)
    tnx = treeheight(tree.d2::iTgbmct)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem(tree.d1::iTgbmct, tna, ll, ϵ)
  else
    tnx = treeheight(tree.d1::iTgbmct)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem(tree.d2::iTgbmct, tna, ll, ϵ)
  end
end




"""
    sum_alone_stem_p(tree::iTgbmct, 
                     tna ::Float64, 
                     ll  ::Float64, 
                     ϵ   ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem_p(tree::iTgbmct, 
                          tna ::Float64, 
                          ll  ::Float64, 
                          ϵ   ::Float64)

  et = e(tree)
  if tna < et
    ll += log(1.0 + ϵ)
  end
  tna -= et

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmct)
    tnx = treeheight(tree.d2::iTgbmct)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem_p(tree.d1::iTgbmct, tna, ll, ϵ)
  else
    tnx = treeheight(tree.d1::iTgbmct)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem_p(tree.d2::iTgbmct, tna, ll, ϵ)
  end
end




"""
    sum_alone_stem(tree::iTgbmbd, 
                   tna ::Float64,
                   exx ::Bool, 
                   ll  ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem(tree::iTgbmbd, 
                        tna ::Float64,
                        ll  ::Float64)

  if istip(tree)
    return ll
  end

  et = e(tree)
  if tna < et
    λv  = lλ(tree)
    l   = lastindex(λv)
    λi  = λv[l]
    μi  = lμ(tree)[l]
    ll += log(exp(λi) + exp(μi)) - λi
  end
  tna -= et

  if isfix(tree.d1::iTgbmbd)
    tnx = treeheight(tree.d2::iTgbmbd)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem(tree.d1::iTgbmbd, tna, ll)
  else
    tnx = treeheight(tree.d1::iTgbmbd)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem(tree.d2::iTgbmbd, tna, ll)
  end
end




"""
    sum_alone_stem_p(tree::iTgbmbd, 
                     tna ::Float64, 
                     ll  ::Float64)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
function sum_alone_stem_p(tree::iTgbmbd, 
                          tna ::Float64, 
                          ll  ::Float64)

  et = e(tree)
  if tna < et
    λv  = lλ(tree)
    l   = lastindex(λv)
    λi  = λv[l]
    μi  = lμ(tree)[l]
    ll += log(exp(λi) + exp(μi)) - λi
  end
  tna -= et

  if istip(tree)
    return ll
  end

  if isfix(tree.d1::iTgbmbd)
    tnx = treeheight(tree.d2::iTgbmbd)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem_p(tree.d1::iTgbmbd, tna, ll)
  else
    tnx = treeheight(tree.d1::iTgbmbd)
    if tnx > tna
      tna = tnx
    end
    sum_alone_stem_p(tree.d2::iTgbmbd, tna, ll)
  end
end
