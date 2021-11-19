#=

birth-death likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#



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
    moments(n::Float64, th::Float64, ϵ::Float64)

Returns `λ` and `μ` under the method of moments given `n` surviving tips,
treeheight `th` and some turnover value `ϵ`.
"""
function moments(n::Float64, th::Float64, ϵ::Float64)
  δ  = 1.0/th * log(n*(1.0 - ϵ) + ϵ)
  λc = δ/(1.0 - ϵ)
  μc = λc - δ
  return λc, μc
end




"""
    stem_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
birth-death process with `λ` and `μ` from stem age.
"""
function stem_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)
  @fastmath begin
    μ += λ === μ ? 1e-14 : 0.0
    log((λ - μ)/(λ - μ*exp(-(λ - μ)*t)))
  end
end




"""
    crown_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
birth-death process with `λ` and `μ` from stem age.
"""
function crown_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)
    μ += λ === μ ? 1e-14 : 0.0
    log(((λ - μ)/(λ - μ*exp(-(λ - μ)*t)))^2)
end




"""
    llik_cbd(tree::sTbd, λ::Float64, μ::Float64)

Log-likelihood up to a constant for constant birth-death 
given a complete `iTree` recursively.
"""
function llik_cbd(tree::sTbd, λ::Float64, μ::Float64)
  if istip(tree) 
    - e(tree)*(λ + μ) + (isextinct(tree) ? log(μ) : 0.0)
  else
    log(λ) - e(tree)*(λ + μ)      +
    llik_cbd(tree.d1::sTbd, λ, μ) +
    llik_cbd(tree.d2::sTbd, λ, μ)
  end
end




"""
    llik_cbd(psi::Vector{sTbd}, 
             idf::Vector{iBffs},
             λ   ::Float64, 
             μ   ::Float64)

Log-likelihood up to a constant for constant birth-death 
given a complete `iTree` for edge trees.
"""
function llik_cbd(psi ::Vector{sTbd}, 
                  idf ::Vector{iBffs},
                  λ    ::Float64, 
                  μ    ::Float64)

  ll = 0.0
  for ψ in psi
    ll += llik_cbd(ψ, λ, μ)
  end

  ll += Float64(lastindex(psi) - 1)/2.0 * log(λ)

  return ll
end




"""
    make_scond(idf::Vector{iBffs}, stem::Bool, ::Type{sTbd})

Return closure for log-likelihood for conditioning
"""
function make_scond(idf::Vector{iBffs}, stem::Bool, ::Type{sTbd})

  b1  = idf[1]
  d1i = d1(b1)
  d2i = d2(b1)

  if stem
    # for whole likelihood
    f = (λ::Float64, μ::Float64, sns::NTuple{3,BitVector}) ->
          cond_ll(λ, μ, sns[1])
    # for new proposal
    f0 = (psi::sTbd, λ::Float64, μ::Float64) -> cond_surv_stem_p(psi, λ, μ)
  else
    # for whole likelihood
    f = (λ::Float64, μ::Float64, sns::NTuple{3,BitVector}) ->
          cond_ll(λ, μ, sns[2]) + cond_ll(λ, μ, sns[3]) + log((λ + μ)/λ)
    # for new proposal
    f0 = function (psi::sTbd, λ::Float64, μ::Float64, ter::Bool)
      if ter
        cond_surv_stem(  psi, λ, μ)
      else
        cond_surv_stem_p(psi, λ, μ)
      end
    end
  end

  return f, f0
end




"""
    cond_ll(λ::Float64, μ::Float64, sn::BitVector)

Condition events when there is only one alive lineage in the crown subtrees 
to only be speciation events.
"""
cond_ll(λ::Float64, μ::Float64, sn::BitVector) =
  Float64(sum(sn)) * log((λ + μ)/λ)



