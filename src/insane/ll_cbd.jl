#=

birth-death likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#



"""
    crown_prob_surv_da(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least two lineage surviving after time `t` for 
birth-death process with `λ` and `μ` for crown age.
"""
function cond_surv_crown(tree::sTbd, λ::Float64, μ::Float64)
  n = sum_alone_stem(tree.d1::sTbd, 0.0, 0.0) +
      sum_alone_stem(tree.d2::sTbd, 0.0, 0.0)

  return n*log((λ + μ)/λ)
end




"""
    stem_prob_surv_da(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
birth-death process with `λ` and `μ` for stem age.
"""
function cond_surv_stem(tree::sTbd, λ::Float64, μ::Float64)
  n = sum_alone_stem(tree, 0.0, 0.0)
  return n*log((λ + μ)/λ)
end




"""
    sum_alone_stem(tree::sTbd, tna::Float64, n::Float64, t::Float64)

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
    if isfix(tree.d2::sTbd)
      return n
    else
      tnx = treeheight(tree.d2::sTbd, 0.0, 0.0)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem(tree.d1::sTbd, tna, n)
    end
  else
    tnx = treeheight(tree.d1::sTbd, 0.0, 0.0)
    tna = tnx > tna ? tnx : tna
    sum_alone_stem(tree.d2::sTbd, tna, n)
  end

end




"""
    stem_prob_surv_da(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
birth-death process with `λ` and `μ` for stem age.
"""
function cond_surv_stem_p(tree::sTbd, λ::Float64, μ::Float64)
  n = sum_alone_stem_p(tree, 0.0, 0.0)
  return n*log((λ + μ)/λ)
end




"""
    sum_alone_stem(tree::sTbd, tna::Float64, n::Float64, t::Float64)

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
    if isfix(tree.d2::sTbd)
      return n
    else
      tnx = treeheight(tree.d2::sTbd, 0.0, 0.0)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem_p(tree.d1::sTbd, tna, n)
    end
  else
    tnx = treeheight(tree.d1::sTbd, 0.0, 0.0)
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
     stree_ll_cbd(tree::sTbd,
                  ll  ::Float64, 
                  λc  ::Float64, 
                  μc  ::Float64,
                  dri ::BitArray{1}, 
                  ldr ::Int64,
                  wpr ::Int64,
                  ix  ::Int64, 
                  px  ::Int64)

Return the Log-likelihood under constant birth-death 
of a grafted subtree determined by `dri`. 
"""
function stree_ll_cbd(tree::sTbd,
                      ll  ::Float64, 
                      λc  ::Float64, 
                      μc  ::Float64,
                      dri ::BitArray{1}, 
                      ldr ::Int64,
                      wpr ::Int64,
                      ix  ::Int64, 
                      px  ::Int64)

  if ix === ldr
    if px === wpr
      if isfix(tree.d1::sTbd)
        ll += llik_cbd(tree.d2::sTbd, λc, μc)
      elseif isfix(tree.d2::sTbd)
        ll += llik_cbd(tree.d1::sTbd, λc, μc)
      end
    else
      px += 1
      if isfix(tree.d1::sTbd)
        ll += 
          stree_ll_cbd(tree.d1::sTbd, ll, λc, μc, dri, ldr, wpr, ix, px)
      else
        ll +=
          stree_ll_cbd(tree.d2::sTbd, ll, λc, μc, dri, ldr, wpr, ix, px)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::sTbd)
    if ifx1 && isfix(tree.d2::sTbd)
      ix += 1
      if dri[ix]
        ll += 
          stree_ll_cbd(tree.d1::sTbd, ll, λc, μc, dri, ldr, wpr, ix, px)
      else
        ll += 
          stree_ll_cbd(tree.d2::sTbd, ll, λc, μc, dri, ldr, wpr, ix, px)
      end
    elseif ifx1
      ll += 
        stree_ll_cbd(tree.d1::sTbd, ll, λc, μc, dri, ldr, wpr, ix, px)
    else
      ll +=
        stree_ll_cbd(tree.d2::sTbd, ll, λc, μc, dri, ldr, wpr, ix, px)
    end
  end

  return ll
end



