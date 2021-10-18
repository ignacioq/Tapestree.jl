#=

fossilized birth-death likelihoods

Jérémy Andréoletti
Adapted from birth-death likelihoods by Ignacio Quintero Mächler

v(^-^v)

Created 11 10 2021
=#



"""
    cond_surv_crown(tree::sTfbd, λ::Float64, μ::Float64)

Log-probability of at least two lineage surviving for fossilized birth-death 
process with `λ` and `μ` for crown age.
"""
function cond_surv_crown(tree::sTfbd, λ::Float64, μ::Float64)
  n = sum_alone_stem(tree.d1::sTfbd, 0.0, 0.0) +
      sum_alone_stem(tree.d2::sTfbd, 0.0, 0.0)

  return n*log((λ + μ)/λ)
end




"""
    cond_surv_stem(tree::sTfbd, λ::Float64, μ::Float64)

Log-probability of at least one lineage surviving for fossilized birth-death 
process with `λ` and `μ` for stem age.
"""
function cond_surv_stem(tree::sTfbd, λ::Float64, μ::Float64)
  n = sum_alone_stem(tree, 0.0, 0.0)
  return n*log((λ + μ)/λ)
end




"""
    sum_alone_stem(tree::sTfbd, tna::Float64, n::Float64)

Count nodes in stem lineage when a diversification event could have 
returned an overall extinction.
"""
function sum_alone_stem(tree::sTfbd, tna::Float64, n::Float64)

  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  # Tip
  if !defd1 && !defd2
    return n
  end

  # Sampled ancestors
  if (defd1 && !defd2) 
    sum_alone_stem(tree.d1::sTfbd, max(0,tna-e(tree)), n)
  elseif (defd2 && !defd1) 
    sum_alone_stem(tree.d2::sTfbd, max(0,tna-e(tree)), n)
  
  else
    if tna < e(tree)
      n += 1.0
    end
    tna -= e(tree)

    if isfix(tree.d1::sTfbd)
      if isfix(tree.d2::sTfbd)
        # Birth of 2 fixed daughters: extinction is now impossible
        return n
      else
        tnx = treeheight(tree.d2::sTfbd)
        tna = tnx > tna ? tnx : tna
        sum_alone_stem(tree.d1::sTfbd, tna, n)
      end
    else
      tnx = treeheight(tree.d1::sTfbd)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem(tree.d2::sTfbd, tna, n)
    end
  end

end




"""
    cond_surv_stem_p(tree::sTfbd, λ::Float64, μ::Float64)

Log-probability of at least one lineage surviving after time `t` for 
fossilized birth-death process with `λ` and `μ` for stem age.
"""
function cond_surv_stem_p(tree::sTfbd, λ::Float64, μ::Float64)
  n = sum_alone_stem_p(tree, 0.0, 0.0)
  return n*log((λ + μ)/λ)
end




"""
    sum_alone_stem(tree::sTfbd, tna::Float64, n::Float64)

Count nodes in stem lineage when a diversification event could have 
returned an overall extinction.
"""
function sum_alone_stem_p(tree::sTfbd, tna::Float64, n::Float64)

  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if tna < e(tree)
    n += 1.0
  end
  tna -= e(tree)

  # Tip
  if !defd1 && !defd2
     return n
  end

  # Sampled ancestors
  if (defd1 && !defd2) 
    sum_alone_stem(tree.d1::sTfbd, max(0,tna-e(tree)), n)
  elseif (defd2 && !defd1) 
    sum_alone_stem(tree.d2::sTfbd, max(0,tna-e(tree)), n)

  else
    if isfix(tree.d1::sTfbd)
      if isfix(tree.d2::sTfbd)
        # Birth of 2 fixed daughters: extinction is now impossible
        return n
      else
        tnx = treeheight(tree.d2::sTfbd)
        tna = tnx > tna ? tnx : tna
        sum_alone_stem_p(tree.d1::sTfbd, tna, n)
      end
    else
      tnx = treeheight(tree.d1::sTfbd)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem_p(tree.d2::sTfbd, tna, n)
    end
  end

end




"""
    stem_prob_surv_cfbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
fossilized birth-death process with `λ` and `μ` from stem age.
"""
function stem_prob_surv_cfbd(λ::Float64, μ::Float64, t::Float64)
  @fastmath begin
    μ += λ === μ ? 1e-14 : 0.0
    log((λ - μ)/(λ - μ*exp(-(λ - μ)*t)))
  end
end




"""
    crown_prob_surv_cfbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
fossilized birth-death process with `λ` and `μ` from stem age.
"""
function crown_prob_surv_cfbd(λ::Float64, μ::Float64, t::Float64)
    μ += λ === μ ? 1e-14 : 0.0
    log(((λ - μ)/(λ - μ*exp(-(λ - μ)*t)))^2)
end




"""
    llik_cfbd(tree::sTfbd, λ::Float64, μ::Float64, ψ::Float64)

Log-likelihood up to a constant for constant fossilized birth-death 
given a complete `iTree` recursively.
"""
function llik_cfbd(tree::sTfbd, λ::Float64, μ::Float64, ψ::Float64)
  if istip(tree)
    - e(tree)*(λ + μ + ψ) + 
        (isextinct(tree) ? log(μ) : 0.0) + 
        (isfossil( tree) ? log(ψ) : 0.0)
  
  elseif issampledancestor(tree)
    - e(tree)*(λ + μ + ψ) + log(ψ) +
        (isdefined(tree, :d1) ? llik_cfbd(tree.d1::sTfbd, λ, μ, ψ) : 0.0) + 
        (isdefined(tree, :d2) ? llik_cfbd(tree.d2::sTfbd, λ, μ, ψ) : 0.0)
  else
    - e(tree)*(λ + μ + ψ) + log(λ) +
        llik_cfbd(tree.d1::sTfbd, λ, μ, ψ) + 
        llik_cfbd(tree.d2::sTfbd, λ, μ, ψ)
  end
end




"""
     stree_ll_cfbd(tree::sTfbd,
                  ll  ::Float64, 
                  λc  ::Float64, 
                  μc  ::Float64,
                  ψc  ::Float64,
                  dri ::BitArray{1}, 
                  ldr ::Int64,
                  wpr ::Int64,
                  ix  ::Int64, 
                  px  ::Int64)

Return the Log-likelihood under constant fossilized birth-death 
of a grafted subtree determined by `dri`. 
"""
function stree_ll_cfbd(tree::sTfbd,
                      ll  ::Float64, 
                      λc  ::Float64, 
                      μc  ::Float64,
                      ψc  ::Float64,
                      dri ::BitArray{1}, 
                      ldr ::Int64,
                      wpr ::Int64,
                      ix  ::Int64, 
                      px  ::Int64)

  if ix === ldr
    if px === wpr
      if isfix(tree.d1::sTfbd)
        ll += llik_cfbd(tree.d2::sTfbd, λc, μc, ψc)
      elseif isfix(tree.d2::sTfbd)
        ll += llik_cfbd(tree.d1::sTfbd, λc, μc, ψc)
      end
    else
      px += 1
      if isfix(tree.d1::sTfbd)
        ll += 
          stree_ll_cfbd(tree.d1::sTfbd, ll, λc, μc, ψc, dri, ldr, wpr, ix, px)
      else
        ll +=
          stree_ll_cfbd(tree.d2::sTfbd, ll, λc, μc, ψc, dri, ldr, wpr, ix, px)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::sTfbd)
    if ifx1 && isfix(tree.d2::sTfbd)
      ix += 1
      if dri[ix]
        ll += 
          stree_ll_cfbd(tree.d1::sTfbd, ll, λc, μc, ψc, dri, ldr, wpr, ix, px)
      else
        ll += 
          stree_ll_cfbd(tree.d2::sTfbd, ll, λc, μc, ψc, dri, ldr, wpr, ix, px)
      end
    elseif ifx1
      ll += 
        stree_ll_cfbd(tree.d1::sTfbd, ll, λc, μc, ψc, dri, ldr, wpr, ix, px)
    else
      ll +=
        stree_ll_cfbd(tree.d2::sTfbd, ll, λc, μc, ψc, dri, ldr, wpr, ix, px)
    end
  end

  return ll
end



