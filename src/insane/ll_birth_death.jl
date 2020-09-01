#=

birth-death likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    stem_prob_surv(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
birth-death process with `λ` and `μ` from stem age.
"""
function stem_prob_surv(λ::Float64, μ::Float64, t::Float64)
  μ += λ == μ ? 1e-14 : 0.0
  log(1.0 - (μ - μ*exp(-(λ - μ)*t))/(λ - μ * exp(-(λ - μ)*t)))
end




"""
    crown_prob_surv(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
birth-death process with `λ` and `μ` from stem age.
"""
function crown_prob_surv(λ::Float64, μ::Float64, t::Float64)
  μ += λ == μ ? 1e-14 : 0.0
  2.0*log(1.0 - (μ - μ*exp(-(λ - μ)*t))/(λ - μ * exp(-(λ - μ)*t))) +
  log(2.0*λ)
end




"""
    llik_cbd(tree::iTree, λ::Float64, μ::Float64)

Log-likelihood up to a constant for constant birth-death 
given a complete `iTree` recursively.
"""
function llik_cbd(tree::iTree, λ::Float64, μ::Float64)
  if istip(tree) 
    - pe(tree)*(λ + μ) + (isextinct(tree) ? log(μ) : 0.0)
  else
    log(2.0*λ) - pe(tree)*(λ + μ) +
    llik_cbd(tree.d1::iTree, λ, μ) + 
    llik_cbd(tree.d2::iTree, λ, μ)
  end
end




"""
    stree_ll_cbd(tree::iTree,
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
function stree_ll_cbd(tree::iTree,
                      ll  ::Float64, 
                      λc  ::Float64, 
                      μc  ::Float64,
                      dri ::BitArray{1}, 
                      ldr ::Int64,
                      wpr ::Int64,
                      ix  ::Int64, 
                      px  ::Int64)

  if ix == ldr
    if px == wpr
      if isfix(tree.d1::iTree)
        ll += llik_cbd(tree.d2::iTree, λc, μc)
      elseif isfix(tree.d2::iTree)
        ll += llik_cbd(tree.d1::iTree, λc, μc)
      end
    else
      px += 1
      if isfix(tree.d1::iTree)
        ll += 
          stree_ll_cbd(tree.d1::iTree, ll, λc, μc, dri, ldr, wpr, ix, px)
      else
        ll +=
          stree_ll_cbd(tree.d2::iTree, ll, λc, μc, dri, ldr, wpr, ix, px)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1::iTree)
    if ifx1 && isfix(tree.d2::iTree)
      ix += 1
      if dri[ix]
        ll += 
          stree_ll_cbd(tree.d1::iTree, ll, λc, μc, dri, ldr, wpr, ix, px)
      else
        ll += 
          stree_ll_cbd(tree.d2::iTree, ll, λc, μc, dri, ldr, wpr, ix, px)
      end
    elseif ifx1
      ll += 
        stree_ll_cbd(tree.d1::iTree, ll, λc, μc, dri, ldr, wpr, ix, px)
    else
      ll +=
        stree_ll_cbd(tree.d2::iTree, ll, λc, μc, dri, ldr, wpr, ix, px)
    end
  end

  return ll
end



