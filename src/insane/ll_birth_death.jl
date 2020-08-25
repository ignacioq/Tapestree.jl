#=

birth-death likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




"""
    llik_cbd(tree::iTree, λ::Float64, μ::Float64)

Log-likelihood up to a constant for constant birth-death 
given a complete `iTree`.
"""
function llik_cbd(tree::iTree, λ::Float64, μ::Float64)

  if istip(tree) 
    - pe(tree)*(λ + μ) + (isextinct(tree) ? log(μ) : 0.0)
  else
    llik_cbd(tree.d1, λ, μ) + 
    llik_cbd(tree.d2, λ, μ) + 
    log(λ) - pe(tree)*(λ + μ)
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
      if isfix(tree.d1)
        ll += llik_cbd(tree.d2, λc, μc)
      elseif isfix(tree.d2)
        ll += llik_cbd(tree.d1, λc, μc)
      end
    else
      px += 1
      if isfix(tree.d1)
        ll += 
          stree_ll_cbd(tree.d1, ll, λc, μc, dri, ldr, wpr, ix, px)
      else
        ll +=
          stree_ll_cbd(tree.d2, ll, λc, μc, dri, ldr, wpr, ix, px)
      end
    end
  elseif ix < ldr
    ifx1 = isfix(tree.d1)
    if ifx1 && isfix(tree.d2)
      ix += 1
      if dri[ix]
        ll += 
          stree_ll_cbd(tree.d1, ll, λc, μc, dri, ldr, wpr, ix, px)
      else
        ll += 
          stree_ll_cbd(tree.d2, ll, λc, μc, dri, ldr, wpr, ix, px)
      end
    elseif ifx1
      ll += 
        stree_ll_cbd(tree.d1, ll, λc, μc, dri, ldr, wpr, ix, px)
    else
      ll +=
        stree_ll_cbd(tree.d2, ll, λc, μc, dri, ldr, wpr, ix, px)
    end
  end

  return ll
end



