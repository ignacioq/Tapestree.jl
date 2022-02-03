#=

fossilized birth-death likelihoods

Jérémy Andréoletti
Adapted from birth-death likelihoods by Ignacio Quintero Mächler

v(^-^v)

Created 11 10 2021
=#



#="""
    cond_surv_crown(tree::sTfbd, λ::Float64, μ::Float64)

Log-probability of at least two lineage surviving for fossilized birth-death 
process with `λ` and `μ` for crown age.
"""
function cond_surv_crown(tree::sTfbd, λ::Float64, μ::Float64)
  survdr1 = survivaldr(tree.d1::sTfbd)
  survdr2 = survivaldr(tree.d2::sTfbd)
  ldr1 = lastindex(survdr1)
  ldr2 = lastindex(survdr2)
  n = sum_alone_stem(treermfossils.d1::sTfbd, 0.0, 0.0, survdr1, ldr2, 0) +
      sum_alone_stem(treermfossils.d2::sTfbd, 0.0, 0.0, survdr1, ldr2, 0)
  return n*log((λ + μ)/λ)
end




"""
    cond_nothing(tree::sTfbd, λ::Float64, μ::Float64)

No conditioning when computing tree likelihood.
"""
cond_nothing(tree::sTfbd, λ::Float64, μ::Float64) = 0




"""
    stem_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
fossilized birth-death process with `λ` and `μ` from stem age.
"""
function stem_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)
  @fastmath begin
    μ += λ === μ ? 1e-14 : 0.0
    - log((λ - μ)/(λ - μ*exp(-(λ - μ)*t)))
  end
end




"""
    crown_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
fossilized birth-death process with `λ` and `μ` from crown age.
"""
function crown_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)
    μ += λ === μ ? 1e-14 : 0.0
    - 2.0 * log((λ - μ)/(λ - μ*exp(-(λ - μ)*t))) - log(λ)
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
                  ll   ::Float64, 
                  λc   ::Float64, 
                  μc   ::Float64,
                  ψc   ::Float64,
                  dri  ::BitArray{1}, 
                  ldr  ::Int64,
                  wpr  ::Int64,
                  ix   ::Int64, 
                  px   ::Int64)

Return the Log-likelihood under constant fossilized birth-death 
of a grafted subtree determined by `dri`. 
"""
function stree_ll_cfbd(tree::sTfbd,
                      ll   ::Float64, 
                      λc   ::Float64, 
                      μc   ::Float64,
                      ψc   ::Float64,
                      dri  ::BitArray{1}, 
                      ldr  ::Int64,
                      wpr  ::Int64,
                      ix   ::Int64, 
                      px   ::Int64)

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



=#