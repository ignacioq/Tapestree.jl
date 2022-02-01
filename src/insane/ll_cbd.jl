#=

birth-death likelihoods

Ignacio Quintero Mächler

t(-_-t)

Created 06 07 2020
=#




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
    - log((λ - μ)/(λ - μ*exp(-(λ - μ)*t)))
  end
end




"""
    crown_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
birth-death process with `λ` and `μ` from stem age.
"""
function crown_prob_surv_cbd(λ::Float64, μ::Float64, t::Float64)
    μ += λ === μ ? 1e-14 : 0.0
    - 2.0 * log((λ - μ)/(λ - μ*exp(-(λ - μ)*t))) - log(λ)
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
    llik_cbd(xi::Vector{sTbd}, 
             λ  ::Float64, 
             μ  ::Float64)

Log-likelihood up to a constant for constant birth-death 
given a complete `iTree` for decoupled trees.
"""
function llik_cbd(xi::Vector{sTbd}, 
                  λ  ::Float64, 
                  μ  ::Float64)

  ll = 0.0
  for ξ in xi
    ll += llik_cbd(ξ, λ, μ)
  end

  ll += Float64(lastindex(xi) - 1)/2.0 * log(λ)

  return ll
end




"""
    make_scond(idf::Vector{iBffs}, stem::Bool, ::Type{sTbd})

Return closure for log-likelihood for conditioning
"""
function make_scond(idf::Vector{iBffs}, stem::Bool, ::Type{sTbd})

  if stem
    # for whole likelihood
    f = (λ::Float64, μ::Float64, sns::NTuple{3,BitVector}) ->
          cond_ll(λ, μ, sns[1])
    # for new proposal
    f0 = (xi::sTbd, λ::Float64, μ::Float64, ter::Bool) -> 
            cond_surv_stem_p(xi, λ, μ)
  else
    # for whole likelihood
    f = (λ::Float64, μ::Float64, sns::NTuple{3,BitVector}) ->
          cond_ll(λ, μ, sns[2]) + 
          cond_ll(λ, μ, sns[3]) + log((λ + μ)/λ)
    # for new proposal
    f0 = function (xi::sTbd, λ::Float64, μ::Float64, ter::Bool)
      if ter
        cond_surv_stem(  xi, λ, μ)
      else
        cond_surv_stem_p(xi, λ, μ)
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



