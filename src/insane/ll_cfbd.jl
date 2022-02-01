#=

fossilized birth-death likelihoods

Jérémy Andréoletti
Adapted from birth-death likelihoods by Ignacio Quintero Mächler

v(^-^v)

Created 16 12 2021
=#




"""
    stem_prob_surv_cfbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
fossilized birth-death process with `λ` and `μ` from stem age.
"""
function stem_prob_surv_cfbd(λ::Float64, μ::Float64, t::Float64)
  @fastmath begin
    μ += λ === μ ? 1e-14 : 0.0
    - log((λ - μ)/(λ - μ*exp(-(λ - μ)*t)))
  end
end




"""
    crown_prob_surv_cfbd(λ::Float64, μ::Float64, t::Float64)

Log-probability of at least one lineage surviving after time `t` for 
fossilized birth-death process with `λ` and `μ` from crown age.
"""
function crown_prob_surv_cfbd(λ::Float64, μ::Float64, t::Float64)
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
    # fossil tips are labelled extinct despite no actual extinction event
    if     isfossil(tree)  return - e(tree)*(λ + μ + ψ) + log(ψ)
    elseif isextinct(tree) return - e(tree)*(λ + μ + ψ) + log(μ)
    else                   return - e(tree)*(λ + μ + ψ) 
    end
  
  elseif issampledancestor(tree)
    return - e(tree)*(λ + μ + ψ) + log(ψ) +
             (isdefined(tree, :d1) ? llik_cfbd(tree.d1::sTfbd, λ, μ, ψ) : 0.0) + 
             (isdefined(tree, :d2) ? llik_cfbd(tree.d2::sTfbd, λ, μ, ψ) : 0.0)
  
  else
    return - e(tree)*(λ + μ + ψ) + log(λ) +
             llik_cfbd(tree.d1::sTfbd, λ, μ, ψ) + 
             llik_cfbd(tree.d2::sTfbd, λ, μ, ψ)
  end
end




"""
    llik_cfbd(Ξ::Vector{sTfbd}, 
              λ::Float64, 
              μ::Float64,
              ψ::Float64)

Log-likelihood up to a constant for constant fossilized birth-death 
given a complete `iTree` for decoupled trees.
"""
function llik_cfbd(Ξ::Vector{sTfbd}, 
                   λ::Float64, 
                   μ::Float64,
                   ψ::Float64)

  ll = 0.0
  nsa = 0.0
  for ξ in Ξ
    if !isextinct(ξ)
      # number of sampled ancestors (fossils tips are labelled extinct)
      nsa += isfossil(ξ)
    end
    ll += llik_cfbd(ξ, λ, μ, ψ)
  end
  
  ll += Float64(lastindex(Ξ) - nsa - 1)/2.0 * log(λ)

  return ll
end




"""
    make_scond(idf::Vector{iBfffs}, stem::Bool, ::Type{sTfbd})

Return closure for log-likelihood for conditioning
"""
function make_scond(idf::Vector{iBfffs}, stem::Bool, ::Type{sTfbd})

  if stem
    # for whole likelihood
    f = (λ::Float64, μ::Float64, sns::NTuple{3,BitVector}) ->
          cond_ll(λ, μ, sns[1])
    # for new proposal
    f0 = (tree::sTfbd, λ::Float64, μ::Float64, ter::Bool) -> 
            cond_surv_stem_p(tree, λ, μ)
  else
    # for whole likelihood
    f = (λ::Float64, μ::Float64, sns::NTuple{3,BitVector}) ->
          cond_ll(λ, μ, sns[2]) + 
          cond_ll(λ, μ, sns[3]) + log((λ + μ)/λ)
    # for new proposal
    f0 = function (tree::sTfbd, λ::Float64, μ::Float64, ter::Bool)
      if ter
        cond_surv_stem(  tree, λ, μ)
      else
        cond_surv_stem_p(tree, λ, μ)
      end
    end
  end

  return f, f0
end



