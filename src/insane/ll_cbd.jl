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
    if isfix(tree.d2::sTbd)
      return n
    else
      tnx = treeheight(tree.d2::sTbd)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem(tree.d1::sTbd, tna, n)
    end
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
    if isfix(tree.d2::sTbd)
      return n
    else
      tnx = treeheight(tree.d2::sTbd)
      tna = tnx > tna ? tnx : tna
      sum_alone_stem_p(tree.d1::sTbd, tna, n)
    end
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
             μ   ::Float64, 
             scond::Function)

Log-likelihood up to a constant for constant birth-death 
given a complete `iTree` for edge trees.
"""
function llik_cbd(psi::Vector{sTbd}, 
                  idf::Vector{iBffs},
                  λ   ::Float64, 
                  μ   ::Float64, 
                  scond::Function)

  ll = 0.0
  for ψ in psi
    ll += llik_cbd(ψ, λ, μ)
  end

  ll += Float64(lastindex(psi) - 1)/2.0 * log(λ) + scond(psi, λ, μ)

  return ll
end




"""
    make_cond(idf::Vector{iBffs}, stem::Bool)

Make closure for conditioning function
"""
function make_cond(idf::Vector{iBffs}, stem::Bool)

  # conditioning
  if stem
    function f(psi::Vector{sTbd}, λ::Float64, μ::Float64)
      cond_surv_stem_p(psi[1], λ, μ)
    end
  else
    b1  = idf[1]
    d1i = d1(b1)
    d2i = d2(b1)
    t1  = iszero(d1(idf[d1i]))
    t2  = iszero(d1(idf[d2i]))

    if t1 
      if t2
        f = let d1i = d1i, d2i = d2i
          (psi::Vector{sTbd}, λ::Float64, μ::Float64) ->
            cond_surv_stem(psi[d1i], λ, μ) + 
            cond_surv_stem(psi[d2i], λ, μ) + log((λ + μ)/λ)
        end
      else
        f = let d1i = d1i, d2i = d2i
          (psi::Vector{sTbd}, λ::Float64, μ::Float64) ->
          cond_surv_stem(  psi[d1i], λ, μ) + 
          cond_surv_stem_p(psi[d2i], λ, μ) + log((λ + μ)/λ)
        end 
      end
    elseif t2
      f = let d1i = d1i, d2i = d2i
        (psi::Vector{sTbd}, λ::Float64, μ::Float64) ->
        cond_surv_stem_p(psi[d1i], λ, μ) + 
        cond_surv_stem(  psi[d2i], λ, μ) + log((λ + μ)/λ)
      end
    else
      f = let d1i = d1i, d2i = d2i
        (psi::Vector{sTbd}, λ::Float64, μ::Float64) ->
        cond_surv_stem_p(psi[d1i], λ, μ) + 
        cond_surv_stem_p(psi[d2i], λ, μ) + log((λ + μ)/λ)
      end
    end
  end

  return f
end


