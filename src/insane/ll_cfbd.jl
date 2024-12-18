#=

fossilized birth-death likelihoods with epochs

Jérémy Andréoletti and Ignacio Quintero Mächler

v(^-^v)
t(-_-t)

Created 16 12 2021
=#




"""
    llik_cfbd(Ξ  ::Vector{sTfbd},
              λ  ::Float64,
              μ  ::Float64,
              ψ  ::Vector{Float64},
              nλ ::Float64,
              ψts::Vector{Float64},
              bst::Vector{Float64},
              eix::Vector{Int64})

Log-likelihood up to a constant for constant fossilized birth-death
given a complete `iTree` for decoupled trees with epochs `ψts`.
"""
function llik_cfbd(Ξ  ::Vector{sTfbd},
                   λ  ::Float64,
                   μ  ::Float64,
                   ψ  ::Vector{Float64},
                   nλ ::Float64,
                   ψts::Vector{Float64},
                   bst::Vector{Float64},
                   eix::Vector{Int64})
  @inbounds begin

    nep = lastindex(ψts) + 1
    ll  = 0.0
    for i in Base.OneTo(lastindex(Ξ))
      ξ   = Ξ[i]
      ll += llik_cfbd(ξ, λ, μ, ψ, bst[i], ψts, eix[i], nep)
    end

    ll += nλ * log(λ)
  end

  return ll
end




"""
    llik_cfbd(tree::sTfbd,
              λ   ::Float64,
              μ   ::Float64,
              ψ   ::Vector{Float64},
              t   ::Float64,
              ψts ::Vector{Float64},
              ix  ::Int64,
              nep ::Int64)

Log-likelihood up to a constant for constant fossilized birth-death
given a complete `iTree` recursively.
"""
function llik_cfbd(tree::sTfbd,
                   λ   ::Float64,
                   μ   ::Float64,
                   ψ   ::Vector{Float64},
                   t   ::Float64,
                   ψts ::Vector{Float64},
                   ix  ::Int64,
                   nep ::Int64)
  @inbounds begin

    ei  = e(tree)
    ll = 0.0
    # if epoch change
    while ix < nep && t - ei < ψts[ix]
      li   = t - ψts[ix]
      ll  -= li*(λ + μ + ψ[ix])
      ei  -= li
      t    = ψts[ix]
      ix  += 1
    end

    ll -= ei*(λ + μ + ψ[ix])
    t  -= ei

    if def1(tree)
      if def2(tree)
        ll += log(λ)                                              +
              llik_cfbd(tree.d1::sTfbd, λ, μ, ψ, t, ψts, ix, nep) +
              llik_cfbd(tree.d2::sTfbd, λ, μ, ψ, t, ψts, ix, nep)
      else
        ll += log(ψ[ix])                                          +
              llik_cfbd(tree.d1::sTfbd, λ, μ, ψ, t, ψts, ix, nep)
      end
    else
      ll += (isextinct(tree) ? log(μ)     : 0.0) +
            (isfossil(tree)  ? log(ψ[ix]) : 0.0)
    end
  end

  return ll
end



