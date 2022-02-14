#=

occurrence birth-death likelihoods

Jérémy Andréoletti

v(^-^v)

Created 11 02 2022
=#




"""
    llik_cobd(tree  ::sTfbd,
              ωtimes::Vector{Float64},
              λ     ::Float64,
              μ     ::Float64,
              ψ     ::Float64,
              ω     ::Float64)

Log-likelihood up to a constant for constant occurrence birth-death 
given a complete `iTree` recursively.
"""
function llik_cobd(tree  ::sTfbd,
                   ωtimes::Vector{Float64},
                   λ     ::Float64,
                   μ     ::Float64,
                   ψ     ::Float64,
                   ω     ::Float64)
  
  se = Float64[]
  ee = Float64[]

  # Compute the FBD tree likelihood and get speciation/extinction event times
  ll = _llik_cfbd_eventimes!(tree::sTfbd, λ, μ, ψ, 0.0, se, ee)

  # Construct the LTT (Lineages Through Time)
  push!(se, 0.0)

  lse = lastindex(se)
  lee = lastindex(ee)

  events = append!(se, ee)
  jumps_order = sortperm(events)
  
  jumps = ones(Int64, lse)
  append!(jumps,  fill(-1, lee))
  jumps = jumps[jumps_order]

  LTTn = cumsum!(jumps,jumps)   # number of lineages at each LTT step
  LTTt = events[jumps_order]    # time of each LTT step
  Δt   = diff(LTTt)             # duration of each LTT step
  #LTT = Ltt(LTTn, LTTt)

  # Update the log-likelihood with the observed occurrences
  kω = 0
  i  = 1
  for ωtime in ωtimes
    # count the number of occurrences between each LTT step
    if ωtime < LTTt[i]
      kω += 1
    else
      ni = LTTn[i]
      Δti = Δt[i]
      # log-probability of observing kω occurrences (Poisson distribution)
      ll += logpdf(Poisson(ni*ω*Δti), kω)
      
      #= ### much slower method when kω>20 ###
      Eω = ni*ω*Δti  # expected number of occurrences during the LTT step i
      ll -= Eω
      if kω != 0
        ll += kω*log(Eω) - log(factorial(kω))
        kω = 0
      end=#

      i += 1
    end
  end

  return ll
end




"""
    llik_cobd(tree  ::sTfbd,
              ωtimes::Vector{Float64},
              LTT   ::Ltt, 
              λ     ::Float64,
              μ     ::Float64,
              ψ     ::Float64,
              ω     ::Float64)

Log-likelihood up to a constant for constant occurrence birth-death 
given a complete `iTree` recursively.
"""
function llik_cobd(tree  ::sTfbd,
                   ωtimes::Vector{Float64},
                   LTT   ::Ltt,
                   λ     ::Float64,
                   μ     ::Float64,
                   ψ     ::Float64,
                   ω     ::Float64)
  
  # Compute the FBD tree likelihood
  ll = llik_cfbd(tree::sTfbd, λ, μ, ψ)

  # Update the log-likelihood with the observed occurrences
  Δt = diff(LTT.t)  # duration of each LTT step
  kω = 0
  i  = 1
  for ωtime in ωtimes
    # count the number of occurrences between each LTT step
    if ωtime < LTT.t[i]
      kω += 1
    else
      ni = LTT.n[i]
      Δti = Δt[i]
      # log-probability of observing kω occurrences (Poisson distribution)
      ll += logpdf(Poisson(ni*ω*Δti), kω)
      i += 1
    end
  end

  return ll
end




"""
    _llik_cfbd_eventimes!(tree  ::sTfbd,
                          λ     ::Float64,
                          μ     ::Float64,
                          ψ     ::Float64,
                          t     ::Float64, 
                          se    ::Vector{Float64}, 
                          ee    ::Vector{Float64})

Log-likelihood up to a constant for constant occurrence birth-death 
given a complete `iTree` recursively.
"""
function _llik_cfbd_eventimes!(tree  ::sTfbd,
                               λ     ::Float64,
                               μ     ::Float64,
                               ψ     ::Float64,
                               t     ::Float64, 
                               se    ::Vector{Float64}, 
                               ee    ::Vector{Float64})
  
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  
  et = e(tree)

  # tip
  if !defd1 && !defd2
    # fossil tips are labelled extinct despite no actual extinction event
    if     isfossil(tree)
      return - e(tree)*(λ + μ + ψ) + log(ψ)
    elseif isextinct(tree)
      push!(ee, t + et)
      return - e(tree)*(λ + μ + ψ) + log(μ)
    else                   
      return - e(tree)*(λ + μ + ψ) 
    end
  
  # sampled ancestor
  elseif defd1 ⊻ defd2
    return - e(tree)*(λ + μ + ψ) + log(ψ) +
           (defd1 ? _llik_cfbd_eventimes!(tree.d1::sTfbd, λ, μ, ψ, t + et, se, ee) : 0.0) + 
           (defd2 ? _llik_cfbd_eventimes!(tree.d2::sTfbd, λ, μ, ψ, t + et, se, ee) : 0.0)
  
  # bifurcation
  else
    push!(se, t + et)
    return - e(tree)*(λ + μ + ψ) + log(λ) +
             _llik_cfbd_eventimes!(tree.d1::sTfbd, λ, μ, ψ, t + et, se, ee) + 
             _llik_cfbd_eventimes!(tree.d2::sTfbd, λ, μ, ψ, t + et, se, ee)
  end
end




"""
    llik_cobd(Ξ::Vector{sTfbd}, 
              λ::Float64, 
              μ::Float64,
              ψ::Float64)

Log-likelihood up to a constant for constant occurrence birth-death 
given a complete `iTree` for decoupled trees.
"""
function llik_cobd(Ξ::Vector{sTfbd}, 
                   λ::Float64, 
                   μ::Float64,
                   ψ::Float64)

  ll = 0.0
  nsa = 0.0    # number of sampled ancestors
  for ξ in Ξ
    nsa += issampledancestor(ξ)
    ll += llik_cobd(ξ, λ, μ, ψ, ω)
  end
  
  ll += Float64(lastindex(Ξ) - nsa - 1)/2.0 * log(λ)

  return ll
end



