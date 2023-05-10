#=

occurrence birth-death likelihoods

Jérémy Andréoletti

v(^-^v)

Created 11 02 2022
=#




"""
    llik_cobd(Ξ     ::Vector{sTfbd}, 
              ωtimes::Vector{Float64},
              LTT   ::Ltt,
              λ     ::Float64, 
              μ     ::Float64,
              ψ     ::Vector{Float64},
              ω     ::Vector{Float64},
              nλ    ::Float64,
              tep   ::Vector{Float64},
              bst   ::Vector{Float64},
              eix   ::Vector{Int64})

Log-likelihood up to a constant for constant occurrence birth-death 
given a complete decoupled `iTree` and occurrence `ωtimes`.
"""
function llik_cobd(Ξ     ::Vector{sTfbd}, 
                   ωtimes::Vector{Float64},
                   LTT   ::Ltt,
                   λ     ::Float64, 
                   μ     ::Float64,
                   ψ     ::Vector{Float64},
                   ω     ::Vector{Float64},
                   nλ    ::Float64,
                   tep   ::Vector{Float64},
                   bst   ::Vector{Float64},
                   eix   ::Vector{Int64})
  @inbounds begin

    # Likelihood of the tree with fossils
    ll = llik_cfbd(Ξ, λ, μ, ψ, nλ, tep, bst, eix)
    
    # Likelihood of the observed occurrences
    ll += ω_llik(ωtimes, ω, tep, LTT)
  
  end
  
  return ll
end




"""
    llik_cobd(tree::sTfbd,
              ωtimes::Vector{Float64},
              λ::Float64,
              μ::Float64,
              ψ::Float64,
              ω::Float64)
    llik_cobd(tree::sTfbd,
              ωtimes::Vector{Float64},
              λ::Float64,
              μ::Float64,
              ψ::Vector{Float64},
              ω::Vector{Float64},
              tep::Vector{Float64})

Log-likelihood up to a constant for constant occurrence birth-death given a 
complete `iTree` recursively.
"""
llik_cobd(tree  ::sTfbd,
          ωtimes::Vector{Float64}, 
          λ     ::Float64, 
          μ     ::Float64, 
          ψ     ::Float64,
          ω     ::Float64) =
  llik_cobd(tree, ωtimes, λ, μ, [ψ], [ω], treeheight(tree), Float64[], 1, 1)
llik_cobd(tree  ::sTfbd,
          ωtimes::Vector{Float64}, 
          λ     ::Float64, 
          μ     ::Float64, 
          ψ     ::Vector{Float64}, 
          ω     ::Vector{Float64},
          tep   ::Vector{Float64}) =
  llik_cobd(tree, ωtimes, λ, μ, ψ, ω, treeheight(tree), tep, 1, lastindex(tep)+1)




"""
    llik_cobd(tree  ::sTfbd,
              ωtimes::Vector{Float64},
              λ     ::Float64,
              μ     ::Float64,
              ψ     ::Vector{Float64},
              ω     ::Vector{Float64},
              tor   ::Float64,
              tep   ::Vector{Float64},
              ix    ::Int64,
              nep   ::Int64)

Log-likelihood up to a constant for constant occurrence birth-death 
given a complete `iTree` recursively.
"""
function llik_cobd(tree  ::sTfbd,
                   ωtimes::Vector{Float64},
                   λ     ::Float64,
                   μ     ::Float64,
                   ψ     ::Vector{Float64},
                   ω     ::Vector{Float64},
                   tor   ::Float64,
                   tep   ::Vector{Float64},
                   ix    ::Int64,
                   nep   ::Int64)
  
  @inbounds begin
  
    se = Float64[]
    ee = Float64[]

    # compute the FBD tree likelihood and get speciation/extinction event times
    ll = _llik_cfbd_eventimes!(tree::sTfbd, λ, μ, ψ, tor, tep, ix, nep, se, ee)

    # construct the LTT (Lineages Through Time)
    push!(se, tor)
    
    lse = lastindex(se)
    lee = lastindex(ee)

    events = append!(se, ee)
    jumps_order = sortperm(events, rev=true)
    
    jumps = ones(Int64,lse)
    append!(jumps,  fill(-1,lee))
    jumps = jumps[jumps_order]

    LTTn = cumsum!(jumps,jumps)    # Number of lineages at each LTT step
    LTTt = events[jumps_order]     # Time of each LTT step

    # add a point at present
    push!(LTTt, 0.0)
    push!(LTTn, last(LTTn))

    LTT = Ltt(LTTn, LTTt)

    # update the log-likelihood with the observed occurrences
    ll += ω_llik(ωtimes, ω, tep, LTT)

  end

  return ll
end




"""
    llik_cobd(tree  ::sTfbd,
              ωtimes::Vector{Float64},
              LTT   ::Ltt,
              λ     ::Float64,
              μ     ::Float64,
              ψ     ::Vector{Float64},
              ω     ::Vector{Float64},
              t     ::Float64,
              tep   ::Vector{Float64},
              ix    ::Int64,
              nep   ::Int64)

Log-likelihood up to a constant for constant occurrence birth-death 
given a complete `iTree` recursively.
"""
function llik_cobd(tree  ::sTfbd,
                   ωtimes::Vector{Float64},
                   LTT   ::Ltt,
                   λ     ::Float64,
                   μ     ::Float64,
                   ψ     ::Vector{Float64},
                   ω     ::Vector{Float64},
                   t     ::Float64,
                   tep   ::Vector{Float64},
                   ix    ::Int64,
                   nep   ::Int64)
  
  @inbounds begin

    # Compute the FBD tree likelihood
    ll = llik_cfbd(tree, λ, μ, ψ, t, tep, ix, nep)

    # Update the log-likelihood with the observed occurrences
    ll += ω_llik(ωtimes, ω, tep, LTT)

  end

  return ll
end




"""
    _llik_cfbd_eventimes!(tree::sTfbd,
                          λ   ::Float64,
                          μ   ::Float64,
                          ψ   ::Vector{Float64},
                          t   ::Float64,
                          ψts ::Vector{Float64},
                          ix  ::Int64,
                          nep ::Int64, 
                          se  ::Vector{Float64}, 
                          ee  ::Vector{Float64})

Log-likelihood up to a constant for piecewise-constant fossilized birth-death
given a complete `iTree` recursively, while recording the times of speciations
and extinctions in order to construct the LTT.
"""
function _llik_cfbd_eventimes!(tree::sTfbd,
                               λ   ::Float64,
                               μ   ::Float64,
                               ψ   ::Vector{Float64},
                               t   ::Float64,
                               ψts ::Vector{Float64},
                               ix  ::Int64,
                               nep ::Int64, 
                               se  ::Vector{Float64}, 
                               ee  ::Vector{Float64})
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
        push!(se, t)
        ll += log(λ)                                              +
              _llik_cfbd_eventimes!(tree.d1::sTfbd, λ, μ, ψ, t, ψts, ix, nep, se, ee) +
              _llik_cfbd_eventimes!(tree.d2::sTfbd, λ, μ, ψ, t, ψts, ix, nep, se, ee)
      else
        ll += log(ψ[ix])                                          +
              _llik_cfbd_eventimes!(tree.d1::sTfbd, λ, μ, ψ, t, ψts, ix, nep, se, ee)
      end
    else
      if isextinct(tree)
        push!(ee, t)
        ll += log(μ)
      end
      ll += (isfossil(tree) ? log(ψ[ix]) : 0.0)
    end
  end

  return ll
end




"""
    ω_llik(ωtimes::Vector{Float64},
           ω     ::Vector{Float64},
           tep   ::Vector{Float64},
           LTT   ::Ltt)

Log-likelihood of the observed occurrences for a given `LTT` 
(Lineages Through Time) trajectory.
"""
function ω_llik(ωtimes::Vector{Float64},
                ω     ::Vector{Float64},
                tep   ::Vector{Float64},
                LTT   ::Ltt)
  
  @inbounds begin
  
    ll = 0.0
    
    lω = lastindex(ωtimes) # total number of occurrences
    kω = 0                 # number of occurrences per interval
    iω = 1                 # occurrence index
    
    ep   = 1               # epoch
    epti = LTT.t[1]        # beggining of current epoch
    eptf = ifelse(isempty(tep), 0.0, tep[ep]) # end of current epoch
    epc  = 0               # whether an iteration is an epoch change

    for (iLTT, ni) in enumerate(LTT.n[1:(end-1)])
      # count the number of occurrences between each LTT step
      while iω <= lω && ωtimes[iω] > LTT.t[iLTT+1]
        if ωtimes[iω] < eptf
          epc = 2
          break
        end
        kω += 1
        iω += 1
      end
      
      # log-probability of observing kω occurrences (Poisson distribution)
      if epc == 2
        Δti  = LTT.t[iLTT] - eptf
        ep  += 1
        epti = eptf
        eptf = ifelse(ep<=lastindex(tep), tep[ep], 0.0)
        epc  = 1
      elseif epc == 1
        Δti = epti - LTT.t[iLTT]
        epc = 0
      else
        Δti = LTT.t[iLTT] - LTT.t[iLTT+1]
      end
      Eω = ni*ω[ep]*Δti  # expected number of occurrences during the LTT step i

      if kω==0
        ll -= Eω
      elseif kω==1
        ll += log(Eω) - Eω
      else
        ll += logpdf(Poisson(Eω), kω)
        
        ### Much slower method when kω>20 ###
        # ll += kω*log(Eω) - Eω - log(factorial(big(kω)))

        ### Very fast Ramanujan's approximation ###
        ### Asymptotic error of 1/(1400*n^3): https://en.wikipedia.org/wiki/Stirling%27s_approximation#Versions_suitable_for_calculators
        # ll += kω*log(Eω) - Eω - (kω*log(kω) - kω + log(kω*(1+4*kω*(1+2*kω))+1/30)/6 + log(pi)/2)

        kω  = 0
      end
    end
  end

  return ll
end




"""
    llrLTT(ξc    ::sTfbd,
           ξp    ::sTfbd,
           bi    ::iBffs,
           ωtimes::Vector{Float64},
           ωvec  ::Vector{Float64},
           LTTc  ::Ltt)

Calculates the difference in log-likelihood between two subtrees `ξc` and `ξp`,
given occurrence sampling rates `ωvec` at `ωtimes` and a current global `LTTc` 
(Lineages Through Time) trajectory.
"""
function llrLTT(ξc    ::sTfbd,
                ξp    ::sTfbd,
                bi    ::iBffs,
                ωtimes::Vector{Float64},
                ωvec  ::Vector{Float64},
                LTTc  ::Ltt)
    
    @show ξc, ξp
    # Get LTT values for ξc and ξp
    tiξ = ti(bi)
    LTTξc = ltt(ξc, tiξ)
    LTTξp = ltt(ξp, tiξ)

    tfξ = min(LTTξc.t[end], LTTξp.t[end])
    tiξ-tfξ > treeheight(ξp) || tiξ-tfξ ≈ treeheight(ξp) || @show LTTξp.t, tiξ-tfξ, treeheight(ξp)

    # ΔLTT  = diff_LTTs(LTTξp, LTTξc)
    ΔLTT  = diff_LTTs_bis(LTTξp, LTTξc)

    # Combine the LTTc and ΔLTT to create the proposal LTT
    # LTTp  = sum_LTTs(LTTc, ΔLTT)   # WARNING: keep constant points?
    LTTp  = sum_LTTs_bis(LTTc, ΔLTT)   # WARNING: keep constant points?
    # @show LTTξc.n
    # @show LTTξc.t
    # @show LTTξp.n
    # @show LTTξp.t
    # @show ΔLTT.n
    # @show ΔLTT.t
    # @show LTTc.n
    # @show LTTc.t
    # @show LTTp.n
    # @show LTTp.t
    @assert all(LTTp.n .> 0)

    # Compute Δt intervals
    idx_min = findlast( x -> ( x >= tiξ || x ≈ tiξ ), LTTp.t)
    idx_max = findfirst(x -> ( x <= tfξ || x ≈ tfξ ), LTTp.t)
    LTTpt = LTTp.t[idx_min:idx_max]
    Δt    = diff(LTTpt)

    # Match the time indices
    ΔLTTidx = [findfirst(x -> (x<=t || x≈t), ΔLTT.t) for t in LTTpt]
    LTTcidx = [findfirst(x -> (x<=t || x≈t), LTTc.t) for t in LTTpt]

    # Initialize the log-likelihood ratio
    llr = 0.0
    kω  = 0
    i   = 1
    liω = lastindex(ωtimes)
    
    j = 1
    for (i, tii) in enumerate(LTTpt[2:end])
      @show i
      
      # count the number of occurrences in the interval
      while j <= liω && ωtimes[j] >= tii
        j  += 1
        kω += 1
      end

      # update the likelihood ratio
      @show j, kω
      if kω == 0
        llr -= ΔLTT.n[ΔLTTidx[i]] * ωvec[j-1] * Δt[i]
      else
        ΔLTTni = ΔLTT.n[ΔLTTidx[i]]
        #@show ΔLTT.n
        #@show ΔLTT.t
        #@show LTTc.n
        #@show LTTc.t
        llr += kω * log(1 + ΔLTTni/LTTc.n[LTTcidx[i]]) - ΔLTTni * ωvec[j-1] * Δt[i]
        kω   = 0
      end
    end

    # Create a mask to remove points with constant n
    li = lastindex(LTTp.n)
    if (li<=2)
      mask = fill(true, li)
    else
      mask = [true ; diff(LTTp.n[1:(end-1)]) .!= 0 ; true]
    end
    return llr, Ltt(LTTp.n[mask], LTTp.t[mask])
end




"""
    sum_LTTs(LTT1::Ltt, LTT2::Ltt)

Calculate the sum of two LTTs with possibly different time points.
"""
function sum_LTTs(LTT1::Ltt, LTT2::Ltt)
  
  # Find the maximum time value between the two LTTs
  ti1, ti2 = LTT1.t[1], LTT2.t[1]
  ti = max(ti1, ti2)
  
  # Prepare the time and lineage vectors for each LTT
  t1, n1 = (ti1 == ti) ? (LTT1.t, LTT1.n) : (vcat(ti, LTT1.t), vcat(0, LTT1.n))
  t2, n2 = (ti2 == ti) ? (LTT2.t, LTT2.n) : (vcat(ti, LTT2.t), vcat(0, LTT2.n))

  lt1, lt2 = lastindex(t1), lastindex(t2)

  # Initialize result vectors with an upper bound on the number of time points
  t = Vector{Float64}(undef, lt1 + lt2)
  n = Vector{Int64}(undef, lt1 + lt2)

  n[1], t[1] = n1[1] + n2[1], ti   # first values
  i1, i2, i_res = 1, 1, 2          # indices
  tii1, tii2 = t1[2], t2[2]        # next time points

  # Iterate through both LTTs and add them
  while i1 < lt1 || i2 < lt2

    # Move to the next time point in one or both LTTs
    if tii1 > tii2
      ti  = tii1
      i1 += 1
      tii1 = i1 < lt1 ? t1[i1+1] : -1.0
    elseif tii1 < tii2
      ti  = tii2
      i2 += 1
      tii2 = i2 < lt2 ? t2[i2+1] : -1.0
    else  # t1[i1+1] == t2[i2+1]
      ti  = tii1
      i1 += 1
      i2 += 1
      tii1 = i1 < lt1 ? t1[i1+1] : -1.0
      tii2 = i2 < lt2 ? t2[i2+1] : -1.0
    end
    
    n_res = n1[i1] + n2[i2]
    
    # Add the time point and lineage count if not constant
    if n_res != n[i_res-1] || (i1==lt1 && i2==lt2)
      n[i_res] = n_res
      t[i_res] = ti
      i_res += 1
    end
  end

  # Resize the result vectors to the actual number of time points
  resize!(t, i_res-1)
  resize!(n, i_res-1)

  # Remove the first element if it's zero
  if n[1] == 0
    popfirst!(t)
    popfirst!(n)
  end

  return Ltt(n, t)
end




"""
    diff_LTTs(LTT1::Ltt, LTT2::Ltt)

Calculate the difference between two LTTs with possibly different time points.
"""
function diff_LTTs(LTT1::Ltt, LTT2::Ltt)
  
  # Find the maximum time value between the two LTTs
  ti1, ti2 = LTT1.t[1], LTT2.t[1]
  ti = max(ti1, ti2)
  
  # Prepare the time and lineage vectors for each LTT
  t1, n1 = (ti1 == ti) ? (LTT1.t, LTT1.n) : (vcat(ti, LTT1.t), vcat(0, LTT1.n))
  t2, n2 = (ti2 == ti) ? (LTT2.t, LTT2.n) : (vcat(ti, LTT2.t), vcat(0, LTT2.n))


  lt1, lt2 = lastindex(t1), lastindex(t2)

  # Initialize result vectors with an upper bound on the number of time points
  t = Vector{Float64}(undef, lt1 + lt2)
  n = Vector{Int64}(undef, lt1 + lt2)

  n[1], t[1] = n1[1] - n2[1], ti   # first values
  i1, i2, i_res = 1, 1, 2          # indices
  tii1, tii2 = t1[2], t2[2]        # next time points

  # Iterate through both LTTs and substract them
  while i1 < lt1 || i2 < lt2

    # Move to the next time point in one or both LTTs
    if tii1 > tii2
      ti  = tii1
      i1 += 1
      tii1 = i1 < lt1 ? t1[i1+1] : -1.0
    elseif tii1 < tii2
      ti  = tii2
      i2 += 1
      tii2 = i2 < lt2 ? t2[i2+1] : -1.0
    else  # t1[i1+1] == t2[i2+1]
      ti  = tii1
      i1 += 1
      i2 += 1
      tii1 = i1 < lt1 ? t1[i1+1] : -1.0
      tii2 = i2 < lt2 ? t2[i2+1] : -1.0
    end
    
    n_res = n1[i1] - n2[i2]
    
    # Add the time point and lineage count if not constant
    if n_res != n[i_res-1] || (i1==lt1 && i2==lt2)
      n[i_res] = n_res
      t[i_res] = ti
      i_res += 1
    end
  end

  # Resize the result vectors to the actual number of time points
  resize!(t, i_res-1)
  resize!(n, i_res-1)

  # Remove the first element if it's zero
  if n[1] == 0
    popfirst!(t)
    popfirst!(n)
  end

  return Ltt(n, t)
end




function sum_LTTs_bis(LTT1::Ltt, LTT2::Ltt)
    # Combine time points from both LTTs and remove duplicates
    all_times = sort(unique(x->abs(round(x, digits=10)), [LTT1.t ; LTT2.t]), rev=true)

    # Interpolate both LTTs on the all_times vector
    n1_interpolated = interpolate_LTT(LTT1, all_times)
    n2_interpolated = interpolate_LTT(LTT2, all_times)

    # Sum the interpolated lineage counts
    n_sum = n1_interpolated .+ n2_interpolated

    return Ltt(n_sum, all_times)

    # # Create a mask to remove points with constant n
    # li = lastindex(n_sum)
    # if (li<=2)
    #   mask = fill(true, li)
    # else
    #   mask = [true ; diff(n_sum[1:(end-1)]) .!= 0 ; true]
    # end

    # return Ltt(n_sum[mask], all_times[mask])
end




function diff_LTTs_bis(LTT1::Ltt, LTT2::Ltt)
    # Combine time points from both LTTs and remove duplicates
    all_times = sort(unique(x->abs(round(x, digits=10)), [LTT1.t ; LTT2.t]), rev=true)

    # Interpolate both LTTs on the all_times vector
    n1_interpolated = interpolate_LTT(LTT1, all_times)
    n2_interpolated = interpolate_LTT(LTT2, all_times)

    # Sum the interpolated lineage counts
    n_sum = n1_interpolated .- n2_interpolated

    return Ltt(n_sum, all_times)

    # # Create a mask to remove points with constant n
    # li = lastindex(n_sum)
    # if (li<=2)
    #   mask = fill(true, li)
    # else
    #   mask = [true ; diff(n_sum[1:(end-1)]) .!= 0 ; true]
    # end

    # return Ltt(n_sum[mask], all_times[mask])
end




"""
    interpolate_LTT(LTT::Ltt, times_union::Vector{Float64})

Interpolate the number of lineages in `LTT` at each time point in `times_union`.
"""
function interpolate_LTT(LTT::Ltt, times_union::Vector{Float64})
  
  # If the LTT is more recent than the time of origin, start with 0 lineage
  tor = times_union[1]
  tv, nv = ifelse(LTT.t[1]≈tor, (LTT.t, LTT.n), ([tor; LTT.t], [0; LTT.n]))
  ltv = lastindex(tv)
  
  n_interp = zeros(length(times_union))
  idx_LTT = 1

  for (i, t) in enumerate(times_union)
      while idx_LTT < ltv && (t <= tv[idx_LTT + 1] || t ≈ tv[idx_LTT + 1])
          idx_LTT += 1
      end

      n_interp[i] = nv[idx_LTT]
  end

  return n_interp
end





# """
#     _llik_cfbd_eventimes!(tree  ::sTfbd,
#                           λ     ::Float64,
#                           μ     ::Float64,
#                           ψ     ::Float64,
#                           t     ::Float64, 
#                           se    ::Vector{Float64}, 
#                           ee    ::Vector{Float64})

# Log-likelihood up to a constant for constant occurrence birth-death 
# given a complete `iTree` recursively.
# """
# function _llik_cfbd_eventimes!(tree  ::sTfbd,
#                                λ     ::Float64,
#                                μ     ::Float64,
#                                ψ     ::Float64,
#                                t     ::Float64, 
#                                se    ::Vector{Float64}, 
#                                ee    ::Vector{Float64})
  
#   defd1 = isdefined(tree, :d1)
#   defd2 = isdefined(tree, :d2)
  
#   et = e(tree)

#   # tip
#   if !defd1 && !defd2
#     # fossil tips are labelled extinct despite no actual extinction event
#     if     isfossil(tree)
#       return - e(tree)*(λ + μ + ψ) + log(ψ)
#     elseif isextinct(tree)
#       push!(ee, t - et)
#       return - e(tree)*(λ + μ + ψ) + log(μ)
#     else                   
#       return - e(tree)*(λ + μ + ψ) 
#     end
  
#   # sampled ancestor
#   elseif defd1 ⊻ defd2
#     return - e(tree)*(λ + μ + ψ) + log(ψ) +
#            (defd1 ? _llik_cfbd_eventimes!(tree.d1::sTfbd, λ, μ, ψ, t - et, se, ee) : 0.0) + 
#            (defd2 ? _llik_cfbd_eventimes!(tree.d2::sTfbd, λ, μ, ψ, t - et, se, ee) : 0.0)
  
#   # bifurcation
#   else
#     push!(se, t - et)
#     return - e(tree)*(λ + μ + ψ) + log(λ) +
#              _llik_cfbd_eventimes!(tree.d1::sTfbd, λ, μ, ψ, t - et, se, ee) + 
#              _llik_cfbd_eventimes!(tree.d2::sTfbd, λ, μ, ψ, t - et, se, ee)
#   end
# end



