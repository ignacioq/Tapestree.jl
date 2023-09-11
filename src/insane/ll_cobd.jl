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
    llik_cobd(tree  ::sTfbd,
              ωtimes::Vector{Float64}, 
              LTT   ::Ltt,
              λ     ::Float64, 
              μ     ::Float64, 
              ψ     ::Float64,
              ω     ::Float64)
    llik_cobd(tree  ::sTfbd,
              ωtimes::Vector{Float64}, 
              LTT   ::Ltt,
              λ     ::Float64, 
              μ     ::Float64, 
              ψ     ::Vector{Float64}, 
              ω     ::Vector{Float64},
              tep   ::Vector{Float64})

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
llik_cobd(tree  ::sTfbd,
          ωtimes::Vector{Float64}, 
          LTT   ::Ltt,
          λ     ::Float64, 
          μ     ::Float64, 
          ψ     ::Float64,
          ω     ::Float64) =
  llik_cobd(tree, ωtimes, LTT, λ, μ, [ψ], [ω], treeheight(tree), Float64[], 1, 1)
llik_cobd(tree  ::sTfbd,
          ωtimes::Vector{Float64}, 
          LTT   ::Ltt,
          λ     ::Float64, 
          μ     ::Float64, 
          ψ     ::Vector{Float64}, 
          ω     ::Vector{Float64},
          tep   ::Vector{Float64}) =
  llik_cobd(tree, ωtimes, LTT, λ, μ, ψ, ω, treeheight(tree), tep, 1, lastindex(tep)+1)




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

    isapprox(LTTt[end], 0.0, atol=1e-15) && (LTTt[end] = 0.0)

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
    
    tep0 = push!(copy(tep), 0.0)
    epti = LTT.t[1]        # beggining of current epoch
    ep   = findfirst(x -> x<=epti, tep0)  # current epoch
    eptf = tep0[ep]        # end of current epoch
    ωep  = ω[ep]           # ω value at current epoch
    ep_change = false      # whether an iteration is an epoch change
    ep_before = false      # whether the interval is before an epoch change

    for (iLTT, ni) in enumerate(LTT.n[1:(end-1)])
      # @show iLTT, ni
      end_reached = false

      # check whether the current epoch ends in this interval
      if LTT.t[iLTT+1] < eptf
        ep_change = true
        ep_before = true
      end

      while !end_reached
        # count the number of occurrences between each LTT step or epoch
        while iω <= lω && ωtimes[iω] >= max(LTT.t[iLTT+1], eptf)
          kω += 1
          iω += 1
        end
        
        # compute the time interval between LTT steps or epoch changes
        # no epoch change during the time interval
        if !ep_change
          Δti = LTT.t[iLTT] - LTT.t[iLTT+1]
          end_reached = true
        # interval before an epoch change
        elseif ep_before
          Δti  = LTT.t[iLTT] - eptf
          ep  += 1
          epti = eptf
          eptf = tep0[ep]
          ep_before = false
        # interval after an epoch change
        else
          ωep = ω[ep]
          if LTT.t[iLTT+1] < eptf
            Δti = epti - eptf
            ep_change = true
            ep_before = true
          else
            Δti = epti - LTT.t[iLTT+1]
            ep_change   = false
            end_reached = true
          end
        end
      
        ll -= ni*ωep*Δti  # expected number of occurrences during the LTT step i
        # @show ni, ωep, Δti, kω
        
        if kω > 0
          ll += kω*log(ni*ωep)
          kω  = 0
        end
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
           ω     ::Vector{Float64},
           tep   ::Vector{Float64},
           LTTc  ::Ltt,
           ixi   ::Int64)

Calculates the difference in log-likelihood between two subtrees `ξc` and `ξp`,
given a record of occurrences `ωtimes` and a current global `LTTc` trajectory.
"""
function llrLTT(ξc    ::sTfbd,
                ξp    ::sTfbd,
                bi    ::iBffs,
                ωtimes::Vector{Float64},
                ω     ::Vector{Float64},
                tep   ::Vector{Float64},
                LTTc  ::Ltt,
                ixi   ::Int64)
    
  # @show ξc, ξp
  nep = lastindex(tep) + 1
  # Get LTT values for ξc and ξp
  tiξ   = ti(bi)
  LTTξc = ltt(ξc, tiξ)
  LTTξp = ltt(ξp, tiξ)

  tfξ = min(LTTξc.t[end], LTTξp.t[end])
  @assert tiξ-tfξ > treeheight(ξp) || tiξ-tfξ ≈ treeheight(ξp) "LTTξp.t, tiξ-tfξ, treeheight(ξp) =", LTTξp.t, tiξ-tfξ, treeheight(ξp)

  ΔLTT  = diff_LTTs(LTTξp, LTTξc)

  # @show ΔLTT.n, ΔLTT.t
  if any(ΔLTT.n .!= 0)
    idx_min = searchsortedlast(LTTc.t, tiξ-1e-12, rev=true)
    idx_max = searchsortedfirst(LTTc.t, tfξ+1e-12, rev=true)
    LTTct = LTTc.t[idx_min:idx_max]
    LTTcn = LTTc.n[idx_min:idx_max]
    
    ΔLTT.t[1] = LTTct[1]
    steps = copy(ΔLTT.t)
    # @show LTTcn, LTTct
    append!(steps, LTTct[2:end])
    sort!(steps, rev=true)
    # @show steps

    # Match the time indices
    ΔLTTidx = [searchsortedlast(ΔLTT.t, t-1e-12, rev=true) for t in steps]
    LTTcidx = [searchsortedlast(LTTct, t-1e-12, rev=true) for t in steps]
    # @show ΔLTTidx
    # @show LTTcidx

    i = searchsortedfirst(ωtimes, steps[1], rev=true)
    if i <= lastindex(ωtimes)
      nsteps  = lastindex(steps)
      kω_ΔLTT = zeros(Int, nsteps)
      step    = 1
      while i <= lastindex(ωtimes) && ωtimes[i]>steps[end]
        # @show ωtimes[i], ΔLTT.t[step]
        while step < nsteps && ωtimes[i] < steps[step+1]
          step += 1
          # @show step, ΔLTT.t[step]
        end
        kω_ΔLTT[step] += 1
        i += 1
      end

      llr_kω = sum(kω_ΔLTT .* log.(1 .+ ΔLTT.n[ΔLTTidx] ./ LTTcn[LTTcidx]))
    else
      llr_kω = 0
    end
    # @show kω_ΔLTT
    # @show ΔLTT.n[ΔLTTidx], ΔLTT.t[ΔLTTidx]
    # @show LTTcn[LTTcidx], LTTct[LTTcidx]

    Lc = zeros(nep)
    Lp = zeros(nep)
    _treelength!(ξc, tiξ, Lc, tep, ixi, nep)
    _treelength!(ξp, tiξ, Lp, tep, ixi, nep)

    llr_L = sum((Lc .- Lp) .* ω)
    
    LTTc  = sum_LTTs(LTTc, ΔLTT)
  else
    llr_kω = llr_L = 0.0
  end

  return llr_kω + llr_L, LTTc
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
  t1, n1 = isapprox(ti1, ti; atol=1e-12) ? (LTT1.t, LTT1.n) : (vcat(ti, LTT1.t), vcat(0, LTT1.n))
  t2, n2 = isapprox(ti2, ti; atol=1e-12) ? (LTT2.t, LTT2.n) : (vcat(ti, LTT2.t), vcat(0, LTT2.n))

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
    if isapprox(tii1, tii2; atol=1e-12)
      ti  = tii1
      i1 += 1
      i2 += 1
      tii1 = i1 < lt1 ? t1[i1+1] : -1.0
      tii2 = i2 < lt2 ? t2[i2+1] : -1.0
    elseif tii1 > tii2
      ti  = tii1
      i1 += 1
      tii1 = i1 < lt1 ? t1[i1+1] : -1.0
    else # tii1 < tii2
      ti  = tii2
      i2 += 1
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
  t1, n1 = isapprox(ti1, ti; atol=1e-12) ? (LTT1.t, LTT1.n) : (vcat(ti, LTT1.t), vcat(0, LTT1.n))
  t2, n2 = isapprox(ti2, ti; atol=1e-12) ? (LTT2.t, LTT2.n) : (vcat(ti, LTT2.t), vcat(0, LTT2.n))

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
    if isapprox(tii1, tii2; atol=1e-12)
      ti  = tii1
      i1 += 1
      i2 += 1
      tii1 = i1 < lt1 ? t1[i1+1] : -1.0
      tii2 = i2 < lt2 ? t2[i2+1] : -1.0
    elseif tii1 > tii2
      ti  = tii1
      i1 += 1
      tii1 = i1 < lt1 ? t1[i1+1] : -1.0
    else # tii1 < tii2
      ti  = tii2
      i2 += 1
      tii2 = i2 < lt2 ? t2[i2+1] : -1.0
    end
    
    n_res = n1[i1] - n2[i2]
    
    # Add the time point and lineage count if not constant
    # if n_res != n[i_res-1] || (i1==lt1 && i2==lt2)
    if n_res != n[i_res-1]
      n[i_res] = n_res
      t[i_res] = ti
      i_res += 1
    end
  end

  # Resize the result vectors to the actual number of time points
  resize!(t, i_res-1)
  resize!(n, i_res-1)

  return Ltt(n, t)
end




# """
#     ω_llr(ωtimes::Vector{Float64},
#           ωc    ::Vector{Float64},
#           ωp    ::Vector{Float64},
#           tep   ::Vector{Float64},
#           LTT   ::Ltt)

# Log-likelihood ratio of the observed occurrences for a given `LTT` 
# (Lineages Through Time) trajectory.
# """
# function ω_llr(ωtimes::Vector{Float64},
#                ωc    ::Vector{Float64},
#                ωp    ::Vector{Float64},
#                tep   ::Vector{Float64},
#                LTT   ::Ltt)
  
#   @inbounds begin
  
#     llr = fill(0.0, lastindex(ωc))
    
#     lω = lastindex(ωtimes) # total number of occurrences
#     kω = 0                 # number of occurrences per interval
#     iω = 1                 # occurrence index
    
#     tep0 = push!(copy(tep), 0.0)
#     epti = LTT.t[1]        # beggining of current epoch
#     ep   = findfirst(x -> x<=epti, tep0)  # current epoch
#     eptf = tep0[ep]        # end of current epoch
#     ωcep  = ωc[ep]         # ωc value at current epoch
#     ωpep  = ωp[ep]         # ωp value at current epoch
#     ep_change = false      # whether an iteration is an epoch change
#     ep_before = false      # whether the interval is before an epoch change

#     for (iLTT, ni) in enumerate(LTT.n[1:(end-1)])
#       # @show iLTT, ni
#       end_reached = false

#       # check whether the current epoch ends in this interval
#       if LTT.t[iLTT+1] < eptf
#         ep_change = true
#         ep_before = true
#       end

#       while !end_reached
#         # count the number of occurrences between each LTT step or epoch
#         while iω <= lω && ωtimes[iω] >= max(LTT.t[iLTT+1], eptf)
#           kω += 1
#           iω += 1
#         end
        
#         # compute the time interval between LTT steps or epoch changes
#         # no epoch change during the time interval
#         if !ep_change
#           Δti = LTT.t[iLTT] - LTT.t[iLTT+1]
#           end_reached = true
#         # interval before an epoch change
#         elseif ep_before
#           Δti  = LTT.t[iLTT] - eptf
#           epti = eptf
#           eptf = tep0[ep + 1]
#         # interval after an epoch change
#         else
#           ωcep = ωc[ep]
#           ωpep = ωp[ep]
#           if LTT.t[iLTT+1] < eptf
#             Δti = epti - eptf
#             ep_change = true
#             ep_before = true
#           else
#             Δti = epti - LTT.t[iLTT+1]
#             ep_change   = false
#             end_reached = true
#           end
#         end

#         llr[ep] -= ni*(ωpep-ωcep)*Δti
        
#         # log-probability ratio of observing kω occurrences
#         if kω > 0
#           llr[ep] += kω*log(ωpep/ωcep)
#           kω  = 0
#         end

#         if ep_before
#           ep += 1
#           ep_before = false
#         end
#       end
#     end
#   end

#   return llr
# end




# """
#     llrLTT(ξc    ::sTfbd,
#            ξp    ::sTfbd,
#            bi    ::iBffs,
#            ωtimes::Vector{Float64},
#            ω     ::Vector{Float64},
#            tep   ::Vector{Float64},
#            LTTc  ::Ltt,
#            ixi   ::Int64)

# Calculates the difference in log-likelihood between two subtrees `ξc` and `ξp`,
# given a record of occurrences `ωtimes` and a current global `LTTc` trajectory.
# """
# function llrLTT(ξc    ::sTfbd,
#                 ξp    ::sTfbd,
#                 bi    ::iBffs,
#                 ωtimes::Vector{Float64},
#                 ω     ::Vector{Float64},
#                 tep   ::Vector{Float64},
#                 LTTc  ::Ltt,
#                 ixi   ::Int64)
    
#   # @show ξc, ξp
#   # Get LTT values for ξc and ξp
#   tiξ   = ti(bi)
#   LTTξc = ltt(ξc, tiξ)
#   LTTξp = ltt(ξp, tiξ)

#   tfξ = min(LTTξc.t[end], LTTξp.t[end])
#   @assert tiξ-tfξ > treeheight(ξp) || tiξ-tfξ ≈ treeheight(ξp) "LTTξp.t, tiξ-tfξ, treeheight(ξp) =", LTTξp.t, tiξ-tfξ, treeheight(ξp)

#   # ΔLTT  = diff_LTTs(LTTξp, LTTξc)
#   ΔLTT  = diff_LTTs_bis(LTTξp, LTTξc)

#   # Combine the LTTc and ΔLTT to create the proposal LTT
#   # LTTp  = sum_LTTs(LTTc, ΔLTT)   # WARNING: keep constant points?
#   LTTp  = sum_LTTs_bis(LTTc, ΔLTT)   # WARNING: remove points?
#   # @show LTTξc.n
#   # @show LTTξc.t
#   # @show LTTξp.n
#   # @show LTTξp.t
#   # @show ΔLTT.n
#   # @show ΔLTT.t
#   # @show LTTc.n
#   # @show LTTc.t
#   # @show LTTp.n
#   # @show LTTp.t
#   # @show sum_LTTs_bis(LTTc, ΔLTT).n
#   # @show sum_LTTs_bis(LTTc, ΔLTT).t
#   @assert all(LTTp.n .> 0)

#   # Compute Δt intervals
#   idx_min = findlast( x -> ( x >= tiξ-1e-12 ), LTTp.t)
#   idx_max = findfirst(x -> ( x <= tfξ+1e-12 ), LTTp.t)
#   LTTpt = LTTp.t[idx_min:idx_max]
#   # @show LTTpt

#   # Match the time indices
#   ΔLTTidx = [findlast(x -> (x >= t-1e-12), ΔLTT.t) for t in LTTpt]
#   LTTcidx = [findlast(x -> (x >= t-1e-12), LTTc.t) for t in LTTpt]
#   # @show ΔLTTidx
#   # @show LTTcidx

#   # Initialize indices
#   lω = lastindex(ωtimes) # total number of occurrences
#   kω = 0                 # number of occurrences per interval
#   iω = 1                 # occurrence index

#   tep0 = push!(copy(tep), 0.0)
#   epti = LTTpt[1]        # beggining of current epoch
#   ep   = findfirst(x -> x<=epti, tep0)  # current epoch
#   eptf = tep0[ep]        # end of current epoch
#   ωep  = ω[ep]           # ω value at current epoch
#   ep_change = false      # whether an iteration is an epoch change
#   ep_before = false      # whether the interval is before an epoch change

#   llr = 0.0              # log-likelihood ratio
  
#   for (iLTT, tii) in enumerate(LTTpt[1:(end-1)])
#     end_reached = false

#     # check whether the current epoch ends in this interval
#     if LTTpt[iLTT+1] < eptf
#       ep_change = true
#       ep_before = true
#     end

#     while !end_reached
#       # count the number of occurrences between each LTT step or epoch
#       while iω <= lω && ωtimes[iω] >= max(LTTpt[iLTT+1], eptf)
#         kω += 1
#         iω += 1
#       end
#       # @show ep_change, ep_before
      
#       # compute the time interval between LTT steps or epoch changes
#       # no epoch change during the time interval
#       if !ep_change
#         Δti = tii - LTTpt[iLTT+1]
#         end_reached = true
#       # interval before an epoch change
#       elseif ep_before
#         Δti  = tii - eptf
#         ep  += 1
#         epti = eptf
#         eptf = tep0[ep]
#         ep_before = false
#       # interval after an epoch change
#       else
#         ωep = ω[ep]
#         if LTTpt[iLTT+1] < eptf
#           Δti = epti - eptf
#           ep_change = true
#           ep_before = true
#         else
#           Δti = epti - LTTpt[iLTT+1]
#           ep_change   = false
#           end_reached = true
#         end
#       end

#       # update the likelihood ratio for kω observed occurrences
#       ΔLTTni = ΔLTT.n[ΔLTTidx[iLTT]]
#       # @show ΔLTTni, ωep, Δti, kω
#       if ΔLTTni != 0
#         llr -= ΔLTTni * ωep * Δti
#         if kω > 0
#           # @show ΔLTT.n
#           # @show ΔLTT.t
#           # @show LTTc.n
#           # @show LTTc.t
#           # @show LTTc.n[LTTcidx[iLTT]]
#           llr += kω * log(1 + ΔLTTni/LTTc.n[LTTcidx[iLTT]])
#         end
#       end

#       kω = 0
#     end
#   end

#   # Create a mask to remove points with constant n
#   li = lastindex(LTTp.n)
#   if (li<=2)
#     mask = fill(true, li)
#   else
#     mask = [true ; diff(LTTp.n[1:(end-1)]) .!= 0 ; true]
#   end

#   return llr, Ltt(LTTp.n[mask], LTTp.t[mask])
# end




# function sum_LTTs_bis(LTT1::Ltt, LTT2::Ltt)
#     # Combine time points from both LTTs and remove duplicates
#     all_times = sort(unique(x->abs(round(x, digits=10)), [LTT1.t ; LTT2.t]), rev=true)

#     # Interpolate both LTTs on the all_times vector
#     n1_interpolated = interpolate_LTT(LTT1, all_times)
#     n2_interpolated = interpolate_LTT(LTT2, all_times)

#     # Sum the interpolated lineage counts
#     n_sum = n1_interpolated .+ n2_interpolated

#     return Ltt(n_sum, all_times)

#     # # Create a mask to remove points with constant n
#     # li = lastindex(n_sum)
#     # if (li<=2)
#     #   mask = fill(true, li)
#     # else
#     #   mask = [true ; diff(n_sum[1:(end-1)]) .!= 0 ; true]
#     # end

#     # return Ltt(n_sum[mask], all_times[mask])
# end




# function diff_LTTs_bis(LTT1::Ltt, LTT2::Ltt)
#     # Combine time points from both LTTs and remove duplicates
#     all_times = sort(unique(x->abs(round(x, digits=10)), [LTT1.t ; LTT2.t]), rev=true)

#     # Interpolate both LTTs on the all_times vector
#     n1_interpolated = interpolate_LTT(LTT1, all_times)
#     n2_interpolated = interpolate_LTT(LTT2, all_times)

#     # Sum the interpolated lineage counts
#     n_sum = n1_interpolated .- n2_interpolated

#     return Ltt(n_sum, all_times)

#     # # Create a mask to remove points with constant n
#     # li = lastindex(n_sum)
#     # if (li<=2)
#     #   mask = fill(true, li)
#     # else
#     #   mask = [true ; diff(n_sum[1:(end-1)]) .!= 0 ; true]
#     # end

#     # return Ltt(n_sum[mask], all_times[mask])
# end




# """
#     interpolate_LTT(LTT::Ltt, times_union::Vector{Float64})

# Interpolate the number of lineages in `LTT` at each time point in `times_union`.
# """
# function interpolate_LTT(LTT::Ltt, times_union::Vector{Float64})
  
#   # If the LTT is more recent than the time of origin, start with 0 lineage
#   tor = times_union[1]
#   tv, nv = ifelse(LTT.t[1]≈tor, (LTT.t, LTT.n), ([tor; LTT.t], [0; LTT.n]))
#   ltv = lastindex(tv)
  
#   n_interp = zeros(length(times_union))
#   idx_LTT = 1

#   for (i, t) in enumerate(times_union)
#       while idx_LTT < ltv && (t <= tv[idx_LTT + 1] || isapprox(t, tv[idx_LTT + 1]; atol=1e-12))
#           idx_LTT += 1
#       end

#       n_interp[i] = nv[idx_LTT]
#   end

#   return n_interp
# end



