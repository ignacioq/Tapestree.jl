#=

constant occurrence birth-death simulation

Jérémy Andréoletti

v(^-^v)

Created 11 02 2022

=#




"""
    cobd_wait(n::Float64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

Sample a waiting time for constant occurrence birth-death when `n` species 
are alive with speciation rate `λ`, extinction rate `μ`, fossil sampling rates 
`ψ` (fossils included in the tree) and `ω`(fossil occurrences).
"""
cobd_wait(n::Float64, λ::Float64, μ::Float64, 
          ψ::Float64, ω::Float64) = rexp(n*(λ + μ + ψ + ω))




"""
    cobd_wait(λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

Sample a per-lineage waiting time for constant occurrence birth-death 
with speciation rate `λ`, extinction rate `μ`, fossil sampling rates 
`ψ` (fossils included in the tree) and `ω`(fossil occurrences).
"""
cobd_wait(λ::Float64, μ::Float64, ψ::Float64, ω::Float64) = rexp(λ + μ + ψ + ω)




"""
    λevent(λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

Return `true` if speciation event.
"""
λevent(λ::Float64, μ::Float64, ψ::Float64, ω::Float64) = (λ/(λ+μ+ψ+ω)) > rand()




"""
    λevent(μ::Float64, ψ::Float64, ω::Float64)

Return `true` if extinction event, conditionned on "not a speciation event".
"""
μevent(μ::Float64, ψ::Float64, ω::Float64) = (μ/(μ + ψ + ω)) > rand()




"""
    ψevent(ψ::Float64, ω::Float64)

Return `true` if the sampled fossil is included in the tree.
"""
ψevent(ψ::Float64, ω::Float64) = (ψ/(ψ + ω)) > rand()




"""
    sim_cobd(t::Float64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64)
    sim_cobd(t::Float64, λ::Float64, μ::Float64, ψ::Vector{Float64}, ω::Vector{Float64})

Simulate a constant occurrence birth-death `iTree` of height `t` with speciation
rate `λ`, extinction rate `μ`, fossil sampling rates `ψ` (fossils included in 
the tree) and `ω`(fossil occurrences).
"""
sim_cobd(t::Float64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64) =
  sim_cobd(t,λ,μ,[ψ],[ω],Float64[],1,1)
sim_cobd(t::Float64, λ::Float64, μ::Float64, ψ::Vector{Float64}, ω::Vector{Float64}, tep::Vector{Float64}) =
  sim_cobd(t,λ,μ,ψ,ω,tep,1,lastindex(tep)+1)




"""
    sim_cobd(t  ::Float64, 
             λ  ::Float64, 
             μ  ::Float64, 
             ψ  ::Vector{Float64},
             ω  ::Vector{Float64},
             tep::Vector{Float64},
             ix ::Int64,
             nep::Int64)

Simulate a constant occurrence birth-death `iTree` of height `t` with speciation 
rate `λ`, extinction rate `μ`, piecewise-constant fossil sampling rates `ψ`
(fossils included in the tree) and `ω`(fossil occurrences), with `nep` epochs 
at times `tep`, starting at index `ix`.
"""
function sim_cobd(t  ::Float64, 
                  λ  ::Float64, 
                  μ  ::Float64, 
                  ψ  ::Vector{Float64},
                  ω  ::Vector{Float64},
                  tep::Vector{Float64},
                  ix ::Int64,
                  nep::Int64)

  @inbounds ψi = ψ[ix]
  @inbounds ωi = ω[ix]

  tw = cobd_wait(λ, μ, ψi, ωi)

  # ψ/ω epoch change
  if ix < nep
    @inbounds tepi = tep[ix]
    if t - tw < tepi
      if tepi < t
        t0, ωtimes = sim_cobd(tepi, λ, μ, ψ, ω, tep, ix + 1, nep)
        adde!(t0, t-tepi)
      else
        t0, ωtimes = sim_cobd(t, λ, μ, ψ, ω, tep, ix + 1, nep)
      end
      return t0, ωtimes
    end
  end

  if tw > t
    return sTfbd(t, false, false, false), Float64[]
  end

  # speciation
  if λevent(λ, μ, ψi, ωi)
    d1, ωtimes  = sim_cobd(t - tw, λ, μ, ψ, ω, tep, ix, nep)
    d2, ωtimes2 = sim_cobd(t - tw, λ, μ, ψ, ω, tep, ix, nep)
    append!(ωtimes, ωtimes2)
    return sTfbd(d1, d2, tw, false, false, false), ωtimes
  
  # extinction
  elseif μevent(μ, ψi, ωi)
    return sTfbd(tw, true, false, false), Float64[]
  
  # fossil sampling (included in the tree)
  elseif ψevent(ψi, ωi)
    t0, ωtimes = sim_cobd(t - tw, λ, μ, ψ, ω, tep, ix, nep)
    return sTfbd(t0, tw, false, true, false), ωtimes
  
  # fossil occurrence sampling (not included in the tree)
  else
    t0, ωtimes = sim_cobd(t - tw, λ, μ, ψ, ω, tep, ix, nep)
    push!(ωtimes, t - tw)
    adde!(t0, tw)
    return t0, ωtimes
  end
end




# """
#     _sim_cobd_t(t     ::Float64,
#                 λ     ::Float64,
#                 μ     ::Float64,
#                 ψ     ::Vector{Float64},
#                 ω     ::Vector{Float64},
#                 tep   ::Vector{Float64},
#                 ωtimes::Vector{Float64},
#                 ix    ::Int64,
#                 nep   ::Int64,
#                 lr    ::Float64,
#                 lU    ::Float64,
#                 Iρi   ::Float64,
#                 na    ::Int64,
#                 nn    ::Int64,
#                 nlim  ::Int64)

# Simulate a constant occurrence birth-death `iTree` of height `t` with speciation
# rate `λ`, extinction rate `μ` and fossilization rate `ψ` for terminal branches,
# conditioned on no fossilizations.
# """
# function _sim_cobd_t(t     ::Float64,
#                      λ     ::Float64,
#                      μ     ::Float64,
#                      ψ     ::Vector{Float64},
#                      ω     ::Vector{Float64},
#                      tep   ::Vector{Float64},
#                      ωtimes::Vector{Float64},
#                      ix    ::Int64,
#                      nep   ::Int64,
#                      lr    ::Float64,
#                      lU    ::Float64,
#                      Iρi   ::Float64,
#                      na    ::Int64,
#                      nn    ::Int64,
#                      nlim  ::Int64)

#   if isfinite(lr) && nn < nlim

#     @inbounds ψi  = ψ[ix]
#     @inbounds ωi  = ω[ix]

#     tw = cobd_wait(λ, μ, ψi, ωi)

#     # ψ epoch change
#     if ix < nep
#       @inbounds tepi = tep[ix]
#       if t - tw < tepi
#         e0 = t - tepi
#         t0, na, nn, lr = 
#           _sim_cobd_t(tepi, λ, μ, ψ, ω, tep, ωtimes, ix + 1, nep, lr, lU, Iρi, na, nn, nlim)
#         sete!(t0, e(t0) + e0)
#         return t0, na, nn, lr
#       end
#     end

#     if tw > t
#       na += 1
#       nlr = lr
#       if na > 1
#         nlr += log(Iρi * Float64(na)/Float64(na-1))
#       end
#       if nlr < lr && lU >= nlr
#         return sTfbd(), na, nn, NaN
#       else
#         return sTfbd(t, false, false, false), na, nn, nlr
#       end
#     end

#     # speciation
#     if λevent(λ, μ, ψi, ωi)
#       nn += 1
#       d1, na, nn, lr =
#         _sim_cobd_t(t - tw, λ, μ, ψ, ω, tep, ωtimes, ix, nep, lr, lU, Iρi, na, nn, nlim)
#       d2, na, nn, lr =
#         _sim_cobd_t(t - tw, λ, μ, ψ, ω, tep, ωtimes, ix, nep, lr, lU, Iρi, na, nn, nlim)

#       return sTfbd(d1, d2, tw, false, false, false), na, nn, lr
    
#     # extinction
#     elseif μevent(μ, ψi, ωi)
#       return sTfbd(tw, true, false, false), na, nn, lr
    
#     # fossil sampling (included in the tree)
#     elseif ψevent(ψi, ωi)
#       return sTfbd(t, false, false, false), na, nn, NaN
  
#     # fossil occurrence sampling (not included in the tree)
#     else
#       t0, na, nn, lr =
#         _sim_cobd_t(t - tw, λ, μ, ψ, ω, tep, ωtimes, ix, nep, lr, lU, Iρi, na, nn, nlim)
#       push!(ωtimes, t - tw)
#       adde!(t0, tw)
#       return t0, na, nn, lr
#     end
#   end

#   return sTfbd(), na, nn, NaN
# end




# """
#     _sim_cobd_i(t     ::Float64,
#                 te    ::Float64,
#                 λ     ::Float64,
#                 μ     ::Float64,
#                 ψ     ::Vector{Float64},
#                 ω     ::Vector{Float64},
#                 tep   ::Vector{Float64},
#                 ωtimes::Vector{Float64},
#                 ix    ::Int64,
#                 nep   ::Int64,
#                 na    ::Int64,
#                 nf    ::Int64,
#                 nn    ::Int64,
#                 nlim  ::Int64)

# Simulate a constant occurrence birth-death `iTree` of height `t` with
# speciation rate `λ`, extinction rate `μ` and fossilization rate `ψ`
# for internal branches, conditioned on no fossilizations.
# """
# function _sim_cobd_i(t     ::Float64,
#                      te    ::Float64,
#                      λ     ::Float64,
#                      μ     ::Float64,
#                      ψ     ::Vector{Float64},
#                      ω     ::Vector{Float64},
#                      tep   ::Vector{Float64},
#                      ωtimes::Vector{Float64},
#                      ix    ::Int64,
#                      nep   ::Int64,
#                      na    ::Int64,
#                      nf    ::Int64,
#                      nn    ::Int64,
#                      nlim  ::Int64)

#   if iszero(nf) && nn < nlim

#     @inbounds ψi  = ψ[ix]
#     @inbounds ωi  = ω[ix]

#     tw = cobd_wait(λ, μ, ψi, ωi)

#     # ψ epoch change
#     if ix < nep
#       @inbounds tepi = tep[ix]
#       if t - tw < tepi > te
#         e0 = t - tepi
#         t0, na, nf, nn  = 
#           _sim_cobd_i(tepi, te, λ, μ, ψ, ω, tep, ωtimes, ix + 1, nep, na, nf, nn, nlim)
#         sete!(t0, e(t0) + e0)
#         return t0, na, nf, nn
#       end
#     end

#     if tw > (t - te)
#       na += 1
#       return sTfbd(t - te, false, false, false), na, nf, nn
#     end

#     # speciation
#     if λevent(λ, μ, ψi, ωi)
#       nn += 1
#       d1, na, nf, nn = 
#         _sim_cobd_i(t - tw, te, λ, μ, ψ, ω, tep, ωtimes, ix, nep, na, nf, nn, nlim)
#       d2, na, nf, nn = 
#         _sim_cobd_i(t - tw, te, λ, μ, ψ, ω, tep, ωtimes, ix, nep, na, nf, nn, nlim)

#       return sTfbd(d1, d2, tw, false, false, false), na, nf, nn
    
#     # extinction
#     elseif μevent(μ, ψi, ωi)
#       return sTfbd(tw, true, false, false), na, nf, nn
    
#     # fossil sampling (included in the tree)
#     elseif ψevent(ψi, ωi)
#       return sTfbd(), na, 1, nn
  
#     # fossil occurrence sampling (not included in the tree)
#     else
#       t0, na, nf, nn =
#         _sim_cobd_i(t - tw, te, λ, μ, ψ, ω, tep, ωtimes, ix, nep, na, nf, nn, nlim)
#       push!(ωtimes, t - tw)
#       adde!(t0, tw)
#       return t0, na, nf, nn
#     end
#   end

#   return sTfbd(), na, nf, nn
# end




# """
#     _sim_cobd_it(t     ::Float64,
#                  λ     ::Float64,
#                  μ     ::Float64,
#                  ψ     ::Vector{Float64},
#                  ω     ::Vector{Float64},
#                  tep   ::Vector{Float64},
#                  ωtimes::Vector{Float64},
#                  ix    ::Int64,
#                  nep   ::Int64,
#                  lr    ::Float64,
#                  Iρi   ::Float64,
#                  na    ::Int64,
#                  nn    ::Int64,
#                  nlim  ::Int64)

# Simulate a constant occurrence birth-death `iTree` of height `t` with
# speciation rate `λ`, extinction rate `μ` and fossilization rate `ψ`
# for continuing internal branches, conditioned on no fossilizations.
# """
# function _sim_cobd_it(t     ::Float64,
#                       λ     ::Float64,
#                       μ     ::Float64,
#                       ψ     ::Vector{Float64},
#                       ω     ::Vector{Float64},
#                       tep   ::Vector{Float64},
#                       ωtimes::Vector{Float64},
#                       ix    ::Int64,
#                       nep   ::Int64,
#                       lr    ::Float64,
#                       Iρi   ::Float64,
#                       na    ::Int64,
#                       nn    ::Int64,
#                       nlim  ::Int64)

#   if isfinite(lr) && nn < nlim

#     @inbounds ψi  = ψ[ix]
#     @inbounds ωi  = ω[ix]

#     tw = cobd_wait(λ, μ, ψi, ωi)

#     # ψ epoch change
#     if ix < nep
#       @inbounds tepi = tep[ix]
#       if t - tw < tepi
#         e0 = t - tepi
#         t0, na, nn, lr = 
#           _sim_cobd_it(tepi, λ, μ, ψ, ω, tep, ωtimes, ix + 1, nep, lr, Iρi, na, nn, nlim)
#         sete!(t0, e(t0) + e0)
#         return t0, na, nn, lr
#       end
#     end

#     if tw > t
#       na += 1
#       lr += log(Iρi)
#       return sTfbd(t, false, false, false), na, nn, lr
#     end

#     # speciation
#     if λevent(λ, μ, ψi, ωi)
#       nn += 1
#       d1, na, nn, lr =
#         _sim_cobd_it(t - tw, λ, μ, ψ, ω, tep, ωtimes, ix, nep, lr, Iρi, na, nn, nlim)
#       d2, na, nn, lr =
#         _sim_cobd_it(t - tw, λ, μ, ψ, ω, tep, ωtimes, ix, nep, lr, Iρi, na, nn, nlim)

#       return sTfbd(d1, d2, tw, false, false, false), na, nn, lr
   
#     # extinction
#     elseif μevent(μ, ψi, ωi)
#       return sTfbd(tw, true, false, false), na, nn, lr
    
#     # fossil sampling (included in the tree)
#     elseif ψevent(ψi, ωi)
#       return sTfbd(), na, nn, NaN
  
#     # fossil occurrence sampling (not included in the tree)
#     else
#       t0, na, nn, lr =
#         _sim_cobd_it(t - tw, λ, μ, ψ, ω, tep, ωtimes, ix, nep, lr, Iρi, na, nn, nlim)
#       push!(ωtimes, t - tw)
#       adde!(t0, tw)
#       return t0, na, nn, lr
#     end
#   end

#   return sTfbd(), na, nn, NaN
# end




# """
#     sim_cobd(t::Float64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64, 
#              na::Int64, nfos::Int64)

# Simulate a constant occurrence birth-death `iTree` of height `t` with speciation 
# rate `λ`, extinction rate `μ`, fossil sampling rates `ψ` (fossils included in the 
# tree) and `ω`(fossil occurrences). `na` initializes the number of alived tips.
# """
# function sim_cobd(t::Float64, 
#                   λ::Float64, 
#                   μ::Float64, 
#                   ψ::Float64,
#                   ω::Float64,
#                   na::Int64,
#                   nfos::Int64)

#   tw = cobd_wait(λ, μ, ψ, ω)
#   ωtimes = Float64[]

#   if tw > t
#     na += 1
#     return sTfbd(t), ωtimes, na, nfos
#   end

#   if λevent(λ, μ, ψ, ω)
#     # speciation
#     d1, ωtimes1, na, nfos = sim_cobd(t - tw, λ, μ, ψ, ω, na, nfos)
#     d2, ωtimes2, na, nfos = sim_cobd(t - tw, λ, μ, ψ, ω, na, nfos)
#     append!(ωtimes1, ωtimes2)
#     return sTfbd(d1, d2, tw), ωtimes1, na, nfos

#   elseif μevent(μ, ψ, ω)
#     # extinction
#     return sTfbd(tw, true), ωtimes, na, nfos

#   elseif ψevent(ψ, ω)
#     # fossil sampling (included in the tree)
#     nfos += 1
#     @show nfos
#     d, ωtimes, na, nfos = sim_cobd(t - tw, λ, μ, ψ, ω, na, nfos)
#     return sTfbd(d, tw, false, true, false), ωtimes, na, nfos
  
#   else
#     # fossil occurrence sampling (not included in the tree)
#     d, ωtimes, na, nfos = sim_cobd(t - tw, λ, μ, ψ, ω, na, nfos)
#     push!(ωtimes, t - tw)
#     adde!(d, tw)
#     return d, ωtimes, na, nfos
#   end
# end




# """
#    sim_cobd_b(n::Int64, λ::Float64, μ::Float64, ψ::Float64, ω::Float64)

# Simulate constant occurrence birth-death in backward time.
# """
# function sim_cobd_b(n::Int64, 
#                     λ::Float64, 
#                     μ::Float64, 
#                     ψ::Float64,
#                     ω::Float64)

#   nF = Float64(n)
#   nI = n
#   ωtimes = Float64[]

#   # disjoint trees vector 
#   tv = sTfbd[]
#   for i in Base.OneTo(nI)
#     push!(tv, sTfbd(0.0))
#   end

#   th = 0.0

#   # start simulation
#   while true
#     w = cobd_wait(nF, λ, μ, ψ, ω)

#     # track backward time
#     th += w

#     for t in tv
#       adde!(t, w)
#     end

#     # if speciation
#     if λevent(λ, μ, ψ, ω)
#       if isone(nI)
#         return tv[nI], ωtimes
#       else
#         j, k = samp2(Base.OneTo(nI))
#         tv[j] = sTfbd(tv[j], tv[k], 0.0)
#         deleteat!(tv,k)
#         nI -= 1
#         nF -= 1.0
#       end
    
#     # if extinction
#     elseif μevent(μ, ψ, ω)
#       nI += 1
#       nF += 1.0
#       push!(tv, sTfbd(0.0, true))
    
#     # if fossil sampling (included in the tree)
#     elseif ψevent(ψ, ω)
#       j = rand(Base.OneTo(nI))
#       tv[j] = sTfbd(tv[j], 0.0, false, true, false)
    
#     # if fossil occurrence sampling (not included in the tree)
#     else
#       push!(ωtimes, th)
#     end
#   end
# end




# """
#     sim_cobd_b(λ::Float64, 
#                μ::Float64, 
#                ψ::Float64,
#                ω::Float64, 
#                mxth::Float64, 
#                maxn::Int64)

# Simulate constant occurrence birth-death in backward time conditioned on 
# 1 survival and not having a greater tree height than `mxth`.
# """
# function sim_cobd_b(λ::Float64, 
#                     μ::Float64, 
#                     ψ::Float64,
#                     ω::Float64, 
#                     mxth::Float64, 
#                     maxn::Int64)

#   nF = 1.0
#   nI = 1
#   ωtimes = Float64[]

#   # disjoint trees vector 
#   tv = [sTfbd(0.0, false)]

#   th = 0.0

#   # start simulation
#   while true
#     w   = cobd_wait(nF, λ, μ, ψ, ω)

#     # track backward time
#     th += w

#     if nI > maxn
#       return tv[nI], ωtimes, (mxth + 0.1)
#     end
    
#     if th > mxth 
#      return tv[nI], ωtimes, th
#     end

#     for t in tv
#       adde!(t, w)
#     end

#     # if speciation
#     if λevent(λ, μ, ψ, ω)
#       if isone(nI)
#         return tv[nI], ωtimes, th
#       else
#         j, k = samp2(Base.OneTo(nI))
#         tv[j] = sTfbd(tv[j], tv[k], 0.0)
#         deleteat!(tv,k)
#         nI -= 1
#         nF -= 1.0
#       end
    
#     # if extinction
#     elseif μevent(μ, ψ, ω)
#       nI += 1
#       nF += 1.0
#       push!(tv, sTfbd(0.0, true))
    
#     # if fossil sampling (included in the tree)
#     elseif ψevent(ψ, ω)
#       j = rand(Base.OneTo(nI))
#       tv[j] = sTfbd(tv[j], 0.0, false, true, false)
    
#     # if fossil occurrence sampling (not included in the tree)
#     else
#       push!(ωtimes, th)
#     end
#   end
# end





