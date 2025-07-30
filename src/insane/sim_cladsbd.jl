#=

clads simulation

Ignacio Quintero Mächler

t(-_-t)

Created 28 07 2025
=#




#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Sample conditional on number of species
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#



# """
#     sim_cladsbd(n       ::Int64;
#               λ0      ::Float64 = 1.0,
#               α       ::Float64 = 0.0,
#               σλ      ::Float64 = 0.1,
#               δt      ::Float64 = 1e-3,
#               init    ::Symbol  = :stem,
#               nstar   ::Int64   = n + 2,
#               p       ::Float64 = 5.0,
#               warnings::Bool    = true,
#               maxt    ::Float64 = δt*1e7)

# Simulate `cTbd` according to a birth-death clads.
# """
# function sim_cladsbd(n       ::Int64;
#                    λ0      ::Float64 = 1.0,
#                    α       ::Float64 = 0.0,
#                    σλ      ::Float64 = 0.1,
#                    δt      ::Float64 = 1e-3,
#                    init    ::Symbol  = :stem,
#                    nstar   ::Int64   = n + 2,
#                    p       ::Float64 = 5.0,
#                    warnings::Bool    = true,
#                    maxt    ::Float64 = δt*1e7)

#   # simulate in non-recursive manner
#   e0, e1, el, λs, ea, na, simt =
#     _sedges_cladsbd(nstar, log(λ0), α, σλ, δt, sqrt(δt), init, maxt)

#   if simt >= maxt
#     warnings && @warn "simulation surpassed maximum time"
#     return cTbd()
#   end

#   # transform to iTree
#   t = cTbd(e0, e1, el, λs, ea, e1[1], 1, δt)

#   # sample a time when species(t) == `n`
#   nt = ltt(t)
#   tn = times_n(n, nt)
#   c  = usample(tn, p)

#   if iszero(c)
#     warnings && @warn "tree not sampled, try increasing `p`"
#     return cTbd()
#   else
#     # cut the tree
#     t = cutbottom(t, simt - c)
#     return t
#   end
# end





# """
#     _sedges_cladsbd(n    ::Int64,
#                   λ0   ::Float64,
#                   α    ::Float64,
#                   σλ   ::Float64,
#                   δt   ::Float64,
#                   srδt ::Float64,
#                   init::Symbol)

# Simulate `cladspb` just until hitting `n` alive species. Note that this is
# a biased sample for a tree conditional on `n` species.
# """
# function _sedges_cladsbd(n    ::Int64,
#                        λ0   ::Float64,
#                        α    ::Float64,
#                        σλ   ::Float64,
#                        δt   ::Float64,
#                        srδt ::Float64,
#                        init::Symbol,
#                        maxt ::Float64)

#   # edges
#   e0 = Int64[]
#   e1 = Int64[]

#   if init == :stem
#     # edges alive
#     ea = [1]
#     # first edge
#     push!(e0, 1)
#     push!(e1, 2)
#     # max index
#     mxi0 = n*2
#     # edge lengths
#     el = [0.0]
#     # lambda vector for each edge
#     λs = [Float64[]]
#     # starting speciation rate
#     push!(λs[1], λ0)
#     # lastindex for each edge
#     li = [1]

#     na = 1 # current number of alive species
#     ne = 2 # current maximum node number

#   elseif init == :crown
#     # edges alive
#     ea = [2, 3]
#     # first edges
#     push!(e0, 1, 2, 2)
#     push!(e1, 2, 3, 4)
#     # max index
#     mxi0 = n*2
#     # edge lengths
#     el = [0.0, 0.0, 0.0]
#     # lambda vector for each edge
#     λs = [Float64[], Float64[], Float64[]]
#     # starting speciation rate
#     push!(λs[1], λ0, λ0)
#     push!(λs[2], λ0)
#     push!(λs[3], λ0)
#     # lastindex for each edge
#     li = [2, 1, 1]

#     na = 2 # current number of alive species
#     ne = 4 # current maximum node number

#   else
#     @error "$init does not match stem or crown"
#   end

#   ieaa = Int64[] # indexes of ea to add
#   iead = Int64[] # indexes of ea to delete

#   # simulation time
#   simt = 0.0

#   @inbounds begin

#     # start simulation
#     while true

#       # keep track of time
#       simt += δt

#       # time guard
#       if simt > maxt
#         return e0, e1, el, λs, ea, na, simt
#       end

#       # one time step for all edges alive `ea`
#       for (i,v) in enumerate(ea)

#         λsi = λs[v]
#         lii = li[v]
#         λt  = λsi[lii]

#         # update edge length
#         el[v] += δt
#         li[v] += 1

#         # sample new speciation rates
#         λt1 = rnorm(λt + α*δt, srδt*σλ)
#         push!(λsi, λt1)
#         λm = exp(0.5*(λt + λt1))

#         # if speciation event
#         if divev(λm, δt)

#           # if reached `n` species
#           if n === na

#             # update λs and δt for other lineages
#             for vi in ea[i+1:end]
#               el[vi] += δt
#               λsi = λs[vi]
#               lvi = li[vi]
#               λt  = λsi[lvi]

#               push!(λsi, rnorm(λt + α*δt, srδt*σλ))
#             end

#             # to add
#             if !isempty(ieaa)
#               append!(ea, ieaa)
#               empty!(ieaa)
#             end

#            # to delete
#             if !isempty(iead)
#               deleteat!(ea, iead)
#               empty!(iead)
#             end

#             return e0, e1, el, λs, ea, na, simt
#           end

#           ### add new edges
#           # start node
#           push!(e0, e1[v], e1[v])

#           # end nodes
#           push!(e1, ne + 1, ne + 2)

#           # push to edge length
#           push!(el, 0.0, 0.0)

#           # push speciation vector
#           push!(λs, [λt1], [λt1])

#           # push length of vector
#           push!(li, 1, 1)

#           # to update living edges
#           push!(iead, i)
#           push!(ieaa, ne, ne + 1)

#           # update `na` and `ne`
#           ne += 2
#           na += 1
#         end
#       end

#       # to add
#       if !isempty(ieaa)
#         append!(ea, ieaa)
#         empty!(ieaa)
#       end

#       # to delete
#       if !isempty(iead)
#         deleteat!(ea, iead)
#         empty!(iead)
#       end
#     end
#   end
# end






#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Sample conditional on time
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    sim_cladsbd(t   ::Float64;
                λ0  ::Float64 = 1.0,
                μ0  ::Float64 = 0.5,
                α   ::Float64 = 0.0,
                σλ  ::Float64 = 0.1,
                σμ  ::Float64 = 0.1,
                nlim::Int64   = 10_000,
                init::Symbol  = :crown)

Simulate `cTbd` according to a birth-death clads
conditional in stopping at time `t`.
"""
function sim_cladsbd(t   ::Float64;
                     λ0  ::Float64 = 1.0,
                     μ0  ::Float64 = 0.5,
                     α   ::Float64 = 0.0,
                     σλ  ::Float64 = 0.1,
                     σμ  ::Float64 = 0.1,
                     nlim::Int64   = 10_000,
                     init::Symbol  = :crown)

  if init === :crown
    lλ0 = log(λ0)
    lμ0 = log(μ0)
    d1, nn = _sim_cladsbd(t, lλ0, lμ0, α, σλ, σμ, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    d2, nn = _sim_cladsbd(t, lλ0, lμ0, α, σλ, σμ, nn + 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    tree = cTbd(d1, d2, 0.0, false, lλ0, lμ0)
   elseif init === :stem
    tree, nn = _sim_cladsbd(t, log(λ0), log(μ0), α, σλ, σμ, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_cladsbd(t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 nn  ::Int64,
                 nlim::Int64)

Simulate `cTbd` according to a birth-death clads.
"""
function _sim_cladsbd(t   ::Float64,
                      λt  ::Float64,
                      μt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      nn  ::Int64,
                      nlim::Int64)

  if nn < nlim

    λi = exp(λt)
    μi = exp(μt)
    tw = cbd_wait(λi, μi)

    if tw > t
      return cTbd(t, false, false, λt), nn
    end

    if λorμ(λi, μi)
      nn += 1
      d1, nn = _sim_cladsbd(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ),
                  α, σλ, σμ, nn, nlim)
      d2, nn = _sim_cladsbd(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ),
                  α, σλ, σμ, nn, nlim)
 
      return cTbd(d1, d2, tw, false, false, λt, μt), nn
    else
      return cTbd(tw, true, false, λt, μt), nn
    end
  end

  return cTbd(), nn
end




"""
    _sim_cladsbd_t(t   ::Float64,
                   λt  ::Float64,
                   μt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64)

Simulate `cTbd` according to a birth-death clads for
terminal branches.
"""
function _sim_cladsbd_t(t   ::Float64,
                        λt  ::Float64,
                        μt  ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64,
                        nlim::Int64)

  if isfinite(lr) && nn < nlim

    λi = exp(λt)
    μi = exp(μt)
    tw = cbd_wait(λi, μi)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(iρi * Float64(na)/Float64(na-1))
      end
      if nlr >= lr || lU < nlr
        return cTbd(t, false, false, λt, μt), na, nn, nlr
      else
        return cTbd(), na, nn, NaN
      end
    end

    if λorμ(λi, μi)
      nn += 1
      td1, na, nn, lr =
        _sim_cladsbd_t(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
          α, σλ, σμ, lr, lU, iρi, na, nn, nlim)
      td2, na, nn, lr =
        _sim_cladsbd_t(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
          α, σλ, σμ, lr, lU, iρi, na, nn, nlim)

      return cTbd(td1, td2, tw, false, false, λt, μt), na, nn, lr
    else
      return cTbd(tw, true, false, λt, μt), na, nn, lr
    end
  end

  return cTbd(), na, nn, NaN
end




"""
    _sim_cladsbd_i(t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64,
                   xfs::Vector{Float64})

Simulate `cTbd` according to a birth-death clads.
"""
function _sim_cladsbd_i(t   ::Float64,
                        λt  ::Float64,
                        μt  ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        na  ::Int64,
                        nn  ::Int64,
                        nlim::Int64,
                        xfs::Vector{Float64})

  if nn < nlim

    λi = exp(λt)
    μi = exp(μt)
    tw = cbd_wait(λi, μi)

    if tw > t
      na += 1
      push!(xfs, λt)
      return cTbd(t, false, false, λt, μt), na, nn
    end

    if λorμ(λi, μi)
      nn += 1
      d1, na, nn = 
        _sim_cladsbd_i(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
          α, σλ, σμ, na, nn, nlim, xfs)
      d2, na, nn = 
        _sim_cladsbd_i(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
          α, σλ, σμ, na, nn, nlim, xfs)

      return cTbd(d1, d2, tw, false, false, λt, μt), na, nn
    else
      return cTbd(tw, true, false, λt, μt), na, nn
    end
  end

  return cTbd(), na, nn
end



"""
    _sim_cladsbd_it(t   ::Float64,
                    λt  ::Float64,
                    μt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    iρi ::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

Simulate `cTbd` according to a birth-death clads for
internal terminal branches.
"""
function _sim_cladsbd_it(t   ::Float64,
                         λt  ::Float64,
                         μt  ::Float64,
                         α   ::Float64,
                         σλ  ::Float64,
                         σμ  ::Float64,
                         lr  ::Float64,
                         lU  ::Float64,
                         iρi ::Float64,
                         na  ::Int64,
                         nn  ::Int64,
                         nlim::Int64)

  if lU < lr && nn < nlim

    λi = exp(λt)
    μi = exp(μt)
    tw = cbd_wait(λi, μi)

   if tw > t
      na += 1
      lr += log(iρi)
      return cTbd(t, false, false, λt, μt), na, nn, lr
    end

    if λorμ(λi, μi)
      nn += 1
      td1, na, nn, lr = 
        _sim_cladsbd_it(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
          α, σλ, σμ, lr, lU, iρi, na, nn, nlim)
      td2, na, nn, lr = 
        _sim_cladsbd_it(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
          α, σλ, σμ, lr, lU, iρi, na, nn, nlim)

      return cTbd(td1, td2, tw, false, false, λt, μt), na, nn, lr
    else
      return cTbd(tw, true, false, λt, μt), na, nn, lr
    end
  end

  return cTbd(), na, nn, NaN
end




"""
    _sim_cladsbd_surv(t   ::Float64,
                      λt  ::Float64,
                      μt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      surv::Bool,
                      nn  ::Int64)

Simulate `cTbd` according to a birth-death clads.
"""
function _sim_cladsbd_surv(t   ::Float64,
                           λt  ::Float64,
                           μt  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           σμ  ::Float64,
                           surv::Bool,
                           nn  ::Int64)

  if !surv && nn < 500

    λi = exp(λt)
    μi = exp(μt)
    tw = cbd_wait(λi, μi)

    if tw > t
      return true, nn
    end

    if λorμ(λi, μi)
      nn += 1
      surv, nn = _sim_cladsbd_surv(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
        α, σλ, σμ, surv, nn)
      surv, nn = _sim_cladsbd_surv(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
        α, σλ, σμ, surv, nn)
 
      return surv, nn
    else
      return surv, nn
    end
  end

  return true, nn
end



