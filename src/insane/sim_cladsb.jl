#=

Anagenetic GBM pure-birth Simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Sample conditional on number of species
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#



# """
#     sim_cladsb(n       ::Int64;
#               λ0      ::Float64 = 1.0,
#               α       ::Float64 = 0.0,
#               σλ      ::Float64 = 0.1,
#               δt      ::Float64 = 1e-3,
#               init    ::Symbol  = :stem,
#               nstar   ::Int64   = n + 2,
#               p       ::Float64 = 5.0,
#               warnings::Bool    = true,
#               maxt    ::Float64 = δt*1e7)

# Simulate `cTb` according to a pure-birth clads.
# """
# function sim_cladsb(n       ::Int64;
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
#     _sedges_cladsb(nstar, log(λ0), α, σλ, δt, sqrt(δt), init, maxt)

#   if simt >= maxt
#     warnings && @warn "simulation surpassed maximum time"
#     return cTb()
#   end

#   # transform to iTree
#   t = cTb(e0, e1, el, λs, ea, e1[1], 1, δt)

#   # sample a time when species(t) == `n`
#   nt = ltt(t)
#   tn = times_n(n, nt)
#   c  = usample(tn, p)

#   if iszero(c)
#     warnings && @warn "tree not sampled, try increasing `p`"
#     return cTb()
#   else
#     # cut the tree
#     t = cutbottom(t, simt - c)
#     return t
#   end
# end





# """
#     _sedges_cladsb(n    ::Int64,
#                   λ0   ::Float64,
#                   α    ::Float64,
#                   σλ   ::Float64,
#                   δt   ::Float64,
#                   srδt ::Float64,
#                   init::Symbol)

# Simulate `cladsb` just until hitting `n` alive species. Note that this is
# a biased sample for a tree conditional on `n` species.
# """
# function _sedges_cladsb(n    ::Int64,
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
    sim_cladsb(t   ::Float64;
                λ0  ::Float64 = 1.0,
                α   ::Float64 = 0.0,
                σλ  ::Float64 = 0.1,
                nlim::Int64   = 10_000,
                init::Symbol  = :crown)

Simulate `cTb` according to a pure-birth clads
conditional in stopping at time `t`.
"""
function sim_cladsb(t   ::Float64;
                     λ0  ::Float64 = 1.0,
                     α   ::Float64 = 0.0,
                     σλ  ::Float64 = 0.1,
                     nlim::Int64   = 10_000,
                     init::Symbol  = :crown)

  if init === :crown
    lλ0 = log(λ0)
    d1, nn = _sim_cladsb(t, lλ0, α, σλ, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    d2, nn = _sim_cladsb(t, lλ0, α, σλ, nn + 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    tree = cTb(d1, d2, 0.0, false, lλ0)
   elseif init === :stem
    tree, nn = _sim_cladsb(t, log(λ0), α, σλ, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_cladsb(t   ::Float64,
                 λt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 nn  ::Int64,
                 nlim::Int64)


Simulate `cTb` according to a pure-birth clads.
"""
function _sim_cladsb(t   ::Float64,
                      λt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      nn  ::Int64,
                      nlim::Int64)

  if nn < nlim

    tw = cb_wait(exp(λt))

    if tw > t
      return cTb(t, false, λt), nn
    end

    nn += 1
    d1, nn = _sim_cladsb(t - tw, rnorm(λt + α, σλ), α, σλ, nn, nlim)
    d2, nn = _sim_cladsb(t - tw, rnorm(λt + α, σλ), α, σλ, nn, nlim)
 
    return cTb(d1, d2, tw, false, λt), nn
  end

  return cTb(), nn
end





"""
    _sim_cladsb_t(t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn ::Int64,
                   nlim::Int64)

Simulate `cTb` according to a pure-birth clads for
terminal branches.
"""
function _sim_cladsb_t(t   ::Float64,
                        λt  ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64,
                        nlim::Int64)

  if isfinite(lr) && nn < nlim

    tw = cb_wait(exp(λt))

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(iρi * Float64(na)/Float64(na-1))
      end
      if nlr >= lr || lU < nlr
        return cTb(t, false, λt), na, nn, nlr
      else
        return cTb(), na, nn, NaN
      end
    end

    nn += 1
    td1, na, nn, lr =
      _sim_cladsb_t(t - tw, rnorm(λt + α, σλ), α, σλ, 
        lr, lU, iρi, na, nn, nlim)
    td2, na, nn, lr =
      _sim_cladsb_t(t - tw, rnorm(λt + α, σλ), α, σλ, 
        lr, lU, iρi, na, nn, nlim)

    return cTb(td1, td2, tw, false, λt), na, nn, lr
  end

  return cTb(), na, nn, NaN
end




"""
    _sim_cladsb_i(t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   nn  ::Int64,
                   nlim::Int64,
                   xfs::Vector{Float64})

Simulate `cTb` according to a pure-birth clads.
"""
function _sim_cladsb_i(t   ::Float64,
                        λt  ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        nn  ::Int64,
                        nlim::Int64,
                        xfs::Vector{Float64})

  if nn < nlim

    tw = cb_wait(exp(λt))

    if tw > t
      push!(xfs, λt)
      return cTb(t, false, λt), nn
    end

    nn += 1
    d1, nn = _sim_cladsb_i(t - tw, rnorm(λt + α, σλ), α, σλ, nn, nlim, xfs)
    d2, nn = _sim_cladsb_i(t - tw, rnorm(λt + α, σλ), α, σλ, nn, nlim, xfs)
 
    return cTb(d1, d2, tw, false, λt), nn
  end

  return cTb(), nn
end




"""
    _sim_cladsb_it(t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    iρi ::Float64,
                    nn  ::Int64,
                    nlim::Int64)

Simulate `cTb` according to a pure-birth clads for
internal terminal branches.
"""
function _sim_cladsb_it(t   ::Float64,
                         λt  ::Float64,
                         α   ::Float64,
                         σλ  ::Float64,
                         lr  ::Float64,
                         lU  ::Float64,
                         iρi ::Float64,
                         nn  ::Int64,
                         nlim::Int64)

  if lU < lr && nn < nlim

    tw = cb_wait(exp(λt))

    if tw > t
      lr += log(iρi)
      return cTb(t, false, λt), nn, lr
    end

    nn += 1
    td1, nn, lr = 
      _sim_cladsb_it(t - tw, rnorm(λt + α, σλ), α, σλ, lr, lU, iρi, nn, nlim)
    td2, nn, lr = 
      _sim_cladsb_it(t - tw, rnorm(λt + α, σλ), α, σλ, lr, lU, iρi, nn, nlim)

    return cTb(td1, td2, tw, false, λt), nn, lr
  end

  return cTb(), nn, NaN
end




