#=

Anagenetic GBM birth-death Simulation

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
#     sim_gbmbd(n       ::Int64;
#               λ0      ::Float64 = 1.0,
#               μ0      ::Float64 = 0.1,
#               α       ::Float64 = 0.0,
#               σλ      ::Float64 = 0.1,
#               σμ      ::Float64 = 0.1,
#               start   ::Symbol  = :stem,
#               δt      ::Float64 = 1e-3,
#               nstar   ::Int64   = 2*n,
#               p       ::Float64 = 5.0,
#               warnings::Bool    = true)

# Simulate `iTbdX` according to geometric Brownian motions for birth and death
# rates.
# """
# function sim_gbmbd(n       ::Int64;
#                    λ0      ::Float64 = 1.0,
#                    μ0      ::Float64 = 0.1,
#                    α       ::Float64 = 0.0,
#                    σλ      ::Float64 = 0.1,
#                    σμ      ::Float64 = 0.1,
#                    start   ::Symbol  = :stem,
#                    δt      ::Float64 = 1e-3,
#                    nstar   ::Int64   = 2*n,
#                    p       ::Float64 = 5.0,
#                    warnings::Bool    = true,
#                    maxt    ::Float64 = δt*1e6)

#   # simulate in non-recursive manner
#   e0, e1, el, λs, μs, ea, ee, na, simt =
#     _sedges_gbmbd(nstar, log(λ0), log(μ0), α, σλ, σμ, xt1, βλ, σx, δt, sqrt(δt), start, maxt)

#   if simt >= maxt
#     warnings && @warn "simulation surpassed maximum time"
#   end

#   # transform to iTree
#   t = iTbdX(e0, e1, el, λs, μs, ea, ee, e1[1], 1, δt)

#   if iszero(ntipsalive(t))
#     warnings && @warn "tree went extinct"
#     return t
#   end

#   # sample a time when species(t) == `n`
#   nt = ltt(t)
#   tn = times_n(n, nt)
#   c  = usample(tn, p)

#   if iszero(c)
#     warnings && @warn "tree not sampled, try increasing `p`"
#     return iTbdX(0.0, 0.0, 0.0, false, false, Float64[], Float64[])
#   else
#     # cut the tree
#     t = cutbottom(t, simt - c)
#     return t
#   end
# end





# """
#     _sedges_gbmbd(n    ::Int64,
#                   λ0   ::Float64,
#                   μ0   ::Float64,
#                   α    ::Float64,
#                   σλ   ::Float64,
#                   σμ   ::Float64,
#                   δt   ::Float64,
#                   srδt ::Float64,
#                   start::Symbol)

# Simulate `gbmbd` just until hitting `n` alive species. Note that this is
# a biased sample for a tree conditional on `n` species.
# """
# function _sedges_gbmbd(n    ::Int64,
#                        λ0   ::Float64,
#                        μ0   ::Float64,
#                        α    ::Float64,
#                        σλ   ::Float64,
#                        σμ   ::Float64,
#                        δt   ::Float64,
#                        srδt ::Float64,
#                        start::Symbol,
#                        maxt ::Float64)

#   # edges
#   e0 = Int64[]
#   e1 = Int64[]
#   # edges extinct
#   ee = Int64[]

#   if start == :stem
#     # edges alive
#     ea = [1]
#     # first edge
#     push!(e0, 1)
#     push!(e1, 2)
#     # max index
#     mxi0 = n*2
#     # edge lengths
#     el = [0.0]
#     # lambda and mu vector for each edge
#     λs = [Float64[]]
#     μs = [Float64[]]
#     # starting speciation rate
#     push!(λs[1], λ0)
#     push!(μs[1], μ0)
#     # lastindex for each edge
#     li = [1]

#     na = 1 # current number of alive species
#     ne = 2 # current maximum node number

#   elseif start == :crown
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
#     μs = [Float64[], Float64[], Float64[]]
#     # starting speciation and extinction rate
#     push!(λs[1], λ0, λ0)
#     push!(λs[2], λ0)
#     push!(λs[3], λ0)
#     push!(μs[1], μ0, μ0)
#     push!(μs[2], μ0)
#     push!(μs[3], μ0)
#     # lastindex for each edge
#     li = [2, 1, 1]

#     na = 2 # current number of alive species
#     ne = 4 # current maximum node number

#   else
#     @error "$start does not match stem or crown"
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
#         return e0, e1, el, λs, ea, ee, na, simt
#       end

#       # one time step for all edges alive `ea`
#       for (i,v) in enumerate(ea)

#         λsi = λs[v]
#         μsi = μs[v]
#         lii = li[v]
#         λt  = λsi[lii]
#         μt  = μsi[lii]

#         # update edge length
#         el[v] += δt
#         li[v] += 1

#         # sample new speciation and extinction rates
#         λt1 = rnorm(λt + α*δt, srδt*σλ)
#         μt1 = rnorm(μt, srδt*σμ)
#         push!(λsi, λt1)
#         push!(μsi, μt1)
#         λm = exp(0.5*(λt + λt1))
#         μm = exp(0.5*(μt + μt1))

#         # if diversification event
#         if divev(λm, μm, δt)

#           #if speciation
#           if λorμ(λm, μm)

#             # if reached `n` species
#             if n === na

#               # update λs and δt for other lineages
#               for vi in ea[i+1:end]
#                 el[vi] += δt
#                 λsi = λs[vi]
#                 μsi = μs[vi]
#                 lvi = li[vi]
#                 λt  = λsi[lvi]
#                 μt  = μsi[lvi]

#                 push!(λsi, rnorm(λt + α*δt, srδt*σλ))
#                 push!(μsi, rnorm(μt, srδt*σμ))
#               end

#               # to add
#               if !isempty(ieaa)
#                 append!(ea, ieaa)
#                 empty!(ieaa)
#               end

#              # to delete
#               if !isempty(iead)
#                 deleteat!(ea, iead)
#                 empty!(iead)
#               end

#               return e0, e1, el, λs, μs, ea, ee, na, simt
#             end

#             ### add new edges
#             # start node
#             push!(e0, e1[v], e1[v])

#             # end nodes
#             push!(e1, ne + 1, ne + 2)

#             # push to edge length
#             push!(el, 0.0, 0.0)

#             # push speciation vector
#             push!(λs, [λt1], [λt1])
#             push!(μs, [μt1], [μt1])

#             # push length of vector
#             push!(li, 1, 1)

#             # to update living edges
#             push!(iead, i)
#             push!(ieaa, ne, ne + 1)

#             # update `na` and `ne`
#             ne += 2
#             na += 1

#           #if extinction
#           else
#             # if tree goes extinct
#             if isone(na)
#               # extinct edges
#               push!(ee, v)
#               # delete from alive lineages
#               deleteat!(ea, i)

#               return e0, e1, el, λs, μs, ea, ee, 0, simt
#             end

#             # extinct edges
#             push!(ee, v)
#             # to update alive lineages
#             push!(iead, i)
#             # update number of alive species
#             na -= 1
#           end
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




# #=
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Sample conditional on time
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =#





# """
#     sim_gbmbd(t   ::Float64;
#               λ0  ::Float64 = 1.0,
#               α   ::Float64 = 0.0,
#               σλ  ::Float64 = 0.1,
#               μ   ::Float64 = 0.2,
#               δt  ::Float64 = 1e-3,
#               nlim::Int64   = 10_000,
#               init::Symbol  = :crown)

# Simulate `iTbdX` according to geometric Brownian motions for birth and death
# rates.
# """
# function sim_gbmbd(t   ::Float64;
#                    λ0  ::Float64 = 1.0,
#                    μ0  ::Float64 = 0.2,
#                    α   ::Float64 = 0.0,
#                    σλ  ::Float64 = 0.1,
#                    σμ  ::Float64 = 0.1,
#                    δt  ::Float64 = 1e-3,
#                    nlim::Int64   = 10_000,
#                    init::Symbol  = :crown)

#   if init === :crown
#     lλ0 = log(λ0)
#     lμ0 = log(μ0)
#     d1, nn = _sim_gbmbd(t, lλ0, lμ0, α, σλ, σμ, xt1, βλ, σx, δt, sqrt(δt), 0, 1, nlim)
#     if nn >= nlim
#       @warn "maximum number of lineages surpassed"
#     end

#     d2, nn = _sim_gbmbd(t, lλ0, lμ0, α, σλ, σμ, xt1, βλ, σx, δt, sqrt(δt), 0, 1, nlim)
#     if nn >= nlim
#       @warn "maximum number of lineages surpassed"
#     end

#     tree = iTbdX(d1, d2, 0.0, δt, 0.0, false, false,
#       Float64[lλ0, lλ0], Float64[lμ0, lμ0])

#   elseif init === :stem
#     tree, nn = _sim_gbmbd(t, lλ0, lμ0, α, σλ, σμ, xt1, βλ, σx, δt, sqrt(δt), 0, 1, nlim)

#     if nn >= nlim
#       @warn "maximum number of lineages surpassed"
#     end
#   else
#     @error string(init, " does not match either crown or stem")
#   end

#   return tree
# end




"""
    _sim_gbmbd(t   ::Float64,
               λt  ::Float64,
               μt  ::Float64,
               α   ::Float64,
               σλ  ::Float64,
               σμ  ::Float64,
               xt  ::Float64,
               βλ  ::Float64,
               σx  ::Float64,
               δt  ::Float64,
               srδt::Float64,
               na  ::Int64,
               nn  ::Int64,
               nlim::Int64)

Simulate `iTbdX` according to geometric Brownian motions for birth and death
rates, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmbd(t   ::Float64,
                    λt  ::Float64,
                    μt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    xt  ::Float64,
                    βλ  ::Float64,
                    σx  ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

  if nn < nlim

    xv = Float64[xt]
    λv = Float64[λt]
    μv = Float64[μt]
    bt = 0.0

    while true

      if t <= δt
        bt += t
        t   = max(0.0,t)
        srt = sqrt(t)
        xt1 = rnorm(xt, srt*σx)
        λt1 = rnorm(λt + (α + βλ*xt)*t, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)
        push!(xv, xt1)
        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, t)
          # if speciation
          if λorμ(λm, μm)
            nn += 1
            na  += 2
            return iTbdX(iTbdX(0.0, δt, 0.0, false, false,
                               Float64[λt1, λt1], Float64[μt1, μt1],
                               Float64[xt1, xt1]),
                         iTbdX(0.0, δt, 0.0, false, false,
                               Float64[λt1, λt1], Float64[μt1, μt1],
                               Float64[xt1, xt1]),
                         bt, δt, t, false, false, λv, μv, xv), na, nn
          # if extinction
          else
            return iTbdX(bt, δt, t, true, false, λv, μv, xv), na, nn
          end
        end

        na += 1
        return iTbdX(bt, δt, t, false, false, λv, μv, xv), na, nn
      end

      t  -= δt
      bt += δt
      xt1 = rnorm(xt, srδt*σx)
      λt1 = rnorm(λt + (α + βλ*xt)*δt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)
      push!(λv, λt1)
      push!(μv, μt1)
      push!(xv, xt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        # if speciation
        if λorμ(λm, μm)
          nn += 1
          td1, na, nn =
            _sim_gbmbd(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt, na, nn, nlim)
          td2, na, nn =
            _sim_gbmbd(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt, na, nn, nlim)

          return iTbdX(td1, td2, bt, δt, δt, false, false, λv, μv, xv), na, nn
        # if extinction
        else
          return iTbdX(bt, δt, δt, true, false, λv, μv, xv), na, nn
        end
      end

      xt = xt1
      λt = λt1
      μt = μt1
    end
  end

  return iTbdX(), na, nn
end




"""
    _sim_gbmbd_t(t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 xt  ::Float64,
                 βλ  ::Float64,
                 σx  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 Iρi ::Float64,
                 na  ::Int64,
                 nn  ::Int64,
                 nlim::Int64,
                 xist::Vector{Float64},
                 xfst::Vector{Float64},
                 est ::Vector{Float64})

Simulate `iTbdX` according to geometric Brownian motions for birth and death
rates, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmbd_t(t   ::Float64,
                      λt  ::Float64,
                      μt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      xt  ::Float64,
                      βλ  ::Float64,
                      σx  ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      lr  ::Float64,
                      lU  ::Float64,
                      Iρi ::Float64,
                      na  ::Int64,
                      nn  ::Int64,
                      nlim::Int64,
                      xist::Vector{Float64},
                      xfst::Vector{Float64},
                      est ::Vector{Float64})

  if isfinite(lr) && nn < nlim

    xv = Float64[xt]
    λv = Float64[λt]
    μv = Float64[μt]
    bt = 0.0

    while true

      if t <= δt
        bt += t
        t   = max(0.0,t)
        srt = sqrt(t)
        xt1 = rnorm(xt, srt*σx)
        λt1 = rnorm(λt + (α + βλ*xt)*t, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)
        push!(xv, xt1)
        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, t)
          # if speciation
          if λorμ(λm, μm)
            nn += 1
            na += 2
            if na === 2
              nlr = lr + log(Iρi*2.0)
            else
              nlr = lr + log(Iρi * Iρi * Float64(na)/Float64(na-2))
            end
            if nlr < lr && lU >= nlr
              return iTbdX(), na, nn, NaN
            else
              push!(xist, xv[1], xv[1])
              push!(xfst, xt1,   xt1)
              push!(est,  bt,    bt)

              return iTbdX(iTbdX(0.0, δt, 0.0, false, false,
                                 Float64[λt1, λt1], Float64[μt1, μt1],
                                 Float64[xt1, xt1]),
                           iTbdX(0.0, δt, 0.0, false, false,
                                Float64[λt1, λt1], Float64[μt1, μt1],
                                Float64[xt1, xt1]),
                           bt, δt, t, false, false, λv, μv, xv), na, nn, nlr
            end
          # if extinction
          else
            return iTbdX(bt, δt, t, true, false, λv, μv, xv), na, nn, lr
          end
        end

        na += 1
        nlr = lr
        if na > 1
          nlr += log(Iρi * Float64(na)/Float64(na-1))
        end
        if nlr < lr && lU >= nlr
          return iTbdX(), na, nn, NaN
        else
          push!(xist, xv[1])
          push!(xfst, xt1)
          push!(est,  bt)

          return iTbdX(bt, δt, t, false, false, λv, μv, xv), na, nn, nlr
        end
      end

      t  -= δt
      bt += δt
      xt1 = rnorm(xt, srδt*σx)
      λt1 = rnorm(λt + (α + βλ*xt)*δt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)
      push!(λv, λt1)
      push!(μv, μt1)
      push!(xv, xt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        # if speciation
        if λorμ(λm, μm)
          nn += 1
          td1, na, nn, lr =
            _sim_gbmbd_t(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
              lr, lU, Iρi, na, nn, nlim, xist, xfst, est)
          td2, na, nn, lr =
            _sim_gbmbd_t(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
              lr, lU, Iρi, na, nn, nlim, xist, xfst, est)

          return iTbdX(td1, td2, bt, δt, δt, false, false, λv, μv, xv),
                 na, nn, lr
        # if extinction
        else
          return iTbdX(bt, δt, δt, true, false, λv, μv, xv), na, nn, lr
        end
      end

      xt = xt1
      λt = λt1
      μt = μt1
    end
  end

  return iTbdX(), na, nn, NaN
end




"""
    _sim_gbmbd_i(t   ::Float64,
                 λt  ::Float64,
                 μt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 xt  ::Float64,
                 βλ  ::Float64,
                 σx  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 nn  ::Int64,
                 nlim::Int64,
                 λsp ::Vector{Float64},
                 μsp ::Vector{Float64},
                 xsp ::Vector{Float64})

Simulate `iTbdX` according to geometric Brownian motions for birth and death
rates, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmbd_i(t   ::Float64,
                      λt  ::Float64,
                      μt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      xt  ::Float64,
                      βλ  ::Float64,
                      σx  ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      nn  ::Int64,
                      nlim::Int64,
                      λsp ::Vector{Float64},
                      μsp ::Vector{Float64},
                      xsp ::Vector{Float64})

  if nn < nlim

    xv = Float64[xt]
    λv = Float64[λt]
    μv = Float64[μt]
    bt = 0.0

    while true

      if t <= δt
        bt += t
        t   = max(0.0,t)
        srt = sqrt(t)
        xt1 = rnorm(xt, srt*σx)
        λt1 = rnorm(λt + (α + βλ*xt)*t, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)
        push!(xv, xt1)
        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, t)
          # if speciation
          if λorμ(λm, μm)
            nn += 1
            push!(λsp, λt1, λt1)
            push!(μsp, μt1, μt1)
            push!(xsp, xt1, xt1)

            return iTbdX(iTbdX(0.0, δt, 0.0, false, false,
                               Float64[λt1, λt1], Float64[μt1, μt1],
                               Float64[xt1, xt1]),
                         iTbdX(0.0, δt, 0.0, false, false,
                               Float64[λt1, λt1], Float64[μt1, μt1],
                               Float64[xt1, xt1]),
                         bt, δt, t, false, false, λv, μv, xv), nn
          # if extinction
          else
            return iTbdX(bt, δt, t, true, false, λv, μv, xv), nn
          end
        end
        push!(λsp, λt1)
        push!(μsp, μt1)
        push!(xsp, xt1)

        return iTbdX(bt, δt, t, false, false, λv, μv, xv), nn
      end

      t  -= δt
      bt += δt
      xt1 = rnorm(xt, srδt*σx)
      λt1 = rnorm(λt + (α + βλ*xt)*δt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)
      push!(λv, λt1)
      push!(μv, μt1)
      push!(xv, xt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        # if speciation
        if λorμ(λm, μm)
          nn += 1
          td1, nn =
            _sim_gbmbd_i(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
              nn, nlim, λsp, μsp, xsp)
          td2, nn =
            _sim_gbmbd_i(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
              nn, nlim, λsp, μsp, xsp)

          return iTbdX(td1, td2, bt, δt, δt, false, false, λv, μv, xv), nn
        # if extinction
        else
          return iTbdX(bt, δt, δt, true, false, λv, μv, xv), nn
        end
      end

      xt = xt1
      λt = λt1
      μt = μt1
    end
  end

  return iTbdX(), nn
end




"""
    _sim_gbmbd_it(nsδt::Float64,
                  t   ::Float64,
                  λt  ::Float64,
                  μt  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  xt  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  lr  ::Float64,
                  lU  ::Float64,
                  Iρi ::Float64,
                  na  ::Int64,
                  nn  ::Int64,
                  nlim::Int64)

Simulate `iTbdX` according to geometric Brownian motions for birth and death
rates, starting with a non-standard δt with a limit in the number of species.
"""
function _sim_gbmbd_it(nsδt::Float64,
                       t   ::Float64,
                       λt  ::Float64,
                       μt  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       xt  ::Float64,
                       βλ  ::Float64,
                       σx  ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       Iρi ::Float64,
                       na  ::Int64,
                       nn  ::Int64,
                       nlim::Int64)

  xv = Float64[xt]
  λv = Float64[λt]
  μv = Float64[μt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt += t
    t   = max(0.0,t)
    srt = sqrt(t)
    xt1 = rnorm(xt, srt*σx)
    λt1 = rnorm(λt + (α + βλ*xt)*t, srt*σλ)
    μt1 = rnorm(μt, srt*σμ)
    push!(xv, xt1)
    push!(λv, λt1)
    push!(μv, μt1)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if divev(λm, μm, t)
      # if speciation
      if λorμ(λm, μm)
        nn += 1
        na += 2
        lr += 2.0*log(Iρi)
        return iTbdX(iTbdX(0.0, δt, 0.0, false, false,
                           Float64[λt1, λt1], Float64[μt1, μt1],
                           Float64[xt1, xt1]),
                     iTbdX(0.0, δt, 0.0, false, false,
                           Float64[λt1, λt1], Float64[μt1, μt1],
                           Float64[xt1, xt1]),
                     bt, δt, t, false, false, λv, μv, xv), na, nn, lr
      # if extinction
      else
        return iTbdX(bt, δt, t, true, false, λv, μv, xv), na, nn, lr
      end
    end

    na += 1
    lr += log(Iρi)
    return iTbdX(bt, δt, t, false, false, λv, μv, xv), na, nn, lr
  end

  t  -= nsδt
  bt += nsδt
  srnsδt = sqrt(nsδt)
  xt1 = rnorm(xt, srnsδt*σx)
  λt1 = rnorm(λt + (α + βλ*xt)*nsδt, srnsδt*σλ)
  μt1 = rnorm(μt, srnsδt*σμ)
  push!(λv, λt1)
  push!(μv, μt1)
  push!(xv, xt1)

  λm = exp(0.5*(λt + λt1))
  μm = exp(0.5*(μt + μt1))

  if divev(λm, μm, nsδt)
    # if speciation
    if λorμ(λm, μm)
      nn += 1
      td1, na, nn, lr =
        _sim_gbmbd_it(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
          lr, lU, Iρi, na, nn, nlim)
      td2, na, nn, lr =
        _sim_gbmbd_it(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
          lr, lU, Iρi, na, nn, nlim)

      return iTbdX(td1, td2, bt, δt, nsδt, false, false, λv, μv, xv), na, nn, lr
    else
    # if extinction
      return iTbdX(bt, δt, nsδt, true, false, λv, μv, xv), na, nn, lr
    end
  end

  xt = xt1
  λt = λt1
  μt = μt1

  if lU < lr && nn < nlim

    ## second: standard δt
    while true

      if t <= δt
        bt += t
        t   = max(0.0,t)
        srt = sqrt(t)
        xt1 = rnorm(xt, srt*σx)
        λt1 = rnorm(λt + (α + βλ*xt)*t, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)
        push!(xv, xt1)
        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, t)
          # if speciation
          if λorμ(λm, μm)
            nn += 1
            na += 2
            lr += 2.0*log(Iρi)
            return iTbdX(iTbdX(0.0, δt, 0.0, false, false,
                             Float64[λt1, λt1], Float64[μt1, μt1],
                             Float64[xt1, xt1]),
                        iTbdX(0.0, δt, 0.0, false, false,
                             Float64[λt1, λt1], Float64[μt1, μt1],
                             Float64[xt1, xt1]),
                        bt, δt, t, false, false, λv, μv, xv), na, nn, lr
          # if extinction
          else
            return iTbdX(bt, δt, t, true, false, λv, μv, xv), na, nn, lr
          end
        end

        na += 1
        lr += log(Iρi)
        return iTbdX(bt, δt, t, false, false, λv, μv, xv), na, nn, lr
      end

      t  -= δt
      bt += δt
      xt1 = rnorm(xt, srδt*σx)
      λt1 = rnorm(λt + (α + βλ*xt)*δt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)
      push!(λv, λt1)
      push!(μv, μt1)
      push!(xv, xt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        # if speciation
        if λorμ(λm, μm)
          nn += 1
          td1, na, nn, lr =
            _sim_gbmbd_it(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
              lr, lU, Iρi, na, nn, nlim)
          td2, na, nn, lr =
            _sim_gbmbd_it(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
              lr, lU, Iρi, na, nn, nlim)

          return iTbdX(td1, td2, bt, δt, δt, false, false, λv, μv, xv), 
                 na, nn, lr
        # if extinction
        else
          return iTbdX(bt, δt, δt, true, false, λv, μv, xv), na, nn, lr
        end
      end

      xt = xt1
      λt = λt1
      μt = μt1
    end
  end

  return iTbdX(), na, nn, NaN
end




"""
    _sim_gbmbd_it(t   ::Float64,
                  λt  ::Float64,
                  μt  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  xt  ::Float64,
                  βλ  ::Float64,
                  σx  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  lr  ::Float64,
                  lU  ::Float64,
                  Iρi ::Float64,
                  na  ::Int64,
                  nn  ::Int64,
                  nlim::Int64)

Simulate `iTbdX` according to geometric Brownian motions for birth and death
rates, starting with a non-standard δt with a limit in the number of species.
"""
function _sim_gbmbd_it(t   ::Float64,
                       λt  ::Float64,
                       μt  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       xt  ::Float64,
                       βλ  ::Float64,
                       σx  ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       Iρi ::Float64,
                       na  ::Int64,
                       nn  ::Int64,
                       nlim::Int64)

  if lU < lr && nn < nlim

    xv = Float64[xt]
    λv = Float64[λt]
    μv = Float64[μt]
    bt = 0.0

    ## second: standard δt
    while true

      if t <= δt
        bt += t
        t   = max(0.0,t)
        srt = sqrt(t)
        xt1 = rnorm(xt, srt*σx)
        λt1 = rnorm(λt + (α + βλ*xt)*t, srt*σλ)
        μt1 = rnorm(μt, srt*σμ)
        push!(xv, xt1)
        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if divev(λm, μm, t)
          # if speciation
          if λorμ(λm, μm)
            nn += 1
            na += 2
            lr += 2.0*log(Iρi)
            return iTbdX(iTbdX(0.0, δt, 0.0, false, false,
                             Float64[λt1, λt1], Float64[μt1, μt1],
                             Float64[xt1, xt1]),
                        iTbdX(0.0, δt, 0.0, false, false,
                             Float64[λt1, λt1], Float64[μt1, μt1],
                             Float64[xt1, xt1]),
                        bt, δt, t, false, false, λv, μv, xv), na, nn, lr
          # if extinction
          else
            return iTbdX(bt, δt, t, true, false, λv, μv, xv), na, nn, lr
          end
        end

        na += 1
        lr += log(Iρi)
        return iTbdX(bt, δt, t, false, false, λv, μv, xv), na, nn, lr
      end

      t  -= δt
      bt += δt
      xt1 = rnorm(xt, srδt*σx)
      λt1 = rnorm(λt + (α + βλ*xt)*δt, srδt*σλ)
      μt1 = rnorm(μt, srδt*σμ)
      push!(λv, λt1)
      push!(μv, μt1)
      push!(xv, xt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        # if speciation
        if λorμ(λm, μm)
          nn += 1
          td1, na, nn, lr =
            _sim_gbmbd_it(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
              lr, lU, Iρi, na, nn, nlim)
          td2, na, nn, lr =
            _sim_gbmbd_it(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt,
              lr, lU, Iρi, na, nn, nlim)

          return iTbdX(td1, td2, bt, δt, δt, false, false, λv, μv, xv), na, nn, lr
        # if extinction
        else
          return iTbdX(bt, δt, δt, true, false, λv, μv, xv), na, nn, lr
        end
      end

      xt = xt1
      λt = λt1
      μt = μt1
    end
  end

  return iTbdX(), na, nn, NaN
end




# """
#     _sim_gbmbd_surv(t   ::Float64,
#                     λt  ::Float64,
#                     μt  ::Float64,
#                     α   ::Float64,
#                     σλ  ::Float64,
#                     σμ  ::Float64,
#                     δt  ::Float64,
#                     srδt::Float64,
#                     surv::Bool,
#                     nn  ::Int64)

# Simulate `iTbdX` according to geometric Brownian motions for birth and death
# rates, with a limit on the number lineages allowed to reach.
# """
# function _sim_gbmbd_surv(t   ::Float64,
#                          λt  ::Float64,
#                          μt  ::Float64,
#                          α   ::Float64,
#                          σλ  ::Float64,
#                          σμ  ::Float64,
#                          δt  ::Float64,
#                          srδt::Float64,
#                          surv::Bool,
#                          nn  ::Int64)

#   if !surv && nn < 200

#     while true

#       if t <= δt

#         t   = max(0.0,t)
#         μt1 = rnorm(μt, sqrt(t)*σμ)
#         μm  = exp(0.5*(μt + μt1))

#         # if extinction
#         if rand() < exp(0.5*(μt + μt1))*t
#           return surv, nn
#         else
#           return true, nn
#         end

#         return true, nn
#       end

#       t  -= δt

#       λt1 = rnorm(λt + α*δt, srδt*σλ)
#       μt1 = rnorm(μt, srδt*σμ)

#       λm = exp(0.5*(λt + λt1))
#       μm = exp(0.5*(μt + μt1))

#       if divev(λm, μm, δt)
#         # if speciation
#         if λorμ(λm, μm)
#           nn += 1
#           surv, nn =
#             _sim_gbmbd_surv(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt, surv, nn)
#           surv, nn =
#             _sim_gbmbd_surv(t, λt1, μt1, α, σλ, σμ, xt1, βλ, σx, δt, srδt, surv, nn)

#           return surv, nn
#         # if extinction
#         else
#           return surv, nn
#         end
#       end

#       xt = xt1
#       λt = λt1
#       μt = μt1
#     end
#   end

#   return true, nn
# end




