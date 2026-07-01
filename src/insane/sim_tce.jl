#=

Anagenetic trait pure-birth simulation

Ignacio Quintero Mächler

t(-_-t)

Created 19 01 2026
=#




# #=
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Sample conditional on number of species
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =#



# """
#     sim_tb(n       ::Int64;
#               λ0      ::Float64 = 1.0,
#               α       ::Float64 = 0.0,
#               σλ      ::Float64 = 0.1,
#               δt      ::Float64 = 1e-3,
#               init    ::Symbol  = :stem,
#               nstar   ::Int64   = n + 2,
#               p       ::Float64 = 5.0,
#               warnings::Bool    = true,
#               maxt    ::Float64 = δt*1e7)

# Simulate `iTxb` according to a pure-birth geometric Brownian motion.
# """
# function sim_tb(n       ::Int64;
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
#     _sedges_tb(nstar, log(λ0), α, σλ, δt, sqrt(δt), init, maxt)

#   if simt >= maxt
#     warnings && @warn "simulation surpassed maximum time"
#     return iTxb()
#   end

#   # transform to iTree
#   t = iTxb(e0, e1, el, λs, ea, e1[1], 1, δt)

#   # sample a time when species(t) == `n`
#   nt = ltt(t)
#   tn = times_n(n, nt)
#   c  = usample(tn, p)

#   if iszero(c)
#     warnings && @warn "tree not sampled, try increasing `p`"
#     return iTxb()
#   else
#     # cut the tree
#     t = cutbottom(t, simt - c)
#     return t
#   end
# end





# """
#     _sedges_tb(n    ::Int64,
#                   λ0   ::Float64,
#                   α    ::Float64,
#                   σλ   ::Float64,
#                   δt   ::Float64,
#                   srδt ::Float64,
#                   init::Symbol)

# Simulate `gbmb` just until hitting `n` alive species. Note that this is
# a biased sample for a tree conditional on `n` species.
# """
# function _sedges_tb(n    ::Int64,
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
    sim_tb(t   ::Float64;
           x0  ::Float64 = 0.0,
           σ20 ::Float64 = 0.1,
           ασ  ::Float64 = 0.0,
           σσ  ::Float64 = 0.1,
           λ0  ::Float64 = 1.0,
           αλ  ::Float64 = 0.0,
           βλ  ::Float64 = 0.0,
           σλ  ::Float64 = 0.1,
           δt  ::Float64 = 1e-3,
           nlim::Int64   = 10_000,
           init::Symbol  = :crown)

Simulate `iTxb` according to a trait dependent pure-birth geometric 
Brownian motion conditional in stopping at time `t`.
"""
function sim_tb(t   ::Float64;
                x0  ::Float64 = 0.0,
                σ20 ::Float64 = 0.1,
                ασ  ::Float64 = 0.0,
                σσ  ::Float64 = 0.1,
                λ0  ::Float64 = 1.0,
                αλ  ::Float64 = 0.0,
                βλ  ::Float64 = 0.0,
                σλ  ::Float64 = 0.1,
                δt  ::Float64 = 1e-3,
                nlim::Int64   = 10_000,
                init::Symbol  = :crown)

  if init === :crown
    lλ0  = log(λ0)
    lσ20 = log(σ20)
    d1, nn = _sim_tb(t, x0, lσ20, ασ, σσ, lλ0, αλ, βλ, σλ, 
               δt, sqrt(δt), 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    d2, nn = _sim_tb(t, x0, lσ20, ασ, σσ, lλ0, αλ, βλ, σλ,  
               δt, sqrt(δt), nn + 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    tree = iTxb(d1, d2, 0.0, δt, 0.0, false, [lλ0, lλ0], [x0, x0], [lσ20, lσ20])

   elseif init === :stem
    tree, nn = _sim_tb(t, x0, lσ20, ασ, σσ, log(λ0), αλ, βλ, σλ, 
                 δt, sqrt(δt), 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_tb(t   ::Float64,
            xt  ::Float64,
            lσ2t::Float64,
            ασ  ::Float64,
            σσ  ::Float64,
            lλt ::Float64,
            αλ  ::Float64,
            βλ  ::Float64,
            σλ  ::Float64,
            δt  ::Float64,
            srδt::Float64,
            nn  ::Int64,
            nlim::Int64)

Simulate `iTxb` according to a trait dependent pure-birth 
geometric Brownian motion.
"""
function _sim_tb(t   ::Float64,
                 xt  ::Float64,
                 lσ2t::Float64,
                 ασ  ::Float64,
                 σσ  ::Float64,
                 lλt ::Float64,
                 αλ  ::Float64,
                 βλ  ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 nn  ::Int64,
                 nlim::Int64)

  if nn < nlim

    lλv = Float64[lλt]
    xv  = Float64[xt]
    lσ2 = Float64[lσ2t]
    bt  = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t

        srt = sqrt(t)

        # draw trait rate
        lσ2t1 = rnorm(lσ2t + ασ*t, srt * σσ)
        push!(lσ2, lσ2t1)

        # draw new trait
        xt1 = rnorm(xt, srt * exp(0.25*(lσ2t + lσ2t1)))
        push!(xv, xt1)

        # draw speciation rates
        lλt1 = rnorm(lλt + αλ*t + βλ*(xt1 - xt), srt*σλ)
        push!(lλv, lλt1)

        λm = exp(0.5*(lλt + lλt1))

        if divev(λm, t)
          nn += 1
          return iTxb(iTxb(0.0, δt, 0.0, false, 
                           [lλt1, lλt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                      iTxb(0.0, δt, 0.0, false, 
                           [lλt1, lλt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                      bt, δt, t, false, lλv, xv, lσ2), nn
        end

        return iTxb(bt, δt, t, false, lλv, xv, lσ2), nn
      end

      t  -= δt
      bt += δt

      # draw trait rate
      lσ2t1 = rnorm(lσ2t + ασ*δt, srδt * σσ)
      push!(lσ2, lσ2t1)

      # draw new trait
      xt1 = rnorm(xt, srδt * exp(0.25*(lσ2t + lσ2t1)))
      push!(xv, xt1)

      # draw speciation rates
      lλt1 = rnorm(lλt + αλ*δt + βλ*(xt1 - xt), srδt*σλ)
      push!(lλv, lλt1)

      λm = exp(0.5*(lλt + lλt1))

      if divev(λm, δt)
        nn += 1
        td1, nn = _sim_tb(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
          δt, srδt, nn, nlim)
        td2, nn = _sim_tb(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
          δt, srδt, nn, nlim)

        return iTxb(td1, td2, bt, δt, δt, false, lλv, xv, lσ2), nn
      end

      lλt  = lλt1
      xt   = xt1
      lσ2t = lσ2t1
    end
  end

  return iTxb(), nn
end




"""
    _sim_tb_t(t   ::Float64,
              xt  ::Float64,
              lσ2t::Float64,
              ασ  ::Float64,
              σσ  ::Float64,
              lλt ::Float64,
              αλ  ::Float64,
              βλ  ::Float64,
              σλ  ::Float64,
              δt  ::Float64,
              srδt::Float64,
              lr  ::Float64,
              lU  ::Float64,
              iρi ::Float64,
              na  ::Int64,
              nn  ::Int64,
              nlim::Int64)

Simulate `iTxb` according to a pure-birth geometric Brownian motion for
terminal branches.
"""
function _sim_tb_t(t   ::Float64,
                   xt  ::Float64,
                   lσ2t::Float64,
                   ασ  ::Float64,
                   σσ  ::Float64,
                   lλt ::Float64,
                   αλ  ::Float64,
                   βλ  ::Float64,
                   σλ  ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64)

  if isfinite(lr) && nn < nlim

    lλv = Float64[lλt]
    xv  = Float64[xt]
    lσ2 = Float64[lσ2t]
    bt  = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t

        srt = sqrt(t)

        # draw trait rate
        lσ2t1 = rnorm(lσ2t + ασ*t, srt * σσ)
        push!(lσ2, lσ2t1)

        # draw new trait
        xt1 = rnorm(xt, srt * exp(0.25*(lσ2t + lσ2t1)))
        push!(xv, xt1)

        # draw speciation rates
        lλt1 = rnorm(lλt + αλ*t + βλ*(xt1 - xt), srt*σλ)
        push!(lλv, lλt1)

        λm = exp(0.5*(lλt + lλt1))

        if divev(λm, t)
          nn += 1
          na += 2
          if na === 2
            nlr = lr + log(iρi*2.0)
          else
            nlr = lr + log(iρi * iρi * Float64(na)/Float64(na-2))
          end
          if nlr < lr && lU >= nlr
            return iTxb(), na, nn, NaN
          else
            return iTxb(iTxb(0.0, δt, 0.0, false, 
                             [λt1, λt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                        iTxb(0.0, δt, 0.0, false, 
                             [λt1, λt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                        bt, δt, t, false, lλv, xv, lσ2), na, nn, nlr
          end
        else
          na += 1
          nlr = lr
          if na > 1
            nlr += log(iρi * Float64(na)/Float64(na-1))
          end
          if nlr >= lr
            return iTxb(bt, δt, t, false, lλv, xv, lσ2), na, nn, nlr
          elseif lU < nlr
            return iTxb(bt, δt, t, false, lλv, xv, lσ2), na, nn, nlr
          else
            return iTxb(), na, nn, NaN
          end
        end
      end

      t  -= δt
      bt += δt

      # draw trait rate
      lσ2t1 = rnorm(lσ2t + ασ*δt, srδt * σσ)
      push!(lσ2, lσ2t1)

      # draw new trait
      xt1 = rnorm(xt, srδt * exp(0.25*(lσ2t + lσ2t1)))
      push!(xv, xt1)

      # draw speciation rates
      lλt1 = rnorm(lλt + αλ*δt + βλ*(xt1 - xt), srδt*σλ)
      push!(lλv, lλt1)

      λm = exp(0.5*(lλt + lλt1))

      if divev(λm, δt)
        nn += 1
        td1, na, nn, lr =
          _sim_tb_t(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
            δt, srδt, lr, lU, iρi, na, nn, nlim)
        td2, na, nn, lr =
          _sim_tb_t(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
            δt, srδt, lr, lU, iρi, na, nn, nlim)

        return iTxb(td1, td2, bt, δt, δt, false, lλv, xv, lσ2), na, nn, lr
      end

      lλt  = lλt1
      xt   = xt1
      lσ2t = lσ2t1
    end
  end

  return iTxb(), na, nn, NaN
end




"""
    _sim_tb_it(nsδt::Float64,
               t   ::Float64,
               xt  ::Float64,
               lσ2t::Float64,
               ασ  ::Float64,
               σσ  ::Float64,
               lλt ::Float64,
               αλ  ::Float64,
               βλ  ::Float64,
               σλ  ::Float64,
               δt  ::Float64,
               srδt::Float64,
               lr  ::Float64,
               lU  ::Float64,
               iρi ::Float64,
               nn  ::Int64,
               nlim::Int64)

Simulate `iTxb` according to a pure-birth geometric Brownian motion,
starting with a non-standard `δt` with a limit in the number of species.
"""
function _sim_tb_it(nsδt::Float64,
                    t   ::Float64,
                    xt  ::Float64,
                    lσ2t::Float64,
                    ασ  ::Float64,
                    σσ  ::Float64,
                    lλt ::Float64,
                    αλ  ::Float64,
                    βλ  ::Float64,
                    σλ  ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    iρi ::Float64,
                    nn  ::Int64,
                    nlim::Int64)

  lλv = Float64[lλt]
  xv  = Float64[xt]
  lσ2 = Float64[lσ2t]
  bt  = 0.0

  ## first: non-standard δt
  if t <= nsδt + accerr
    t   = isapprox(t, 0.0) ? 0.0 : isapprox(t, nsδt) ? nsδt : t
    bt += t

    srt = sqrt(t)

    # draw trait rate
    lσ2t1 = rnorm(lσ2t + ασ*t, srt * σσ)
    push!(lσ2, lσ2t1)

    # draw new trait
    xt1 = rnorm(xt, srt * exp(0.25*(lσ2t + lσ2t1)))
    push!(xv, xt1)

    # draw speciation rates
    lλt1 = rnorm(lλt + αλ*t + βλ*(xt1 - xt), srt*σλ)
    push!(lλv, lλt1)

    λm = exp(0.5*(lλt + lλt1))

    if divev(λm, t)
      nn += 1
      lr += 2.0*log(iρi)
      return iTxb(iTxb(0.0, δt, 0.0, false, 
                       [lλt1, lλt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                  iTxb(0.0, δt, 0.0, false, 
                       [lλt1, lλt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                  bt, δt, t, false, lλv, xv, lσ2), nn, lr
    else
      lr += log(iρi)
      return iTxb(bt, δt, t, false, lλv, xv, lσ2), nn, lr
    end
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

 # draw trait rate
  lσ2t1 = rnorm(lσ2t + ασ*nsδt, srnsδt * σσ)
  push!(lσ2, lσ2t1)

  # draw new trait
  xt1 = rnorm(xt, srnsδt * exp(0.25*(lσ2t + lσ2t1)))
  push!(xv, xt1)

  # draw speciation rates
  lλt1 = rnorm(lλt + αλ*nsδt + βλ*(xt1 - xt), srnsδt*σλ)
  push!(lλv, lλt1)

  λm = exp(0.5*(lλt + lλt1))

  if divev(λm, nsδt)
    nn += 1
    td1, nn, lr =
      _sim_tb_it(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
        δt, srδt, lr, lU, iρi, nn, nlim)
    td2, nn, lr =
      _sim_tb_it(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
        δt, srδt, lr, lU, iρi, nn, nlim)

    return iTxb(td1, td2, bt, δt, nsδt, false, lλv, xv, lσ2), nn, lr
  end

  lλt  = lλt1
  xt   = xt1
  lσ2t = lσ2t1

  if lU < lr && nn < nlim

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t

        srt = sqrt(t)

        # draw trait rate
        lσ2t1 = rnorm(lσ2t + ασ*t, srt * σσ)
        push!(lσ2, lσ2t1)

        # draw new trait
        xt1 = rnorm(xt, srt * exp(0.25*(lσ2t + lσ2t1)))
        push!(xv, xt1)

        # draw speciation rates
        lλt1 = rnorm(lλt + αλ*t + βλ*(xt1 - xt), srt*σλ)
        push!(lλv, lλt1)

        λm = exp(0.5*(lλt + lλt1))

        if divev(λm, t)
          nn += 1
          lr  += 2.0*log(iρi)
          return iTxb(iTxb(0.0, δt, 0.0, false, 
                           [lλt1, lλt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                      iTxb(0.0, δt, 0.0, false, 
                           [lλt1, lλt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                      bt, δt, t, false, lλv, xv, lσ2), nn, lr
        else
          lr += log(iρi)
          return iTxb(bt, δt, t, false, lλv, xv, lσ2), nn, lr
        end
      end

      t  -= δt
      bt += δt

      # draw trait rate
      lσ2t1 = rnorm(lσ2t + ασ*δt, srδt * σσ)
      push!(lσ2, lσ2t1)

      # draw new trait
      xt1 = rnorm(xt, srδt * exp(0.25*(lσ2t + lσ2t1)))
      push!(xv, xt1)

      # draw speciation rates
      lλt1 = rnorm(lλt + αλ*δt + βλ*(xt1 - xt), srδt*σλ)
      push!(lλv, lλt1)

      λm = exp(0.5*(lλt + lλt1))

      if divev(λm, δt)
        nn += 1
        td1, nn, lr =
          _sim_tb_it(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
            δt, srδt, lr, lU, iρi, nn, nlim)
        td2, nn, lr =
          _sim_tb_it(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
            δt, srδt, lr, lU, iρi, nn, nlim)

        return iTxb(td1, td2, bt, δt, δt, false, lλv, xv, lσ2), nn, lr
      end

      lλt  = lλt1
      xt   = xt1
      lσ2t = lσ2t1
    end
  end

  return iTxb(), nn, NaN
end




"""
    _sim_tb_it(t   ::Float64,
               xt  ::Float64,
               lσ2t::Float64,
               ασ  ::Float64,
               σσ  ::Float64,
               lλt ::Float64,
               αλ  ::Float64,
               βλ  ::Float64,
               σλ  ::Float64,
               δt  ::Float64,
               srδt::Float64,
               lr  ::Float64,
               lU  ::Float64,
               iρi ::Float64,
               nn  ::Int64,
               nlim::Int64)

Simulate `iTxb` according to a pure-birth geometric Brownian motion for
terminal branches.
"""
function _sim_tb_it(t   ::Float64,
                    xt  ::Float64,
                    lσ2t::Float64,
                    ασ  ::Float64,
                    σσ  ::Float64,
                    lλt ::Float64,
                    αλ  ::Float64,
                    βλ  ::Float64,
                    σλ  ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    iρi ::Float64,
                    nn  ::Int64,
                    nlim::Int64)

  if lU < lr && nn < nlim

    lλv = Float64[lλt]
    xv  = Float64[xt]
    lσ2 = Float64[lσ2t]
    bt  = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t

        srt = sqrt(t)

        # draw trait rate
        lσ2t1 = rnorm(lσ2t + ασ*t, srt * σσ)
        push!(lσ2, lσ2t1)

        # draw new trait
        xt1 = rnorm(xt, srt * exp(0.25*(lσ2t + lσ2t1)))
        push!(xv, xt1)

        # draw speciation rates
        lλt1 = rnorm(lλt + αλ*t + βλ*(xt1 - xt), srt*σλ)
        push!(lλv, lλt1)

        λm = exp(0.5*(lλt + lλt1))

        if divev(λm, t)
          nn += 1
          lr += 2.0*log(iρi)
          return iTxb(iTxb(0.0, δt, 0.0, false, 
                      [lλt1, lλt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                      iTxb(0.0, δt, 0.0, false, 
                      [lλt1, lλt1], [xt1, xt1], [lσ2t1, lσ2t1]),
                      bt, δt, t, false, lλv, xv, lσ2), nn, lr
        end

        lr += log(iρi)
        return iTxb(bt, δt, t, false, lλv, xv, lσ2), nn, lr
      end

      t  -= δt
      bt += δt

      # draw trait rate
      lσ2t1 = rnorm(lσ2t + ασ*δt, srδt * σσ)
      push!(lσ2, lσ2t1)

      # draw new trait
      xt1 = rnorm(xt, srδt * exp(0.25*(lσ2t + lσ2t1)))
      push!(xv, xt1)

      # draw speciation rates
      lλt1 = rnorm(lλt + αλ*δt + βλ*(xt1 - xt), srδt*σλ)
      push!(lλv, lλt1)

      λm = exp(0.5*(lλt + lλt1))

      if divev(λm, δt)
        nn += 1
        td1, nn, lr =
          _sim_tb_it(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
            δt, srδt, lr, lU, iρi, nn, nlim)
        td2, nn, lr =
          _sim_tb_it(t, xt1, lσ2t1, ασ, σσ, lλt1, αλ, βλ, σλ, 
            δt, srδt, lr, lU, iρi, nn, nlim)

        return iTxb(td1, td2, bt, δt, δt, false, lλv, xv, lσ2), nn, lr
      end

      lλt  = lλt1
      xt   = xt1
      lσ2t = lσ2t1
    end
  end

  return iTxb(), nn, NaN
end



