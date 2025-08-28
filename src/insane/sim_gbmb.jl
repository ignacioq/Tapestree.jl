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



"""
    sim_gbmb(n       ::Int64;
              λ0      ::Float64 = 1.0,
              α       ::Float64 = 0.0,
              σλ      ::Float64 = 0.1,
              δt      ::Float64 = 1e-3,
              init    ::Symbol  = :stem,
              nstar   ::Int64   = n + 2,
              p       ::Float64 = 5.0,
              warnings::Bool    = true,
              maxt    ::Float64 = δt*1e7)

Simulate `iTb` according to a pure-birth geometric Brownian motion.
"""
function sim_gbmb(n       ::Int64;
                   λ0      ::Float64 = 1.0,
                   α       ::Float64 = 0.0,
                   σλ      ::Float64 = 0.1,
                   δt      ::Float64 = 1e-3,
                   init    ::Symbol  = :stem,
                   nstar   ::Int64   = n + 2,
                   p       ::Float64 = 5.0,
                   warnings::Bool    = true,
                   maxt    ::Float64 = δt*1e7)

  # simulate in non-recursive manner
  e0, e1, el, λs, ea, na, simt =
    _sedges_gbmb(nstar, log(λ0), α, σλ, δt, sqrt(δt), init, maxt)

  if simt >= maxt
    warnings && @warn "simulation surpassed maximum time"
    return iTb()
  end

  # transform to iTree
  t = iTb(e0, e1, el, λs, ea, e1[1], 1, δt)

  # sample a time when species(t) == `n`
  nt = ltt(t)
  tn = times_n(n, nt)
  c  = usample(tn, p)

  if iszero(c)
    warnings && @warn "tree not sampled, try increasing `p`"
    return iTb()
  else
    # cut the tree
    t = cutbottom(t, simt - c)
    return t
  end
end





"""
    _sedges_gbmb(n    ::Int64,
                  λ0   ::Float64,
                  α    ::Float64,
                  σλ   ::Float64,
                  δt   ::Float64,
                  srδt ::Float64,
                  init::Symbol)

Simulate `gbmb` just until hitting `n` alive species. Note that this is
a biased sample for a tree conditional on `n` species.
"""
function _sedges_gbmb(n    ::Int64,
                       λ0   ::Float64,
                       α    ::Float64,
                       σλ   ::Float64,
                       δt   ::Float64,
                       srδt ::Float64,
                       init::Symbol,
                       maxt ::Float64)

  # edges
  e0 = Int64[]
  e1 = Int64[]

  if init == :stem
    # edges alive
    ea = [1]
    # first edge
    push!(e0, 1)
    push!(e1, 2)
    # max index
    mxi0 = n*2
    # edge lengths
    el = [0.0]
    # lambda vector for each edge
    λs = [Float64[]]
    # starting speciation rate
    push!(λs[1], λ0)
    # lastindex for each edge
    li = [1]

    na = 1 # current number of alive species
    ne = 2 # current maximum node number

  elseif init == :crown
    # edges alive
    ea = [2, 3]
    # first edges
    push!(e0, 1, 2, 2)
    push!(e1, 2, 3, 4)
    # max index
    mxi0 = n*2
    # edge lengths
    el = [0.0, 0.0, 0.0]
    # lambda vector for each edge
    λs = [Float64[], Float64[], Float64[]]
    # starting speciation rate
    push!(λs[1], λ0, λ0)
    push!(λs[2], λ0)
    push!(λs[3], λ0)
    # lastindex for each edge
    li = [2, 1, 1]

    na = 2 # current number of alive species
    ne = 4 # current maximum node number

  else
    @error "$init does not match stem or crown"
  end

  ieaa = Int64[] # indexes of ea to add
  iead = Int64[] # indexes of ea to delete

  # simulation time
  simt = 0.0

  @inbounds begin

    # start simulation
    while true

      # keep track of time
      simt += δt

      # time guard
      if simt > maxt
        return e0, e1, el, λs, ea, na, simt
      end

      # one time step for all edges alive `ea`
      for (i,v) in enumerate(ea)

        λsi = λs[v]
        lii = li[v]
        λt  = λsi[lii]

        # update edge length
        el[v] += δt
        li[v] += 1

        # sample new speciation rates
        λt1 = rnorm(λt + α*δt, srδt*σλ)
        push!(λsi, λt1)
        λm = exp(0.5*(λt + λt1))

        # if speciation event
        if divev(λm, δt)

          # if reached `n` species
          if n === na

            # update λs and δt for other lineages
            for vi in ea[i+1:end]
              el[vi] += δt
              λsi = λs[vi]
              lvi = li[vi]
              λt  = λsi[lvi]

              push!(λsi, rnorm(λt + α*δt, srδt*σλ))
            end

            # to add
            if !isempty(ieaa)
              append!(ea, ieaa)
              empty!(ieaa)
            end

           # to delete
            if !isempty(iead)
              deleteat!(ea, iead)
              empty!(iead)
            end

            return e0, e1, el, λs, ea, na, simt
          end

          ### add new edges
          # start node
          push!(e0, e1[v], e1[v])

          # end nodes
          push!(e1, ne + 1, ne + 2)

          # push to edge length
          push!(el, 0.0, 0.0)

          # push speciation vector
          push!(λs, [λt1], [λt1])

          # push length of vector
          push!(li, 1, 1)

          # to update living edges
          push!(iead, i)
          push!(ieaa, ne, ne + 1)

          # update `na` and `ne`
          ne += 2
          na += 1
        end
      end

      # to add
      if !isempty(ieaa)
        append!(ea, ieaa)
        empty!(ieaa)
      end

      # to delete
      if !isempty(iead)
        deleteat!(ea, iead)
        empty!(iead)
      end
    end
  end
end






#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Sample conditional on time
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    sim_gbmb(t   ::Float64;
              λ0  ::Float64 = 1.0,
              α   ::Float64 = 0.0,
              σλ  ::Float64 = 0.1,
              δt  ::Float64 = 1e-3,
              nlim::Int64   = 10_000,
              init::Symbol  = :crown)

Simulate `iTb` according to a pure-birth geometric Brownian motion
conditional in stopping at time `t`.
"""
function sim_gbmb(t   ::Float64;
                   λ0  ::Float64 = 1.0,
                   α   ::Float64 = 0.0,
                   σλ  ::Float64 = 0.1,
                   δt  ::Float64 = 1e-3,
                   nlim::Int64   = 10_000,
                   init::Symbol  = :crown)

  if init === :crown
    lλ0 = log(λ0)
    d1, nn = _sim_gbmb(t, lλ0, α, σλ, δt, sqrt(δt), 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    d2, nn = _sim_gbmb(t, lλ0, α, σλ, δt, sqrt(δt), 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    tree = iTb(d1, d2, 0.0, δt, 0.0, false, [lλ0, lλ0])
   elseif init === :stem
    tree, nn = _sim_gbmb(t, log(λ0), α, σλ, δt, sqrt(δt), nn + 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_gbmb(t   ::Float64,
               λt  ::Float64,
               α   ::Float64,
               σλ  ::Float64,
               δt  ::Float64,
               srδt::Float64,
               nn ::Int64,
               nlim::Int64)

Simulate `iTb` according to a pure-birth geometric Brownian motion.
"""
function _sim_gbmb(t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    nn ::Int64,
                    nlim::Int64)

  if nn < nlim

    λv = Float64[λt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, t)
          nn += 1
          return iTb(iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                      iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                      bt, δt, t, false, λv), nn
        end

        return iTb(bt, δt, t, false, λv), nn
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nn += 1
        td1, nn = _sim_gbmb(t, λt1, α, σλ, δt, srδt, nn, nlim)
        td2, nn = _sim_gbmb(t, λt1, α, σλ, δt, srδt, nn, nlim)

        return iTb(td1, td2, bt, δt, δt, false, λv), nn
      end

      λt = λt1
    end
  end

  return iTb(), nn
end





"""
    _sim_gbmb_t(t   ::Float64,
                 λt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 Iρi ::Float64,
                 na  ::Int64,
                 nn ::Int64,
                 nlim::Int64)

Simulate `iTb` according to a pure-birth geometric Brownian motion for
terminal branches.
"""
function _sim_gbmb_t(t   ::Float64,
                      λt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      lr  ::Float64,
                      lU  ::Float64,
                      Iρi ::Float64,
                      na  ::Int64,
                      nn ::Int64,
                      nlim::Int64)

  if isfinite(lr) && nn < nlim

    λv = Float64[λt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, t)
          nn += 1
          na += 2
          if na === 2
            nlr = lr + log(Iρi*2.0)
          else
            nlr = lr + log(Iρi * Iρi * Float64(na)/Float64(na-2))
          end
          if nlr < lr && lU >= nlr
            return iTb(), na, nn, NaN
          else
            return iTb(iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                        iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                        bt, δt, t, false, λv), na, nn, nlr
          end
        else
          na += 1
          nlr = lr
          if na > 1
            nlr += log(Iρi * Float64(na)/Float64(na-1))
          end
          if nlr >= lr
            return iTb(bt, δt, t, false, λv), na, nn, nlr
          elseif lU < nlr
            return iTb(bt, δt, t, false, λv), na, nn, nlr
          else
            return iTb(), na, nn, NaN
          end
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nn += 1
        td1, na, nn, lr =
          _sim_gbmb_t(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, na, nn, nlim)
        td2, na, nn, lr =
          _sim_gbmb_t(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, na, nn, nlim)

        return iTb(td1, td2, bt, δt, δt, false, λv), na, nn, lr
      end

      λt = λt1
    end
  end

  return iTb(), na, nn, NaN
end




"""
    _sim_gbmb_it(nsδt::Float64,
                  t   ::Float64,
                  λt  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  lr  ::Float64,
                  lU  ::Float64,
                  Iρi ::Float64,
                  na  ::Int64,
                  nn  ::Int64,
                  nlim::Int64)
Simulate `iTb` according to a pure-birth geometric Brownian motion,
starting with a non-standard δt with a limit in the number of species.
"""
function _sim_gbmb_it(nsδt::Float64,
                       t   ::Float64,
                       λt  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       Iρi ::Float64,
                       nn  ::Int64,
                       nlim::Int64)

  λv = Float64[λt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt + accerr
    t   = isapprox(t, 0.0) ? 0.0 : isapprox(t, nsδt) ? nsδt : t
    bt += t
    λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
    λm  = exp(0.5*(λt + λt1))
    push!(λv, λt1)

    if divev(λm, t)
      nn += 1
      lr += 2.0*log(Iρi)
      return iTb(iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                  iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                  bt, δt, t, false, λv), nn, lr
    else
      lr += log(Iρi)
      return iTb(bt, δt, t, false, λv), nn, lr
    end
  end

  t  -= nsδt
  bt += nsδt

  λt1 = rnorm(λt + α*nsδt, sqrt(nsδt)*σλ)
  λm  = exp(0.5*(λt + λt1))
  push!(λv, λt1)

  if divev(λm, nsδt)
    nn += 1
    td1, nn, lr =
      _sim_gbmb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)
    td2, nn, lr =
      _sim_gbmb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)

    return iTb(td1, td2, bt, δt, nsδt, false, λv), nn, lr
  end

  λt = λt1

  if lU < lr && nn < nlim

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, t)
          nn += 1
          lr  += 2.0*log(Iρi)
          return iTb(iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                      iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                      bt, δt, t, false, λv), nn, lr
        else
          lr += log(Iρi)
          return iTb(bt, δt, t, false, λv), nn, lr
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nn += 1
        td1, nn, lr =
          _sim_gbmb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)
        td2, nn, lr =
          _sim_gbmb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)

        return iTb(td1, td2, bt, δt, δt, false, λv), nn, lr
      end

      λt = λt1
    end
  end

  return iTb(), nn, NaN
end




"""
    _sim_gbmb_it(t   ::Float64,
                 λt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 Iρi ::Float64,
                 nn ::Int64,
                 nlim::Int64)

Simulate `iTb` according to a pure-birth geometric Brownian motion for
terminal branches.
"""
function _sim_gbmb_it(t   ::Float64,
                       λt  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       Iρi ::Float64,
                       nn ::Int64,
                       nlim::Int64)

  if lU < lr && nn < nlim

    λv = Float64[λt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, t)
          nn += 1
          lr  += 2.0*log(Iρi)
          return iTb(iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                      iTb(0.0, δt, 0.0, false, [λt1, λt1]),
                      bt, δt, t, false, λv), nn, lr
        end

        lr += log(Iρi)
        return iTb(bt, δt, t, false, λv), nn, lr
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, δt)
        nn += 1
        td1, nn, lr =
          _sim_gbmb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)
        td2, nn, lr =
          _sim_gbmb_it(t, λt1, α, σλ, δt, srδt, lr, lU, Iρi, nn, nlim)

        return iTb(td1, td2, bt, δt, δt, false, λv), nn, lr
      end

      λt = λt1
    end
  end

  return iTb(), nn, NaN
end




"""
    divev(λ::Float64, δt::Float64)

Return true if diversification event.
"""
divev(λ::Float64, δt::Float64) = @fastmath rand() < λ*δt


