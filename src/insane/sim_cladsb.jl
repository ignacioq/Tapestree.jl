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
    sim_cladsb(n       ::Int64;
              λ0      ::Float64 = 1.0,
              α       ::Float64 = 0.0,
              σλ      ::Float64 = 0.1,
              init    ::Symbol  = :stem,
              nstar   ::Int64   = n + 2,
              p       ::Float64 = 5.0,
              warnings::Bool    = true,
              maxt    ::Float64 = δt*1e7)

Simulate `cTb` according to a pure-birth clads.
"""
function sim_cladsb(n       ::Int64;
                    λ0      ::Float64 = 1.0,
                    α       ::Float64 = 0.0,
                    σλ      ::Float64 = 0.1,
                    init    ::Symbol  = :stem,
                    nstar   ::Int64   = n + 1,
                    p       ::Float64 = 5.0,
                    warnings::Bool    = true,
                    maxt    ::Float64 = 1e7)


  # simulate in non-recursive manner
  e0, e1, el, λs, ea, na, simt =
    _sedges_cladsb(nstar, log(λ0), α, σλ, init, maxt)

  if simt >= maxt
    warnings && @warn "simulation surpassed maximum time"
    return cTb()
  end

  # transform to iTree
  t = cTb(e0, e1, el, λs, ea, e1[1], 1)

  # sample a time when species(t) == `n`
  ntt = ltt(t)
  tn  = times_n(n, ntt)
  c   = usample(tn, p)

  if iszero(c)
    warnings && @warn "tree not sampled, try increasing `p`"
    return cTb()
  else
    # cut the tree
    t = cutbottom(t, simt - c)
    return t
  end
end




"""
    _sedges_cladsb(n   ::Int64,
                   λ0  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   init::Symbol,
                   maxt::Float64)

Simulate `cladsb` just until hitting `n` alive species.
"""
function _sedges_cladsb(n   ::Int64,
                        λ0  ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        init::Symbol,
                        maxt::Float64)

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
    # starting speciation rate
    λs = Float64[λ0]

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
    el = zeros(Float64, 3)
    # starting speciation rates
    λs = [λ0, rnorm(λ0 + α, σλ), rnorm(λ0 + α, σλ)]

    na = 2 # current number of alive species
    ne = 4 # current maximum node number

  else
    @error "$init does not match stem or crown"
  end

  # simulation time
  simt = 0.0

  @inbounds begin

    # start simulation
    while true

      λa  = λs[ea]
      eλa = exp.(λa)

      # combined event rate
      Λ = sum(eλa)

      # time to next speciation event
      tw = cb_wait(Λ)

      # keep track of time
      simt += tw

      # time guard
      if simt > maxt
        return e0, e1, el, λs, ea, na, simt
      end

      # assign to an alive lineage according to their speciation rate
      si = sample(eλa)
      λi = λa[si]
      wl = ea[si]

      # one time step for all edges alive `ea`
      @simd for i in ea
        el[i] += tw
      end

      ### add new edges
      # start node
      push!(e0, e1[wl], e1[wl])

      # end nodes
      push!(e1, ne + 1, ne + 2)

      # push to edge length
      push!(el, 0.0, 0.0)

      # push speciation vector
      push!(λs, rnorm(λi + α, σλ), rnorm(λi + α, σλ))

      # update living edges
      deleteat!(ea, si)
      push!(ea, ne, ne + 1)

      # update `na` and `ne`
      ne += 2
      na += 1

      # if reached `n` species
      if n === na

        Λ = sum(exp, λs[ea])

        # time before to next speciation event
        tw = rand()*cb_wait(Λ)

        # one time step for all edges alive `ea`
        @simd for i in ea
          el[i] += tw
        end

        return e0, e1, el, λs, ea, na, simt
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




