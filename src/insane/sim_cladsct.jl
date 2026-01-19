#=

cladogenetic constant turnover simulation

Ignacio Quintero Mächler

t(-_-t)

Created 24 07 2025
=#




#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Sample conditional on number of species
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    sim_cladsct(n       ::Int64;
                λ0      ::Float64 = 1.0,
                α       ::Float64 = 0.0,
                σλ      ::Float64 = 0.1,
                ϵ       ::Float64 = 0.5,
                init    ::Symbol  = :stem,
                nstar   ::Int64   = 2*n,
                p       ::Float64 = 5.0,
                warnings::Bool    = true,
                maxt    ::Float64 = 1e7)

Simulate `cTct` according to a constant turnover clads.
"""
function sim_cladsct(n       ::Int64;
                     λ0      ::Float64 = 1.0,
                     α       ::Float64 = 0.0,
                     σλ      ::Float64 = 0.1,
                     ϵ       ::Float64 = 0.5,
                     init    ::Symbol  = :stem,
                     nstar   ::Int64   = 2*n,
                     p       ::Float64 = 5.0,
                     warnings::Bool    = true,
                     maxt    ::Float64 = 1e7)

  # simulate in non-recursive manner
  e0, e1, el, λs, ea, ee, na, simt =
    _sedges_cladsct(nstar, log(λ0), α, σλ, ϵ, init, maxt)

  if simt >= maxt
    warnings && @warn "simulation surpassed maximum time"
    return cTct()
  end

  # transform to iTree
  t = cTct(e0, e1, el, λs, ea, ee, e1[1], 1)

  # sample a time when species(t) == `n`
  ntt = ltt(t)
  tn  = times_n(n, ntt)
  c   = usample(tn, p)

  if iszero(c)
    warnings && @warn "tree not sampled, try increasing `p`"
    return cTct()
  else
    # cut the tree
    t = cutbottom(t, simt - c)
    return t
  end
end




"""
    _sedges_cladsct(n   ::Int64,
                    λ0  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    init::Symbol,
                    maxt::Float64)

Simulate `cladsb` just until hitting `n` alive species. Note that this is
a biased sample for a tree conditional on `n` species.
"""
function _sedges_cladsct(n   ::Int64,
                         λ0  ::Float64,
                         α   ::Float64,
                         σλ  ::Float64,
                         ϵ   ::Float64,
                         init::Symbol,
                         maxt::Float64)

  # edges
  e0 = Int64[]
  e1 = Int64[]
  # edges extinct
  ee = Int64[]

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
    @error string(init, " does not match stem or crown")
  end

  # simulation time
  simt = 0.0

  @inbounds begin

    # start simulation
    while true

      λa  = λs[ea]
      eλa = exp.(λa)

      # combined event rate
      Λ = sum(eλa) * (1.0 + ϵ)

      # time to next speciation event
      tw = cb_wait(Λ)

      # keep track of time
      simt += tw

      # time guard
      if simt > maxt
        return e0, e1, el, λs, ea, ee, na, simt
      end

      # assign to an alive lineage according to their speciation rate
      si = sample(eλa)
      λi = λa[si]
      wl = ea[si]

      # one time step for all edges alive `ea`
      @simd for i in ea
        el[i] += tw
      end

      # if speciation
      if λorμ(ϵ)

        ## add new edges
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

          Λ = sum(exp, λs[ea]) * (1.0 + ϵ)

          # time before to next speciation event
          tw = rand()*cb_wait(Λ)

          # one time step for all edges alive `ea`
          @simd for i in ea
            el[i] += tw
          end

          return e0, e1, el, λs, ea, ee, na, simt
        end

      # if extinction
      else
        # if tree goes extinct
        if isone(na)
          # extinct edges
          push!(ee, wl)
          # delete from alive lineages
          deleteat!(ea, si)

          return e0, e1, el, λs, ea, ee, 0, simt
        end

        # extinct edges
        push!(ee, wl)
        # delete from alive lineages
        deleteat!(ea, si)
        # update number of alive species
        na -= 1
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
    sim_cladsct(t   ::Float64;
                λ0  ::Float64 = 1.0,
                α   ::Float64 = 0.0,
                σλ  ::Float64 = 0.1,
                ϵ   ::Float64 = 0.5,
                nlim::Int64   = 10_000,
                init::Symbol  = :crown)

Simulate `cTct` according to a constant turnover clads
conditional in stopping at time `t`.
"""
function sim_cladsct(t   ::Float64;
                     λ0  ::Float64 = 1.0,
                     α   ::Float64 = 0.0,
                     σλ  ::Float64 = 0.1,
                     ϵ   ::Float64 = 0.5,
                     nlim::Int64   = 10_000,
                     init::Symbol  = :crown)

  if init === :crown
    lλ0 = log(λ0)
    d1, nn = _sim_cladsct(t, lλ0, α, σλ, ϵ, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    d2, nn = _sim_cladsct(t, lλ0, α, σλ, ϵ, nn + 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    tree = cTct(d1, d2, 0.0, false, lλ0)
   elseif init === :stem
    tree, nn = _sim_cladsct(t, log(λ0), α, σλ, ϵ, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_cladsct(t   ::Float64,
                 λt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 ϵ   ::Float64,
                 nn  ::Int64,
                 nlim::Int64)


Simulate `cTct` according to a constant turnover clads.
"""
function _sim_cladsct(t   ::Float64,
                      λt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      ϵ   ::Float64,
                      nn  ::Int64,
                      nlim::Int64)

  if nn < nlim

    λi = exp(λt)
    tw = ct_wait(λi, ϵ)

    if tw > t
      return cTct(t, false, false, λt), nn
    end

    if λorμ(ϵ)
      nn += 1
      d1, nn = _sim_cladsct(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ, nn, nlim)
      d2, nn = _sim_cladsct(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ, nn, nlim)
 
      return cTct(d1, d2, tw, false, false, λt), nn
    else
      return cTct(tw, true, false, λt), nn
    end
  end

  return cTct(), nn
end




"""
    _sim_cladsct_t(t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   ϵ   ::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64)

Simulate `cTct` according to a constant turnover clads for
terminal branches.
"""
function _sim_cladsct_t(t   ::Float64,
                        λt  ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        ϵ   ::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64,
                        nlim::Int64)

  if isfinite(lr) && nn < nlim

    λi = exp(λt)
    tw = ct_wait(λi, ϵ)

    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(iρi * Float64(na)/Float64(na-1))
      end
      if nlr >= lr || lU < nlr
        return cTct(t, false, false, λt), na, nn, nlr
      else
        return cTct(), na, nn, NaN
      end
    end

    if λorμ(ϵ)
      nn += 1
      td1, na, nn, lr =
        _sim_cladsct_t(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ,
          lr, lU, iρi, na, nn, nlim)
      td2, na, nn, lr =
        _sim_cladsct_t(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ,
          lr, lU, iρi, na, nn, nlim)

      return cTct(td1, td2, tw, false, false, λt), na, nn, lr
    else
      return cTct(tw, true, false, λt), na, nn, lr
    end
  end

  return cTct(), na, nn, NaN
end




"""
    _sim_cladsct_i(t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   ϵ   ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64,
                   xfs::Vector{Float64})

Simulate `cTct` according to a constant turnover clads.
"""
function _sim_cladsct_i(t   ::Float64,
                        λt  ::Float64,
                        α   ::Float64,
                        σλ  ::Float64,
                        ϵ   ::Float64,
                        na  ::Int64,
                        nn  ::Int64,
                        nlim::Int64,
                        xfs::Vector{Float64})

  if nn < nlim

    λi = exp(λt)
    tw = ct_wait(λi, ϵ)

    if tw > t
      na += 1
      push!(xfs, λt)
      return cTct(t, false, false, λt), na, nn
    end

    if λorμ(ϵ)
      nn += 1
      d1, na, nn = 
        _sim_cladsct_i(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ, na, nn, nlim, xfs)
      d2, na, nn = 
        _sim_cladsct_i(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ, na, nn, nlim, xfs)

      return cTct(d1, d2, tw, false, false, λt), na, nn
    else
      return cTct(tw, true, false, λt), na, nn
    end
  end

  return cTct(), na, nn
end



"""
    _sim_cladsct_it(t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    lr  ::Float64,
                    lU  ::Float64,
                    iρi ::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

Simulate `cTct` according to a constant turnover clads for
internal terminal branches.
"""
function _sim_cladsct_it(t   ::Float64,
                         λt  ::Float64,
                         α   ::Float64,
                         σλ  ::Float64,
                         ϵ   ::Float64,
                         lr  ::Float64,
                         lU  ::Float64,
                         iρi ::Float64,
                         na  ::Int64,
                         nn  ::Int64,
                         nlim::Int64)
  if lU < lr && nn < nlim

    λi = exp(λt)
    tw = ct_wait(λi, ϵ)

    if tw > t
      na += 1
      lr += log(iρi)
      return cTct(t, false, false, λt), na, nn, lr
    end

    if λorμ(ϵ)
      nn += 1
      td1, na, nn, lr = _sim_cladsct_it(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ,
        lr, lU, iρi, na, nn, nlim)
      td2, na, nn, lr = _sim_cladsct_it(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ,
        lr, lU, iρi, na, nn, nlim)

      return cTct(td1, td2, tw, false, false, λt), na, nn, lr
    else
      return cTct(tw, true, false, λt), na, nn, lr
    end
  end

  return cTct(), na, nn, NaN
end




"""
    _sim_cladsct_surv(t   ::Float64,
                      λt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      ϵ   ::Float64,
                      surv::Bool,
                      nn  ::Int64)

Simulate `cTct` according to a constant turnover clads.
"""
function _sim_cladsct_surv(t   ::Float64,
                           λt  ::Float64,
                           α   ::Float64,
                           σλ  ::Float64,
                           ϵ   ::Float64,
                           surv::Bool,
                           nn  ::Int64)

  if !surv && nn < 500

    λi = exp(λt)
    tw = ct_wait(λi, ϵ)

    if tw > t
      return true, nn
    end

    if λorμ(ϵ)
      nn += 1
      surv, nn = 
        _sim_cladsct_surv(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ, surv, nn)
      surv, nn = 
        _sim_cladsct_surv(t - tw, rnorm(λt + α, σλ), α, σλ, ϵ, surv, nn)
 
      return surv, nn
    else
      return surv, nn
    end
  end

  return true, nn
end





"""
    ct_wait(λ::Float64, ϵ::Float64)

Sample a per-lineage waiting time for birth-death with constant extinction with 
speciation rate `λ` and turnover rate `ϵ`.
"""
ct_wait(λ::Float64, ϵ::Float64) = rexp(λ*(1.0 + ϵ))



