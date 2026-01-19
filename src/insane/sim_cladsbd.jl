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



"""
    sim_cladsbd(n       ::Int64;
                λ0      ::Float64 = 1.0,
                α       ::Float64 = 0.0,
                σλ      ::Float64 = 0.1,
                μ       ::Float64 = 0.5,
                init    ::Symbol  = :stem,
                nstar   ::Int64   = 2*n,
                p       ::Float64 = 5.0,
                warnings::Bool    = true,
                maxt    ::Float64 = 1e7)

Simulate `cTbd` according to a birth-death clads.
"""
function sim_cladsbd(n       ::Int64;
                     λ0      ::Float64 = 1.0,
                     μ0      ::Float64 = 0.5,
                     α       ::Float64 = 0.0,
                     σλ      ::Float64 = 0.1,
                     σμ      ::Float64 = 0.1,
                     init    ::Symbol  = :stem,
                     nstar   ::Int64   = 2*n,
                     p       ::Float64 = 5.0,
                     warnings::Bool    = true,
                     maxt    ::Float64 = 1e7)

  # simulate in non-recursive manner
  e0, e1, el, λs, μs, ea, ee, na, simt =
    _sedges_cladsbd(nstar, log(λ0), log(μ0), α, σλ, σμ, init, maxt)

  if simt >= maxt
    warnings && @warn "simulation surpassed maximum time"
    return cTbd()
  end

  # transform to iTree
  t = cTbd(e0, e1, el, λs, μs, ea, ee, e1[1], 1)

  # sample a time when species(t) == `n`
  ntt = ltt(t)
  tn  = times_n(n, ntt)
  c   = usample(tn, p)

  if iszero(c)
    warnings && @warn "tree not sampled, try increasing `p`"
    return cTbd()
  else
    # cut the tree
    t = cutbottom(t, simt - c)
    return t
  end
end




"""
    _sedges_cladsbd(n   ::Int64,
                    λ0  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    μ   ::Float64,
                    init::Symbol,
                    maxt::Float64)

Simulate `cladsb` just until hitting `n` alive species. Note that this is
a biased sample for a tree conditional on `n` species.
"""
function _sedges_cladsbd(n   ::Int64,
                         λ0  ::Float64,
                         μ0   ::Float64,
                         α    ::Float64,
                         σλ   ::Float64,
                         σμ   ::Float64,
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
    # starting speciation and extinction rate
    λs = Float64[λ0]
    μs = Float64[μ0]

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
    # starting speciation and extinction rates
    λs = [λ0, rnorm(λ0 + α, σλ), rnorm(λ0 + α, σλ)]
    μs = [μ0, rnorm(μ0, σμ),     rnorm(μ0, σμ)]

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
      μa  = μs[ea]
      eλa = exp.(λa)
      eμa = exp.(μa)
      sλμ = eλa .+ eμa

      # combined event rate
      Λ = sum(sλμ)

      # time to next speciation event
      tw = cb_wait(Λ)

      # keep track of time
      simt += tw

      # time guard
      if simt > maxt
        return e0, e1, el, λs, μs, ea, ee, na, simt
      end

      # assign to an alive lineage according to their speciation rate
      si = sample(sλμ)
      λi = λa[si]
      μi = μa[si]
      wl = ea[si]

      # one time step for all edges alive `ea`
      @simd for i in ea
        el[i] += tw
      end

      # if speciation
      if λorμ(eλa[si], eμa[si])

        ## add new edges
        # start node
        push!(e0, e1[wl], e1[wl])

        # end nodes
        push!(e1, ne + 1, ne + 2)

        # push to edge length
        push!(el, 0.0, 0.0)

        # push speciation vector
        push!(λs, rnorm(λi + α, σλ), rnorm(λi + α, σλ))

        # push extinction vector
        push!(μs, rnorm(μi, σμ), rnorm(μi, σμ))

        # update living edges
        deleteat!(ea, si)
        push!(ea, ne, ne + 1)

        # update `na` and `ne`
        ne += 2
        na += 1

        # if reached `n` species
        if n === na

          Λ = sum(exp, λs[ea]) + sum(exp, μs[ea])

          # time before to next speciation event
          tw = rand()*cb_wait(Λ)

          # one time step for all edges alive `ea`
          @simd for i in ea
            el[i] += tw
          end

          return e0, e1, el, λs, μs, ea, ee, na, simt
        end

      # if extinction
      else
        # if tree goes extinct
        if isone(na)
          # extinct edges
          push!(ee, wl)
          # delete from alive lineages
          deleteat!(ea, si)

          return e0, e1, el, λs, μs, ea, ee, 0, simt
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
      return cTbd(t, false, false, λt, μt), nn
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
                   μt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64,
                   λfs ::Vector{Float64},
                   μfs ::Vector{Float64},)

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
                        λfs ::Vector{Float64},
                        μfs ::Vector{Float64})

  if nn < nlim

    λi = exp(λt)
    μi = exp(μt)
    tw = cbd_wait(λi, μi)

    if tw > t
      na += 1
      push!(λfs, λt)
      push!(μfs, μt)
      return cTbd(t, false, false, λt, μt), na, nn
    end

    if λorμ(λi, μi)
      nn += 1
      d1, na, nn = 
        _sim_cladsbd_i(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
          α, σλ, σμ, na, nn, nlim, λfs, μfs)
      d2, na, nn = 
        _sim_cladsbd_i(t - tw, rnorm(λt + α, σλ), rnorm(μt, σμ), 
          α, σλ, σμ, na, nn, nlim, λfs, μfs)

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



