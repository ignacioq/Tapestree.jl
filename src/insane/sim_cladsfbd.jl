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
    sim_cladsfbd(n       ::Int64;
                 λ0      ::Float64         = 1.0,
                 μ0      ::Float64         = 0.5,
                 αλ      ::Float64         = 0.0,
                 αμ      ::Float64         = 0.0,
                 σλ      ::Float64         = 0.1,
                 σμ      ::Float64         = 0.1,
                 ψ       ::Vector{Float64} = [0.1],
                 ψts     ::Vector{Float64} = Float64[],
                 init    ::Symbol          = :stem,
                 nstar   ::Int64           = 2*n,
                 p       ::Float64         = 5.0,
                 warnings::Bool            = true,
                 maxt    ::Float64         = 1e7)

Simulate `cTfbd` according to a birth-death clads.
"""
function sim_cladsfbd(n       ::Int64;
                      λ0      ::Float64         = 1.0,
                      μ0      ::Float64         = 0.5,
                      αλ      ::Float64         = 0.0,
                      αμ      ::Float64         = 0.0,
                      σλ      ::Float64         = 0.1,
                      σμ      ::Float64         = 0.1,
                      ψ       ::Vector{Float64} = [0.1],
                      ψ_epoch ::Vector{Float64} = Float64[],
                      init    ::Symbol          = :stem,
                      nstar   ::Int64           = 2*n,
                      p       ::Float64         = 5.0,
                      warnings::Bool            = true,
                      maxt    ::Float64         = 1e7)

  # simulate in non-recursive manner
  e0, e1, el, λs, μs, ea, ee, ef, na, simt =
    _sedges_cladsfbd(nstar, log(λ0), log(μ0), αλ, αμ, σλ, σμ, ψ, ψ_epoch, 
      init, maxt)

  if simt >= maxt
    warnings && @warn "simulation surpassed maximum time"
    return cTfbd()
  end

  # transform to iTree
  t = cTfbd(e0, e1, el, λs, μs, ea, ee, ef, e1[1], 1)

  # sample a time when species(t) == `n`
  ntt = ltt(t)
  tn  = times_n(n, ntt)
  c   = usample(tn, p)

  if iszero(c)
    warnings && @warn "tree not sampled, try increasing `p`"
    return cTfbd()
  else
    # cut the tree
    t = cutbottom(t, simt - c)
    return t
  end
end




"""
    _sedges_cladsfbd(n   ::Int64,
                     λ0  ::Float64,
                     μ0  ::Float64,
                     αλ  ::Float64,
                     αμ  ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     ψ   ::Vector{Float64},
                     ψts ::Vector{Float64},
                     init::Symbol,
                     maxt::Float64)

Simulate `cladsfbd` just until hitting `n` alive species. Note that this is
a biased sample for a tree conditional on `n` species.
"""
function _sedges_cladsfbd(n   ::Int64,
                          λ0  ::Float64,
                          μ0  ::Float64,
                          αλ  ::Float64,
                          αμ  ::Float64,
                          σλ  ::Float64,
                          σμ  ::Float64,
                          ψ   ::Vector{Float64},
                          ψts ::Vector{Float64},
                          init::Symbol,
                          maxt::Float64)

  # edges
  e0 = Int64[]
  e1 = Int64[]
  ee = Int64[] # edges extinct
  ef = Int64[] # edges fossil

  if init == :stem
    # edges alive
    ea = [1]
    # first edge
    push!(e0, 1)
    push!(e1, 2)
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
    # edge lengths
    el = zeros(Float64, 3)
    # starting speciation and extinction rates
    λs = [λ0, rnorm(λ0 + αλ, σλ), rnorm(λ0 + αλ, σλ)]
    μs = [μ0, rnorm(μ0 + αμ, σμ), rnorm(μ0 + αμ, σμ)]

    na = 2 # current number of alive species
    ne = 4 # current maximum node number

  else
    @error string(init, " does not match stem or crown")
  end

  # simulation time
  simt = 0.0
  nep  = lastindex(ψts) + 1        # n epochs
  ix   = 1                         # epoch index ix
  ψi   = ψ[ix]                     # epoch rate
  et   = nep > 1 ? ψts[1] : Inf   # current epoch

  @inbounds begin

    # start simulation
    while true

      λa  = λs[ea]
      μa  = μs[ea]
      eλa = exp.(λa)
      eμa = exp.(μa)
      sλμ = eλa .+ eμa

      # combined event rate
      Λ = sum(sλμ) + Float64(na)*ψi

      # time to next speciation event
      tw = cb_wait(Λ)

      # if change of epoch is needed
      while simt + tw >= et
        simt = et
        ix  += 1
        ψi  = ψ[ix]
        et  = ix < nep ? ψts[ix] : Inf

        # time to next speciation event
        Λ = sum(sλμ) + Float64(na)*ψi
        tw = cb_wait(Λ)
      end

      # keep track of time
      simt += tw

      # time guard
      if simt > maxt
        return e0, e1, el, λs, μs, ea, ee, ef, na, simt
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
      if λevent(eλa[si], eμa[si], ψi)

        ## add new edges
        # start node
        push!(e0, e1[wl], e1[wl])

        # end nodes
        push!(e1, ne + 1, ne + 2)

        # push to edge length
        push!(el, 0.0, 0.0)

        # push speciation
        push!(λs, rnorm(λi + αλ, σλ), rnorm(λi + αλ, σλ))

        # push extinction
        push!(μs, rnorm(μi + αμ, σμ), rnorm(μi + αμ, σμ))

        # update living edges
        deleteat!(ea, si)
        push!(ea, ne, ne + 1)

        # update `na` and `ne`
        ne += 2
        na += 1

        # if reached `n` species
        if n === na

          Λ = sum(exp, λs[ea]) .+ sum(exp, μs[ea]) + Float64(na)*ψi

          # time before to next speciation event
          tw = rand()*cb_wait(Λ)

          # one time step for all edges alive `ea`
          @simd for i in ea
            el[i] += tw
          end

          return e0, e1, el, λs, μs, ea, ee, ef, na, simt
        end

      # if extinction
      elseif μevent(eμa[si], ψi)
        # if tree goes extinct
        if isone(na)
          # extinct edges
          push!(ee, wl)
          # delete from alive lineages
          deleteat!(ea, si)

          return e0, e1, el, λs, μs, ea, ee, ef, 0, simt
        end

        # extinct edges
        push!(ee, wl)
        # delete from alive lineages
        deleteat!(ea, si)
        # update number of alive species
        na -= 1

      # if fossilization
      else
        # start node
        push!(e0, e1[wl])

        # end nodes
        push!(e1, ne + 1)

        # push to edge length
        push!(el, 0.0)

        # push speciation and extinction
        push!(λs, λi)
        push!(μs, μi)

        # to update fossil edges
        push!(ef, wl)

        # update living edges
        deleteat!(ea, si)
        push!(ea, ne)

        # update `na` and `ne`
        ne += 1
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
    sim_cladsfbd(t   ::Float64;
                 λ0  ::Float64         = 1.0,
                 μ0  ::Float64         = 0.5,
                 αλ  ::Float64         = 0.0,
                 αμ  ::Float64         = 0.0,
                 σλ  ::Float64         = 0.1,
                 σμ  ::Float64         = 0.1,
                 ψ   ::Vector{Float64} = [0.1],
                 ψts ::Vector{Float64} = Float64[],
                 nlim::Int64           = 10_000,
                 init::Symbol          = :crown)

Simulate `cTfbd` according to a birth-death clads
conditional in stopping at time `t`.
"""
function sim_cladsfbd(t   ::Float64;
                      λ0  ::Float64         = 1.0,
                      μ0  ::Float64         = 0.5,
                      αλ  ::Float64         = 0.0,
                      αμ  ::Float64         = 0.0,
                      σλ  ::Float64         = 0.1,
                      σμ  ::Float64         = 0.1,
                      ψ   ::Vector{Float64} = [0.1],
                      ψts ::Vector{Float64} = Float64[],
                      nlim::Int64           = 10_000,
                      init::Symbol          = :crown)

  # only include epochs where the tree occurs
  tix = findfirst(x -> x < t, ψts)
  if !isnothing(tix)
    ψ   = ψ[tix:end]
    ψts = ψts[tix:end]
  end
  nep  = lastindex(ψts) + 1

  if init === :crown
    lλ0 = log(λ0)
    lμ0 = log(μ0)
    d1, nn = _sim_cladsfbd(t, lλ0, lμ0, αλ, αμ, σλ, σμ, ψ, ψts, 
      1, nep, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    d2, nn = _sim_cladsfbd(t, lλ0, lμ0, αλ, αμ, σλ, σμ, ψ, ψts, 
      1, nep, nn + 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    tree = cTfbd(d1, d2, 0.0, false, false, false, lλ0, lμ0)
  elseif init === :stem
    tree, nn = _sim_cladsfbd(t, log(λ0), log(μ0), αλ, αμ, σλ, σμ, ψ, ψts, 
      1, nep, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_cladsfbd(t   ::Float64,
                  λt  ::Float64,
                  μt  ::Float64,
                  αλ  ::Float64,
                  αμ  ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  ψ   ::Vector{Float64},
                  ψts ::Vector{Float64},
                  ix  ::Int64,
                  nep ::Int64,
                  nn  ::Int64,
                  nlim::Int64)

Simulate `cTfbd` according to a birth-death clads.
"""
function _sim_cladsfbd(t   ::Float64,
                       λt  ::Float64,
                       μt  ::Float64,
                       αλ  ::Float64,
                       αμ  ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       ψ   ::Vector{Float64},
                       ψts ::Vector{Float64},
                       ix  ::Int64,
                       nep ::Int64,
                       nn  ::Int64,
                       nlim::Int64)

  if nn < nlim

    @inbounds ψi = ψ[ix]

    λi = exp(λt)
    μi = exp(μt)
    tw = cfbd_wait(λi, μi, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        t0, nn = _sim_cladsfbd(ψti, λt, μt, αλ, αμ, σλ, σμ, ψ, ψts, 
                   ix + 1, nep, nn, nlim)
        sete!(t0, e(t0) + e0)
        return t0, nn
      end
    end

    if tw > t
      return cTfbd(t, false, false, false, λt, μt), nn
    end

    # if speciation
    if λevent(λi, μi, ψi)
      nn += 1
      d1, nn = _sim_cladsfbd(t - tw, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ),
                  αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, nn, nlim)
      d2, nn = _sim_cladsfbd(t - tw, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ),
                  αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, nn, nlim)
 
      return cTfbd(d1, d2, tw, false, false, false, λt, μt), nn
    # if extinction
    elseif μevent(μi, ψi)

      return cTfbd(tw, true, false, false, λt, μt), nn
    # if fossil sampling
    else
      d1, nn = _sim_cladsfbd(t - tw, λt, μt, αλ, αμ, σλ, σμ, ψ, ψts, 
                 ix, nep, nn, nlim)

      return cTfbd(d1, tw, false, true, false, λt, μt), nn
    end
  end

  return cTfbd(), nn
end




"""
    _sim_cladsfbd_t(t   ::Float64,
                    λt  ::Float64,
                    μt  ::Float64,
                    αλ  ::Float64,
                    αμ  ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Vector{Float64},
                    ψts ::Vector{Float64},
                    ix  ::Int64,
                    nep ::Int64,
                    lr  ::Float64,
                    lU  ::Float64,
                    iρi ::Float64,
                    na  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

Simulate `cTfbd` according to a birth-death clads for
terminal branches.
"""
function _sim_cladsfbd_t(t   ::Float64,
                         λt  ::Float64,
                         μt  ::Float64,
                         αλ  ::Float64,
                         αμ  ::Float64,
                         σλ  ::Float64,
                         σμ  ::Float64,
                         ψ   ::Vector{Float64},
                         ψts ::Vector{Float64},
                         ix  ::Int64,
                         nep ::Int64,
                         lr  ::Float64,
                         lU  ::Float64,
                         iρi ::Float64,
                         na  ::Int64,
                         nn  ::Int64,
                         nlim::Int64)

  if isfinite(lr) && nn < nlim

    @inbounds ψi = ψ[ix]

    λi = exp(λt)
    μi = exp(μt)
    tw = cfbd_wait(λi, μi, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        t0, na, nn, lr = _sim_cladsfbd_t(ψti, λt, μt, αλ, αμ, σλ, σμ, 
           ψ, ψts, ix + 1, nep, lr, lU, iρi, na, nn, nlim)
        sete!(t0, e(t0) + e0)
        return t0, na, nn, lr
      end
    end

    # if survives to present
    if tw > t
      na += 1
      nlr = lr
      if na > 1
        nlr += log(iρi * Float64(na)/Float64(na-1))
      end
      if nlr >= lr || lU < nlr
        return cTfbd(t, false, false, false, λt, μt), na, nn, nlr
      else
        return cTfbd(), na, nn, NaN
      end
    end

    # if speciation
    if λevent(λi, μi, ψi)
      nn += 1
      d1, na, nn, lr =
        _sim_cladsfbd_t(t - tw, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, lr, lU, iρi, na, nn, nlim)
      d2, na, nn, lr =
        _sim_cladsfbd_t(t - tw, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, lr, lU, iρi, na, nn, nlim)

      return cTfbd(d1, d2, tw, false, false, false, λt, μt), na, nn, lr
    # if extinction
    elseif μevent(μi, ψi)

      return cTfbd(tw, true, false, false, λt, μt), na, nn, lr
    # if fossil sampling
    else
      return cTfbd(), na, nn, NaN
    end
  end

  return cTfbd(), na, nn, NaN
end




"""
    _sim_cladsfbd_i(t   ::Float64,
                    te  ::Float64,
                    λt  ::Float64,
                    μt  ::Float64,
                    αλ  ::Float64,
                    αμ  ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Vector{Float64},
                    ψts ::Vector{Float64},
                    ix  ::Int64,
                    nep ::Int64,
                    na  ::Int64,
                    nf  ::Int64,
                    nn  ::Int64,
                    nlim::Int64,
                    λfs ::Vector{Float64},
                    μfs ::Vector{Float64})

Simulate `cTfbd` according to a birth-death clads.
"""
function _sim_cladsfbd_i(t   ::Float64,
                         te  ::Float64,
                         λt  ::Float64,
                         μt  ::Float64,
                         αλ  ::Float64,
                         αμ  ::Float64,
                         σλ  ::Float64,
                         σμ  ::Float64,
                         ψ   ::Vector{Float64},
                         ψts ::Vector{Float64},
                         ix  ::Int64,
                         nep ::Int64,
                         na  ::Int64,
                         af  ::Bool,
                         nn  ::Int64,
                         nlim::Int64,
                         λfs ::Vector{Float64},
                         μfs ::Vector{Float64})

  if !af && nn < nlim

    @inbounds ψi = ψ[ix]

    λi = exp(λt)
    μi = exp(μt)
    tw = cfbd_wait(λi, μi, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        t0, na, af, nn = _sim_cladsfbd_i(ψti, te, λt, μt, αλ, αμ, σλ, σμ, 
          ψ, ψts, ix + 1, nep, na, af, nn, nlim, λfs, μfs)
        sete!(t0, e(t0) + e0)
        return t0, na, af, nn
      end
    end

    if tw > (t - te)
      na += 1
      push!(λfs, λt)
      push!(μfs, μt)
      return cTfbd(t - te, false, false, false, λt, μt), na, af, nn
    end

    # if speciation
    if λevent(λi, μi, ψi)
      nn += 1
      d1, na, af, nn = 
        _sim_cladsfbd_i(t - tw, te, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, na, af, nn, nlim, λfs, μfs)
      d2, na, af, nn = 
        _sim_cladsfbd_i(t - tw, te, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, na, af, nn, nlim, λfs, μfs)

      return cTfbd(d1, d2, tw, false, false, false, λt, μt), na, af, nn
    # if extinction
    elseif μevent(μi, ψi)

      return cTfbd(tw, true, false, false, λt, μt), na, af, nn
    # if fossil sampling
    else
      return cTfbd(), na, true, nn
    end
  end

  return cTfbd(), na, af, nn
end




"""
    _sim_cladsfbd_i(t   ::Float64,
                    te  ::Float64,
                    λt  ::Float64,
                    μt  ::Float64,
                    αλ  ::Float64,
                    αμ  ::Float64,
                    σλ  ::Float64,
                    σμ  ::Float64,
                    ψ   ::Vector{Float64},
                    ψts ::Vector{Float64},
                    ix  ::Int64,
                    nep ::Int64,
                    na  ::Int64,
                    nf  ::Int64,
                    nn  ::Int64,
                    nlim::Int64)

Simulate `cTfbd` according to a birth-death clads.
"""
function _sim_cladsfbd_i(t   ::Float64,
                         te  ::Float64,
                         λt  ::Float64,
                         μt  ::Float64,
                         αλ  ::Float64,
                         αμ  ::Float64,
                         σλ  ::Float64,
                         σμ  ::Float64,
                         ψ   ::Vector{Float64},
                         ψts ::Vector{Float64},
                         ix  ::Int64,
                         nep ::Int64,
                         na  ::Int64,
                         af  ::Bool,
                         nn  ::Int64,
                         nlim::Int64)

  if !af && nn < nlim

    @inbounds ψi = ψ[ix]

    λi = exp(λt)
    μi = exp(μt)
    tw = cfbd_wait(λi, μi, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        t0, na, af, nn = _sim_cladsfbd_i(ψti, te, λt, μt, αλ, αμ, σλ, σμ, 
          ψ, ψts, ix + 1, nep, na, af, nn, nlim)
        sete!(t0, e(t0) + e0)
        return t0, na, af, nn
      end
    end

    if tw > (t - te)
      na += 1
      return cTfbd(t - te, false, false, false, λt, μt), na, af, nn
    end

    # if speciation
    if λevent(λi, μi, ψi)
      nn += 1
      d1, na, af, nn = 
        _sim_cladsfbd_i(t - tw, te, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, na, af, nn, nlim)
      d2, na, af, nn = 
        _sim_cladsfbd_i(t - tw, te, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, na, af, nn, nlim)

      return cTfbd(d1, d2, tw, false, false, false, λt, μt), na, af, nn
    # if extinction
    elseif μevent(μi, ψi)

      return cTfbd(tw, true, false, false, λt, μt), na, af, nn
    # if fossil sampling
    else
      return cTfbd(), na, true, nn
    end
  end

  return cTfbd(), na, af, nn
end




"""
    _sim_cladsfbd_it(t   ::Float64,
                     λt  ::Float64,
                     μt  ::Float64,
                     αλ  ::Float64, 
                     αμ  ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     ψ   ::Vector{Float64},
                     ψts ::Vector{Float64},
                     ix  ::Int64,
                     nep ::Int64,
                     lr  ::Float64,
                     lU  ::Float64,
                     iρi ::Float64,
                     na  ::Int64,
                     nn  ::Int64,
                     nlim::Int64)

Simulate `cTfbd` according to a birth-death clads for
internal terminal branches.
"""
function _sim_cladsfbd_it(t   ::Float64,
                          λt  ::Float64,
                          μt  ::Float64,
                          αλ  ::Float64, 
                          αμ  ::Float64,
                          σλ  ::Float64,
                          σμ  ::Float64,
                          ψ   ::Vector{Float64},
                          ψts ::Vector{Float64},
                          ix  ::Int64,
                          nep ::Int64,
                          lr  ::Float64,
                          lU  ::Float64,
                          iρi ::Float64,
                          na  ::Int64,
                          nn  ::Int64,
                          nlim::Int64)

  if lU < lr && nn < nlim

    @inbounds ψi = ψ[ix]

    λi = exp(λt)
    μi = exp(μt)
    tw = cfbd_wait(λi, μi, ψi)

    # ψ epoch change
    if ix < nep
      @inbounds ψti = ψts[ix]
      if t - tw < ψti
        e0 = t - ψti
        t0, na, nn, lr = _sim_cladsfbd_it(ψti, λt, μt, αλ, αμ, σλ, σμ, 
          ψ, ψts, ix + 1, nep, lr, lU, iρi, na, nn, nlim)
        sete!(t0, e(t0) + e0)
        return t0, na, nn, lr
      end
    end

    if tw > t
      na += 1
      lr += log(iρi)
      return cTfbd(t, false, false, false, λt, μt), na, nn, lr
    end

    # if speciation
    if λevent(λi, μi, ψi)
      nn += 1
      d1, na, nn, lr = 
        _sim_cladsfbd_it(t - tw, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, lr, lU, iρi, na, nn, nlim)
      d2, na, nn, lr = 
        _sim_cladsfbd_it(t - tw, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, lr, lU, iρi, na, nn, nlim)

      return cTfbd(d1, d2, tw, false, false, false, λt, μt), na, nn, lr
    # if extinction
    elseif μevent(μi, ψi)

      return cTfbd(tw, true, false, false, λt, μt), na, nn, lr
    # if fossil sampling
    else

      return cTfbd(), na, nn, NaN
    end

  end

  return cTfbd(), na, nn, NaN
end




"""
    _sim_cladsfbd_surv(t   ::Float64,
                      λt  ::Float64,
                      μt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      surv::Bool,
                      nn  ::Int64)

Simulate `cTfbd` according to a birth-death clads.
"""
function _sim_cladsfbd_surv(t   ::Float64,
                            λt  ::Float64,
                            μt  ::Float64,
                            αλ  ::Float64,
                            αμ  ::Float64,
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
      surv, nn = 
        _sim_cladsfbd_surv(t - tw, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, surv, nn)
      surv, nn =
        _sim_cladsfbd_surv(t - tw, rnorm(λt + αλ, σλ), rnorm(μt + αμ, σμ), 
          αλ, αμ, σλ, σμ, surv, nn)

      return surv, nn
    else
      return surv, nn
    end
  end

  return true, nn
end


