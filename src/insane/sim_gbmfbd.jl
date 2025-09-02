#=

Fossilized birth-death diffusion simulation

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#


"""
    event(λ::Float64, μ::Float64, ψ::Float64, δt::Float64)

Return true if an event for a `ifbd`.
"""
event(λ::Float64, μ::Float64, ψ::Float64, δt::Float64) =
  rand() < (λ + μ + ψ)*δt






#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Sample conditional on number of species
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    sim_gbmfbd(n       ::Int64;
               λ0      ::Float64         = 1.0,
               μ0      ::Float64         = 0.2,
               αλ      ::Float64         = 0.0,
               αμ      ::Float64         = 0.0,
               σλ      ::Float64         = 0.1,
               σμ      ::Float64         = 0.1,
               ψ       ::Vector{Float64} = [0.1],
               ψts     ::Vector{Float64} = Float64[],
               init    ::Symbol          = :stem,
               δt      ::Float64         = 1e-3,
               nstar   ::Int64           = 2*n,
               p       ::Float64         = 5.0,
               warnings::Bool            = true,
               maxt    ::Float64         = δt*1e6)

Simulate `iTfbd` according to geometric Brownian motions for birth and death
rates. Note that it will not necessarily evolve along all of the epochs.
"""
function sim_gbmfbd(n       ::Int64;
                    λ0      ::Float64         = 1.0,
                    μ0      ::Float64         = 0.2,
                    αλ      ::Float64         = 0.0,
                    αμ      ::Float64         = 0.0,
                    σλ      ::Float64         = 0.1,
                    σμ      ::Float64         = 0.1,
                    ψ       ::Vector{Float64} = [0.1],
                    ψts     ::Vector{Float64} = Float64[],
                    init    ::Symbol          = :stem,
                    δt      ::Float64         = 1e-3,
                    nstar   ::Int64           = 2*n,
                    p       ::Float64         = 5.0,
                    warnings::Bool            = true,
                    maxt    ::Float64         = δt*1e6)

  # simulate in non-recursive manner
  e0, e1, el, λs, μs, ea, ee, ef, na, simt =
    _sedges_gbmfbd(nstar, log(λ0), log(μ0), αλ, αμ, σλ, σμ, ψ, ψts, 
      δt, sqrt(δt), init, maxt)

  if simt >= maxt
    warnings && @warn "simulation surpassed maximum time"
  end

  # transform to iTree
  t = iTfbd(e0, e1, el, λs, μs, ea, ee, ef, e1[1], 1, δt)

  if iszero(ntipsalive(t))
    warnings && @warn "tree went extinct"
    return t
  end

  # sample a time when species(t) == `n`
  nts = ltt(t)
  tn  = times_n(n, nts)
  c   = usample(tn, p)

  if iszero(c)
    warnings && @warn "tree not sampled, try increasing `p`"
    return iTfbd()
  else
    # cut the tree
    t = cutbottom(t, simt - c)
    return t
  end
end




"""
    _sedges_gbmfbd(n   ::Int64,
                   λ0  ::Float64,
                   μ0  ::Float64,
                   αλ  ::Float64,
                   αμ  ::Float64,
                   σλ  ::Float64,
                   σμ  ::Float64,
                   ψ   ::Vector{Float64},
                   ψts ::Vector{Float64},
                   δt  ::Float64,
                   srδt::Float64,
                   init::Symbol,
                   maxt ::Float

Simulate `gbmfbd` just until hitting `n` alive species. Note that this is
a biased sample for a tree conditional on `n` species.
"""
function _sedges_gbmfbd(n   ::Int64,
                        λ0  ::Float64,
                        μ0  ::Float64,
                        αλ  ::Float64,
                        αμ  ::Float64,
                        σλ  ::Float64,
                        σμ  ::Float64,
                        ψ   ::Vector{Float64},
                        ψts ::Vector{Float64},
                        δt  ::Float64,
                        srδt::Float64,
                        init::Symbol,
                        maxt ::Float64)
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
    # max index
    mxi0 = n*2
    # edge lengths
    el = [0.0]
    # lambda and mu vector for each edge
    λs = [Float64[]]
    μs = [Float64[]]
    # starting speciation rate
    push!(λs[1], λ0)
    push!(μs[1], μ0)
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
    μs = [Float64[], Float64[], Float64[]]
    # starting speciation and extinction rate
    push!(λs[1], λ0, λ0)
    push!(λs[2], λ0)
    push!(λs[3], λ0)
    push!(μs[1], μ0, μ0)
    push!(μs[2], μ0)
    push!(μs[3], μ0)
    # lastindex for each edge
    li = [2, 1, 1]

    na = 2 # current number of alive species
    ne = 4 # current maximum node number

  else
    @error string(init, "does not match stem or crown")
  end

  ieaa = Int64[]                   # indexes of ea to add
  iead = Int64[]                   # indexes of ea to delete
  simt = 0.0                       # simulation time
  nep  = lastindex(ψts) + 1        # n epochs
  ix   = 1                         # epoch index ix
  ψi   = ψ[ix]                     # epoch rate
  et   = nep > 1 ? ψts[1] : -Inf   # current epoch

  @inbounds begin

    # start simulation
    while true

      # keep track of time
      simt += δt

      # time guard
      if simt > maxt
        return e0, e1, el, λs, μs, ea, ee, ef, na, simt
      end

      # select epoch's fossilization rate
      if simt - δt < et <= simt
        ix += 1
        ψi = ψ[ix]
        et = ix < nep ? ψts[ix] : -Inf
      end

      # one time step for all edges alive `ea`
      for (i,v) in enumerate(ea)

        λsi = λs[v]
        μsi = μs[v]
        lii = li[v]
        λt  = λsi[lii]
        μt  = μsi[lii]

        # update edge length
        el[v] += δt
        li[v] += 1

        # sample new speciation and extinction rates
        λt1 = rnorm(λt + αλ*δt, srδt*σλ)
        μt1 = rnorm(μt + αμ*δt, srδt*σμ)
        push!(λsi, λt1)
        push!(μsi, μt1)
        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        # if event
        if event(λm, μm, ψi, δt)

          #if speciation
          if λevent(λm, μm, ψi)

            # if reached `n` species
            if n === na

              # update λs and δt for other lineages
              for vi in ea[i+1:end]
                el[vi] += δt
                λsi = λs[vi]
                μsi = μs[vi]
                lvi = li[vi]
                λt  = λsi[lvi]
                μt  = μsi[lvi]

                push!(λsi, rnorm(λt + αλ*δt, srδt*σλ))
                push!(μsi, rnorm(μt + αμ*δt, srδt*σμ))
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

              return e0, e1, el, λs, μs, ea, ee, ef, na, simt
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
            push!(μs, [μt1], [μt1])

            # push length of vector
            push!(li, 1, 1)

            # to update living edges
            push!(iead, i)
            push!(ieaa, ne, ne + 1)

            # update `na` and `ne`
            ne += 2
            na += 1

          #if extinction
          elseif μevent(μm, ψi)
            # if tree goes extinct
            if isone(na)
              # extinct edges
              push!(ee, v)
              # delete from alive lineages
              deleteat!(ea, i)

              return e0, e1, el, λs, μs, ea, ee, ef, 0, simt
            end

            # extinct edges
            push!(ee, v)
            # to update alive lineages
            push!(iead, i)
            # update number of alive species
            na -= 1

          # if fossilization
          else
            ### add new edges
            # start node
            push!(e0, e1[v])

            # end nodes
            push!(e1, ne + 1)

            # push to edge length
            push!(el, 0.0)

            # push speciation vector
            push!(λs, [λt1])
            push!(μs, [μt1])

            # push length of vector
            push!(li, 1)

            # to update living edges
            push!(iead, i)
            push!(ieaa, ne)

            # to update fossil edges
            push!(ef, v)

            # update `na` and `ne`
            ne += 1
          end
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




# #=
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Sample conditional on time
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =#





"""
    sim_gbmfbd(t   ::Float64;
               λ0  ::Float64         = 1.0,
               μ0  ::Float64         = 0.2,
               αλ  ::Float64         = 0.0,
               αμ  ::Float64         = 0.0,
               σλ  ::Float64         = 0.1,
               σμ  ::Float64         = 0.1,
               ψ   ::Vector{Float64} = [0.1],
               ψts ::Vector{Float64} = Float64[],
               δt  ::Float64         = 1e-3,
               nlim::Int64           = 10_000,
               init::Symbol          = :crown)

Simulate `iTfbd` according to geometric Brownian motions for birth and death
rates.
"""
function sim_gbmfbd(t   ::Float64;
                    λ0  ::Float64         = 1.0,
                    μ0  ::Float64         = 0.2,
                    αλ  ::Float64         = 0.0,
                    αμ  ::Float64         = 0.0,
                    σλ  ::Float64         = 0.1,
                    σμ  ::Float64         = 0.1,
                    ψ   ::Vector{Float64} = [0.1],
                    ψts ::Vector{Float64} = Float64[],
                    δt  ::Float64         = 1e-3,
                    nlim::Int64           = 10_000,
                    init::Symbol          = :crown)

  # only include epochs where the tree occurs
  tix = findfirst(x -> x < t, ψts)
  if !isnothing(tix)
    ψ   = ψ[tix:end]
    ψts = ψts[tix:end]
  end
  nep  = lastindex(ψts) + 1

  lλ0 = log(λ0)
  lμ0 = log(μ0)

  if init === :crown
    d1, nn = _sim_gbmfbd(t, lλ0, lμ0, αλ, αμ, σλ, σμ, ψ, ψts, 1, nep, 
               δt, sqrt(δt), 1, nlim)
    nn >= nlim && @warn "maximum number of lineages surpassed"

    d2, nn = _sim_gbmfbd(t, lλ0, lμ0, αλ, αμ, σλ, σμ, ψ, ψts, 1, nep, 
               δt, sqrt(δt), nn + 1, nlim)
    nn >= nlim && @warn "maximum number of lineages surpassed"

    tree = iTfbd(d1, d2, 0.0, δt, 0.0, false, false, false,
      Float64[lλ0, lλ0], Float64[lμ0, lμ0])

  elseif init === :stem
    tree, nn = _sim_gbmfbd(t, lλ0, lμ0, αλ, αμ, σλ, σμ, ψ, ψts, 1, nep, 
                 δt, sqrt(δt), 1, nlim)
    nn >= nlim && @warn "maximum number of lineages surpassed"

  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_gbmfbd(t   ::Float64,
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
                δt  ::Float64,
                srδt::Float64,
                nn  ::Int64,
                nlim::Int64)

Simulate `iTfbd` according to geometric Brownian motions for birth and death
rates, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmfbd(t   ::Float64,
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
                     δt  ::Float64,
                     srδt::Float64,
                     nn  ::Int64,
                     nlim::Int64)

  if nn < nlim

    @inbounds begin
      ψi = ψ[ix]
      et = ix < nep ? ψts[ix] : -Inf
    end

    λv = [λt]
    μv = [μt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        t   = max(0.0, t)
        bt += t
        srt = sqrt(t)
        λt1 = rnorm(λt + αλ*t, srt*σλ)
        μt1 = rnorm(μt + αμ*t, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if 0.0 < et <= t
          @inbounds ψi = ψ[ix+1]
        end

        if event(λm, μm, ψi, t)
          # if speciation
          if λevent(λm, μm, ψi)
            nn += 1
            return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         bt, δt, t, false, false, false, λv, μv), nn
          # if extinction
          elseif μevent(μm, ψi)
            return iTfbd(bt, δt, t, true, false, false, λv, μv), nn
          # fossil sampling
          else
            return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         bt, δt, t, false, true, false, λv, μv), nn
          end
        end

        return iTfbd(bt, δt, t, false, false, false, λv, μv), nn
      end

      if t - δt < et <= t
        ix += 1
        @inbounds begin
          ψi = ψ[ix]
          et = ix < nep ? ψts[ix] : -Inf
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + αλ*δt, srδt*σλ)
      μt1 = rnorm(μt + αμ*δt, srδt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if event(λm, μm, ψi, δt)
        # if speciation
        if λevent(λm, μm, ψi)
          nn += 1
          td1, nn =
            _sim_gbmfbd(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
              nn, nlim)
          td2, nn =
            _sim_gbmfbd(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
              nn, nlim)

          return iTfbd(td1, td2, bt, δt, δt, false, false, false, λv, μv),
                 nn
        # if extinction
        elseif μevent(μm, ψi)

          return iTfbd(bt, δt, δt, true, false, false, λv, μv), nn
        else
          td1, nn =
            _sim_gbmfbd(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
              nn, nlim)

          return iTfbd(td1, bt, δt, δt, false, true, false, λv, μv), nn
        end
      end

      λt = λt1
      μt = μt1
    end
  end

  return iTfbd(), nn
end




"""
    _sim_gbmfbd_t(t   ::Float64,
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
                  δt  ::Float64,
                  srδt::Float64,
                  lr  ::Float64,
                  lU  ::Float64,
                  iρi ::Float64,
                  na  ::Int64,
                  nn  ::Int64,
                  nlim::Int64)

Simulate `iTfbd` according to geometric Brownian motions for birth and death
rates, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmfbd_t(t   ::Float64,
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
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       iρi ::Float64,
                       na  ::Int64,
                       nn  ::Int64,
                       nlim::Int64)

  if isfinite(lr) && nn < nlim

    @inbounds begin
      ψi = ψ[ix]
      et = ix < nep ? ψts[ix] : -Inf
    end

    λv = [λt]
    μv = [μt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        t   = max(0.0, t)
        bt += t
        srt = sqrt(t)
        λt1 = rnorm(λt + αλ*t, srt*σλ)
        μt1 = rnorm(μt + αμ*t, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if 0.0 < et <= t
          @inbounds ψi = ψ[ix+1]
        end

        if event(λm, μm, ψi, t)
          # if speciation
          if λevent(λm, μm, ψi)
            nn += 1
            na += 2
            if na === 2
              nlr = lr + log(iρi*2.0)
            else
              nlr = lr + log(iρi * iρi * Float64(na)/Float64(na-2))
            end
            if nlr < lr && lU >= nlr
              return iTfbd(), na, nn, NaN
            else
              return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                           iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                           bt, δt, t, false, false, false, λv, μv), na, nn, nlr
            end
          # if extinction
          elseif μevent(μm, ψi)
            return iTfbd(bt, δt, t, true, false, false, λv, μv), na, nn, lr
          # fossil sampling
          else
            return iTfbd(), na, nn, NaN
          end
        end

        na += 1
        nlr = lr
        if na > 1
          nlr += log(iρi * Float64(na)/Float64(na-1))
        end
        if nlr < lr && lU >= nlr
          return iTfbd(), na, nn, NaN
        else
          return iTfbd(bt, δt, t, false, false, false, λv, μv), na, nn, nlr
        end
      end

      if t - δt < et <= t
        ix += 1
        @inbounds begin
          ψi = ψ[ix]
          et = ix < nep ? ψts[ix] : -Inf
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + αλ*δt, srδt*σλ)
      μt1 = rnorm(μt + αμ*δt, srδt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if event(λm, μm, ψi, δt)
        # if speciation
        if λevent(λm, μm, ψi)
          nn += 1
          td1, na, nn, lr =
            _sim_gbmfbd_t(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
              lr, lU, iρi, na, nn, nlim)
          td2, na, nn, lr =
            _sim_gbmfbd_t(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
              lr, lU, iρi, na, nn, nlim)

          return iTfbd(td1, td2, bt, δt, δt, false, false, false, λv, μv),
                 na, nn, lr
        # if extinction
        elseif μevent(μm, ψi)
          return iTfbd(bt, δt, δt, true, false, false, λv, μv), na, nn, lr
        else
          return iTfbd(), na, nn, NaN
        end
      end

      λt = λt1
      μt = μt1
    end
  end

  return iTfbd(), na, nn, NaN
end




"""
    _sim_gbmfbd_i(t   ::Float64,
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
                  δt  ::Float64,
                  srδt::Float64,
                  na  ::Int64,
                  nf  ::Int64,
                  nn  ::Int64,
                  nlim::Int64)

Simulate `iTfbd` according to geometric Brownian motions for birth and death
rates, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmfbd_i(t   ::Float64,
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
                       δt  ::Float64,
                       srδt::Float64,
                       na  ::Int64,
                       nf  ::Int64,
                       nn  ::Int64,
                       nlim::Int64)

  if iszero(nf) && nn < nlim

    @inbounds begin
      ψi = ψ[ix]
      et = ix < nep ? ψts[ix] : -Inf
    end

    λv = [λt]
    μv = [μt]
    bt = 0.0

    while true

      if t - te <= δt + accerr
        dtf = t - te
        dtf = isapprox(dtf, 0.0) ? 0.0 : isapprox(dtf, δt) ? δt : dtf
        dtf = max(0.0, dtf)

        bt += dtf
        srt = sqrt(dtf)
        λt1 = rnorm(λt + αλ*dtf, srt*σλ)
        μt1 = rnorm(μt + αμ*dtf, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if t - dtf < et <= t
          @inbounds ψi = ψ[ix+1]
        end

        if event(λm, μm, ψi, dtf)
          # if speciation
          if λevent(λm, μm, ψi)
            na += 2
            nn += 1
            return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         bt, δt, dtf, false, false, false, λv, μv), na, nf, nn
          # if extinction
          elseif μevent(μm, ψi)
            return iTfbd(bt, δt, dtf, true, false, false, λv, μv), na, nf, nn
          # fossil sampling
          else
            return iTfbd(), na, 1, nn
          end
        end

        na += 1
        return iTfbd(bt, δt, dtf, false, false, false, λv, μv), na, nf, nn
      end

      if t - δt < et <= t
        ix += 1
        @inbounds begin
          ψi = ψ[ix]
          et = ix < nep ? ψts[ix] : -Inf
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + αλ*δt, srδt*σλ)
      μt1 = rnorm(μt + αμ*δt, srδt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if event(λm, μm, ψi, δt)
        # if speciation
        if λevent(λm, μm, ψi)
          nn += 1
          td1, na, nf, nn =
            _sim_gbmfbd_i(t, te, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, 
              δt, srδt, na, nf, nn, nlim)
          td2, na, nf, nn =
            _sim_gbmfbd_i(t, te, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, 
              δt, srδt, na, nf, nn, nlim)

          return iTfbd(td1, td2, bt, δt, δt, false, false, false, λv, μv),
                 na, nf, nn
        # if extinction
        elseif μevent(μm, ψi)

          return iTfbd(bt, δt, δt, true, false, false, λv, μv), na, nf, nn
        else
          return iTfbd(), na, 1, nn
        end
      end

      λt = λt1
      μt = μt1
    end
  end

  return iTfbd(), na, 1, nn
end




"""
    _sim_gbmfbd_it(nsδt::Float64,
                   t   ::Float64,
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
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64)

Simulate `iTfbd` according to geometric Brownian motions for birth and death
rates, starting with a non-standard δt with a limit in the number of species.
"""
function _sim_gbmfbd_it(nsδt::Float64,
                        t   ::Float64,
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
                        δt  ::Float64,
                        srδt::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64,
                        nlim::Int64)

  λv = [λt]
  μv = [μt]
  bt = 0.0

  @inbounds begin
    ψi = ψ[ix]
    et = ix < nep ? ψts[ix] : -Inf
  end

 ## first: non-standard δt
  if t <= nsδt + accerr
    t   = isapprox(t, 0.0) ? 0.0 : isapprox(t, nsδt) ? nsδt : t
    t   = max(0.0, t)
    bt += t
    srt = sqrt(t)
    λt1 = rnorm(λt + αλ*t, srt*σλ)
    μt1 = rnorm(μt + αμ*t, srt*σμ)

    push!(λv, λt1)
    push!(μv, μt1)

    λm = exp(0.5*(λt + λt1))
    μm = exp(0.5*(μt + μt1))

    if 0.0 < et <= t
      @inbounds ψi = ψ[ix+1]
    end

    if event(λm, μm, ψi, t)
      # if speciation
      if λevent(λm, μm, ψi)
        nn += 1
        na += 2
        lr += 2.0*log(iρi)
        return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                           [λt1, λt1], [μt1, μt1]),
                     iTfbd(0.0, δt, 0.0, false, false, false,
                           [λt1, λt1], [μt1, μt1]),
                     bt, δt, t, false, false, false, λv, μv), na, nn, lr
      # if extinction
      elseif μevent(μm, ψi)
        return iTfbd(bt, δt, t, true, false, false, λv, μv), na, nn, lr
      # if fossil sampling
      else
        return iTfbd(), na, nn, NaN
      end
    end

    na += 1
    lr += log(iρi)
    return iTfbd(bt, δt, t, false, false, false, λv, μv), na, nn, lr
  end

  if t - nsδt < et <= t
    ix += 1
    @inbounds begin
      ψi = ψ[ix]
      et = ix < nep ? ψts[ix] : -Inf
    end
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt + αλ*nsδt, srnsδt*σλ)
  μt1 = rnorm(μt + αμ*nsδt, srnsδt*σμ)

  push!(λv, λt1)
  push!(μv, μt1)

  λm = exp(0.5*(λt + λt1))
  μm = exp(0.5*(μt + μt1))

  if event(λm, μm, ψi, nsδt)
    # if speciation
    if λevent(λm, μm, ψi)
      nn += 1
      td1, na, nn, lr =
        _sim_gbmfbd_it(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, iρi, na, nn, nlim)
      td2, na, nn, lr =
        _sim_gbmfbd_it(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, δt, srδt,
          lr, lU, iρi, na, nn, nlim)

      return iTfbd(td1, td2, bt, δt, nsδt, false, false, false, λv, μv),
             na, nn, lr
    # if extinction
    elseif μevent(μm, ψi)
      return iTfbd(bt, δt, nsδt, true, false, false, λv, μv), na, nn, lr
    # if fossil sampling
    else
      return iTfbd(), na, nn, NaN
    end
  end

  λt = λt1
  μt = μt1

  if lU < lr &&  nn < nlim

    ## second: standard δt
    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        t   = max(0.0, t)
        bt += t
        srt = sqrt(t)
        λt1 = rnorm(λt + αλ*t, srt*σλ)
        μt1 = rnorm(μt + αμ*t, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if 0.0 < et <= t
          @inbounds ψi = ψ[ix+1]
        end

        if event(λm, μm, ψi, t)
          # if speciation
          if λevent(λm, μm, ψi)
            nn += 1
            na += 2
            lr += 2.0*log(iρi)
            return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         bt, δt, t, false, false, false, λv, μv), na, nn, lr
          # if extinction
          elseif μevent(μm, ψi)
            return iTfbd(bt, δt, t, true, false, false, λv, μv), na, nn, lr
          # if fossil sampling
          else
            return iTfbd(), na, nn, NaN
          end
        end

        na += 1
        lr += log(iρi)
        return iTfbd(bt, δt, t, false, false, false, λv, μv), na, nn, lr
      end

      if t - δt < et <= t
        ix += 1
        @inbounds begin
          ψi = ψ[ix]
          et = ix < nep ? ψts[ix] : -Inf
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + αλ*δt, srδt*σλ)
      μt1 = rnorm(μt + αμ*δt, srδt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if event(λm, μm, ψi, δt)
        # if speciation
        if λevent(λm, μm, ψi)
          nn += 1
          td1, na, nn, lr =
            _sim_gbmfbd_it(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, 
              δt, srδt, lr, lU, iρi, na, nn, nlim)
          td2, na, nn, lr =
            _sim_gbmfbd_it(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, 
              δt, srδt, lr, lU, iρi, na, nn, nlim)

          return iTfbd(td1, td2, bt, δt, δt, false, false, false, λv, μv),
                 na, nn, lr
        # if extinction
        elseif μevent(μm, ψi)
          return iTfbd(bt, δt, δt, true, false, false, λv, μv), na, nn, lr
        else
          return iTfbd(), na, nn, NaN
        end
      end

      λt = λt1
      μt = μt1
    end
  end

  return iTfbd(), na, nn, NaN
end




"""
    _sim_gbmfbd_it(t   ::Float64,
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
                   δt  ::Float64,
                   srδt::Float64,
                   lr  ::Float64,
                   lU  ::Float64,
                   iρi ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64)

Simulate `iTfbd` according to geometric Brownian motions for birth and death
rates, starting with a non-standard δt with a limit in the number of species.
"""
function _sim_gbmfbd_it(t   ::Float64,
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
                        δt  ::Float64,
                        srδt::Float64,
                        lr  ::Float64,
                        lU  ::Float64,
                        iρi ::Float64,
                        na  ::Int64,
                        nn  ::Int64,
                        nlim::Int64)

  if lU < lr &&  nn < nlim

    @inbounds begin
      ψi = ψ[ix]
      et = ix < nep ? ψts[ix] : -Inf
    end

    λv = [λt]
    μv = [μt]
    bt = 0.0

    ## second: standard δt
    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        t   = max(0.0, t)
        bt += t
        srt = sqrt(t)
        λt1 = rnorm(λt + αλ*t, srt*σλ)
        μt1 = rnorm(μt + αμ*t, srt*σμ)

        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if event(λm, μm, ψi, t)
          # if speciation
          if λevent(λm, μm, ψi)
            nn += 1
            na += 2
            lr += 2.0*log(iρi)
            return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         bt, δt, t, false, false, false, λv, μv), na, nn, lr
          # if extinction
          elseif μevent(μm, ψi)
            return iTfbd(bt, δt, t, true, false, false, λv, μv), na, nn, lr
          # if fossil sampling
          else
            return iTfbd(), na, nn, NaN
          end
        end

        na += 1
        lr += log(iρi)
        return iTfbd(bt, δt, t, false, false, false, λv, μv), na, nn, lr
      end

      if t - δt < et <= t
        ix += 1
        @inbounds begin
          ψi = ψ[ix]
          et = ix < nep ? ψts[ix] : -Inf
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + αλ*δt, srδt*σλ)
      μt1 = rnorm(μt + αμ*δt, srδt*σμ)

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if event(λm, μm, ψi, δt)
        # if speciation
        if λevent(λm, μm, ψi)
          nn += 1
          td1, na, nn, lr =
            _sim_gbmfbd_it(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, 
              δt, srδt, lr, lU, iρi, na, nn, nlim)
          td2, na, nn, lr =
            _sim_gbmfbd_it(t, λt1, μt1, αλ, αμ, σλ, σμ, ψ, ψts, ix, nep, 
              δt, srδt, lr, lU, iρi, na, nn, nlim)

          return iTfbd(td1, td2, bt, δt, δt, false, false, false, λv, μv),
                 na, nn, lr
        # if extinction
        elseif μevent(μm, ψi)
          return iTfbd(bt, δt, δt, true, false, false, λv, μv), na, nn, lr
        else
          return iTfbd(), na, nn, NaN
        end
      end

      λt = λt1
      μt = μt1
    end
  end

  return iTfbd(), na, nn, NaN
end





"""
    _sim_gbmfbd_surv(t   ::Float64,
                     λt  ::Float64,
                     μt  ::Float64,
                     αλ  ::Float64,
                     αμ  ::Float64,
                     σλ  ::Float64,
                     σμ  ::Float64,
                     δt  ::Float64,
                     srδt::Float64,
                     surv::Bool,
                     nn  ::Int64)

Returns if survived under the `fbdd` model.
"""
function _sim_gbmfbd_surv(t   ::Float64,
                          λt  ::Float64,
                          μt  ::Float64,
                          αλ  ::Float64,
                          αμ  ::Float64,
                          σλ  ::Float64,
                          σμ  ::Float64,
                          δt  ::Float64,
                          srδt::Float64,
                          surv::Bool,
                          nn  ::Int64)

  if !surv && nn < 200

    while true

      if t <= δt
        t   = max(0.0, t)
        μt1 = rnorm(μt + αμ*δt, sqrt(t)*σμ)

        # if extinction
        if rand() < exp(0.5*(μt + μt1))*t
          return surv, nn
        else
          return true, nn
        end

        return true, nn
      end

      t  -= δt
      λt1 = rnorm(λt + αλ*δt, srδt*σλ)
      μt1 = rnorm(μt + αμ*δt, srδt*σμ)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if divev(λm, μm, δt)
        # if speciation
        if λorμ(λm, μm)
          nn += 1
          surv, nn =
            _sim_gbmfbd_surv(t, λt1, μt1, αλ, αμ, σλ, σμ, δt, srδt, surv, nn)
          surv, nn =
            _sim_gbmfbd_surv(t, λt1, μt1, αλ, αμ, σλ, σμ, δt, srδt, surv, nn)

          return surv, nn
        # if extinction
        else
          return surv, nn
        end
      end

      λt = λt1
      μt = μt1
    end
  end

  return true, nn
end




"""
    _sim_gbmfbd_fx(t   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   ix  ::Int64,
                   tz  ::Vector{Float64},
                   zλ  ::Vector{Float64},
                   zμ  ::Vector{Float64},
                   ψ   ::Float64,
                   na  ::Int64,
                   nn  ::Int64,
                   nlim::Int64)

Simulate `iTfbd` where `λ(t)` & `μ(t)` follow `zλ` and `zμ`.
"""
function _sim_gbmfbd_fx(t   ::Float64,
                        δt  ::Float64,
                        srδt::Float64,
                        ix  ::Int64,
                        tz  ::Vector{Float64},
                        zλ  ::Vector{Float64},
                        zμ  ::Vector{Float64},
                        ψ   ::Float64,
                        na  ::Int64,
                        nn  ::Int64,
                        nlim::Int64)

  if nn < nlim

    λt = linpred(t, tz[ix], tz[ix+1], zλ[ix], zλ[ix+1])
    μt = linpred(t, tz[ix], tz[ix+1], zμ[ix], zμ[ix+1])
    λv = Float64[λt]
    μv = Float64[μt]
    bt = 0.0

    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        t   = max(0.0,t)
        bt += t
        srt = sqrt(t)

        while 0.0 < tz[ix]
          ix += 1
        end
        ix -= 1
        λt1 = linpred(t, tz[ix], tz[ix+1], zλ[ix], zλ[ix+1])
        μt1 = linpred(t, tz[ix], tz[ix+1], zμ[ix], zμ[ix+1])

        push!(λv, λt1)
        push!(μv, μt1)

        λm = exp(0.5*(λt + λt1))
        μm = exp(0.5*(μt + μt1))

        if event(λm, μm, ψ, t)
          # if speciation
          if λevent(λm, μm, ψ)
            nn += 1
            na += 2
            return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         bt, δt, t, false, false, false, λv, μv), na, nn
          # if extinction
          elseif μevent(μm, ψ)
            return iTfbd(bt, δt, t, true, false, false, λv, μv), na, nn
          # fossil sampling
          else
            return iTfbd(iTfbd(0.0, δt, 0.0, false, false, false,
                               [λt1, λt1], [μt1, μt1]),
                         bt, δt, t, false, true, false, λv, μv), na, nn
          end
        end

        na += 1
        return iTfbd(bt, δt, t, false, false, false, λv, μv), na, nn
      end

      t  -= δt
      bt += δt

      while t < tz[ix]
        ix += 1
      end
      ix -= 1
      λt1 = linpred(t, tz[ix], tz[ix+1], zλ[ix], zλ[ix+1])
      μt1 = linpred(t, tz[ix], tz[ix+1], zμ[ix], zμ[ix+1])

      push!(λv, λt1)
      push!(μv, μt1)

      λm = exp(0.5*(λt + λt1))
      μm = exp(0.5*(μt + μt1))

      if event(λm, μm, ψ, δt)
        # if speciation
        if λevent(λm, μm, ψ)
          nn += 1
          td1, na, nn =
            _sim_gbmfbd_fx(t, δt, srδt, ix, tz, zλ, zμ, ψ, na, nn, nlim)
          td2, na, nn =
            _sim_gbmfbd_fx(t, δt, srδt, ix, tz, zλ, zμ, ψ, na, nn, nlim)

          return iTfbd(td1, td2, bt, δt, δt, false, false, false, λv, μv),
                 na, nn
        # if extinction
        elseif μevent(μm, ψ)

          return iTfbd(bt, δt, δt, true, false, false, λv, μv), na, nn
        else
          td1, na, nn =
            _sim_gbmfbd_fx(t, δt, srδt, ix, tz, zλ, zμ, ψ, na, nn, nlim)

          return iTfbd(td1, bt, δt, δt, false, true, false, λv, μv), na, nn
        end
      end

      λt = λt1
      μt = μt1
    end
  end

  return iTfbd(), na, nn
end

