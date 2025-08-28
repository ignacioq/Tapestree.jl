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




"""
    sim_gbmct(n       ::Int64;
              λ0      ::Float64 = 1.0,
              α       ::Float64 = 0.0,
              σλ      ::Float64 = 0.1,
              ϵ       ::Float64 = 0.0,
              δt      ::Float64 = 1e-3,
              init    ::Symbol  = :stem,
              nstar   ::Int64   = 2*n,
              p       ::Float64 = 5.0,
              warnings::Bool    = true,
              maxt    ::Float64 = δt*1e7)

Simulate `iTct` according to a geometric Brownian motion for birth rates and
constant turnover.
"""
function sim_gbmct(n       ::Int64;
                   λ0      ::Float64 = 1.0,
                   α       ::Float64 = 0.0,
                   σλ      ::Float64 = 0.1,
                   ϵ       ::Float64 = 0.0,
                   δt      ::Float64 = 1e-3,
                   init    ::Symbol  = :stem,
                   nstar   ::Int64   = 2*n,
                   p       ::Float64 = 5.0,
                   warnings::Bool    = true,
                   maxt    ::Float64 = δt*1e7)

  # simulate in non-recursive manner
  e0, e1, el, λs, ea, ee, na, simt =
    _sedges_gbmct(nstar, log(λ0), α, σλ, ϵ, δt, sqrt(δt), init, maxt)

  if simt >= maxt
    warnings && @warn "simulation surpassed maximum time"
    return iTct(0.0, 0.0, 0.0, false, false, Float64[])
  end

  # transform to iTree
  t = iTct(e0, e1, el, λs, ea, ee, e1[1], 1, δt)

  if ntipsalive(t) < 2
    warnings && @warn "tree went extinct"
    return t
  end

  # sample a time when species(t) == `n`
  nts = ltt(t)
  tn = times_n(n, nts)
  c  = usample(tn, p)

  if iszero(c)
    warnings && @warn "tree not sampled, try increasing `p`"
    return iTct(0.0, 0.0, 0.0, false, false, Float64[])
  else
    # cut the tree
    t = cutbottom(t, simt - c)
    return t
  end
end





"""
    _sedges_gbmct(n    ::Int64,
                  λ0   ::Float64,
                  α    ::Float64,
                  σλ   ::Float64,
                  ϵ    ::Float64,
                  δt   ::Float64,
                  srδt ::Float64,
                  init::Symbol,
                  maxt ::Float64)

Simulate `gbmct` just until hitting `n` alive species. Note that this is
a biased sample for a tree conditional on `n` species.
"""
function _sedges_gbmct(n    ::Int64,
                       λ0   ::Float64,
                       α    ::Float64,
                       σλ   ::Float64,
                       ϵ    ::Float64,
                       δt   ::Float64,
                       srδt ::Float64,
                       init::Symbol,
                       maxt ::Float64)

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
        return e0, e1, el, λs, ea, ee, na, simt
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

        # if diversification event
        if divevϵ(λm, ϵ, δt)

          #if speciation
          if λorμ(ϵ)

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

              return e0, e1, el, λs, ea, ee, na, simt
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

          #if extinction
          else
            # if tree goes extinct
            if isone(na)
              # extinct edges
              push!(ee, v)
              # delete from alive lineages
              deleteat!(ea, i)

              return e0, e1, el, λs, ea, ee, 0, simt
            end

            # extinct edges
            push!(ee, v)
            # to update alive lineages
            push!(iead, i)
            # update number of alive species
            na -= 1
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




#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Sample conditional on time
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    sim_gbmct(t   ::Float64;
              λ0  ::Float64 = 1.0,
              α   ::Float64 = 0.0,
              σλ  ::Float64 = 0.1,
              ϵ   ::Float64 = 0.2,
              δt  ::Float64 = 1e-3,
              nlim::Int64   = 10_000,
              init::Symbol  = :crown)

Simulate `iTct` according to a geometric Brownian motion for birth rates and
constant turnover.
"""
function sim_gbmct(t   ::Float64;
                   λ0  ::Float64 = 1.0,
                   α   ::Float64 = 0.0,
                   σλ  ::Float64 = 0.1,
                   ϵ   ::Float64 = 0.2,
                   δt  ::Float64 = 1e-3,
                   nlim::Int64   = 10_000,
                   init::Symbol  = :crown)

  if init === :crown
    lλ0 = log(λ0)
    d1, nn = _sim_gbmct(t, lλ0, α, σλ, ϵ, δt, sqrt(δt), 0, 1, nlim)
    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    d2, nn = _sim_gbmct(t, lλ0, α, σλ, ϵ, δt, sqrt(δt), 0, nn + 1, nlim)
    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end

    tree = iTct(d1, d2, 0.0, δt, 0.0, false, false, Float64[lλ0, lλ0])
  elseif init === :stem
    tree, nn = _sim_gbmct(t, log(λ0), α, σλ, ϵ, δt, sqrt(δt), 0, 1, nlim)

    if nn >= nlim
      @warn "maximum number of lineages surpassed"
    end
  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_gbmct(t   ::Float64,
               λt  ::Float64,
               α   ::Float64,
               σλ  ::Float64,
               ϵ   ::Float64,
               δt  ::Float64,
               srδt::Float64,
               na  ::Int64,
               nn ::Int64,
               nlim::Int64)

Simulate `iTct` according to a geometric Brownian motion for birth rates and
constant turnover, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmct(t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    na  ::Int64,
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

        if divevϵ(λm, ϵ, t)
          # if speciation
          if λorμ(ϵ)
            nn += 1
            na += 2
            return iTct(
                    iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                    iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                    bt, δt, t, false, false, λv), na, nn
          # if extinction
          else
            return iTct(bt, δt, t, true, false, λv), na, nn
          end
        end

        na  += 1
        return iTct(bt, δt, t, false, false, λv), na, nn
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)
      push!(λv, λt1)
      λm  = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(ϵ)
          nn += 1
          td1, na, nn = _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt, na, nn, nlim)
          td2, na, nn = _sim_gbmct(t, λt1, α, σλ, ϵ, δt, srδt, na, nn, nlim)

          return iTct(td1, td2, bt, δt, δt, false, false, λv), na, nn
        # if extinction
        else
          return iTct(bt, δt, δt, true, false, λv), na, nn
        end
      end

      λt = λt1
    end
  end

  return iTct(), na, nn
end




"""
    _sim_gbmct_t(t   ::Float64,
                 λt  ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 ϵ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 lr  ::Float64,
                 lU  ::Float64,
                 Iρi ::Float64,
                 na  ::Int64,
                 nn  ::Int64,
                 nlim::Int64)

Simulate `iTct` according to a geometric Brownian motion for birth rates and
constant turnover, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmct_t(t   ::Float64,
                      λt  ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      ϵ   ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      lr  ::Float64,
                      lU  ::Float64,
                      Iρi ::Float64,
                      na  ::Int64,
                      nn  ::Int64,
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

        if divevϵ(λm, ϵ, t)
          # if speciation
          if λorμ(ϵ)
            nn += 1
            na += 2
            if na === 2
              nlr = lr + log(Iρi*2.0)
            else
              nlr = lr + log(Iρi * Iρi * Float64(na)/Float64(na-2))
            end
            if nlr < lr && lU >= nlr
              return iTct(), na, nn, NaN
            else
              return iTct(iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                          iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                          bt, δt, t, false, false, λv), na, nn, nlr
            end
          # if extinction
          else
            return iTct(bt, δt, t, true, false, λv), na, nn, lr
          end
        end

        na += 1
        nlr = lr
        if na > 1
          nlr += log(Iρi * Float64(na)/Float64(na-1))
        end
        if nlr < lr && lU >= nlr
          return iTct(), na, nn, NaN
        else
          return iTct(bt, δt, t, false, false, λv), na, nn, nlr
        end
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)
      push!(λv, λt1)
      λm  = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(ϵ)
          nn += 1
          td1, na, nn, lr = 
            _sim_gbmct_t(t, λt1, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn, nlim)
          td2, na, nn, lr = 
            _sim_gbmct_t(t, λt1, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn, nlim)

          return iTct(td1, td2, bt, δt, δt, false, false, λv), na, nn, lr
        # if extinction
        else
          return iTct(bt, δt, δt, true, false, λv), na, nn, lr
        end
      end

      λt = λt1
    end
  end

  return iTct(), na, nn, NaN
end




"""
    _sim_gbmct_it(nsδt::Float64,
                  t   ::Float64,
                  λt  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  lr  ::Float64,
                  lU  ::Float64,
                  Iρi ::Float64,
                  na  ::Int64,
                  nn  ::Int64,
                  nlim::Int64)

Simulate `iTct` according to a geometric Brownian motion for birth rates and
constant turnover, starting with a non-standard δt with a limit in the number
of species.
"""
function _sim_gbmct_it(nsδt::Float64,
                       t   ::Float64,
                       λt  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       ϵ   ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       Iρi ::Float64,
                       na  ::Int64,
                       nn  ::Int64,
                       nlim::Int64)

  λv = Float64[λt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt + accerr
    t   = isapprox(t, 0.0) ? 0.0 : isapprox(t, nsδt) ? nsδt : t
    bt += t
    λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divevϵ(λm, ϵ, t)
      # if speciation
      if λorμ(ϵ)
        nn += 1
        na += 2
        lr += 2.0*log(Iρi)
        return iTct(iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                    iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                    bt, δt, t, false, false, λv), na, nn, lr
      # if extinction
      else
        return iTct(bt, δt, t, true, false, λv), na, nn, lr
      end
    end

    na += 1
    lr += log(Iρi)
   return iTct(bt, δt, t, false, false, λv), na, nn, lr
  end

  t  -= nsδt
  bt += nsδt

  λt1 = rnorm(λt + α*nsδt, sqrt(nsδt)*σλ)
  λm  = exp(0.5*(λt + λt1))
  push!(λv, λt1)

  if divevϵ(λm, ϵ, nsδt)
    # if speciation
    if λorμ(ϵ)
      nn += 1
      td1, na, nn, lr = 
        _sim_gbmct_it(t, λt1, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn, nlim)
      td2, na, nn, lr = 
        _sim_gbmct_it(t, λt1, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn, nlim)

      return iTct(td1, td2, bt, δt, nsδt, false, false, λv), na, nn, lr
    else
    # if extinction
      return iTct(bt, δt, nsδt, true, false, λv), na, nn, lr
    end
  end

  λt = λt1

  if lU < lr && nn < nlim

    ## second: standard δt
    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        λm  = exp(0.5*(λt + λt1))
        push!(λv, λt1)

        if divevϵ(λm, ϵ, t)
          # if speciation
          if λorμ(ϵ)
            nn += 1
            na += 2
            lr += 2.0*log(Iρi)
            return iTct(
                      iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                      iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                           bt, δt, t, false, false, λv), na, nn, lr
          # if extinction
          else
            return iTct(bt, δt, t, true, false, λv), na, nn, lr
          end
        end
        na += 1
        lr += log(Iρi)
        return iTct(bt, δt, t, false, false, λv), na, nn, lr
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)
      λm  = exp(0.5*(λt + λt1))
      push!(λv, λt1)

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(ϵ)
          nn += 1
          td1, na, nn, lr = 
            _sim_gbmct_it(t, λt1, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn, nlim)
          td2, na, nn, lr = 
            _sim_gbmct_it(t, λt1, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn, nlim)

          return iTct(td1, td2, bt, δt, δt, false, false, λv), na, nn, lr
        # if extinction
        else
          return iTct(bt, δt, δt, true, false, λv), na, nn, lr
        end
      end

      λt = λt1
    end
  end

  return iTct(), na, nn, NaN
end




"""
    _sim_gbmct_it(t   ::Float64,
                  λt  ::Float64,
                  α   ::Float64,
                  σλ  ::Float64,
                  ϵ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  lr  ::Float64,
                  lU  ::Float64,
                  Iρi ::Float64,
                  na  ::Int64,
                  nn  ::Int64,
                  nlim::Int64)

Simulate `iTct` according to a geometric Brownian motion for birth rates and
constant turnover, starting with a non-standard δt with a limit in the number
of species.
"""
function _sim_gbmct_it(t   ::Float64,
                       λt  ::Float64,
                       α   ::Float64,
                       σλ  ::Float64,
                       ϵ   ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       lr  ::Float64,
                       lU  ::Float64,
                       Iρi ::Float64,
                       na  ::Int64,
                       nn  ::Int64,
                       nlim::Int64)

  if lU < lr && nn < nlim

    λv = Float64[λt]
    bt = 0.0

    ## second: standard δt
    while true

      if t <= δt + accerr
        t   = isapprox(t, δt) ? δt : isapprox(t, 0.0) ? 0.0 : t
        bt += t
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)
        λm  = exp(0.5*(λt + λt1))
        push!(λv, λt1)

        if divevϵ(λm, ϵ, t)
          # if speciation
          if λorμ(ϵ)
            nn += 1
            na += 2
            lr += 2.0*log(Iρi)
            return iTct(
                      iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                      iTct(0.0, δt, 0.0, false, false, Float64[λt1, λt1]),
                           bt, δt, t, false, false, λv), na, nn, lr
          # if extinction
          else
            return iTct(bt, δt, t, true, false, λv), na, nn, lr
          end
        end
        na += 1
        lr += log(Iρi)
        return iTct(bt, δt, t, false, false, λv), na, nn, lr
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)
      λm  = exp(0.5*(λt + λt1))
      push!(λv, λt1)

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(ϵ)
          nn += 1
          td1, na, nn, lr = 
            _sim_gbmct_it(t, λt1, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn, nlim)
          td2, na, nn, lr = 
            _sim_gbmct_it(t, λt1, α, σλ, ϵ, δt, srδt, lr, lU, Iρi, na, nn, nlim)

          return iTct(td1, td2, bt, δt, δt, false, false, λv), na, nn, lr
        # if extinction
        else
          return iTct(bt, δt, δt, true, false, λv), na, nn, lr
        end
      end

      λt = λt1
    end
  end

  return iTct(), na, nn, NaN
end




"""
    _sim_gbmct_surv(t   ::Float64,
                    λt  ::Float64,
                    α   ::Float64,
                    σλ  ::Float64,
                    ϵ   ::Float64,
                    δt  ::Float64,
                    srδt::Float64,
                    surv::Bool,
                    nn ::Int64)

Simulate `iTct` according to a geometric Brownian motion for birth rates and
constant turnover, with a limit on the number lineages allowed to reach.
"""
function _sim_gbmct_surv(t   ::Float64,
                         λt  ::Float64,
                         α   ::Float64,
                         σλ  ::Float64,
                         ϵ   ::Float64,
                         δt  ::Float64,
                         srδt::Float64,
                         surv::Bool,
                         nn  ::Int64)

  if !surv && nn < 200

    while true

      if t <= δt

        t   = max(0.0,t)
        λt1 = rnorm(λt + α*t, sqrt(t)*σλ)

        # if extinction
        if rand() < ϵ*exp(0.5*(λt + λt1))*t
          return surv, nn
        end

        return true, nn
      end

      t  -= δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)
      λm  = exp(0.5*(λt + λt1))

      if divevϵ(λm, ϵ, δt)
        # if speciation
        if λorμ(ϵ)
          nn += 1
          surv, nn = _sim_gbmct_surv(t, λt1, α, σλ, ϵ, δt, srδt, surv, nn)
          surv, nn = _sim_gbmct_surv(t, λt1, α, σλ, ϵ, δt, srδt, surv, nn)

          return surv, nn
        # if extinction
        else
          return surv, nn
        end
      end

      λt = λt1
    end
  end

  return true, nn
end




"""
    divevϵ(λ::Float64, ϵ::Float64, δt::Float64)

Return true if diversification event for `ϵ` parametization.
"""
divevϵ(λ::Float64, ϵ::Float64, δt::Float64) =
  @fastmath rand() < ((1.0 + ϵ)*λ*δt)




"""
    λorμ(ϵ::Float64)

Return true if speciation event for `ϵ` parametization.
"""
λorμ(ϵ::Float64) = @fastmath rand() < (1.0/(1.0 + ϵ))



