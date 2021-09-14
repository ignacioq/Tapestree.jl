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
    sim_gbmce(n    ::Int64;
              λ0   ::Float64 = 1.0, 
              α    ::Float64 = 0.0, 
              σλ   ::Float64 = 0.1, 
              μ    ::Float64 = 0.0,
              δt   ::Float64 = 1e-3,
              init ::Symbol  = :crown,
              nstar::Int64   = 2*n,
              p    ::Float64 = 5.0)

Simulate `iTgbmpb` according to a pure-birth geometric Brownian motion.
"""
function sim_gbmce(n    ::Int64;
                   λ0   ::Float64 = 1.0, 
                   α    ::Float64 = 0.0, 
                   σλ   ::Float64 = 0.1, 
                   μ    ::Float64 = 0.0,
                   δt   ::Float64 = 1e-3,
                   init ::Symbol  = :crown,
                   nstar::Int64   = 2*n,
                   p    ::Float64 = 5.0)

  # simulate in non-recursive manner
  e0, e1, el, λs, ea, ee, na, simt = 
    _sedges_gbmce(nstar, log(λ0), α, σλ, μ, δt, sqrt(δt))

  # transform to iTree
  t = iTgbmce(e0, e1, el, λs, ea, ee, e1[1], 1, δt)

  if iszero(snan(t, 0))
    @warn "tree went extinct"
    return t
  end

  # sample a time when species(t) == `n`
  nt = ltt(t)
  tn = times_n(n, nt)
  c  = usample(tn, p)

  if iszero(c)
    @warn "tree not sampled, try increasing `p`"
    return iTgbmce()
  else
    # cut the tree
    t = cutbottom(t, c)
    return t
  end
end





"""
    _sedges_gbmce(n   ::Int64, 
                  λ0  ::Float64, 
                  α   ::Float64, 
                  σλ  ::Float64, 
                  μ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

Simulate `gbmce` just until hitting `n` alive species. Note that this is 
a biased sample for a tree conditional on `n` species.
"""
function _sedges_gbmce(n   ::Int64, 
                       λ0  ::Float64, 
                       α   ::Float64, 
                       σλ  ::Float64, 
                       μ   ::Float64,
                       δt  ::Float64,
                       srδt::Float64)

  # edges
  e0 = Int64[]
  e1 = Int64[]
  # edges extinct
  ee = Int64[]
  # edges alive
  ea = [1]
  # first edge
  push!(e0,1)
  push!(e1,2)
  # max index
  mxi0 = n*2
  # edge lengths
  el = [0.0]
  # lambda vector for each edge
  λs = [Float64[]]

  na = 1 # current number of alive species
  ne = 2 # current maximum node number
  ieaa = Int64[] # indexes of ea to add
  iead = Int64[] # indexes of ea to delete

  # starting speciation rate 
  push!(λs[1], λ0)
  # lastindex for each edge
  li = [1]

  # simulation time
  simt = 0.0

  @inbounds begin

    # start simulation
    while true

      # keep track of time
      simt += δt

      # one time step for all edges alive `ea`
      for (i,v) in enumerate(ea)

        λsi = λs[v]
        lii = li[v]
        λt  = λsi[lii]

        # update edge length
        el[v] += δt
        li[v] += 1

        # sample new speciation
        λt1 = rnorm(λt + α*δt, srδt*σλ)
        push!(λsi, λt1)
        λm = exp(0.5*(λt + λt1))

        # if diversification event
        if divev(λm, μ, δt)

          #if speciation 
          if λorμ(λm, μ)

            # if reached `n` species
            if n === na

              # in case of events at same time in different lineages
              if !isempty(ieaa)
                append!(ea, ieaa)
                empty!(ieaa)
              end
              if !isempty(iead)
                deleteat!(ea, iead)
                empty!(iead)
              end

              # update λs and dt for other lineages
              for ii in ea[i+1]:ea[end]
                el[ii] += δt
                λsi = λs[ii]
                lii = li[ii]
                λt  = λsi[lii]

                push!(λsi, rnorm(λt + α*δt, srδt*σλ))
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
              push!(ee, v)      # extinct edges
              return e0, e1, el, λs, ea, ee, 0, simt
            end

            push!(ee, v)      # extinct edges
            push!(iead, i)    # to update alive lineages

            na -= 1            # update number of alive species
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
    sim_gbmce(t   ::Float64;
              λ0  ::Float64 = 1.0,
              α   ::Float64 = 0.0,
              σλ  ::Float64 = 0.1,
              μ   ::Float64 = 0.2,
              δt  ::Float64 = 1e-3,
              nlim::Int64   = 10_000,
              init::Symbol  = :crown)

Simulate `iTgbmpb` according to a pure-birth geometric Brownian motion.
"""
function sim_gbmce(t   ::Float64;
                   λ0  ::Float64 = 1.0,
                   α   ::Float64 = 0.0,
                   σλ  ::Float64 = 0.1,
                   μ   ::Float64 = 0.2,
                   δt  ::Float64 = 1e-3,
                   nlim::Int64   = 10_000,
                   init::Symbol  = :crown)

  if init === :crown
    lλ0 = log(λ0)
    d1, nsp = _sim_gbmce(t, lλ0, α, σλ, μ, δt, sqrt(δt), 1, nlim)
    if nsp >= nlim 
      @warn "maximum number of lineages surpassed"
    end

    d2, nsp = _sim_gbmce(t, lλ0, α, σλ, μ, δt, sqrt(δt), 1, nlim)
    if nsp >= nlim 
      @warn "maximum number of lineages surpassed"
    end

    tree = iTgbmce(d1, d2, 0.0, δt, 0.0, false, false, Float64[lλ0, lλ0])
  elseif init === :stem
    tree, nsp = _sim_gbmce(t, log(λ0), α, σλ, μ, δt, sqrt(δt), 1, nlim)

    if nsp >= nlim 
      @warn "maximum number of lineages surpassed"
    end
  else
    @error string(init, " does not match either crown or stem")
  end

  return tree
end




"""
    _sim_gbmce(t   ::Float64,
              λt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              δt  ::Float64,
              srδt::Float64)

Simulate `iTgbmce` according to a geometric Brownian motion.
"""
function _sim_gbmce(t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

  λv = Float64[λt]
  bt = 0.0

  while true

    if t <= δt
      bt  += t

      t = max(0.0,t)
      srt = sqrt(t)
      λt1 = rnorm(λt + α*t, srt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, μ, t)
        # if speciation
        if λorμ(λm, μ)
          return iTgbmce(iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                         iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                  bt, δt, t, false, false, λv)
        # if extinction
        else
          return iTgbmce(bt, δt, t, true, false, λv)
        end
      end

      return iTgbmce(bt, δt, t, false, false, λv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt + α*δt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divev(λm, μ, δt)
      # if speciation
      if λorμ(λm, μ)
        return iTgbmce(_sim_gbmce(t, λt1, α, σλ, μ, δt, srδt), 
                       _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt), 
                bt, δt, δt, false, false, λv)
      # if extinction
      else
        return iTgbmce(bt, δt, δt, true, false, λv)
      end
    end

    λt = λt1
  end
end




"""
    _sim_gbmce(t   ::Float64,
              λt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              δt  ::Float64,
              srδt::Float64,
              nsp ::Int64,
              nlim::Int64)

Simulate `iTgbmce` according to a geometric Brownian motion with a limit
on the number lineages allowed to reach.
"""
function _sim_gbmce(t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64,
                   nsp ::Int64,
                   nlim::Int64)

  if nsp < nlim

    λv = Float64[λt]
    bt = 0.0

    while true

      if t <= δt
        bt  += t

        t  = max(0.0,t)
        srt = sqrt(t)
        λt1 = rnorm(λt + α*t, srt*σλ)

        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, μ, t)
          # if speciation
          if λorμ(λm, μ)
            nsp += 1
            return iTgbmce(
                     iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                     iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                     bt, δt, t, false, false, λv), nsp
          # if extinction
          else
            return iTgbmce(bt, δt, t, true, false, λv), nsp
          end
        end

        return iTgbmce(bt, δt, t, false, false, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, μ, δt)
        # if speciation
        if λorμ(λm, μ)
          nsp += 1
          td1, nsp = _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt, nsp, nlim)
          td2, nsp = _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt, nsp, nlim)

          return iTgbmce(td1, td2, bt, δt, δt, false, false, λv), nsp
        # if extinction
        else
          return iTgbmce(bt, δt, δt, true, false, λv), nsp
        end
      end

      λt = λt1
    end

  else
    return iTgbmce(), nsp
  end
end





"""
    _sim_gbmce(nsδt::Float64,
              t   ::Float64,
              λt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              δt  ::Float64,
              srδt::Float64)

Simulate `iTgbmce` according to a geometric Brownian motion starting 
with a non-standard δt.
"""
function _sim_gbmce(nsδt::Float64,
                   t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64)

  λv = Float64[λt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt  += t

    t   = max(0.0,t)
    srt = sqrt(t)
    λt1 = rnorm(λt + α*t, srt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divev(λm, μ, t)
      # if speciation
      if λorμ(λm, μ)
        return iTgbmce(
                 iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                 iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                 bt, δt, t, false, false, λv)
      # if extinction
      else
        return iTgbmce(bt, δt, t, true, false, λv)
      end
    end

    return iTgbmce(bt, δt, t, false, false, λv)
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt + α*nsδt, srnsδt*σλ)

  push!(λv, λt1)

  λm = exp(0.5*(λt + λt1))

  if divev(λm, μ, nsδt)
    # if speciation
    if λorμ(λm, μ)
      return iTgbmce(_sim_gbmce(t, λt1, α, σλ, μ, δt, srδt), 
                     _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt), 
              bt, δt, nsδt, false, false, λv)
    # if extinction
    else
      return iTgbmce(bt, δt, nsδt, true, false, λv)
    end
  end

  λt = λt1

  ## second: standard δt
  while true

    if t <= δt
      bt  += t

      t   = max(0.0,t)
      srt = sqrt(t)
      λt1 = rnorm(λt + α*t, srt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, μ, t)
        # if speciation
        if λorμ(λm, μ)
          return iTgbmce(
                   iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                   iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                   bt, δt, t, false, false, λv)
        # if extinction
        else
          return iTgbmce(bt, δt, t, true, false, λv)
        end
      end

      return iTgbmce(bt, δt, t, false, false, λv)
    end

    t  -= δt
    bt += δt

    λt1 = rnorm(λt + α*δt, srδt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divev(λm, μ, δt)
      # if speciation
      if λorμ(λm, μ)
        return iTgbmce(_sim_gbmce(t, λt1, α, σλ, μ, δt, srδt), 
                       _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt), 
                bt, δt, δt, false, false, λv)
      # if extinction
      else
        return iTgbmce(bt, δt, δt, true, false, λv)
      end
    end

    λt = λt1
  end
end





"""
    _sim_gbmce(nsδt::Float64,
              t   ::Float64,
              λt  ::Float64,
              α   ::Float64,
              σλ  ::Float64,
              μ   ::Float64,
              δt  ::Float64,
              srδt::Float64, 
              nsp ::Int64,
              nlim::Int64)

Simulate `iTgbmce` according to a geometric Brownian motion starting 
with a non-standard δt with a limit in the number of species.
"""
function _sim_gbmce(nsδt::Float64,
                   t   ::Float64,
                   λt  ::Float64,
                   α   ::Float64,
                   σλ  ::Float64,
                   μ   ::Float64,
                   δt  ::Float64,
                   srδt::Float64, 
                   nsp ::Int64,
                   nlim::Int64)

  λv = Float64[λt]
  bt = 0.0

  ## first: non-standard δt
  if t <= nsδt
    bt  += t

    t   = max(0.0, t)
    srt = sqrt(t)
    λt1 = rnorm(λt + α*t, srt*σλ)

    push!(λv, λt1)

    λm = exp(0.5*(λt + λt1))

    if divev(λm, μ, t)
      # if speciation
      if λorμ(λm, μ)
        nsp += 1
        return iTgbmce(
                 iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                 iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                 bt, δt, t, false, false, λv), nsp
      # if extinction
      else
        return iTgbmce(bt, δt, t, true, false, λv), nsp
      end
    end

    return iTgbmce(bt, δt, t, false, false, λv), nsp
  end

  t  -= nsδt
  bt += nsδt

  srnsδt = sqrt(nsδt)

  λt1 = rnorm(λt + α*nsδt, srnsδt*σλ)

  push!(λv, λt1)

  λm = exp(0.5*(λt + λt1))

  if divev(λm, μ, nsδt)
    # if speciation
    if λorμ(λm, μ)
      nsp += 1
      td1, nsp = _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt, nsp, nlim)
      td2, nsp = _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt, nsp, nlim)

      return iTgbmce(td1, td2, bt, δt, nsδt, false, false, λv), nsp
    else
    # if extinction
      return iTgbmce(bt, δt, nsδt, true, false, λv), nsp
    end
  end

  λt = λt1

  if nsp < nlim

    ## second: standard δt
    while true

      if t <= δt
        bt  += t

        t   = max(0.0,t)
        srt = sqrt(t)
        λt1 = rnorm(λt + α*t, srt*σλ)

        push!(λv, λt1)

        λm = exp(0.5*(λt + λt1))

        if divev(λm, μ, t)
          # if speciation
          if λorμ(λm, μ)
            nsp += 1

            return iTgbmce(
                     iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                     iTgbmce(0.0, δt, 0.0, false, false, Float64[λt1, λt1]), 
                     bt, δt, t, false, false, λv), nsp
          # if extinction
          else
            return iTgbmce(bt, δt, t, true, false, λv), nsp
          end
        end

        return iTgbmce(bt, δt, t, false, false, λv), nsp
      end

      t  -= δt
      bt += δt

      λt1 = rnorm(λt + α*δt, srδt*σλ)

      push!(λv, λt1)

      λm = exp(0.5*(λt + λt1))

      if divev(λm, μ, δt)
        # if speciation
        if λorμ(λm, μ)
          nsp += 1
          td1, nsp = _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt, nsp, nlim)
          td2, nsp = _sim_gbmce(t, λt1, α, σλ, μ, δt, srδt, nsp, nlim)

          return iTgbmce(td1, td2, bt, δt, δt, false, false, λv), nsp
        # if extinction
        else
          return iTgbmce(bt, δt, δt, true, false, λv), nsp
        end
      end

      λt = λt1
    end
  else
    return iTgbmce(), nsp
  end
end



