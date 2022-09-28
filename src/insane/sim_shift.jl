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
    sim_shift(n       ::Int64;
              λ0      ::Float64 = 1.0,
              μ0      ::Float64 = 0.1,
              pshift  ::Float64 = 0.1,
              σλ      ::Float64 = 0.4,
              start   ::Symbol  = :stem,
              δt      ::Float64 = 1e-3,
              nstar   ::Int64   = 2*n,
              p       ::Float64 = 5.0,
              warnings::Bool    = true,
              maxt    ::Float64 = δt*1e7)

Simulate `iTbd` according to a shift model in speciation rates.
"""
function sim_shift(n       ::Int64;
                   λ0      ::Float64 = 1.0,
                   μ0      ::Float64 = 0.1,
                   pshift  ::Float64 = 1e-3,
                   σλ      ::Float64 = 1.0,
                   start   ::Symbol  = :stem,
                   δt      ::Float64 = 1e-3,
                   nstar   ::Int64   = 2*n,
                   p       ::Float64 = 5.0,
                   warnings::Bool    = true,
                   maxt    ::Float64 = δt*1e7)

  # simulate in non-recursive manner
  e0, e1, el, λs, μs, ea, ee, na, simt =
    _sedges_shift(nstar, log(λ0), log(μ0), pshift, σλ, δt, sqrt(δt), start, maxt)

  if simt >= maxt
    warnings && @warn "simulation surpassed maximum time"
    return iTbd(0.0, 0.0, 0.0, false, false, Float64[], Float64[])
  end

  # transform to iTree
  t = iTbd(e0, e1, el, λs, μs, ea, ee, e1[1], 1, δt)

  if iszero(ntipsalive(t))
    warnings && @warn "tree went extinct"
    return t
  end

  # sample a time when species(t) == `n`
  nt = ltt(t)
  tn = times_n(n, nt)
  c  = usample(tn, p)

  if iszero(c)
    warnings && @warn "tree not sampled, try increasing `p`"
    return iTbd(0.0, 0.0, 0.0, false, false, Float64[], Float64[])
  else
    # cut the tree
    t = cutbottom(t, simt - c)
    return t
  end
end





"""
    _sedges_shift(n     ::Int64,
                  λ0    ::Float64,
                  μ0    ::Float64,
                  pshift::Float64,
                  σλ    ::Float64,
                  δt    ::Float64,
                  srδt  ::Float64,
                  start ::Symbol,
                  maxt  ::Float64)

Simulate shift model just until hitting `n` alive species. Note that this is
a biased sample for a tree conditional on `n` species.
"""
function _sedges_shift(n     ::Int64,
                       λ0    ::Float64,
                       μ0    ::Float64,
                       pshift::Float64,
                       σλ    ::Float64,
                       δt    ::Float64,
                       srδt  ::Float64,
                       start ::Symbol,
                       maxt  ::Float64)

  # edges
  e0 = Int64[]
  e1 = Int64[]
  # edges extinct
  ee = Int64[]

  if start == :stem
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

  elseif start == :crown
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
    @error "$start does not match stem or crown"
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
        return e0, e1, el, λs, μs, ea, ee, na, simt
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
        if rand() < pshift
          λt1 = λt + randn()*σλ
        else
          λt1 = λt
        end
        μt1 = μt
        push!(λsi, λt1)
        push!(μsi, μt1)
        λm = exp(λt1)
        μm = exp(μt1)

        # if diversification event
        if divev(λm, μm, δt)

          #if speciation
          if λorμ(λm, μm)

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

                push!(λsi, λt)
                push!(μsi, μt)
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

              return e0, e1, el, λs, μs, ea, ee, na, simt
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
          else
            # if tree goes extinct
            if isone(na)
              # extinct edges
              push!(ee, v)
              # delete from alive lineages
              deleteat!(ea, i)

              return e0, e1, el, λs, μs, ea, ee, 0, simt
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





