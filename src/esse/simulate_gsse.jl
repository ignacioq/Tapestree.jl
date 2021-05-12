#=
Simulation for sse_g

Ignacio Quintero Mächler

t(-_-t)

January 12 2017
=#




"""
    simulate_sse(λ       ::Array{Float64,1},
                 μ       ::Array{Float64,1},
                 l       ::Array{Float64,1},
                 g       ::Array{Float64,1},
                 q       ::Array{Float64,1},
                 t       ::Float64;
                 δt      ::Float64 = 1e-4,
                 ast     ::Int64   = 0,
                 nspp_min::Int64   = 1,
                 nspp_max::Int64   = 200_000,
                 retry_ext::Bool   = true,
                 rejectel0::Bool   = false)

Simulate tree according to `SSE.g`.
"""
function simulate_sse(λ       ::Array{Float64,1},
                      μ       ::Array{Float64,1},
                      l       ::Array{Float64,1},
                      g       ::Array{Float64,1},
                      q       ::Array{Float64,1},
                      t       ::Float64;
                      δt      ::Float64 = 1e-4,
                      ast     ::Int64   = 0,
                      nspp_min::Int64   = 1,
                      nspp_max::Int64   = 200_000,
                      retry_ext::Bool   = true,
                      rejectel0::Bool   = true,
                      verbose  ::Bool   = true)

  # make simulation
  ed, el, st, n, S, k = 
    simulate_edges(λ, μ, l, g, q, t, δt, ast, nspp_max)

  if verbose

    @info "Tree with $n extant species successfully simulated"

    if iszero(n)
      @warn "\n
      What would you do if an endangered animal is eating an endangered plant? \n 
      Sometimes nature is too cruel..."
      printstyled("tree went extinct... \n", color=:light_red)
    end

    if n > nspp_max
      @warn string("Simulation surpassed the maximum of lineages allowed : ", nspp_max)
    end

    if rejectel0 && in(0.0, el)
      @warn "Bad Luck! a lineage speciated at time 0.0... \n 
      rerun simulation"
    end
  end

  if retry_ext 
    while iszero(size(ed,1)) || (rejectel0 && in(0.0, el)) ||
          n < nspp_min       || n > nspp_max 
      ed, el, st, n, S, k = 
        simulate_edges(λ, μ, l, g, q, t, δt, ast, nspp_max)

      if verbose
        @info "Tree with $n extant species successfully simulated"

        if iszero(n)
          @warn "\n
          What would you do if an endangered animal is eating an endangered plant? \n 
          Sometimes nature is too cruel..."
          printstyled("tree went extinct... \n", color=:light_red)

          @info "But, don't worry, will rerun the simulation..."
        end

        if n > nspp_max
          @warn string("Simulation surpassed the maximum of lineages allowed : ", nspp_max)

          @info "But, don't worry, will rerun the simulation..."
        end

        if rejectel0 && in(0.0, el)
          @warn "Bad Luck! a lineage speciated at time 0.0... \n 
          rerun simulation"

          @info "But, don't worry, will rerun the simulation..."
        end
      end
    end
  else 
    if iszero(size(ed,1))
      return Dict{Int64, Vector{Float64}}(), ed, el
    end
  end

  # organize in postorder
  ed     = numberedges(ed, n)
  ed, el = postorderedges(ed, el, n)

  ## round branch lengths
  # find out order of simulation
  i = 1
  while !isone(δt*Float64(10^i))
    i += 1
  end

  el = map(x -> round(x; digits = i+1), el)::Array{Float64,1}

  # organize states
  tip_val = Dict(i => st[i] for i = 1:n)

  tip_val = states_to_values(tip_val, S, k)

  return tip_val, ed, el
end




"""
    simulate_edges(λ       ::Array{Float64,1},
                   μ       ::Array{Float64,1},
                   l       ::Array{Float64,1},
                   g       ::Array{Float64,1},
                   q       ::Array{Float64,1},
                   simt    ::Float64, 
                   δt      ::Float64,
                   si      ::Int64,
                   nspp_max::Int64)

Simulate edges tree according to `sse_g`.
"""
function simulate_edges(λ       ::Array{Float64,1},
                        μ       ::Array{Float64,1},
                        l       ::Array{Float64,1},
                        g       ::Array{Float64,1},
                        q       ::Array{Float64,1},
                        simt    ::Float64, 
                        δt      ::Float64,
                        si      ::Int64,
                        nspp_max::Int64)

  h = div(isqrt(length(q)*4 + 1) + 1, 2)
  k = div(length(l), h)

  # create states
  S = create_states(k, h)

  # number of states
  ns = length(S)

  # areas & hidden states
  as = 1:k
  hs = 0:(h-1)

  # if initial state not assigned
  if iszero(si)
    # vector of initial state probabilities
    isp = Array{Float64,1}(undef, length(S))

    # assign initial state probabilities
    init_states_pr!(isp, λ, μ, l, g, q, S, k, h, as, hs)

    # preallocate vector of individual area probabilities 
    spr = Array{Float64,1}(undef, ns)

    # sample initial state
    si = prop_sample(spr, isp, ns)
  end

  # n = 2 starting lineages (assume no stem branch)
  n = 2

  # edges alive
  ea = [1, 2]

  # edge array
  ed = zeros(Int64, nspp_max*2, 2)
  ed[ea,:] = [1 2;
              1 3]

  # edge lengths
  el = zeros(nspp_max*2)

  # make probability vectors
  λpr, μpr, gpr, qpr, Sλpr, Sμpr, Sgpr, Sqpr = 
   event_probs(λ, μ, l, g, q, S, k, h, δt, ns, as, hs)

  # preallocate states probabilities for specific lengths of each event
  svλ = Array{Float64,1}[]
  svμ = Array{Float64,1}[]
  svg = Array{Float64,1}[]
  svq = Array{Float64,1}[]
  for i in Base.OneTo(ns)
    push!(svλ, zeros(length(λpr[i])))
    push!(svμ, zeros(length(μpr[i])))
    push!(svg, zeros(length(gpr[i])))
    push!(svq, zeros(length(qpr[i])))
  end

  # make state change vectors
  gtos, μtos, qtos, λtos = 
    makecorresschg(gpr, μpr, qpr, λpr, S, as, hs, k, ns)

  # model first speciation event for daughter inheritance
  st = λtos[si][prop_sample(svλ[si], λpr[si], length(svλ[si]))]

  # random assignment to daughters
  if rand() < 0.5
    reverse!(st)
  end

  # start simulation
  while true

    # keep track of time
    simt -= δt

    # one time step
    for i in Base.OneTo(n)

      # update edge length
      el[ea[i]] += δt

      sti = st[i]

      #=
        gain event?
      =#
      if rand() < Sgpr[sti]
        # sample from the transition probabilities
        @inbounds st[i] = 
          gtos[sti][prop_sample(svg[sti], gpr[sti], length(svg[sti]))]

        # no more events at this time
        continue
      end

      #=
        μ or loss event?
      =#
      if rand() < Sμpr[sti]

        # if extinction
        if isone(length(S[sti].g))

          # update time in other extant lineages
          el[ea[(i+1):n]] .+= δt

          # node to remove
          nod = ed[ea[i],1]

          # top or bottom edge of node
          top = nod === ed[ea[i]+1, 1] 

          # remove edge (only preserving persisting species)
          ned = findfirst(isequal(0), ed)[1]

          # if goes extinct
          if ned === 3
            return zeros(Int64,0,2), Float64[], Int64[], 0, Sgh[], 0
          end

          @views ed[ea[i]:(ned-1),:] = ed[(ea[i]+1):ned,:]

          # remove edge length
          @views el[ea[i]:(ned-1)] = el[(ea[i]+1):ned]

          # update alive lineages
          deleteat!(ea,i)
          ea[i:end] .-= 1

          # remove node and extra edge
          @views pr = findfirst(isequal(nod), ed[:,1])::Int64
          
          @views so = findfirst(isequal(nod), ed[:,2])

          if isnothing(so)
            da = 0
          else
            da = so::Int64
          end

          if !iszero(da)
            ed[da,2] = ed[pr,2]
            el[da]  += el[pr]
          end

          # remove intervening edge
          @views ed[pr:(ned-2),:] = ed[(pr+1):(ned-1),:]

          # remove edge lengths of intervening node
          @views el[pr:(ned-2)] = el[(pr+1):(ned-1)]

          ## update alive lineages
          # is the remaining node terminal?
          @views if iszero(da) || in(ed[da,2], ed[:,1])
            ea[i:end] .-= 1
          else 
            ea[findfirst(isequal(pr),ea)::Int64] = da
            sort!(ea)
            if top
              ea[(i+1):end] .-= 1
            else
              ea[i:end]     .-= 1
            end
          end

          # update living species states
          deleteat!(st,i)

          # update n species
          n = lastindex(ea)

          # break loop
          break

        # if local extinction (state change)
        else
          @inbounds st[i] = 
            μtos[sti][prop_sample(svμ[sti], μpr[sti], length(svμ[sti]))]

          # no more events at this time
          continue
        end

      end

      #=
          q event?
      =#
      if rand() < Sqpr[sti]
        @inbounds st[i] = 
          qtos[sti][prop_sample(svq[sti], qpr[sti], length(svq[sti]))]

        # no more events at this time
        continue
      end

      #=
          λ event?
      =#
      if rand() < Sλpr[sti]

        # update time in other extant lineages
        el[ea[(i+1):n]] .+= δt

        ### add new edges
        # start node
        ed[ea[end] + 1,1] = ed[ea[end] + 2,1] = ed[ea[i],2]

        # end nodes
        ed[ea[end] + 1,2] = maximum(ed)+1
        ed[ea[end] + 2,2] = maximum(ed)+1

        # update living edges
        push!(ea, ea[end]+1, ea[end]+2)
        deleteat!(ea, i)

        # update living states according to speciation event
        nst1, nst2 = 
          λtos[sti][prop_sample(svλ[sti], λpr[sti], length(svλ[sti]))]

        # random assignment to daughters
        if rand() < 0.5
          push!(st,nst1,nst2)
        else
          push!(st,nst2,nst1)
        end

        deleteat!(st,i)

        # update number of species alive
        n = lastindex(ea)

        # no more events for any of the remaining lineages
        break
      end

    end

    simt < 0.0 && break

    if n > nspp_max 
      return zeros(Int64,0,2), Float64[], Int64[], n, Sgh[], 0
    end
  end

  # remove 0s
  ed = ed[1:(2n-2),:]
  el = el[1:(2n-2)]

  return ed, el, st, n, S, k
end





"""
    init_states_pr!(isp  ::Array{Float64,1},
                    λ    ::Array{Float64,1},
                    μ    ::Array{Float64,1},
                    l    ::Array{Float64,1},
                    g    ::Array{Float64,1},
                    q    ::Array{Float64,1},
                    S    ::Array{Sgh,1},
                    k    ::Int64,
                    h    ::Int64,
                    as   ::UnitRange{Int64},
                    hs   ::UnitRange{Int64})

 Estimate starting transition probabilities for initial states
"""
function init_states_pr!(isp  ::Array{Float64,1},
                         λ    ::Array{Float64,1},
                         μ    ::Array{Float64,1},
                         l    ::Array{Float64,1},
                         g    ::Array{Float64,1},
                         q    ::Array{Float64,1},
                         S    ::Array{Sgh,1},
                         k    ::Int64,
                         h    ::Int64,
                         as   ::UnitRange{Int64},
                         hs   ::UnitRange{Int64})

  for si in Base.OneTo(length(S))

    s = S[si]

    # allopatric speciation rates
    sr = 0.0
    if length(s.g) != k
      sr += λ[(s.h+1)*(k+1)]
    end

    # loss rates
    lr = 0.0
    for i = setdiff(as, s.g)
      lr += l[s.h*k+i]
    end

    # gain rates
    gr = 0.0
    if length(s.g) > 1
      for ta = s.g, fa = setdiff(s.g, ta)
        gr += g[s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)]
      end
    end

    # hidden states
    hr = 0.0
    for fa = setdiff(hs, s.h)
        hr += q[fa*(h-1) + (s.h > fa ? s.h : s.h+1)]
    end

    # add all and hidden rates
    isp[si] = sr + lr + gr + hr
  end

  return nothing
end





"""
    event_probs(λ    ::Array{Float64,1},
                μ    ::Array{Float64,1},
                l    ::Array{Float64,1},
                g    ::Array{Float64,1},
                q    ::Array{Float64,1},
                β    ::Array{Float64,1},
                r    ::Array{Float64,1},
                S    ::Array{Sgh,1},
                k    ::Int64,
                h    ::Int64,
                ny   ::Int64,
                model::NTuple{3,Bool},
                md   ::Bool, 
                δt   ::Float64,
                ns   ::Int64,
                as   ::UnitRange{Int64},
                hs   ::UnitRange{Int64})

Create event probabilities for each state given the input parameters.
"""
function event_probs(λ    ::Array{Float64,1},
                     μ    ::Array{Float64,1},
                     l    ::Array{Float64,1},
                     g    ::Array{Float64,1},
                     q    ::Array{Float64,1},
                     S    ::Array{Sgh,1},
                     k    ::Int64,
                     h    ::Int64,
                     δt   ::Float64,
                     ns   ::Int64,
                     as   ::UnitRange{Int64},
                     hs   ::UnitRange{Int64})

  @inbounds begin
    ### make fixed vectors of approximate probabilities for each state
    ## λ
    λpr = Array{Float64,1}[]
    for s in S
      push!(λpr, zeros(length(s.g) + div(length(vicsubsets(s.g)),2)))
    end
    # fill it up
    for si in Base.OneTo(ns)
      s = S[si]
      na = 0
      for a in s.g
        na += 1
        λpr[si][na] = λ[(k+1)*s.h + a]*δt
      end
      # between-area speciation
      if na > 1
        λpr[si][na+1:end] .= λ[k+1 + (k+1)*s.h]*δt
      end
    end

    ## μ and loss
    μpr = Array{Float64,1}[]
    for s in S
      push!(μpr, zeros(length(s.g)))
    end
    # fill it up
    for si in Base.OneTo(ns)
      s = S[si]
      if isone(length(s.g))
        for a in s.g
          μpr[si][1] = μ[k*s.h + a]*δt
        end
      else
        na = 0
        for a in s.g
          na += 1
          μpr[si][na] = l[k*s.h + a]*δt
        end
      end
    end

    ## gain
    gpr = Array{Float64,1}[]
    for s in S
      push!(gpr, zeros(length(setdiff(as,s.g))))
    end
    # fill it up
    for si in Base.OneTo(ns)
      s = S[si]
      for (i,ta) = enumerate(setdiff(as,s.g)), fa = s.g
        gpr[si][i] += 
          g[s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)]*δt
      end
    end

    ## hidden states
    qpr = Array{Float64,1}[]
    for s in S
      push!(qpr, zeros(h-1))
    end
    # fill it up
    for si in Base.OneTo(ns)
      s = S[si]
       for (i, th) = enumerate(setdiff(hs, s.h))
        qpr[si][i] = q[s.h*(h-1) + (s.h > th ? s.h : s.h+1)]*δt
      end
    end

    # make λ, μ and g probability sums
    Sλpr = zeros(ns)
    Sμpr = zeros(ns)
    Sgpr = zeros(ns)
    Sqpr = zeros(ns)
    for i in Base.OneTo(ns)
      Sλpr[i] = sum(λpr[i])
      Sμpr[i] = sum(μpr[i])
      Sgpr[i] = sum(gpr[i])
      Sqpr[i] = sum(qpr[i])
    end
  end

  return λpr, μpr, gpr, qpr, Sλpr, Sμpr, Sgpr, Sqpr
end




