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

Simulate tree according to `sse_g`.
"""
function simulate_sse(λ       ::Array{Float64,1},
                      μ       ::Array{Float64,1},
                      l       ::Array{Float64,1},
                      g       ::Array{Float64,1},
                      q       ::Array{Float64,1},
                      t       ::Float64;
                      δt      ::Float64 = 1e-4,
                      ast     ::Int64   = 0,
                      nspp_max::Int64   = 200_000,
                      retry_ext::Bool   = true,
                      rejectel0::Bool   = true,
                      verbose  ::Bool   = true,
                      rm_ext   ::Bool   = true)

  # make simulation
  ed, el, st, ea, ee, n, S, k = 
    simulate_edges(λ, μ, l, g, q, t, δt, ast, nspp_max)

  ne = lastindex(ee)

  if verbose
    @info "Tree with $n extant and $ne extinct species successfully simulated"

    if iszero(n)
      @warn "\n
      What would you do if an endangered animal is eating an endangered plant? \n 
      Sometimes nature is too cruel..."
      printstyled("tree went extinct... \n", color=:light_red)
    end

    if rejectel0 && in(0.0, el)
      @warn "Bad Luck! a lineage speciated at time 0.0... \n 
      rerun simulation"
    end
  end

  if n > nspp_max
    if verbose
      @warn string("Simulation surpassed the maximum of lineages allowed : ", nspp_max)
    end
  end

  if retry_ext 
    while iszero(size(ed,1)) || (rejectel0 && in(0.0, el))

      ed, el, st, ea, ee, n, S, k = 
        simulate_edges(λ, μ, l, g, q, t, δt, ast, nspp_max)

      if verbose
        @info "Tree with $n extant and $ne extinct species successfully simulated"

        if iszero(n)
          @warn "\n
          What would you do if an endangered animal is eating an endangered plant? \n 
          Sometimes nature is too cruel..."
          printstyled("tree went extinct... \n", color=:light_red)

          @info "But, don't worry, will rerun the simulation..."
        end

        if rejectel0 && in(0.0, el)
          @warn "Bad Luck! a lineage speciated at time 0.0... \n 
          rerun simulation"

          @info "But, don't worry, will rerun the simulation..."
        end
      end

      if n > nspp_max
        if verbose
          @warn string("Simulation surpassed the maximum of lineages allowed : ", nspp_max)
        end
      end
    end
  else 
    if iszero(size(ed,1))
      return Dict{Int64, Vector{Float64}}(), ed, el
    end
  end

  # tip numbers
  tN = ed[ea,2]

  # tip states
  tS = st[ea]

  # remove extinct
  if rm_ext 
    ed, el = remove_extinct(ed, el, ee)
    nt = n
  else
    nt = n + ne
  end

  # organize in postorder
  ed, tv = numberedges(ed, tN, tS)
  ed, el = postorderedges(ed, el, nt)

  ## round branch lengths
  # find out order of simulation
  i = 1
  while !isone(δt*Float64(10^i))
    i += 1
  end
  el = map(x -> round(x; digits = i+1), el)::Array{Float64,1}

  tv = states_to_values(tv, S, k)

  return tv, ed, el
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

  # edges alive
  ea = [1, 2]

  # edges extinct
  ee = Int64[]

  # edge array
  ed = zeros(Int64, nspp_max*2, 2)
  ed[ea,:] = [1 2;
              1 3]

  el = zeros(nspp_max*2)            # edge lengths
  st = zeros(Int64,nspp_max*2)      # state for each edge

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

  # make BitVector of being one area (global extinction) or not
  i1S = falses(length(S))
  for (i,s) in enumerate(S)
    if isone(length(s.g))
      i1S[i] = true
    end
  end

  # model first speciation event for daughter inheritance
  st[ea] = λtos[si][prop_sample(svλ[si], λpr[si], length(svλ[si]))]

  i0 = 3 # current first edge with 0
  n  = 2 # current number of species
  mx = 3 # current maximum node number

  ieaa = Int64[] # indexes of ea to add
  iead = Int64[] # indexes of ea to delete

  @inbounds begin

    # start simulation
    while true

      # keep track of time
      simt -= δt

      simt < 0.0 && break

      # one time step for all edges alive `ea`
      for (i,v) in enumerate(ea)

        el[v] += δt         # update edge length
        sti    = st[v]      # get state of living edge

        #=
        speciation
        =#
        if rand() < Sλpr[sti]

          n  += 1

          if n >= nspp_max 
            ed = ed[1:(i0-1),:]
            el = el[1:(i0-1)]
            st = st[1:(i0-1)]

            return ed, el, st, ea, ee, n, S, k
          end

          ### add new edges
          # start node
          ed[i0,1] = ed[i0+1,1] = ed[v,2]

          # end nodes
          ed[i0,    2] = mx + 1
          ed[i0 + 1,2] = mx + 2

          # to update living edges
          push!(iead, i)
          push!(ieaa, i0, i0 + 1)

          # update living states according to speciation event
          s1, s2 = 
            λtos[sti][prop_sample(svλ[sti], λpr[sti], length(svλ[sti]))]

          st[i0]   = s1
          st[i0+1] = s2

          # update `i0`, `n` and `mx`
          i0 += 2
          mx += 2

        #=
          extinction
        =#
        elseif rand() < Sμpr[sti]

          # if global extinction
          if i1S[sti]

            # if tree goes extinct
            if isone(n)
              return zeros(Int64,0,2), Float64[], Int64[], Int64[], Int64[], 0, Sgh[], 0
            end

            push!(ee, v)      # extinct edges
            push!(iead, i)    # to update alive lineages

            n -= 1            # update number of alive species

          # if local extinction (state change)
          else

            st[v] = 
              μtos[sti][prop_sample(svμ[sti], μpr[sti], length(svμ[sti]))]

          end

        #=
          gain event
        =#
        elseif rand() < Sgpr[sti]

          st[v] = 
            gtos[sti][prop_sample(svg[sti], gpr[sti], length(svg[sti]))]

        #=
            q event?
        =#
        elseif rand() < Sqpr[sti]

          st[v] = 
            qtos[sti][prop_sample(svq[sti], qpr[sti], length(svq[sti]))]

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

    # remove 0s
    ed = ed[1:(i0-1),:]
    el = el[1:(i0-1)]
    st = st[1:(i0-1)]
  end

  return ed, el, st, ea, ee, n, S, k
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





"""
    save_esse_sim(tv::Dict{Int64, Vector{Float64}},
                  ed::Array{Int64,2},
                  el::Array{Float64,1},
                  out_file::String)

Write esse simulation output, tree and states.
"""
function save_esse_sim(tv      ::Dict{Int64, Vector{Float64}},
                       ed      ::Array{Int64,2},
                       el      ::Array{Float64,1},
                       out_file::String)

  n    = length(tv)
  nnod = n - 1

  wt = ed[:,2] .<= n 

  tlab = Dict(i => string("t", i) for i in ed[wt,2])

  sv  = zeros(n, length(tv[1]))
  lbs = String[]
  for i in Base.OneTo(length(tv))
    push!(lbs, tlab[i])
    sv[i,:] = tv[i]
  end

  @rput ed el nnod sv lbs

  str = reval("""
              library(ape)
              t <- list(edge = ed, edge.length = el, Nnode = nnod, tip.label = lbs)
              class(t) <- "phylo"
              write.tree(t, file = "$out_file.tre")
              write.table(data.frame(lbs,sv), file = "$out_file.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
            """)
end

