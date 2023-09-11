#=
Simulation for sse_g

Ignacio Quintero Mächler

t(-_-t)

January 12 2017
=#




"""
    simulate_sse(λ          ::Array{Float64,1},
                 μ          ::Array{Float64,1},
                 l          ::Array{Float64,1},
                 g          ::Array{Float64,1},
                 q          ::Array{Float64,1},
                 t          ::Float64;
                 δt         ::Float64 = 1e-4,
                 ast        ::Int64   = 0,
                 nspp_max   ::Int64   = 100_000,
                 retry_ext  ::Bool    = true,
                 rejectel0  ::Bool    = true,
                 verbose    ::Bool    = true,
                 rm_ext     ::Bool    = true,
                 states_only::Bool    = false, 
                 start      ::Symbol  = :crown)

Simulate tree according to the geographic **esse** model. The number of areas 
and hidden states is inferred from the parameter vectors, but they must be 
consistent between them and with the covariates. See tutorial for an example.

...
# Arguments
- `λ::Array{Float64,1}`: rates for within-area and between-area speciation.
- `μ::Array{Float64,1}`: per-area extinction rates when it leads to global 
extinction.
- `l::Array{Float64,1}`: per-area extinction rates when it leads to local 
extinction.
- `g::Array{Float64,1}`: colonization rates between areas. 
- `q::Array{Float64,1}`: transition rates between hidden states.
- `t::Float64`: simulation time.
- `δt::Float64 = 1e-4`: time step to perform simulations. Smaller more precise 
but more computationally expensive.
- `ast::Int64 = 0`: initial state. `0` specifies random sampling based on the 
input parameters.
- `nspp_max::Int64 = 100_000`: maximum number of species allowed to stop 
simulation.
- `retry_ext::Bool = true`: automatically restart simulation if simulation goes
extinct.
- `rejectel0::Bool = true`: reject simulations where there are edges of 0 
length.
- `verbose::Bool = true`: print messages.
- `rm_ext::Bool = true`: remove extinct taxa from output.
- `states_only::Bool = false`: if only return tip states (faster).
- `start::Symbol  = :crown`: if `crown`, starts after a speciation event with 
two lineages, if `stem`, starts with one lineage.
...

...
# Returned values
- Dictionary with tip number and corresponding state.
- Array with parent -> daughter edges.
- Array with edge lengths.
- Number of maximum species.
...
"""
function simulate_sse(λ          ::Array{Float64,1},
                      μ          ::Array{Float64,1},
                      l          ::Array{Float64,1},
                      g          ::Array{Float64,1},
                      q          ::Array{Float64,1},
                      t          ::Float64;
                      δt         ::Float64 = 1e-4,
                      ast        ::Int64   = 0,
                      nspp_max   ::Int64   = 100_000,
                      retry_ext  ::Bool    = true,
                      rejectel0  ::Bool    = true,
                      verbose    ::Bool    = true,
                      rm_ext     ::Bool    = true,
                      states_only::Bool    = false, 
                      start      ::Symbol  = :crown)

  # make simulation
  ed, el, st, ea, ee, n, S, k = 
    simulate_edges(λ, μ, l, g, q, t, δt, ast, nspp_max, start)

  ne = lastindex(ee)

  maxsp = false

  if verbose
    @info "Tree with $n extant and $ne extinct species successfully simulated"

    if n < 1
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

  if (n + ne) >= (nspp_max - (start == :crown ? 0 : 1))
    if verbose
      @warn string("Simulation surpassed the maximum of lineages allowed : ", nspp_max)
    end
    maxsp = true
  end

  if retry_ext 
    while n < 1 || (rejectel0 && in(0.0, el))

      ed, el, st, ea, ee, n, S, k = 
        simulate_edges(λ, μ, l, g, q, t, δt, ast, nspp_max, start)

      ne = lastindex(ee)

      if verbose
        @info "Tree with $n extant and $ne extinct species successfully simulated"

        if n < 1
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

      if (n + ne) > (nspp_max - (start == :crown ? 0 : 1))
        if verbose
          @warn string("Simulation surpassed the maximum of lineages allowed : ", nspp_max)
        end
        maxsp = true
      end
    end
  else 
    if n < 1
      return Dict{Int64, Vector{Float64}}(), ed, el, false
    end
  end

  # tip numbers
  tN = ed[ea,2]

  tS = st[ea]
  # tip states
  if states_only
    tv = tip_dictionary(tS)
    tv = states_to_values(tv, S, k)
    return tv, ed, el, maxsp
  end

  # remove extinct
  if rm_ext 
    if n === 1
      return Dict{Int64, Vector{Float64}}(), ed, el, false
    else
      ed, el = remove_extinct(ed, el, ee)
      nt = n
    end
  else
    nt = n + ne
  end

  # organize in postorder
  ed, tv = numberedges(ed, tN, tS)
  ed, el = postorderedges(ed, el, nt)

  ## round branch lengths
  map!(x -> round(x; digits = abs(ceil(Int64, log10(δt)))+1), el, el)

  tv = states_to_values(tv, S, k)

  return tv, ed, el, maxsp
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
                        nspp_max::Int64,
                        start   ::Symbol)

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

  # edges 
  ed = zeros(Int64, nspp_max*2, 2)
  # edges extinct
  ee = Int64[]

  if start == :crown
    # edges alive
    ea = [1, 2]

    # edge array
    ed = zeros(Int64, nspp_max*2, 2)
    ed[ea,:] = [1 2;
                1 3]
    mxi0 = (nspp_max*2 + 1)
  else
    # edges alive
    ea = [1]

    # edge array
    ed[ea,:] = [1 2]
    mxi0 = (nspp_max*2)
  end

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

  # start states
  if start == :crown
    # model first speciation event for daughter inheritance
    st[ea] = λtos[si][prop_sample(svλ[si], λpr[si], length(svλ[si]))]
    i0 = 3 # current first edge with 0
    n  = 2 # current number of species
    mx = 3 # current maximum node number
  else
    st[ea] .= si
    i0 = 2 # current first edge with 0
    n  = 1 # current number of species
    mx = 2 # current maximum node number
  end

  atol = 10.0^(floor(Int64, log10(δt)) - 2)

  ieaa = Int64[] # indexes of ea to add
  iead = Int64[] # indexes of ea to delete

  @inbounds begin

    # start simulation
    while true

      # keep track of time
      simt -= δt

      # one time step for all edges alive `ea`
      for (i,v) in enumerate(ea)

        el[v] += δt         # update edge length
        sti    = st[v]      # get state of living edge

        #=
        speciation
        =#
        if rand() < Sλpr[sti]

          if i0 === mxi0
            ed = ed[1:(i0-1),:]
            el = el[1:(i0-1)]
            st = st[1:(i0-1)]

            # in case of events at same time in different lineages
            if !isempty(ieaa)
              append!(ea, ieaa)
              empty!(ieaa)
            end
            if !isempty(iead)
              deleteat!(ea, iead)
              empty!(iead)
            end

            # update dt for other lineages
            el[ea[i+1]:ea[end]] .+= δt

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
          n  += 1

        #=
          extinction
        =#
        elseif rand() < Sμpr[sti]

          # if global extinction
          if i1S[sti]

            # if tree goes extinct
            if isone(n)
              push!(ee, v)      # extinct edges
              return ed, el, st, Int64[], ee, 0, S, k
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

      # time stop trigger
      if isapprox(simt, 0.0, atol = atol) || simt < 0.0
        break
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

