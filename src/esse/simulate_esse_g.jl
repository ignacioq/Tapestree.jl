#=
Simulation for ESSE

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
                 β          ::Array{Float64,1},
                 t          ::Float64,
                 x          ::Array{Float64,1},
                 y          ::Array{Float64,N},
                 cov_mod    ::Tuple{Vararg{String}};
                 δt         ::Float64 = 1e-4,
                 ast        ::Int64   = 0,
                 nspp_max   ::Int64   = 100_000,
                 retry_ext  ::Bool    = true,
                 rejectel0  ::Bool    = true,
                 verbose    ::Bool    = true,
                 rm_ext     ::Bool    = true,
                 states_only::Bool    = false,
                 start      ::Symbol  = :crown) where {N}

Simulate tree according to the geographic **sse** model. The number of areas 
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
- `β::Array{Float64,1}`: per-area effect of covariates.
- `t::Float64`: simulation time.
- `x::Array{Float64,1}`: times where the covariate `y` is sampled.
- `y::Array{Float64,N}`: value of covariates, i.e., `f(x)`. Can be multivariate.
- `cov_mod::Tuple{Vararg{String}}`: specifies which rates are affected by 
covariates: `s` for speciation, `e` for extinction, and `g` for colonization.
More than 1 is possible. 
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
                      β          ::Array{Float64,1},
                      t          ::Float64,
                      x          ::Array{Float64,1},
                      y          ::Array{Float64,N},
                      cov_mod    ::Tuple{Vararg{String}};
                      δt         ::Float64 = 1e-4,
                      ast        ::Int64   = 0,
                      nspp_max   ::Int64   = 100_000,
                      retry_ext  ::Bool   = true,
                      rejectel0  ::Bool   = false,
                      verbose    ::Bool   = true,
                      rm_ext     ::Bool   = true,
                      states_only::Bool   = false,
                      start      ::Symbol = :crown) where {N}

  # make simulation
  ed, el, st, ea, ee, n, S, k = 
    simulate_edges(λ, μ, l, g, q, β, x, y, ast, t, δt, cov_mod, nspp_max, start)

  ne = lastindex(ee)

  maxsp = false

  if verbose
    @info "Tree with $n extant and $ne extinct species successfully simulated"

    if n < 2
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

  if (n + ne) > (nspp_max - (start == :crown ? 0 : 1))
    if verbose
      @warn string("Simulation surpassed the maximum of lineages allowed : ", nspp_max)
    end
    maxsp = true
  end

  if retry_ext 
    while n < 2 || (rejectel0 && in(0.0, el))

      ed, el, st, ea, ee, n, S, k = 
        simulate_edges(λ, μ, l, g, q, β, x, y, ast, t, δt, cov_mod, nspp_max, start)

      ne = lastindex(ee)

      if verbose
        @info "Tree with $n extant and $ne extinct species successfully simulated"

        if n < 2
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
    if n < 2
      return Dict{Int64, Vector{Float64}}(), ed, el, false
    end
  end

  # tip numbers
  tN = ed[ea,2]

  # tip states
  tS = st[ea]

  # tip states
  if states_only
    tv = tip_dictionary(tS)
    tv = states_to_values(tv, S, k)
    return tv, ed, el, maxsp
  end

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
                   β       ::Array{Float64,1},
                   x       ::Array{Float64,1},
                   y       ::Array{Float64,N},
                   si      ::Int64,
                   simt    ::Float64,
                   δt      ::Float64,
                   cov_mod ::Tuple{Vararg{String}},
                   nspp_max::Int64,
                   start   ::Symbol) where {N}

Simulate edges tree according to EGeoHiSSE.
"""
function simulate_edges(λ       ::Array{Float64,1},
                        μ       ::Array{Float64,1},
                        l       ::Array{Float64,1},
                        g       ::Array{Float64,1},
                        q       ::Array{Float64,1},
                        β       ::Array{Float64,1},
                        x       ::Array{Float64,1},
                        y       ::Array{Float64,N},
                        si      ::Int64,
                        simt    ::Float64,
                        δt      ::Float64,
                        cov_mod ::Tuple{Vararg{String}},
                        nspp_max::Int64,
                        start   ::Symbol) where {N}

  h  = div(isqrt(length(q)*4 + 1) + 1, 2)
  k  = div(length(l), h)
  ny = size(y, 2)

  # if multidimensional
  md = !isone(N)

  # check number of parameters and data
  model = id_mod(cov_mod, k, h, ny, λ, μ, l, g, q, β)

  # make approximate time function
  af! = make_af(x, y, Val(size(y,2)))

  # preallocate af! result vector 
  r = Array{Float64,1}(undef, size(y,2))

  # get function estimates at `simt`
  af!(simt, r)

  # create states
  S = create_states(k, h)

  ns = length(S)

  # areas & hidden states
  as = 1:k
  hs = 0:(h-1)

  # if initial state not assigned
  if iszero(si)

    # vector of initial state probabilities
    isp = Array{Float64,1}(undef, length(S))

    # assign initial state probabilities
    init_states_pr!(isp, λ, μ, l, g, q, β, S, r, k, h, ny, model, md, as, hs)

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
    event_probs(λ, μ, l, g, q, β, r, S, k, h, ny, model, md, δt, ns, as, hs)

  # make λ, μ, g & q probability functions
  updλpr! = make_updλpr!(λ, β, λpr, Sλpr, δt, k, model, md, S, ny)
  updμpr! = make_updμpr!(μ, l, β, μpr, Sμpr, δt, k, h, model, md, S, ny)
  updgpr! = make_updgpr!(g, β, gpr, Sgpr, δt, k, h, model, md, as, S, ny)
  updqpr! = make_updqpr!(Sqpr)

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
  gtos, μtos, qtos, λtos = makecorresschg(gpr, μpr, qpr, λpr, S, as, hs, k, ns)

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

      # estimate environment
      af!(simt, r)

      # update probabilities
      updλpr!(r)
      updμpr!(r)
      updgpr!(r)

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
    id_mod(cov_mod::Tuple{Vararg{String}}, 
           k      ::Int64, 
           h      ::Int64, 
           ny     ::Int64,
           λ      ::Array{Float64,1}, 
           μ      ::Array{Float64,1}, 
           l      ::Array{Float64,1},
           g      ::Array{Float64,1},
           q      ::Array{Float64,1},
           β      ::Array{Float64,1})

Check if number of parameters and `z(t)` functions 
is consistent with specified model.
"""
function id_mod(cov_mod::Tuple{Vararg{String}}, 
                k      ::Int64, 
                h      ::Int64, 
                ny     ::Int64,
                λ      ::Array{Float64,1}, 
                μ      ::Array{Float64,1}, 
                l      ::Array{Float64,1},
                g      ::Array{Float64,1},
                q      ::Array{Float64,1},
                β      ::Array{Float64,1})

  all_ok = true

  # check consistent length of parameter vectors proposed
   if (lastindex(λ) != (k+1)*h   || 
       lastindex(μ) != k*h       ||
       lastindex(l) != k*h       ||
       lastindex(g) != k*(k-1)*h ||
       lastindex(q) != h*(h-1))
    @info string(k, " areas with ", h, " hidden states should have the following parameter lengths:
       λ = ", ((k+1)*h), " parameters
       μ = ", (k*h), " parameters
       l = ", (k*h), " parameters
       g = ", (k*(k-1)*h), " parameter matrix
       q = ", (h*(h-1)), " parameters")
    @error "Parameter vector lengths not consistent with the number of states"
  end

  model = [false, false, false]

  for m in cov_mod
    # if speciation
    if occursin(r"^[s|S][A-za-z]*", m) 
      model[1] = true
    end
    # if extinction
    if occursin(r"^[e|E][A-za-z]*", m) 
      model[2] = true
    end
    # if dispersal
    if occursin(r"^[t|T|r|R|q|Q][A-za-z]*", m)
      model[3] = true
    end
  end

  # beta length should be
  yppar = ny == 1 ? 1 : div(ny,
         model[1]*k +
         model[2]*k +
         model[3]*k*(k-1))

  betaE = model[1]*k*h*yppar + model[2]*k*h*yppar + model[3]*k*(k-1)*h*yppar

  if lastindex(β) == betaE
    @info "Number of parameters are consistent with specified model"
  else
    all_ok = false
    @error string("this covariate model with ", k," areas, ",ny," covariates and ", h, " hidden states should have a β vector length of ", betaE)
  end

  mexp = "$(model[1] ? "speciation," : "")$(model[2] ? "extinction," : "")$(model[3] ? "transition," : "")"
  mexp = replace(mexp, "," => ", ")
  mexp = mexp[1:(end-2)]

  if all_ok 
    @info string("Simulating ", mexp, " covariate SSE model with ", k," states, ",h, " hidden states and ", ny, " covariates") 
  else
    @error "Parameter and z(t) function number not consistent with specified model"
  end

  return tuple(model...)::Tuple{Bool, Bool, Bool}
end





"""
    init_states_pr!(isp  ::Array{Float64,1},
                    λ    ::Array{Float64,1},
                    μ    ::Array{Float64,1},
                    l    ::Array{Float64,1},
                    g    ::Array{Float64,1},
                    q    ::Array{Float64,1},
                    β    ::Array{Float64,1},
                    S    ::Array{Sgh,1},
                    r    ::Array{Float64,1},
                    k    ::Int64,
                    h    ::Int64,
                    ny   ::Int64,
                    model::NTuple{3,Bool},
                    md   ::Bool,
                    as   ::UnitRange{Int64},
                    hs   ::UnitRange{Int64})

Estimate starting transition probabilities for initial states.
"""
function init_states_pr!(isp  ::Array{Float64,1},
                         λ    ::Array{Float64,1},
                         μ    ::Array{Float64,1},
                         l    ::Array{Float64,1},
                         g    ::Array{Float64,1},
                         q    ::Array{Float64,1},
                         β    ::Array{Float64,1},
                         S    ::Array{Sgh,1},
                         r    ::Array{Float64,1},
                         k    ::Int64,
                         h    ::Int64,
                         ny   ::Int64,
                         model::NTuple{3,Bool},
                         md   ::Bool,
                         as   ::UnitRange{Int64},
                         hs   ::UnitRange{Int64})

  # expected number of covariates
  ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

  # y per parameter
  yppar = isone(ny) ? 1 : div(ny,ncov)

  # starting values for models
  m2s = model[1]*k*h
  m3s = m2s + model[2]*k*h
  y2s = model[1]*yppar*k
  y3s = y2s + model[2]*yppar*k

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
    if model[3]
      gr = 0.0
      if length(s.g) > 1
        for ta = s.g, fa = setdiff(s.g, ta)
          gr += expf(g[s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)], 
                     β[m3s + s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)], 
                     md ? r[y3s + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)] : r[1])
        end
      end
    else
      gr = 0.0
      if length(s.g) > 1
        for ta = s.g, fa = setdiff(s.g, ta)
          gr += g[s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)]
        end
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

  # expected number of covariates
  ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

  # y per parameter
  yppar = isone(ny) ? 1 : div(ny,ncov)

  # starting values for models
  m2s = model[1]*k
  m3s = m2s + model[2]*k
  y2s = model[1]*yppar*k
  y3s = y2s + model[2]*yppar*k

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
      # within-area speciation
      if model[1]
        for a in s.g
          na += 1
          λpr[si][na] = expf(λ[(k+1)*s.h + a], β[s.h*ncov + a], md ? r[a] : r[1])*δt
        end
      else
        for a in s.g
          na += 1
          λpr[si][na] = λ[(k+1)*s.h + a]*δt
        end
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
        if model[2]
          for a in s.g
            μpr[si][1] = expf(μ[k*s.h + a], β[m2s + s.h*ncov + a], md ? r[y2s + a] : r[1])*δt
          end
        else
          for a in s.g
            μpr[si][1] = μ[k*s.h + a]*δt
          end
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
      if model[3]
        for (i,ta) = enumerate(setdiff(as,s.g)), fa = s.g
          gpr[si][i] += 
            expf(q[s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)], 
                 β[m3s + s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)], 
                 md ? r[y3s + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)] : r[1])*δt
        end
      else
        for (i,ta) = enumerate(setdiff(as,s.g)), fa = s.g
          gpr[si][i] += 
            g[s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)]*δt
        end
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
    makecorresschg(gpr::Array{Array{Float64,1},1},
                   μpr::Array{Array{Float64,1},1},
                   qpr::Array{Array{Float64,1},1},
                   λpr::Array{Array{Float64,1},1},
                   S  ::Array{Sgh, 1},
                   as::UnitRange{Int64},
                   hs::UnitRange{Int64},
                   k ::Int64,
                   ns::Int64)

Make vector of vectors of corresponding states changes 
for **gains** and **looses** and **hidden states**.
"""
function makecorresschg(gpr::Array{Array{Float64,1},1},
                        μpr::Array{Array{Float64,1},1},
                        qpr::Array{Array{Float64,1},1},
                        λpr::Array{Array{Float64,1},1},
                        S  ::Array{Sgh, 1},
                        as ::UnitRange{Int64},
                        hs ::UnitRange{Int64},
                        k  ::Int64,
                        ns ::Int64)
  # for *gains*
  gtos = Array{Int64,1}[]
  for i in Base.OneTo(ns)
    push!(gtos, zeros(Int64,length(gpr[i])))
  end
  for si in Base.OneTo(ns)
    s = S[si]
    for (i,ta) = enumerate(setdiff(as,s.g))
      gtos[si][i] = findfirst(r -> isSghequal(r,Sgh(union(s.g, ta),s.h)), S)
    end
  end

  # # for *losses*
  μtos = Array{Int64,1}[]
  for i in Base.OneTo(ns)
    push!(μtos, zeros(Int64,length(μpr[i])))
  end
  for si in Base.OneTo(ns)
    s = S[si]
    isone(length(s.g)) && continue
    na = 0
    for a in s.g
      na += 1
      μtos[si][na] = findfirst(r -> isSghequal(r,Sgh(setdiff(s.g, a),s.h)), S)
    end
  end

  # # for *hidden states*
  qtos = Array{Int64,1}[]
  for i in Base.OneTo(ns)
    push!(qtos, zeros(Int64,length(qpr[i])))
  end
  for si in Base.OneTo(ns)
    s = S[si]
     for (i, th) = enumerate(setdiff(hs, s.h))
      qtos[si][i] = th*(2^k-1) + si - s.h*(2^k-1)
    end
  end

  # # for *speciation*
  λtos = Array{Array{Int64,1},1}[]
  for i in Base.OneTo(ns)
    push!(λtos, fill([0,0],length(λpr[i])))
  end
  for si in Base.OneTo(ns)
    s = S[si]
    na = 0
    if isone(length(s.g))
      for a in s.g
        na += 1
        λtos[si][na] = [si,si]
      end
    else
      # within-area speciation
      for a in s.g
        na += 1
        λtos[si][na] = [a + (2^k-1)*s.h, si]
      end
      # between-area speciation
      vs = vicsubsets(s.g)[1:div(end,2)]
      for v in vs
        na += 1
        λtos[si][na] = 
          Int64[findfirst(r -> isSghequal(r,Sgh(v[1],s.h)), S),
           findfirst(r -> isSghequal(r,Sgh(v[2],s.h)), S)]
      end
    end
  end

  return gtos, μtos, qtos, λtos
end





"""
    make_updλpr!(λ    ::Array{Float64,1},
                 β    ::Array{Float64,1},
                 λpr  ::Array{Array{Float64,1},1},
                 Sλpr ::Array{Float64,1},
                 δt   ::Float64,
                 k    ::Int64,
                 model::NTuple{3,Bool},
                 md   ::Bool,
                 ny   ::Int64)

Make λ instantaneous probability function for each state.
"""
function make_updλpr!(λ    ::Array{Float64,1},
                      β    ::Array{Float64,1},
                      λpr  ::Array{Array{Float64,1},1},
                      Sλpr ::Array{Float64,1},
                      δt   ::Float64,
                      k    ::Int64,
                      model::NTuple{3,Bool},
                      md   ::Bool,
                      S  ::Array{Sgh, 1},
                      ny   ::Int64)

  # expected number of covariates
  ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

  # y per parameter
  yppar = isone(ny) ? 1 : div(ny,ncov)

  # starting values for models
  m2s = model[1]*k
  m3s = m2s + model[2]*k
  y2s = model[1]*yppar*k
  y3s = y2s + model[2]*yppar*k

  lS = length(S)

  # if dependent on `z(t)` and multidimensional
  function f1(r::Array{Float64,1})

    @inbounds begin
      for si in Base.OneTo(lS)

        s = S[si]

        na = 0
        # within-area speciation
        for a in s.g
          na += 1
          λpr[si][na] = expf(λ[(k+1)*s.h + a], β[s.h*ncov + a], md ? r[a] : r[1])*δt
        end
        # between-area speciation
        if na > 1
          λpr[si][k+1:end] .= λ[k+1 + (k+1)*s.h]*δt
        end

        Sλpr[si] = sum(λpr[si]) 

      end
    end
  end

  # if dependent on `z(t)` and unidimensional
  function f2(r::Array{Float64,1})
  end

  if model[1]
    return f1
  else
    return f2
  end
end





"""
    make_updμpr!(μ    ::Array{Float64,1},
                 l    ::Array{Float64,1},
                 β    ::Array{Float64,1},
                 μpr  ::Array{Array{Float64,1},1},
                 Sμpr ::Array{Float64,1},
                 δt   ::Float64,
                 k    ::Int64,
                 h    ::Int64,
                 model::NTuple{3,Bool},
                 md   ::Bool,
                 ny   ::Int64)

Make `μ` instantaneous probability function
"""
function make_updμpr!(μ    ::Array{Float64,1},
                      l    ::Array{Float64,1},
                      β    ::Array{Float64,1},
                      μpr  ::Array{Array{Float64,1},1},
                      Sμpr ::Array{Float64,1},
                      δt   ::Float64,
                      k    ::Int64,
                      h    ::Int64,
                      model::NTuple{3,Bool},
                      md   ::Bool,
                      S  ::Array{Sgh, 1},
                      ny   ::Int64)
  
  # expected number of covariates
  ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

  # y per parameter
  yppar = isone(ny) ? 1 : div(ny,ncov)

  # starting values for models
  m2s = model[1]*k
  m3s = m2s + model[2]*k
  y2s = model[1]*yppar*k
  y3s = y2s + model[2]*yppar*k

  lS = length(S)

  function f1(r ::Array{Float64,1})

    @inbounds begin
      for si in Base.OneTo(lS)

        s = S[si]

        if isone(length(s.g))
          for a in s.g
            μpr[si][1] = 
              expf(μ[k*s.h + a], β[m2s + s.h*ncov + a], md ? r[y2s + a] : r[1])*δt
          end
        else
          for a in s.g
            μpr[si][a] =
              expf(μ[k*s.h + a], β[m2s + s.h*ncov + a], md ? r[y2s + a] : r[1])*δt
          end
        end

        Sμpr[si] = sum(μpr[si]) 
      end
    end
  end

  # if dependent on `z(t)` and unidimensional
  function f2(r::Array{Float64,1})
  end

  if model[2]
    return f1
  else
    return f2
  end
end





"""
    make_updgpr!(g    ::Array{Float64,1},
                 β    ::Array{Float64,1},
                 gpr  ::Array{Array{Float64,1},1},
                 Sgpr ::Array{Float64,1},
                 δt   ::Float64,
                 k    ::Int64,
                 h    ::Int64,
                 model::NTuple{3,Bool},
                 md   ::Bool,
                 as   ::UnitRange{Int64},
                 S    ::Array{Sgh,1},
                 ny   ::Int64)

Make gain instantaneous probability function
"""
function make_updgpr!(g    ::Array{Float64,1},
                      β    ::Array{Float64,1},
                      gpr  ::Array{Array{Float64,1},1},
                      Sgpr ::Array{Float64,1},
                      δt   ::Float64,
                      k    ::Int64,
                      h    ::Int64,
                      model::NTuple{3,Bool},
                      md   ::Bool,
                      as   ::UnitRange{Int64},
                      S    ::Array{Sgh,1},
                      ny   ::Int64)

  # expected number of covariates
  ncov = model[1]*k + model[2]*k + model[3]*k*(k-1)

  # y per parameter
  yppar = isone(ny) ? 1 : div(ny,ncov)

  # starting values for models
  m3s = model[1]*k + model[2]*k
  y3s = model[1]*yppar*k + model[2]*yppar*k

  # do setdiff before
  sdf = Array{Int64,1}[]
  for s = S
    push!(sdf, setdiff(as,s.g))
  end

  lS = length(S)

  function f1(r::Array{Float64,1})

    @inbounds begin
      for si in Base.OneTo(lS)

        s = S[si]

        for (i,ta) = enumerate(sdf[si])
          gpr[si][i] = 0.0
          for fa = s.g
            gpr[si][i] += 
              expf(g[s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)], 
                   β[m3s + s.h*k*(k-1) + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)], 
                   md ? r[y3s + (fa-1)*(k-1) + (ta > fa ? ta - 1 : ta)] : r[1])*δt
          end
        end

        Sgpr[si] = sum(gpr[si]) 

      end
    end
  end

  function f2(r::Array{Float64,1})
  end

  if model[3]
    return f1
  else
    return f2
  end
end





"""
    make_updqpr!(Sqpr::Array{Float64,1})

Make hidden states instantaneous probability function
"""
function make_updqpr!(Sqpr::Array{Float64,1})

  function f(si ::Int64)
    return Sqpr[si]
  end

  return f
end





"""
    expf(α::Float64, β::Float64, x::Float64)

Exponential regression function for rates with base rate `α`, 
coefficient `β` and covariate `x`.
"""
expf(α::Float64, β::Float64, x::Float64) = α * exp(β * x)





"""
    prop_sample(s::Array{Float64,1}, prv::Array{Float64,1}, ns::Int64)

Sample one state with the respective (relative) probabilities 
given `prv`.
"""
function prop_sample(s::Array{Float64,1}, prv::Array{Float64,1}, ns::Int64)

  cumsum!(s, prv)
  @simd for i in Base.OneTo(ns)
    s[i] /= s[ns]  
  end

  return searchsortedfirst(s, rand())
end



