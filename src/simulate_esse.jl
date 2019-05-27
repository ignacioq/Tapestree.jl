#=
Simulation for ESSE

Ignacio Quintero Mächler

t(-_-t)

January 12 2017
=#





"""
    simesse(λ       ::Array{Float64,1}, 
            μ       ::Array{Float64,1},
            Q       ::Array{Float64,2},
            x       ::Array{Float64,1},
            y       ::Array{Float64,1},
            esse_mod::String;
            nspp    ::Int64   = 20,
            δt      ::Float64 = 1e-3)

Simulate tree according to ESSE.
"""
function simesse(λ       ::Array{Float64,1}, 
                 μ       ::Array{Float64,1},
                 Q       ::Array{Float64,2},
                 β       ::Array{Float64},
                 x       ::Array{Float64,1},
                 y       ::Array{Float64};
                 esse_mod::String  = "speciation",
                 nspp    ::Int64   = 20,
                 δt      ::Float64 = 1e-3)

  # make simulation
  ed, el, st, simt, af, r = simedges(λ, μ, Q, β, x, y, nspp, δt, esse_mod)

  # organize in postorder
  ed     = numberedges(ed, nspp)
  ed, el = postorderedges(ed, el, nspp)

  tip_val = Dict(i => st[i] for i = 1:nspp)

  ##
  # organize new x and y for inference
  th = tree_height(el, ed, nspp)

  # is y multidimensional
  md = size(y,2) > 1

  # if simulation time is equal to simulated tree height
  if isapprox(th, simt, atol = 1e-9)

    idx = searchsortedlast(x, simt)

    nx = x[1:idx]

    push!(nx, th)

    if md
      ny = y[1:idx, :]
      af(th, r)
      ny = vcat(ny, reshape(r,(1,:)))
    else
      ny = y[1:idx]
      push!(ny, af(th))
    end

    ny = reverse(ny, dims = 1)

  # if simulation time is not equal to simulated tree height
  else
    idxs = searchsortedlast(x, simt)
    idxt = searchsortedlast(x, th)

    # new start
    nzero = simt - th

    idxz = searchsortedfirst(x, nzero)

    nx = x[idxz:idxs]
    ny = y[idxz:idxs, :]

    # push to time
    pushfirst!(nx, nzero)
    push!(nx, simt)

    # make x start at 0
    nx .-= nzero

    if md
      ny = y[idxz:idxs, :]
      af(nzero, r)
      ny = vcat(ny, reshape(r,(1,:)))
    else
      ny = y[idxz:idxs]
      push!(ny, af(nzero))
    end

    if md
      af(th, r)
      ny = vcat(ny, reshape(r,(1,:)))
    else
      push!(ny, af(th))
    end

    ny = reverse(ny, dims = 1)
  end

  tip_val = states_to_values(tip_val, lastindex(λ))

  return tip_val, ed, el, nx, ny
end





"""
    simedges(λ       ::Array{Float64,1}, 
             μ       ::Array{Float64,1},
             Q       ::Array{Float64,2},
             β       ::Array{Float64,1},
             af      ::Function,
             nspp    ::Int64,
             δt      ::Float64,
             esse_mod::String)

Simulate edges tree according to ESSE.
"""
function simedges(λ       ::Array{Float64,1}, 
                  μ       ::Array{Float64,1},
                  Q       ::Array{Float64,2},
                  β       ::Array{Float64,1},
                  x       ::Array{Float64,1},
                  y       ::Array{Float64},
                  nspp    ::Int64,
                  δt      ::Float64,
                  esse_mod::String)

  k = length(λ)

  # if multidimensional
  md = size(y, 2) > 1

  # check number of parameters and data
  ws, we, wq = id_mod(esse_mod, k, md, size(y, 2), λ, μ, β, Q)

  # total simulation time
  simt = 0.0

  # make approximate time function
  af = make_approxf(x, y)

  if md 
    if wq
      r = Array{Float64,1}(undef, k*k - k)
    else
      r = Array{Float64,1}(undef, k)
    end
  end

  # get function estimates at `0.0`
  if md
    af(0.0, r)
  else
    r = af(0.0)
  end

  # preallocate vector of transitions probabilities after a transition has
  # been selected
  qpr = Array{Float64,1}(undef,k)
  spr = Array{Float64,1}(undef,k)

  # estimate starting q probabilities for initial states
  if wq
    # non-diagonal Q indices
    Qndi = setdiff(1:(k*k), 1:(k+1):(k*k))

    # Q column indices without diagonals
    Qci = Array{Int64,1}[]
    for i in Base.OneTo(k)
      push!(Qci, setdiff(1+k*(i-1): k+k*(i-1), 1:(k+1):(k*k)))
    end

    # β indices to match the ones for Q
    βi = UnitRange{Int64}[]
    for j in Base.OneTo(k)
      push!(βi, (j-1)k + (2-j):(j*k-j))
    end

    if md
      qsum = sum(l -> expf(Q[Qndi[l]], β[l], r[l]), 1:lastindex(Qndi))
      for j in Base.OneTo(k)
        qpr[j] = 
          sum(i -> expf(Q[Qci[j][i]], β[βi[j][i]], r[βi[j][i]]), 1:(k-1))/qsum
      end
    else 
      qsum = sum(l -> expf(Q[Qndi[l]], β[l], r), 1:lastindex(Qndi))
      for j in Base.OneTo(k)
        qpr[j] = 
          sum(i -> expf(Q[Qci[j][i]], β[βi[j][i]], r), 1:(k-1))/qsum
      end
    end
  else
    qsum = sum(Q)
    sum!(qpr, transpose(Q))./qsum
  end

  # sample initial state
  si = prop_sample(spr, qpr, k)

  # n = 2 starting lineages (assume no stem branch)
  n = 2

  # edges alive
  ea = [1, 2]

  # state vector
  st = fill(si, 2)

  # edge array
  ed = zeros(Int64, 2nspp, 2)
  ed[ea,:] = [1 2;
              1 3]

  # edge lengths
  el = zeros(2nspp)

  # make λ, μ & Q probability functions
  λpr = make_λpr(λ, β, δt, md, ws)
  μpr = make_μpr(μ, β, δt, md, we)
  Qpr = make_Qpr(Q, β, δt, md, wq)

  # make function for estimating transition probabilites
  rqpr! = make_rqpr(Q, β, δt, md, wq)

  # start simulation
  while true

    # keep track of time
    simt += δt

    # estimate `z(simt)` or `z_i(simt)`
    if md
      af(simt, r)
    else
      r = af(simt)
    end

    # one time step
    for i in Base.OneTo(n)

      # update edge length
      el[ea[i]] += δt

      #=
          Q event?
      =#
      if rand() < Qpr(st[i], r)
        # get `out` transition probabilities
        rqpr!(st[i], qpr, r)
        # sample from the transition probabilities
        st[i] = prop_sample(spr, qpr, k)
      end

      #=
          μ event?
      =#
      if rand() < μpr(st[i], r)

        # update time in other extant lineages
        el[ea[(i+1):n]] .+= δt

        # node to remove
        nod = ed[ea[i],1]

        # top or bottom edge of node?
        top = ea[i] == findfirst(isequal(nod), ed[:,1])

        # remove edge (only preserving persisting species)
        ned = findfirst(isequal(0), ed)[1]
        
        if ned == 3
          printstyled("What would you do if an endangenered animal is eating an endangered plant? Sometimes nature is too harsh... \n", 
            color=:red)
          error("tree went extinct... rerun simulation")

        end

        @views ed[ea[i]:(ned-1),:] = ed[(ea[i]+1):ned,:]

        # remove edge length
        @views el[ea[i]:(ned-1)] = el[(ea[i]+1):ned]

        # update alive lineages
        deleteat!(ea,i)
        ea[i:end] .-= 1

        # remove node and extra edge
        @views pr = findfirst(isequal(nod), ed[:,1])
        @views da = findfirst(isequal(nod), ed[:,2])

        if da == nothing
          da = 0
        end

        if da != 0
          ed[da,2] = ed[pr,2]
          el[da]  += el[pr]
        end

        # remove intervening edge
        @views ed[pr:(ned-2),:] = ed[(pr+1):(ned-1),:]

        # remove edge lengths of intervening node
        @views el[pr:(ned-2)] = el[(pr+1):(ned-1)]

        ## update alive lineages
        # is the remaining node terminal?
        if da == 0 || in(ed[da,2], ed[:,1])
          ea[i:end] .-= 1
        else 
          ea[findfirst(isequal(pr),ea)] = da
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

        # no more events at this time
        break
      end

      #=
          λ event?
      =#
      if rand() < λpr(st[i], r)

        ### add new edges
        # start node
        ed[ea[end] + 1,1] = ed[ea[i],2]
        ed[ea[end] + 2,1] = ed[ea[i],2]

        # end nodes
        ed[ea[end] + 1,2] = maximum(ed)+1
        ed[ea[end] + 2,2] = maximum(ed)+1

        # update time in other extant lineages
        el[ea[(i+1):n]] .+= δt

        # update living edges
        push!(ea, ea[end]+1, ea[end]+2)
        deleteat!(ea, i)

        # update living states
        push!(st,st[i],st[i])
        deleteat!(st,i)

        # update number of species alive
        n = lastindex(ea)

        # no more events at this time
        break
      end

    end

    n == nspp + 1 && break

  end

  ## remove last species and organize
  # identify daughter node and delete last two edges
  @views da = findfirst(isequal(ed[end,1]), ed[:,2])

  ed = ed[1:(2nspp-2),:]
  el = el[1:(2nspp-2)]

  # organize extant species
  pop!(ea)
  pop!(ea)
  push!(ea,da)

  sort!(ea)

  # organize extant states
  stf = st[end]
  pop!(st)
  pop!(st)

  insert!(st,findfirst(isequal(da), ea),stf)

  return ed, el, st, simt, af, r
end





"""
    id_mod(esse_mod::String, k::Int64,  md::Bool, nzt::Int64)

Check if number of parameters and `z(t)` functions 
is consistent with specified model.
"""
function id_mod(esse_mod::String, 
                k  ::Int64, 
                md ::Bool, 
                nzt::Int64,
                λ  ::Array{Float64,1}, 
                μ  ::Array{Float64,1}, 
                β  ::Array{Float64}, 
                Q  ::Array{Float64})

  allok = true

  # if speciation, extinction or transition model
  if occursin(r"^[s|S][A-za-z]*", esse_mod)

    mdd = "speciation" 
    ws, we, wq = true, false, false

    if lastindex(λ) == lastindex(μ) == lastindex(β) == size(Q,1) == size(Q,2)
      printstyled("Number of parameters are consistent with specified model \n", 
        color=:green)
    else
      allok = false
      printstyled("Parameter number for speciation ESSE with $k states:
                   λ = $(k) parameters
                   μ = $(k) parameters
                   Q = $k x $k parameter matrix
                   β = $(k) parameters", 
        color=:red)
    end

  elseif occursin(r"^[e|E][A-za-z]*", esse_mod)

    if lastindex(λ) == lastindex(μ) == lastindex(β) == size(Q,1) == size(Q,2)
      printstyled("Number of parameters are consistent with specified model \n", 
        color=:green)    
    else
      allok = false
      printstyled("Parameter number for extinction ESSE with $k states:
                   λ = $(k) parameters
                   μ = $(k) parameters
                   Q = $k x $k parameter matrix
                   β = $(k) parameters", 
        color=:red)
      error("Parameter number not consistent with specified model")
    end

    mdd = "extinction" 
    ws, we, wq = false, true, false

  elseif occursin(r"^[t|T|r|R|q|Q][A-za-z]*", esse_mod)

    if lastindex(λ) == lastindex(μ) == size(Q,1) == size(Q,2) && 
       lastindex(β) == (k*k - k)
      printstyled("Number of parameters are consistent with specified model \n", 
        color=:green)    
    else
      allok = false
      printstyled("Parameter number for extinction ESSE with $k states:
                   λ = $(k) parameters
                   μ = $(k) parameters
                   Q = $k x $k parameter matrix
                   β = $(k*k - k) parameters", 
        color=:red)
    end
    mdd = "transition" 
    ws, we, wq = false, false, true
  else
    ws, we, wq = false, false, false
    warning("esse_mod does not match any ESSE model, running MuSSE \n")
  end

  if md 
    if (ws || we) && nzt != k
      allok = false
      printstyled("For character-specific $mdd ESSE model for $k states:
               z(t) = $(k) functions",
        color=:red)
    end
    if wq && nzt != (k*k - k)
      allok = false
      printstyled("For character-specific $mdd ESSE model for $k states:
               z(t) = $(k*k - k) functions",
        color=:red)
    end
  end

  if allok 
    printstyled("Simulating $mdd ESSE model with $k states \n", 
      color=:green)
  else
    error("Parameter and z(t) function number not consistent with specified model")
  end

  return ws, we, wq
end




"""
    make_λpr(λ ::Array{Float64,1},
             β ::Array{Float64}, 
             δt::Float64,
             md::Bool,
             ws::Bool)

Make λ instantaneous probability function
"""
function make_λpr(λ ::Array{Float64,1},
                  β ::Array{Float64,1}, 
                  δt::Float64,
                  md::Bool,
                  ws::Bool)

  # if dependent on `z(t)` and multidimensional
  function f1(st ::Int64,
              aft::Array{Float64,1})
    return expf(λ[st], β[st], aft[st])*δt
  end

  # if dependent on `z(t)` and unidimensional
  function f2(st ::Int64,
              aft::Float64)
    return expf(λ[st], β[st], aft)*δt
  end

  # if **not** dependent on `z(t)` and multidimensional
  function f3(st ::Int64,
              aft::Array{Float64,1})
    return λ[st]*δt
  end

  # if **not** dependent on `z(t)` and unidimensional
  function f4(st ::Int64,
              aft::Float64)
    return λ[st]*δt
  end

  if ws
    if md
      return f1
    else
      return f2
    end
  else
    if md
      return f3
    else
      return f4
    end
  end
end





"""
    make_μpr(μ ::Array{Float64,1},
             β ::Array{Float64}, 
             δt::Float64,
             md::Bool,
             we::Bool)

Make `μ` instantaneous probability function
"""
function make_μpr(μ ::Array{Float64,1},
                  β ::Array{Float64,1}, 
                  δt::Float64,
                  md::Bool,
                  we::Bool)

  # if dependent on `z(t)` and multidimensional
  function f1(st ::Int64,
              aft::Array{Float64,1})
    return expf(μ[st], β[st], aft[st])*δt
  end

  # if dependent on `z(t)` and unidimensional
  function f2(st ::Int64,
              aft::Float64)
    return expf(μ[st], β[st], aft)*δt
  end

  # if **not** dependent on `z(t)` and multidimensional
  function f3(st ::Int64,
              aft::Array{Float64,1})
    return μ[st]*δt
  end

  # if **not** dependent on `z(t)` and unidimensional
  function f4(st ::Int64,
              aft::Float64)
    return μ[st]*δt
  end

  if we
    if md
      return f1
    else
      return f2
    end
  else
    if md
      return f3
    else
      return f4
    end
  end
end





"""
    make_Qpr(Q ::Array{Float64,2},
             β ::Array{Float64}, 
             δt::Float64,
             md::Bool,
             wq::Bool)

Make `Q` instantaneous probability function
"""
function make_Qpr(Q ::Array{Float64,2},
                  β ::Array{Float64,1}, 
                  δt::Float64,
                  md::Bool,
                  wq::Bool)
  k = size(Q,1)

  # Q row indices without diagonals
  Qri = Array{Int64,1}[]
  for i in Base.OneTo(k)
    push!(Qri, setdiff(i:k:(k*k), 1:(k+1):(k*k)))
  end

  # β indices to match the ones for Q
  βi = Array{Int64,1}[]
  for j in Base.OneTo(k)
    tv = Int64[]
    for i in Base.OneTo(k-1)
      push!(tv, Qri[j][i] - i)
    end
    push!(βi, tv)
  end

  ## Functions
  # if dependent on `z(t)` and multidimensional
  function f1(st ::Int64,
              aft::Array{Float64,1})
    rsum = 0.0
    for i in Base.OneTo(k-1)
      rsum += expf(Q[Qri[st][i]], β[βi[st][i]], aft[βi[st][i]])
    end
    return rsum*δt
  end

  # if dependent on `z(t)` and unidimensional
  function f2(st ::Int64,
              aft::Float64)
    rsum = 0.0
    for i in Base.OneTo(k-1)
      rsum += expf(Q[Qri[st][i]], β[βi[st][i]], aft)
    end
    return rsum*δt
  end

  # if **not** dependent on `z(t)` and multidimensional
  function f3(st ::Int64,
              aft::Array{Float64,1})
    return @views sum(Q[st,:])*δt
  end

  # if **not** dependent on `z(t)` and unidimensional
  function f4(st ::Int64,
              aft::Float64)
    return @views sum(Q[st,:])*δt
  end

  if wq
    if md
      return f1
    else
      return f2
    end
  else
    if md
      return f3
    else
      return f4
    end
  end
end





"""
    make_rqpr(Q ::Array{Float64,2},
              β ::Array{Float64}, 
              δt::Float64,
              md::Bool,
              wq::Bool)

Make cumulative probability for instantaneous change in 
a given row in `Q` function.
"""
function make_rqpr(Q ::Array{Float64,2},
                   β ::Array{Float64}, 
                   δt::Float64,
                   md::Bool,
                   wq::Bool)
  k = size(Q,1)

  # Q row indices without diagonals
  Qri = Array{Int64,1}[]
  for i in Base.OneTo(k)
    push!(Qri, setdiff(i:k:(k*k), 1:(k+1):(k*k)))
  end

  # β indices to match the ones for Q
  βi = Array{Int64,1}[]
  for j in Base.OneTo(k)
    tv = Int64[]
    for i in Base.OneTo(k-1)
      push!(tv, Qri[j][i] - i)
    end
    push!(βi, tv)
  end

  ## Functions
  # if dependent on `z(t)` and multidimensional
  function f1(sti ::Int64,
              qpr::Array{Float64,1},
              aft::Array{Float64,1})
    l = 1
    for i in Base.OneTo(k)
      if i == sti
        qpr[i] = 0.0
      else
        qpr[i] = expf(Q[Qri[sti][l]], β[βi[sti][l]], aft[βi[sti][l]])
        l += 1
      end
    end

    return nothing
  end

  # if dependent on `z(t)` and unidimensional
  function f2(sti ::Int64,
              qpr::Array{Float64,1},
              aft::Float64)

    l = 1
    for i in Base.OneTo(k)
      if i == sti
        qpr[i] = 0.0
      else
        qpr[i] = expf(Q[Qri[sti][l]], β[βi[sti][l]], aft)
        l += 1
      end
    end

    return nothing
  end

  # if **not** dependent on `z(t)` and multidimensional
  function f3(sti ::Int64,
              qpr::Array{Float64,1},
              aft::Array{Float64,1})

    for i in Base.OneTo(k)
      qpr[i] = Q[sti,i]
    end

    return nothing
  end

  # if **not** dependent on `z(t)` and unidimensional
  function f4(sti ::Int64,
              qpr::Array{Float64,1},
              aft::Float64)
    for i in Base.OneTo(k)
      qpr[i] = Q[sti,i]
    end

    return nothing
  end

  if wq
    if md
      return f1
    else
      return f2
    end
  else
    if md
      return f3
    else
      return f4
    end
  end
end







