#=

insane tree plot

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#




"""
    b(tree::T) 
    d(tree::T) 
    lb(tree::T)
    ld(tree::T)
    t(tree::T) 
    nd(tree::T)
    dμ(tree::iT)

Predefined functions for plotting: 
  `b`  speciation rates
  `d`  extinction rates
  `lb` log speciation rates
  `lb` log speciation rates
  `t`  turnover
  `nd` net diversification
  `dμ` change in speciation rates
"""
b(tree::iT)  = exp.(lλ(tree))
d(tree::iT)  = exp.(lμ(tree))
lb(tree::iT) = lλ(tree)
ld(tree::iT) = lμ(tree)
t(tree::iT)  = exp.(lμ(tree)) ./ exp.(lλ(tree))
lt(tree::iT) = log.(exp.(lμ(tree)) ./ exp.(lλ(tree)))
nd(tree::iT) = exp.(lλ(tree)) .- exp.(lμ(tree))
function dλ(tree::iT)
  dd = diff(exp.(lλ(tree)))
  return append!(dd, dd[end])
end
function dμ(tree::iT)
  dd = diff(exp.(lμ(tree)))
  return append!(dd, dd[end])
end
function dλc(tree::iT)
  lv = lλ(tree)
  fill(exp(lv[end]) - exp(lv[1]), lastindex(lv))
end
function dμc(x)
  lv = lμ(x)
  fill(exp(lv[end]) - exp(lv[1]), lastindex(lv))
end



"""
    function f(tree::T;
               zf  ::Function,
               shownodes  = (T <: iTf),
               tip        = false,
               speciation = false,
               extinct    = false,
               fossil     = true,
               type       = :phylogram) where {T <: iT}

Recipe for plotting a Type `iT`.
"""
@recipe function f(tree::T,
                   zf  ::Function;
                   shownodes  = (T <: iTf),
                   tip        = false,
                   speciation = false,
                   extinct    = false,
                   fossil     = true,
                   type       = :phylogram,
                   simple     = false) where {T <: iT}

  x = Float64[]
  y = Float64[]
  z = Float64[]

  th  = treeheight(tree)
  nts = ntips(tree)

  if type === :lengthrates
    _rplottree_lr!(tree, 0.0, 1:nts, zf, x, y, z)
  else
    _rplottree!(tree,     th, 1:nts, zf, x, y, z, simple)
  end

  ntF = Float64(nts)

  # plot defaults
  legend              --> :none
  colorbar_fontfamily --> :Helvetica
  grid                --> :off
  fontfamily          --> :Helvetica

  if type === :radial
    x, y, z = append_forradial(x, y, z, 50)
    polar_coords!(x, y, 360.0/ntF, th)

    ylims           --> (-th*1.05, th*1.05)
    xlims           --> (-th*1.05, th*1.05)
    xticks          --> (nothing)
    xshowaxis       --> false
    colorbar        --> true
    yshowaxis       --> false
    yticks          --> nothing

  elseif type === :lengthrates
    xguide          --> "cumulative rates"
    colorbar        --> true
    yshowaxis       --> false
    yticks        --> nothing

  else
    xlims           --> (-th*0.05, th*1.05)
    xguide          --> "time"
    xflip           --> true
    tickfontfamily  --> :Helvetica
    tickfontsize    --> 8
    xtick_direction --> :out

    if type === :phylogram
      ylims         --> (1.0-(0.05*ntF), ntF+(0.05*ntF))
      colorbar      --> true
      yshowaxis     --> false
      yticks        --> nothing

    else
      colorbar  --> false
      yshowaxis --> true
    end
  end

  if shownodes

    xN = Float64[]
    yN = Float64[]

    _rplottree!(tree, th, 1:nts, xN, yN)

    shape = Symbol[:circle]
    col   = Symbol[:pink]
    alpha =
      Float64[(0.5+0.5*(!isdefined(tree, :fx) || isfix(tree))) *
              Float64(speciation)]
    _nodeproperties!(tree, shape, col, alpha,
      Float64(tip), Float64(speciation), Float64(extinct), Float64(fossil))

    @series begin
      markershape --> shape
      markercolor --> col
      markeralpha --> alpha
      markersize  --> 2.0
      if type === :phylogram
        xN, yN
      elseif type === :radial
        polar_coords!(xN, yN, 360.0/ntF, th)
        xN, yN
      end
    end
  end

  line_z     --> z
  linecolor  --> cgrad(:roma, rev = true)

  if type === :rates
    return x, z
  else
    return x, y
  end
end





"""
    _rplottree!(tree  ::T,
                xc    ::Float64,
                yr    ::UnitRange{Int64},
                zf    ::Function,
                x     ::Array{Float64,1},
                y     ::Array{Float64,1},
                z     ::Array{Float64,1},
                simple::Bool) where {T <: iT}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplottree!(tree  ::T,
                     xc    ::Float64,
                     yr    ::UnitRange{Int64},
                     zf    ::Function,
                     x     ::Array{Float64,1},
                     y     ::Array{Float64,1},
                     z     ::Array{Float64,1},
                     simple::Bool) where {T <: iT}

  # tree δt and nsδt
  δt = dt(tree)

  # add horizontal lines
  yc = Float64(yr[1] + yr[end])*0.5

  if simple
    zv = zf(tree)
    l  = lastindex(zv)
    push!(z, mean(zf(tree)))
    push!(y, yc)
    push!(x, xc)
  else
    # plot function
    zv = copy(zf(tree))
    l  = lastindex(zv)
    @simd for i in Base.OneTo(l-1)
      push!(x, xc - Float64(i-1)*δt)
      push!(y, yc)
      push!(z, zv[i])
    end
  end

  zc = last(zv)
  push!(x, xc - (Float64(l-2)*δt + fdt(tree)), NaN)
  push!(y, yc, NaN)
  push!(z, zc, NaN)

  xc -= e(tree)

  if def1(tree)
    if def2(tree)
      ntip1 = ntips(tree.d1)
      ntip2 = ntips(tree.d2)

      # add vertical lines
      push!(x, xc, xc, NaN)

      yr1 = yr[1:ntip1]
      yr2 = yr[(ntip1+1):(ntip1+ntip2)]
      push!(y, Float64(yr1[1] + yr1[end])*0.5,
               Float64(yr2[1] + yr2[end])*0.5,
               NaN)

      push!(z, zc, zc, NaN)

      _rplottree!(tree.d1, xc, yr1, zf, x, y, z, simple)
      _rplottree!(tree.d2, xc, yr2, zf, x, y, z, simple)
    else
      _rplottree!(tree.d1, xc, yr, zf, x, y, z, simple)
    end
  end
end




"""
    _rplottree_lr!(tree::T,
                   xc  ::Float64,
                   yr  ::UnitRange{Int64},
                   zf  ::Function,
                   x   ::Array{Float64,1},
                   y   ::Array{Float64,1},
                   z   ::Array{Float64,1}) where {T <: iT}

Returns `x` and `y` coordinates in order to plot a tree of type `iT`
where branch lengths reflect the cumulative from function `zf`.
"""
function _rplottree_lr!(tree::T,
                        xc  ::Float64,
                        yr  ::UnitRange{Int64},
                        zf  ::Function,
                        x   ::Array{Float64,1},
                        y   ::Array{Float64,1},
                        z   ::Array{Float64,1}) where {T <: iT}

  # tree δt and nsδt
  δt = dt(tree)

  # add horizontal lines
  yc = Float64(yr[1] + yr[end])*0.5

  # plot function
  zv = copy(zf(tree))

  # append
  append!(y, fill(yc, lastindex(zv)))
  append!(z, zv)
  zc     = last(zv)
  zv[1]  = xc
  append!(x, cumsum!(zv, zv))
  xc = last(x)

  push!(x, NaN)
  push!(y, NaN)
  push!(z, NaN)

  if def1(tree)
    if def2(tree)
      ntip1 = ntips(tree.d1)
      ntip2 = ntips(tree.d2)

      # add vertical lines
      push!(x, xc, xc, NaN)

      yr1 = yr[1:ntip1]
      yr2 = yr[(ntip1+1):(ntip1+ntip2)]
      push!(y, Float64(yr1[1] + yr1[end])*0.5,
               Float64(yr2[1] + yr2[end])*0.5,
               NaN)

      push!(z, zc, zc, NaN)

      _rplottree_lr!(tree.d1, xc, yr1, zf, x, y, z)
      _rplottree_lr!(tree.d2, xc, yr2, zf, x, y, z)
    else
      _rplottree_lr!(tree.d1, xc, yr, zf, x, y, z)
    end
  end
end






"""
    _rplottree!(tree::T,
                xc  ::Float64,
                yr  ::UnitRange{Int64},
                x   ::Array{Float64,1},
                y   ::Array{Float64,1}) where {T <: iTree}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplottree!(tree::T,
                     xc  ::Float64,
                     yr  ::UnitRange{Int64},
                     x   ::Array{Float64,1},
                     y   ::Array{Float64,1}) where {T <: iTree}

  # add horizontal lines
  push!(x, xc)
  xc -= e(tree)
  push!(x, xc, NaN)
  yc = (yr[1] + yr[end])*0.5
  push!(y, yc, yc, NaN)

  if def1(tree)
    if def2(tree)
      ntip1 = ntips(tree.d1)
      ntip2 = ntips(tree.d2)

      # add vertical lines
      push!(x, xc, xc, NaN)

      yr1 = yr[1:ntip1]
      yr2 = yr[(ntip1+1):(ntip1+ntip2)]

      push!(y, Float64(yr1[1] + yr1[end])*0.5,
               Float64(yr2[1] + yr2[end])*0.5,
               NaN)
      _rplottree!(tree.d1, xc, yr1, x, y)
      _rplottree!(tree.d2, xc, yr2, x, y)
    else
      _rplottree!(tree.d1, xc, yr, x, y)
    end
  end
end




"""
    f(tree::T;
      shownodes  = (T <: iTf),
      showlabels = (T <: Tlabel),
      tip        = false,
      speciation = false,
      extinct    = false,
      fossil     = true,
      textsize   = 8,
      type       = :phylogram) where {T <: iTree}

Recipe for plotting a Type `iTree`. Displays type-specific nodes if `shownodes
== true`. True by default for `sTf` trees to make sampled ancestors visible.
"""
@recipe function f(tree::T;
                   shownodes  = (T <: iTf),
                   showlabels = (T <: Tlabel),
                   tip        = false,
                   speciation = false,
                   extinct    = false,
                   fossil     = true,
                   textsize   = 8,
                   type       = :phylogram) where {T <: iTree}

  x = Float64[]
  y = Float64[]

  th  = treeheight(tree)
  nts = ntips(tree)

  _rplottree!(tree, th, 1:nts, x, y)

  ntF = Float64(nts)

  if type === :phylogram

    xlims           --> (-th*0.05, th*1.05)
    ylims           --> (1.0-(0.05*ntF), ntF+(0.05*ntF))
    xguide          --> "time"
    xflip           --> true
    fontfamily      --> :Helvetica
    tickfontfamily  --> :Helvetica
    tickfontsize    --> 8
    xtick_direction --> :out

 elseif type === :radial

    x, y = append_forradial(x, y, 50)
    polar_coords!(x, y, 360.0/ntF, th)

    xlims           --> (-th*1.05, th*1.05)
    ylims           --> (-th*1.05, th*1.05)
    xticks          --> (nothing)
    xshowaxis       --> false
  else
    @error "$type must be either phylogram of radial"
  end

  # plot defaults
  legend          --> false
  seriescolor     --> :black
  grid            --> :off
  yticks          --> (nothing)
  yshowaxis       --> false

  if shownodes

    xN = Float64[]
    yN = Float64[]

    th = treeheight(tree)
    nt = ntips(tree)

    _rplottree!(tree, th, 1:nt, xN, yN)

    shape = Symbol[:circle]
    col   = Symbol[:pink]
    alpha =
      Float64[(0.5+0.5*(!isdefined(tree, :fx) || isfix(tree))) *
              Float64(speciation)]
    _nodeproperties!(tree, shape, col, alpha,
      Float64(tip), Float64(speciation), Float64(extinct), Float64(fossil))

    @series begin
      markershape --> shape
      markercolor --> col
      markeralpha --> alpha
      markersize  --> 2.0
      if type === :phylogram
        xN, yN
      elseif type === :radial
        polar_coords!(xN, yN, 360.0/ntF, th)
        xN, yN
      end
    end
  end

  if showlabels
    labels = String[]
    _tiplabels!(tree, labels)

    @series begin
      seriestype         := :scatter
      primary            := false
      markercolor        := :black
      markershape        := :circle
      markersize         := 0
      markeralpha        := fill(0.0,nts)
      series_annotations := map(x -> (x, :Helvetica, :left, textsize, :black), labels)

      xa = fill(0.0 - 0.02*th, nts)
      ya = collect(1.0:1.0:ntF)

      if type === :phylogram
        xa, ya
      elseif type === :radial
        polar_coords!(xa, ya, 360.0/ntF, th)
        xa, ya
      end
    end
  end

  return x, y
end




"""
    append_forradial(x::Vector{Float64}, y::Vector{Float64}, n::Int64)

Appends `n` new data for making circle into `x` and `y`.
"""
function append_forradial(x::Vector{Float64},
                          y::Vector{Float64}, n::Int64)

  nnan = div(lastindex(y),3)
  nani = 1
  i1   = 1
  i2   = 2

  while nani != nnan
    y1 = y[i1]
    y2 = y[i2]

    if y1 === y2
      i1   += 3
      i2   += 3
      nani += 1
      continue
    else
      lr = collect(LinRange(y1, y2, n))
      y = append!(y[1:(i1-1)], lr,          y[(i2+1):end])
      x = append!(x[1:i1], fill(x[i1], n-2), x[i2:end])
      i1   += n+1
      i2   += n+1
      nani += 1
    end
  end
  return x, y
end




"""
    append_forradial(x::Vector{Float64},
                     y::Vector{Float64},
                     z::Vector{Float64},
                     n::Int64)

Appends `n` new data for making circle into `x` and `y`.
"""
function append_forradial(x::Vector{Float64},
                          y::Vector{Float64},
                          z::Vector{Float64},
                          n::Int64)

  nnan = count(isnan, y)
  nani = 1
  i1   = 1
  i2   = 2

  while nani != nnan
    y1 = y[i1]
    y2 = y[i2]

    if y1 === y2
      i1    = findnext(isnan, y, i2+1) + 1
      i2    = findnext(isnan, y, i1+2) - 1
      nani += 1
      continue
    else
      n = ceil(Int64, 2.0*abs(y1 - y2))
      lr = collect(LinRange(y1, y2, n))
      y = append!(y[1:(i1-1)], lr,          y[(i2+1):end])
      x = append!(x[1:i1], fill(x[i1], n-2), x[i2:end])
      z = append!(z[1:i1], fill(z[i1], n-2), z[i2:end])

      i1   += n+1
      i2   += n+1
      nani += 1
    end
  end

  return x, y, z
end




"""
    polar_coords!(x::Vector{Float64}, y::Vector{Float64}, α::Float64)

Transform `x` and `y` cartesian coordinates into polar coordinates.
"""
function polar_coords!(x ::Vector{Float64},
                       y ::Vector{Float64},
                       α ::Float64,
                       th::Float64)

  @simd for i in Base.OneTo(lastindex(x))
    x[i]  = th - x[i]
    r     = x[i]
    a     = α * y[i]
    x[i] *= cos(a*π/180.0)
    y[i]  = r * sin(a*π/180.0)
  end
end




"""
    _nodeproperties!(tree      ::T,
                     shape     ::Vector{Symbol},
                     col       ::Vector{Symbol},
                     alpha     ::Vector{Float64},
                     tip       ::Float64,
                     speciation::Float64,
                     extinct   ::Float64,
                     fossil    ::Float64) where {T <: iTree}

Completes the lists of node shapes, colors and alphas according to their
properties.
"""
function _nodeproperties!(tree      ::T,
                          shape     ::Vector{Symbol},
                          col       ::Vector{Symbol},
                          alpha     ::Vector{Float64},
                          tip       ::Float64,
                          speciation::Float64,
                          extinct   ::Float64,
                          fossil    ::Float64) where {T <: iTree}

  fx = !isdefined(tree, :fx) || isfix(tree)

  if def1(tree)
    if def2(tree)
      # speciation event
      push!(shape, :circle, fill(:none,5)...)
      push!(col, :gray, fill(:white,5)...)
      push!(alpha, (0.5+0.5*fx)*speciation, 0, 0, 0, 0, 0)
      _nodeproperties!(tree.d1, shape, col, alpha,
        tip, speciation, extinct, fossil)
      push!(shape, :none, :none)
      push!(col, :white, :white)
      push!(alpha, 0, 0)
      _nodeproperties!(tree.d2, shape, col, alpha,
        tip, speciation, extinct, fossil)
    else
      push!(shape, :none, :none, :square)
      push!(col, :white, :white, :purple)
      push!(alpha, 0, 0, (0.5+0.5*fx)*fossil)
      _nodeproperties!(tree.d1, shape, col, alpha,
        tip, speciation, extinct, fossil)
    end
  else
    # tip
    if isfossil(tree)
      push!(shape, :square); push!(col, :purple)
      push!(alpha, (0.5+0.5*fx)*fossil)
    elseif isextinct(tree)
      push!(shape, :circle); push!(col, :blue)
      push!(alpha, (0.5+0.5*fx)*extinct)
    else
      push!(shape, :circle); push!(col, :blue)
      push!(alpha, (0.5+0.5*fx)*tip)
    end
  end

  return nothing
end




#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# diversity through time
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




"""
    f(nt::Ltt)

Recipe for plotting lineage through time plots of type `Ltt`.
"""
@recipe function f(nt::Ltt)

  x = nt.t
  y = nt.n

  # plot defaults
  legend          --> false
  xguide          --> "time"
  yguide          --> "N lineages"
  xflip           --> true
  seriescolor     --> :black
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  tick_direction  --> :out
  seriestype      --> :steppost
  if maximum(y) >= 10.0
    yscale         --> :log10
  end

  return  x, y
end




"""
    f(nts::Vector{Ltt})

Recipe for plotting lineage through time plots of type `Ltt`.
"""
@recipe function f(nts::Vector{Ltt})

  x = Float64[]
  y = Float64[]

  for nt in nts
    append!(x, nt.t, NaN)
    append!(y, nt.n, NaN)
  end

  # plot defaults
  legend          --> false
  xguide          --> "time"
  yguide          --> "Diversity"
  xflip           --> true
  seriescolor     --> :black
  seriesalpha     --> min(1.0, 10.0/Float64(lastindex(nts)))
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  tick_direction  --> :out
  seriestype      --> :steppost
  if maximum(x -> isnan(x) ? 1.0 : x, y) >= 10.0
    yscale         --> :log10
  end

  return  x, y
end




"""
    f(nts::Vector{Ltt}, tdt::Float64)

Recipe for plotting lineage through time plots of type `Ltt`.
"""
@recipe function f(nts::Vector{Ltt}, 
                   tdt::Float64, 
                   q0 = [0.025, 0.975],
                   q1 = [0.25,  0.75],
                   q2 = Float64[])

  n  = lastindex(nts)
  th = maximum(map(x -> maximum(x.t), nts))

  # make time vector (present = 0.0)
  ts  = [0.0:tdt:th...]
  lts = length(ts)

  q = zeros(Int64, lts, n)
  for j in Base.OneTo(n), i in Base.OneTo(lts)
    q[i,j] = nspt(nts[j], ts[i])
  end

  if !isempty(q0) 
    Q0 = Array{Float64}(undef, lts, 2)
  end
  if !isempty(q1) 
    Q1 = Array{Float64}(undef, lts, 2)
  end
  if !isempty(q2) 
    Q2 = Array{Float64}(undef, lts, 2)
  end

  M = Array{Float64}(undef, lts)
  for i in Base.OneTo(lts)
    qi = q[i,:]
    filter!(!isnan, qi)
    if !isempty(q0) 
      Q0[i,:] = quantile(qi, q0)
    end
    if !isempty(q1) 
      Q1[i,:] = quantile(qi, q1)
    end
    if !isempty(q2) 
      Q2[i,:] = quantile(qi, q2)
    end
    M[i] = exp.(mean(log, qi))
  end

  # plot defaults
  legend          --> false
  xguide          --> "time"
  yguide          --> "Diversity"
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  tick_direction  --> :out
  fillcolor       --> :orange
  fillalpha       --> 0.3

  if maximum(x -> isnan(x) ? 1.0 : x, M) >= 10.0
    yscale         --> :log10
  end

  if !isempty(q0)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh0 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh0, (ts[i], Q0[i,1]))
      end
      for i in lts:-1:1
        push!(sh0, (ts[i], Q0[i,2]))
      end

      # Shape(sh0)
      sh0
    end
  end

  if !isempty(q1)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh1 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh1, (ts[i], Q1[i,1]))
      end
      for i in lts:-1:1
        push!(sh1, (ts[i], Q1[i,2]))
      end

      # Shape(sh1)
      sh1
    end
  end

  if !isempty(q2)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh2 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh2, (ts[i], Q2[i,1]))
      end
      for i in lts:-1:1
        push!(sh2, (ts[i], Q2[i,2]))
      end

      # Shape(sh2)
      sh2
    end
  end

  # midline
  @series begin
    seriestype := :line
    linecolor --> "#00213a99"
    linewidth --> 1.4

    ts, M
  end
end



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# rates through time
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




"""
    f(tree::T,
      f  ::Function,
      dt  ::Float64;
      q0 = [0.025, 0.975],
      q1 = [0.25,  0.75],
      q2 = Float64[])  where {T <: iT}

Recipe for plotting values given by `f` through time for a `iT`.
"""
@recipe function f(tree::T,
                   f   ::Function,
                   dt  ::Float64;
                   q0 = [0.025, 0.975],
                   q1 = [0.25,  0.75],
                   q2 = Float64[]) where {T <: iT}

  # prepare data
  ts, r = time_rate(tree, dt, f)
  lts = lastindex(ts)
  m   = time_quantile(r, [0.5])

  # common shape plot defaults
  legend          --> :none
  xguide          --> "time"
  yguide          --> string(f)[2:end]*"(t)"
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out
  ytick_direction --> :out
  fillcolor       --> :orange
  fillalpha       --> 0.3

  if !isempty(q0)
    qr0 = time_quantile(r, q0)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh0 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh0, (ts[i], qr0[i,1]))
      end
      for i in lts:-1:1
        push!(sh0, (ts[i], qr0[i,2]))
      end

      # Shape(sh0)
      sh0
    end
  end

  if !isempty(q1)
    qr1 = time_quantile(r, q1)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh0 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh0, (ts[i], qr1[i,1]))
      end
      for i in lts:-1:1
        push!(sh0, (ts[i], qr1[i,2]))
      end

      # Shape(sh0)
      sh0
    end
  end

  if !isempty(q2)
    qr2 = time_quantile(r, q2)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh0 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh0, (ts[i], qr2[i,1]))
      end
      for i in lts:-1:1
        push!(sh0, (ts[i], qr2[i,2]))
      end

      # Shape(sh0)
      sh0
    end
  end

  # midline
  @series begin
    seriestype := :line
    linecolor --> "#00304999"
    linewidth --> 1.4

    ts, m
  end
end




"""
    function f(trees::Vector{T},
               f    ::Function,
               tdt  ::Float64;
               af = x -> quantile(x, 0.5),
               q0 = [0.025, 0.975],
               q1 = [0.25,  0.75],
               q2 = Float64[]) where {T <: iT}

Recipe for plotting values given by `f` through time for a `iT`.
"""
@recipe function f(trees::Vector{T},
                   f    ::Function,
                   tdt  ::Float64;
                   af = x -> quantile(x, 0.5),
                   q0 = [0.025, 0.975],
                   q1 = [0.25,  0.75],
                   q2 = Float64[]) where {T <: iT}

  ntrees = lastindex(trees)
  riv = Vector{Float64}[]

  lts = 0
  ts  = Float64[]
  for t in trees

    # tree extracting function
    tsi, ri = time_rate(t, tdt, f)

    # aggregating function
    ri = map(af, ri)

    if lastindex(tsi) > lts
      ts  = tsi
      reverse!(ts)
      lts = lastindex(tsi)
    end

    reverse!(ri)
    push!(riv, ri)
  end

  # estimate quantiles
  q = fill(NaN, lts, ntrees)
  # estimate quantiles
  for i in Base.OneTo(ntrees)
    ri = riv[i]
    lr = lastindex(ri)
    q[1:lr,i] = ri
  end

  if !isempty(q0) 
    Q0 = Array{Float64}(undef, lts, 2)
  end
  if !isempty(q1) 
    Q1 = Array{Float64}(undef, lts, 2)
  end
  if !isempty(q2) 
    Q2 = Array{Float64}(undef, lts, 2)
  end

  M = Array{Float64}(undef, lts)
  for i in Base.OneTo(lts)
    qi = q[i,:]
    filter!(!isnan, qi)
    if !isempty(q0) 
      Q0[i,:] = quantile(qi, q0)
    end
    if !isempty(q1) 
      Q1[i,:] = quantile(qi, q1)
    end
    if !isempty(q2) 
      Q2[i,:] = quantile(qi, q2)
    end
    M[i] = quantile(qi, 0.5)
  end

  # common shape plot defaults
  legend          --> :none
  xguide          --> "time"
  yguide          --> string(f)[2:end]*"(t)"
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out
  ytick_direction --> :out
  fillcolor       --> :orange
  fillalpha       --> 0.3

  if !isempty(q0)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh0 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh0, (ts[i], Q0[i,1]))
      end
      for i in lts:-1:1
        push!(sh0, (ts[i], Q0[i,2]))
      end

      # Shape(sh0)
      sh0
    end
  end

  if !isempty(q1)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh1 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh1, (ts[i], Q1[i,1]))
      end
      for i in lts:-1:1
        push!(sh1, (ts[i], Q1[i,2]))
      end

      # Shape(sh1)
      sh1
    end
  end

  if !isempty(q2)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      sh2 = Tuple{Float64,Float64}[]
      for i in Base.OneTo(lts)
        push!(sh2, (ts[i], Q2[i,1]))
      end
      for i in lts:-1:1
        push!(sh2, (ts[i], Q2[i,2]))
      end

      # Shape(sh2)
      sh2
    end
  end


  # midline
  @series begin
    seriestype := :line
    linecolor --> "#00213a99"
    linewidth --> 1.4

    ts, M
  end
end




#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# traits
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




"""
    function f(tree::T; type::Symbol = :trait)

Recipe for plotting with the tree or the trait evolutions for `sTX`.
"""
@recipe function f(tree::T; type = :trait) where {T <: sTX}

  x = Float64[]
  y = Float64[]

  if type === :tree

    th = treeheight(tree)
    nt = ntips(tree)

    _rplottree!(tree, th, 1:nt, x, y)

    # plot defaults
    legend          --> false
    xguide          --> "time"
    seriescolor     --> :black
    xlims           --> (-th*0.05, th*1.05)
    ylims           --> (1.0-(0.05*Float64(nt)), nt+(0.05*Float64(nt)))
    xflip           --> true
    fontfamily      --> :Helvetica
    tickfontfamily  --> :Helvetica
    tickfontsize    --> 8
    grid            --> :off
    xtick_direction --> :out
    yticks          --> (nothing)
    yshowaxis       --> false

    return x, y

  elseif type === :trait

    th = treeheight(tree)

    _rplottrait!(tree, th, x, y)

    yfilt = filter(x -> !isnan(x), y)

    # plot defaults
    legend          --> false
    xguide          --> "time"
    yguide          --> "trait"
    seriescolor     --> :purple
    xlims           --> (-th*0.05, th*1.05)
    ylims           --> (minimum(yfilt), maximum(yfilt))
    xflip           --> true
    fontfamily      --> :Helvetica
    tickfontfamily  --> :Helvetica
    tickfontsize    --> 8
    grid            --> :off
    xtick_direction --> :out

    return x, y
  else

    @warn "$type neither tree nor trait"
  end

end





"""
    function f(tree::T; type::Symbol = :trait)

Recipe for plotting with the tree or the trait evolutions for `sTX`.
"""
@recipe function f(trees::Vector{T}) where {T <: sTX}

  x = Float64[]
  y = Float64[]

  th = treeheight(trees[1])
  n  = lastindex(trees)

  for t in trees
    _rplottrait!(t, treeheight(t), x, y)
  end

  yfilt = filter(x -> !isnan(x), y)

  # plot defaults
  legend          --> false
  xguide          --> "time"
  yguide          --> "trait"
  seriescolor     --> :purple
  seriesalpha     --> min(1.0, 10.0/Float64(n))
  xlims           --> (-th*0.05, th*1.05)
  ylims           --> (minimum(yfilt), maximum(yfilt))
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out

  return x, y
end




"""
    _rplottree!(tree::T,
                xc  ::Float64,
                x   ::Array{Float64,1},
                y   ::Array{Float64,1}) where {T <: sTX}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplottrait!(tree::T,
                      xc  ::Float64,
                      x   ::Array{Float64,1},
                      y   ::Array{Float64,1}) where {T <: sTX}

  # add horizontal lines
  push!(x, xc)
  xc -= e(tree)
  push!(x, xc, NaN)
  push!(y, xi(tree), xf(tree), NaN)

  if def1(tree)
    _rplottrait!(tree.d1, xc, x, y)
    if def2(tree)
      _rplottrait!(tree.d2, xc, x, y)
    end
  end
end




"""
    function f(tree::T; type::Symbol = :trait)

Recipe for plotting with the tree or the trait evolutions for `sTX`.
"""
@recipe function f(tree::T) where {T <: iTX}

  x = Float64[]
  y = Float64[]

  th = treeheight(tree)

  _rplottrait!(tree, th, x, y)

  yfilt = filter(x -> !isnan(x), y)

  # plot defaults
  legend          --> false
  xguide          --> "time"
  yguide          --> "trait"
  seriescolor     --> :purple
  xlims           --> (-th*0.05, th*1.05)
  ylims           --> (minimum(yfilt), maximum(yfilt))
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out

  return x, y
end




"""
    function f(tree::T; type::Symbol = :trait)

Recipe for plotting with the tree or the trait evolutions for `sTX`.
"""
@recipe function f(trees::Vector{T}) where {T <: iTX}

  x = Float64[]
  y = Float64[]

  th = treeheight(trees[1])
  n  = lastindex(trees)

  for t in trees
    _rplottrait!(t, treeheight(t), x, y)
  end

  yfilt = filter(x -> !isnan(x), y)

  # plot defaults
  legend          --> false
  xguide          --> "time"
  yguide          --> "trait"
  seriescolor     --> :purple
  seriesalpha     --> min(1.0, 10.0/Float64(n))
  xlims           --> (-th*0.05, th*1.05)
  ylims           --> (minimum(yfilt), maximum(yfilt))
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out

  return x, y
end



"""
    _rplottree!(tree::T,
                xc  ::Float64,
                x   ::Array{Float64,1},
                y   ::Array{Float64,1}) where {T <: iTX}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplottrait!(tree::T,
                      xc  ::Float64,
                      x   ::Array{Float64,1},
                      y   ::Array{Float64,1}) where {T <: iTX}

  # tree δt and nsδt
  δt = dt(tree)

  # add horizontal lines
  xvi = xv(tree)
  l   = lastindex(xvi)
  @simd for i in Base.OneTo(l-1)
    push!(x, xc - Float64(i-1)*δt)
    push!(y, xvi[i])
  end

  push!(x, xc - (Float64(l-2)*δt + fdt(tree)), NaN)
  push!(y, xvi[l], NaN)

  xc -= e(tree)

  if def1(tree)
    _rplottrait!(tree.d1, xc, x, y)
    if def2(tree)
      _rplottrait!(tree.d2, xc, x, y)
    end
  end
end

