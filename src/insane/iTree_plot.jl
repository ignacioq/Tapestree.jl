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
  `b`  birth (speciation or initiation) rates
  `d`  death (extinction) rates
  `c`  completion rates [for protracted models]
  `lb` log birth (speciation or initiation) rates
  `ld` log death (extinction) rates
  `lc` log completion rates [for protracted models]
  `t`  turnover
  `lt` log turnover
  `nd` net diversification
  `dμ` change in speciation rates
"""
b(tree::iT)     = exp.(lλ(tree))
d(tree::iT)     = exp.(lμ(tree))
b(tree::iTpbd)  = exp.(lb(tree))
c(tree::iTpbd)  = exp.(lλ(tree))
lb(tree::iT)    = lλ(tree) # lb(tree::iTpbd) is defined independently
ld(tree::iT)    = lμ(tree)
lc(tree::iTpbd) = lλ(tree)
t(tree::iT)     = d(tree) ./ b(tree)
lt(tree::iT)    = log.(d(tree) ./ b(tree))
nd(tree::iT)    = b(tree) .- d(tree)
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
    f(tree::T;
      labsize    = 8,
      type       = :phylogram,
      showlabels = (T <: Tlabel),
      shownodes  = (false, false, T <: Union{iTf, iTpbd}),
      shapes     = [:none, :none, :square],
      colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
      shsizes    = [3.0, 3.0, 3.0],
      showda     = false,
      col_da     = ["#a9a9a9", :black]) where {T <: iTree}

Recipe for plotting a Type `iTree`. Displays type-specific nodes if `shownodes
= true`. True by default for `sTf` trees to make sampled ancestors visible.
"""
@recipe function f(tree::T;
                   labsize    = 8,
                   type       = :phylogram,
                   showlabels = (T <: Tlabel),
                   shownodes  = (false, false, T <: Union{iTf, iTpbd}),
                   shapes     = [:none, :none, :square],
                   colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
                   shsizes    = [3.0, 3.0, 3.0],
                   showda     = false,
                   col_da     = ["#a9a9a9", :black]) where {T <: iTree}

  x  = Float64[]
  y  = Float64[]
  z  = Float64[]
  nodet = Int64[]   # 0 = speciation, 1 = extinction, 2 = fossilization
  xnode = Float64[]
  ynode = Float64[]
  zstyle = Symbol[]

  th  = treeheight(tree)
  nts = ntips(tree)

  _rplottree!(tree, th, 0, x, y, z, zstyle, nodet, xnode, ynode, shownodes)
  # pop!(zstyle)
  # @show x
  # @show y
  # @show z
  # @show zstyle

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

  if typeof(tree) == iTpbd
    linestyle     --> zstyle
  end

  @series begin
    seriestype  := :path

    if showda
      line_z     --> z
      linecolor  --> palette(col_da, 2)
      return x, y, z
    end

    return x, y
  end

  if any(shownodes)

    @series begin
      seriestype  := :scatter
      markershape -->       shapes[nodet]
      markercolor -->       colors[nodet]
      markerstrokecolor --> colors[nodet]
      markersize  -->       shsizes[nodet]
      if type === :phylogram
        xnode, ynode
      elseif type === :radial
        polar_coords!(xnode, ynode, 360.0/ntF, th)
        xnode, ynode
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
      series_annotations := map(x -> (x, :Helvetica, :left, labsize, :black), labels)

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
end




"""
    _rplottree!(tree ::T,
                xc   ::Float64,
                nn   ::Int64,
                nx   ::Int64,
                x    ::Array{Float64,1},
                y    ::Array{Float64,1},
                z    ::Array{Float64,1},
                zstyle ::Array{Symbol,1},
                nodet::Array{Int64,1},
                xnode::Array{Float64,1},
                ynode::Array{Float64,1},
                show ::NTuple{3,Bool}) where {T <: iTree}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree` and 
`z` vector differentiating fixed `1` from data augmented `0` components.
"""
function _rplottree!(tree ::T,
                     xc   ::Float64,
                     i    ::Int64,
                     x    ::Array{Float64,1},
                     y    ::Array{Float64,1},
                     z    ::Array{Float64,1},
                     zstyle ::Array{Symbol,1},
                     nodet::Array{Int64,1},
                     xnode::Array{Float64,1},
                     ynode::Array{Float64,1},
                     show ::NTuple{3,Bool}) where {T <: iTree}

  xe = xc - e(tree)

  if def1(tree)
    if def2(tree)

      y1, i = _rplottree!(tree.d1, xe, i, x, y, z, zstyle, nodet, xnode, ynode, show)
      y2, i = _rplottree!(tree.d2, xe, i, x, y, z, zstyle, nodet, xnode, ynode, show)

      yc = (y1 + y2)*0.5

      # add vertical lines
      push!(x, xe, xe, NaN, xe, xe, NaN)
      push!(y, y1, yc, NaN, yc, y2, NaN)

      z1 = Float64(isfix(tree.d1))
      z2 = Float64(isfix(tree.d2))
      push!(z, z1, z1, NaN, z2, z2, NaN)
      
      zstyle1 = ifelse(isgood(tree.d1), :solid, :dot)
      zstyle2 = ifelse(isgood(tree.d2), :solid, :dot)
      push!(zstyle, zstyle1, zstyle1, :solid, zstyle2, zstyle2, :solid)

      # nodes
      if show[1]
        push!(nodet, 1)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    else
      yc, i = _rplottree!(tree.d1, xe, i, x, y, z, zstyle, nodet, xnode, ynode, show)

      if show[3]
        push!(nodet, 3)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    end
  else
    i += 1
    yc = Float64(i)
    if isextinct(tree)
      if show[2]
        push!(nodet, 2)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    elseif isfossil(tree)
      if show[3]
        push!(nodet, 3)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    end
  end

  # add horizontal lines
  push!(x, xc, xe, NaN)
  push!(y, yc, yc, NaN)
  zc = Float64(isfix(tree))
  push!(z, zc, zc, NaN)    
  zstylec = ifelse(isgood(tree), :solid, :dot)
  push!(zstyle, zstylec, zstylec, :solid)

  return yc, i
end




"""
    append_forradial(x::Vector{Float64}, y::Vector{Float64}, n::Int64)

Appends `n` new data for making circle into `x` and `y`.
"""
function append_forradial(x::Vector{Float64},
                          y::Vector{Float64}, 
                          n::Int64)

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
    f(tree::T,
      zf  ::Function;
      type       = :phylogram,
      showlabels = (T <: Tlabel),
      shownodes  = (false, false, (T <: Union{iTf, iTpbd})),
      shapes     = [:circle, :circle, :square],
      colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
      shsizes    = [0.0, 0.0, 3.0],
      simple     = false) where {T <: iT}

Recipe for plotting a Type `iT`.
"""
@recipe function f(tree::T,
                   zf  ::Function;
                   type       = :phylogram,
                   showlabels = (T <: Tlabel),
                   shownodes  = (false, false, T <: Union{iTf, iTpbd}),
                   shapes     = [:none, :none, :square],
                   colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
                   shsizes    = [0.0, 0.0, 3.0],
                   simple     = false) where {T <: iT}

  x = Float64[]
  y = Float64[]
  z = Float64[]
  nodet = Int64[]   # 0 = speciation, 1 = extinction, 2 = fossilization
  xnode = Float64[]
  ynode = Float64[]

  th  = treeheight(tree)
  nts = ntips(tree)

  _rplottree!(tree, zf, th, 0, x, y, z, nodet, xnode, ynode, shownodes, simple)

  ntF = Float64(nts)

  # plot defaults
  fontfamily              --> :Helvetica
  legend                  --> :none
  colorbar_tickfontfamily --> :Helvetica
  grid                    --> :off

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

  @series begin
    seriestype  := :path
    line_z     --> z
    linecolor  --> cgrad(:roma, rev = true)

    if type === :rates
      return x, z
    else
      return x, y
    end
  end

  if type != :rates && any(shownodes)
    @series begin
      seriestype  := :scatter
      markershape -->       shapes[nodet]
      markercolor -->       colors[nodet]
      markerstrokecolor --> colors[nodet]
      markersize  -->       shsizes[nodet]
      if type === :phylogram
        xnode, ynode
      elseif type === :radial
        polar_coords!(xnode, ynode, 360.0/ntF, th)
        xnode, ynode
      end
    end
  end
end




"""
    _rplottree!(tree  ::T,
                zf    ::Function,
                xc    ::Float64,
                i     ::Int64,
                x     ::Array{Float64,1},
                y     ::Array{Float64,1},
                z     ::Array{Float64,1},
                nodet ::Array{Int64,1},
                xnode ::Array{Float64,1},
                ynode ::Array{Float64,1},
                show  ::NTuple{3,Bool},
                simple::Bool) where {T <: iT}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplottree!(tree  ::T,
                     zf    ::Function,
                     xc    ::Float64,
                     i     ::Int64,
                     x     ::Array{Float64,1},
                     y     ::Array{Float64,1},
                     z     ::Array{Float64,1},
                     nodet ::Array{Int64,1},
                     xnode ::Array{Float64,1},
                     ynode ::Array{Float64,1},
                     show  ::NTuple{3,Bool},
                     simple::Bool) where {T <: iT}

  xe = xc - e(tree)

  if def1(tree)
     if def2(tree)

      y1, i = _rplottree!(tree.d1, zf, xe, i, x, y, z, nodet, xnode, ynode, 
        show, simple)
      y2, i = _rplottree!(tree.d2, zf, xe, i, x, y, z, nodet, xnode, ynode,
        show, simple)

      yc = (y1 + y2)*0.5
      zc = last(zf(tree))

      # add vertical lines
      push!(x, xe, xe, NaN)
      push!(y, y1, y2, NaN)
      push!(z, zc, zc, NaN)

      # nodes
      if show[1]
        push!(nodet, 1)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    else

      yc, i = _rplottree!(tree.d1, zf, xe, i, x, y, z, nodet, xnode, ynode, 
        show, simple)

      if show[3]
        push!(nodet, 3)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    end
  else
    i += 1
    yc = Float64(i)
    if isextinct(tree)
      if show[2]
        push!(nodet, 2)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    elseif isfossil(tree)
      if show[3]
        push!(nodet, 3)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    end
  end

  # tree δt and nsδt
  δt = dt(tree)

  # add horizontal lines
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

  return yc, i
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
function polar_coords!(x  ::Vector{Float64},
                       y  ::Vector{Float64},
                       α  ::Float64,
                       tor::Float64)

  @simd for i in Base.OneTo(lastindex(x))
    x[i]  = tor - x[i]
    r     = x[i]
    a     = α * y[i]
    x[i] *= cos(a*π/180.0)
    y[i]  = r * sin(a*π/180.0)
  end
end




"""
    plotω(tree       ::T,
          ωtimes     ::Vector{Float64};
          labsize    = 8,
          type       = :phylogram,
          showlabels = (T <: Tlabel),
          shownodes  = (false, false, T <: Union{iTf, iTpbd}),
          shapes     = [:none, :none, :square],
          colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
          showda     = false,
          col_da     = ["#a9a9a9", :black],
          yω         = 0.98+(1.1/50-0.05)*ntips(tree)) where {T <: iTree}

Recipe for plotting a tree with fossil occurrences.
"""
function plotω(tree       ::T,
               ωtimes     ::Vector{Float64};
               labsize    = 8,
               type       = :phylogram,
               showlabels = (T <: Tlabel),
               shownodes  = (false, false, T <: Union{iTf, iTpbd}),
               shapes     = [:none, :none, :square],
               colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
               showda     = false,
               col_da     = ["#a9a9a9", :black],
               yω         = 0.98+(1.1/50-0.05)*ntips(tree)) where {T <: iTree}

  plot(tree, labsize=labsize, type=type, showlabels=showlabels, shownodes=shownodes, shapes=shapes,
       colors=colors, showda=showda, col_da=col_da)
  ymin = 1.0-0.05*ntips(tree)
  scatter!(ωtimes, [max(rnorm(yω, (yω-ymin)/5), ymin) for ωt in ωtimes], label="occurrences", mc=:grey, ms=2, ma=0.5)

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
  xflip           --> true
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

  return x, y
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
    # M[i] = exp.(mean(log, qi))
    M[i] = median(qi)
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




"""
    plotω(LTT        ::Ltt,
          ωtimes     ::Vector{Float64})

    plotω(LTT        ::Vector{Ltt},
          ωtimes     ::Vector{Float64})
    
    plotω(LTT        ::Vector{Ltt},
          ωtimes     ::Vector{Float64}, 
          tdt        ::Float64, 
          q0         = [0.025, 0.975],
          q1         = [0.25,  0.75],
          q2         = Float64[])

Recipe for plotting lineage through time plots of type `Ltt`, together with fossil occurrences.
"""
function plotω(LTT        ::Ltt,
               ωtimes     ::Vector{Float64})

  plot(LTT)
  scatter!(ωtimes, [1+abs(rnorm(0.0, maximum(LTT.n)/500)) for ωt in ωtimes], label="occurrences", mc=:grey, ms=1.5, ma=0.5)

end

function plotω(LTT        ::Vector{Ltt},
               ωtimes     ::Vector{Float64})

  plot(LTT)
  scatter!(ωtimes, [1+abs(rnorm(0.0, maximum([maximum(LTTi.n) for LTTi in LTT])/500)) for ωt in ωtimes], label="occurrences", mc=:grey, ms=1.5, ma=0.5)

end

function plotω(LTT        ::Vector{Ltt},
               ωtimes     ::Vector{Float64}, 
               tdt        ::Float64, 
               q0         = [0.025, 0.975],
               q1         = [0.25,  0.75],
               q2         = Float64[])

  plot(LTT, tdt, q0, q1, q2)
  scatter!(ωtimes, [1+abs(rnorm(0.0, maximum([maximum(LTTi.n) for LTTi in LTT])/500)) for ωt in ωtimes], label="occurrences", mc=:grey, ms=1.5, ma=0.5)

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
      f   ::Function,
      dt  ::Float64;
      t_af = mean,
      fillcolor = :orange,
      linecolor = "#00304999",
      q0   = [0.025, 0.975],
      q1   = [0.25,  0.75],
      q2   = Float64[]) where {T <: iT}

Recipe for plotting values given by `f` through time for a `iT`.
"""
@recipe function f(tree::T,
                   f   ::Function,
                   δt  ::Float64;
                   t_af = mean,
                   fillcolor = :orange,
                   linecolor = "#00304999",
                   q0   = [0.025, 0.975],
                   q1   = [0.25,  0.75],
                   q2   = Float64[]) where {T <: iT}

  # prepare data
  ts, r = time_rate(tree, f, δt)
  lts = lastindex(ts)
  m   = map(t_af, r)

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
  fillcolor       --> fillcolor
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
    linecolor --> linecolor
    linewidth --> 1.4

    ts, m
  end
end




"""
    function f(trees::Vector{T},
               f    ::Function,
               δt   ::Float64;
               t_af  = mean,
               fillcolor = :orange,
               linecolor = "#00213a99",
               tv_af = x -> quantile(x, 0.5),
               q0    = [0.025, 0.975],
               q1    = [0.25,  0.75],
               q2    = Float64[]) where {T <: iT}

Recipe for plotting values given by `f` through time for a `iT`.
"""
@recipe function f(trees::Vector{T},
                   f    ::Function,
                   δt   ::Float64;
                   t_af  = mean,
                   fillcolor = :orange,
                   linecolor = "#00213a99",
                   tv_af = x -> quantile(x, 0.5),
                   q0    = [0.025, 0.975],
                   q1    = [0.25,  0.75],
                   q2    = Float64[]) where {T <: iT}

  ntrees = lastindex(trees)
  riv = Vector{Float64}[]

  lts = 0
  ts  = Float64[]
  for t in trees

    # tree extracting function
    tsi, ri = time_rate(t, f, δt)

    # aggregating function
    ri = map(t_af, ri)

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
    isempty(qi) && continue

    if !isempty(q0) 
      Q0[i,:] = quantile(qi, q0)
    end
    if !isempty(q1) 
      Q1[i,:] = quantile(qi, q1)
    end
    if !isempty(q2) 
      Q2[i,:] = quantile(qi, q2)
    end
    M[i] = tv_af(qi)
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
  fillcolor       --> fillcolor
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
    linecolor --> linecolor
    linewidth --> 1.4

    ts, M
  end
end




"""
    f(rates::Vector{Float64},
      tor::Float64;
      fillcolor = :orange,
      linecolor = "#00304999",
      tv_af = x -> quantile(x, 0.5),
      q0 = [0.025, 0.975],
      q1 = [0.25,  0.75],
      q2 = Float64[])

Recipe for plotting constant rates through time.
"""
@recipe function f(rates::Vector{Float64},
                   tor::Float64;
                   fillcolor = :orange,
                   linecolor = "#00304999",
                   tv_af = x -> quantile(x, 0.5),
                   q0 = [0.025, 0.975],
                   q1 = [0.25,  0.75],
                   q2 = Float64[])

  # Compute quantiles for uncertainty bands
  if !isempty(q0)
    Q0 = quantile(rates, q0)
  end
  if !isempty(q1)
    Q1 = quantile(rates, q1)
  end
  if !isempty(q2)
    Q2 = quantile(rates, q2)
  end

  M = tv_af(rates)

  # Common plot defaults
  legend          --> :none
  xguide          --> "time"
  yguide          --> "rate(t)"
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out
  ytick_direction --> :out
  fillcolor       --> fillcolor
  fillalpha       --> 0.3

  # Plot uncertainty bands as rectangles over [0, tor]

  if !isempty(q0)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      x_values = [0.0, tor, tor, 0.0]
      y_values = [Q0[1], Q0[1], Q0[2], Q0[2]]
      x_values, y_values
    end
  end

  if !isempty(q1)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      x_values = [0.0, tor, tor, 0.0]
      y_values = [Q1[1], Q1[1], Q1[2], Q1[2]]
      x_values, y_values
    end
  end

  if !isempty(q2)
    @series begin
      seriestype := :shape
      linecolor  := nothing

      x_values = [0.0, tor, tor, 0.0]
      y_values = [Q2[1], Q2[1], Q2[2], Q2[2]]
      x_values, y_values
    end
  end

  # Midline for the constant rate
  @series begin
    seriestype := :line
    linecolor --> linecolor
    linewidth --> 1.4
    x_values = [0.0, tor]
    y_values = [M, M]
    x_values, y_values
  end
end




"""
    f(rates::Vector{Vector{Float64}},
      tor::Float64,
      ψω_epoch::Vector{Float64};
      fillcolor = :saddlebrown,
      linecolor = :saddlebrown,
      tv_af = x -> quantile(x, 0.5),
      q0 = [0.025, 0.975],
      q1 = [0.25,  0.75],
      q2 = Float64[])

Recipe for plotting piecewise-constant rates through time.
"""
@recipe function f(rates::Vector{Vector{Float64}},
                   tor::Float64,
                   ψω_epoch::Vector{Float64};
                   fillcolor = :saddlebrown,
                   linecolor = :saddlebrown,
                   tv_af = x -> quantile(x, 0.5),
                   q0 = [0.025, 0.975],
                   q1 = [0.25,  0.75],
                   q2 = Float64[])

  # Ensure breakpoints start at tor and end at 0
  bp = vcat([tor], ψω_epoch, [0.0])

  # Common plot defaults
  legend          --> :none
  xguide          --> "time"
  yguide          --> "rate(t)"
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out
  ytick_direction --> :out
  fillcolor       --> fillcolor
  fillalpha       --> 0.3

  # Loop over each interval and plot the rate and uncertainty bands
  for i in 1:length(rates)
    t_start = bp[i]
    t_end = bp[i+1]
    rate_samples = rates[i]  # Posterior samples for interval i

    # Compute central tendency and quantiles
    M = tv_af(rate_samples)

    if !isempty(q0)
      Q0 = quantile(rate_samples, q0)
    end
    if !isempty(q1)
      Q1 = quantile(rate_samples, q1)
    end
    if !isempty(q2)
      Q2 = quantile(rate_samples, q2)
    end

    # Plot uncertainty bands as rectangles for the current interval
    if !isempty(q0)
      @series begin
        seriestype := :shape
        linecolor  := nothing

        x_values = [t_start, t_end, t_end, t_start]
        y_values = [Q0[1], Q0[1], Q0[2], Q0[2]]
        x_values, y_values
      end
    end

    if !isempty(q1)
      @series begin
        seriestype := :shape
        linecolor  := nothing

        x_values = [t_start, t_end, t_end, t_start]
        y_values = [Q1[1], Q1[1], Q1[2], Q1[2]]
        x_values, y_values
      end
    end

    if !isempty(q2)
      @series begin
        seriestype := :shape
        linecolor  := nothing

        x_values = [t_start, t_end, t_end, t_start]
        y_values = [Q2[1], Q2[1], Q2[2], Q2[2]]
        x_values, y_values
      end
    end

    # Midline for the current interval
    @series begin
      seriestype := :line
      linecolor --> linecolor
      linewidth --> 1.4
      x_values = [t_start, t_end]
      y_values = [M, M]
      x_values, y_values
    end
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

Returns `x` and `y` coordinates in order to plot a tree of type `sTX`.
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

Returns `x` and `y` coordinates in order to plot a tree of type `sTX`.
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

