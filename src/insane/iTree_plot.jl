#=

insane tree plot

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#




"""
    f(tree::T, lv::Function, dt::Float64, e::Bool) where {T <: iTgbm}

Recipe for plotting values given by `lv` through time for a `iTgbm`.
"""
@recipe function f(tree::T, 
                   lv  ::Function, 
                   dt  ::Float64,
                   e   ::Bool) where {T <: iTgbm}

  # prepare data
  ts, r = time_rate(tree, dt, lv)

  if e
    fx = exp.(time_quantile(r, [0.0, 0.25, 0.5, 0.75, 1.0]))
  else
    fx = time_quantile(r, [0.0, 0.25, 0.5, 0.75, 1.0])
  end

  lf = size(fx,1)

  # common shape plot defaults
  legend          --> :none
  xguide          --> "time"
  yguide          --> string(lv)[2:end]*"(t)"
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out
  ytick_direction --> :out
  fillcolor       --> :orange
  fillalpha       --> 0.3

  # range shape
  @series begin
    seriestype := :shape
    linecolor  := nothing

    sh0 = Tuple{Float64,Float64}[]
    for i in Base.OneTo(lf)
      push!(sh0, (ts[i], fx[i,1]))
    end
    for i in lf:-1:1
      push!(sh0, (ts[i], fx[i,5]))
    end

    # Shape(sh0)
    sh0
  end

  # [0.25, 0.75] quantile range shape
  @series begin
    seriestype := :shape
    linecolor  := nothing

    sh1 = Tuple{Float64,Float64}[]
    for i in Base.OneTo(lf)
      push!(sh1, (ts[i], fx[i,2]))
    end
    for i in lf:-1:1
      push!(sh1, (ts[i], fx[i,4]))
    end

    # Shape(sh1)
    sh1
  end

  # midline
  @series begin
    seriestype := :line
    linecolor --> "#00304999"
    linewidth --> 1.4

    ts, fx[:,3]
  end

end




"""
    function f(tree::T, zfun::Function) where {T <: iTgbm}

Recipe for plotting a Type `iTgbm`.
"""
@recipe function f(tree::T, zfun::Function) where {T <: iTgbm}

  x = Float64[]
  y = Float64[]
  z = Float64[]

  th = treeheight(tree)
  nt = ntips(tree)

  _rplottree!(tree, th, 1:nt, zfun, x, y, z)

  # plot defaults
  line_z          --> z
  linecolor       --> :inferno
  legend          --> :none
  colorbar        --> true
  xguide          --> "time"
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
end




"""
    _rplottree!(tree::T, 
                xc  ::Float64, 
                yr  ::UnitRange{Int64},
                zfun::Function,
                x   ::Array{Float64,1}, 
                y   ::Array{Float64,1},
                z   ::Array{Float64,1}) where {T <: iTgbm}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplottree!(tree::T, 
                     xc  ::Float64, 
                     yr  ::UnitRange{Int64},
                     zfun::Function,
                     x   ::Array{Float64,1}, 
                     y   ::Array{Float64,1},
                     z   ::Array{Float64,1}) where {T <: iTgbm}

  # tree δt and nsδt
  δt = dt(tree)

  # add horizontal lines
  yc = Float64(yr[1] + yr[end])*0.5
  zv = exp.(zfun(tree))
  l  = lastindex(zv)
  @simd for i in Base.OneTo(l-1)
    push!(x, xc - Float64(i-1)*δt)
    push!(y, yc)
    push!(z, zv[i])
  end

  push!(x, xc - (Float64(l-2)*δt + fdt(tree)), NaN)
  push!(y, yc, NaN)
  push!(z, zv[l], NaN)

  xc -= e(tree)
  
  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if defd1 && defd2

    ntip1 = ntips(tree.d1)
    ntip2 = ntips(tree.d2)

    # add vertical lines
    push!(x, xc, xc, NaN)

    yr1 = yr[1:ntip1]
    yr2 = yr[(ntip1+1):(ntip1+ntip2)]
    push!(y, Float64(yr1[1] + yr1[end])*0.5,
             Float64(yr2[1] + yr2[end])*0.5,
             NaN)

    push!(z, z[end-1], z[end-1], NaN)

    _rplottree!(tree.d1, xc, yr1, zfun, x, y, z)
    _rplottree!(tree.d2, xc, yr2, zfun, x, y, z)
  
  elseif defd1  _rplottree!(tree.d1, xc, yr, zfun, x, y, z)
  elseif defd2  _rplottree!(tree.d2, xc, yr, zfun, x, y, z)
  end

end




"""
    function f(tree::iTgbmct, zfun::Function, ϵ::Float64)

Recipe for plotting extinction on a `iTgbmct` given `ϵ`.
"""
@recipe function f(tree::iTgbmct, zfun::Function, ϵ::Float64)

  x = Float64[]
  y = Float64[]
  z = Float64[]

  th = treeheight(tree)
  nt = ntips(tree)

  _rplottree!(tree, th, 1:nt, zfun, x, y, z)

  @simd for i in Base.OneTo(lastindex(z))
    z[i] *= ϵ
  end

  # plot defaults
  line_z          --> z
  linecolor       --> :inferno
  legend          --> :none
  colorbar        --> true
  xguide          --> "time"
  xlims           --> (-th*0.05, th*1.1)
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

  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if defd1 && defd2

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
  elseif defd1
    _rplottree!(tree.d1, xc, yr, x, y)
  elseif defd2
    _rplottree!(tree.d2, xc, yr, x, y)
  end
end




"""
    f(tree::T; shownodes  = (T <: sTf),
               showlabels = (T === sT_label || T === sTf_label)) 
               where {T <: iTree}

Recipe for plotting a Type `iTree`. Displays type-specific nodes if `shownodes 
== true`. True by default for `sTf` trees to make sampled ancestors visible.
"""
@recipe function f(tree::T; 
                   shownodes  = (T <: sTf),
                   showlabels = (T === sT_label || T === sTf_label)) where {
                                                                    T <: iTree}

  x = Float64[]
  y = Float64[]

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

  if shownodes
    shape = Symbol[:circle]
    col   = Symbol[:pink]
    alpha = Float64[0.5+0.5*(!isdefined(tree, :fx) || isfix(tree))]
    _nodeproperties!(tree, shape, col, alpha)

    markershape       --> shape
    markercolor       --> col
    markeralpha       --> alpha
  end

  if showlabels
    labels = String[]
    _tiplabels!(tree, labels)

    txt = [(0.0, i, labels[i]) for i in 1:nt]

    @series begin
      seriestype         := :scatter
      primary            := false
      markercolor        := :black
      markershape        := :circle
      markersize         := 0
      markeralpha        := fill(0.0,nt)
      series_annotations := map(x -> (x, :Helvetica, :left, 8, :black), labels)
      fill(0.0 - 0.02*th, nt), 1:nt
    end
  end

  return x, y
end




"""
    _nodeproperties!(tree ::T, 
                     shape::Vector{Symbol}, 
                     col  ::Vector{Symbol}, 
                     alpha::Vector{Float64}) where {T <: iTree}

Completes the lists of node shapes, colors and alphas according to their 
properties.
"""
function _nodeproperties!(tree ::T, 
                          shape::Vector{Symbol}, 
                          col  ::Vector{Symbol}, 
                          alpha::Vector{Float64}) where {T <: iTree}

  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)
  fx = !isdefined(tree, :fx) || isfix(tree)

  if defd1 || defd2
    # not a tip

    if defd1 && defd2
      # speciation event
      push!(shape, :circle, fill(:none,5)...)
      push!(col, :gray, fill(:white,5)...)
      push!(alpha, 0.5+0.5*fx, 0, 0, 0, 0, 0)
      _nodeproperties!(tree.d1, shape, col, alpha)
      push!(shape, :none, :none)
      push!(col, :white, :white)
      push!(alpha, 0, 0)
      _nodeproperties!(tree.d2, shape, col, alpha)
    elseif isfossil(tree)
      # sampled ancestor
      push!(shape, :none, :none, :square)
      push!(col, :white, :white, :red)
      push!(alpha, 0, 0, 0.5+0.5*fx)
      _nodeproperties!(defd1 ? tree.d1 : tree.d2, shape, col, alpha)
    else
      # error
      @warn "Star node incorrectly defined (neither tip, speciation or fossil)"
      push!(shape, :none, :none, :star)
      push!(col, :white, :white, :yellow)
      push!(alpha, 0, 0, 0.5+0.5*fx)
      _nodeproperties!(defd1 ? tree.d1 : tree.d2, shape, col, alpha)
    end
  else
    # tip
    isfossil(tree) ? (push!(shape, :square); push!(col, :red)) : 
                     (push!(shape, :circle); push!(col, :blue))
    push!(alpha, 0.5+0.5*fx)
  end
  
  return nothing
end




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
  if maximum(y)>=10
    yaxis         --> :log
  end

  return  x, y
end




"""
    function f(tree::T; type::Symbol = :trait)

Recipe for plotting with the tree or the trait evolutions for `iTreeX`.
"""
@recipe function f(tree::T; type = :trait) where {T <: iTreeX}

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

Recipe for plotting with the tree or the trait evolutions for `iTreeX`.
"""
@recipe function f(trees::Vector{T}) where {T <: iTreeX}

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
                y   ::Array{Float64,1}) where {T <: iTreeX}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplottrait!(tree::T, 
                      xc  ::Float64, 
                      x   ::Array{Float64,1}, 
                      y   ::Array{Float64,1}) where {T <: iTreeX} 

  # add horizontal lines
  push!(x, xc)
  xc -= e(tree)
  push!(x, xc, NaN)
  push!(y, xi(tree), xf(tree), NaN)

  defd1 = isdefined(tree, :d1)
  defd2 = isdefined(tree, :d2)

  if defd1 && defd2
    _rplottrait!(tree.d1, xc, x, y)
    _rplottrait!(tree.d2, xc, x, y)
  elseif defd1
    _rplottrait!(tree.d1, xc, x, y)
  elseif defd2
    _rplottrait!(tree.d2, xc, x, y)
  end
end



