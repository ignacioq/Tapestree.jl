#=

insane tree plot

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#


# make this into a recipe


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
  fontfamily      --> font(2, "Helvetica")
  xflip           --> true
  xtickfont       --> font(8, "Helvetica")
  grid            --> :off
  xtick_direction --> :out
  ytick_direction --> :out
  fillcolor       --> plot_color(:orange, 0.3)

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

    Shape(sh0)
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

    Shape(sh1)
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
    _rplottree!(tree::iTgbm, 
                xc  ::Float64, 
                yr  ::UnitRange{Int64},
                zfun::Function,
                x   ::Array{Float64,1}, 
                y   ::Array{Float64,1},
                z   ::Array{Float64,1})

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

  if isdefined(tree, :d1)
    ntip1 = sntn(tree.d1, 0)
    ntip2 = sntn(tree.d2, 0)

    yr1 = yr[1:ntip1]
    yr2 = yr[(ntip1+1):(ntip1+ntip2)]

    xcmpe = xc - e(tree)
    # add vertical lines
    push!(x, xcmpe, xcmpe)
    push!(y, Float64(yr1[1] + yr1[end])*0.5, 
             Float64(yr2[1] + yr2[end])*0.5)
    push!(z, z[end-1])
    push!(z, z[end-2])

    push!(x, NaN)
    push!(y, NaN)
    push!(z, NaN)

    _rplottree!(tree.d1, xcmpe, yr1, zfun, x, y, z)
    _rplottree!(tree.d2, xcmpe, yr2, zfun, x, y, z)
  end

end




"""
    function f(tree::T, zfun::Function)

Recipe for plotting a Type `iTgbm`.
"""
@recipe function f(tree::T, zfun::Function) where {T <: iTgbm}

  x = Float64[]
  y = Float64[]
  z = Float64[]

  th = treeheight(tree)
  nt = sntn(tree, 0)

  _rplottree!(tree, th, 1:nt, zfun, x, y, z)

  # plot defaults
  line_z          --> z
  linecolor       --> :inferno
  legend          --> :none
  colorbar        --> true
  xguide          --> "time"
  fontfamily      --> font(2, "Helvetica")
  xlims           --> (0, th)
  ylims           --> (0, nt+1)
  xflip           --> true
  xtickfont       --> font(8, "Helvetica")
  grid            --> :off
  xtick_direction --> :out
  yticks          --> (nothing)
  yshowaxis       --> false

  return x, y
end




"""
    function f(tree::T, zfun::Function, ϵ::Float64)

Recipe for plotting a Type `iTgbmct` given `ϵ`.
"""
@recipe function f(tree::iTgbmct, zfun::Function, ϵ::Float64)

  x = Float64[]
  y = Float64[]
  z = Float64[]

  th = treeheight(tree)
  nt = sntn(tree, 0)

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
  fontfamily      --> font(2, "Helvetica")
  xlims           --> (0, th)
  ylims           --> (0, nt+1)
  xflip           --> true
  xtickfont       --> font(8, "Helvetica")
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
  xc  -= e(tree)
  push!(x, xc, NaN)
  yc = (yr[1] + yr[end])*0.5
  push!(y, yc, yc, NaN)

  if isdefined(tree, :d1)
    ntip1 = sntn(tree.d1, 0)
    ntip2 = sntn(tree.d2, 0)

    yr1 = yr[1:ntip1]
    yr2 = yr[(ntip1+1):(ntip1+ntip2)]

    # add vertical lines
    push!(x, xc, xc, NaN)
    push!(y, Float64(yr1[1] + yr1[end])*0.5,
             Float64(yr2[1] + yr2[end])*0.5,
             NaN)

    _rplottree!(tree.d1, xc, yr1, x, y)
    _rplottree!(tree.d2, xc, yr2, x, y)
  end

end




"""
    function f(tree::T) where {T <: iTree}
Recipe for plotting a Type `iTree`.
"""
@recipe function f(tree::T) where {T <: iTree}

  x = Float64[]
  y = Float64[]

  th = treeheight(tree)
  nt = sntn(tree, 0)

  _rplottree!(tree, th, 1:nt, x, y)

  # plot defaults
  legend          --> false
  xguide          --> "time"
  fontfamily      --> font(2, "Helvetica")
  seriescolor     --> :black
  xlims           --> (0, th)
  ylims           --> (0, nt+1)
  xflip           --> true
  xtickfont       --> font(8, "Helvetica")
  grid            --> :off
  xtick_direction --> :out
  yticks          --> (nothing)
  yshowaxis       --> false

  return x, y
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
  fontfamily      --> font(2, "Helvetica")
  seriescolor     --> :black
  tickfont        --> font(8, "Helvetica")
  grid            --> :off
  tick_direction  --> :out
  seriestype      --> :steppost
  if maximum(y)>=10
    yaxis         --> :log
  end

  return  x, y
end
