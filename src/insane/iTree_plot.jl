#=

insane tree plot

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#





"""
    rplottree!(tree::iTgbm, 
               xc  ::Float64, 
               yr  ::UnitRange{Int64},
               zfun::Function,
               x   ::Array{Float64,1}, 
               y   ::Array{Float64,1},
               z   ::Array{Float64,1})

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function rplottree!(tree::T, 
                    xc  ::Float64, 
                    yr  ::UnitRange{Int64},
                    zfun::Function,
                    x   ::Array{Float64,1}, 
                    y   ::Array{Float64,1},
                    z   ::Array{Float64,1}) where {T <: iTgbm}

  # add horizontal lines
  xv  = ts(tree)
  yc  = Float64(yr[1] + yr[end])*0.5
  zv  = exp.(zfun(tree))
  lts = lastindex(xv)
  for i in Base.OneTo(lts)
    push!(x, xc - xv[i])
    push!(y, yc)
    push!(z, zv[i])
  end

  # push!(x, NaN)
  # push!(y, NaN)
  # push!(z, NaN)

  if !istip(tree)
    ntip1 = sntn(tree.d1)
    ntip2 = sntn(tree.d2)

    yr1 = yr[1:ntip1]
    yr2 = yr[(ntip1+1):(ntip1+ntip2)]

    xcmpe = xc - pe(tree)
    # add vertical lines
    push!(x, xcmpe, xcmpe)
    push!(y, Float64(yr1[1] + yr1[end])*0.5, 
             Float64(yr2[1] + yr2[end])*0.5)
    push!(z, z[end-1])
    push!(z, z[end-2])

    # push!(x, NaN)
    # push!(y, NaN)
    # push!(z, NaN)

    rplottree!(tree.d1, xcmpe, yr1, zfun, x, y, z)
    rplottree!(tree.d2, xcmpe, yr2, zfun, x, y, z)
  end

end




"""
    function f(tree::iTgbm, cfun::Function)
Recipe for plotting a Type `iTgbm`.
"""
@recipe function f(tree::T, zfun::Function) where {T <: iTgbm}

  x = Float64[]
  y = Float64[]
  z = Float64[]
  rplottree!(tree, treeheight(tree), 1:sntn(tree), zfun, x, y, z)

  # plot defaults
  line_z          --> z
  linecolor       --> :inferno
  legend          --> :none
  colorbar        --> true
  xguide          --> "time"
  fontfamily      --> font(2, "Helvetica")
  xlims           --> (0, treeheight(tree))
  ylims           --> (0, sntn(tree)+1)
  xflip           --> true
  xtickfont       --> font(8, "Helvetica")
  grid            --> :off
  xtick_direction --> :out
  yticks          --> (nothing)
  yshowaxis       --> false

  # nan_inds = findall(x -> isnan(x),z)
  # for i in eachindex(nan_inds)
  #     start = 1 + (i == 1 ? 0 : nan_inds[i - 1])
  #     stop = nan_inds[i] - 1
  #     @series begin
  #         line_z --> z[start:stop]
  #         x[start:stop], y[start:stop]
  #     end
  # end
  return x, y
end




"""
    rplottree!(tree::T, 
              xc  ::Float64, 
              yr  ::UnitRange{Int64},
              x   ::Array{Float64,1}, 
              y   ::Array{Float64,1}) where {T <: iTree}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function rplottree!(tree::T, 
                    xc  ::Float64, 
                    yr  ::UnitRange{Int64},
                    x   ::Array{Float64,1}, 
                    y   ::Array{Float64,1}) where {T <: iTree}

  # add horizontal lines
  push!(x, xc)
  xc  -= pe(tree)
  push!(x, xc, NaN)
  yc = (yr[1] + yr[end])*0.5
  push!(y, yc, yc, NaN)

  if !istip(tree)
    ntip1 = sntn(tree.d1)
    ntip2 = sntn(tree.d2)

    yr1 = yr[1:ntip1]
    yr2 = yr[(ntip1+1):(ntip1+ntip2)]

    # add vertical lines
    push!(x, xc, xc, NaN)
    push!(y, Float64(yr1[1] + yr1[end])*0.5,
             Float64(yr2[1] + yr2[end])*0.5,
             NaN)

    rplottree!(tree.d1, xc, yr1, x, y)
    rplottree!(tree.d2, xc, yr2, x, y)
  end

end




"""
    function f(tree::T) where {T <: iTree}
Recipe for plotting a Type `iTree`.
"""
@recipe function f(tree::T) where {T <: iTree}

  x = Float64[]
  y = Float64[]
  rplottree!(tree, treeheight(tree), 1:sntn(tree), x, y)

  # plot defaults
  legend          --> false
  xguide          --> "time"
  fontfamily      --> font(2, "Helvetica")
  seriescolor     --> :black
  xlims           --> (0, treeheight(tree))
  ylims           --> (0, sntn(tree)+1)
  xflip           --> true
  xtickfont       --> font(8, "Helvetica")
  grid            --> :off
  xtick_direction --> :out
  yticks          --> (nothing)
  yshowaxis       --> false

  return x, y
end


