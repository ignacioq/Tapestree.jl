#=

insane tree plot

Ignacio Quintero Mächler

t(-_-t)

Created 07 07 2020
=#




"""
    rplottree(tree::T, 
              xc  ::Float64, 
              yr  ::UnitRange{Int64},
              x   ::Array{Float64,1}, 
              y   ::Array{Float64,1}) where {T <: iTree}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function rplottree(tree::T, 
                   xc  ::Float64, 
                   yr  ::UnitRange{Int64},
                   x   ::Array{Float64,1}, 
                   y   ::Array{Float64,1}) where {T <: iTree}

  # add horizontal lines
  push!(x, xc)
  xc  -= pe(tree)
  push!(x, xc)
  yc = (yr[1] + yr[end])/2.0
  push!(y, yc, yc)
  push!(x, NaN)
  push!(y, NaN)

  if !istip(tree)
    ntip1 = sntn(tree.d1)
    ntip2 = sntn(tree.d2)

    yr1 = yr[1:ntip1]
    yr2 = yr[(ntip1+1):(ntip1+ntip2)]

    # add vertical lines
    push!(x, xc, xc)
    push!(y, (yr1[1] + yr1[end])/2.0)
    push!(y, (yr2[1] + yr2[end])/2.0)
    push!(x, NaN)
    push!(y, NaN)

    rplottree(tree.d1, xc, yr1, x, y)
    rplottree(tree.d2, xc, yr2, x, y)
  end

end




"""
    function f(tree::T) where {T <: iTree}

Recipe for plotting a Type `iTree`.
"""
@recipe function f(tree::T) where {T <: iTree}

  x = Float64[]
  y = Float64[]
  rplottree(tree, treeheight(tree), 1:sntn(tree), x, y)
  
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
  yshowaxis       --> false

  return x, y
end


