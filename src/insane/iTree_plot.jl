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




#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Tree plot
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#


"""
    f(tree::T;
      labsize    = 8,
      type       = :phylogram,
      showlabels = (T <: Tlabel),
      shownodes  = (false, false, (T <: iTf)),
      shapes     = [:circle, :circle, :square],
      colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
      shsizes    = [0.0, 0.0, 2.0],
      showda     = false,
      col_da     = ["#a9a9a9", :black]) where {T <: iTree}

Recipe for plotting a Type `iTree`. Displays type-specific nodes if `shownodes
= true`. True by default for `sTf` trees to make sampled ancestors visible.
"""
@recipe function f(tree::T;
                   labsize    = 8,
                   type       = :phylogram,
                   showlabels = (T <: Tlabel),
                   shownodes  = (false, false, T <: iTf),
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

  th  = treeheight(tree)
  nts = ntips(tree)

  _rplottree!(tree, th, 0, x, y, z, nodet, xnode, ynode, shownodes)

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
      seriestype         := :scatter
      markershape       --> shapes[nodet]
      markercolor       --> colors[nodet]
      markerstrokecolor --> colors[nodet]
      markersize        --> shsizes[nodet]
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
      ylims             --> (-th*1.05, th*1.15)

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
                     nodet::Array{Int64,1},
                     xnode::Array{Float64,1},
                     ynode::Array{Float64,1},
                     show ::NTuple{3,Bool}) where {T <: iTree}

  xe = xc - e(tree)

  if def1(tree)
    if def2(tree)

      y1, i = _rplottree!(tree.d1, xe, i, x, y, z, nodet, xnode, ynode, show)
      y2, i = _rplottree!(tree.d2, xe, i, x, y, z, nodet, xnode, ynode, show)

      yc = (y1 + y2)*0.5

      # add vertical lines
      push!(x, xe, xe, NaN, xe, xe, NaN)
      push!(y, y1, yc, NaN, yc, y2, NaN)

      z1 = Float64(isfix(tree.d1))
      z2 = Float64(isfix(tree.d2))
      push!(z, z1, z1, NaN, z2, z2, NaN)

      # nodes
      if show[1]
        push!(nodet, 1)
        push!(xnode, xe)
        push!(ynode, yc)
      end
    else
      yc, i = _rplottree!(tree.d1, xe, i, x, y, z, nodet, xnode, ynode, show)

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




#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Tree plot painted by function  f
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#



"""
    f(tree::T,
      f   ::Function;
      type       = :phylogram,
      showlabels = (T <: Tlabel),
      shownodes  = (false, false, (T <: iTf)),
      shapes     = [:circle, :circle, :square],
      colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
      shsizes    = [0.0, 0.0, 3.0],
      simple     = false) where {T <: iT}

Recipe for plotting a Type `iT`.
"""
@recipe function f(tree::T,
                   f   ::Function;
                   type       = :phylogram,
                   showlabels = (T <: Tlabel),
                   shownodes  = (false, false, T <: iTf),
                   shapes     = [:none, :none, :square],
                   colors     = ["#BACBDB", "#DA6A00", "#4D8FC3"],
                   shsizes    = [0.0, 0.0, 3.0],
                   simple     = false) where {T <: iTree}

  x = Float64[]
  y = Float64[]
  z = Float64[]
  nodet = Int64[]   # 0 = speciation, 1 = extinction, 2 = fossilization
  xnode = Float64[]
  ynode = Float64[]

  th  = treeheight(tree)
  nts = ntips(tree)

  _rplottree!(tree, f, th, 0, x, y, z, nodet, xnode, ynode, shownodes, simple)

  ntF = Float64(nts)

  # plot defaults
  fontfamily              --> :Helvetica
  legend                  --> :none
  colorbar                --> true
  colorbar_tickfontfamily --> :Helvetica
  grid                    --> :off
  yshowaxis               --> false
  yticks                  --> nothing
  tickfontfamily          --> :Helvetica
  tickfontsize            --> 8

  if type === :radial
    x, y, z = append_forradial(x, y, z, 50)
    polar_coords!(x, y, 360.0/ntF, th)
    ylims           --> (-th*1.05, th*1.05)
    xlims           --> (-th*1.05, th*1.05)
    xticks          --> nothing
    xshowaxis       --> false
  else
    xlims           --> (-th*0.05, th*1.05)
    ylims           --> (1.0-(0.05*ntF), ntF+(0.05*ntF))
    xguide          --> "time"
    xflip           --> true
    xtick_direction --> :out
    colorbar        --> true
  end

  @series begin
    seriestype  := :path
    line_z     --> z
    linecolor  --> cgrad(:roma, rev = true)

    return x, y
  end

  if any(shownodes)
    @series begin
      seriestype         := :scatter
      markershape       --> shapes[nodet]
      markercolor       --> colors[nodet]
      markerstrokecolor --> colors[nodet]
      markersize        --> shsizes[nodet]
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
                f     ::Function,
                xc    ::Float64,
                i     ::Int64,
                x     ::Array{Float64,1},
                y     ::Array{Float64,1},
                z     ::Array{Float64,1},
                nodet ::Array{Int64,1},
                xnode ::Array{Float64,1},
                ynode ::Array{Float64,1},
                show  ::NTuple{3,Bool},
                simple::Bool) where {T <: iTree}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplottree!(tree  ::T,
                     f     ::Function,
                     xc    ::Float64,
                     i     ::Int64,
                     x     ::Array{Float64,1},
                     y     ::Array{Float64,1},
                     z     ::Array{Float64,1},
                     nodet ::Array{Int64,1},
                     xnode ::Array{Float64,1},
                     ynode ::Array{Float64,1},
                     show  ::NTuple{3,Bool},
                     simple::Bool) where {T <: iTree}

  xe = xc - e(tree)

  if def1(tree)
     if def2(tree)

      y1, i = _rplottree!(tree.d1, f, xe, i, x, y, z, nodet, xnode, ynode, 
        show, simple)
      y2, i = _rplottree!(tree.d2, f, xe, i, x, y, z, nodet, xnode, ynode,
        show, simple)

      yc = (y1 + y2)*0.5
      zc = last(f(tree))

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

      yc, i = _rplottree!(tree.d1, f, xe, i, x, y, z, nodet, xnode, ynode, 
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
    zv = f(tree)
    l  = lastindex(zv)
    push!(z, mean(f(tree)))
    push!(y, yc)
    push!(x, xc)
  else
    # plot function
    zv = copy(f(tree))
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




#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Function f painted by function zf
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    f(f::Function, tree::T; zf = f) where {T <: iTree}

Recipe for plotting function f painted by function zf.
"""
@recipe function f(f::Function, tree::T; zf = f) where {T <: iTree}

  th = treeheight(tree)
  x = Float64[]
  y = Float64[]
  if zf === f
    _rplotf!(tree, th, x, y, f)
    z = y
  else
    z = Float64[]
    _rplotf!(tree, th, x, y, z, f, zf)
  end
  yfilt = filter(x -> !isnan(x), y)
  ymn, ymx = extrema(yfilt)
  pe = 0.02*(ymx - ymn)

  # plot defaults
  legend           --> false
  xguide           --> "time"
  yguide           --> string(f,"(t)")
  xlims            --> (-th*0.05, th*1.05)
  xflip            --> true
  fontfamily       --> :Helvetica
  tickfontfamily   --> :Helvetica
  tickfontsize     --> 8
  grid             --> :off
  xtick_direction  --> :out
  linecolor        --> cgrad(:roma, rev = true)
  line_z           --> z
  ylims            --> (ymn - pe, ymx + pe)
  colorbar         --> true

  return x, y
end





"""
    f(f::Function, tree::Vector{T}; zf = f) where {T <: iTree}

Recipe for plotting function f painted by function zf across multiple trees.
"""
@recipe function f(f::Function, trees::Vector{T}; zf = f) where {T <: iTree}

  th = treeheight(trees[1])
  n  = lastindex(trees)

  x = Float64[]
  y = Float64[]
  if zf === f
    for t in trees
      _rplotf!(t, th, x, y, f)
    end
    z = y
  else
    z = Float64[]
    for t in trees
      _rplotf!(t, th, x, y, z, f, zf)
    end
  end
  yfilt = filter(x -> !isnan(x), y)
  ymn, ymx = extrema(yfilt)
  pe = 0.02*(ymx - ymn)

  # plot defaults
  legend           --> false
  xguide           --> "time"
  yguide           --> string(f,"(t)")
  xlims            --> (-th*0.05, th*1.05)
  xflip            --> true
  fontfamily       --> :Helvetica
  tickfontfamily   --> :Helvetica
  tickfontsize     --> 8
  grid             --> :off
  xtick_direction  --> :out
  line_z           --> z
  linecolor        --> cgrad(:roma, rev = true)
  ylims            --> (ymn - pe, ymx + pe)
  colorbar         --> true
  seriesalpha      --> min(1.0, 10.0/Float64(n))

  return x, y
end




"""
    _rplotf!(tree::T,
             xc  ::Float64,
             x   ::Array{Float64,1},
             y   ::Array{Float64,1},
             f   ::Function) where {T <: Txs}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplotf!(tree::T,
                  xc  ::Float64,
                  x   ::Array{Float64,1},
                  y   ::Array{Float64,1},
                  f   ::Function) where {T <: iTree}

  # tree δt and nsδt
  δt = dt(tree)

  # add horizontal lines
  l = lastindex(f(tree))
  @simd for i in Base.OneTo(l-1)
    push!(x, xc - Float64(i-1)*δt)
  end
  push!(x, xc - (Float64(l-2)*δt + fdt(tree)), NaN)

  append!(y, f(tree))
  push!(y, NaN)

  xc -= e(tree)

  if def1(tree)
    _rplotf!(tree.d1, xc, x, y, f)
    if def2(tree)
      _rplotf!(tree.d2, xc, x, y, f)
    end
  end
end




"""
    _rplotf!(tree::T,
                 xc  ::Float64,
                 x   ::Array{Float64,1},
                 y   ::Array{Float64,1},
                 z   ::Array{Float64,1}
                 yf  ::Function,
                 zf  ::Function) where {T <: Txs}

Returns `x` and `y` coordinates in order to plot a tree of type `iTree`.
"""
function _rplotf!(tree::T,
                  xc  ::Float64,
                  x   ::Array{Float64,1},
                  y   ::Array{Float64,1},
                  z   ::Array{Float64,1},
                  yf  ::Function,
                  zf  ::Function) where {T <: iTree}

  # tree δt and nsδt
  δt = dt(tree)

  # add horizontal lines
  l   = lastindex(yf(tree))
  @simd for i in Base.OneTo(l-1)
    push!(x, xc - Float64(i-1)*δt)
  end
  push!(x, xc - (Float64(l-2)*δt + fdt(tree)), NaN)

  append!(y, yf(tree))
  push!(y, NaN)
  append!(z, zf(tree))
  push!(z, NaN)

  xc -= e(tree)

  if def1(tree)
    _rplotf!(tree.d1, xc, x, y, z, yf, zf)
    if def2(tree)
      _rplotf!(tree.d2, xc, x, y, z, yf, zf)
    end
  end
end





"""
    function f(tree::T; type::Symbol = :trait)

Recipe for plotting punctuated equilibrium trees.
"""
@recipe function f(tree::sTpe)

  x = Float64[]
  y = Float64[]

  th = treeheight(tree)

  _rplottrait!(tree, th, xi(tree), x, y)

  yfilt = filter(x -> !isnan(x), y)
  ymn = minimum(yfilt)
  ymx = maximum(yfilt)
  rng = ymx - ymn

  # plot defaults
  legend          --> false
  xguide          --> "time"
  yguide          --> "trait"
  seriescolor     --> :purple
  xlims           --> (-th*0.05, th*1.05)
  ylims           --> (ymn - 0.05*rng, ymx + 0.05*rng)
  xflip           --> true
  fontfamily      --> :Helvetica
  tickfontfamily  --> :Helvetica
  tickfontsize    --> 8
  grid            --> :off
  xtick_direction --> :out

  return x, y

end



"""
    _rplottrait!(tree::sTpe,
                 xc  ::Float64,
                 x   ::Array{Float64,1},
                 y   ::Array{Float64,1})

Returns `x` and `y` coordinates in order to plot a tree of type `sTpe`.
"""
function _rplottrait!(tree::sTpe,
                      xc  ::Float64,
                      yc  ::Float64,
                      x   ::Array{Float64,1},
                      y   ::Array{Float64,1})

  xii = xi(tree)

  # add vertical lines
  if yc != xii
    push!(x, xc, xc,  NaN)
    push!(y, yc, xii, NaN)
  end

  ei = e(tree)

  # add horizontal lines
  xfi = xf(tree)
  push!(x, xc, xc - ei, NaN)
  push!(y, xii, xfi, NaN)

  xc -= ei

  if def1(tree)
    _rplottrait!(tree.d1, xc, xfi, x, y)
    if def2(tree)
      _rplottrait!(tree.d2, xc, xfi, x, y)
    end
  end
end






#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Function f aggregated by af through time
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    f(f   ::Function,
      δt  ::Float64,
      tree::T;
      af = mean,
      q0 = Float64[],
      q1 = Float64[],
      q2 = Float64[]) where {T <: iTree}

Recipe for plotting values given by `f` through time for a `iT`.
"""
@recipe function f(f   ::Function,
                   δt  ::Float64,
                   tree::T;
                   af = mean,
                   q0 = Float64[],
                   q1 = Float64[],
                   q2 = Float64[]) where {T <: iTree}

  # prepare data
  ts, r = time_rate(tree, f, δt)
  lts = lastindex(ts)
  m   = map(af, r)

  # common shape plot defaults
  legend          --> :none
  xguide          --> "time"
  yguide          --> string(f,"(t)")
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
               δt   ::Float64;
               af  = mean,
               vaf = x -> quantile(x, 0.5),
               q0  = [0.025, 0.975],
               q1  = [0.25,  0.75],
               q2  = Float64[]) where {T <: iTree}

Recipe for plotting values given by `f` through time for a `iT`.
"""
@recipe function f(f    ::Function,
                   δt   ::Float64,
                   trees::Vector{T};
                   af  = mean,
                   vaf = x -> quantile(x, 0.5),
                   q0  = [0.025, 0.975],
                   q1  = [0.25,  0.75],
                   q2  = Float64[]) where {T <: iTree}

  ntrees = lastindex(trees)
  riv = Vector{Float64}[]

  lts = 0
  ts  = Float64[]
  for t in trees

    # tree extracting function
    tsi, ri = time_rate(t, f, δt)

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

  M = fill(NaN, lts)
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
    M[i] = vaf(qi)
  end

  # common shape plot defaults
  legend          --> :none
  xguide          --> "time"
  yguide          --> string(f,"(t)")
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





#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Diversity through time
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#



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






