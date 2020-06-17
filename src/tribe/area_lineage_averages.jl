#=

functions to estimate area and lineage trait averages 
and lineage specific differences

Ignacio Quintero Mächler

t(-_-t)

May 15 2017

=#





"""
    deltaX!(δX   ::Array{Float64,3}, 
            X    ::Array{Float64,2},
            m    ::Int64,
            ntip ::Int64,
            narea::Int64)

Estimate pairwise trait and range distances between all lineages
conditional on a given `δY`.
"""
function deltaX!(δX   ::Array{Float64,3}, 
                 X    ::Array{Float64,2},
                 wcol ::Array{Array{Int64,1},1},
                 m    ::Int64,
                 ntip ::Int64,
                 narea::Int64)

  @inbounds begin

    for i = Base.OneTo(m), l = wcol[i], j = wcol[i]
      l == j && continue
      δX[j,l,i] = X[i,j] - X[i,l]
    end

  end

  return nothing
end





"""
    deltaY!(δY   ::Array{Float64,3}, 
            Y    ::Array{Int64,3},
            m    ::Int64,
            ntip ::Int64,
            narea::Int64)

Estimate pairwise trait and range distances between all lineages
conditional on a given `δY`.
"""
function deltaY!(δY   ::Array{Float64,3}, 
                 Y    ::Array{Int64,3},
                 wcol ::Array{Array{Int64,1},1},
                 m    ::Int64,
                 ntip ::Int64,
                 narea::Int64)

  @inbounds begin

    for i = Base.OneTo(m), l = wcol[i], j = wcol[i]

      j == l && continue

      sl        = 0.0
      δY[l,j,i] = 0.0
      @simd for k = Base.OneTo(narea)
        if Y[i,j,k] == 1
          sl        += 1.0
          δY[l,j,i] += Float64(Y[i,l,k])
        end
      end
      δY[l,j,i] /= sl
    end

  end

  return nothing
end





"""
    deltaXY!(δX   ::Array{Float64,3}, 
             δY   ::Array{Float64,3},
             X    ::Array{Float64,2},
             Y    ::Array{Int64,3},
             m    ::Int64,
             ntip ::Int64,
             narea::Int64)

Estimate pairwise trait and range distances between all lineages.
"""
function deltaXY!(δX   ::Array{Float64,3}, 
                  δY   ::Array{Float64,3},
                  X    ::Array{Float64,2},
                  Y    ::Array{Int64,3},
                  wcol ::Array{Array{Int64,1},1},
                  m    ::Int64,
                  ntip ::Int64,
                  narea::Int64)

  @inbounds begin

    for i = Base.OneTo(m), l = wcol[i], j = wcol[i]

      j == l && continue

      # X differences
      δX[j,l,i] = X[i,j] - X[i,l]

      # Y area overlap
      sl        = 0.0
      δY[l,j,i] = 0.0
      @simd for k = Base.OneTo(narea)
        if Y[i,j,k] == 1
          sl        += 1.0
          δY[l,j,i] += Float64(Y[i,l,k])
        end
      end
      δY[l,j,i] /= sl
    end

  end

  return nothing
end





"""
    sde!(LA   ::Array{Float64,2},
         δX   ::Array{Float64,3}, 
         δY   ::Array{Float64,3},
         m    ::Int64,
         ntip ::Int64)

Estimate the lineage averages for the SDE of trait evolution
"""
function sde!(LAp   ::Array{Float64,2},
              LAn   ::Array{Float64,2},
              δX   ::Array{Float64,3}, 
              δY   ::Array{Float64,3},
              wcol ::Array{Array{Int64,1},1},
              m    ::Int64,
              ntip ::Int64)

  @inbounds begin

    for i = Base.OneTo(m), j = wcol[i]
      LAp[i,j] = 0.0
      LAn[i,j] = 0.0
      for l = wcol[i]
        l == j && continue
        y = δY[l,j,i]
        iszero(y) && continue
        x = δX[l,j,i]
        LAp[i,j] += x*y
        LAn[i,j] += sign(x) * y * exp(-abs(x))
      end
    end
  end

  return nothing
end





"""
    lindiff!(LD   ::Array{Float64,3},
             δX   ::Array{Float64,3},
             Y    ::Array{Int64,3},
             m    ::Int64,
             ntip ::Int64,
             narea::Int64)

Estimate area-specific distance averages 
for colonization and extinction rates.
"""
function lindiff!(LD   ::Array{Float64,3},
                  δX   ::Array{Float64,3},
                  Y    ::Array{Int64,3},
                  wcol ::Array{Array{Int64,1},1},
                  m    ::Int64,
                  ntip ::Int64,
                  narea::Int64)
  @inbounds begin

    for k = Base.OneTo(narea), i = Base.OneTo(m), l = wcol[i]

      xmin = Inf
      for j = wcol[i]
        j == l && continue
        if Y[i,j,k] == 1
          x = abs(δX[j,l,i])
          iszero(x) && continue
          xmin = x < xmin ? x : xmin
        end
      end
      LD[i,l,k] = isinf(xmin) ? 0.0 : xmin
    end
  end

  return nothing
end






"""
    Xupd_linavg!(δxi  ::Array{Float64,2},
                 lapi ::Array{Float64,1},
                 lani ::Array{Float64,1},
                 ldi  ::Array{Float64,2},
                 wci  ::Array{Int64,1},
                 xpi  ::Array{Float64,1},
                 xi   ::Int64,
                 xj   ::Int64,
                 Y    ::Array{Int64,3},
                 δyi  ::Array{Float64,2},
                 narea::Int64)

Re-estimate lineage specific means 
for a node update when `ωx >= 0.0`
"""
function Xupd_linavg!(δxi  ::Array{Float64,2},
                      lapi ::Array{Float64,1},
                      lani ::Array{Float64,1},
                      ldi  ::Array{Float64,2},
                      wcol ::Array{Array{Int64,1},1},
                      xpi  ::Array{Float64,1},
                      xi   ::Int64,
                      xj   ::Int64,
                      Y    ::Array{Int64,3},
                      δY   ::Array{Float64,3},
                      narea::Int64)

  @inbounds begin

    # estimate pairwise distances
    wci = wcol[xi]
    for l = wci, j = wci
      l == j && continue
      δxi[j,l] = xpi[j] - xpi[l]
    end

    # estimate lineage averages
    for l = wci 
      lapi[l] = 0.0
      lani[l] = 0.0
      for j = wci
        j == l && continue
        y = δY[j,l,xi]
        iszero(y) && continue
        x = δxi[j,l]
        lapi[l] += x * y
        lani[l] += sign(x) * y * exp(-abs(x))
      end
    end

    # estimate lineage sum of distances
    for k = Base.OneTo(narea), l = wci
      xmin = Inf
      for j = wci
        j == l && continue
        if Y[xi,j,k] == 1
          x = abs(δxi[j,l])
          iszero(x) && continue
          xmin = x < xmin ? x : xmin
        end
      end
      ldi[l,k] = xmin == Inf ? 0.0 : xmin
    end
  end

  return nothing
end


