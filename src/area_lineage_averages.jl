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
          sl += 1.0
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
function sde!(LA   ::Array{Float64,2},
              δX   ::Array{Float64,3}, 
              δY   ::Array{Float64,3},
              wcol ::Array{Array{Int64,1},1},
              m    ::Int64,
              ntip ::Int64)

  @inbounds begin

    for i = Base.OneTo(m), j = wcol[i]
      LA[i,j] = 0.0
      for l = wcol[i]
        l == j && continue
        y = δY[l,j,i]
        iszero(y) && continue
        x = δX[l,j,i]
        LA[i,j] += sign(x) * y * exp(-abs(x))
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

    for i = Base.OneTo(m), k = Base.OneTo(narea), l = wcol[i]

      LD[i,l,k] = 0.0
      sj        = 0.0

      for j = wcol[i]
        j == l && continue
        y          = Float64(Y[i,j,k])
        LD[i,l,k] += abs(δX[j,l,i])*y
        sj        += y
      end
      LD[i,l,k] /= (iszero(sj) ? 1.0 : sj)
    
    end
  end

  return nothing
end





"""
    Xupd_linavg!(δxi  ::Array{Float64,2},
                 lai  ::Array{Float64,1},
                 ldi  ::Array{Float64,2},
                 wci  ::Array{Int64,1},
                 xpi  ::Array{Float64,1},
                 xi   ::Int64,
                 xj   ::Int64,
                 Y    ::Array{Int64,3},
                 δyi  ::Array{Float64,2},
                 narea::Int64)

Re-estimate lineage specific means 
for a node update.
"""
function Xupd_linavg!(δxi  ::Array{Float64,2},
                      lai  ::Array{Float64,1},
                      ldi  ::Array{Float64,2},
                      wci  ::Array{Int64,1},
                      xpi  ::Array{Float64,1},
                      xi   ::Int64,
                      xj   ::Int64,
                      Y    ::SubArray{Int64,2,Array{Int64,3},Tuple{Int64,Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}}},true},
                      δyi  ::SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true},
                      narea::Int64)

  @inbounds begin

    # estimate pairwise distances
    for l = wci, j = wci
      l == j && continue
      δxi[j,l] = xpi[j] - xpi[l]
    end

    # estimate lineage averages
    lai[:] = 0.0
    for l = wci, j = wci
        j == l && continue
        lai[l] += sign(δxi[j,l]) * δyi[j,l] * exp(-abs(δxi[j,l]))
    end

    # estimate lineage sum of distances
    ldi[:] = 0.0
    for k = Base.OneTo(narea), l = wci
      sj = 0.0
      for j = wci
        j == l && continue
        y         = Float64(Y[j,k])
        ldi[l,k] += abs(δxi[j,l])*y
        sj       += y
      end
      ldi[l,k] /= (iszero(sj) ? 1.0 : sj)
    end
  end

  return nothing
end




