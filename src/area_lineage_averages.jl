#=

Functions to estimate area and lineage trait averages 
and lineage specific differences

Ignacio Quintero Mächler

t(-_-t)

May 15 2017

=#





"""
    area_lineage_means!(AA::Array{Float64,2}, LA::Array{Float64,2}, X::Array{Float64,2}, Y::Array{Int64,3}, wcol::Array{Array{Int64,1},1}, m::Int64)

Estimate area means according to presence
absence of species and linage means according
to area averages.
"""
function area_lineage_means!(AA  ::Array{Float64,2}, 
                             LA   ::Array{Float64,2},
                             AO   ::Array{Int64,2},
                             X    ::Array{Float64,2}, 
                             Y    ::Array{Int64,3}, 
                             wcol ::Array{Array{Int64,1},1},
                             m    ::Int64,
                             narea::Int64)

  @inbounds begin

    for k in Base.OneTo(m)

      # area averages
      for j in Base.OneTo(narea)
        AA[k,j] = 0.0::Float64
        sumY    = 0.0::Float64
        AO[k,j] = 0::Int64
        for i in wcol[k]
          if Y[k,i,j] == 1
            AA[k,j] += X[k,i]::Float64
            sumY    += 1.0::Float64
            AO[k,j]  = 1::Int64
          end
        end

        AA[k,j] /= (sumY == 0.0 ? 1.0 : sumY)::Float64
      end

      # lineage average
      for i = wcol[k]
        LA[k,i] = 0.0::Float64
        sden    = 0.0::Float64
        for j = Base.OneTo(narea) 
          if Y[k,i,j] == 1
            LA[k,i] += AA[k,j]::Float64
            sden    += 1.0::Float64
          end
        end
        
        LA[k,i] /= sden::Float64
      end
    end

  end

  nothing
end





"""
    linarea_diff!(LD::Array{Float64,3}, X::Array{Float64,2}, AA::Array{Float64,2}, narea::Int64, ntip::Int64, m::Int64)

Create multi-dimensional array with 
lineage and area averages differences in X.
"""
function linarea_diff!(LD   ::Array{Float64,3},
                       X    ::Array{Float64,2},
                       AA   ::Array{Float64,2},
                       AO   ::Array{Int64,2},
                       narea::Int64,
                       ntip ::Int64,
                       m    ::Int64)
  @inbounds begin

    for j = Base.OneTo(narea), n = Base.OneTo(ntip), i = Base.OneTo(m)

      if AO[i,j] == 0
        setindex!(LD, 0.0, i, n, j)
      else
        setindex!(LD, abs(X[i,n] - AA[i,j]), i, n, j)
      end
    end

  end

  nothing
end





"""
    linarea_branch_avg!(avg_Δx ::Array{Float64,1}, LD::Array{Float64,3}, bridx_a::Array{Array{Array{Int64,1},1},1}, narea::Int64, nedge::Int64)

Estimate the branch average of lineage differences in each specific area.
"""
function linarea_branch_avg!(avg_Δx ::Array{Float64,2},
                             LD     ::Array{Float64,3},
                             bridx_a::Array{Array{UnitRange{Int64},1},1},
                             narea  ::Int64,
                             nedge  ::Int64)
  @inbounds begin

    for j = Base.OneTo(narea), i = Base.OneTo(nedge - 1)
      setindex!(avg_Δx, mean(LD[bridx_a[j][i]]), i, j)
    end

  end

  nothing
end





"""
    Xupd_linavg(k::Int64, wck::Array{Int64,1}, X::Array{Float64,2}, Y::Array{Int64,3}, narea::Int64)

Re-estimate lineage specific means 
for a branch update.
"""
function Xupd_linavg(k    ::Int64, 
                     wck  ::Array{Int64,1},
                     X    ::Array{Float64,2},
                     Y    ::Array{Int64,3},
                     narea::Int64)

  Sk, Sx = symp_traits(X, Y, wck, k)

  # new area averages
  aa, ao = area_averages(Sx, Sk, narea)

  # new lineage averages
  la = lineage_averages(Sk, aa)

  # new lineage differences
  ld = linarea_difference(k, X, aa, ao, wck, narea)

  return aa, la, ld
end




"""
    symp_traits(X::Array{Float64,2}, Y::Array{Int64,3}, wck::Array{Int64,1}, k::Int64)

Return the trait per occupied area per lineage.
"""
function symp_traits(X  ::Array{Float64,2}, 
                     Y  ::Array{Int64,3}, 
                     wck::Array{Int64,1}, 
                     k  ::Int64)

  @inbounds begin
    
    const Sk  = Y[k,wck,:]::Array{Int64,2}
    const xx  = X[k,wck]::Array{Float64,1}
    
    const nrow, ncol = size(Sk)
    const Sx = zeros(Float64,nrow,ncol)

    for j = Base.OneTo(ncol), i = Base.OneTo(nrow)
      if Sk[i,j] == 1
        Sx[i,j] = xx[i]::Float64
      end
    end
  
  end
  
  return (Sk, Sx)::Tuple{Array{Int64,2}, Array{Float64,2}}
end





"""
    area_averages(Sx::Array{Float64,2}, Sk::Array{Int64,2}, narea::Int64)

Estimate area averages based on per area per
lineage values.
"""
function area_averages(Sx   ::Array{Float64,2}, 
                       Sk   ::Array{Int64,2}, 
                       narea::Int64)
  @inbounds begin

    const Sxind = indices(Sx,1)
    const aa    = zeros(Float64,narea)
    const ao    = zeros(Int64,narea)

    for j in Base.OneTo(narea)
      aa[j] = 0.0::Float64
      sumY  = 0.0::Float64
      for i in Sxind
        if Sk[i,j] == 1
          aa[j] += Sx[i,j]::Float64
          sumY  += 1.0::Float64
          ao[j]  = 1::Int64
        end
      end

      aa[j] /= (sumY == 0.0 ? 1.0 : sumY)::Float64
    end

  end

  return (aa, ao)::Tuple{Array{Float64,1},Array{Int64,1}}
end





"""
    lineage_averages(Sk::Array{Int64,2}, AA::Vector{Float64})

Estimate lineage specific competition averages 
based on sympatry and area averages.
"""
function lineage_averages(Sk::Array{Int64,2}, 
                          AA::Vector{Float64})
  @inbounds begin

    const nrow, ncol = size(Sk)::Tuple{Int64,Int64}
    const la = zeros(Float64,nrow)
        
    for i = Base.OneTo(nrow)
      la[i] = 0.0
      sden  = 0.0
      for j = Base.OneTo(ncol)
        if Sk[i,j] == 1
          la[i] += AA[j]
          sden  += 1.0
        end
      end
      la[i] /= sden::Float64
    end

  end

  return la::Array{Float64,1}
end




"""
    linarea_difference(k::Int64, X::Array{Float64,2}, AA::Vector{Float64}, wck ::Array{Int64,1}, narea::Int64)

Create a multi-dimensional array with 
lineage and area averages differences in X.
"""
function linarea_difference(k    ::Int64,
                            X    ::Array{Float64,2},
                            aa   ::Array{Float64,1},
                            ao   ::Array{Int64,1},
                            wck  ::Array{Int64,1},
                            narea::Int64)
  @inbounds begin

    const nsp = endof(wck)::Int64
    const ld = Array{Float64}(nsp, narea)

    for j = Base.OneTo(narea), i = Base.OneTo(nsp)
      if ao[j] == 0
        setindex!(ld, 0.0, i, j)
      else
        setindex!(ld, abs(X[k,wck[i]] - aa[j]), i, j)
      end
    end

  end

  return ld::Array{Float64,2}
end


