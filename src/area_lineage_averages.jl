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
                             LA  ::Array{Float64,2},
                             X   ::Array{Float64,2}, 
                             Y   ::Array{Int64,3}, 
                             wcol::Array{Array{Int64,1},1},
                             m   ::Int64)

  @inbounds begin

    for k in Base.OneTo(m) 
      
      wck = wcol[k]
      Sk  = Y[k,wck,:]
      xx  = X[k,wck]
    
      nrow, ncol = size(Sk)
  
      Sx = Array{Float64}(nrow,ncol)
      for j = Base.OneTo(ncol), i = Base.OneTo(nrow)
        Sx[i,j] = Sk[i,j] * xx[i]
      end
      
      # area averages
      for j in Base.OneTo(ncol)
        sumX = 0.0
        sumY = 0

        for i in Base.OneTo(nrow)
          sumX += Sx[i,j]
          sumY += Sk[i,j]
        end

        AA[k,j] = sumX/(sumY == 0 ? 1 : sumY) 
      end

      # lineage average
      for i = Base.OneTo(nrow)
        snum = 0.0
        sden = 0

        for j = Base.OneTo(ncol) 
          snum += Sk[i,j] * AA[k,j]
          sden += Sk[i,j]
        end
        
        LA[k,wck[i]] = snum/sden
      end
    end

  end
end




"""
    linarea_diff!(LD   ::Array{Float64,3}, X::Array{Float64,2}, AA::Array{Float64,2}, narea::Int64, ntip::Int64, m::Int64)

Create multi-dimensional array with 
lineage and area averages differences in X.
"""
function linarea_diff!(LD   ::Array{Float64,3},
                       X    ::Array{Float64,2},
                       AA   ::Array{Float64,2},
                       narea::Int64,
                       ntip ::Int64,
                       m    ::Int64)
  @inbounds begin

    for j = Base.OneTo(narea), n = Base.OneTo(ntip), i = Base.OneTo(m)
      if AA[i,j] == 0.0
        continue
      end
      setindex!(LD, abs(X[i,n] - AA[i,j]), i, n, j)
    end

  end

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
  aa = area_averages(Sx, Sk, narea)

  # new lineage averages
  la = lineage_averages(Sk, aa)

  # new lineage differences
  ld = linarea_difference(k, X, aa, wck, narea)

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
    
    Sk  = Y[k,wck,:]
    xx  = X[k,wck]
    
    nrow, ncol = size(Sk)
    Sx = Array{Float64}(nrow,ncol)

    for j = Base.OneTo(ncol), i = Base.OneTo(nrow)
      Sx[i,j] = Sk[i,j] * xx[i]
    end
  
  end
  
  Sk, Sx
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

    Sxind = indices(Sx,1)
    aa    = zeros(narea)

    for j in Base.OneTo(narea)
      sumX::Float64 = 0.0
      sumY::Int64   = 0

      for i in Sxind
        sumX += Sx[i,j]
        sumY += Sk[i,j]
      end

      aa[j] = sumX/(sumY == 0 ? 1 : sumY)
    end

  end

  aa
end





"""
    lineage_averages(Sk::Array{Int64,2}, AA::Vector{Float64})

Estimate lineage specific competition averages 
based on sympatry and area averages.
"""
function lineage_averages(Sk::Array{Int64,2}, 
                          AA::Vector{Float64})
  @inbounds begin

    nrow, ncol = size(Sk)
    la = Vector{Float64}(nrow)
        
    for i = Base.OneTo(nrow)
      snum = 0.0
      sden = 0
      for j = Base.OneTo(ncol) 
        snum += Sk[i,j] * AA[j]
        sden += Sk[i,j]
      end
      la[i] = snum/sden
    end

  end

  la
end




"""
    linarea_difference(k::Int64, X::Array{Float64,2}, AA::Vector{Float64}, wck ::Array{Int64,1}, narea::Int64)

Create a multi-dimensional array with 
lineage and area averages differences in X.
"""
function linarea_difference(k    ::Int64,
                            X    ::Array{Float64,2},
                            AA   ::Vector{Float64},
                            wck  ::Array{Int64,1},
                            narea::Int64)
  @inbounds begin

    nsp::Int64 = endof(wck)

    ld = Array{Float64}(nsp, narea)

    for j = Base.OneTo(narea), i = Base.OneTo(nsp)
      if AA[j] == 0.0
        setindex!(ld, 0.0, i, j)
      else
        setindex!(ld, abs(X[k,wck[i]] - AA[j]), i, j)
      end
    end

  end

  ld
end


