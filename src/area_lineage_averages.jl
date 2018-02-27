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
function area_lineage_means!(AA   ::Array{Float64,2}, 
                             LA   ::Array{Float64,2},
                             AO   ::Array{Int64,2},
                             X    ::Array{Float64,2}, 
                             Y    ::Array{Int64,3}, 
                             wcol ::Array{Array{Int64,1},1},
                             m    ::Int64,
                             narea::Int64)

  @inbounds begin

    for i = Base.OneTo(m)

      # area averages
      for k = Base.OneTo(narea)
        AA[i,k] = 0.0::Float64
        sumY    = 0.0::Float64
        AO[i,k] = 0::Int64
        @simd for j = wcol[i]::Array{Int64,1}
          if Y[i,j,k]::Int64 == 1
            AA[i,k] += X[i,j]::Float64
            sumY    += 1.0::Float64
            AO[i,k]  = 1::Int64
          end
        end
        if sumY != 0.0
          AA[i,k] /= sumY::Float64
        end
      end

      # lineage average
      for j = wcol[i]::Array{Int64,1}
        LA[i,j] = 0.0::Float64
        sumY    = 0.0::Float64
        @simd for k = Base.OneTo(narea) 
          if Y[i,j,k]::Int64 == 1
            LA[i,j] += AA[i,k]::Float64
            sumY    += 1.0::Float64
          end
        end
        LA[i,j] /= sumY::Float64
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

    for k = Base.OneTo(narea), j = Base.OneTo(ntip), i = Base.OneTo(m)
      if isnan(X[i,j])
        continue
      end
      setindex!(LD, 
                (AO[i,k] == 0 ? 0.0 : abs(X[i,j] - AA[i,k]))::Float64, 
                i, j, k)
    end

  end
  nothing
end





"""
    make_la_branch_avg(bridx_a, lY::Int64, narea::Int64, nedge::Int64)

Make function to estimate the branch average of lineage differences in each specific area.
"""
function make_la_branch_avg(bridx_a::Array{Array{UnitRange{Int64},1},1},
                            lY     ::Int64,
                            m      ::Int64,
                            narea  ::Int64,
                            nedge  ::Int64)

  # create branch indexes removing tips
  const bridx_a_nt = deepcopy(bridx_a)
  for k=Base.OneTo(narea), j=Base.OneTo(nedge)
    sdi = setdiff(bridx_a[k][j], m:m:lY)
    bridx_a_nt[k][j] = sdi[1]:sdi[end]
  end

  function f(avg_Δx::Array{Float64,2},
             LD    ::Array{Float64,3})

    @inbounds begin
    
      for k = Base.OneTo(narea), i = Base.OneTo(nedge - 1)
        setindex!(avg_Δx, mean(LD[bridx_a_nt[k][i]]), i, k)
      end

    end

    nothing
  end

  return f
end




"""
    Xupd_linavg(k::Int64, wci::Array{Int64,1}, X::Array{Float64,2}, Y::Array{Int64,3}, narea::Int64)

Re-estimate lineage specific means 
for a branch update.
"""
function Xupd_linavg!(aa   ::Array{Float64,1},
                      la   ::Array{Float64,1},
                      ld   ::Array{Float64,2},
                      ao   ::Array{Int64,2},
                      i    ::Int64, 
                      wci  ::Array{Int64,1},
                      xi   ::Array{Float64,1},
                      Y    ::Array{Int64,3},
                      narea::Int64)
  @inbounds begin

    for k = Base.OneTo(narea)
      sumY  = 0.0::Float64
      aa[k] = 0.0::Float64
      @simd for j = wci
        if Y[i,j,k] == 1
          aa[k] += xi[j]::Float64
          sumY  += 1.0
        end
      end
      if sumY != 0.0
        aa[k] /= sumY::Float64
      end
    end

    fill!(la, NaN)
    for j = wci
      la[j] = 0.0
      sumY  = 0.0::Float64
      @simd for k = Base.OneTo(narea)
        if Y[i,j,k]::Int64 == 1
          la[j] += aa[k]::Float64
          sumY  += 1.0
        end
      end
      la[j] /= sumY::Float64
    end

    fill!(ld, NaN)
    for k = Base.OneTo(narea), j = wci
      setindex!(ld, (ao[i,k] == 0 ? 0.0 : abs(xi[j] - aa[k])), j, k)
    end

  end

  return nothing 
end

