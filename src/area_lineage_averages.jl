#=

functions to estimate area and lineage trait averages 
and lineage specific differences

Ignacio Quintero Mächler

t(-_-t)

May 15 2017

=#





"""
    deltaXY!(ΔX   ::Array{Float64,3}, 
             ΔY   ::Array{Float64,3},
             X    ::Array{Float64,2},
             Y    ::Array{Int64,3},
             m    ::Int64,
             ntip ::Int64,
             narea::Int64)

Estimate pairwise trait and range distances between all lineages.
"""
function deltaXY!(ΔX   ::Array{Float64,3}, 
                  ΔY   ::Array{Float64,3},
                  X    ::Array{Float64,2},
                  Y    ::Array{Int64,3},
                  m    ::Int64,
                  ntip ::Int64,
                  narea::Int64)

  @inbounds begin

    for i = Base.OneTo(m), l = Base.OneTo(ntip)

      isnan(X[i,l]) && continue

      for j = Base.OneTo(ntip)

        isnan(X[i,j]) && continue

        if j == l 
          ΔX[j,l,i] = 0.0
          ΔY[j,l,i] = 1.0
        else
          # X differences
          ΔX[j,l,i] = X[i,j] - X[i,l]

          # Y area overlap
          sl        = 0.0
          ΔY[l,j,i] = 0.0
          @simd for k = Base.OneTo(narea)
            if Y[i,j,k] == 1
              sl += 1.0
              ΔY[l,j,i] += Float64(Y[i,l,k])
            end
          end
          ΔY[l,j,i] /= sl
        end
      end
    end

  end

  return nothing
end





"""
    sde!(LA   ::Array{Float64,2},
         ΔX   ::Array{Float64,3}, 
         ΔY   ::Array{Float64,3},
         m    ::Int64,
         ntip ::Int64)

Estimate the lineage averages for the SDE of trait evolution
"""
function sde!(LA   ::Array{Float64,2},
              ΔX   ::Array{Float64,3}, 
              ΔY   ::Array{Float64,3},
              m    ::Int64,
              ntip ::Int64)

  @inbounds begin

    for i = Base.OneTo(m), j = Base.OneTo(ntip)

      LA[i,j] = 0.0

      for l = Base.OneTo(ntip)
        l == j && continue
        x = ΔX[l,j,i]
        isnan(x) && continue
        y = ΔY[l,j,i]
        iszero(y) && continue
        @fastmath LA[i,j] += sign(x) * y * exp(-abs(x))
      end
    end
  end

  return nothing
end





"""
    lindiff!(LD   ::Array{Float64,3},
             ΔX   ::Array{Float64,3},
             Y    ::Array{Int64,3},
             m    ::Int64,
             ntip ::Int64,
             narea::Int64)

Estimate area-specific distance averages 
for colonization and extinction rates.
"""
function lindiff!(LD   ::Array{Float64,3},
                  ΔX   ::Array{Float64,3},
                  Y    ::Array{Int64,3},
                  m    ::Int64,
                  ntip ::Int64,
                  narea::Int64)
  @inbounds begin

    for i = Base.OneTo(m), k = Base.OneTo(narea), l = Base.OneTo(ntip)

      LD[i,l,k] = 0.0
      sk        = 0.0

      for j = Base.OneTo(ntip)
        j == l && continue
        x = ΔX[j,l,i]
        isnan(x) && continue
        y          = Float64(Y[i,j,k])
        LD[i,l,k] += abs(x)*y
        sk        += y
      end

      LD[i,l,k] /= (iszero(sk) ? 1.0 : sk)
    end
  end

  return nothing
end





"""
    Xupd_linavg(k::Int64, 
                wci::Array{Int64,1},
                X::Array{Float64,2},
                Y::Array{Int64,3},
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
                      Y    ::Array{Int64,3},
                      δyi  ::Array{Float64,2},
                      narea::Int64)

  @inbounds begin

    @simd for l = wci
      δxi[xj,l] = xpi[xj] - xpi[l]
      δxi[l,xj] = xpi[l]  - xpi[xj]
    end

    for l = wci
      lai[l] = 0.0
      for j = wci
        j == l && continue
        lai[l] += sign(δxi[j,l]) * δyi[j,l] * exp(-abs(δxi[j,l]))
      end
    end

    for k = Base.OneTo(narea), l = wci
      ldi[l,k] = 0.0
      sk       = 0.0
      for j = wci
        j == l && continue
        y         = Float64(Y[xi,j,k])
        ldi[l,k] += abs(δxi[j,l])*y
        sk       += y
      end
      ldi[l,k] /= (iszero(sk) ? 1.0 : sk)
    end
  end

  return nothing
end





#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# """
#     area_lineage_means!(AA::Array{Float64,2}, LA::Array{Float64,2}, X::Array{Float64,2}, Y::Array{Int64,3}, wcol::Array{Array{Int64,1},1}, m::Int64)

# Estimate area means according to presence
# absence of species and linage means according
# to area averages.
# """
# function area_lineage_means!(AA   ::Array{Float64,2}, 
#                              LA   ::Array{Float64,2},
#                              AO   ::Array{Int64,2},
#                              X    ::Array{Float64,2}, 
#                              Y    ::Array{Int64,3}, 
#                              wcol ::Array{Array{Int64,1},1},
#                              m    ::Int64,
#                              narea::Int64)

#   @inbounds begin

#     for i = Base.OneTo(m)

#       # area averages
#       for k = Base.OneTo(narea)
#         AA[i,k] = 0.0::Float64
#         sumY    = 0.0::Float64
#         AO[i,k] = 0::Int64
#         @simd for j = wcol[i]::Array{Int64,1}
#           if Y[i,j,k]::Int64 == 1
#             AA[i,k] += X[i,j]::Float64
#             sumY    += 1.0::Float64
#             AO[i,k]  = 1::Int64
#           end
#         end
#         if sumY != 0.0
#           AA[i,k] /= sumY::Float64
#         end
#       end

#       # lineage average
#       for j = wcol[i]::Array{Int64,1}
#         LA[i,j] = 0.0::Float64
#         sumY    = 0.0::Float64
#         @simd for k = Base.OneTo(narea) 
#           if Y[i,j,k]::Int64 == 1
#             LA[i,j] += AA[i,k]::Float64
#             sumY    += 1.0::Float64
#           end
#         end
#         LA[i,j] /= sumY::Float64
#       end
#     end
#   end

#   return nothing
# end




# """
#     linarea_diff!(LD::Array{Float64,3}, X::Array{Float64,2}, AA::Array{Float64,2}, narea::Int64, ntip::Int64, m::Int64)

# Create multi-dimensional array with 
# lineage and area averages differences in X.
# """
# function linarea_diff!(LD   ::Array{Float64,3},
#                        X    ::Array{Float64,2},
#                        AA   ::Array{Float64,2},
#                        AO   ::Array{Int64,2},
#                        narea::Int64,
#                        ntip ::Int64,
#                        m    ::Int64)
#   @inbounds begin

#     for k = Base.OneTo(narea), j = Base.OneTo(ntip), i = Base.OneTo(m)
#       if isnan(X[i,j])
#         continue
#       end
#       setindex!(LD, 
#                 (iszero(AO[i,k]) ? 0.0 : abs(X[i,j] - AA[i,k]))::Float64, 
#                 i, j, k)
#     end
#   end

#   return nothing
# end




# """
#     Xupd_linavg(k::Int64, wci::Array{Int64,1}, X::Array{Float64,2}, Y::Array{Int64,3}, narea::Int64)

# Re-estimate lineage specific means 
# for a branch update.
# """
# function Xupd_linavg!(aa   ::Array{Float64,1},
#                       la   ::Array{Float64,1},
#                       ld   ::Array{Float64,2},
#                       ao   ::Array{Int64,2},
#                       i    ::Int64, 
#                       wci  ::Array{Int64,1},
#                       xi   ::Array{Float64,1},
#                       Y    ::Array{Int64,3},
#                       narea::Int64)
#   @inbounds begin

#     for k = Base.OneTo(narea)
#       sumY  = 0.0::Float64
#       aa[k] = 0.0::Float64
#       @simd for j = wci
#         if Y[i,j,k] == 1
#           aa[k] += xi[j]::Float64
#           sumY  += 1.0
#         end
#       end
#       if sumY != 0.0
#         aa[k] /= sumY::Float64
#       end
#     end

#     fill!(la, NaN)
#     for j = wci
#       la[j] = 0.0
#       sumY  = 0.0::Float64
#       @simd for k = Base.OneTo(narea)
#         if Y[i,j,k]::Int64 == 1
#           la[j] += aa[k]::Float64
#           sumY  += 1.0
#         end
#       end
#       la[j] /= sumY::Float64
#     end

#     fill!(ld, NaN)
#     for k = Base.OneTo(narea), j = wci
#       setindex!(ld, (iszero(ao[i,k]) ? 0.0 : abs(xi[j] - aa[k])), j, k)
#     end

#   end

#   return nothing 
# end

