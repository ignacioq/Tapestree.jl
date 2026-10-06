#=

Survival conditioning using simulation proposals.

Ignacio Quintero Mächler

t(-_-t)

Created 16 11 2021
=#



"""
    m_survival(f::Function, ntry::Int64, surv::Int64, args...)

Sample the total number of `m` trials until both simulations survive
for birth-death model in `f`.
"""
function m_survival(f::Function, ntry::Int64, surv::Int64, args...)

  ntries = 1
  m      = 1.0

  # if survival of process with 1 lineage
  if isone(surv)

    while true
      s1, n1 = f(args..., false, 1)

      s1 && break
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end

  # if survival of process with 2 lineages
  elseif surv === 2

    while true
      s1, n1 = f(args..., false, 1)

      if s1
        s2, n2 = f(args..., false, 1)
        s2 && break
      end
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end
  end

  return m
end



