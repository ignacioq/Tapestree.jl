#=

Survival nodes for conditioning

Ignacio Quintero Mächler

t(-_-t)

Created 16 11 2021
=#




"""
    m_surv_cbd(t::Float64, λ::Float64, μ::Float64, ntry::Int64, surv::Int64)

Sample the total number of `m` trials until both simulations survive
for constant birth-death.
"""
function m_surv_cbd(t::Float64, λ::Float64, μ::Float64, ntry::Int64, surv::Int64)

  ntries = 1
  m      = 1.0

  # if survival of process with 1 lineage
  if isone(surv)

    while true
      s1, n1 = sim_cbd_surv(t, λ, μ, false, 1)

      s1 && break
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end

  # if survival of process with 2 lineages
  elseif surv === 2

    while true
      s1, n1 = sim_cbd_surv(t, λ, μ, false, 1)

      if s1
        s2, n2 = sim_cbd_surv(t, λ, μ, false, 1)
        s2 && break
      end
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end
  end

  return m
end




"""
    m_surv_gbmce(t   ::Float64,
                 lλ0 ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 μ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 ntry::Int64,
                 surv::Int64)

Sample the total number of `m` trials until both simulations survive
for `gbmce`.
"""
function m_surv_gbmce(t   ::Float64,
                      lλ0 ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      μ   ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ntry::Int64,
                      surv::Int64)
  ntries = 1
  m      = 1.0

  # if survival of process with 1 lineage
  if isone(surv)

    while true
      s1, n1 = _sim_gbmce_surv(t, lλ0, α, σλ, μ, δt, srδt, false, 1)

      s1 && break
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end

  # if survival of process with 2 lineages
  elseif surv === 2

    while true
      s1, n1 = _sim_gbmce_surv(t, lλ0, α, σλ, μ, δt, srδt, false, 1)
      if s1
        s2, n2 = _sim_gbmce_surv(t, lλ0, α, σλ, μ, δt, srδt, false, 1)
        s2 && break
      end
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end
  end
  return m
end




"""
    m_surv_gbmct(t   ::Float64,
                 lλ0 ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 ϵ   ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 ntry::Int64,
                 surv::Int64)

Sample the total number of `m` trials until both simulations survive
for `gbmct`.
"""
function m_surv_gbmct(t   ::Float64,
                      lλ0 ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      ϵ   ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ntry::Int64,
                      surv::Int64)
  ntries = 1
  m      = 1.0

  # if survival of process with 1 lineage
  if isone(surv)

    while true
      s1, n1 = _sim_gbmct_surv(t, lλ0, α, σλ, ϵ, δt, srδt, false, 1)

      s1 && break
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end

  # if survival of process with 2 lineages
  elseif surv === 2

    while true

      s1, n1 = _sim_gbmct_surv(t, lλ0, α, σλ, ϵ, δt, srδt, false, 1)

      if s1
        s2, n2 = _sim_gbmct_surv(t, lλ0, α, σλ, ϵ, δt, srδt, false, 1)
        s2 && break
      end
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end
  end

  return m
end




"""
    m_surv_gbmbd(t   ::Float64,
                 lλ0 ::Float64,
                 lμ0 ::Float64,
                 α   ::Float64,
                 σλ  ::Float64,
                 σμ  ::Float64,
                 δt  ::Float64,
                 srδt::Float64,
                 ntry::Int64,
                 c   ::Int64)

Sample the total number of `m` trials until both simulations survive
for `gbmbd`.
"""
function m_surv_gbmbd(t   ::Float64,
                      lλ0 ::Float64,
                      lμ0 ::Float64,
                      α   ::Float64,
                      σλ  ::Float64,
                      σμ  ::Float64,
                      δt  ::Float64,
                      srδt::Float64,
                      ntry::Int64,
                      surv::Int64)
  ntries = 1
  m      = 1.0

  # if survival of process with 1 lineage
  if isone(surv)

    while true
      s1, n1 = _sim_gbmbd_surv(t, lλ0, lμ0, α, σλ, σμ, δt, srδt, false, 1)

      s1 && break
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end

  # if survival of process with 2 lineages
  elseif surv === 2

    while true
      s1, n1 = _sim_gbmbd_surv(t, lλ0, lμ0, α, σλ, σμ, δt, srδt, false, 1)

      if s1
        s2, n2 = _sim_gbmbd_surv(t, lλ0, lμ0, α, σλ, σμ, δt, srδt, false, 1)
        s2 && break
      end
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end

  # no conditioning
  end

  return m
end




"""
    m_surv_gbmfbd(t   ::Float64,
                  λ0  ::Float64,
                  μ0  ::Float64,
                  αλ  ::Float64,
                  αμ  ::Float64,
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  ntry::Int64,
                  c   ::Int64)

Sample the total number of `m` trials until both simulations survive
for `fbdd`.
"""
function m_surv_gbmfbd(t   ::Float64,
                       λ0  ::Float64,
                       μ0  ::Float64,
                       αλ  ::Float64,
                       αμ  ::Float64,
                       σλ  ::Float64,
                       σμ  ::Float64,
                       δt  ::Float64,
                       srδt::Float64,
                       ntry::Int64,
                       c   ::Int64)

  ntries = 1
  m      = 1.0

  # if survival of process with 1 lineage
  if isone(c)

    while true
      s1, n1 = _sim_gbmfbd_surv(t, λ0, μ0, αλ, αμ, σλ, σμ, δt, srδt, false, 1)

      s1 && break
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end

  # if survival of process with 2 lineages
  elseif c === 2

    while true
      s1, n1 = _sim_gbmfbd_surv(t, λ0, μ0, αλ, αμ, σλ, σμ, δt, srδt, false, 1)

      if s1
        s2, n2 = _sim_gbmfbd_surv(t, λ0, μ0, αλ, αμ, σλ, σμ, δt, srδt, false, 1)
        s2 && break
      end
      ntries === ntry && break

      m      += 1.0
      ntries += 1
    end

  # no conditioning
  end

  return m
end

