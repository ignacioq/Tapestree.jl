#=

Diffused Brownian motion simulation on fixed trees

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#




"""
    sim_dbm(tree::iTree, 
            x0  ::Float64,
            σ0  ::Float64,
            α   ::Float64,
            γ   ::Float64,
            δt  ::Float64)

Simulate a diffused Brownian motion given starting values.
"""
function sim_dbm(tree::iTree, 
                 x0  ::Float64,
                 σ20  ::Float64,
                 α   ::Float64,
                 γ   ::Float64,
                 δt  ::Float64)

  _sim_dbm(tree, x0, log(σ20), α, γ, δt, sqrt(δt))
end




"""
    _sim_dbm(tree::iTree, 
             x0  ::Float64,
             lσ0 ::Float64,
             α   ::Float64,
             γ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Simulate a diffused Brownian motion given starting values recursively.
"""
function _sim_dbm(tree::iTree, 
                  x0  ::Float64,
                  lσ0 ::Float64,
                  α   ::Float64,
                  γ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  et = e(tree)

  # simulate dbm
  if iszero(et)
    xv   = Float64[x0,  x0]
    lσ   = Float64[lσ0, lσ0]
    fdti = 0.0
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    xv, lσ  = dbm(x0, lσ0, α, γ, fdti, δt, srδt, nt)
  end

  if def1(tree)
    if def2(tree)
      x0  = xv[end]
      lσ0 = lσ[end]
      sTxs(_sim_dbm(tree.d1, x0, lσ0, α, γ, δt, srδt), 
           _sim_dbm(tree.d2, x0, lσ0, α, γ, δt, srδt), 
           et, δt, fdti, xv, lσ)
    else
      sTxs(_sim_dbm(tree.d1, xv[end], lσ[end], α, γ, δt, srδt), 
           et, δt, fdti, xv, lσ)
    end
  else
    sTxs(et, δt, fdti, xv, lσ)
  end
end




"""
    sim_dbm(tree::Tlabel, 
            x0  ::Float64,
            σ0  ::Float64,
            α   ::Float64,
            γ   ::Float64,
            δt  ::Float64)

Simulate a diffused Brownian motion given starting values.
"""
function sim_dbm(tree::Tlabel, 
                 x0  ::Float64,
                 σ0  ::Float64,
                 α   ::Float64,
                 γ   ::Float64,
                 δt  ::Float64)
  
  xs = Dict{String, Float64}()
  tr = _sim_dbm(tree, x0, log(σ0), α, γ, δt, sqrt(δt), xs)
  return tr, xs
end




"""
    _sim_dbm(tree::Tlabel, 
             x0  ::Float64,
             lσ0 ::Float64,
             α   ::Float64,
             γ   ::Float64,
             δt  ::Float64,
             srδt::Float64,
             xs  ::Dict{String, Float64})

Simulate a diffused Brownian motion given starting values recursively.
"""
function _sim_dbm(tree::Tlabel, 
                  x0  ::Float64,
                  lσ0 ::Float64,
                  α   ::Float64,
                  γ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  xs  ::Dict{String, Float64})

  et = e(tree)

  # simulate dbm
  if iszero(et)
    xv   = Float64[x0,  x0]
    lσ   = Float64[lσ0, lσ0]
    fdti = 0.0
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    xv, lσ = dbm(x0, lσ0, α, γ, δt, fdti, srδt, nt)
  end

  if def1(tree)
    if def2(tree)
      x0  = xv[end]
      lσ0 = lσ[end]
      sTxs(_sim_dbm(tree.d1, x0, lσ0, α, γ, δt, srδt, xs), 
           _sim_dbm(tree.d2, x0, lσ0, α, γ, δt, srδt, xs), 
           et, δt, fdti, xv, lσ)
    else
      xs[label(tree)] = xv[end]
      sTxs(_sim_dbm(tree.d1, xv[end], lσ[end], α, γ, δt, srδt, xs), 
           et, δt, fdti, xv, lσ)
    end
  else
    xs[label(tree)] = xv[end]
    sTxs(et, δt, fdti, xv, lσ)
  end
end



