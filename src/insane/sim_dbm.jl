#=

Diffused Brownian motion simulation on fixed trees

Ignacio Quintero Mächler

t(-_-t)

Created 25 01 2024
=#




"""
    sim_dbm(tree::iTree; 
            x0  ::Float64 = 0.0,
            αx  ::Float64 = 0.0,
            σ20 ::Float64 = 0.1,
            ασ  ::Float64 = 0.0,
            γ   ::Float64 = 0.1,
            δt  ::Float64 = 1e-3)

Simulate a diffused Brownian motion given starting values.
"""
function sim_dbm(tree::iTree; 
                 x0  ::Float64 = 0.0,
                 αx  ::Float64 = 0.0,
                 σ20 ::Float64 = 0.1,
                 ασ  ::Float64 = 0.0,
                 γ   ::Float64 = 0.1,
                 δt  ::Float64 = 1e-3)

  _sim_dbm(tree, x0, αx, log(σ20), ασ, γ, δt, sqrt(δt))
end




"""
    _sim_dbm(tree::iTree, 
             x0  ::Float64,
             αx  ::Float64,
             lσ20::Float64,
             ασ  ::Float64,
             γ   ::Float64,
             δt  ::Float64,
             srδt::Float64)

Simulate a diffused Brownian motion given starting values recursively.
"""
function _sim_dbm(tree::iTree, 
                  x0  ::Float64,
                  αx  ::Float64,
                  lσ20::Float64,
                  ασ  ::Float64,
                  γ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  et = e(tree)

  # simulate dbm
  if iszero(et)
    xv   = Float64[x0,  x0]
    lσ2   = Float64[lσ20, lσ20]
    fdti = 0.0
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    xv, lσ2  = dbm(x0, αx, lσ20, ασ, γ, fdti, δt, srδt, nt)
  end

  if def1(tree)
    if def2(tree)
      x0  = xv[end]
      lσ20 = lσ2[end]
      sTxs(_sim_dbm(tree.d1, x0, αx, lσ20, ασ, γ, δt, srδt), 
           _sim_dbm(tree.d2, x0, αx, lσ20, ασ, γ, δt, srδt), 
           et, δt, fdti, xv, lσ2)
    else
      sTxs(_sim_dbm(tree.d1, xv[end], αx, lσ2[end], ασ, γ, δt, srδt), 
           et, δt, fdti, xv, lσ2)
    end
  else
    sTxs(et, δt, fdti, xv, lσ2)
  end
end




"""
    sim_dbm(tree::Tlabel, 
            x0  ::Float64,
            αx  ::Float64,
            σ20 ::Float64,
            ασ  ::Float64,
            γ   ::Float64,
            δt  ::Float64)

Simulate a diffused Brownian motion given starting values.
"""
function sim_dbm(tree::Tlabel, 
                 x0  ::Float64,
                 αx  ::Float64,
                 σ20 ::Float64,
                 ασ  ::Float64,
                 γ   ::Float64,
                 δt  ::Float64)
  
  xs = Dict{String, Float64}()
  tr = _sim_dbm(tree, x0, αx, log(σ20), ασ, γ, δt, sqrt(δt), xs)
  return tr, xs
end




"""
    _sim_dbm(tree::Tlabel, 
             x0  ::Float64,
             αx  ::Float64,
             lσ20::Float64,
             ασ  ::Float64,
             γ   ::Float64,
             δt  ::Float64,
             srδt::Float64,
             xs  ::Dict{String, Float64})

Simulate a diffused Brownian motion given starting values recursively.
"""
function _sim_dbm(tree::Tlabel, 
                  x0  ::Float64,
                  αx  ::Float64,
                  lσ20::Float64,
                  ασ  ::Float64,
                  γ   ::Float64,
                  δt  ::Float64,
                  srδt::Float64,
                  xs  ::Dict{String, Float64})

  et = e(tree)

  # simulate dbm
  if iszero(et)
    xv   = Float64[x0,  x0]
    lσ2   = Float64[lσ20, lσ20]
    fdti = 0.0
  else
    nt, fdti = divrem(et, δt, RoundDown)
    nt = Int64(nt)

    if iszero(fdti)
      fdti = δt
      nt  -= 1
    end

    xv, lσ2 = dbm(x0, αx, lσ20, ασ, γ, δt, fdti, srδt, nt)
  end

  if def1(tree)
    if def2(tree)
      x0  = xv[end]
      lσ20 = lσ2[end]
      sTxs(_sim_dbm(tree.d1, x0, αx, lσ20, ασ, γ, δt, srδt, xs), 
           _sim_dbm(tree.d2, x0, αx, lσ20, ασ, γ, δt, srδt, xs), 
           et, δt, fdti, xv, lσ2)
    else
      xs[label(tree)] = xv[end]
      sTxs(_sim_dbm(tree.d1, xv[end], αx, lσ2[end], ασ, γ, δt, srδt, xs), 
           et, δt, fdti, xv, lσ2)
    end
  else
    xs[label(tree)] = xv[end]
    sTxs(et, δt, fdti, xv, lσ2)
  end
end



