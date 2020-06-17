#=

MuSSE equations

Ignacio Quintero Mächler

t(-_-t)

16 02 2020

=#





"""
    fbisse(du::Array{Float64,1}, 
           u::Array{Float64,1}, 
           p::Array{Float64,1}, 
           t::Float64)

BiSSE differential equations.
"""
function fbisse(du::Array{Float64,1}, 
                u::Array{Float64,1}, 
                p::Array{Float64,1}, 
                t::Float64)
  @inbounds begin
    du[1]= -1.0 * (p[1] + p[3] + p[5]) * u[1] + p[5] * u[2] + 2.0 * p[1] * u[3] * u[1]
    du[2]= -1.0 * (p[2] + p[4] + p[6]) * u[2] + p[6] * u[1] + 2.0 * p[2] * u[4] * u[2]
    du[3]= -1.0 * (p[1] + p[3] + p[5]) * u[3] + p[3] + p[5] * u[4] + p[1] * u[3]^2
    du[4]= -1.0 * (p[2] + p[4] + p[6]) * u[4] + p[4] + p[6] * u[3] + p[2] * u[4]^2
  end
  
  return nothing
end





"""
    fbisseE(du::Array{Float64,1}, 
            u::Array{Float64,1}, 
            p::Array{Float64,1}, 
            t::Float64)

BiSSE extinction differential equations.
"""
function fbisseE(du::Array{Float64,1}, 
                 u::Array{Float64,1}, 
                 p::Array{Float64,1}, 
                 t::Float64)
  @inbounds begin
    du[1] = -1.0 * (p[1] + p[3] + p[5]) * u[1] + p[3] + p[5] * u[2] + p[1] * u[1]^2
    du[2] = -1.0 * (p[2] + p[4] + p[6]) * u[2] + p[4] + p[6] * u[1] + p[2] * u[2]^2
  end
  
  return nothing
end





"""
    make_fbisseM(afE!, idxl::Int64)

Make function for BiSSE Matrix form of differential equations for flow 
algorithm and given an Extinction approximation function.
"""
function make_fbisseM(afE!, idxl::Int64)

  rE  = Array{Float64,1}(undef,2)

  function f(du::Array{Float64,2}, 
             u ::Array{Float64,2}, 
             p ::Array{Array{Float64,1},1}, 
             t ::Float64)
    @inbounds begin
      par = p[idxl]
      afE!(t, rE, p)
      du[1,1] = -1.0 * (par[1] + par[3] + par[5]) * u[1,1] + par[5] * u[2,1] + 2.0 * par[1] * rE[1] * u[1,1]
      du[2,1] = -1.0 * (par[2] + par[4] + par[6]) * u[2,1] + par[6] * u[1,1] + 2.0 * par[2] * rE[2] * u[2,1]
      du[1,2] = -1.0 * (par[1] + par[3] + par[5]) * u[1,2] + par[5] * u[2,2] + 2.0 * par[1] * rE[1] * u[1,2]
      du[2,2] = -1.0 * (par[2] + par[4] + par[6]) * u[2,2] + par[6] * u[1,2] + 2.0 * par[2] * rE[2] * u[2,2]
    end

    return nothing
  end
end

