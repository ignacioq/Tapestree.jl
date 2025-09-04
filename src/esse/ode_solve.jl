#=

Create ODE numerical integration function

Ignacio Quintero Mächler

t(-_-t)

September 19 2017

=#







"""
    solvef(int::OrdinaryDiffEqCore.ODEIntegrator,
           u  ::Array{Float64,1},
           ti ::Float64,
           tf ::Float64)

Solve an IDE integrator for new `u` and and `ti` and `tf`.
"""
function solvef(int::OrdinaryDiffEqCore.ODEIntegrator,
                u  ::Array{Float64,1},
                ti ::Float64,
                tf ::Float64)
  @inbounds begin
    reinit!(int, u, t0 = ti, tf = tf)
    return solve!(int).u[1]
  end
end





# """
#     make_solver(odef, p0::Array{Float64,1})
# Make **ODE** solver for `odef`.
# """
# function make_solver(odef, p0::Array{Float64,1}, u0::Array{Float64,1})

#   prob = ODEProblem(odef, u0, (0.0,1.0), p0)

#   int = init(prob,
#              Tsit5(),
#              save_everystep  = false, 
#              calck           = false,
#              force_dtmin     = true,
#              save_start      = false,
#              initialize_save = false,
#              maxiters        = 100_000_000,
#              verbose         = false)

#   function f(u ::Array{Float64,1}, 
#              p ::Array{Float64,1},
#              ti::Float64,
#              tf::Float64)
#     int.p = p
#     reinit!(int, u, t0 = ti, tf = tf)
#     solve!(int).u[1]
#   end
#   return f
# end





