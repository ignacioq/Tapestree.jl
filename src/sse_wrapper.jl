#=

Wrapper

Ignacio Quintero Mächler

t(-_-t)

September 26 2017

=#


"""
    runSSE(tip_val ::Dict{Int64,Array{Float64,1}},
           edges   ::Array{Float64,2},
           edlen   ::Array{Float64,1},
           out_file::String;
           niter   ::Int64  = 10_000,
           nthin   ::Int64  = 10,
           model   ::String = "musse")

Run a SSE model.
"""
function runSSE(tip_val ::Dict{Int64,Array{Float64,1}},
                edges   ::Array{Int64,2},
                edlen   ::Array{Float64,1},
                out_file::String;
                niter   ::Int64   = 10_000,
                nthin   ::Int64   = 10,
                model   ::String  = "musse",
                λpriors ::Float64 = .1,
                μpriors ::Float64 = .1,
                qpriors ::Float64 = .1)

end
