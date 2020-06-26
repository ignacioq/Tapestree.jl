#=

insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 25 06 2020
=#


struct itree
    g::SimpleDiGraph{Int64}
    el::Dict{LightGraphs.SimpleGraphs.SimpleEdge{Int64}, Float64}
end