#=

Abstract insane tree structure

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#



"""
    iTree

An abstract type for all composite recursive types 
representing a binary phylogenetic tree for `insane` use
"""
abstract type iTree end




"""
    iTree

An abstract type for all composite recursive types 
representing a simple binary phylogenetic tree for `insane` use
"""
abstract type iT <: iTree end




"""
    iTgbm

An abstract type for all composite recursive types 
representing a binary phylogenetic tree with Geometric
Brownian motion rates for `insane` use
"""
abstract type iTgbm <: iTree end
