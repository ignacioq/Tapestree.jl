using Tapestree
using Test

# write your own tests here
@test 2 == 2

# read tree
tree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_5.tre"))

@test isa(tree, sT_label)