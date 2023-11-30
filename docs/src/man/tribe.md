# TRIBE

## Reference

Quintero, I., & Landis, M. J. (2020). Interdependent phenotypic and biogeographic evolution driven by biotic interactions. Systematic biology, 69(4), 739-755. [https://doi.org/10.1093/sysbio/syz082](https://doi.org/10.1093/sysbio/syz082)

## Example

Open Julia and load the Tapestree package: 
```julia
using Tapestree
```

Specify the path to the phylogenetic tree (in a format that `ape::read.tree()` can read):
```julia
finches_tree_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "finches_rescaled.tre")
```

Specify data. Data should be a `.txt` file where each row is a species, first 
the species name that matches the tree tip labels, second the phenotypic data 
and then the species presence in each area (`0` if absent and `1` if present). 
Open `finches_pca1.txt` in the data folder to see an example.
```julia
finches_data_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "finches_pca1.txt")
```

Specify output file (`homedir()` is an alias to your home folder)
```julia
out_file = *(homedir(),"...")
```

Run the `tribe()` (TRIBE: Trait and Range Interspecific Biogeographic Evolution) model:
```julia
tribe(finches_tree_file, finches_data_file, out_file)
```

For optional (keyword) arguments, see below

## Function Documentation
```@docs
tribe
simulate_tribe
```
