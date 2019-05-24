# Tapestree

[![Travis](https://travis-ci.org/ignacioq/Tapestree.jl.svg?branch=master)](https://travis-ci.org/ignacioq/Tapestree.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/ks4wkrv6d4qw5wn1?svg=true)](https://ci.appveyor.com/project/ignacioq/tapestree-jl)
[![Coverage Status](https://coveralls.io/repos/ignacioq/Tapestree.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ignacioq/Tapestry.jl?branch=master)
[![codecov.io](http://codecov.io/github/ignacioq/Tapestree.jl/coverage.svg?branch=master)](http://codecov.io/github/ignacioq/Tapestry.jl?branch=master)

## Reference

Package that implements a joint model of trait evolution and biogeographic history described in: 
Quintero, I. and Landis, Michael J. Interdependent Phenotypic and Biogeographic Evolution Driven by Biotic Interactions. bioRxiv 560912 [link](https://doi.org/10.1101/560912).

## Usage

### Requirements:
  * Julia v1.1.x
  * `Tapestree` Package installed along with `RCall`, `Optim` and `ProgressMeter`. Install packages by typing `]` in the julia prompt and typing `add <package_name>`. For example, for Tapestree: `add Tapestree`.
  * R installed
  * R `ape` package installed.

### Inference


1. Open Julia v1.1.x

2. Load Tapestree package: 
```julia
using Tapestree
```

3. Specify the path to the phylogenetic tree (in a format that `ape` can read):
```julia
finches_tree_file = "/directory_where_Tapestree_was_cloned/Tapestree/data/finches_rescaled.tre"
```

4. Specify data. Data should be a `.txt` file where each row is a species, first the species name that matches the tree tip labels, second the phenotypic data and then the species presence in each area (`0` if absent and `1` if present) . Open `finches_pca1.txt` in the data folder to see an example.
```julia
finches_data_file = "/directory_where_Tapestree_was_cloned/Tapestree/data/finches_pca1.txt"
```

5. Specify output file (`homedir()` is an alias to your home folder)
```julia
out_file  = *(homedir(),"...")
```

6. Run the `tribe()` (TRIBE: Trait and Range Interspecific Biogeographic Evolution) model:
```julia
tribe(finches_tree_file, finches_data_file, out_file)
```

7. Further options for `tribe()` are
```julia
min_dt  = 0.01                       # a float describing the percentage of tree height allowed for discretization (lower values are more precise but take longer).
niter   = 10_000                     # an integer for the number of iterations.
nburn   = 5_000                      # an integer for the number of iterations in the adaptive burn-in phase.
nthin   = 100                        # an integer for the iteration sampling frequency.
saveXY  = (true, 1_000)              # a tuple of length 2: first is a boolean to save (or not) data augmented histories, second an integer for sampling frequency.
saveDM  = (true, 1_000)              # a tuple of length 2: first is a boolean to save (or not) data augmented deterministic effects, second an integer for sampling frequency.
ωxprior = (0.,10.)                   # a tuple of length 2 for the normal prior of ωx, first the mean, second the variance.
ω1prior = (0.,10.)                   # a tuple of length 2 for the normal prior of ω1, first the mean, second the variance.
ω0prior = (0.,10.)                   # a tuple of length 2 for the normal prior of ω0, first the mean, second the variance.
σ²prior = 1e-1                       # a float for the mean of the exponential prior for σ².
λprior  = 1e-1                       # a float for the mean of the exponential prior for both λs.
weight  = (0.15,0.05,0.02,0.02,5e-3) # a tuple of length 5 specifying the probabilities to update σ², ωx, ω1 & ω0, and λ1 & λ0 respectively.
λ1i     = 1.0                        # a float for the starting value for λ1.
λ0i     = 0.5                        # a float for the starting value for λ0.
ωxi     = 0.0                        # a float for the starting value for ωx.
ω1i     = 0.0                        # a float for the starting value for ω1.
ω0i     = 0.0                        # a float for the starting value for ω0.
fix_ωx  = false                      # a boolean to make inference without ωx.
fix_ω1  = false                      # a boolean to make inference without ω1.
fix_ω0  = false                      # a boolean to make inference without ω0.
```

8. The output is a `.log` file with the results of the MCMC chain, and optionally (if `saveXY = (true, k)`), an R data file (`.Rdata`) with the augmented data histories. R code to manipulate and visualize this output are provided upon request.


### Simulation

1. Specify the path to the phylogenetic tree (in a format that `ape` can read):
```julia
finches_tree_file = "/directory_where_Tapestree_was_cloned/Tapestree/data/finches_rescaled.tre"
```

2. Perform simulation (here with 0.0 as the inital trait value and 6 areas on the finches tree)
```julia
x_init  = 0.0
n_areas = 6
tip_values, tip_areas, tree, bts = simulate_tribe(x_init, n_areas, finches_tree_file)
```

3. Further options for `simulate_tribe()` are
```julia
ωx       = 0.0   # a float for simulated value of ωx.
σ²       = 0.5   # a float for simulated value of σ².
ω1       = 0.0   # a float for simulated value of ω1.
ω0       = 0.0   # a float for simulated value of ω0.
λ1       = 0.5   # a float for simulated value of λ1.
λ0       = 0.2   # a float for simulated value of λ0.
const_δt = 1e-4  # a float for the delta t used to approximate the simulation (lower values are more accurate but at a slight computation cost).
```

4. Specify output file (`homedir()` is an alias to your home folder)
```julia
out_file  = *(homedir(),"...")
```

5. Run the `tribe()` (optional parameters are the same as with inference):
```julia
tribe(tip_values, tip_areas, tree, bts, out_file)
```



