# ESSE

## Reference

Quintero, I., Landis, M. J., Jetz, W., & Morlon, H. (2022). The build-up of the present-day tropical diversity of tetrapods. Proceedings of the National Academy of Sciences, 2023. 120 (20) e2220672120. [https://doi.org/10.1073/pnas.2220672120](https://doi.org/10.1073/pnas.2220672120)


## Example

In this example we run some random data for 50 species across 2 areas, where
each area has a specific covariate that affects speciation and we also allow
for 2 hidden states.

Open Julia and load the Tapestree package: 
```julia
using Tapestree
```

Specify the path to the phylogenetic tree (in a format that `ape::read.tree()` can read):
```julia
tree_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_50.tre")
```

Specify state data. Data should be a `.txt` file where each row is a species, 
first the species name that matches the tree tip labels and each subsequent
column specifying presence or absence in a given area with `0` or `1`, 
respecitvely.
Open `st2_data.txt` in the data folder to see an example for 2 areas.
```julia
states_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "st2_data.txt")
```

Specify covariate data ``y = f(x)``. Data should be a `.txt` file where the 
first column is time ``x`` in backward fashion (the present is ``0`` and the 
past is ``> 0``), and the subsequent columns are the respective time covariates
for each area ``f(x)``. If there is only one covariate, the same is used across
all areas, if not, the number of covariates should match the specific model. 
More than one covariate per area is allowed, and in the case of covariates 
affecting colonization rates, they should match the number of possible 
colonization parameters between all areas.
Open `env_data_2.txt` in the data folder to see an example for 2 covariates for
2 areas.
```julia
states_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "st2_data.txt")
```

Specify output MCMC file (`homedir()` is an alias to your home folder)
```julia
out_file = *(homedir(),"...")
```

Specify the (optional) output file (`homedir()` is an alias to your home 
folder) for the node marginal probabilities. 
```julia
out_states = *(homedir(),"...")
```


Run the `esse()` (ESSE: Environmental and State dependent Speciation and 
Extinction) model, with covariates affecting speciation rates:
```julia
esse(tree_file, out_file, 2, envdata_file = envdata_file, 
  states_file = states_file, out_states = out_states, cov_mod = ("s",))
```

If one would like to make the covariates affect other rates, such as dispersal,
in addition to speciation rates, one would specify the following covariate
model `cov_mod = ("s","g")`. Note however that this has not been validated. 
Moreover, covariate effect on extinction is non retrievable from extant-only
phylogenetic trees.

### Parallel MC3 

It is encouraged to use Metropolis coupled MCMC (MC3) for more robust 
convergence (the posterior surfaced is highly peaked).

Load the Distributed package, set the number of processors for Julia, and
make Tapestree available to all (see the [Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/#man-distributed) 
package for more information). Below we add 3 processors.
```julia
using Distributed
addprocs(3)
@everywhere using Tapestree
```

### Set parameter constraints

To constrain parameters to be equal to one another or to fix them to be `0`, 
it is necessary to create equalities between parameters and pass them as an 
argument. In general parameters are specified as follows: 
  "<parameter name>\_<area>\_<hidden state>"
- Speciation is "lambda" (e.g. speciation rate for area A and hidden state 0: `lambda_A_0`)
- Local extinction is "loss" (e.g. local extinction rate for area B and hidden state 1: `loss_B_1`)
- colonization is "gain" (e.g. gain rate from A -> B and hidden state 1: `gain_AB_1`)
- the effect of the covariate is "beta\_<effect parameter name>" (e.g. effect of first covariate on speciation in area A and hidden state 0: `beta_lambda_1_A_0` - here the 1 after the lambda is because is the first covariate)
- transition between hidden states is "q" (e.g. transition rate from hidden state 0 -> 1 and hidden state 1: `q_01`). 
You have to make sure that the given constraints apply to the specification of 
model. In the following example, we constraint local and global extinction rates
to be the same across 2 areas and 2 hidden states , and constrain the hidden
states transition rates.
```julia
cpar = ("q_01 = q_10",
        "loss_A_0 = mu_A_0",
        "loss_B_0 = mu_B_0",
        "loss_A_1 = mu_A_1",
        "loss_B_1 = mu_B_1")
```

We now run esse, using MC3 with Metropolis-Hastings MCMC (here using 3 parallel chains).
```julia
esse(tree_file, out_file, 2,
     envdata_file = envdata_file,
     states_file  = states_file, 
     out_states   = out_states,
     constraints  = cpar,
     cov_mod      = ("s",),
     ncch         = 3,
     parallel     = true,
     niter        = 5_000,
     nthin        = 100,
     dt           = 0.8,
     nburn        = 1_000, 
     mc           = "mh",
     node_ps      = (true, 100))
```

For optional (keyword) arguments, see below

## Function Documentation
```@docs
esse
simulate_sse
```
