# Constant fossilized birth-death process (CFBD)

We now consider that species have an instantaneous rate of being sampled as fossils, ``\psi``, besides the constant rate of speciation ``\lambda`` and extinction ``\mu``.

## Simulations

To simulate under the CFBD model, say, for a period of ``10``, with speciation rate of ``0.5``, extinction rate of ``0.3``, and fossilization rate of ``0.4``, one can use
```julia
tr = sim_cfbd(10.0, 0.5, 0.3, 0.4)
```

## Inference

To perform inference on the CFBD model, we can use
```julia
r, tv = insane_cfbd(tree,
                    nburn    = 1_000,
                    niter    = 50_000,
                    nthin    = 50, 
                    ofile    = "<directory>",
                    λ_prior  = (1.0, 1.0),
                    μ_prior  = (1.0, 1.0),
                    ψ_prior  = (1.0, 1.0),
                    survival = true,
                    tρ       = Dict("" => 1.0))
```
where we now have a Gamma prior for the speciation ``\lambda``, extinction ``\mu`` and fossilization rates ``\psi``, and we can also specify if we want to condition on survival of the process with `survival`.

!!! note
    the `tree` object must be of type `sTf_label` [Insane tree and model input/output](@ref), which is the automatic type when reading a tree with `read_newick`.


# Episodic fossilized birth-death process (eFBD)

We now consider that the fossilization rate can vary through time in piece-wise constant fashion, for instance, for different stratigraphic epochs. The periods are not inferred but must instead be _pre_-defined.


## Simulations

To simulate under the eFBD model, for a period of, say, ``10``, with speciation rate of ``0.5``, extinction rate of ``0.3``, and fossilization rates of ``\{0.4, 0.1, 0.8\}`` at periods ``(10, 7)``, ``(7, 3)`` and ``(3, 0)``:

```julia
tr = sim_cfbd(10.0, 0.5, 0.3, [0.4, 0.1, 0.8], [7.0, 3.0])
```

## Inference

To perform inference on the eFBD model, we input the times defining the periods for different fossilization rates `ψ_epoch`. If we define the periods used above for simulation, we can use

```julia
r, tv = insane_cfbd(tree,
                    nburn    = 1_000,
                    niter    = 50_000,
                    nthin    = 50, 
                    ofile    = "<directory>",
                    λ_prior  = (1.0, 1.0),
                    μ_prior  = (1.0, 1.0),
                    ψ_prior  = (1.0, 1.0),
                    ψ_epoch  = [7.0, 3.0],
                    survival = true,
                    tρ       = Dict("" => 1.0))
```

# Adding external fossil occurrences

Finally, one can incorporate known fossil occurrences that are not included in the empirical tree. For instance, we might know of ``4`` fossil occurrences for a given species, but only one of them is represented in the empirical tree. Adding the other ``3`` occurrences might be desirable to better inform the fossilization rates, which has cascading effects on inferred speciation and extinction.

Given the eFBD model assumptions, this _additional_ fossil occurrences act as sampled ancestors, and their position on the tree does not matter (since we assume that all lineages at a given period share the same fossilization rate). Thus, we only need to input an Integer vector specifying the number of _additional_ fossil occurrences per period. Following the example above, suppose we have ``3``, ``2`` and ``5`` _additional_ occurrences per period, which we would specify in the `f_epoch` argument, as follows

```julia
r, tv = insane_cfbd(tree,
                    nburn    = 1_000,
                    niter    = 50_000,
                    nthin    = 50, 
                    ofile    = "<directory>",
                    λ_prior  = (1.0, 1.0),
                    μ_prior  = (1.0, 1.0),
                    ψ_prior  = (1.0, 1.0),
                    ψ_epoch  = [7.0, 3.0],
                    f_epoch  = [3, 2, 5],
                    survival = true,
                    tρ       = Dict("" => 1.0))
```

Full documentation
```@docs
sim_cfbd
insane_cfbd
```





