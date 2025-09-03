# Constant pure-birth (Yule) process (CPB)

The simplest diversification model assumes no extinction (``\mu = 0``) and a constant speciation rate ``\lambda``, also known as, the pure-birth or Yule model.

## Simulations

To simulate a pure-birth tree one can use `sim_cb`. For instance, for a period of ``10`` time units and a speciation rate of ``\lambda = 0.5``:
```julia
tr = sim_cb(10.0, 0.5)
```

## Inference

To perform inference on a tree (of type `sT_label`), we can use the `insane_cb` function (cb = constant pure-birth).
```julia
r, tv = insane_cb(tree,
                   nburn  = 500,
                   niter  = 1_000,
                   nthin  = 2,
                   nflush = nthin,
                   ofile  = "<directory>")
```

!!! note 
    that For this specific CPB model, where the sampling fraction, ``\rho``, is ``1``, there are no unobserved components since we assume no extinction and that all species have been sampled. In this case, all the trees in the tree vector will be exactly the same.

In the following example, we now specify a global sampling fraction of `0.8`.
```julia
r, tv = insane_cb(tree,
                   nburn  = 500,
                   niter  = 1_000,
                   nthin  = 2,
                   nflush = nthin,
                   ofile  = "<directory>",
                   tρ     = Dict("" => 0.8))
```

This time, the DA trees are different from one another since they have data augmented lineages that represent that proportion of species not included in the tree. The position change from tree to tree because we are integrating over their unknown placement. You can check this by plotting the trees (check [Insane plots](@ref)). For instance, to plot the first tree in the vector:
```julia
plot(tv[1])
```

!!! note
    Insane uses Gibbs sampling across most parameter updates in INSANE. So, the prior for speciation in the CPB is a Gamma prior, and its parameters can be specified with argument `λ_prior`, The default is `λ_prior = (1.0, 1.0)`.


# Constant birth-death process (CBD)

We now consider that species have an instantaneous constant rate of extinction, ``\mu``, besides the constant rate of speciation ``\lambda``.

## Simulations

To simulate under the CBD model, say, for a period of ``10``, with speciation rate of ``0.5``, and extinction rate of ``0.3``, a one can use
```julia
tr = sim_cbd(10.0, 0.5, 0.3)
```

## Inference

To perform inference on the CBD model, we can use
```julia
r, tv = insane_cbd(tree,
                   nburn    = 1_000,
                   niter    = 100_000,
                   nthin    = 100, 
                   ofile    = "<out file directory>",
                   λ_prior  = (1.0, 1.0),
                   μ_prior  = (1.0, 1.0),
                   survival = true,
                   tρ       = Dict("" => 1.0))
```
where we now have a Gamma prior for the extinction ``\mu``, and we can also specify if we want to condition on survival of the process with `survival`.


Full documentation
```@docs
sim_cb
sim_cbd
insane_cb
insane_cbd
```

