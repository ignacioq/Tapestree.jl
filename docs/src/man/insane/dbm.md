# Diffused Brownian motion process (DBM)

This is a model of trait evolution where an univariate trait ``x(t)`` evolves under a diffused Brownian motion with an underlying evolutionary rate ``\sigma^2(t)`` that is also itself evolving separately according to Geometric Brownian motion:
```math
dx(t) = \alpha_x dt + \sigma(t) d W(t), \\
d(\text{ln}(\sigma^2(t)) = \alpha_{\sigma} dt + \gamma d W(t)
```
where ``\alpha_x`` is the trait drift (general trait tendency to increase or decrease), ``\alpha_{\sigma}`` is the drift in evolutionary rates, and ``\gamma`` represents the heterogeneity in evolutionary rates.


## Simulations

In the DBM, we are not modeling the realization of the tree, but rather a process that evolved along the tree. Thus, to simulate, we need a tree to be used as template so that we can simulate on top. So if we have a `tree` object of type `sT_label` or `sTf_label`, we can simulate a trait with a starting trait value of `x0` (``x(t = 0)``), with drift `αx` (``\alpha_x``), undergoing a starting rate of `σ20` (``\sigma^2(t = 0)``) with drift `ασ` (``\alpha_{\sigma}``), with a maximal discretization time step of `δt`:

```julia
sim_dbm(tree, x0 = 0.0, αx = 0.0, σ20 = 0.1, ασ = 0.0, γ = 0.1, δt = 1e-3)
```

which returns a tree of type `sTxs`, holding the simulation.

It might also be useful to return a `Dictionary` with the final trait values at each sampled species. For this, use the same function, in the same argument order but not using named (keyword) arguments. For instance, the same simulation above, but returning both a `sTxs` tree named `tr` and a dictionary of species values named `xd`:

```julia
tr, xd = sim_dbm(tree, 0.0, 0.0, 0.1, 0.0, 0.1, 1e-3)
```

We can plot the resulting tree using Tapestree's plot recipes ([Insane plots](@ref)). For example to plot the trait evolution colored by the logarithmic rates:
```julia
plot(xv, tr, zf = lσ2)
```

## Inference


For a given `sT_label` or `sTf_label` type `tree` object and a Dictionary `xavg` with a `String` key pointing to a `Float64` number (_i.e._, `Dict{String, Float64}`), where matching tip labels point to the trait value. If species labels are not included in the dictionary, the trait value is assumed missing. 

```julia
r, td = insane_dbm(tree, 
                   xavg,
                   γ_prior  = (0.05, 0.05),
                   αx_prior = (0.0, 10.0),
                   ασ_prior = (0.0, 10.0),
                   nburn    = 10_000,
                   niter    = 100_000,
                   nthin    = 1_000,
                   nflush   = 1_000,
                   ofile    = "<...directory...>",
                   δt       = 1e-3)
```

Finally, error or uncertainty around trait values can be included (assuming Normal variance) by setting another Dictionary, called say `xsv` (also a `Dict{String, Float64}`) where tip values point to the variance around trait values. Again, if tip labels are not included in this dictionary, it is assumed that there is no error around tip values.
Then you specify this dictionary on the argument `xs`:
```julia
r, td = insane_dbm(tree, 
                   xavg, 
                   xs       = xsv,
                   γ_prior  = (0.05, 0.05),
                   αx_prior = (0.0, 10.0),
                   ασ_prior = (0.0, 10.0),
                   nburn    = 10_000,
                   niter    = 100_000,
                   nthin    = 1_000,
                   nflush   = 1_000,
                   ofile    = "<...directory...>",
                   δt       = 1e-3)
```

Full documentation
```@docs
sim_dbm
insane_dbm
```
