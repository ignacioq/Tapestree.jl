# Constant occurrence birth-death process (COBD)

We now consider that species have an instantaneous rate of being sampled as fossil occurrences that are not placed in the phylogeny, ``\omega``, besides the constant rate of speciation ``\lambda``, extinction ``\mu`` and sampling of phylogenetic fossils ``\psi``.

## Simulations

To simulate under the COBD model, say, for a period of ``10``, with speciation rate of ``0.5``, extinction rate of ``0.3``, fossilization rate of ``0.4``, and occurrence sampling rate of ``2.0``, one can use
```julia
tr, ωtimes = sim_cobd(10.0, 0.5, 0.3, 0.4, 2.0)
```

where `ωtimes` represents the record of occurrences, i.e. a vector of sampling times.

## Inference

To perform inference on the COBD model, we can use
```julia
out, trees = insane_cobd(tree,
                         ωtimes,
                         nburn    = 1_000,
                         niter    = 50_000,
                         nthin    = 50, 
                         ofile    = "<directory>",
                         λ_prior  = (1.5, 1.0),
                         μ_prior  = (1.5, 1.0),
                         ψ_prior  = (1.0, 1.0),
                         ω_prior  = (1.0, 0.2),
                         survival = true,
                         tρ       = Dict("" => 1.0))
```
where we now have a Gamma prior for the speciation ``\lambda``, extinction ``\mu``, fossilization ``\psi`` and occurrence sampling rates ``\omega``, and we can also specify if we want to condition on survival of the process with `survival`.

!!! note
    the `tree` object must be of type `sTf_label` (see [Insane tree and model input/output](@ref)), which is the automatic type when reading a tree with `read_newick`.


# Episodic occurrence birth-death process (eOBD)

We now consider that the fossilization rate can vary through time in piece-wise constant fashion, for instance, for different stratigraphic epochs. The periods are not inferred but must instead be _pre_-defined.


## Simulations

To simulate under the eOBD model, for a period of, say, ``10``, with speciation rate of ``0.5``, extinction rate of ``0.3``, and fossilization rates of ``\{0.4, 0.1, 0.8\}`` and occurrence sampling rates ``\{2.4, 1.1, 1.8\}`` at periods ``(10, 7)``, ``(7, 3)`` and ``(3, 0)``:

```julia
tr, ωtimes = sim_cobd(10.0, 0.5, 0.3, [0.4, 0.1, 0.8], [2.4, 1.1, 1.8], [7.0, 3.0])
```

## Inference

To perform inference on the eOBD model, we input the times defining the periods for different fossilization and occurrence sampling rates `ψω_epoch`. If we define the periods used above for simulation, we can use

```julia
out, trees = insane_cobd(tree,
                         ωtimes,
                         nburn    = 1_000,
                         niter    = 50_000,
                         nthin    = 50, 
                         ofile    = "<directory>",
                         λ_prior  = (1.5, 1.0),
                         μ_prior  = (1.5, 1.0),
                         ψ_prior  = (1.0, 1.0),
                         ω_prior  = (1.0, 0.2),
                         ψω_epoch = [7.0, 3.0],
                         survival = true,
                         tρ       = Dict("" => 1.0))
```

# Occurrence birth-death diffusion process

This model follows from the eOBD model, but relaxing the assumption of homogeneous rates for instead allowing per-lineage instantaneous speciation rates ``\lambda(t)`` and extinction rates ``\mu(t)`` to follow separate Geometric Brownian motions (GBM), such that

```math
d(\text{ln}(\lambda(t)) = \alpha_{\lambda} dt + \sigma_{\lambda} d W(t), \\
d(\text{ln}(\mu(t)) = \alpha_{\mu} dt + \sigma_{\mu} d W(t),
```
where ``W(t)`` is the Wiener process (_i.e._, standard Brownian motion), ``\alpha_{\lambda}`` and ``\alpha_{\mu}`` are the drift in speciation and extinction, and ``\sigma_{\lambda}`` and ``\sigma_{\mu}`` are the diffusion for speciation and extinction, respectively. 

## Simulations

One can specify a time or a number of lineages to simulate under the OBDD. For example, for ``10`` species
```julia
sim_gbmobd(10, λ0 = 1.0, μ0 = 0.2, αλ = 0.0, αμ = 0.0, σλ = 0.1, σμ = 0.1, ψ = [0.1], ω = [1.0])
```

Or one can simulate for a given amount of time, say, ``5`` time units
```julia
sim_gbmobd(5.0, λ0 = 1.0, μ0 = 0.2)
```

## Inference

One can perform inference using:
```julia
out, trees = insane_gbmobd(tree,
                           ωtimes,
                           nburn    = 1_000,
                           niter    = 50_000,
                           nthin    = 50, 
                           ofile    = "<directory>",
                           λa_prior = (1.5, 1.0),
                           μa_prior = (1.5, 1.0),
                           αλ_prior = (0.0, 1.0),
                           αμ_prior = (0.0, 1.0),
                           σλ_prior = (3.0, 0.5),
                           σμ_prior = (3.0, 0.5),
                           ψ_prior  = (1.0, 1.0),
                           ω_prior  = (1.0, 1.0),
                           tρ       = Dict("" => 1.0))
```
where we have Gamma priors on the initial (root) speciation and extinction rates as well as on the sampling rates, Normal priors on the drift for speciation and extinction, and Inverse Gamma prior for the speciation and extinction diffusion rates. Usually, an informative prior such as `σμ_prior = (3.0, 0.1)` or more is often needed.

# Using fossil occurrences to inform ``\psi``

If the phylogenetic fossils included in `tree` have been selectively sampled to include only one occurrence per morphospecies, and if additional occurrences of these morphospecies exist, some of them can be added as specified in [Adding external fossil occurrences](@ref) to inform ``\psi``.
This option is complementary to the use of occurrences in the OBD process:
- Occurrences informing ``\psi``: assumes sampling from observed lineages in the reconstructed tree only, including more occurrences that necessary to correct the sampling bias can lead to underestimate past diversity and rates ;
- Occurrences informing ``\omega``: assumes sampling from all observed and unsampled lineages, no constraint to the number of occurrences that can be included.

Full documentation
```@docs
sim_gbmobd
insane_gbmobd
```
