# Fossilized birth-death diffusion process

This model follows from the episodic fossilized birth-death model, where we relax the assumption of constant rates for instead allowing the per-lineage instantaneous speciation rates ``\lambda(t)`` and extinction rates ``\mu(t)`` to follow separate Geometric Brownian motions (GBM), such that

```math
d(\text{ln}(\lambda(t)) = \alpha_{\lambda} dt + \sigma_{\lambda} d W(t), \\
d(\text{ln}(\mu(t)) = \alpha_{\mu} dt + \sigma_{\mu} d W(t),
```
where ``W(t)`` is the Wiener process (_i.e._, standard Brownian motion), ``\alpha_{\lambda}`` and ``\alpha_{\mu}`` are the drift in speciation and extinction, and ``\sigma_{\lambda}`` and ``\sigma_{\mu}`` are the diffusion for speciation and extinction, respectively. 


## Simulations

One can specify a time or a number of lineages to simulate under the FBDD. For example, for ``50`` species
```julia
sim_gbmfbd(50, λ0 = 1.0, μ0 = 0.9, αλ = 0.0, αμ = 0.0, σλ = 0.1, σμ = 0.1),
```

Or one can simulate for a given amount of time, say, ``10`` time units
```julia
sim_gbmfbd(10.0, 1.0, 0.9, 0.0, 0.0, 0.1, 0.1),
```

## Inference

One can perform inference using:
```julia
r, tv = insane_gbmfbd(tree,
                      nburn    = 1_000,
                      niter    = 50_000,
                      nthin    = 50, 
                      ofile    = "<directory>",
                      λ0_prior = (0.05, 148.41),
                      μ0_prior = (0.05, 148.41),
                      αλ_prior = (0.0, 1.0),
                      αμ_prior = (0.0, 1.0),
                      σλ_prior = (0.05, 0.05),
                      σμ_prior = (3.0, 0.1),
                      tρ       = Dict("" => 1.0))
```
where we have log-normal priors on the initial (root) values for ``\lambda_0`` and ``\mu_0``, normal priors on the drift for speciation and extinction, and Inverse Gamma prior for the speciation and extinction diffusion rates. Usually, an informative prior such as `σμ_prior = (3.0, 0.1)` or more is often needed.

If _additional_ occurrences exist, they can be added as specified in [Adding external fossil occurrences](@ref)


Full documentation
```@docs
sim_gbmfbd
insane_gbmfbd
```
