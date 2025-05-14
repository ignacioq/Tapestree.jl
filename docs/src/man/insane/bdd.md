# Birth-death diffusion (BDD) processes

These models all assume that the per-lineage instantaneous speciation rates ``\lambda(t)`` follow a Geometric Brownian motion (GBM), such that

```math
d(\text{ln}(\lambda(t)) = \alpha dt + \sigma_{\lambda} d W(t),
```
where ``\alpha`` is the drift, ``\sigma_{\lambda}`` is the diffusion for speciation, and ``W(t)`` is the Wiener process (_i.e._, standard Brownian motion), and holding different assumptions on extinction.

The likelihood is approximated using the Euler method to solve Stochastic Differential Equations, and require to set a minimum time step for the discretization of the augmented diffusion paths. In INSANE `δt` (by default ``10^{-3}``) sets the time step to perform the Euler approximation, and denotes the proportion with respect to the tree height (the duration of the tree), so the time step will be ``δt = th \times 10^-3``, where ``th`` is the tree height.

## Pure birth diffusion (``\mu(t) = 0``)

In the simplest BDD model, the pure-birth diffusion, we assume there is no extinction, that is ``\mu(t) = 0``.

### Simulations

For all the BDD models, we have the possibility to simulate conditioned on total simulation time or the number of species.

To simulate conditioned on some number of species, say, ``20``, with starting speciation rate of ``\lambda_0 = 1.0``, drift of ``\alpha = 0`` and diffusion of ``\sigma_{\lambda} = 0.1``, we can use:

```julia
 sim_gbmpb(20, λ0 = 1.0, α = 0.0, σλ = 0.1)
```

Similarly, to simulate conditioned on time, say, ``10`` time units, with the same parameters, we can use:

```julia
 sim_gbmpb(10.0, λ0 = 1.0, α = 0.0, σλ = 0.1)
```

Other options are available, such as `δt` which controls the time step of the simulation, which uses the Euler approximation. Similarly, `init` can be `:crown` or `:stem` to simulate starting with ``1`` or ``2`` lineages. 

!!! note
    The simulations conditioned on time must input a `Float64` as first (and mandatory) argument while those conditioned on number of species must rater input a `Int64`.


### Inference

To perform inference under this model we can use:

```julia
r, tv = insane_gbmpb(tree,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     α_prior  = (0.0, 10.0),
                     σλ_prior = (0.05, 0.05),
                     tρ    = Dict("" => 1.0))
```

Here, the prior for the drift `α_prior` is a Normal distribution, with the first element representing the mean and the second the standard deviation, and the prior for the diffusion in speciation rates. `σλ_prior`, is an Inverse Gamma. 


## Birth-death diffusion process with constant extinction (``\mu(t) = \mu``)

Here we assume that there is extinction, but that is constant across lineages and time.

### Simulations

As with the pure-birth diffusion, we can simulate conditioned on reaching some number of lineages or some time. Now, however, we need to specify the extinction rate as well, for instance to generate a tree with ``20`` species and some parameter values:

```julia
 sim_gbmce(20, λ0 = 1.0, α = 0.0, σλ = 0.1, μ = 0.7)
```


### Inference

To perform inference we can use
```julia
r, tv = insane_gbmce(tree,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     α_prior  = (0.0, 10.0),
                     σλ_prior = (0.05, 0.05),
                     μ_prior  = (1.0, 1.0),
                     tρ       = Dict("" => 1.0))
```
where we now specify the Gamma prior `μ_prior` for the extinction.


## Birth-death diffusion process with constant turnover (``\mu(t) = \epsilon \lambda(t)``)

One can also assume that turnover, ``\epsilon``, that is, the ration of extinction over speciation (_i.e._, ``\frac{\mu}{\lambda}``), is constant.

### Simulations

As with the other BDD models, we can simulate both conditioned on the number of species or in time. For instance, here we condition on time with a turnover rate of ``0.4``.

```julia
 sim_gbmct(10.0, λ0 = 1.0, α = 0.0, σλ = 0.1, ϵ = 0.4)
```



### Inference

Performing inference can then be done using:
```julia
r, tv = insane_gbmct(tree,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     α_prior  = (0.0, 10.0),
                     σλ_prior = (0.05, 0.05),
                     ϵ_prior  = (0.0, 100.0),
                     tρ       = Dict("" => 1.0))
```
where we now specify a Uniform prior `ϵ_prior` for the turnover. Also, Tapestree resorts to simple MH updates for ``\epsilon``, so an increased burnin is recommended. 


## Birth-death diffusion process

Finally, we arrive at the most complex model, where extinction also follows a GBM, that is,
```math
d(\text{ln}(\mu(t)) = \sigma_{\mu} d W(t),
```
where ``\sigma_{\mu}`` is the diffusion for extinction.

!!! note 
    This is unidentifiable, unless we specify strong priors on the extinction diffusion coefficient ``\sigma_{\mu}``, and still do not do well in recovering simulated values, at least for small trees.


### Simulations

Again, one can specify a time or a number of lineages to simulate under the BDD. For example, for ``50`` species

```julia
 sim_gbmbd(50, λ0 = 1.0, μ0 = 1.0, α = 0.0, σλ = 0.1, σμ = 0.1)
```

### Inference

One can perform inference using:
```julia
r, tv = insane_gbmbd(tree,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     λ0_prior = (0.05, 148.41),
                     μ0_prior = (0.05, 148.41),
                     α_prior  = (0.0, 10.0),
                     σλ_prior = (0.05, 0.05),
                     σμ_prior = (3.0, 0.1),
                     tρ       = Dict("" => 1.0))
```
where we have log-normal priors on the initial (root) values for ``\lambda_0`` and ``\mu_0``, and Inverse Gamma prior for sigma. Usually, an informative prior such as `σμ_prior = (3.0, 0.1)` or more is needed.



## Birth-death diffusion process with informed extinction

One can also define a branch-specific fixed extinction function and make inference on the GBM speciation. For example, in Quintero _et al._ (2024) _Science_, we estimated extinction rates based on fossil occurrences and then used this as input.

### Inference

To make inference under this model, we need two vectors where each element links to a given branch in the reconstructed tree, in the same order as the phylogenetic tree. Each element of these vectors are themselves a vector of type `Float64`, and correspond to the time at which extinction is recovered in time backward order (_e.g._, [1.1, 1.0, ..., 0.1, 0.08, 0.0]) and to the extinction at those times. INSANE uses a linear approximation function between the input sampled points.

To make inference under this model, we now input these vectors. Lets call them `time_vector` (of type `Vector{Vector{Float64}}`) and `extinction_vector` (also of type `Vector{Vector{Float64}}`), then:
```julia
r, tv = insane_gbmbd(tree,
                     time_vector,
                     extinction_vector,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     tρ       = Dict("" => 1.0))
```

The only thing one has to make sure is that the order of `time_vector` and `extinction_vector` correspond to the order of the branches in `tree`, which can be a bit involved. However, this can be achieved by understanding that the trees are recursive and are ordered either by left or right branch, and using some of the tree utility functions. 

For example, one can use the `subclade` function, to extract a subclade given a set of tip labels. One can also use the `make_idf` function, which splits the tree into a vector of branches in the same order that INSANE will process the tree. This way one can see the ordering of branches, and then order the `time_vector` and the `extinction_vector` in the same order.

If only one extinction function is used for all the tree, then one has to only create a vector of the time and extinction vector of equal size as the number of branches.

To estimate the number of branches in your empirical tree (saved in `<...tree directory...>`), we can do:

```julia
# read tree
tre = read_newick("<...tree directory...>")

## make extinction function vectors for each branch
# set a dummy smapling fraction (you can set the real one afterwards! - this is just to estimate the number of branches)
tρ  = Dict(li => 1.0 for li in tiplabels(tre))
idf = make_idf(tre, tρ, Inf)

# number of branches
nb = length(idf)
```

If we have a global curve with extinction `[0.06, 0.02, 0.05]` at times `[0.5, 0.3, 0.1]`, then we would have to create the following input objects for `time_vector` and `extinction_vector`:

```julia
extinction_vector = fill([0.06, 0.02, 0.05], nb)
time_vector = fill([0.5, 0.3, 0.1], nb)
```


Full documentation
```@docs
sim_gbmpb
sim_gbmce
sim_gbmct
sim_gbmbd
insane_gbmpb
insane_gbmce
insane_gbmct
insane_gbmbd
```