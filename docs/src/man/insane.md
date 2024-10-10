# INSANE

## Reference

Quintero, I., Lartillot, N., Morlon, H. (in prep). Imbalanced speciation pulses sustain the radiation of mammals. 


## Insane Bayesian data augmentation

INSANE uses Bayesian data augmentation (DA) to perform inference on a number of evolutionary models on phylogenetic trees. As such, performing inference will output posterior samples for the governing parameters as well as _complete_ or _data augmented_ trees, that is, trees that include probable configurations of unobserved variables such as the lineages that went extinct in the past or the underlying (latent) speciation rate. 

## Insane tree and model input/output

### Insane trees

INSANE uses different types of trees, but all are recursive structures. All trees have at least three fields: 
  * `d1` which corresponds to the left daughter (either another tree structure or `nothing` if it is a tip).
  * `d2` which corresponds to the right daughter (also either another tree structure or `nothing` if it is a tip).
  * `e` a decimal number specifying the edge length of the tree.

For all tree output from simulations and inference (all expect labelled trees):
  * `fx` a Boolean stating if tree is fixed or not. That is, if it corresponds to a sampled (`fx = true`) or unsampled (`fx = false`) lineage.

For all trees that allow for extinction:
  * `iμ` a Boolean stating if tree is extinct (`iμ = true`) or not (`iμ = false`).


### Insane simple trees (`sT`)

1. `sT_label`: is a simple labelled tree, and is the one required as input to perform inference. Note that all other model specific tree types can be converted to `sT_label` by using the latter as a function: `tree = sT_label(tree_of_other_type)`.
  * It has the additional field `l` for label.
2. `sTb`: is a simple pure-birth tree for the constant pure-birth (Yule) model.
3. `sTbd`: is a simple birth-death tree for the constant birth_death model.

### Insane BDD trees (`iT`)

All `iT` trees specify to birth-death diffusion trees with different extinction assumptions and all hold the following extra fields:

  * `dt` a decimal specifying the time step of the GBM discretization.
  * `fdt` a decimal specifying the final time step of the GBM discretization.
  * `lλ` an array specifying the Brownian motion evolution of log-speciation rates.

And have the following concrete types:

1. `iTb`, `iTce` and `iTct`: are BDD trees with no-extinction (pure-birth), constant extinction and constant turnover, respectively.
2. `iTbd`: is the full BDD tree.
  * It has the additional field `lμ`:  an array specifying the Brownian motion evolution of log-extinction rates.


### Reading and saving newick trees

Tapestree can read files in the simple newick format using the `read_newick` function:
```julia
tree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_50.tre"))
```

Note that the tree has type `sT_label`, which stands for simple labelled tree. You can check this using
```julia
typeof(tree)
```

Similarly, it can also write trees using `write_newick`
```julia
write_newick(tree, "<directory>")
```

### Reading and saving model output

All models can return and/or save the output by writing directly to a file on the fly and within the Julia session when the model finishes. There are two outputs:

1. The governing parameters trace, which is saved as a `.log` file and returned as the first object once the algorithm finished. These can be conveniently read in the Tracer software (https://github.com/beast-dev/tracer/releases/tag/v1.7.2)[https://github.com/beast-dev/tracer/releases/tag/v1.7.2].
2. The DA trees as a tree vector which are saved as a `.txt` file (of the same name as the `.log`) and returned as the second object once the algorithm finished. 

The DA trees written in the insane-specific `.txt` file can be read using the `iread()` function, which only needs the specific file directory as input, but also accepts an optional (keyword) argument `ix` that indicates the specific tree iterations to read, as an `OrdinalRange` object. For instance, to read only the first ``50`` trees, one can use `iread("<directory to txt>", ix = 1:50`. To read only the trees every ``10`` iterations from the ``100`` to ``400`` sampled, one can use `ix = 100:10:400`, and so on. This can be helpful to avoid high computation costs of reading very large files.

Similarly, you can save any individual or vector of insane trees using the `iwrite()` function: `iwrite(trees, "<directory to txt>")`.


## Insane models

### Tree input

All inference functions require a phylogenetic tree of type `sT_label`, that is, a simple labelled tree. This is the default object type when using the `read_newick` function. However, when using simulations from models or `iread`, the resulting type is specific to the model (for computational efficiency). One can easily create a tree of type `sT_label` from any other tree by using `sT_label(tree)`. The output can then be used to perform inference.


### Common (keyword) arguments across all insane inference models

* `nburn`: specifies the number of iterations to discard as burn-in.
* `niter`: specifies the number of MCMC iterations. 
* `nthin`: specifies the iteration frequency at which to save the parameters **in the julia session** (_i.e._,`nthin = 2` specifies saving every 2 iterations), 
* `nflush`: specifies the frequency at which to save **to file**. 
* `ofile`: specifies the directory where the results will be written. 
* `tρ`: controls the sampling fraction and receives a `Dictionary` as input, with a `String` key pointing to a `Float64` number (_i.e._, `Dict{String, Float64}`). If the dictionary is of length 1 with an empty string, then the insane sets this as a the global sampling fraction. For example, to set a sampling fraction of `0.6`, one show input `tρ = Dict("" => 0.6)`. Most times, however, sampling fraction is not uniform across the tree, but rather some part so the tree is more heavily sampled than others, to accomodate these variability, you can input a dictionary of the same length as the number of tips in the tree, where the dictionary key string is the tip label pointing to the specific sampling fraction value. For example, for two tips, named `tip_1` and `tip_2`, one could input `tρ = Dict("tip_1" => 0.5, "tip_2" => 0.3)`.
* `prints`: specifies the number of seconds to refresh the progress meter.
* `survival`: For those modesl with extinction, `true` or `false` if to condition the likelihood on survival.


### Constant rate models

#### Constant pure-birth (Yule) process (CB)

The simplest diversification model assumes no extinction and a constant speciation rate ``\lambda``, also known as, the pure-birth or Yule model.

##### Simulations

To simulate a pure-birth tree one can use `sim_cb`. For instance, for a period of ``10`` time units and a speciation rate of ``\lambda = 0.5``:
```julia
tr = sim_cb(10.0, 0.5)
```

Full documentation
```@docs
sim_cb
```

##### Inference

To perform inference on a tree (of type `sT_label`), we can use the `insane_cb` function (cb = constant birth).
```julia
r, tv = insane_cb(tree,
                  nburn  = 500,
                  niter  = 1_000,
                  nthin  = 2,
                  nflush = nthin,
                  ofile  = "<directory>")
```

Note that for this specific CB model, where the sampling fraction, ``\rho``, is ``1``, there are no unobserved components since we assume no extinction and that all species have been sampled. In this case, all the trees in the tree vector will be exactly the same.

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

Finally, it is important to note that we use Gibbs sampling across most parameter updates in INSANE. So, the prior for speciation in the CB is a Gamma prior, and its parameters can be specified with argument `λ_prior`, The default is `λ_prior = (1.0, 1.0)`.

Full documentation
```@docs
insane_cb
```

#### Constant birth-death process (CBD)

We now consider that species have an instantaneous constant rate of extinction, ``\mu``, besides the constant rate of speciation ``\lambda``.

##### Simulations

To simulate under the CBD model, say, for a period of ``10``, with speciation rate of ``0.5``, and extinction rate of ``0.3``, a one can use
```julia
tr = sim_cbd(10.0, 0.5, 0.3)
```

Full documentation
```@docs
sim_cbd
```

##### Inference

To perform  inference on the CBD model, we can use
```julia
r, tv = insane_cbd(tree,
                   nburn    = 1_000,
                   niter    = 50_000,
                   nthin    = 50, 
                   ofile    = "<directory>",
                   λ_prior  = (1.0, 1.0),
                   μ_prior  = (1.0, 1.0),
                   survival = true,
                   tρ       = Dict("" => 1.0))
```
where we now have a Gamma prior for the extinction ``\mu``, and we can also specify if we want to condition on survival of the process with `survival`.

Full documentation
```@docs
insane_cbd
```


### Birth-death diffusion (BDD) processes

These models all assume that the per-lineage instantaneous speciation rates ``\lambda(t)`` follow a Geometric Brownian motion (GBM), such that

```math
d(\text{ln}(\lambda(t)) = \alpha dt + \sigma_{\lambda} d W(t),
```
where ``\alpha`` is the drift, ``\sigma_{\lambda}`` is the diffusion for speciation, and ``W(t)`` is the Wiener process, and holding different assumptions on extinction.

The likelihood is approximated using the Euler method to solve Stochastic Differential Equations, and require to set a minimum time step for the discretization of the augmented diffusion paths. In INSANE `δt` (by default ``10^{-3}``) sets the time step to perform the Euler approximation, and denotes the proportion with respect to the tree height (the duration of the tree), so the time step will be ``δt = th \times 10^-3``, where ``th`` is the tree height.

#### Pure birth diffusion (``\mu(t) = 0``)

In the simplest BDD model, the pure-birth diffusion, we assume there is no extinction, that is ``\mu(t) = 0``.

##### Simulations

For all the BDD models, we have the possibility to simulate conditioned on total simulation time or the number of species.

To simulate conditioned on some number of species, say, ``20``, with starting speciation rate of ``\lambda_0 = 1.0``, drift of ``\alpha = 0`` and diffusion of ``\sigma_{\lambda} = 0.1``, we can use:

```julia
 sim_gbmb(20, λ0 = 1.0, α = 0.0, σλ = 0.1)
```

Similarly, to simulate conditioned on time, say, ``10`` time units, with the same parameters, we can use:

```julia
 sim_gbmb(10.0, λ0 = 1.0, α = 0.0, σλ = 0.1)
```

Other options are available, such as `δt` which controls the time step of the simulation, which uses the Euler approximation. Similarly, `init` can be `:crown` or `:stem` to simulate starting with ``1`` or ``2`` lineages. 

Note then that the simulations conditioned on time must input a `Float64` as first (and mandatory) argument while those conditioned on number of species must rater input a `Int64`.

Full documentation
```@docs
sim_gbmb
```

##### Inference

To perform inference under this model we can use:

```julia
r, tv = insane_gbmb(tree,
                    nburn    = 1_000,
                    niter    = 50_000,
                    nthin    = 50, 
                    ofile    = "<directory>",
                    α_prior  = (0.0, 10.0),
                    σλ_prior = (0.05, 0.5),
                    tρ    = Dict("" => 1.0))
```

Here, the prior for the drift `α_prior` is a Normal distribution, with the first element representing the mean and the second the standard deviation, and the prior for the diffusion of speciation is an Inverse Gamma. 

Full documentation
```@docs
insane_gbmb
```

#### Birth-death diffusion process with constant extinction (``\mu(t) = \mu``)

Here we assume that there is extinction, but that is constant across lineages and time.

##### Simulations

As with the pure-birth diffusion, we can simulate conditioned on reaching some number of lineages or some time. Now, however, we need to specify the extinction rate as well, for instance to generate a tree with ``20`` species and some parameter values:

```julia
 sim_gbmce(20, λ0 = 1.0, α = 0.0, σλ = 0.1, μ = 0.7)
```

Full documentation
```@docs
sim_gbmce
```

##### Inference

To perform inference we can use
```julia
r, tv = insane_gbmce(tree,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     α_prior  = (0.0, 10.0),
                     σλ_prior = (0.05, 0.5),
                     μ_prior  = (1.0, 1.0),
                     tρ       = Dict("" => 1.0))
```
where we now specify the Gamma prior `μ_prior` for the extinction.

Full documentation
```@docs
insane_gbmce
```


#### Birth-death diffusion process with constant turnover (``\mu(t) = \epsilon \lambda(t)``)

One can also assume that turnover, ``\epsilon``, that is, the ration of extinction over speciation (_i.e._, ``\frac{\mu}{\lambda}``), is constant.

##### Simulations

As with the other BDD models, we can simulate both conditioned on the number of species or in time. For instance, here we condition on time with a turnover rate of ``0.4``.

```julia
 sim_gbmct(10.0, λ0 = 1.0, α = 0.0, σλ = 0.1, ϵ = 0.4)
```


Full documentation
```@docs
sim_gbmct
```

##### Inference

Performing inference can then be done using:
```julia
r, tv = insane_gbmct(tree,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     α_prior  = (0.0, 10.0),
                     σλ_prior = (0.05, 0.5),
                     ϵ_prior  = (0.0, 100.0),
                     tρ       = Dict("" => 1.0))
```
where we now specify a Uniform prior `ϵ_prior` for the turnover. Also, INSANE resorts to simple MH updates for ``\epsilon``, so an increased burnin is recommended. 

Full documentation
```@docs
insane_gbmct
```

#### Birth-death diffusion process

Finally, we arrive at the most complex model, where extinction also follows a GBM, that is,
```math
d(\text{ln}(\mu(t)) = \sigma_{\mu} d W(t),
```
where ``\sigma_{\mu}`` is the diffusion for extinction.

Note that this is unidentifiable, unless we specify strong priors on the extinction diffusion coefficient ``\sigma_{\mu}``, and still do not do well in recovering simulated values, at least for small trees.


##### Simulations

Again, one can speify a time or a number of lineages to simulate under the BDD. For example, for ``50`` species

```julia
 sim_gbmbd(50, λ0 = 1.0, μ0 = 1.0, α = 0.0, σλ = 0.1, σμ = 0.1)
```

Full documentation
```@docs
sim_gbmbd
```

##### Inference

One can perform inference using:
```julia
r, tv = insane_gbmbd(tree,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     λa_prior = (0.0, 100.0),
                     μa_prior = (0.0, 100.0),
                     α_prior  = (0.0, 10.0),
                     σλ_prior = (0.5, 0.05),
                     σμ_prior = (3.0, 0.5),
                     tρ       = Dict("" => 1.0))
```
where we have uniform priors on the initial (root) values for ``\lambda_0`` and ``\mu_0``, and Inverse Gamma prior for sigma. Usually, an informative prior such as `σμ_prior = (3.0, 0.5)` or more is needed.

Full documentation
```@docs
insane_gbmbd
```


#### Birth-death diffusion process with informed extinction

One can also define a branch-specific fixed extinction function and make inference on the GBM speciation. For example, in Quintero et al. (in prep), we estimated extinction rates based on fossil occurrences and then used this as input.

##### Inference

To make inference under this model, we need two vectors where each element links to a given branch in the reconstructed tree, in the same order as the phylogenetic tree. Each element of these vectors are themselves a vector of type `Float64`, and correspond to the time at which extinction is recovered in time backward order (_e.g._, [1.1, 1.0, ..., 0.1, 0.08, 0.0]) and to the extinction at those times. INSANE uses a linear approximation function between the input sampled points.

To make inference under this model, we now input these vectors. Lets call them time_vector (of type `Vector{Vector{Float64}}`) and extinction_vector (also of type `Vector{Vector{Float64}}`), then:
```julia
r, tv = insane_gbmbd(tree,
                     time_vector,
                     extinction_vector,
                     nburn    = 1_000,
                     niter    = 50_000,
                     nthin    = 50, 
                     ofile    = "<directory>",
                     λa_prior = (0.0, 100.0),
                     μa_prior = (0.0, 100.0),
                     α_prior  = (0.0, 10.0),
                     σλ_prior = (0.5, 0.05),
                     σμ_prior = (3.0, 0.5),
                     tρ       = Dict("" => 1.0))
```

The only thing one has to make sure is that the order of `time_vector` and `extinction_vector` correspond to the order of the branches in `tree`, which can be a bit involved. However, this can be achieved by understanding that the trees are recursive and are ordered either by left or right branch, and using some of the tree utility functions.

For example, one can use the `subclade` function, to extract a subclade given a set of tip labels. One can also use the `make_idf` function, which splits the tree into a vector of branches in the same order that INSANE will process the tree. This way one can see the ordering of branches, and then order the `time_vector` and the `extinction_vector` in the same order.


## Insane tree data access and processing functions

### Simple tree information

Simple information about the tree, such as the number tips, the number of extinct tips, the tree height (duration of the tree) and the tree length (sum of all branch lengths) can be performed using, respectively:
```julia
ntips(tree)
ntipsextinct(tree)
treeheight(tree)
treelength(tree)
```

Julia makes it simple to look at statistics across a vector of trees. For example, using the package `Statistics`, we can estimate the average number of extinct species on tree vector `tv` by simply:
```julia
using Statistics 

mean(ntipsextinct, tv)
```
See the documentation of `mean` for more details, but basically, `mean`, and many other functions in Julia allow to perform a undefined function on each element before calculating the mean. In this case, we are estimating the number of extinct tips in each tree in `tv`, and then averaging over them.

For labelled trees, one can extract the tip labels using `tiplabels`. Moreover we can create subclades based on a vector of tips, where the subclade will be the minimum tree that has all the tips. For example, for a vector `tip_vector` holding Strings that correspond to the tip labels in the tree of type sT_label, we can use
```julia
subclade(tree, tip_vector)
```
However, most time we want to extract subclades of other types of trees. Since these do not hold label information but should be ordered in the same order as the `sT_label` tree, one has to use both. Note, if you change the order of either the `sT_label` tree or the single of vector of trees of other types, this will not work. Thus, if we want a the subclades that have the tips in `tip_vector`, we can use
```julia
subclade(tv, tree, tip_vector, true)
```
where the last argument states if returning the stem or crown tree.


One can also estimate the Lineage Through Time (LTT), or, for a data augmented tree (or a vector of trees), the Diversity Through Time (DTT) using 
```julia
ltt(tree)
```
To be clear, the LTT is usually used to describe the accumulation of reconstructed lineages (those that have been sampled) while DTT is used to describe estimated diversity (sampled **and unsampled** lineages). Thus, the result is either LTT or DTT simply depending on the tree you use as input.

We can also estimate the ltt for a tree vector `tv` using
```julia
ltt(tv)
```

### BDD information

#### Estimating posterior average rates along the tree

Of particular interest is the estimation of posterior average rates along the reconstructed tree. Note that since the data augmented (unsampled) lineages change between different iterations of the algorithm, we obtain lineage-specific instantaneous rate distributions only for the reconstructed (observed) tree (the tree we used as input). Consequently, we first need to remove the data augmented lineages from all the trees in the posterior tree vector:
```julia
tv0 = remove_unsampled(tv)
```

We can then estimate the average tree using
```julia
tm = imean(tv0)
```

We can also estimate any quantile tree, for instance, for the ``0.25`` quantile tree:
```julia
t025 = iquantile(tv0, 0.25)
```

Clearly, these resulting trees can then be further scrutinized as with any other tree in INSANE.


#### Other data access and averaging functions

If one wants to obtain the range (_i.e._, extrema) of the output of function `f` on `tree`, for example, the maximum and minimum speciation rates:
```julia
irange(tree, b)
```

If one wants to sample, recursively, some function at regular intervals along a tree, one can use `sample`. For example if we want to sample speciation rates every ``0.1`` time units, we can use
```julia
sample(tv, b, 0.1)
```
Note that here we are sampling along each branch of the tree in recursive order, not sampling across lineages through time. If we would like to extract an array across lineages in a given tree of the output of function `f`, we would use `time_rate`. For example, if we want the cross-lineage extinction rates of a tree of type `iTbd` sampled every ``0.5`` time units, we would use
```julia
time_rate(tv, d, 0.5)
```
which returns a vector of vectors, where each element is a time holding the rates (in this case extinction rates) of all contemporary lineages at that time.

Finally, a convenience wrapper to extract information recursively from a tree is `trextract`. For example, if we want all branch lengths for a tree, we can use
```julia
trextract(tree, e)
```

Below are some functions to obtain data from trees.

Full documentation
```@docs
tiplabels
ntips
ntipsalive
ntipsextinct
treeheight
treelength
ltt
iscrowntree
irange
time_rate
trextract
subclade
lλ
lμ
e 
```

## Insane tree manipulation functions

Two important manipulation functions are, first to be able to remove extinct lineages, which can be performed on a tree or a tree vector using
```julia
remove_extinct(tree)
```

Similarly, as shown above, one can remove the unsampled lineages (all the data augmented lineages) on a single or vector of trees using
```julia
remove_unsampled(tree)
```
Note that these two are different. First, when performing simulations, the tree is not `fixed`, which means that if you run `remove_unsampled`, you will remove the tree. You would have to fix the tree before, which can be done using `fixtree!(tree)`. Also, if sampling fraction is not $1$, `remove_unsampled` will also remove lineages alive that were not sampled, while `remove_extinct` will only remove those lineages extinct.


Full documentation
```@docs
reorder!
rm_stem!
fixtree!
remove_extinct
remove_unsampled
```

## Insane plots

Tapestree holds recipes to plot phylogenetic trees, results and perform aggregating functions. The output will rely in the number and type of arguments you provide.

First, however, we need to load the `Plots` package
```julia
using Plots
```

### Insane plots across all tree types

#### Simple tree plot

The most basic plot function simply plots the tree:
```julia
plot(tree)
```

One can reorder the tree according to balance (have one daughter always have the largest number of tips) by using `reorder!(tree)`, which orders in place the tree and helps in visualization.



If one has a tree vector, we could, for example, sample 4 of them at random and simply plot them together using

```julia
ti = rand(tv,4)

p0 = plot(ti[1])
p1 = plot(ti[2])
p2 = plot(ti[3])
p3 = plot(ti[4])

plot(p0, p1, p2, p3)
```
If the tree is of type `sT_label`, labels will be shown automatically, but you can toggle this off with `showlabels = false`.

One can also plot the tree radially (as a fan) using
```julia
plot(tree, type = :radial)
```

#### LTT and DTT plots

We can plot the LTT or DTT by using the `ltt` result (of type `Ltt`, check (check [ Simple tree information](@ref))) as input:

```julia
plot(ltt(tree), linewidth = 2.0)
```

Moreover, if we input a vector of `Ltt` we will plot each LTT individually, or, better, if we add a decimal number argument, it will use it as sampling frequency through time and return the mean and desired quantiles of lineage or diversity through time using the arguments `q0` and `q1` (by default `q0 = [0.025, 0.975]` and `q1 = [0.25, 0.75]`):

```julia
lttv = ltt(tv)
plot(lttv, 0.1)
```


### Insane BDD plots

These plotting functions are specific to BDD type trees (_i.e._, of `iT` supertype).

#### Plot rates on the phylogram

To "paint" the tree with the instantaneous lineage-specific rates of speciation ``\lambda(t)``, we can use:
```julia
plot(tree, b)
```
Here, `b`, which stands for "birth rates" is a convenience wrapper around `exp.(lλ(tree))`: it extracts the log-speciation vector from a give `iT` tree using `lλ`, and then returns the exponential.

This plotting function also allows to plot the death rates (only where extinction is also a diffusion, _i.e._, `iTbd`) using
```julia
plot(tree, d)
```
where `d`, which stands for "death rates" is a wrapper around `exp.(lμ(tree))`.

In general, this plotting recipe receives a tree and a function that is applied recursively to paint the tree. Thus, we can use any custom made function that extracts information from the tree. Some predefined ones are:
* `lb`: log-(speciation) birth rates
* `ld`: log-(extinction) death rates 
* `t`: turnover (extinction/speciation) rates
* `lt`: log turnover rates
* `nd`: net diversification (speciation - extinction) rates

We can also plot these trees radially using the `type = :radial`.

#### Plot the underlying rates along a tree

To plot how rates evolve across time, that is, to plot the rates in the y axis, one can use the `type = :rates` argument:
```julia
plot(tree, b, type = :rates)
```

Similarly, we can plot the average for a tree (or other aggregating function as median, geometric mean, _etc._) and custom quantiles (as in [LTT and DTT plots](@ref)) for a given tree by adding a decimal number argument representing the sampling frequency through time.
```julia
plot(tree, b, 0.1)
```
To change the aggregating function, we modify the function `t_af` (by default `t_af = mean`), to the desired one.

#### Plot the rates across tree vectors

Often we would like to plot the average rates across a series of data augmented trees. This can be done by adding a decimal number argument (and, using a tree vector as input).
For instance, to estimate average speciation rates (using wrapping function `b`) through time across tree vector `tv`, every ``0.1`` time units, we use:
```julia
plot(tv, b, 0.1)
```

We can choose the function to aggregate rates across lineages for each single tree using `_af` (by default `t_af = mean`, and then to aggregate these tree averages using the `tv_af` function (by default `tv_af = x -> quantile(x, 0.5)`, that is, the median).

