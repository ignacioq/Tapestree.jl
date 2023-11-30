# INSANE

## Reference

Quintero, I., Lartillot, N., Morlon, H. (in prep). Imbalanced speciation pulses sustain the radiation of mammals. 


## Insane Bayesian data augmentation

INSANE uses Bayesian data augmentation (DA) to perform inference on a number of evolutionary models on phylogenetic trees. As such, performing inference will output posterior samples for the governing parameters as well as _complete_ or _data augmented_ trees, that is, trees that include probable configurations of unobserved variables such as the lineages that went extinct in the past or the underlying (latent) speciation rate. 

## Insane tree and model input/output


### Reading and saving newick trees

Tapestree can read files in the simple newick format using the `read_newick` function:
```jldoctest
using Tapestree
tree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_50.tre"))
ntips(tree)

# output

50
```

```

Note that the tree has type `sT_label`, which stands for simple labelled tree. YOu can check this using
```julia
typeof(tree)
```

Similarly, it can also write trees using `write_newick`
```julia
write_newick(tree, "<directory>")
```

### Reading and saving model output

All models can return and/or save the output by writing directly to a file on the fly and in julia when the model stops. There are two outputs:

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


### Constant rates

#### Constant pure-birth (Yule) process (CPB)

##### Simulations

You can simulate a pure-birth tree using `sim_cpb`. For instance, for a period of ``10`` time units and a speciation rate of ``\lambda = 0.5``:
```julia
tr = sim_cpb(10.0, 0.5)
```

##### Inference

The simplest diversification model assumes no extinction and a constant speciation rate ``\lambda``, also known as, the pure-birth or Yule model. To perform inference on a tree (of type `sT_label`), we can use the `insane_cpb` function (cpb = constant pure-birth).
```julia
r, tv = insane_cpb(tree,
                   nburn  = 500,
                   niter  = 1_000,
                   nthin  = 2,
                   nflush = nthin,
                   ofile  = "<directory>")
```

Note that for this specific CPB model, where the sampling fraction, ``\rho``, is $1$, there are no unobserved components since we assume no extinction and that all species have been sampled. In this case, all the trees in the tree vector will be exactly the same.

In the following example, we now specify a global sampling fraction of `0.8`.
```julia
r, tv = insane_cpb(tree,
                   nburn  = 500,
                   niter  = 1_000,
                   nthin  = 2,
                   nflush = nthin,
                   ofile  = "<directory>",
                   tρ     = Dict("" => 0.8))
```

This time, the DA trees are different from one another since they have data augmented lineages that represent that proportion of species not included in the tree. The position change from tree to tree because we are integrating over their unknown placement. You can check this by plotting the trees (check [Insane plots](@ref)). For instance, to plot the first tree:
```julia
plot(tv[1])
```

##### Full documentation
```@docs
insane_cpb
```


#### Constant birth-death process (CBD)

##### Simulations

##### Inference


### Birth-death diffusion processes

#### Pure birth diffusion (``\mu(t) = 0``)

##### Simulations

##### Inference

#### Birth-death diffusion process with constant extinction (``\mu(t) = \mu``)

##### Simulations

##### Inference

#### Birth-death diffusion process with constant turnover (``\mu(t) = \epsilon \lambda(t)``)

##### Simulations

##### Inference



#### Birth-death diffusion process

##### Simulations

##### Inference



#### Birth-death diffusion process informed by fossils

##### Simulations

##### Inference



## Insane plots

Tapestree holds many recipes to plot phylogenetic trees, model results and aggregates


```julia
using Plots
```

## Insane tree functions

<!-- 

3. How many tips does the tree have? (use the `ntips` function)

```julia
ntips(tr)
```

4. How many of these tips are extinct? (use `ntipsextinct`)

```julia
ntipsextinct(tr)
```

5. Check that the tree results from a process of $10$ time units using the `treeheight` function.

```julia
treeheight(tr)
```

6. Plot the resulting tree using `plot(tr)` (note, you have to load the `Plots` package).

```julia
using Plots

plot(tr)
```

7. Estimate the MLE for the speciation rate (note: here we start with one lineage, so no need to condition on observing the tree). You can obtain the tree length (sum of all branch lengths) using `treelength`.

```julia
mle = (Float64(ntips(tr))-1.0)/treelength(tr)
```


8. Now let's run a Bayesian analysis. Because of internal workings we need to change the type of the tree to `sT_label`, this can be easily done by using the `sT_label` function on the tree. Then, use the `insane_cpb` function to run analysis. The only required argument is the tree, but better to specify other keywords, such as `ofile` that specifies the directory to save the results, and `niter` and `nthin` which specify the number of iterations and the thinning. Finally `λ_prior` specifies the parameters of the Gamma prior on the speciation rate.

```julia
tr = sT_label(tr)

r = insane_cpb(tr,
               niter = 1_000,
               nthin = 2, 
               ofile = homedir()*"/repos/tscience_pcm/qmd/yule")
```

9. Explore the resulting MCMC chain in Tracer and compare the posterior distribution for $\lambda$ with the MLE.


10. Read the tree in the Dropbox folder for session 04 using the funciton `read_newick` ("newick" is a basic representation file for a phylogenetic tree). Note that the type os already `sT_label`.

```julia
tr = read_newick("/Users/quintero/Library/CloudStorage/Dropbox/231109 PCM Julia - Instructor/231109 PCM Julia - participants/session_4/tree50.tre")

typeof(tr)
```

11. Make inference under this tree two times ($1000$ or $2000$ iterations should be more than fine with a thinning of $2$), one assuming that the sampling fraction is $\rho = 1$ and another were $\rho = 0.8$. To specify a global sampling fraction (we assume the sampling proportion is uniformly distributed across the tree), you can use the keyword `tρ`, which requires q dictionary, for instance, to specify $\rho = 0.8$, use `tρ = Dict("" => 0.8)`. _NOTE:_ We could also specify a tip specific sampling fraction by creating a dictionary where each entry is `"<tip label>" => rho_i)`. 

```julia
r = insane_cpb(tr,
               niter = 2_000,
               nthin = 5, 
               ofile = homedir()*"/repos/tscience_pcm/qmd/yule_rho1",
               tρ    = Dict("" => 1.0))

r = insane_cpb(tr,
               niter = 2_000,
               nthin = 5, 
               ofile = homedir()*"/repos/tscience_pcm/qmd/yule_rho0.8",
               tρ    = Dict("" => 0.8))
```


# Constant birth-death process


12. Make inference on the `tree50.tre` from above, but assuming a constant birth-death model and $\rho = 1$, to do this use the `insane_cbd` function. For this we might need more iterations, perhaps $50000$ and sampling every $50$ should be fine. All of insane models return an array with the parameters as well as a vector of all the data augmented trees (so you can use `r, tv = insane_cbd(..)`, and r will hold the MCMC run and tv the data augmented trees). Note that, as with the yule process, insane writes a `.txt` file aside from the `.log` file. This `.txt` file writes the data augmented trees; which you can read using the function `iread`, or write using the function `iwrite`.

```julia
r, tv = insane_cbd(tr,
                   niter = 50_000,
                   nthin = 50, 
                   ofile = homedir()*"/repos/tscience_pcm/qmd/bd",
                   tρ    = Dict("" => 1.0))
```

13. Compare the results from assuming no extinction (Yule process) to assuming a birth-death model.

14. Plot four different data augmented trees

```julia
using Plots

ti = rand(tv,4)

p0 = plot(ti[1])
p1 = plot(ti[2])
p2 = plot(ti[3])
p3 = plot(ti[4])

plot(p0, p1, p2, p3)
```

15. What is the average number of extinct lineages?

```julia
using Statistics 

mean(ntipsextinct, tv)
```

16. What is the average tree length?

```julia
mean(treelength, tv)
```

17. Plot the Lineages Through Time (LTT) for the reconstructed tree and, on top, for 5 randomly selected data augmented tree. You can estimate the LTT for each tree using `ltt(tree)`, which returns an object that can be plotted using `Plots`. Does the reconstructed tree seem like having originated by a birth-death process?

```julia
plot(ltt(tr), linewidth = 2.0)

for ti in rand(tv, 5)
  plot!(ltt(ti), linecolor = :orange)
end

plot!()
```

18. Plot the reconstructed Diversity Through Time across all the distribution of data augmented trees (you can do this by using `plot(ltt(<vector of trees>), 0.1)`, here `0.1` specifies how often to sample the diversity, here, every $0.1$ time units.).

```julia
plot(ltt(tv), 0.1)
```


# Birth-death diffusion (BDD) process


19. Let us now use the same tree `tree50.tre` and make inference under the BDD process without extinction. For this use the function `insane_gbmpb`.

```julia
r, tv = insane_gbmpb(tr,
                     niter = 50_000,
                     nthin = 50, 
                     ofile = homedir()*"/repos/tscience_pcm/qmd/bdd_pb",
                     tρ    = Dict("" => 1.0))
```


20. Plot four of the data augmented trees, note that here, to add colors based on the rates you have to type `plot(tree, b)`, where `b` is birth rates.

```julia
tvi = rand(tv, 4)

p0 = plot(tvi[1], b)
p1 = plot(tvi[2], b)
p2 = plot(tvi[3], b)
p3 = plot(tvi[4], b)

plot(p0, p1, p2, p3)
```

21. Estimate the average speciation rate through time across the whole distribution of trees. For this, use `plot(<tree vector>, b, 0.1)`, where `b` again denotes that we want the birth (speciation) rate, and `0.1` ithat we sample every `0.1` time units.

```julia
plot(tv, b, 0.1)
```

22. Estimate the average posterior rates across all of the data augmented trees and plot them. Note that all the unobserved parts of the DA trees change, so we can only get a posterior distribution of rates on the observed lineages, that is, on the reconstructed tree. So, first we need to remove the unsampled part of the tree, we can do this with the function `remove_unsampled`. We can then estimate the mean speciation rates by using the function `imean`.

```julia
tv0 = remove_unsampled(tv)

tm = imean(tv0)

plot(tm, b)
```

23. Plot again the average speciation rates, but plot them in the y axis, this can be achieved by specifying the keyword `type = :rates` in the plot function.

```julia
plot(tm, b, type = :rates)
```

24. Now make inference under the BDD process assuming constant extinction. For this use the function `insane_gbmce`, compare the parameter results using Tracer, and repeat the plots and analyses from Exercises 20 to 23.

```julia
r, tv = insane_gbmce(tr,
                     niter = 50_000,
                     nthin = 50, 
                     ofile = homedir()*"/repos/tscience_pcm/qmd/bdd_ce",
                     tρ    = Dict("" => 1.0))
```

```julia
tvi = rand(tv, 4)

p0 = plot(tvi[1], b)
p1 = plot(tvi[2], b)
p2 = plot(tvi[3], b)
p3 = plot(tvi[4], b)

plot(p0, p1, p2, p3)
```

```julia
plot(tv, b, 0.1)
```

```julia
tv0 = remove_unsampled(tv)

tm = imean(tv0)

plot(tm, b)
```

```julia
plot(tm, b, type = :rates)
```

25. Finally, let's repeat this analyses but assuming that extinction rates also follow a Geometric Brownian motion. Note, however, that this becomes identifiable unless we specify an information prior for the diffusion of extinction rates, $\sigma_{\mu}$. Priors for the diffusion rate of speciation and extinction rates, $\sigma_{\lambda}$ and $\sigma_{\mu}$, are assumed to come from an Inverse-Gamma distribution (conjugate prior for variances). Thus, we can specify $\sigma_{\mu} = \Gamma^{-1}(3.0, 0.1)$, in the keyword argument `σλ_prior = (3.0, 0.1)`. Note that we might need some extra iterations.

```julia
r, tv = insane_gbmbd(tr,
                     niter = 200_000,
                     nthin = 200, 
                     ofile = homedir()*"/repos/tscience_pcm/qmd/bdd_bd",
                     tρ    = Dict("" => 1.0),
                     σλ_prior = (0.05, 0.05),
                     σμ_prior = (3.0, 0.5))
```

26. Compare these results with the previous ones in Tracer, and now plot both speciation rates and extinction rates as done before. For extinction rates change the `b` for `d` (for death rates).

```julia
using LaTeXStrings

p0 = plot(tv, b, 0.1, yguide = L"\lambda")
p1 = plot(tv, d, 0.1, yguide = L"\mu")

plot(p0, p1)
```


```julia
tv0 = remove_unsampled(tv)

tm = imean(tv0)

p0 = plot(tm, b)
p1 = plot(tm, d)

plot(p0, p1)
```

```julia
p0 = plot(tm, b, type = :rates)
p1 = plot(tm, d, type = :rates)

plot(p0, p1)
```

27. Simulate a crown tree (starting with a speciation event) of $20$ species under the BDD model using `sim_gbmbd` and specifying the parameters of your preference (I recommend you increase the starting extinction rate, $\mu_0$. If the tree is not sampled, increase the `p` argument (e.g., just for this example, you can use `p = 1e10`. To simulate a crown tree, use `start = :crown` rather than `start = :stem` (the default). 

```julia
tr = sim_gbmbd(20, μ0 = 0.5, start = :crown, p = 1e10)
```

28. Make sure both crown lineages (the first two lineages survive). Then save only what the reconstructed tree would look like (remove extinct lineages) as a newick file using `write_newick`.

```julia
while iszero(ntipsalive(tr.d1)) || iszero(ntipsalive(tr.d2))
  tr = sim_gbmbd(20, μ0 = 0.5, start = :crown, p = 1e10)
end

tr0 = remove_extinct(tr)

write_newick(tr0, homedir()*"/repos/tscience_pcm/qmd/simtree")
```
 -->