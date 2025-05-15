# Insane input & structures

## Insane tree and model input/output

### Insane trees

INSANE uses different types of trees, but all are recursive structures. All trees have at least three fields: 
  * `d1` which corresponds to the left daughter (either another tree structure or `nothing` if it is a tip).
  * `d2` which corresponds to the right daughter (also either another tree structure or `nothing` if it is not a bifurcation).
  * `e` a decimal number specifying the edge length of the tree.

For all tree output from simulations and inference (all expect labelled trees):
  * `fx` a Boolean stating if tree is fixed or not. That is, if it corresponds to a sampled (`fx = true`) or unsampled (`fx = false`) lineage.

For all trees that allow for extinction:
  * `iμ` a Boolean stating if tree is extinct (`iμ = true`) or not (`iμ = false`).

For fossil trees:
  * `iψ` a Boolean stating if tree is fossil (`iψ = true`) or not (`iψ = false`).


### Insane simple trees (`sT` & `sTf`)

1. `sT_label`: is a simple labelled tree, and is the one required as input to perform inference. 
2. `sTf_label`: is a simple fossil labelled tree, and is the one required as input to perform inference for fossilized birth-death models. 
  * They have the additional field `l` for label.
!!! note 
    All other model specific tree types can be converted to `sT_label` by using the latter as a function: `tree = sT_label(tree_of_other_type)`. For fossil trees as well by using `tree = sTf_label(fossil_tree_of_other_type)`.
3. `sTpb`: is a simple pure-birth tree for the constant pure-birth (Yule) model.
4. `sTbd`: is a simple birth-death tree for the constant birth-death model.
5. `sTfbd`: is a simple fossil birth-death tree for the constant and episodic fossilized birth-death models.


### Insane BDD & FBDD trees (`iT`)

All `iT` trees specify to birth-death diffusion trees with different extinction assumptions and all hold the following extra fields:

  * `dt` a decimal specifying the time step of the GBM discretization.
  * `fdt` a decimal specifying the final time step of the GBM discretization.
  * `lλ` an array specifying the Brownian motion evolution of log-speciation rates.

And have the following concrete types:

1. `iTpb`, `iTce` and `iTct`: are BDD trees with no-extinction (pure-birth), constant extinction and constant turnover, respectively.
2. `iTbd`: the full BDD tree.
3. `iTfbd`: for the fossil BDD tree.
  * These last two have the additional field `lμ`:  an array specifying the Brownian motion evolution of log-extinction rates.

### Insane DBM trees (`sTxs`)

  * `dt` a decimal specifying the time step of the GBM discretization.
  * `fdt` a decimal specifying the final time step of the GBM discretization.
  * `xv` an array specifying the Brownian motion evolution of traits.
  * `lσ2` an array specifying the Geometric Brownian motion evolution of rates.

### Reading and saving newick trees

#### Extant only trees

Tapestree can read files in the simple newick format using the `read_newick` function:
```julia
tree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_5.tre"))
```

!!! note
   This tree has type `sT_label`, which stands for simple labelled tree. You can check this using `typeof(tree)`

#### Trees with fossils

Tapestree reads fossilized birth-death trees in newick format using the `read_newick` function, as above, but specifying as second argument `true`:
```julia
tree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_6.tre"), true)
```

In addition, because of rounding errors, a very recent fossil tip might actually be labelled as an extant tip or vice-versa. To control for the time from which extant species are differentiated from fossil tip ones, the argument `ne` can be used:
```julia
tree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_6.tre"), true, ne = 0.1)
```

Tapestree can also write trees using `write_newick`
```julia
write_newick(tree, "<directory>")
```

!!! note
    Only `*_label` trees (_i.e._, `sT_label` & `sTf_label`) have labels. So if you want to save a DA tree with the original tip labels plus new names to the data augmented trees, you first will have to create a labelled tree from the `DA_tree` and the loaded labelled `tree` and then save it:
    `write_newick(sT_label(DA_tree, tree), "<directory>")`.

### Reading and saving model output

All models can return and/or save the output by writing directly to a file on the fly and within the Julia session when the model finishes. There are two outputs:

1. The governing parameters trace, which is saved as a `.log` file and returned as the first object once the algorithm finished. These can be conveniently read in the Tracer software [https://github.com/beast-dev/tracer/releases/tag/v1.7.2](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).
2. The DA trees as a tree vector which are saved as a `.txt` file (of the same name as the `.log`) and returned as the second object once the algorithm finished. 

The DA trees written in the insane-specific `.txt` file can be read using the `iread()` function, which only needs the specific file directory as input, but also accepts an optional (keyword) argument `ix` that indicates the specific tree iterations to read, as an `OrdinalRange` object. For instance, to read only the first ``50`` trees, one can use `iread("<directory to txt>", ix = 1:50`. To read only the trees every ``10`` iterations from the ``100`` to ``400`` sampled, one can use `ix = 100:10:400`, and so on. This can be helpful to avoid high computation costs of reading very large files.

You can save any individual or vector of insane trees using the `iwrite()` function: `iwrite(trees, "<directory>")`, which can be read by `Tapestree`. For portability, you can also save the DA trees as nexus files where nodes have been annotated with the diffusion vectors using `write_nexus(trees, tree, "<directory>")`, where trees is the tree vector and tree is the original labelled tree (_i.e._, the one you used as input into an insane model).


## Insane inference

### Tree input

All inference functions require a phylogenetic tree of type `sT_label` or `sTf_label`, that is, a simple labelled tree. This is the default object type when using the `read_newick` function. However, when using simulations from models or `iread`, the resulting type is specific to the model (for computational efficiency). One can easily create a tree of type `sT_label` or `sTf_label` from any other tree by using `sT_label(tree)` and `sTf_label(tree)`, respectively. The output can then be used to perform inference.


### Common (keyword) arguments across all insane inference models

* `nburn`: specifies the number of iterations to discard as burn-in.
* `niter`: specifies the number of MCMC iterations. 
* `nthin`: specifies the iteration frequency at which to save the parameters **in the julia session** (_i.e._,`nthin = 2` specifies saving every 2 iterations), 
* `nflush`: specifies the frequency at which to save **to file**. 
* `ofile`: specifies the directory where the results will be written. 
* `tρ`: controls the sampling fraction and receives a `Dictionary` as input, with a `String` key pointing to a `Float64` number (_i.e._, `Dict{String, Float64}`). If the dictionary is of length 1 with an empty string, then the insane sets this as a the global sampling fraction. For example, to set a sampling fraction of `0.6`, one show input `tρ = Dict("" => 0.6)`. Most times, however, sampling fraction is not uniform across the tree, but rather some part so the tree is more heavily sampled than others, to accommodate these variability, you can input a dictionary of the same length as the number of tips in the tree, where the dictionary key string is the tip label pointing to the specific sampling fraction value. For example, for two tips, named `tip_1` and `tip_2`, one could input `tρ = Dict("tip_1" => 0.5, "tip_2" => 0.3)`. Make sure to specify all tips when assigning different sampling fractions across the tips, even ones with `1.0`. 
* `prints`: specifies the number of seconds to refresh the progress meter.
* `survival`: For those models with extinction, `true` or `false` if to condition the likelihood on survival.


Full documentation
```@docs
read_newick
write_newick
iread
iwrite
write_nexus
```
