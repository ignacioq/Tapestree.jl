# Insane tree data access and processing functions

## Basic tree information

Basic information about the tree, such as the number tips, the number of extinct tips, the number of fossils, the tree height (duration of the tree) and the tree length (sum of all branch lengths) can be performed using, respectively:
```julia
ntips(tree)
ntipsextinct(tree)
nfossils(tree)
treeheight(tree)
treelength(tree)
```

## Tree vector statistics

Julia makes it simple to look at statistics across a vector of trees. For example, using the package `Statistics`, we can estimate the average number of extinct species on tree vector `tv` by simply:
```julia
using Statistics 

mean(ntipsextinct, tv)
```
See the documentation of `mean` for more details, but basically, `mean`, and many other functions in Julia allow to perform a undefined function on each element before calculating the mean. In this case, we are estimating the number of extinct tips in each tree in `tv`, and then averaging over them.


## Tree labels and obtaining subtrees

For labelled trees, one can extract the tip labels using `tiplabels`. Moreover we can create subclades based on a vector of tips, where the subclade will be the minimum tree that has all the tips. For example, for a vector `tip_vector` holding Strings that correspond to the tip labels in the tree of type sT_label, we can use
```julia
subclade(tree, tip_vector)
```
However, most time we want to extract subclades of other types of trees. Since these do not hold label information but should be ordered in the same order as the `sT_label` tree, one has to use both. 

!!! warning
    if you change the order of either the `sT_label` tree or the single of vector of trees of other types, this will not work. 

Thus, if we want a the subclades that have the tips in `tip_vector`, we can use
```julia
subclade(tv, tree, tip_vector, true)
```
where the last argument states if returning the stem or crown tree.


## Lineage and Diversity through time (LTT & DTT)

One can also estimate the Lineage Through Time (LTT), or, for a data augmented tree (or a vector of trees), the Diversity Through Time (DTT) using 
```julia
ltt(tree)
```
To be clear, the LTT is usually used to describe the accumulation of reconstructed lineages (those that have been sampled) while DTT is used to describe estimated diversity (sampled **and unsampled** lineages). Thus, the result is either LTT or DTT simply depending on the tree you use as input.

We can also estimate the ltt for a tree vector `tv` using
```julia
ltt(tv)
```


## Trees with diffusion information (_e.g._, BDD, FBDD, DBM)

### Estimating posterior average rates along the tree

Of particular interest is the estimation of posterior average rates along the reconstructed tree. Since the data augmented (unsampled) lineages change between different iterations of the algorithm, we obtain lineage-specific instantaneous rate distributions only for the reconstructed (observed) part of the trees (the tree we used as input). Consequently, we first need to remove the data augmented lineages from all the trees in the posterior tree vector:
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


### Other data access and averaging functions

If one wants to obtain the range (_i.e._, extrema) of the output of function `f` on `tree`, for example, the maximum and minimum speciation rates:
```julia
irange(tree, b)
```

If one wants to sample, recursively, some function at regular intervals along a tree, one can use `sample`. For example if we want to sample speciation rates every ``0.1`` time units, we can use
```julia
sample(tv, b, 0.1)
```
!!! note 
    Here we are sampling along each branch of the tree in recursive order, not sampling across lineages through time. 

If we would like to extract an array across lineages in a given tree of the output of function `f`, we would use `time_rate`. For example, if we want the cross-lineage extinction rates of a tree of type `iTbd` sampled every ``0.5`` time units, we would use
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

!!! note
    `remove_extinct` and `remove_unsampled` are different. First, when performing simulations, the tree is not `fixed`, which means that if you run `remove_unsampled`, you will remove the tree. You would have to fix the tree before, which can be done using `fixtree!(tree)`. Also, if sampling fraction is not $1$, `remove_unsampled` will also remove lineages alive that were not sampled, while `remove_extinct` will only remove those lineages extinct.

For fossil trees, one can remove all fossils using
```julia
remove_fossils(tree)
```
or make a given tree a fossil by using
```julia
fossilize!(tree)
```
which will only make fossil that specific tree (not the recursive daughters).



Full documentation
```@docs
reorder!
rm_stem!
fixtree!
remove_extinct
remove_unsampled
```


