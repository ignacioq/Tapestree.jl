# INSANE

## Insane Bayesian data augmentation

INSANE uses Bayesian data augmentation (DA) to perform inference on a number of evolutionary models on phylogenetic trees. As such, performing inference will output posterior samples for the governing parameters as well as _complete_ or _data augmented_ trees, that is, trees that include probable configurations of unobserved variables such as the lineages that went extinct in the past or the underlying (latent) speciation rates. 

### References

Constant birth-death (CBD) and birth-death diffusions (BDD): Quintero, I., Lartillot, N., Morlon, H. (2024). Imbalanced speciation pulses sustain the radiation of mammals. _Science_, **384**: 1007-1012. [link](https://doi.org/10.1126/science.adj2793)

Constant fossilized birth-death (CFBD) and fossilized birth-death diffusion (FBDD):
Quintero, I., Andréoletti, J., Silvestro, D., & Morlon, H. (2025). "The rise, decline and fall of clades". BioRxiv  [link](https://doi.org/10.1101/2025.03.20.644316)

Diffused Brownian motion (DBM): Quintero, I. (2025). The diffused evolutionary dynamics of morphological novelty. _Proceedings of the National Academy of Sciences_, U.S.A. **122** (18) e2425573122, [link](https://doi.org/10.1073/pnas.2425573122).


## Insane input and structures

```@contents
Pages = [
  "man/insane/io.md"
]
Depth = 2
```

## Insane models

```@contents
Pages = [
  "man/insane/cbd.md",
  "man/insane/bdd.md",
  "man/insane/cfbd.md",
  "man/insane/fbdd.md",
  "man/insane/dbm.md",
]
Depth = 2
```

## Insane plots and utilities

```@contents
Pages = [
  "man/insane/processing.md"
  "man/insane/iplots.md"
]
Depth = 2
```