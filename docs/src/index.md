
# Tapestree.jl

**Unraveling the evolutionary tapestry**:

[Tapestree] (https://github.com/ignacioq/Tapestree.jl) is a 
[Julia](http://julialang.org) package of phylogenetic analyses of 
diversification, trait and biogeographic dynamics.


## Package Features

Tapestree currently holds the following phylogenetic models:

### INSANE models using Bayesian Data augmentation (DA)

#### Diversification (birth-death) models

  - "Constant Birth-Death" (CBD): implements a birth-death model of constant speciation rates and extinction rates under a Bayesian data augmentation algorithm, as described in Quintero, I., Lartillot, N., & Morlon, H. (2024). "Imbalanced speciation pulses sustain the radiation of mammals". _Science_, **384** (6699) 1007-1012  [https://doi.org/10.1126/science.adj2793](https://doi.org/10.1126/science.adj2793)

  - "Birth-Death Diffusion" (BDD): implements birth-death models where speciation rates follow a geometric Brownian motion with either no extinction, constant extinction, constant turnover, informed extinction or extinction diffusion under a Bayesian data augmentation algorithm, as described in Quintero, I., Lartillot, N., & Morlon, H. (2024). "Imbalanced speciation pulses sustain the radiation of mammals". _Science_, **384** (6699) 1007-1012  [https://doi.org/10.1126/science.adj2793](https://doi.org/10.1126/science.adj2793)

  - "Constant Fossilized Birth-Death" (CFBD): implements a birth-death model of constant speciation rates and extinction rates with piece-wise constant temporal changes in fossilization rates under a Bayesian data augmentation algorithm, as described in Quintero, I., Andréoletti, J., Silvestro, D., & Morlon, H. (2025). "The rise, decline and fall of clades". BioRxiv  [https://doi.org/10.1101/2025.03.20.644316](https://doi.org/10.1101/2025.03.20.644316)

  - "Fossilized Birth-Death Diffusion" (FBDD): implements a fossilized birth-death model where speciation and extinction rates follow separate geometric Brownian motions under a Bayesian data augmentation algorithm, as described in Quintero, I., Andréoletti, J., Silvestro, D., & Morlon, H. (2025). "The rise, decline and fall of clades". BioRxiv  [https://doi.org/10.1101/2025.03.20.644316](https://doi.org/10.1101/2025.03.20.644316)

  "Occurrence Birth-Death" (OBD): implements a (fossilized or not) birth-death model of constant speciation rates and extinction rates with piece-wise constant temporal changes in fossilization rates and external occurrences under a Bayesian data augmentation algorithm, as described in Andréoletti J., Quintero, I., Morlon, H. (2025) The Occurrence Birth-Death Diffusion Process: Unraveling Diversification Histories with Fossils and Heterogeneous Rates. bioRxiv [https://doi.org/10.1101/2025.08.26.672414](https://doi.org/10.1101/2025.08.26.672414)

  "Occurrence Birth-Death Diffusion" (OBDD): implements a (fossilized or not) birth-death model where speciation and extinction rates follow separate geometric Brownian motions with piece-wise constant temporal changes in fossilization rates and external occurrences under a Bayesian data augmentation algorithm, as described in Andréoletti J., Quintero, I., Morlon, H. (2025) The Occurrence Birth-Death Diffusion Process: Unraveling Diversification Histories with Fossils and Heterogeneous Rates. bioRxiv [https://doi.org/10.1101/2025.08.26.672414](https://doi.org/10.1101/2025.08.26.672414)



#### Diffused Brownian motion (DBM) model

  - DBM implements a trait evolution model where trait and its underlying rates follow separate geometric Brownian motions under a Bayesian data augmentation algorithm, as described in Quintero, I. (2025). The diffused evolutionary dynamics of morphological novelty. _Proceedings of the National Academy of Sciences_, U.S.A. **122** (18) e2425573122, [https://doi.org/10.1073/pnas.2425573122](https://doi.org/10.1073/pnas.2425573122).

### Environmental and State dependent Speciation and Extinction (ESSE) model

  - ESSE implements a joint geographic and environmental dependent model of diversification as described in Quintero, I., Landis, M. J., Jetz, W., & Morlon, H. (2022). The build-up of the present-day tropical diversity of tetrapods. _Proceedings of the National Academy of Sciences_, 2023. **120** (20) e2220672120. [https://doi.org/10.1073/pnas.2220672120](https://doi.org/10.1073/pnas.2220672120)

### Trait and Range Interspecific Biogeographic Evolution (TRIBE) model

  - TRIBE implements a joint model of trait evolution and biogeographic history shaped by biotic interactions under a Bayesian data augmentation algorithm, as described in Quintero, I., & Landis, M. J. (2020). Interdependent phenotypic and biogeographic evolution driven by biotic interactions. _Systematic biology_, **69** (4), 739-755. [https://doi.org/10.1093/sysbio/syz082](https://doi.org/10.1093/sysbio/syz082)


## Quick Start

```@contents
Pages = [
    "installation.md",
]
Depth = 1
```

```@contents
Pages = [
    "quick_start.md",
]
Depth = 4
```


## Manual

```@contents
Pages = [
    "man/insane/contents.md",
    "man/tribe.md",
    "man/esse.md",
]
Depth = 2
```

## Issues

Submit them to [Issues](https://github.com/ignacioq/Tapestree.jl/issues)
