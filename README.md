# BellsSubgroupCentrality.jl

Implementatino of Bell's subgroup centrality measure. Allows for the calculation
of key measures of subgroup centrality, including:

* **Overall (O)**
* **Global (G)**
* **Local (L)**
* **Boundary (B)**

See `demo_calculations.ipynb` for a reproduction of Betweenness Centrality as calculated in Appendix B of Bell (2014) using a Dolphin network (Figure 4).

##  Use 


### Requirements

`BellsSubgroupCentrality.jl` requires `Graphs.jl`. The optional `BellsSubgroupCentralityExt` module requires Graphs, CSV, XLSX, DataFrames, and DiscreteGraphAlgorithms. 


### Calculating


The quickest way to access one of the 4 subgroup centralities is to

1. Load a `Graphs.jl` graph
2. Define vertices in a group
3. Use `get_sugroup_centrality` to calculate the type of subgroup centrality desired.

E.g., 

```
get_subgroup_centrality(
    graph,
    TYPE,
    group
)[group]
```

where `TYPE` is any of the following symbols:

* `:boundary` or `:b`,
* `:global` or `:g`
* `:local` or `:l`
* `:overall` or `:o`


##  References

Bell, JR. Subgroup centrality measures. Network Science. 2014;2(2):277-297. [doi:10.1017/nws.2014.15](https://www.cambridge.org/core/journals/network-science/article/abs/subgroup-centrality-measures/875D8E7EBF4E33008CB45C6A417C3C66)