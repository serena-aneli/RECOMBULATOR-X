---
sort: 3
---

# preprocess_families

This function identifies the informative families for the estimation of recombination rates. 
For recombination, informative subfamilies are either those with:

* a phased mother and at least one son or phased daughter, called type I families,
* an unphased mother and at least two between sons and phased daughters, called type II families.

Notably, females can be phased when their father is available: in this way, they will be virtually transformed into males, thus being allowed to take part to informative families.

```python
help(recombulatorx.preprocess_families)
```

```text
Help on function preprocess_families in module recombulatorx.families:

preprocess_families(family_graphs)
    Preprocess family graphs for rates estimation.
    
    Does the following preprocessing steps:
    - checks the consistency of the parsed family graphs,
    - finds all the informative subfamilies from each family,
    - phases all females using their father when available,
    - packs the inherited maternal haplotypes into a matrix
    
    family_graphs: sequence of tuples (str, networkx.DiGraph)
        sequence of tuples of family IDs and family graphs
    
    Returns: list of ProcessedFamily
        the preprocessed families
```   
