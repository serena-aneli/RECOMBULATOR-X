---
sort: 1
---

# ped2graph

```python
help(recombulatorx.ped2graph)
```
```text
Help on function ped2graph in module recombulatorx.io:

ped2graph(path)
    This function takes a pedigree file as input and encodes it as graphs.
    
    The pedigree file can be a table in .ped format (headerless and space-separated),
    in .tsv format (with header and tab-separated) or in .xlsx format (experimental).
    The format is inferred by the extension.
    
    The graphs are DiGraphs objects from the networkx module, where each node is an
    individual and edges go from parent to child. The nodes include individual information
    (sex and genotypes) as attributes.
    
    path: str
        path to the pedigree file
    
    Returns: tuple
        a list of tuples (family ID, DiGraphs) for each family ID in the pedigree file and 
        a list of marker names.
    ```   
