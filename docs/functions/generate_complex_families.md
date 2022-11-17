---
sort: 5
---

# generate_complex_families

```python
help(recombulatorx.testing.generate_complex_families)
```
```text
Help on function generate_complex_families in module recombulatorx.testing:

generate_complex_families(n_fam_I: int, n_fam_II: int, recombination_rates, mutation_rates)
    Generate randomly simulated individuals from a given number of type I and II families.
    
    For each family, generate a random number of children (from a Poisson distribution with a
    mean of 2.0), with a minimum of 1 or 2 for type I or II families, respectively. A random 
    number of children (from a binomial distribution with a mean of 0.5) are set to be male and 
    the rest female. For each daughter, a father and his paternal haplotype are generated randomly.
    
    n_fam_I: int
        number of type I families
    n_fam_II: int
        number of type II families
    recombination_rates: array of shape (n_markers - 1, )
        recombination rates between adjacent markers used to generate children
    mutation_rates: array of shape (n_markers, )
        mutation rates of markers used to generate children
    
    Returns:
        an iterator to the generated individuals
    ```
