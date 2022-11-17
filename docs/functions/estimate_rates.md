---
sort: 1
---

# estimate_rates

The estimation of recombination and mutation rates can be launched with the function estimate_rates, which estimates recombination and optionally mutation rates from a set of families. The function takes the following parameters:

* the families,
* the initial recombination rate,
* which mutation rate needs to be estimated (no: no mutation rate estimation, one: just one mutation rate for all markers, all: a mutation rate for each marker),
* the type of implementation (the default implementation is the one using dynamic programming).

**Important**: Genetic markers (from the 7th column on) must be provided according to their physical genomic position. Indeed, the algorithm will infer the recombination rate between M1 and M2, M2 and M3 and so on.

```python
help(recombulatorx.estimate_rates)
```

```text
Help on function estimate_rates in module recombulatorx.estimate:

estimate_rates(families, starting_recombination_rates=0.05, starting_mutation_rates=0.001, estimate_mutation_rates='no', implementation=None, optimization_method='L-BFGS-B', maxiter=1000)
    Estimate recombination and optionally mutation rates from a set of families.
    
    Estimation is done by expectation maximization, that is by finding recombination 
    and mutation rates that maximize the likelihood of observing the given families.
    Implemented using scipy.optimize.minimize to find recombination and mutation rates 
    that minimize the negative log likelihood.
    
    families: sequence of families
        the preprocessed genetic data of all the families
    starting_recombination_rates: array of shape (n_markers - 1, )
        starting probability of recombination between adjacent markers as initial 
    starting_mutation_rates: array of shape (n_markers, )
        probability of unit mutation (+1 or -1) for each marker
    estimate_mutation_rates: string
        one of 'no', 'one', 'all'. If 'no' 
    implementation: string
        the name of the implementation to use for the likelihood computation
    optimization_method: string
        the name of the method to be passed to scipy.optimize.minimize
    maxiter: int
        the maximum number of optimizing iterations (passed to scipy.optimize.minimize)
      
    Returns: arrays
        an array of estimated recombination rates if estimate_mutation_rates is 'no',
        a tuple of an array of estimated recombination rates and a scalar mutation rate if estimate_mutation_rates is 'one',
        a tuple of an array of estimated recombination rates and an array of mutation rates if estimate_mutation_rates is 'all'
```
