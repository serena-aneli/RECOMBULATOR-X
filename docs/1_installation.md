---
sort: 1
permalink: /installation
---

# Installation

You can install recombulator-x via the **pip** command from the standard PyPI repository:

```bash
pip install recombulator-x
```

# Software dependencies

recombulator-x needs the following Python modules and versions: 

- [numpy](https://numpy.org/)>=1.14
- [pandas](https://pandas.pydata.org/)>=0.23
- [scipy](https://scipy.org/)>=1.0
- [networkx](https://networkx.org/)>=2.0. 

Notably, it is not mandatory to have Numba installed: if [`Numba`](https://numba.pydata.org/) it is not installed, just the implementations which do not use Numba will be available. 
