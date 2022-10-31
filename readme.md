# RECOMBULATOR-X

<p align="center">
  <img align="left" width="300" height="300"  src="docs/assets/images/LOGO.png">
</p>
<br/>
<br/>
recombulator-x is a Python module and command line tool for computing the recombination rates between short tandem repeats (STRs) markers along the X chromosome starting from pedigree data in forensic genetics.
<br/>
<br/>
<br/>
<br/>

:open_book: [Documentation](https://serena-aneli.github.io/recombulator-x/)

:page_facing_up: Please cite [Paper]()

<br/>
<br/>

recombulator-x is written in Python (3.6) and can be used either as a module or as a command-line tool. 
It is the first open source implementation of the estimation method introduced in [Nothnagel et al., 2012](https://www.sciencedirect.com/science/article/pii/S1872497312000713?via%3Dihub), which is the gold-standard for the estimation of recombination rates for X-chromosomal markers. We designed recombulator-x to solve some practical issues with the original R implementation. Its main advantages are:

* performance: much faster than the original implementation, thanks to algorithmic improvements (dynamic programming)
* open source: full source code and documentation available from github
* input parsing: reads pedigree data in standard (PED) format
* user friendly: easy installation (via pip) and usage with a simple command-line tool

We thank Prof. Michael Nothnagel for kindly sharing the original R implementation with us, which was an important reference for the development.

## :open_book: Documentation
Full documentation is available online at the :open_book: [dedicated website](https://serena-aneli.github.io/recombulator-x/), or in this repository under ```docs```.


## :wrench: Installation

You can install recombulator-x via the **pip** command from the standard PyPI repository:

```bash
pip install recombulator-x
```

## :mortar_board: Overview

STRs located on the X chromosome are a valuable resource for solving complex kinship cases in forensic genetics thanks to their peculiar inheritance mode. At the same time, the usage of multiple markers linked along the same chromsome, while increasing the evidential weight, also requires proper consideration of the recombination rates between markers in the biostatistical evaluation of kinship.

For more details on X-STR kinship analyses in forensic see [Gomes et al., 2020](https://www.frontiersin.org/articles/10.3389/fgene.2020.00926/full) and [Tillmar et al., 2017](https://www.sciencedirect.com/science/article/pii/S1872497317301126?via%3Dihub).

In the case of forensic X-STRs, recombination rates have been either inferred from population samples through high-density multi-point single nucleotide polymorphism (SNP) data or directly estimated in large pedigree-based studies.

The main statistical approach for the estimation of recombination rates from pedigrees computes the likelihood of kinship by taking into account all possible recombinations within the maternal haplotype, thus resorting to the exponential complexity of the underlying algorithm (see [Nothnagel et al., 2012](https://www.sciencedirect.com/science/article/pii/S1872497312000713?via%3Dihub) for a thorough description of the likelihood computation). Despite a computational update in C++ allowing multi-core parallelization, this approach is expected to be unsuitable when panels of more than 15 X-STRs are considered ([Diegoli et al., 2016](https://www.sciencedirect.com/science/article/pii/S1872497316301247?via%3Dihub)).

We developed recombulator-x to overcome this issue. Built upon the same statistical framework of the previous work (Nothnagel et al., 2012), recombulator-x uses a new computational strategy to infer recombination rates for X-STRs, while taking also the probability of mutation into account. 

## :boom: Additional features

- data consistency checks
- informative family extraction and preprocessing
- multiple implementations included
- numba
- Mutation rates can estimated for each marker separately or as a unique parameter 
- Simulation of pedigrees typed with STRs

## :rocket: Benchmark

- DIRE QUI CHE FUNZIONA PER FAMIGLIE DI TIPO I E NON II? far exceed the computational speed of the previous approach at least for type I families NEL BENCHMARK
Performance has been the main focus of recombulator-x: in a test with simulated data of the same size as 
This is due to algorithmic improvement and
Performance in particular is much improved, 

## :computer: Usage

recombulator-x uses the PED files based on [PLINK](https://www.cog-genomics.org/plink/) pedigree files as input. The PED file format stores sample pedigree information (i.e., the familial relationships between samples) and the genotypes. More information on the input file can be found in this repository under ```docs/3_usage.md```.

The program can be used both as a Python module or a command-line tool. A detailed notebook for the Python module can be found [here]. 
