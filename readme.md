# RECOMBULATOR-X

<p align="center">
  <img align="left" width="350" height="350" src="docs/assets/images/LOGO.png">
</p>
<br/>
<br/>
<br/>
RECOMBULATOR-X is a *Python module* and a *command line tool* for computing the recombination rates between short tandem repeats (STRs) markers along the X chromosome starting from pedigree data.
<br/>
<br/>
<br/>
<br/>

:octocat: [GitHub](https://github.com/serena-aneli/RECOMBULATOR-X)

:open_book: [Documentation](https://serena-aneli.github.io/RECOMBULATOR-X/)

:page_facing_up: Please cite [Paper]()

<br/>
<br/>
<br/>
<br/>


The program is written in Python (3.6) and can be used as a Python module or a command-line tool. 

We designed RECOMBULATOR-X to solve some issues in the computation of recombination rates for X-chromosomal markers in forensics:

* something
* bla
* bla
* bla


## :open_book: Documentation
Full documentation is available online at the :open_book: [dedicated website](https://serena-aneli.github.io/RECOMBULATOR-X/), or in this repository under ```docs```.


## :wrench: Installation

You can install RECOMBULATOR-X via **pip**

```bash
pip install RECOMBULATOR-X
```


## :mortar_board: Overview

STRs located on the X chromosome are a valuable resource for solving complex kinship cases in forensic genetics thanks to their peculiar inheritance mode. At the same time, the usage of multiple markers linked along the same chromsome, while increasing the evidential weight, also requires proper consideration of the recombination rates between markers in the biostatistical evaluation of kinship.

For more details on X-STR kinship analyses in forensic see [Gomes et al., 2020](https://www.frontiersin.org/articles/10.3389/fgene.2020.00926/full) and [Tillmar et al., 2017](https://www.sciencedirect.com/science/article/pii/S1872497317301126?via%3Dihub).

In the case of forensic X-STRs, recombination rates have been either inferred from population samples through high-density multi-point single nucleotide polymorphism (SNP) data or directly estimated in large pedigree-based studies.

The main statistical approach for the estimation of recombination rates from pedigrees computes the likelihood of kinship by taking into account all possible recombinations within the maternal haplotype, thus resorting to the exponential complexity of the underlying algorithm (see [Nothnagel et al., 2012](https://www.sciencedirect.com/science/article/pii/S1872497312000713?via%3Dihub) for a thorough description of the likelihood computation). Despite a computational update in C++ allowing multi-core parallelization, this approach is expected to be unsuitable when panels of more than 15 X-STRs are considered ([Diegoli et al., 2016](https://www.sciencedirect.com/science/article/pii/S1872497316301247?via%3Dihub)).

We developed RECOMBULATOR-X to overcome this issue. Built upon the same statistical framework of the previous work (Nothnagel et al., 2012), RECOMBULATOR-X uses a new computational strategy to infer recombination rates for X-STRs, while taking also the probability of mutation into account. 

## :boom: Features

- *Inference of recombination rate faster for many markers*
- Mutation rates can estimated for each marker separately or as a unique parameter 
- Simulation of pedigrees typed with STRs

## :rocket: Benchmark

- *Inference of recombination rate faster for many markers*
- Mutation rates can estimated for each marker separately or as a unique parameter 
- Simulation of pedigrees typed with STRs

## :computer: Usage

RECOMBULATOR-X uses the PED files based on [PLINK](https://www.cog-genomics.org/plink/) pedigree files as input. The PED file format stores sample pedigree information (i.e., the familial relationships between samples) and the genotypes. More information on the input file can be found in this repository under ```docs/3_usage.md```.

- pacchetto python
- tool linea di comando
