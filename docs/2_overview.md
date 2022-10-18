---
sort: 2
permalink: /overview
---

# Overview

STRs located on the X chromosome are a valuable resource for solving complex kinship cases in forensic genetics thanks to their peculiar inheritance mode. At the same time, the usage of multiple markers linked along the same chromsome, while increasing the evidential weight, also requires proper consideration of the recombination rates between markers in the biostatistical evaluation of kinship.

For more details on X-STR kinship analyses in forensic see [Gomes et al., 2020](https://www.frontiersin.org/articles/10.3389/fgene.2020.00926/full) and [Tillmar et al., 2017](https://www.sciencedirect.com/science/article/pii/S1872497317301126?via%3Dihub).

In the case of forensic X-STRs, recombination rates have been either inferred from population samples through high-density multi-point single nucleotide polymorphism (SNP) data or directly estimated in large pedigree-based studies.

## Statistical framework

Recombination between X-chromosomal markers only happens in female meiosis. This entails that only females can provide information on recombination events, while haploid males can be used to phase their mother/offspring. Ideal linkage-informative families in pedigree studies are therefore three-generation families, including maternal grandfather, mother and one or more sons (“Type I”), or two-generation families of mother and two or more sons (“Type II”). The main statistical approach for the estimation of recombination rates from pedigrees computes the likelihood of kinship by taking into account all possible recombinations within the maternal haplotype, thus resorting to the exponential complexity of the underlying algorithm (see [Nothnagel et al., 2012](https://www.sciencedirect.com/science/article/pii/S1872497312000713?via%3Dihub) for a thorough description of the likelihood computation). Despite a computational update in C++ allowing multi-core parallelization, this approach is expected to be unsuitable when panels of more than 15 X-STRs are considered ([Diegoli et al., 2016](https://www.sciencedirect.com/science/article/pii/S1872497316301247?via%3Dihub)).

We developed RECOMBULATOR-X to overcome this issue. Built upon the same statistical framework of the previous work (Nothnagel et al., 2012), RECOMBULATOR-X uses a new computational strategy, based on dynamic programming, to infer recombination rates for X-STRs, while taking also the probability of mutation into account. Dynamic programming is an optimization technique used in computer science to solve computational problems with overlapping subproblems, where the solution to each subproblem is computed and then efficiently stored in order to avoid recomputing it during the treatment of other overlapping subproblems. In our case, dynamic programming algorithm for the likelihood estimation is implemented explicitly with nested loops and can compiled with Numba, if available, for a big speedup. 
Besides the fast dynamic programming based implementations, we also added the possibility of using the direct approach for the likelihood computations (a *for loop* as in Nothnagel et al., 2012 and a NumPy-based implementation). 

In detail, RECOMBULATOR-X performs the following steps to infer recombination and mutation rates: 

- computation of the likelihood of a son given his mother's genotypes, inheritance vector, recombination and mutation rate;
- for type I families, we sum the likelihoods computed over all possible inheritance vectors (these are disjoint events, since just one inheritance vector will be true). For type II families, we do not know the mother phasing, thus, we repeat the likelihood computation of step 1 by trying over all possible mother phasings;
- multiplication of the likelihoods of all sons within a family in order to obtain the likelihood of the family;
- multiplication of the likelihoods of all families in order to obtain the likelihood of the entire dataset;
- actually, we are not so interested in the likelihoods themselves, but in their parameters, which are the recombination and mutation rates. Hence, using EM, RECOMBULATOR-X finds the parameters maximizing the overall likelihood.


## Features

- Inference of recombination rate faster for many markers
- Mutation rates can estimated for each marker separately or as a unique parameter 
- Simulation of pedigrees typed with STRs

## Assumptions

- Males must be haploid for all the markers: given that our tool is designed for X-chromosomal markers, males have just one copy of the X-chromosome. 
- For the current version of RECOMBULATOR-X, markers must be short tandem repeats (STR).
- Unit mutations: STR fractionary mutations are not allowed and mutation of more than one repeat are assumed to have zero probability
- Genetic markers on the PED files must be provided according to their physical genomic position.

