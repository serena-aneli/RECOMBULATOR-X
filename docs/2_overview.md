---
sort: 2
permalink: /overview
---

# Overview

STRs located on the X chromosome are a valuable resource for solving complex kinship cases in forensic genetics thanks to their peculiar inheritance mode. At the same time, the usage of multiple markers linked along the same chromsome, while increasing the evidential weight, also requires proper consideration of the recombination rates between markers in the biostatistical evaluation of kinship.

For more details on X-STR kinship analyses in forensic see [Gomes et al., 2020](https://www.frontiersin.org/articles/10.3389/fgene.2020.00926/full) and [Tillmar et al., 2017](https://www.sciencedirect.com/science/article/pii/S1872497317301126?via%3Dihub).

In the case of forensic X-STRs, recombination rates have been either inferred from population samples through high-density multi-point single nucleotide polymorphism (SNP) data or directly estimated in large pedigree-based studies.

### Statistical framework

Recombination between X-chromosomal markers only happens in female meiosis. This entails that only females can provide information on recombination events, while haploid males can be used to phase their mother/offspring. Ideal linkage-informative families in pedigree studies are therefore three-generation families, including maternal grandfather, mother and one or more sons (“Type I”), or two-generation families of mother and two or more sons (“Type II”). The main statistical approach for the estimation of recombination rates from pedigrees computes the likelihood of kinship by taking into account all possible recombinations within the maternal haplotype, thus resorting to the exponential complexity of the underlying algorithm (see [Nothnagel et al., 2012](https://www.sciencedirect.com/science/article/pii/S1872497312000713?via%3Dihub) for a thorough description of the likelihood computation). Despite a computational update in C++ allowing multi-core parallelization, this approach is expected to be unsuitable when panels of more than 15 X-STRs are considered ([Diegoli et al., 2016](https://www.sciencedirect.com/science/article/pii/S1872497316301247?via%3Dihub)).

We developed RECOMBULATOR-X to overcome this issue. Built upon the same statistical framework of the previous work (Nothnagel et al., 2012), RECOMBULATOR-X uses a new computational strategy to infer recombination rates for X-STRs, while taking also the probability of mutation into account. 

*more on the new computational strategy*
