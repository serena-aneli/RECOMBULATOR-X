# RECOMBULATOR-X

<p align="center">
  <img align="left" width="350" height="350" src="docs/assets/images/LOGO.png">
</p>
<br/>
<br/>
<br/>
RECOMBULATOR-X is a *Python module* and a *command line tool* for computing the recombination rates between short tandem repeats (STRs) markers along the X chromosome staring from pedigree data.
<br/>
<br/>
<br/>
<br/>
:octocat: [GitHub site](https://github.com/serena-aneli/RECOMBULATOR-X)

:page_facing_up: Please cite [Paper]()

<br/>
<br/>
<br/>
<br/>

---


## Installation

You can install RECOMBULATOR-X via **pip**

```bash
pip install RECOMBULATOR-X
```

---

## Overview

STRs located on the X chromosome are a valuable resource for solving complex kinship cases in forensic genetics thanks to their peculiar inheritance mode. At the same time, the usage of multiple markers linked along the same chromsome, while increasing the evidential weight, also requires proper consideration of the recombination rates between markers in the biostatistical evaluation of kinship.

For more details on X-STR kinship analyses in forensic see [Gomes et al., 2020](https://www.frontiersin.org/articles/10.3389/fgene.2020.00926/full) and [Tillmar et al., 2017](https://www.sciencedirect.com/science/article/pii/S1872497317301126?via%3Dihub).

In the case of forensic X-STRs, recombination rates have been either inferred from population samples through high-density multi-point single nucleotide polymorphism (SNP) data or directly estimated in large pedigree-based studies.

### Statistical framework

Recombination between X-chromosomal markers only happens in female meiosis. This entails that only females can provide information on recombination events, while haploid males can be used to phase their mother/offspring. Ideal linkage-informative families in pedigree studies are therefore three-generation families, including maternal grandfather, mother and one or more sons (“Type I”), or two-generation families of mother and two or more sons (“Type II”). The main statistical approach for the estimation of recombination rates from pedigrees computes the likelihood of kinship by taking into account all possible recombinations within the maternal haplotype, thus resorting to the exponential complexity of the underlying algorithm (see [Nothnagel et al., 2012](https://www.sciencedirect.com/science/article/pii/S1872497312000713?via%3Dihub) for a thorough description of the likelihood computation). Despite a computational update in C++ allowing multi-core parallelization, this approach is expected to be unsuitable when panels of more than 15 X-STRs are considered ([Diegoli et al., 2016](https://www.sciencedirect.com/science/article/pii/S1872497316301247?via%3Dihub)).

We developed RECOMBULATOR-X to overcome this issue. Built upon the same statistical framework of the previous work (Nothnagel et al., 2012), RECOMBULATOR-X uses a new computational strategy to infer recombination rates for X-STRs, while taking also the probability of mutation into account. 

*more on the new computational strategy*

---

### Features

- *Inference of recombination rate faster for many markers*
- Mutation rates can estimated for each marker separately or as a unique parameter 
- Simulation of pedigrees typed with STRs

---

## Usage

### File formats

RECOMBULATOR-X uses the PED files based on [PLINK](https://www.cog-genomics.org/plink/) pedigree files as input. The PED file format stores sample pedigree information (i.e., the familial relationships between samples) and the genotypes.
In particular, the first 6 mandatory columns contain: 

* Family ID
* Individual ID
* Paternal ID
* Maternal ID
* Sex
* Phenotype

The "Sex" field may be coded as: 1=male/2=female; XY=male/XX=female; M=male/F=female; MALE=male/FEMALE=female. If you are using STR markers, you can store within this column the Amelogenin marker. 

The "Phenotype" field refers to the use of PED files in medical research. In non-medical application, it may be -9 which means "unknown".  ATTENZIONE PERCHE' NELLA TABELLA E' -1
        
From the 7th column on, there are the markers genotypes (two columns for a genetic marker, each of the two storing an allele). In case of STRs, the columns contain  numbers, which correspond to the STR repeats or "0" when missing. 

#### Example

Here is a family (each row is an individual):

| **FID** |   **IID**   |   **PAT**   | **MAT** | **SEX** | **PHENO** | **STR-A1** | **STR-A2** | **STR-A1** | **STR-A2** | **STR-A1** | **STR-A2** |
|:-------:|:-----------:|:-----------:|:-------:|:-------:|:---------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|
| FAM_I | GRANDFATHER | 0           | 0       | 1       | -1        | 12         | 0          | 29         | 0          | 39         | 0          |
| FAM_I | MOTHER      | GRANDFATHER | 0       | 2       | -1        | 12         | 16         | 27         | 29         | 34         | 39         |
| FAM_I | SON_1       | 0           | MOTHER  | 1       | -1        | 12         | 0          | 29         | 0          | 34         | 0          |
| FAM_I | FATHER_1    | 0           | 0       | 1       | -1        | 14         | 0          | 21         | 0          | 37         | 0          |
| FAM_I | DAUGHTER_1  | FATHER_1    | MOTHER  | 2       | -1        | 14         | 16         | 21         | 27         | 34         | 37         |
| FAM_I | FATHER_2    | 0           | 0       | 1       | -1        | 18         | 0          | 25         | 0          | 36         | 0          |
| FAM_I | DAUGHTER_2  | FATHER_2    | MOTHER  | 2       | -1        | 12         | 18         | 25         | 29         | 36         | 39         |

---


- pacchetto python
- tool linea di comando







You can use the [editor on GitHub](https://github.com/bioserenina/bioserenina.github.io/edit/main/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/bioserenina/bioserenina.github.io/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
