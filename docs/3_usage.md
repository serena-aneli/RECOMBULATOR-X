---
sort: 3
permalink: /usage
---

# Usage

## File formats

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

### Example

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


