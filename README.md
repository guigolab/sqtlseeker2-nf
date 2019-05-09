# sqtlseeker2-nf

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.0-blue.svg)](http://nextflow.io)
[![Build Status](https://travis-ci.org/dgarrimar/sqtlseeker2-nf.svg?branch=master)](https://travis-ci.org/dgarrimar/sqtlseeker2-nf)

A pipeline for splicing quantitative trait loci (sQTL) mapping written in the Nextflow DSL.

The pipeline performs the following analysis steps:

* Index the genotype file
* Preprocess the transcript expression data
* Test for association between splicing ratios and genetic variants in *cis* (nominal pass)
* Obtain an empirical P-value for each phenotype (permutation pass, optional)
* Control for multiple testing 

For details on each step, please read [sQTLseekeR2](https://github.com/dgarrimar/sQTLseekeR2) documentation.

## Quickstart

Install nextflow with the following command:
```
curl -fsSL get.nextflow.io | bash
```

Pull the docker image:
```
docker pull dgarrimar/sqtlseeker2-nf@sha256:9ddae31aaf8f70f02cd24d3447f4fa3517494da87e12674484d25fe7cf3dc16b
```

Launch the test pipeline with the following command:
```
./nextflow run dgarrimar/sqtseeker2-nf
```

## Pipeline usage

Launching the pipeline with the `--help ` parameter shows the help message:

```
nextflow run sqtlseeker2-nf --help
```

```
N E X T F L O W  ~  version 0.27.2
Launching `sqtlseeker2.nf` [admiring_lichterman] - revision: 28c86caf1c

sqtlseeker2-nf ~ A pipeline for splicing QTL mapping
----------------------------------------------------
Run sQTLseekeR2 on a set of data.

Usage: 
    sqtlseeker2-nf [options]

Options:
--genotype GENOTYPE_FILE    the genotype file
--trexp EXPRESSION_FILE     the transcript expression file
--metadata METADATA_FILE    the metadata file
--genes GENES_FILE          the gene location file
--dir DIRECTORY             the output directory
--mode MODE                 the run mode: nominal or permuted (default: nominal)
--covariates COVARIATES     include covariates in the model (default: false)
--fdr FDR                   false discovery rate level (default: 0.05)
--min_md MIN_MD             minimum effect size reported (default: 0.05)
--svqtl SVQTLS              report svQTLs (default: false)

Additional parameters for mode = nominal:
--ld LD                     threshold for LD-based variant clustering (default: 0, no clustering)
--kn KN                     number of genes per batch in nominal pass (default: 10)

Additional parameters for mode = permuted:
--kp KP                     number of genes per batch in permuted pass (default: 10)
--max_perm MAX_PERM         maximum number of permutations (default: 1000)
```

## Input files and format

`sqtlseeker2-nf` takes as input files the following:

* **Genotype file.**
Contains the genotype of each sample, coded as follows: 0 for REF/REF, 1 for REF/ALT, 2 for ALT/ALT, -1 for missing value.
The first four columns should be: `chr`, `start`, `end` and `snpId`. This file needs to be sorted by coordinate.

* **Transcript expression file.**
Contains the expression of each transcript in each sample (e.g. read counts, RPKM, TPM)
It is not recommended to use transformed (log, quantile, or any non-linear transformation) expression.
Column `trId` and `geneId`, corresponding to the transcript and gene IDs are required. 

* **Metadata file.** Contains the covariate information for each sample. 
In addition, it defines the groups or conditions for which sQTL mapping will be performed.
The first columns should be: `indId`, `sampleId`, `group`, followed by the covariates.
This file defines which samples will be tested.

* **Gene location file.**
Contains the location of each gene. Columns `chr`, `start`, `end` and `geneId` are required. 
This file defines which genes will be tested.

Example [data](data) is available for the test run.

## Pipeline results

sQTL mapping results are saved into the folder specified with the `--dir` parameter. By default it is the `result` directory within the current working folder.

Output files are organinzed into subfolders corresponding to the different `groups` specified in the metadata file: 

```
result
└── groups
    ├── group1                            
    │   ├── all-tests.nominal.tsv          
    │   ├── all-tests.permuted.tsv         
    │   ├── sqtls-${level}fdr.nominal.tsv      
    │   └── sqtls-${level}fdr.permuted.tsv     
    ├── group2
   ...
```

Note: if only a nominal pass was run, files `*.permuted.tsv` will not be present.

Output files contain the following information:

`all-tests.nominal.tsv`

* geneId: gene name	
* snpId: variant name
* F: test statistic
* nb.groups: number of genotype groups
* md: maximum difference in relative expression between genotype groups (sQTL effect size)
* tr.first/tr.second: the transcript IDs of the two transcripts that change the most, in opposite directions
* info: number of individuals in each genotype group, including missing values (-1,0,1,2)
* pv: nominal P-value

if `--svqtl true`
* F.svQTL: svQTL test statistic
* nb.perms.svQTL: number of permutations for svQTL test
* pv.svQTL: svQTL nominal P-value 

if `--ld ${r2}`
* LD: other snps in linkage disequilibrium with *snpId* above a given r<sup>2</sup> threshold > 0

`sqtls-${level}fdr.nominal.tsv` (in addition to the previous)

* fdr: false discovery rate (computed across all nominal tests)
* fdr.svQTL: svQTL FDR

`all-tests.permuted.tsv`

* geneId: gene name
* variants.cis: number of variants tested in *cis*
* LD: median linkage disequilibrium in the region (r<sup>2</sup>)
* best.snp: ID of the top variant
* best.nominal.pv: P-value of the top variant
* shape1: first parameter value of the fitted beta distribution
* shape2: second parameter value of the fitted beta distribution (effective number of independent tests in the region)
* nb.perm: number of permutations
* pv.emp.perm: empirical P-value, computed based on permutations
* pv.emp.beta: empirical P-value, computed based on the fitted beta distribution
* runtime: run time in minutes

`sqtls-${level}fdr.nominal.tsv` (in addition to the previous)

* fdr: false discovery rate (computed across empirical P-values)
* p_tn: gene-level threshold for nominal P-values


## Requirements

`sqtlseeker2-nf` is configured to run using the [Docker](https://www.docker.com/) container engine by default. See the included 
[Dockerfile](docker/Dockerfile) for the configuration details. [Singularity](https://www.sylabs.io/singularity/) is also 
supported, but changes in the Nextflow [configuration](nextflow.config) are required.

In order to run the pipeline with Docker the following dependencies have to be met:

* Java 7/8
* [Nextflow](https://www.nextflow.io) 0.27.0 (or higher)
* [Docker](https://www.docker.com/) 1.12.0 (or higher) 

If you use Singularity:

* [Singularity](https://www.sylabs.io/singularity/) 2.5.0 (or higher)

The pipeline can also be used without Docker/Singularity by installing [sQTLseekeR2](https://github.com/dgarrimar/sQTLseekeR2) on your system.
