# sqtlseeker2-nf

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.0-blue.svg)](http://nextflow.io)
[![Build Status](https://travis-ci.org/dgarrimar/sqtlseeker2-nf.svg?branch=master)](https://travis-ci.org/dgarrimar/sqtlseeker2-nf)

A pipeline for splicing quantitative trait loci (sQTL) mapping written in the Nextflow DSL.

The pipeline performs the following analyses steps:

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
--min_md MIN_MD             minimum effect size reported (default: 0.05)
--svqtl SVQTLS              report svQTLs (default: false)

Additional parameters for mode = nominal:
--ld LD                     threshold for LD-based variant clustering (default: 0, no clustering)
--kn KN                     number of genes per batch in nominal pass (default: 10)
--fdr FDR                   False Discovery Rate (default: 0.05)

Additional parameters for mode = nominal:
--kp KP                     number of genes per batch in permuted pass (default: 10)
--max_perm MAX_PERM         maximum number of permutations (default: 1000)
--fdr FDR                   False Discovery Rate (default: 0.05)
```

## Input format

sqtlseeker2-nf takes as input files the following:

* Genotype file
* Transcript expression file
* Metadata file
* Gene location file

## Pipeline results

sQTL mapping results are saved into the folder specified with the `--dir` parameter. By default it is the `result` directory within the current working folder.

Output files are organinzed into subfolders corresponding to the different `groups` specified in the metadata file: 

```
result
└── groups
    ├── group1                            
    │   ├── all-tests.nominal.tsv          
    │   ├── all-tests.permuted.tsv         
    │   ├── sqtls-0.05fdr.nominal.tsv      
    │   └── sqtls-0.05fdr.permuted.tsv     
    ├── group2
   ...
```

Note: if only a nominal pass was run, files `*.permuted.tsv` will not be present.

## Requirements

sqtlseeker2-nf is configured to run using the [Docker](https://www.docker.com/) container engine by default. See the included 
[Dockerfile](docker/Dockerfile) for the configuration details. [Singularity](https://www.sylabs.io/singularity/) is also 
supported, but changes in the nextflow configuration [file](nextflow.config) are required.

In order to run the pipeline with Docker the following dependencies have to be met:

* Java 7/8
* [Nextflow](https://www.nextflow.io) 0.27.0 (or higher)
* [Docker](https://www.docker.com/) 1.12.0 (or higher) 

If you use Singularity:

* [Singularity](https://www.sylabs.io/singularity/) 2.5.0 (or higher)

The pipeline can also be used without Docker/Singularity by installing [sQTLseekeR2](https://github.com/dgarrimar/sQTLseekeR2) on your system.
