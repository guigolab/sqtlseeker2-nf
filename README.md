# sqtlseeker2.int-nf

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.0-blue.svg)](http://nextflow.io)

A pipeline for condition-biased splicing quantitative trait loci (sQTL) mapping, modified from sqtlseeker2-nf.

The pipeline performs the following analysis steps:

* Index the genotype file
* Preprocess the transcript expression data
* Test for association: splicing ratios ~ genotype + condition + genotype:condition
* Control for multiple testing via [eigenMT](https://github.com/joed3/eigenMT)

In contrast to `sqtlseeker2-nf`, the model assessed here involves a second main effect (a 2-level factor such as gender, disease status, etc.)
in addition to the genotype, and the interaction between the two of them.

For details on each step, please read [sQTLseekeR2.int](https://github.com/guigolab/sQTLseekeR2/tree/interaction) documentation.

The pipeline uses [Nextflow](http://www.nextflow.io) as the execution backend. Please check [Nextflow documentation](http://www.nextflow.io/docs/latest/index.html) for more information.
