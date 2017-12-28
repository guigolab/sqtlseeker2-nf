#!/usr/bin/env Rscript

#### Index genotype file

## 1. Load libraries and arguments

library(optparse)
library(sQTLseekeR2)

option_list <- list(
    make_option(c("-g", "--genotype_file"), type = "character",
                help = "012 encoded genotype file", metavar = "GENOTYPES")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## 2. Input file: genotype information

genotype.f <- opt$genotype_file

if ( is.null(genotype.f) ){
    print_help(opt_parser)
    stop("File not found: ", genotype.f , call. = FALSE)
}

## 3. Index the genotype file

index.genotype(genotype.f)

####  END
