#!/usr/bin/env Rscript

#### Run FDR after nominal pass

## 1. Load libraries and arguments

library(optparse)
library(sQTLseekeR2)
library(data.table)

option_list = list(
  make_option(c("-a", "--all_tests"), type = "character",
              help = "All nominal tests", metavar = "FILE"),
  make_option(c("-f", "--fdr"), type = "numeric", help = "FDR level", 
              metavar = "NUMERIC", default = 0.05),
  make_option(c("-o", "--output"), type = "character",
              help = "Output file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## 2. Input file: all tests  

res.df <- as.data.frame(fread(opt$all_tests, header = TRUE, sep = "\t"), stringsAsFactors = F) # All tests

## 3. Get significant sQTLs and output result

sqtls.df <- sqtls(res.df, FDR = opt$fdr, FDR.svQTL = opt$fdr)
write.table(sqtls.df, file = opt$output, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
#### END