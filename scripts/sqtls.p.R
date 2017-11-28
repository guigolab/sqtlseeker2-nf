#!/usr/bin/env Rscript

#### Run FDR after permuted and nominal passes

## 1. Load libraries and arguments

library(optparse)
library(sQTLseekeR2)
library(data.table)

option_list = list(
  make_option(c("-n", "--nominal"), type = "character",
              help = "All nominal tests", metavar = "FILE"),
  make_option(c("-p", "--permuted"), type = "character",
              help = "All permuted tests", metavar = "FILE"),
  make_option(c("-f", "--fdr"), type = "numeric", help = "FDR level", 
              metavar = "NUMERIC", default = 0.05),
  make_option(c("-t", "--type_fdr"), type = "character", help = "FDR method (BH or qvalue)", 
              metavar = "CHARACTER", default = "BH"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

fdr.level <- opt$fdr
fdr.method <- opt$type_fdr
nominal.f <- opt$nominal
permuted.f <- opt$permuted
output.f <- opt$output

## 2. Input file: all tests  (nominal and permuted)

res.nominal.df <- as.data.frame(fread(nominal.f, header = TRUE, sep = "\t"), stringsAsFactors = FALSE)    # All tests
res.permuted.df <- as.data.frame(fread(permuted.f, header = TRUE, sep = "\t"), stringsAsFactors = FALSE)  # All tests
  
## 3. Get significant sQTLs and output result
  
sqtls.df <- sqtls.p(res.nominal.df, res.permuted.df, method = fdr.method, FDR = fdr.level)
write.table(sqtls.df, file = output.f, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

##### END
