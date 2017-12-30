#!/usr/bin/env Rscript

#### Run FDR after permuted and nominal passes

## 1. Load libraries and arguments

library(optparse)
library(sQTLseekeR2)
library(data.table)

option_list <- list(
    make_option(c("-n", "--nominal"), type = "character",
                help = "All nominal tests", metavar = "FILE"),
    make_option(c("-p", "--permuted"), type = "character",
                help = "All permuted tests", metavar = "FILE"),
    make_option(c("-f", "--fdr"), type = "numeric", help = "FDR level", 
                metavar = "NUMERIC", default = 0.05),
    make_option(c("-t", "--type_fdr"), type = "character", 
                help = "FDR method (BH or qvalue)", metavar = "CHARACTER", 
                default = "BH"),
    make_option(c("-m", "--md_min"), type = "numeric", 
                help = "sQTLs with MD value below the threshold are not reported", 
                metavar = "NUMERIC", default = 0.05),
    make_option(c("-r", "--rm_svqtl"), action = "store_true", 
                help = "significant sQTLs that are also significant svQTLs are not reported", 
                default = FALSE),
    make_option(c("-o", "--output"), type = "character",
                help = "Output file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

nominal.f <- opt$nominal
permuted.f <- opt$permuted
output.f <- opt$output

if ( is.null(nominal.f) || is.null(permuted.f) || is.null(output.f) ){
    print_help(opt_parser)
    stop("Missing/not found input files", call.= FALSE)
}

## 2. Input file: all tests  (nominal and permuted)

res.nominal.df <- as.data.frame(fread(nominal.f, header = TRUE, sep = "\t"), 
                                stringsAsFactors = FALSE)  # All tests
res.permuted.df <- as.data.frame(fread(permuted.f, header = TRUE, sep = "\t"), 
                                 stringsAsFactors = FALSE)  # All tests
  
## 3. Get significant sQTLs and output result
  
sqtls.df <- sqtls.p(res.nominal.df, res.permuted.df, FDR = opt$fdr,
                    method = opt$type_fdr, md.min = opt$md_min,
                    svQTL.removal = opt$rm_svqtl, FDR.svQTL = opt$fdr)

write.table(sqtls.df, file = output.f, quote = FALSE, row.names = FALSE, 
            col.names = TRUE, sep = "\t")

#### END
