#!/usr/bin/env Rscript

#### Run FDR after nominal pass

## 1. Load libraries and arguments

library(optparse)
library(sQTLseekeR2)
library(data.table)

option_list <- list(
    make_option(c("-n", "--nominal"), type = "character",
                help = "all nominal tests", metavar = "FILE"),
    make_option(c("-f", "--fdr"), type = "numeric", help = "FDR level", 
                metavar = "NUMERIC", default = 0.05),
    make_option(c("-r", "--rm_svqtl"), action = "store_true", 
                help = "significant sQTLs that are also significant svQTLs are not reported", 
                default = FALSE),
    make_option(c("-m", "--md_min"), type = "numeric", 
                help = "sQTLs with MD value below the threshold are not reported", 
                metavar = "NUMERIC", default = 0.05),
    make_option(c("-t", "--type_fdr"), type = "character", 
                help = "FDR method (BH or qvalue)", metavar = "CHARACTER", 
                default = "BH"),
    make_option(c("-p", "--plot_pdf"), type = "character",
                help = "output PDF with the P-values distribution and a semi-volcano plot.", 
                metavar = "FILE", default = NULL),
    make_option(c("-o", "--output"), type = "character",
                help = "output file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

nominal.f <- opt$nominal
output.f <- opt$output

if ( is.null(nominal.f) ){
  print_help(opt_parser)
  stop("Missing/not found input file", call.= FALSE)
}

## 2. Input file: all tests  

res.df <- as.data.frame(fread(nominal.f, header = TRUE, sep = "\t"), 
                        stringsAsFactors = F)  # All tests

## 3. Get significant sQTLs and output result

sqtls.df <- sqtls(res.df, FDR = opt$fdr, FDR.svQTL = opt$fdr, 
                  method = opt$type_fdr, md.min = opt$md_min, 
                  out.pdf = opt$plot_pdf, svQTL.removal = opt$rm_svqtl) # Same FDR level for sQTLs and svQTLs

write.table(sqtls.df, file = output.f, quote = FALSE, row.names = FALSE, 
            col.names = TRUE, sep = "\t")
  
#### END