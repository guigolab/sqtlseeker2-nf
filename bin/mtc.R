#!/usr/bin/env Rscript

#### Perform multiple testing correction (eigenMT + FDR) 

## 1. Load libraries and arguments

library(optparse)
library(data.table)

option_list <- list(
    make_option(c("-n", "--nominal"), type = "character",
                help = "all nominal tests", metavar = "FILE"),
    make_option(c("-e", "--eigenMT"), type = "character",
		help = "eigenMT output with colnames snpId, geneId, pv, BF, Meff", 
                metavar = "FILE"),
    make_option(c("-f", "--fdr"), type = "numeric", 
                help = "FDR level", 
                metavar = "NUMERIC", default = 0.05),
    make_option(c("-b", "--bonferroni"), type = "numeric",
		help = "Bonferroni level",
	        metavar = "NUMERIC", default = 0.05),	
    make_option(c("-o", "--output"), type = "character",
                help = "output file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

nominal.f <- opt$nominal
eigenMT.f <- opt$eigenMT
output.f <- opt$output

if ( is.null(nominal.f) || is.null(eigenMT.f) || is.null(output.f) ){
    print_help(opt_parser)
    stop("Required I/O files not found", call.= FALSE)
}

## 2. Read input files  

nominal.df <- fread(nominal.f, header = TRUE, sep = "\t", data.table = F)
eigenMT.df <- fread(eigenMT.f, header = TRUE, sep = "\t", data.table = F)

## 3. Get significant sQTLs and output result

eigenMT.df$fdr <- p.adjust(eigenMT.df$BF, "fdr")
eigenMT.df <- subset(eigenMT.df, subset = fdr < opt$fdr, select = c("geneId", "Meff"))

merged.df <- merge(eigenMT.df, nominal.df, by = "geneId")

corrected.df <- subset(merged.df, subset = pv_snpXcond*Meff < opt$bonferroni) 

write.table(corrected.df, file = output.f, quote = FALSE, row.names = FALSE,
            col.names = TRUE, sep = "\t")

#### END
