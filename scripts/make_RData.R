#!/usr/bin/env Rscript

#### Make RData given an input file

## 1. Load libraries and arguments

library(optparse)
library(data.table)

option_list = list(
  make_option(c("-i", "--input"), type = "character",
              help = "input file", metavar = "INPUT"),
  make_option(c("-o", "--output"), type = "character",
              help = "output RData file", metavar = "RDATA"),
  make_option(c("-H", "--Header"), action = "store_true",
              help = "header is present", default = FALSE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input.f <- opt$input                   
output.f <- opt$output                    

if ( is.null(input.f) || is.null(output.f) ){
  print_help(opt_parser)
  stop("Missing/not found input files", call.= FALSE)
}

## 2. Read input file and write RData

if( grepl("\\.gz$", input.f) ){ input.f <- paste0("zcat < '", input.f, "'") }

te.df <- as.data.frame(fread(input = input.f, header = opt$Header, sep = "\t"))

save(te.df, file = output.f)

####  END
