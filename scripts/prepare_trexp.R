#!/usr/bin/env Rscript

#### Prepare transcript expression file 

## 1. Load libraries and arguments

library(optparse)
library(sQTLseekeR2)
  
option_list <- list(
    make_option(c("-G", "--Group"), type = "character", 
                help = "select sampleIds belonging to Group", metavar = "SELECTION"),
    make_option(c("-t", "--transcript_expr"), type = "character", 
                help = "transcript expression file", metavar = "FILE"),
    make_option(c("-s", "--sample_groups"), type = "character",
                help = "file defining groups of samples (sampleId, group)", metavar = "FILE"),
    make_option(c("-g", "--gene_location"), type = "character",
                help = "gene location file (chr, start, end, id)", metavar = "FILE"),
    make_option(c("-m", "--min_gene_expr"), type = "numeric", default = 1,
                help = "minimum gene expression [default %default]", metavar = "NUMERIC"),
    make_option(c("-i", "--min_transcript_expr"), type = "numeric", default = 0.1,
                help = "minimum transcript expression. [default %default]", metavar = "NUMERIC"),
    make_option(c("-o1", "--output_tre"), type = "character",
                help = "preprocessed transcript expression file", metavar = "FILE"),
    make_option(c("-o2", "--output_gene"), type = "character",
                help = "preprocessed transcript expression file", metavar = "FILE"),
    make_option(c("-S", "--Seed"), type = "numeric", help = "Set seed for random processess",
                metavar = "NUMERIC", default = 123),
    make_option(c("-v", "--verbose"), action = "store_true", 
                help = "print genes and transcripts filtered out [default %default]", 
                default = FALSE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## 2. Input files: transcript expression, sample groups, gene location
  
trans.expr.f <- opt$transcript_expr                        
sample.groups.f <- opt$sample_groups                     
genes.bed.f <- opt$gene_location
sel.group <- opt$Group
out.tre.f <- opt$output_tre
out.gene.f <- opt$output_gene

if ( is.null(trans.expr.f) || is.null(sample.groups.f) || is.null(genes.bed.f) ){
    print_help(opt_parser)
    stop("Missing/not found input files", call.= FALSE)
}
  
genes.bed <- read.table(genes.bed.f, header = TRUE, as.is = TRUE, sep = "\t")
  
## 3. Getting the IDs of the samples in the group of interest
  
sample.groups <- read.table(sample.groups.f, header = TRUE, 
                            as.is = TRUE, sep = "\t")
subset.samples <- subset(sample.groups, group == sel.group)$sampleId                # Select samples of interest            

## 4. Prepare transcript expression

load(trans.expr.f) 
colnames(te.df)[1:2] <- c("trId", "geneId")                                         # Proper 1,2 colnames
subset.samples <- subset.samples[subset.samples%in%colnames(te.df)]                 # Get samples that have quantifications 
te.df <- te.df[, c("trId", "geneId", subset.samples)]                               # Select subgroup of samples = the ones from the group of interest
te.df <- subset(te.df, geneId %in% genes.bed$geneId)                                # Remove from te.df all genes that are not in genes.bed        

set.seed(opt$Seed)
tre.df <- prepare.trans.exp(te.df, min.gene.exp = opt$min_gene_expr, 
                                   min.transcript.exp = opt$min_transcript_expr, 
                                   verbose = opt$verbose)                           # Run

## 5. Save result

save(tre.df, file = out.tre.f)                                                      # Save tre.df as RData

## 6. Preprocess for split

genes.bed <- subset(genes.bed, geneId %in% tre.df$geneId)                           # Remove from gene.bed all genes that are not in tre.df
write.table(genes.bed, file = out.gene.f, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t")                       # Write gene list

#### END
