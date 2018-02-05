#!/usr/bin/env Rscript

#### Prepare transcript expression file 

## 1. Load libraries and arguments

library(optparse)
library(data.table)
library(sQTLseekeR2)
  
option_list <- list(
    make_option(c("-g", "--group"), type = "character", 
                help = "select sampleIds belonging to group", metavar = "CHARACTER"),
    make_option(c("-t", "--transcript_expr"), type = "character", 
                help = "transcript expression file", metavar = "FILE"),
    make_option(c("-m", "--metadata"), type = "character",
                help = "file defining sample groups (sampleId, indId, group, [covariates])", 
                metavar = "FILE"),
    make_option(c("-c", "--covariates"), action = "store_true", 
                help = "prepare covariate file [default %default]", 
                default = FALSE),
    make_option(c("-l", "--gene_location"), type = "character",
                help = "gene location file (chr, start, end, id)", metavar = "FILE"),
    make_option(c("-e", "--min_gene_expr"), type = "numeric", default = 1,
                help = "minimum gene expression [default %default]", metavar = "NUMERIC"),
    make_option(c("-p", "--min_proportion"), type = "numeric", default = 0.8,
                help = "minimum proportion of samples with gene expression above --min_gene_expr. [default %default]", 
                metavar = "NUMERIC"),
    make_option(c("-i", "--min_transcript_expr"), type = "numeric", default = 0.1,
                help = "minimum transcript expression. [default %default]", metavar = "NUMERIC"),
    make_option(c("-d", "--min_dispersion"), type = "numeric", default = 0.1,
                help = "minimum dispersion of transcript relative expression. [default %default]", 
                metavar = "NUMERIC"),
    make_option(c("-o1", "--output_tre"), type = "character",
                help = "preprocessed transcript expression file", metavar = "FILE"),
    make_option(c("-o2", "--output_gene"), type = "character",
                help = "updated gene location file", metavar = "FILE"),
    make_option(c("-o3", "--output_cov"), type = "character",
                help = "prepared covariate file", metavar = "FILE"),
    make_option(c("-s", "--seed"), type = "numeric", help = "Set seed for random processess",
                metavar = "NUMERIC", default = 123),
    make_option(c("-v", "--verbose"), action = "store_true", 
                help = "print genes and transcripts filtered out [default %default]", 
                default = TRUE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## 2. Input files: transcript expression, sample groups, gene location
  
trans.expr.f <- opt$transcript_expr                        
metadata.f <- opt$metadata
genes.bed.f <- opt$gene_location
sel.group <- opt$group
out.tre.f <- opt$output_tre
out.gene.f <- opt$output_gene
out.cov.f <- opt$output_cov

if ( is.null(trans.expr.f) || is.null(metadata.f) || is.null(genes.bed.f) || 
     is.null(out.tre.f) || is.null (out.gene.f) ){
    print_help(opt_parser)
    stop("Missing/not found input files", call.= FALSE)
}
  
genes.bed <- read.table(genes.bed.f, header = TRUE, as.is = TRUE, sep = "\t")
  
## 3. Getting the IDs of the samples in the group of interest
  
metadata <- read.table(metadata.f, header = TRUE,
                       as.is = TRUE, sep = "\t")

subset.df <- subset(metadata, group == sel.group)                                   # Select samples of interest            
subset.samples <- subset.df$sampleId

## 4. Prepare transcript expression

if( grepl("\\.gz$", trans.expr.f) ){
    trans.expr.f <- paste0("zcat < '", trans.expr.f, "'") 
}

te.df <- as.data.frame(fread(input = trans.expr.f, header = TRUE, sep = "\t"))
colnames(te.df)[1:2] <- c("trId", "geneId")                                         # Proper 1,2 colnames
subset.samples <- subset.samples[subset.samples%in%colnames(te.df)]                 # Get samples that have quantifications 
te.df <- te.df[, c("trId", "geneId", subset.samples)]                               # Select subgroup of samples = the ones from the group of interest
te.df <- subset(te.df, geneId %in% genes.bed$geneId)                                # Remove from te.df all genes that are not in genes.bed        

set.seed(opt$seed)
tre.df <- prepare.trans.exp(te.df, min.gene.exp = opt$min_gene_expr, 
                            min.transcript.exp = opt$min_transcript_expr, 
                            min.dispersion = opt$min_dispersion,
                            min.prop = opt$min_proportion,
                            verbose = opt$verbose)                                  # Run

## 5. Prepare covariate file

if(opt$covariates) {
    if ( is.null(out.cov.f) ){
        print_help(opt_parser)
        stop("Missing output covariate file", call.= FALSE)
    }
    covariates.df <- subset.df[, setdiff(colnames(subset.df), 
                                         c("sampleId", "indId", "group"))]
    rownames(covariates.df) <- subset.df$indId
    for(i in 1:ncol(covariates.df)){
        typ <- class(covariates.df[, i])
        if(typ == "character"){
            covariates.df[, i] <- as.factor(covariates.df[, i])
        } else if (typ == "numeric" || typ == "integer"){
            next
        } else {
            stop ("Covariates should be either 'numeric' or 'character'")
        }
    }
    types <- unlist(lapply(covariates.df, class))
    if (opt$verbose) {
        message("Covariate types:\n", 
                paste(names(types), types, sep=": ", collapse = ", "))
    }
    save(covariates.df, file = out.cov.f) 
} else {
    covariates.df <- NULL
    save(covariates.df, file = out.cov.f) 
}                                                                                   

colnames(tre.df)[-c(1:2)] <- subset.df$indId                                        # Rename colnames to individual ID
                                                        
## 6. Save result

save(tre.df, file = out.tre.f)                                                      # Save tre.df as RData

## 7. Preprocess for split

genes.bed <- subset(genes.bed, geneId %in% tre.df$geneId)                           # Remove from gene.bed all genes that are not in tre.df
write.table(genes.bed, file = out.gene.f, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t")                       # Write gene list

#### END
