#!/usr/bin/env Rscript

#### Test gene/SNP association (permutation pass)

## 1. Load libraries and arguments

library(optparse)
library(sQTLseekeR2)

option_list <- list(
    make_option(c("-t", "--transcript_expr"), type = "character",
                help = "Prepared transcript expression RData file", metavar = "FILE"),
    make_option(c("-i", "--indexed_geno"), type = "character",
                help = "Indexed genotype file", metavar = "FILE"),
    make_option(c("-l", "--gene_location"), type = "character",
                help = "gene location chunk file", metavar = "FILE"),
    make_option(c("-c", "--covariates"), type = "character", 
                help = "file of covariates to regress out before testing", 
                metavar = "FILE"),
    make_option(c("-o", "--output_file"), type = "character", help = "output file", 
                metavar = "FILE"),
    make_option(c("-w", "--window"), type = "numeric", 
                help = "genic window in bp [default %default]", 
                metavar = "NUMERIC", default = 5000),
    make_option(c("-e", "--min_nb_ext_scores"), type = "numeric", 
                help = "minimum number of external scores [default %default]", 
                metavar = "NUMERIC", default = 100),
    make_option(c("-M", "--nb_perm_max"), type = "numeric", 
                help = "maximum number of permutations [default %default]", 
                metavar = "NUMERIC", default = 1000),
    make_option(c("-m", "--nb_perm_min"), type = "numeric", 
                help = "minimum number of permutations [default %default]", 
                metavar = "NUMERIC", default = 100),
    make_option(c("-n", "--min_nb_ind_geno"), type = "numeric", 
                help = "minimum number of individuals per genotype group [default %default]", 
                metavar = "NUMERIC", default = 10),
    make_option(c("-s", "--seed"), type = "numeric", help = "set seed for random processess", 
                metavar = "NUMERIC", default = 123),
    make_option(c("-v", "--verbose"), action = "store_true", 
                help = "print progress (geneId, SNP suitability, etc.) [default %default]",
                default = TRUE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## 2. Input files: prepared transcript expression, genotypes (indexed) and gene location (chunk) 

trans.expr.p.f <- opt$transcript_expr
indexed.geno.f <- opt$indexed_geno
gene.loc.chunk <- opt$gene_location 
covariates.f <- opt$covariates  
output.f <- opt$output_file


if ( is.null(trans.expr.p.f) || is.null (indexed.geno.f) || 
     is.null(gene.loc.chunk) || is.null (output.f) ){
    print_help(opt_parser)
    stop("Missing/not found input files", call.= FALSE)
}

load(trans.expr.p.f)                                                            # Load tre.df
genes.bed <- read.table(gene.loc.chunk, header = FALSE, as.is = TRUE)           # Load chunk
colnames(genes.bed) <- c("chr", "start", "end", "geneId")                       # Name chunk
genes <- genes.bed$geneId                                                       # Get gene names
tre.df <- subset(tre.df, geneId %in% genes)                                     # Subset tre.df
load(covariates.f)                                                              # Load covariates.df

## 3. Run association test (permuted) & write result

set.seed(opt$seed)

res.df <- sqtl.seeker.p(tre.df, indexed.geno.f, genes.bed, 
                        covariates = covariates.df,
                        genic.window = opt$window,
                        min.nb.ext.scores = opt$min_nb_ext_scores,
                        nb.perm.max = opt$nb_perm_max,
                        nb.perm.min = opt$nb_perm_min,
                        min.nb.ind.geno = opt$min_nb_ind_geno,
                        verbose = opt$verbose)

if (output.f == "permuted_out.1"){
    if(is.null(res.df)){
        colnms <- c("geneId", "variants.cis", "LD", "best.snp", "best.nominal.pv", 
                    "shape1", "shape2", "nb.perms", "pv.emp", "pv.emp.beta", "runtime")
        write.table(t(colnms), file = output.f, col.names = FALSE, row.names = FALSE, 
                    quote = FALSE, sep = "\t")
    }else{
        write.table(res.df, file = output.f, quote = FALSE, row.names = FALSE, 
                    col.names = TRUE, sep = "\t")
    }
} else {
     write.table(res.df, file = output.f, quote = FALSE, row.names = FALSE, 
                 col.names = FALSE, sep = "\t")
}


#### END
