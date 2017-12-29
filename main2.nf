/*
 * sQTL mapping pipeline
 * vs 0.1
 * Diego Garrido Martin 
 */

/*
 * 0.a Define parameters
 */

params.genotype = "data/snps-012coded.tsv.gz"
params.trexp = "data/transExpression.tsv.gz"
params.metadata = "data/metadata.tsv"
params.genes = "data/genes.bed"
params.covariates = false
params.kn = 2
params.kp = 1
params.fdr = 0.05
params.mode = "nominal"

println """\
         =============================
         s Q T L s e e k e R - N F   
         =============================
         genotype: ${params.genotype}
         trexp: ${params.trexp}
         metadata: ${params.metadata}
         genes: ${params.genes}
         covariates: ${params.covariates}
         kn: ${params.kn}
         kp: ${params.kp}
         fdr: ${params.fdr}
         mode: ${params.mode}"""
         .stripIndent()

/*
 *  0.b Create file objects given parameters
 */

genotype_file = file(params.genotype)
trexp_file = file(params.trexp)
metadata_file = file(params.metadata)
genes_file = file(params.genes)

/*
 *  0.c Generate the 'groups' list
 */

def groups = []
myReader = metadata_file.newReader()
String line
while( line = myReader.readLine() ) {
    def (sampleId, indId, group, covariates) = line.tokenize('\t')
    if( group != 'group' ) {
    	groups += group
    }
}
myReader.close()
groups.unique()
String show = groups.join(", ")

println "groups: $show\n"

/*
 *  1. Index genotype file
 */

process index {

    input:
    file genotype from genotype_file

    output:    
    set file("${genotype.baseName}.bgz"), file("${genotype.baseName}.bgz.tbi") into index_ch
    
    script:
    """
    index_geno.R --genotype_file $genotype 

    """
}

index_ch.into{index2nominal_ch; index2permuted_ch}

/*
 *  2. Preprocess input data
 *     a. make RData from transcript quantifications for faster access
 *     b. preprocess transcript expression matrix
 */

process make_RData {

    input:
    file trexp from trexp_file

    output:
    file 'te.df.RData' into te_ch

    script:
    """
    make_RData.R -i $trexp -o te.df.RData -H

    """
}

process prepare {

    publishDir "result/groups/$group" 
    tag { group }

    input:
    val group from Channel.from(groups)
    file te_rdata from te_ch 
    file metadata from metadata_file
    file genes from genes_file

    output:
    set val(group), file('tre.df.RData') into tre_ch
    set val(group), file('genes.ss.bed') into genes_ch
    set val(group), file('covariates.df.RData') into cov_ch
    
    script:
    if (params.covariates == true)
    """
    prepare_trexp.R -g $group -t $te_rdata -m $metadata --gene_location $genes --covariates --output_tre tre.df.RData --output_gene genes.ss.bed --output_cov covariates.df.RData

    """
    else
    """
    prepare_trexp.R -g $group -t $te_rdata -m $metadata --gene_location $genes --output_tre tre.df.RData --output_gene genes.ss.bed --output_cov covariates.df.RData

    """
}

tre_ch.into {tre2nominal_ch; tre2permuted_ch}
genes_ch.into {genes2nominal_ch; genes2permuted_ch} 
cov_ch.into {cov2nominal_ch; cov2permuted_ch}

/*
 *  3. Run sQTLseekeR (nominal)
 */

tre2nominal_ch.join(cov2nominal_ch).combine(genes2nominal_ch.splitText( by: params.kn, file: "nominal_in" ), by: 0).set{nominal_in_ch}

process nominal_test {

    tag {"$group, $chunk"}
 
    input:
    set val(group), file(tre_rdata), file(cov_rdata), file (chunk) from nominal_in_ch
    set file(indexed_geno), file(tbi) from index2nominal_ch 

    output:
    set val(group), file('nominal_out.*') into nominal_out_ch 
 
    script: 
    """
    res=\$(echo $chunk | sed 's/_in/_out/')
    sqtlseeker.R -t $tre_rdata -i $indexed_geno -l $chunk -c $cov_rdata --asympt -o \$res

    """
}

