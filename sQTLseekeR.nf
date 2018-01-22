/*
 * sQTL mapping pipeline
 * vs 1.0
 * Diego Garrido Martin 
 */

/*
 * 0.a Define parameters
 */

params.genotype = "data/snps-012coded.tsv.gz"
params.trexp = "data/transExpression.tsv.gz"
params.metadata = "data/metadata.tsv"
params.genes = "data/genes.bed"
params.mode = "nominal"
params.covariates = false
params.kn = 2
params.kp = 1
params.fdr = 0.05
params.svqtl = false
params.ld = 0 
params.min_md = 0.05
params.max_perm = 1000 

println """\
         =============================
         s Q T L s e e k e R - N F   
         =============================
         genotype: ${params.genotype}
         trexp: ${params.trexp}
         metadata: ${params.metadata}
         genes: ${params.genes}
         mode: ${params.mode}         
         covariates: ${params.covariates}
         kn: ${params.kn}
         kp: ${params.kp}
         fdr: ${params.fdr}
         svqtl: ${params.svqtl}
         ld: ${params.ld}
         min_md: ${params.min_md}
         max_perm: ${params.max_perm}"""
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
    set file("${genotype.baseName}.*bgz"), file("${genotype.baseName}.*bgz.tbi") into index_ch
    
    script:
    """
    index_geno.R --genotype_file $genotype 

    """
}

index_ch.into{index2nominal_ch; index2permuted_ch}

/*
 *  2. Preprocess input data
 */

process prepare {

    publishDir "result/groups/$group" 
    tag { group }

    input:
    val group from Channel.from(groups)
    file te from trexp_file 
    file metadata from metadata_file
    file genes from genes_file

    output:
    set val(group), file('tre.df.RData') into tre_ch
    set val(group), file('genes.ss.bed') into genes_ch
    set val(group), file('covariates.df.RData') into cov_ch
    
    script:
    if (params.covariates == true)
    """
    prepare_trexp.R --group $group -t $te -m $metadata --gene_location $genes --covariates --output_tre tre.df.RData --output_gene genes.ss.bed --output_cov covariates.df.RData

    """
    else
    """
    prepare_trexp.R --group $group -t $te -m $metadata --gene_location $genes --output_tre tre.df.RData --output_gene genes.ss.bed --output_cov covariates.df.RData

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
    if (params.svqtl == false)
    """
    res=\$(echo $chunk | sed 's/_in/_out/')
    sqtlseeker.R -t $tre_rdata -i $indexed_geno -l $chunk -c $cov_rdata --asympt --ld ${params.ld} -o \$res 

    """
    else
    """
    res=\$(echo $chunk | sed 's/_in/_out/')
    sqtlseeker.R -t $tre_rdata -i $indexed_geno -l $chunk -c $cov_rdata --asympt --svqtl --ld ${params.ld} -o \$res

    """
}

nominal_out_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.into{all_nominal_tests_ch1; all_nominal_tests_ch2}
  
process nominal_mtc {

    publishDir "result/groups/$group"
    tag { group }

    input: 
    set val(group), file('all-tests.nominal.tsv') from all_nominal_tests_ch1

    output:
    set val(group), file('all-tests.nominal.tsv'), file ("sqtls-${params.fdr}fdr.nominal.tsv") into nominal_end_ch
 
    script:
    """
    sqtls.R -n all-tests.nominal.tsv -f ${params.fdr} --rm_svqtl --md_min ${params.min_md} -o sqtls-${params.fdr}fdr.nominal.tsv

    """       
}

/*
 *  4. Run sQTLseekeR (permuted)
 */
 
if (params.mode == "permuted") {

    tre2permuted_ch.join(cov2permuted_ch).combine(genes2permuted_ch.splitText( by: params.kp, file: "permuted_in" ), by: 0).set{permuted_in_ch}

    process permuted_test {

        tag {"$group, $chunk"}

        input:
	set val(group), file(tre_rdata), file(cov_rdata), file (chunk) from permuted_in_ch
        set file(indexed_geno), file(tbi) from index2permuted_ch

        output:
        set val(group), file('permuted_out.*') into permuted_out_ch

        script:
        """
        res=\$(echo $chunk | sed 's/_in/_out/')
        sqtlseeker.p.R -t $tre_rdata -i $indexed_geno -l $chunk -c $cov_rdata -M ${params.max_perm} -o \$res

        """
    }    

    permuted_out_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.set{all_permuted_tests_ch}

    all_nominal_tests_ch2.join(all_permuted_tests_ch).set{all_tests_ch}

    process permuted_mtc {

        publishDir "result/groups/$group"
        tag { group }

        input:
        
        set val(group), file('all-tests.nominal.tsv'), file('all-tests.permuted.tsv') from all_tests_ch

        output:
        set val(group), file('all-tests.permuted.tsv'), file ("sqtls-${params.fdr}fdr.permuted.tsv") into permuted_end_ch
 
        script:
        """
        sqtls.p.R -n all-tests.nominal.tsv -p all-tests.permuted.tsv -f ${params.fdr} --rm_svqtl --md_min ${params.min_md} -o sqtls-${params.fdr}fdr.permuted.tsv

        """
    }

}


