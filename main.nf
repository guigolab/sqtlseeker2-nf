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
params.sgroups = "data/sample.groups.tsv"
params.genes = "data/genes.bed"
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
         sgroups: ${params.sgroups}
         genes: ${params.genes}
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
sgroups_file = file(params.sgroups)
genes_file = file(params.genes)

/*
 *  0.c Generate the 'groups' list
 */

def groups = []
myReader = sgroups_file.newReader()
String line
while( line = myReader.readLine() ) {
    def (sampleID, group) = line.tokenize('\t')
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
    file sgroups from sgroups_file
    file genes from genes_file

    output:
    set val(group), file('tre.df.RData') into tre_ch
    set val(group), file('genes.ss.bed') into genes_ch

    script:
    """
    prepare_trexp.R -G $group -t $te_rdata -s $sgroups --gene_location $genes --output_tre tre.df.RData --output_gene genes.ss.bed

    """
}

tre_ch.into {tre2nominal_ch; tre2permuted_ch}
genes_ch.into {genes2nominal_ch; genes2permuted_ch} 

/*
 *  3. Run sQTLseekeR (nominal)
 */

tre2nominal_ch.combine(genes2nominal_ch.splitText( by: params.kn, file: "nominal_in" ), by: 0).set{nominal_in_ch}

process nominal_test {

    tag {"$group, $chunk"}
 
    input:
    set val(group), file(tre_rdata), file (chunk) from nominal_in_ch
    set file(indexed_geno), file(tbi) from index2nominal_ch 

    output:
    set val(group), file('nominal_out.*') into nominal_out_ch 
 
    script: 
    """
    res=\$(echo $chunk | sed 's/_in/_out/')
    sqtlseeker.R -t $tre_rdata -i $indexed_geno -g $chunk -o \$res

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
    sqtls.R -a all-tests.nominal.tsv -f ${params.fdr} -o sqtls-${params.fdr}fdr.nominal.tsv

    """       
}

/*
 *  4. Run sQTLseekeR (permuted)
 */
 
if (params.mode == "permuted") {

    tre2permuted_ch.combine(genes2permuted_ch.splitText( by: params.kp, file: "permuted_in" ), by: 0).set{permuted_in_ch}

    process permuted_test {

        tag {"$group, $chunk"}

        input:
	set val(group), file(tre_rdata), file (chunk) from permuted_in_ch
        set file(indexed_geno), file(tbi) from index2permuted_ch

        output:
        set val(group), file('permuted_out.*') into permuted_out_ch

        script:
        """
        res=\$(echo $chunk | sed 's/_in/_out/')
        sqtlseeker.p.R -t $tre_rdata -i $indexed_geno -g $chunk -o \$res

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
        sqtls.p.R -n all-tests.nominal.tsv -p all-tests.permuted.tsv -f ${params.fdr} -o sqtls-${params.fdr}fdr.permuted.tsv

        """
    }

}


