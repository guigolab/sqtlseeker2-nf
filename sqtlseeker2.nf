/*
 * Copyright (c) 2019, Centre for Genomic Regulation (CRG)
 *
 * Copyright (c) 2019, Diego Garrido-Martín
 * 
 * This file is part of 'sqtlseeker2-nf': 
 * sQTLseekeR2 in Nextflow, a pipeline for splicing QTL mapping
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


// Define parameters

params.genotype = null
params.trexp = null
params.metadata = null
params.genes = null
params.dir = "result"
params.mode = "nominal"
params.win = 5000
params.covariates = false
params.kn = 10
params.kp = 1
params.fdr = 0.05
params.svqtl = true
params.ld = 0
params.min_md = 0.05
params.max_perm = 1000 
params.help = false


// Print usage

if (params.help) {
  log.info ''
  log.info 'sqtlseeker2-nf ~ A pipeline for splicing QTL mapping'
  log.info '----------------------------------------------------'
  log.info 'Run sQTLseekeR2 on a set of data.'
  log.info ''
  log.info 'Usage: '
  log.info "    ${workflow.projectDir.baseName} [options]"
  log.info ''
  log.info 'Options:'
  log.info '--genotype GENOTYPE_FILE    the genotype file'
  log.info '--trexp EXPRESSION_FILE     the transcript expression file'
  log.info '--metadata METADATA_FILE    the metadata file'
  log.info '--genes GENES_FILE          the gene location file' 
  log.info '--dir DIRECTORY             the output directory'
  log.info '--mode MODE                 the run mode: nominal or permuted (default: nominal)'
  log.info '--win WINDOW		the cis window in bp (default: 5000)'
  log.info '--covariates COVARIATES     include covariates in the model (default: false)'
  log.info '--fdr FDR                   false discovery rate level (default: 0.05)'
  log.info '--min_md MIN_MD             minimum effect size reported (default: 0.05)'
  log.info '--svqtl SVQTLS              test for svQTLs (default: false)'
  log.info ''
  log.info 'Additional parameters for mode = nominal:'
  log.info '--ld LD                     threshold for LD-based variant clustering (default: 0, no clustering)'
  log.info '--kn KN                     number of genes per batch in nominal pass (default: 10)'
  log.info ''
  log.info 'Additional parameters for mode = permuted:'
  log.info '--kp KP                     number of genes per batch in permuted pass (default: 10)'
  log.info '--max_perm MAX_PERM         maximum number of permutations (default: 1000)'
  log.info ''
  exit 1
}


// Check mandatory options

if (!params.genotype) {
    exit 1, "Genotype file not specified."
} else if (!params.trexp){
    exit 1, "Transcript expression file not specified."
} else if (!params.metadata){
    exit 1, "Metadata file not specified."
} else if (!params.genes){
    exit 1, "Gene location file not specified."
}

 
// Print selected options

log.info ""
log.info "sqtlseeker2-nf ~ A pipeline for splicing QTL mapping"
log.info ""
log.info "General parameters"
log.info '------------------'
log.info "Genotype file                      : ${params.genotype}"
log.info "Transcript expression file         : ${params.trexp}"
log.info "Metadata file                      : ${params.metadata}"
log.info "Gene location file                 : ${params.genes}"
log.info "Output directory                   : ${params.dir}"
log.info "Run mode                           : ${params.mode}"
log.info "Cis window                         : ${params.win}"
log.info "Covariates                         : ${params.covariates}"
log.info "FDR level                          : ${params.fdr}"
log.info "Min. effect size                   : ${params.min_md}"
log.info "Test for svQTLs                    : ${params.svqtl}"
log.info ""

if(params.mode == "nominal"){
  log.info 'Additional parameters for mode = nominal'
  log.info '----------------------------------------'
  log.info "LD-based clustering threshold      : ${params.ld}"
  log.info "Genes/batch in nominal pass        : ${params.kn}"
  log.info ""
} else if(params.mode == "permuted"){
  log.info 'Additional parameters for mode = permuted'
  log.info '-----------------------------------------'
  log.info "Genes/batch in permuted pass       : ${params.kp}"
  log.info "Max. number of permutations        : ${params.max_perm}"
  log.info ""
}


// Create file objects given parameters

genotype_file = file(params.genotype)
trexp_file = file(params.trexp)
metadata_file = file(params.metadata)
genes_file = file(params.genes)


// Obtain the 'groups' list

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


// Index genotype file

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


// Preprocess input data

process prepare {

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
    prepare_trexp.R --group "$group" -t $te -m $metadata --gene_location $genes --covariates --output_tre tre.df.RData --output_gene genes.ss.bed --output_cov covariates.df.RData

    """
    else
    """
    prepare_trexp.R --group "$group" -t $te -m $metadata --gene_location $genes --output_tre tre.df.RData --output_gene genes.ss.bed --output_cov covariates.df.RData

    """
}

tre_ch.into {tre2nominal_ch; tre2permuted_ch}
genes_ch.into {genes2nominal_ch; genes2permuted_ch} 
cov_ch.into {cov2nominal_ch; cov2permuted_ch}


// Run sQTLseekeR2 (nominal)

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
    sqtlseeker.R -t $tre_rdata -i $indexed_geno -l $chunk -c $cov_rdata --asympt --svqtl --ld ${params.ld} --window ${params.win} -o \$res

    """
}

nominal_out_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.into{all_nominal_tests_ch1; all_nominal_tests_ch2}

process nominal_mtc {

    publishDir "${params.dir}/groups/$group"
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


// Run sQTLseekeR2 (permuted)
 
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
        sqtlseeker.p.R -t $tre_rdata -i $indexed_geno -l $chunk -c $cov_rdata -M ${params.max_perm} --window ${params.win} -o \$res

        """
    }    

    permuted_out_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.set{all_permuted_tests_ch}

    all_nominal_tests_ch2.join(all_permuted_tests_ch).set{all_tests_ch}

    process permuted_mtc {

        publishDir "${params.dir}/groups/$group"
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


