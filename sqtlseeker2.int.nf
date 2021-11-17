/*
 * Copyright (c) 2019, Centre for Genomic Regulation (CRG)
 *
 * Copyright (c) 2019, Diego Garrido-Mart<C3><AD>n
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
params.covariates = true
params.condition  = "-"
params.win = 5000
params.kn = 10
params.svqtl = false
params.ld = 0 
params.min_nb_ind = 10
params.min_gene_expr = 1
params.min_transcript_expr = 0.1
params.fdr = 0.1
params.bonferroni = 0.05
params.dir = "result"
params.help = false


// Print usage

if (params.help) {
  log.info ''
  log.info 'sqtlseeker2.int-nf ~ A pipeline for condition-biased splicing QTL mapping'
  log.info '----------------------------------------------------'
  log.info 'Run sQTLseekeR2.int on a set of data.'
  log.info ''
  log.info 'Usage: '
  log.info "    ${workflow.projectDir.baseName} [options]"
  log.info ''
  log.info 'Options:'
  log.info '--genotype GENOTYPE_FILE    the genotype file (variant ID format should be "^\$chr_.*")'
  log.info '--trexp EXPRESSION_FILE     the transcript expression file'
  log.info '--metadata METADATA_FILE    the metadata file'
  log.info '--genes GENES_FILE          the gene location file' 
  log.info '--covariates COVARIATES     include covariates in the model (default: false)'
  log.info '--condition CONDITION       a 2-level factor to include in the model (default: -)'
  log.info '--win WINDOW                the cis window in bp (default: 5000)'
  log.info '--kn KN                     number of genes per batch in nominal pass (default: 10)'
  log.info '--svqtl SVQTLS              test for svQTLs (default: false)'
  log.info '--ld LD                     threshold for LD-based variant clustering (default: 0, no clustering)'
  log.info '--min_nb_ind MNI            minimum number of individuals per genotype (default: 0.05)'
  log.info '--min_gene_expr MGE         minimum gene expression (default: 1)'
  log.info '--min_transcript_expr MTE   minimum transcript expression (default: 0.1)'
  log.info '--fdr FDR                   false discovery rate level (default: 0.05)'
  log.info '--bonferroni BONFERRONI     Bonferroni correction level (default: 0.05)' 
  log.info '--dir DIRECTORY             the output directory'
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
log.info "sqtlseeker2.int-nf ~ A pipeline for condition-biased splicing QTL mapping"
log.info ""
log.info "General parameters"
log.info '------------------'
log.info "Genotype file                      : ${params.genotype}"
log.info "Transcript expression file         : ${params.trexp}"
log.info "Metadata file                      : ${params.metadata}"
log.info "Gene location file                 : ${params.genes}"
log.info "Covariates                         : ${params.covariates}"
log.info "Condition (2-level factor)         : ${params.condition}"
log.info "Cis window                         : ${params.win}"
log.info "Genes/batch in nominal pass        : ${params.kn}"
log.info "Test for svQTLs                    : ${params.svqtl}"
log.info "LD-based clustering threshold      : ${params.ld}" 
log.info "Min. #individuals/GT               : ${params.min_nb_ind}"
log.info "Min. gene expression               : ${params.min_gene_expr}"
log.info "Min. transcript expression         : ${params.min_transcript_expr}"
log.info "FDR level                          : ${params.fdr}"
log.info "Bonferroni level                   : ${params.bonferroni}"
log.info "Output directory                   : ${params.dir}"
log.info ""


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


// Preprocess input data

process prepare {

    publishDir "${params.dir}/groups/$group" 
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
    prepare_trexp.R --group "$group" -t $te -m $metadata --gene_location $genes --covariates --output_tre tre.df.RData --output_gene genes.ss.bed --output_cov covariates.df.RData --min_gene_expr ${params.min_gene_expr} --min_transcript_expr ${params.min_transcript_expr}

    """
    else
    """
    prepare_trexp.R --group "$group" -t $te -m $metadata --gene_location $genes --output_tre tre.df.RData --output_gene genes.ss.bed --output_cov covariates.df.RData --min_gene_expr ${params.min_gene_expr} --min_transcript_expr ${params.min_transcript_expr}

    """
}


// Run sQTLseekeR.int

tre_ch.join(cov_ch).combine(genes_ch.splitText( by: params.kn, file: "nominal_in" ), by: 0).set{nominal_in_ch}

process nominal_test {

    tag {"$group, $chunk"}
 
    input:

    set val(group), file(tre_rdata), file(cov_rdata), file (chunk) from nominal_in_ch
    set file(indexed_geno), file(tbi) from index_ch 

    output:
    set val(group), file('nominal_out.*') into nominal_out_ch 
 
    script: 
    if (params.svqtl == false)
    """
    res=\$(echo $chunk | sed 's/_in/_out/')
    sqtlseeker.R -t $tre_rdata -i $indexed_geno -l $chunk -c $cov_rdata --condition ${params.condition} --win ${params.win} --min_nb_ind_geno ${params.min_nb_ind} --asympt --ld ${params.ld} -o \$res 
    """
    else
    """
    res=\$(echo $chunk | sed 's/_in/_out/')
    sqtlseeker.R -t $tre_rdata -i $indexed_geno -l $chunk -c $cov_rdata --condition ${params.condition} --win ${params.win} --min_nb_ind_geno ${params.min_nb_ind} --asympt --svqtl --ld ${params.ld} -o \$res
    """
}

nominal_out_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.into{uncorrected1_ch; uncorrected2_ch}


//  Prepare eigenMT input files 

process eigenMT_prepare {

    input:
    file genotype from genotype_file
    file genes from genes_file    

    output:
    file ('eigenMT.*') into eigenMT_prep_ch

    script:
    """
    # Parse genotypes
    if [[ \$(echo $genotype | grep '.gz\$') ]]; then
        zcat $genotype | cut -f1-3 --complement | sed 1's/snpId/ID/' > eigenMT.GEN    
        zcat $genotype | awk 'BEGIN{OFS=FS="\t"}{print \$4,\$1,\$2}' > eigenMT.GENPOS
    else
        cat $genotype | cut -f1-3 --complement | sed 1's/snpId/ID/' > eigenMT.GEN
        cat $genotype | awk 'BEGIN{OFS=FS="\t"}{print \$4,\$1,\$2}' > eigenMT.GENPOS
    fi

    # Parse genes
    awk 'BEGIN{FS=OFS="\t"}{print \$4,\$1,\$2,\$2}' $genes > eigenMT.PHEPOS    
    """
}


// Run eigenMT multiple testing correction

process eigenMT {
    
    tag { "$group, chr$chr" }


    input:
    set val(group), file(uncorrected) from uncorrected1_ch 
    set file(gen), file(genpos), file(phepos) from eigenMT_prep_ch
    each chr from Channel.of(1..22, 'X') 

    output:
    set val(group), file("eigenMT.$chr") into eigenMT_ch  

    script:
    """
    if [[ \$(wc -l $uncorrected | cut -d' ' -f1) > 1 ]]; then    
        awk 'NR>1{print \$2"\t"\$1"\t"\$8}' $uncorrected | sed "1 s,^,SNP\tgene\tp-value\\n," > qtl
        if [[ \$(grep "^${chr}_" qtl | wc -l) > 0 ]]; then  
            eigenMT.py --CHROM $chr --GEN <(cat <(head -1 $gen) <(grep "^${chr}_" $gen)) --GENPOS $genpos --PHEPOS $phepos --QTL qtl --OUT eigenMT.$chr
            sed -i '1d' eigenMT.$chr        
        else
            touch eigenMT.$chr
        fi
    else
        touch eigenMT.$chr
    fi
    """
}

eigenMT_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.set{mtc_ch}

process mtc {

    publishDir "${params.dir}/groups/$group"
    tag { group }

    input:
    set val(group), file('all-tests.nominal.tsv') from uncorrected2_ch
    set val(group), file('eigenMT.tsv') from mtc_ch

    output:
    set val(group), file('all-tests.nominal.tsv'), file("sqtls-${params.fdr}fdr.eigenMT.tsv"), file('eigenMT.tsv') into end_ch

    script:
    """
    sed -i '1 s,^,snpId\tgeneId\tpv\tBF\tMeff\\n,' eigenMT.tsv
    if [[ \$(wc -l eigenMT.tsv | cut -d' ' -f1) > 1 ]]; then
        mtc.R -n all-tests.nominal.tsv -e eigenMT.tsv -f ${params.fdr} -b ${params.bonferroni} -o sqtls-${params.fdr}fdr.eigenMT.tsv
    else
        echo "geneId\tMeff\t\$(head -1 all-tests.nominal.tsv | cut -f2-)" > sqtls-${params.fdr}fdr.eigenMT.tsv
    fi 
    """
}
