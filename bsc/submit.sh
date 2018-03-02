#!/bin/bash
#SBATCH --job-name=sqtlseeker-nf
#SBATCH --output=sqtlseeker-nf.log
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=48
#SBATCH --tasks-per-node=1
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
export NXF_TEMP=$TMPDIR
module load java
module load intel/2018.1
module load singularity
srun nextflow run sQTLseekeR.nf -w ~/scratch/work -with-report report.html -with-trace -with-timeline time.html -with-dag flowchart.dot --mode permuted --genotype gtex/snps.tsv --trexp gtex/tq.rsem.tsv.gz --genes gtex/genes.bed --svqtl true --metadata gtex/metadata.tsv --kn 10 -with-mpi
