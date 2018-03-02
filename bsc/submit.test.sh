#!/bin/bash
#SBATCH --job-name=sqtlseeker-nf
#SBATCH --output=sqtlseeker-nf.test.log.txt
#SBATCH --ntasks=2
#SBATCH --qos=debug
#SBATCH --cpus-per-task=4
#SBATCH --tasks-per-node=1
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
export NXF_TEMP=$TMPDIR
module load java
module load intel/2018.1
module load singularity
srun nextflow run sQTLseekeR.nf --kn 1 -w ~/scratch/work -with-trace -with-report report.html -with-timeline time.html -with-dag flowchart.html -with-mpi
