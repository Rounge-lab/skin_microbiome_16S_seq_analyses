#!/bin/bash
#SBATCH --job-name=quality_check
#SBATCH --account=p2005_tsd
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1
#SBATCH --output=fastqc-check%j.out

echo "Hello " $USER
echo "my submit directory is"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID

echo "Today is:"
date

module purge
set -o errexit
module load fastqc/0.11.8

echo "Running fastqc.."
fastqc -f fastq $1 $2 -o "/ess/p2005/cluster/skin-metagenomics/data/raw/seqdata/16S/230323_M07166.Project_Glenna-16S1-2023-03-15/QC/"

echo "I've done at "
date



