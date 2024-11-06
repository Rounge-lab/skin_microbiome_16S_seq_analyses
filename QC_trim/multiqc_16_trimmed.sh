#!/bin/bash
# Job name:
#SBATCH --job-name=multiqc-trim-16s
# Project:
#SBATCH --account=p2005_tsd
# Wall clock limit:
#SBATCH --time=1:00:00
# Memory usage:
#SBATCH --mem-per-cpu=1G
## Cpu
#SBATCH --cpus-per-task=2
##Computers
#SBATCH --nodes=1

#SBATCH --output=multiqc_q25_l150%j.out

##Useful lines to know where and when the job starts
echo "Hello " $USER 
echo "my submit directory is"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date


## Set up job environment:
module purge   # clear any inherited modules
set -o errexit # exit on errors
module load MultiQC/1.7-foss-2018b-Python-3.6.6

#Go to folder where trimmed QC files are
cd /ess/p2005/cluster/skin-metagenomics/data/trimmed_files/16S/noprimers_nextera_q25_l150/QC_trimmed/

#MultiQC on all log files in this folder
multiqc --verbose -n "MultiQC_q25_l150" .


echo "I was done at"
date

