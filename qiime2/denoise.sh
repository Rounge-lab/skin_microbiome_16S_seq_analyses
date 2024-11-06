#!/bin/bash
#SBATCH --job-name=denoising_dada2
#SBATCH --account=p2005_tsd
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --output=denoise_above_5000_%j.out

echo "Hello " $USER
echo "my submit directory is"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID

echo "Today is:"
date

echo "Denoising the trimmed 16S sequences with dada2..."

qiime dada2 denoise-paired --p-n-threads 4 --i-demultiplexed-seqs /ess/p2005/cluster/skin-metagenomics/qiime2/all_samples/all_above_5000_reads/16s_trimmed_above_5000_reads.qza --p-trunc-len-f 270 --p-trunc-len-r 210 --o-representative-sequences rep-seqs-all-above5000_f270_r210.qza --o-table asv_all_above_5000_reads_f270_r210.qza --o-denoising-stats denoising-stats_above_5000_f270_r210.qza --output-dir /ess/p2005/cluster/skin-metagenomics/qiime2/all_samples/all_above_5000_reads/DADA2_above_5000_reads_f270_r210 --p-max-ee-f 2 --p-max-ee-r 2 --verbose 

echo "I've done at "
date



