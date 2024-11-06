#!/bin/bash
#SBATCH --job-name=classify-sklearn
#SBATCH --account=p2005_tsd
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --nodes=1
#SBATCH --output=classification_above_5000_f270_r210_%j.out

echo "Hello " $USER
echo "my submit directory is"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID

echo "Today is:"
date

echo "Assigning taxonomy to features with qiime feature-classifier classify-sklearn (machine learning approach).."

#/usr/bin/time -v 
qiime feature-classifier classify-sklearn --p-n-jobs 8 --i-classifier /ess/p2005/cluster/skin-metagenomics/data/silva_reference/silva138_AB_V3-V4_classifier.qza --i-reads rep-seqs-all-above5000_f270_r210.qza --o-classification taxonomy_all_above_5000_reads_f270_r210.qza --output-dir taxonomy_all_above_5000_reads_f270_r210 --verbose

echo "I've done at "
date



