#!/bin/bash
#SBATCH --job-name=import_to_qiime
#SBATCH --account=p2005_tsd
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1
#SBATCH --output=import_%j.out

echo "Hello " $USER
echo "my submit directory is"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID

echo "Today is:"

echo "Creating qiime import artifact..."

#/usr/bin/time -v 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_above_5000_reads_sorted.tsv \
  --output-path 16s_trimmed_above_5000_reads.qza \
  --input-format PairedEndFastqManifestPhred33V2


echo "I've done at "
date



