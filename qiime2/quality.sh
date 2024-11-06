#!/bin/bash
#SBATCH --job-name=qualitycheck
#SBATCH --account=p2005_tsd
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1
#SBATCH --output=quality_soldiers%j.out

echo "Hello " $USER
echo "my submit directory is"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID

echo "Today is:"
date

echo "Summarizing read quality using qiime demux summarize.."

#/usr/bin/time -v 
qiime demux summarize \
  --i-data 16s_trimmed_soldiers.qza \
  --o-visualization 16s_trimmed_soldiers_QC.qzv


echo "I've done at "
date



