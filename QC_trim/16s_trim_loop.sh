#!/bin/bash
# Job name:
#SBATCH --job-name=trim-16S-loop
# Project:
#SBATCH --account=p2005_tsd
# Wall clock limit:
#SBATCH --time=12:00:00
# Memory usage:
#SBATCH --mem-per-cpu=1G
## Cpu
#SBATCH --cpus-per-task=6
##Computers
#SBATCH --nodes=1

#SBATCH --output=trimmed-loop-q20-l150%j.out

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
module load fastqc/0.11.8
module load cutadapt/4.2-GCCcore-11.3.0
module load Trim_Galore/0.6.6-GCCcore-9.3.0-Python-3.8.2

#Script that first does QC, then primer removal with cutadapt, then trimgalore to remove adapters and low-quality reads (including another fastQC).

#Go to folder where raw fastq files are
cd /ess/p2005/cluster/skin-metagenomics/data/raw/seqdata/16S/230323_M07166.Project_Glenna-16S1-2023-03-15

# Define primer sequences (here 16S V3-V4 primers)
forward_primer="CCTACGGGNGGCWGCAG"
reverse_primer="GACTACHVGGGTATCTAATCC"

#Define output folder for where trimmed files should end up
trimmed_folder="/cluster/projects/p2005/skin-metagenomics/data/trimmed_files/16S/noprimers_nextera_q20_l150"
cutadapt_output="/cluster/projects/p2005/skin-metagenomics/data/trimmed_files/16S/noprimers_nextera_q25_l250/cutadapt_output"

# Create directory to store (intermediate) cutadapt output
mkdir -p "${trimmed_folder}/cutadapt_output"

#LOOP through fastq files
for file in *R1_001.fastq.gz
do
    # Set name of sample (same for both reads)
    name=$(basename "${file}" | cut -d '_' -f 1) #sample name. Separates on underspace and chooses the first substring before "_"
    # Perform fastQC on raw seq files
    fastqc ${file} ${file/R1_001/R2_001}
    # Trim primers using Cutadapt
    cutadapt -g $forward_primer -G $reverse_primer --error-rate 0.1 --overlap 3 --info-file $trimmed_folder/cutadapt_output/${name}_cutadapt_report.txt -o $trimmed_folder/cutadapt_output/${name}_R1_trimmed.fq.gz -p $trimmed_folder/cutadapt_output/${name}_R2_trimmed.fq.gz ${file} ${file/R1_001/R2_001}
    # Trim the adapters and low-quality reads from the paired-end reads (results from cutadapt)
    trim_galore --paired --q 25 --cores 6 --length 150 --nextera --fastqc $trimmed_folder/cutadapt_output/${name}_R1_trimmed.fq.gz $trimmed_folder/cutadapt_output/${name}_R2_trimmed.fq.gz -o $trimmed_folder
done


echo "I was done at"
date

