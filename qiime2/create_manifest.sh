#!/bin/bash
#SBATCH --job-name=create_qiime_manifest
#SBATCH --account=p2005_tsd
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1M
#SBATCH --nodes=1
#SBATCH --output=manifest%j.out

echo "Hello " $USER
echo "my submit directory is"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID

echo "Today is:"
date

echo "Creating manifest file for qiime import..."
folder="/ess/p2005/cluster/skin-metagenomics/data/trimmed_files/16S/noprimers_nextera_q25_l150/tokeep" #using the trimmed 16S files, with q>25 and min_length=150, tokeep is all samples with >5000 trimmed reads
manifest_file="manifest_above_5000_reads.tsv" #output file to be created

#header of manifest file
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$manifest_file"

#create sample ID, fwd and rev filepaths for all fastq files in the folder
for file in "$folder"/*; do
    if [[ -f "$file" ]]; then
        filename=$(basename "$file")
        sample_id=$(echo "$filename" | grep -o '^[0-9]\+')
        
        if [[ $filename == *"R1"* ]]; then
        forward_filepath=$(realpath "$file")
        elif [[ $filename == *"R2"* ]]; then
        reverse_filepath=$(realpath "$file")
        echo -e "$sample_id\t$forward_filepath\t$reverse_filepath" #print info for the specific file to a new line
        fi
    fi
done >> "$manifest_file" #send all the printed lines to the manifest file 


echo "I've done at "
date



