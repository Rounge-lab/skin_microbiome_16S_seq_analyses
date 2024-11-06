#!/bin/bash
#SBATCH --job-name=filter-features
#SBATCH --account=p2005_tsd
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10M
#SBATCH --nodes=1
#SBATCH --output=filtering%j.out

echo "Hello " $USER
echo "my submit directory is"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID

echo "Today is:"
date

echo "Filtering ASV table and representative sequences based on features.."

# Assign input arguments to variables
asv_input="$1"
taxonomy_input="$2"
sequence_input="$3"
output_folder="$4"

#Create temporary directory
temp_dir=$(mktemp -d)

# Step 1: Filter ASV table based on taxonomy
qiime taxa filter-table --i-table "$asv_input" --i-taxonomy "$taxonomy_input" --p-mode contains --p-exclude mitochondria,chloroplast --o-filtered-table "$temp_dir/asv_table.qza" --p-verbose

#Step 2: Filter features based on frequency, and remove singletons
qiime feature-table filter-features --i-table "$temp_dir/asv_table.qza" --p-min-frequency 10 --p-min-samples 2 --o-filtered-table "$output_folder/asv_filtered.qza" --p-verbose

#Step 3: Filter sequences based on the filtered ASV table
qiime feature-table filter-seqs --i-data "$sequence_input" --i-table "$output_folder/asv_filtered.qza" --o-filtered-data "$output_folder/rep-seqs-filtered.qza" --p-verbose

#Remove temporary directory
rm -r "$temp_dir"

echo "Filtering complete. Output files are in the '$output_folder' directory."




echo "I was done at:"

date



