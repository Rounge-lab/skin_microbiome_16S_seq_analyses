#Script used to load files from qiime2 into a phyloseq object in R

#Install packages if not already installed
#install.packages("Biostrings")
#install.packages("phyloseq")
#Load necessary packages
library(Biostrings)
library(phyloseq)
library(readr)

## First need to have results from Qiime2 in the correct input format for phyloseq: feature table (.tsv), representative sequences (.fasta), taxonomy file (.tsv), metadata file (.tsv)

#Load ASV feature table
asv_table = read.table("/PATH/asv_filtered_f270_r210.tsv", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
count_table = as.data.frame(asv_table)

#match ASV IDs between count table and taxonomy file
asv_ids = rownames(count_table)
taxonomy = read.table("/PATH/taxonomy_clean.tsv", header = TRUE, row.names=1, sep = "\t") #load pre-cleaned taxonomy file
tax.clean = as.data.frame(taxonomy)
tax_matched = tax.clean[asv_ids, ] #keep only taxa present in count table

#import metadata .tsv file
meta = import_qiime_sample_data("/PATH/metadata_n141.tsv")

#change order of sampleids for count_table to match with metadata
matching_sample_ids = intersect(rownames(meta), colnames(count_table))
count_table_new = count_table[, matching_sample_ids]

#import all representative sequences (including negatives)
rep.seqs = Biostrings::readDNAStringSet("/PATH/sequences.fasta", format="fasta") #representative sequences

#create full phyloseq object for all study samples to be used for Decontam()
asv.table = phyloseq::otu_table(count_table_new, taxa_are_rows = TRUE)
ps_all = phyloseq(asv.table, tax_table(as.matrix(tax_matched)), rep.seqs, sample_data(meta))
#taxa_names(ps_all) = paste0("ASV", seq(ntaxa(ps_all))) #rename from ASV IDs to no ASV1, ASV2 etc
save(ps_all, file="/PATH/ps_f270_r210_non_normalized_allsamples.RData") #save phyloseq object

#Check physeq object
nsamples(ps_all) #samples
ntaxa(ps_all) #ASVs 
sample_names(ps_all)[1:5] #sample names - first 5
rank_names(ps_all) #taxonmic ranks: "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
sample_variables(ps_all) #metadata variables
otu_table(ps_all)[1:5, 1:5] # first 5 ASVs and first 5 samples
tax_table(ps_all)[1:5, 4:7] #first 5 samples and taxonomy order, family, genus and species
top10 = names(sort(taxa_sums(ps_all), decreasing = TRUE)[1:10]) #top 10 ASVs
tax_table(ps_all)[top10][1:10,c(2,6,7)] #the top 10 ASVs, showing phylum, genus and species


  