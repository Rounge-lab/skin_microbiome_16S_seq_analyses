#Script using DESeq2 to find differentially abundant taxa on hands and forearms, contrasting each pairs of sample rounds.

#Loading libraries
library(readr)
library(dplyr)
library(tidyr)
library(phyloseq) 
library(vegan)
library(DESeq2)
library(ggplot2)
library(ggpubr)
require(grid)

#load non-normalized, clean phyloseq object
load("PATH/ps_sol_paper1_clean.RData")

#Hands
hands = subset_samples(ps_soldiers_paper1_clean, skinsite == "hands") %>% 
  prune_taxa(taxa_sums(.) != 0, .)
sample_sums(hands) #total number of ASVs per sample 

#Forearms
forearms = subset_samples(ps_soldiers_paper1_clean, skinsite == "forearm") %>% 
  prune_taxa(taxa_sums(.) != 0, .)
sample_sums(forearms) #total number of ASVs per sample 

#Filter ASV tables to fewer features before differential abundance calculation
hands_filtered = phyloseq::filter_taxa(hands, 
                                       function(x) sum(x >= 10) > (0.1*length(x)), prune=TRUE) #retains all taxa with count >=10 in at least 10 % of hand samples
forearms_filtered = phyloseq::filter_taxa(forearms, 
                                          function(x) sum(x >= 10) > (0.1*length(x)), prune=TRUE) #retains all taxa with count >=10 in at least 10 % of forearm samples
#30 best forearm samples
forearm_30 = prune_samples(sample_sums(forearms_filtered) > 7000, forearms_filtered) #keep only the samples that are kept after rarefaction

#Counts of each feature
summary(rowSums(otu_table(hands))) #hands before
summary(rowSums(otu_table(forearms))) #forearms before
summary(rowSums(otu_table(hands_filtered))) #hands after
summary(rowSums(otu_table(forearm_30))) #forearms after



##DESeq2 analysis##
# Ensure 'sample_round' is a factor and prepare contrasts
contrasts = list(
  c("sample_round",2,1),
  c("sample_round",3,2),
  c("sample_round",3,1)
)

#calculate geometric means prior to estimating size factors, since there are many zeros in the ASV tables
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x>0]), na.rm = na.rm) / length(x))
  } 

deseq_analysis = function(ps, contrasts) {
  results_list = list()
  #Iterate over the contrasts to get the DA results
  for (i in 1:length(contrasts)) {
    contrast = contrasts[[i]]
    print(paste("Processing contrast:", paste(contrast, collapse = " "))) #check
    sample_data(ps)$sample_round <- relevel(sample_data(ps)$sample_round, ref = as.character(contrast[3])) # Relevel sample_round based on the reference level specified in contrast
    #print(paste("Baseline level is now:", levels(sample_data(ps)$sample_round)[1])) #check
    dds = phyloseq_to_deseq2(ps, ~ subjectid + extraction_batch + sample_round) #convert phyloseq object to DESeqDataSet, comparative condition must be placed last in design formula
    geoMeans = apply(counts(dds), 1, gm_mean)
    dds = estimateSizeFactors(dds, geoMeans = geoMeans) 
    dds = DESeq(dds, test = "Wald", fitType="local", sfType="poscounts") #sfType: size factors for accurate normalization, using only positive counts, fitType: type of dispersion estimation, local is for each ASV individually to account for variable dispersions
    
    #Extract results and shrink log2FC-values
    coef_names = resultsNames(dds) #names of all the model variables (coefficients)
    #print(paste("Coefficient names:", paste(coef_names, collapse = ", "))) #if unsure of the coefficients
    res = results(dds, contrast = contrast, pAdjustMethod = "BH")
    coef_index <- which(coef_names == paste0("sample_round_", contrast[2], "_vs_", contrast[3]))
    res = lfcShrink(dds, coef = coef_index, res=res, type="apeglm") #apeglm shrinkage of log2FC
    final_table = cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))
    results_list[[i]] = final_table
  }
  return(results_list)
}

#Use DESeq2 function for each phyloseq object (hands and forearms)
shrunk_results_hands = deseq_analysis(hands_filtered, subjectid + extraction_batch + sample_round)
shrunk_results_forearms = deseq_analysis(forearm_30, subjectid + extraction_batch + sample_round)

# Define results for both hands and forearms and the corresponding filenames
shrunk_results_all <- list(hands = shrunk_results_hands, forearms = shrunk_results_forearms)
file_names_all <- list(
  hands = c("h_2v1_shrunk.rds", "h_3v2_shrunk.rds", "h_3v1_shrunk.rds"),
  forearms = c("f_2v1_shrunk.rds", "f_3v2_shrunk.rds", "f_3v1_shrunk.rds")
)
# Loop through the combined lists and save each result to the corresponding filename
for (site in names(shrunk_results_all)) {
  shrunk_results <- shrunk_results_all[[site]]
  file_names <- file_names_all[[site]]
  
  for (i in seq_along(shrunk_results)) {
    saveRDS(shrunk_results[[i]], file = file_names[i])
  }
}

 
sessionInfo()