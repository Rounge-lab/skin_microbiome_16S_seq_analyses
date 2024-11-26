#Loading libraries
library(dplyr)
library(tidyr)
library(phyloseq) 
library(ggplot2)
library(decontam)

#load phyloseq object
load("/PATH/ps_f270_r210_non_normalized_allsamples.RData") #non-rarefied phyloseq object with all samples, including negative controls
ps_all = subset_samples(ps_all, type!="mock") %>% #remove mocks, they have 100% prevalence, not relevant to include
  prune_taxa(taxa_sums(.) != 0, .) #keep only taxa with non-zero sums

#decontam - find contaminants from all negative controls, threshold 0.5 removes all contaminant ASVS more prevalent in controls than in biological samples
all_contam = isContaminant(ps_all, method = "prevalence", neg = sample_data(ps_all)$group=="negative", 
  batch="extraction_batch", threshold = 0.5, normalize=TRUE, detailed = TRUE) 
table(all_contam$contaminant) #number of non-contaminants vs contaminants

#Contaminants and their relative abundance in each sample
ps_contam = prune_taxa(all_contam$contaminant, ps_all) %>% # phyloseq object of contaminants
  prune_taxa(taxa_sums(.) != 0, .) %>%
  subset_samples(type == "sample") #exclude the negatives
tax_table(ps_contam)[, c(2,6,7)] #taxonomy of contaminants
ps_contam %>% transform_sample_counts(function (x) 100*x/sum(x)) %>% psmelt() %>% 
  ggplot(aes(x=Sample, y=Abundance, fill=Genus)) + geom_bar(stat="identity") + facet_grid(OTU ~ skinsite, scales="free") 

#Distribution score plot (Fig S3A)
all_contam$prev_factor = cut(all_contam$prev, breaks=c(1.9, 3,6,9, Inf), labels=c("2", "3-6", "7-9", "10+"), include.lowest = F) #split prevalence into chosen ranges for coloring
score_plot = all_contam %>% filter(!is.na(p.prev)) %>% 
  ggplot(aes(x=p.prev, fill=prev_factor)) + geom_histogram(bins=20) + 
  scale_fill_manual(values=c("darkseagreen2", "palegreen3", "mediumseagreen", "seagreen", "darkgreen")) + 
  labs(x="Score Statistic", y="Frequency", fill="Prevalence") + theme_bw()

#Prevalence plot (Fig S3B) -  Transform counts to presence/absence and split samples
ps.pa = transform_sample_counts(ps_all, function(abund) ifelse(abund > 0, 1, 0)) #present/absent
ps.pa.neg = prune_samples(sample_data(ps.pa)$group == "negative", ps.pa) #negative control samples
ps.pa.pos = prune_samples(sample_data(ps.pa)$group != "negative", ps.pa) #true biological samples
#dataframe of prevalence in positive and negative samples
df.pa = data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=all_contam$contaminant)
prevalence_plot = ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative controls)") + 
  ylab("Prevalence (Biological samples)") + scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9)) + theme_bw() + labs(color="Contaminant")

#Combine score plot with prevalence plot (Fig S3)
decontam_plot = ggarrange(score_plot, prevalence_plot, labels=c("A", "B"))
ggsave(filename="FigureS3", plot=decontam_plot, width=8, height=3.2)


#Clean phyloseq object (Remove identified contaminants)
ps_soldiers = subset_samples(ps_all, type=="sample") %>% #soldiers only
  prune_taxa(taxa_sums(.) != 0, .) 
ps_soldiers_paper1_clean = prune_taxa(!taxa_names(ps_soldiers) %in% taxa_names(ps_contam), ps_soldiers) #contaminants removed

save(ps_soldiers_paper1_clean, file="ps_sol_paper1_clean.RData") #to be used for rarefaction, taxonomy plotting, DA etc.


sessionInfo()
