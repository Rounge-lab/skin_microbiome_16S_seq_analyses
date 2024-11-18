# Loading libraries ---------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(vegan)
library(dplyr)
library(gtools)
# ----------------------------------------------------------------------------------------------------------------------------------------------------


#load data and extract ASV table for biological and negative samples
load("ps_f270_r210_non_normalized_paper1.RData") #ps_paper1 containing negative controls
load("/PATH/ps_sol_paper1_clean.RData") #clean soldier phyloseq object (after decontam)
asv_soldiers=otu_table(ps_soldiers_paper1_clean)
ps_negatives = subset_samples(ps_paper1, group=="negative") %>%
  prune_taxa(taxa_sums(.) != 0, .)
asv_neg = otu_table(ps_negatives)

#Goodness of fit plot to assess rarefaction depth
goodness_of_fit <- function(asv_table, rarefaction_depth) {
  S <- specnumber(t(asv_table)) #number of species
  asv_t <- t(asv_table) #Transpose the ASV table for rarefaction
  set.seed(111) #Set random seed for reproducibility
  rarefied <- rarefy(asv_t, rarefaction_depth) # Perform rarefaction
  fit <- lm(rarefied ~ S) #Fit linear model
  # Generate plot
  plot(S, rarefied, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "Goodness of Fit")
  abline(fit, col = "red")
  rsq <- summary(fit)$r.squared #Compute R squared
  text(220,50,paste0("R² = ", round(rsq, 2))) #Display R-squared on the plot
  
  return(list(
    model = fit,
    rsquared = rsq
  ))
}
#Run the goodness of fit-plotting function with chosen rarefaction depth
goodness_of_fit_result <- goodness_of_fit(asv_table, rarefaction_depth = 7000) #change ASV-table as desired


#function to calculate rarefaction curves
calculate_rarecurve = function(asv_table, color) {
  asv = as.data.frame(asv_table)
  asv = asv[, mixedsort(colnames(asv))] #sort columns from sample 1,2,3 etc
  asv_t = t(asv)
  S = specnumber(asv_t)
  min_seq_t=min(rowSums(asv_t))
  rarecurve(asv_t, step=100, col=color, sample = 7000, xlim=c(0,30000), tidy=TRUE)
  #abline(v=7000) #sampling depth rarefied 
}

#Individual rarefaction curves
rare_neg = calculate_rarecurve(asv_neg) #negative samples
rare_bio = calculate_rarecurve(asv_soldiers) #biological samples

#Combined curves for biological and negative samples
combined_rarecurve_data = dplyr::bind_rows(
  mutate(rare_bio, Samples = "Biological samples"), 
  mutate(rare_neg, Samples = "Negative controls")
)

#Plot rarefaction curves
rareplot_combined = ggplot(combined_rarecurve_data, aes(x=Sample, y = Species, group=Sample, color=Samples)) +
  geom_point(aes(color=Samples), size=0.3) +
  #xlim(0,30000) +
  theme_minimal() +
  geom_vline(xintercept=7000, linetype="dashed") +
  labs(y="ASVs", x="Reads") + 
  guides(colour=guide_legend(override.aes = list(size=2))) +
  scale_color_manual(
    values = c("Biological samples" = "black", "Negative controls" = "lightgrey"), name = "Sample type") +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=15), legend.title = element_text(size=13), legend.text = element_text(size=11))

#save plot
ggsave(plot=rareplot_combined, file="FigureS4.pdf", width=12)


### RAREFY clean phyloseq object ---------------------------------------------------------------------------------------------------------------------------
set.seed(111) #keep result reproducible
ps_sol_paper1_clean_rarefied = rarefy_even_depth(ps_soldiers_paper1_clean, sample.size = 7000, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

# Phyloseq objects by skin site
hands_clean_rare = subset_samples(ps_sol_paper1_clean_rarefied, skinsite=="hands") %>%
  prune_taxa(taxa_sums(.) != 0, .) %>%
save(hands_clean_rare, file="hands_paper1_clean_rarefied.RData")
f_clean_rare = subset_samples(ps_sol_paper1_clean_rarefied, skinsite=="forearm") %>%
  prune_taxa(taxa_sums(.) != 0, .) %>%
save(f_clean_rare, file="forearms_paper1_clean_rarefied.RData")

# -------------------------------------------------------------------------------------------------------------------------------------------------------

sessionInfo()

