# Load libraries ---------------------------------------------------------------------------------------------------------------------------------------------
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(phyloseq) 
library(grid)
library(gridExtra)
library(ggpubr)
library(vegan)
library(nlme)
library(emmeans)
library(plyr)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------

# Import already cleaned and rarefied phyloseq objects for each site
load("hands_paper1_clean_rarefied.RData")
load("forearms_paper1_clean_rarefied.RData")

## ALPHA DIVERSITY -------------------------------------------------------------------------------------------------------------------------------------------

# List of diversity metrics
richness_metrics <- c("Observed", "Shannon", "InvSimpson")

# Function to create a data frame with diversity metrics
create_alphadiv_df <- function(physeq_obj, metrics) {
  richness_data <- lapply(metrics, function(m) estimate_richness(physeq_obj, measures = m))
  richness_df <- data.frame(
    sampleid = sample_data(physeq_obj)$sampleid,
    subjectid = sample_data(physeq_obj)$subjectid,
    sample_round = factor(sample_data(physeq_obj)$sample_round),
    batch = factor(sample_data(physeq_obj)$extraction_batch)
  )
  # Add richness metrics to the data frame
  alphadiv_df <- cbind(richness_df, setNames(richness_data, metrics))
  return(alphadiv_df)
}

# Create data frames with diversity metrics
all_rounds_hands <- create_alphadiv_df(hands_clean_rare, richness_metrics) #hands alpha div
all_rounds_forearms <- create_alphadiv_df(f_clean_rare, richness_metrics) #forearms alpha div

## STATISTICS ## 

# Function to fit linear mixed effects models for a given richness metric and design variable
fit_mixed_effects_models <- function(df, metrics, design_var) {
  models <- list()
  emmeans_list <- list()
  
  for (metric in metrics) {
    response_var <- if (metric == "InvSimpson") "Inverse_Simpson" else metric  # Handle naming inconsistency
    formula <- as.formula(paste(response_var, "~", design_var, "+ batch"))
    model_name <- paste0("model_h_", tolower(substr(response_var, 1, 3)))
    emmeans_name <- paste0("emmeans_", tolower(substr(response_var, 1, 3)))

    models[[model_name]] <- lme(formula, random = ~ 1 + get(design_var) | subjectid, data = df, method = "REML") #fit model, including subjectid as random effect
    emmeans_list[[emmeans_name]] <- emmeans(models[[model_name]], specs = pairwise ~ design_var) #estimated marginal means and contrast between each sample round pair in the linear model
    
    # Print summary of the model
    print(paste("Summary for model:", model_name))
    print(summary(models[[model_name]]))
  }
  
  # Return models and emmeans
  return(list(models = models, emmeans_list = emmeans_list))
}

# Get the linear mixed effects models and emmeans for each skin site
lme_hands <- fit_mixed_effects_models(all_rounds_hands, richness_metrics, "sample_round")
lme_forearms <- fit_mixed_effects_models(all_rounds_forearms, richness_metrics, "sample_round")
models_hands <- lme_hands$models
emmeans_hands <- lme_hands$emmeans_list
models_forearms <- lme_forearms$models
emmeans_forearms <- lme_forearms$emmeans_list

                       
## PLOTTING ## 

# Function to create alpha diversity boxplots for each sample round
create_alpha_div_plot <- function(data, scales, legend_pos) {
  data %>% 
  gather(key=metric, value=value, c("Shannon", "InvSimpson")) %>%
  mutate(metric=factor(metric, levels = c("Shannon", "InvSimpson"))) %>%
  ggplot(aes(x=sample_round, y=value), fill = sample_round) +
  geom_boxplot(outlier.color="grey", aes(group = sample_round, fill= sample_round)) + 
  geom_jitter(height=0, width=0) +
  labs(title="", x="", y = "") + 
  facet_wrap(~ metric, scales = scales, labeller = labeller(metric = c("InvSimpson"="Inverse Simpson"))) +
  theme_bw() + 
  theme(
    legend.position=legend_pos, 
    axis.text.x = element_blank(), 
    legend.text = element_text(size=14), 
    legend.title = element_text(size=14), 
    strip.text=element_text(size=12), 
    axis.text = element_text(size=10)
  ) + 
  scale_fill_manual(
    values=c("grey80", "deepskyblue4", "darkgoldenrod3"), 
    name = "Sample round", 
    labels = c("Baseline", "Post exercise", "3W Post exercise")
  ) +
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))
}

  
# Alpha div figures for each skinsite over sample rounds
fig_all_hands <- create_alpha_div_plot(all_round_hands, scales="free", legend_pos="none") 
fig_all_forearm <- create_alpha_div_plot(all_round_forearms, scales="free_y", legend_pos="bottom")

# Combine plot for hands and forearm (Fig. 3)
alphadiv_sha_inv <- ggarrange(fig_all_hands, fig_all_forearm, labels=c("A", "B"), ncol=1, nrow=2, common.legend = TRUE, legend="bottom", align = "hv") + 
                      theme(legend.position = "bottom")
ggsave(alphadiv_sha_inv, file="Figure3.pdf", width=210, units="mm", height=150) #A4 width


## SUPPLEMENTARY Fig. S7 A

# Alpha div hands vs forearms data frame
load("ps_sol_paper1_clean_rarefied.RData") #phyloseq object with both skin sites, if not already loaded
alpha_both_sites = data.frame("sampleid" = phyloseq::sample_data(ps_sol_paper1_clean_rarefied)$sampleid, 
                       "subjectid" = phyloseq::sample_data(ps_sol_paper1_clean_rarefied)$subjectid, 
                       "Observed" = phyloseq::estimate_richness(ps_sol_paper1_clean_rarefied, measures = "Observed"),
                       "Shannon" = phyloseq::estimate_richness(ps_sol_paper1_clean_rarefied, measures = "Shannon"),
                       "Inverse_Simpson" = phyloseq::estimate_richness(ps_sol_paper1_clean_rarefied, measures = "InvSimpson"), 
                       "sample_round" = factor(phyloseq::sample_data(ps_sol_paper1_clean_rarefied)$sample_round),
                       "batch" = factor(phyloseq::sample_data(ps_sol_paper1_clean_rarefied)$extraction_batch),
                       "skinsite" = factor(phyloseq::sample_data(ps_sol_paper1_clean_rarefied)$skinsite))

# Plot baseline alpha div hands vs forearms
alpha_both_sites %>% filter(sample_round==1) -> baseline_data
baseline_data %>% 
  gather(key=metric, value=value, c("Observed", "Shannon", "InvSimpson")) %>%
  mutate(metric=factor(metric, levels = c("Observed", "Shannon", "InvSimpson"))) %>%
  ggplot(aes(x=skinsite, y=value), fill = skinsite) + 
  geom_boxplot(outlier.color="grey", aes(group = skinsite, fill= skinsite)) +
  geom_jitter(height=0, width=0) + 
  facet_wrap(~ metric, scales = "free") +
  labs(fill="Skin site", y="Alpha diversity", x="") +
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_discrete(labels=c("hands", "forearms")) + 
  scale_y_continuous(expand=expansion(mult=c(0.05,0.05))) -> alpha_skinsite_base
ggsave(alpha_skinsite_base, file="FigureS7A.pdf")

# Testing alpha div difference between sites for each metric
for metric in richness_metrics {
  # Fit the linear mixed effects model using dynamically created formula
  model_skinsite <- lme(as.formula(paste(metric, "~ skinsite + batch")), random = ~ 1 | subjectid, data = baseline_data, method = "REML")
  
  # Extract and print the model summary (p-values)
  summary_res <- summary(model_skinsite)$tTable
  cat("Summary for metric:", metric, "\n")
  print(summary_res)
  
  # Get and print estimated marginal means
  emmeans_res <- emmeans(model_skinsite, specs = pairwise ~ skinsite)
  cat("\nEstimated marginal means for metric:", metric, "\n")
  print(emmeans_res)
  cat("\n") # Add a new line for better readability
}


# --------------------------------------------------------------------------------------------------------------------------------------------------------------


## BETA DIVERSITY ----------------------------------------------------------------------------------------------------------------------------------------------

# Function to calculate Bray-Curtis distance matrix
calculate_bray_distance <- function(physeq_obj) {
  phyloseq::distance(physeq_obj, method = "bray")
}

# Function to perform PERMANOVA
perform_permanova <- function(distance_matrix, formula, data, by="terms") {
  test.adonis <- adonis2(distance_matrix ~ get(formula), data = data, permutations = 999, by = by) # by="margin" tests each term against each other
  adonis_adjusted <- p.adjust(test.adonis$`Pr(>F)`, method = "BH") # Benjamini-Hochberg method for p-value adjustment
  list(test.adonis = test.adonis, adonis_adjusted = adonis_adjusted) #R² is the percentage of variance explained by the rounds
}
                          
# Calculate distance matrices
bray_all <- calculate_bray_distance(ps_sol_paper1_clean_rarefied) # full matrix
bray_hands <- calculate_bray_distance(hands_clean_rare) # hand samples only
bray_forearms <- calculate_bray_distance(f_clean_rare) # forearm samples only

# Perform PERMANOVA
adonis_all <- perform_permanova(bray_all, 'subjectid + sample_round + skinsite + extraction_batch', data.frame(sample_data(ps_sol_paper1_clean_rarefied)), by="margin") # Test all factors, how much they influence microbiome composition
adonis_hands <- perform_permanova(bray_hands, 'subjectid + sample_round + extraction_batch', data.frame(sample_data(hands_clean_rare)), by="margin")
adonis_forearms <- perform_permanova(bray_forearms, 'subjectid + sample_round + extraction_batch', data.frame(sample_data(f_clean_rare)), by="margin")

# Test for batch effects using permanova
adonis_batch_hands <- perform_permanova(bray_hands, 'extraction_batch', data.frame(sample_data(hands_clean_rare))) #batch effects hands
adonis_batch_forearms <- perform_permanova(bray_forearms, 'extraction_batch', data.frame(sample_data(f_clean_rare))) #batch effects forearms
                          

## PCOA plotting

# Function to create PCoA plot
create_betaplot <- function(ps, pcoa_data, annotation_text) {
  plot_ordination(ps, pcoa_data, color = "sample_round") + 
    geom_line(aes(group=subjectid), color="darkgray", lty="dashed") +
    stat_ellipse() + 
    scale_color_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name = "Sample round", labels = c("Baseline", "Post Exercise", "3W Post Exercise")) + 
    theme_bw() + 
    geom_point(size=2) +
    annotation_custom(grobTree(textGrob(annotation_text, x=0.02, y=0.96, hjust=0, gp=gpar(col="black", fontsize=8))))
}

# Calculate PCoA results
# pcoa_all = ordinate(ps_sol_paper1_clean_rarefied, method="PCoA", distance=bray, formula = ~ sample_round + skinsite + extraction_batch)
pcoa_hands = ordinate(hands_clean_rare, method="PCoA", distance=bray_hands, formula = ~ sample_round + extraction_batch) 
pcoa_forearms = ordinate(f_clean_rare, method="PCoA", distance=bray_forearms, formula = ~ sample_round + extraction_batch)

# Create PCoA plot for each skin site, #statistics found using permanova above
beta_hands <- create_betaplot(hands_clean_rare, pcoa_hands, "R = 0.429, p < 0.001") 
beta_forearms <- create_betaplot(f_clean_rare, pcoa_forearms, "R = 0.458, p < 0.001")

# Combined beta div plot (Fig. 4)
beta_comb = ggarrange(beta_hands, beta_forearms, labels=c("A", "B"), ncol=2, nrow=1, common.legend=TRUE, legend="bottom", align="hv") + geom_point(size=1)
ggsave(beta_comb, file="Figure4.pdf", width=7.5, height=4)


## Batch comparison

# Function to create PCOA plot for batch comparison
create_beta_batch_plot <- function(clean_data, pcoa_data, annotation_text) {
  plot_ordination(clean_data, pcoa_data, color = "extraction_batch") + 
    stat_ellipse() + 
    labs(color = "Batch") +
    scale_color_brewer(palette = "Dark2") + 
    geom_point(size = 2) + 
    theme_bw() + 
    annotation_custom(grobTree(textGrob(annotation_text, x = 0.03, y = 0.98, hjust = 0, gp = gpar(col = "black", fontsize = 8))))
}

# Create beta plot for batch comparison for hands and forearms (Fig. S6)
beta_hands_batch <- create_beta_batch_plot(hands_clean_rare, pcoa_hands, "R² = 0.05568, p = 0.022")
beta_forearm_batch <- create_beta_batch_plot(f_clean_rare, pcoa_forearms, "R² = 0.0755, p = 0.187")
# Combine plots
beta_comb_batch <- ggarrange(beta_hands_batch, beta_forearm_batch, labels = c("A", "B"), ncol = 2, nrow = 1, common.legend = TRUE, legend = "right", align = "hv")
ggsave(beta_comb_batch, file="FigureS6.pdf", width=7.5, height=4)
                          


# SUPPLEMENTARY Beta div hands vs forearms (Fig. S7)

baseline_all = subset_samples(ps_sol_paper1_clean_rarefied, sample_round==1) %>% prune_taxa(taxa_sums(.)>0, .)
bray_base = phyloseq::distance(baseline_all, method="bray")
pcoa_base = ordinate(baseline_all, method="PCoA", distance=bray_base, formula = ~skinsite+extraction_batch)
bray_base_skinsite <- plot_ordination(baseline_all, pcoa_base, color="skinsite") + 
                          stat_ellipse() + theme_bw() + geom_point(size=3) + labs(color="Skin site") + scale_color_discrete(labels=c("hands", "forearms")) 

# permanova test for statistical difference
adonis_skinsite <- perform_permanova(bray_base, 'skinsite + subjectid + extraction_batch', data.frame(sample_data(baseline_all)), by="margin")
bray_base_skinsite <- bray_base_skinsite + annotation_custom(grobTree(textGrob("R² = 0.08313, p = 0.003", x=0.03, y=0.06, hjust=0, gp=gpar(col="black", fontsize=10))))

# Combined alpha and beta div for skin site comparison baseline (Fig. S7)
ggarrange(alpha_skinsite_base, bray_base_skinsite, labels=c("A", "B"), ncol=1, nrow=2) -> base_hand_v_forearm
ggsave(base_hand_v_forearm, filename = "FigureS7", height=5, width=6)             


## Bray-Curtis distances within (intra-individual) vs between (inter-individual)
                          
# Convert Bray-Curtis distance to matrix and dataframe
bray_mat <- as.data.frame(as.matrix(bray_all))
rownames(bray_mat) <- as.numeric(rownames(bray_mat))

# Extract metadata
meta_sol_paper1 <- data.frame(sample_data(ps_sol_paper1_clean_rarefied)) %>%
  mutate(sampleid = as.character(sampleid))

# Function to compute within vs. between distances between all pairs of samples
compute_within_between <- function(bray_matrix, meta_data) {
  bray_matrix %>% 
    rownames_to_column("sample_id_1") %>%  # Convert rownames to a column for use in joins
    pivot_longer(-sample_id_1, names_to = "sample_id_2", values_to = "dist") %>% 
    left_join(meta_sol_paper1 %>% select( # Join metadata to the first sample_id
    sample_id_1 = sampleid, subject_id_1 = subjectid, sample_round_1 = sample_round, skinsite_1 = skinsite)) %>%
    left_join(meta_sol_paper1 %>% select( # Join metadata to the second sample_id
    sample_id_2 = sampleid, subject_id_2 = subjectid, sample_round_2 = sample_round, skinsite_2 = skinsite)) %>% 
    mutate(
      same_individual = ifelse(subject_id_1 == subject_id_2, "Within", "Between"), # Determine if samples are from the same individual
      same_site = ifelse(skinsite_1 == skinsite_2, "Same_skinsite", "Hands vs Forearms"), # Determine if samples are from the same skin site
      same_time = ifelse(sample_round_1 == sample_round_2, "Same_time", "Over time") # Determine if samples are from the same time point
    ) %>%
    filter(sample_id_1 != sample_id_2) %>% # Remove distances to the same sample
    distinct(pair_id = paste(pmin(sample_id_1, sample_id_2), pmax(sample_id_1, sample_id_2), sep = "_"), .keep_all = TRUE) # Remove duplicate pairs of samples
}

# Data frame with computed within/between BC-distances                          
inter_intra_distances <- compute_within_between(bray_mat, meta_sol_paper1)
saveRDS(inter_intra_distances, file="inter_intra_distances.rds") #for later use
                          
# Function to create Bray-Curtis boxplots
create_BC_boxplot <- function(data) {
  data %>%
    group_by(same_individual, same_site, subject_id_1, skinsite_1) %>%  # Group by key variables
    summarize(mean_dist = mean(dist)) %>%  # Calculate mean distances
    ggplot(aes(x = reorder(same_individual, mean_dist), y = mean_dist, fill = same_individual)) +
    geom_boxplot() +  
    labs(y = "Bray-Curtis distance") +  # y-axis label
    scale_fill_manual(values = c("#9e4f4a", "#ecb775")) +  # Define fill colors
    theme_bw() + # Use a clean theme
    theme(text = element_text(size = 16), legend.position = "none")
}
                          
# Same skin site plot
intra_inter_samesite_plot <- create_BC_boxplot(inter_intra_distances) + 
  facet_wrap(~ skinsite_1, scales = "free", labeller = labeller(
    skinsite_1 = function(var) tools::toTitleCase(as.character(var)))) #facet labels are capitalized
                                                                
# Hands vs forearms plot
bray_h_vs_f <- inter_intra_distances %>% filter(skinsite_1 != skinsite_2)
intra_inter_h_vs_f_plot <- create_BC_boxplot(bray_h_vs_f) + 
  facet_wrap(~ same_site) +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

# Combine the plots side by side
intra_inter_beta_plot <- ggarrange(
  intra_inter_samesite_plot + rremove("xlab"),
  intra_inter_h_vs_f_plot + rremove("ylab") + rremove("xlab"),
  widths = c(2, 1)
)

ggsave(intra_inter_beta_plot, file="intra_inter_betaplot.pdf", width=11, height=5)


## Figure 5 - inter- and intraindividual variation
#load("inter_intra_distances.rds") #load or read data set if not already loaded
#inter_intra_distances <- readRDS("PATH/inter_intra_distances.rds")

#  Bray-Curtis distance between samples from the same individual over time (intra-individual), comparing skin sites (Fig 5A)
fig5_a <- inter_intra_distances %>% as.data.frame() %>%
  filter(same_site == "Same_skinsite" & same_time == "Over time" & same_individual == "Within") %>%
  group_by(sample_round_1) %>%
  ggplot(aes(x = reorder(skinsite_1, dist), y = dist, fill = skinsite_1)) + 
  geom_boxplot() + 
  scale_fill_discrete(name = "Skin site", labels = c("hands", "forearms")) + 
  labs(x = "", y = "Bray-Curtis distance (intra-individual)") + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  ylim(0.35, 0.97)

# Bray-Curtis distance between hands and forearms (H-to-F) within each individual (Fig 5B)
fig5_b <- inter_intra_distances %>% as.data.frame() %>%
  filter(same_site == "Hands vs Forearms" & same_time == "Same_time" & same_individual == "Within") %>%
  group_by(sample_round_1) %>%
  ggplot(aes(x = sample_round_1, y = dist, fill = sample_round_1)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name = "Sample round", labels = c("Baseline", "Post Exercise", "3W Post Exercise")) + 
  labs(x="", y="Bray-Curtis distance (H-to-F)") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  ylim(0.35, 0.97)

# Between-distance over time, for each skinsite (Fig 5C)
fig5_c <- inter_intra_distances %>% as.data.frame() %>% 
    filter(same_site=="Same_skinsite" & same_individual=="Between" & same_time=="Same_time") %>% 
    ggplot(aes(x=sample_round_1, y=dist, fill=sample_round_1)) + 
    geom_boxplot() + 
    scale_fill_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name="Sample round", labels=c("Baseline", "Post exercise", "3W Post exercise")) + 
    labs(x="", y="Bray-Curtis distance (inter-individual)") + 
    theme_bw() + 
    facet_wrap(~skinsite_1, labeller=labeller(skinsite_1 = c(hands="hands", forearm="forearms"))) + 
    ylim(0.35, 0.97) + 
    theme(axis.text.x = element_blank())
saveRDS(fig5_c, file="inter_dist_boxes.rds")
#fig5_c <- readRDS("PATH/inter_dist_boxes.rds") 


## Beta dispersion (soldier similarity over time)

analyze_dispersion_and_anosim <- function(distance_matrix, sample_data_frame, site_name) {
    sample_round <- sample_data_frame$sample_round
    dispersion <- betadisper(distance_matrix, sample_round, type = "centroid") # Compute dispersion  
    permutest_result <- permutest(dispersion, pairwise = TRUE, permutations = 999) # Perform permutest
    
    # Perform ANOSIM
    anosim_result <- anosim(distance_matrix, sample_round, permutations = 999)
    anosim_stat <- anosim_result$statistic #R value range between -1 and 1. Positive values suggest that there is more similarity within groups than between groups.
    anosim_signif <- anosim_result$signif
    
    # ANOVA and TukeyHSD
    anova_test <- anova(dispersion)
    tukey_result <- TukeyHSD(dispersion)
    
    # Visualization of dispersion (optional)
    plot(dispersion)
    boxplot(dispersion, xlab = "Sample round") #distribution of distances to centroid for each round - how similar are the soldiers to each other?
    
    # Prepare a data frame for betadisp distances
    beta_disp_df <- data.frame(
        site = site_name, 
        round = sample_round,
        beta_disp = dispersion$distances
    )
    
    return(list(
        dispersion = dispersion, 
        permutest_result = permutest_result, 
        anosim_stat = anosim_stat, 
        anosim_signif = anosim_signif, 
        anova_test = anova_test, 
        tukey_result = tukey_result, 
        beta_disp_df = beta_disp_df
    ))
}

# Beta dispersion for forearms
forearms_results <- analyze_dispersion_and_anosim(bray_forearms, data.frame(sample_data(f_clean_rare)), "forearms")
#print(forearms_results$anosim_stat)
#print(forearms_results$anosim_signif)
#print(forearms_results$anova_test)
#print(forearms_results$tukey_result)

# Beta dispersion for forearms hands
hands_results <- analyze_dispersion_and_anosim(bray_hands, data.frame(sample_data(hands_clean_rare)), "hands")

# Combine data frames
disp_combined <- rbind(forearms_results$beta_disp_df, hands_results$beta_disp_df)

# Beta dispersion (Fig 5D)
fig5_d <- ggplot(disp_combined, aes(x=round, y=beta_disp, fill=round)) +
  geom_boxplot() +
  facet_wrap(~reorder(site,beta_disp)) +
  labs(x="", y="Distance to centroid") + 
  theme_bw() +
  scale_fill_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name = "Sample round", labels = c("Baseline", "Post Exercise", "3W Post Exercise")) +
  ylim(0.35, 0.97) + 
  theme(axis.text.x = element_blank())
ggsave(fig5_d, filename="beta_disp_boxplot.pdf")
saveRDS(betadisp_boxes, file="beta_disp_boxplot.rds")
#fig5_d <- readRDS("PATH/beta_disp_boxplot.rds")

# Figure 5 - Combine all plots to figure
fig5_all <- ggarrange(fig5_a + theme(legend.position = "right", legend.text = element_text(size=11), legend.title=element_text(size=12)), 
                     fig5_b + theme(legend.position = "right", legend.text = element_text(size=11), legend.title=element_text(size=12)), 
                     fig5_c + theme(legend.position = "none", strip.text=element_text(size=11)), 
                     fig5_d + theme(legend.position = "none", strip.text=element_text(size=11)), 
                     ncol=2, nrow=2, labels=c("A", "B", "C", "D"))
ggsave(fig5_all, dpi=300, filename = "Figure5.pdf", useDingbats=FALSE, width=7, height=7, units="in")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

sessionInfo()

