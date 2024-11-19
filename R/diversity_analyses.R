# Load libraries ---------------------------------------------------------------------------------------------------------------------------------------------
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(phyloseq) 
library(grid)
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
alphadiv_sha_inv <- ggarrange(fig_all_hands, fig_all_forearm, labels=c("A", "B"), ncol=1, nrow=2, common.legend = TRUE, legend="bottom", align = "hv") + theme(legend.position = "bottom")
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

# Test all factors, how much they influence microbiome composition (PERMANOVA)
bray = phyloseq::distance(ps_sol_paper1_clean_rarefied, method="bray") #dissimilarity matrix
test.adonis = adonis2(bray ~ subjectid + sample_round + skinsite + extraction_batch, data = data.frame(
  sample_data(ps_sol_paper1_clean_rarefied)), permutations = 999, by="margin") #by="margin" would test each term against each other
test.adonis
adonis_adjusted = p.adjust(test.adonis$`Pr(>F)`, method="BH")

bray_hands = phyloseq::distance(hands_clean_rare, method="bray") 
bray_forearms = phyloseq::distance(f_clean_rare, method="bray") 
#PERMANOVA
test.adonis = adonis2(bray_hands ~ subjectid + sample_round + extraction_batch, data = data.frame(sample_data(hands_clean_rare)), permutations = 999, by="margin") #by="margin" would test each term against each other
test.adonis #R² is the percentage of variance explained by the rounds
adonis_adjusted = p.adjust(test.adonis$`Pr(>F)`, method="BH")

#PCOA plot of beta diversity, including stat results
pcoa_hands = ordinate(hands_clean_rare, method="PCoA", distance=bray_hands, formula = ~ sample_round+extraction_batch) #PCoA results
pcoa_forearms = ordinate(f_clean_rare, method="PCoA", distance=bray_forearms, formula = ~ sample_round+extraction_batch) #PCoA results
pcoa_all = ordinate(ps_sol_paper1_clean_rarefied, method="PCoA", distance=bray, formula = ~sample_round + skinsite+extraction_batch)

betaplot_f = plot_ordination(f_clean_rare, pcoa_forearms, color = "sample_round") + 
  geom_line(aes(group=subjectid), color="darkgray", lty="dashed") +
  stat_ellipse() + 
  scale_color_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name = "Sample round", labels = c("Baseline", "Post Exercise", "3W Post Exercise")) + theme_bw() + geom_point(size=2)
betaplot_h = plot_ordination(hands_clean_rare, pcoa_hands, color = "sample_round") + 
  geom_line(aes(group=subjectid), color="darkgray", lty="dashed") +
  stat_ellipse() + 
  scale_color_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name = "Sample round", labels = c("Baseline", "Post Exercise", "3W Post Exercise")) + theme_bw() + geom_point(size=2)
#+ theme(strip.background = element_blank())
beta_hands = betaplot_h + annotation_custom(grobTree(textGrob("R = 0.429, p < 0.001", x=0.02, y=0.96, hjust=0, gp=gpar(col="black", fontsize=8))))
beta_forearms = betaplot_f + annotation_custom(grobTree(textGrob("R = 0.458, p < 0.001", x=0.02, y=0.96, hjust=0, gp=gpar(col="black", fontsize=8))))

#batch comparison
betaplot_batch_f = plot_ordination(f_clean_rare, pcoa_forearms, color = "extraction_batch") + 
  stat_ellipse() + 
  labs(color="Batch") +
  scale_color_brewer(palette="Dark2") + geom_point(size=2) + theme_bw()
test.adonis = adonis2(bray_forearms ~ extraction_batch, data = data.frame(sample_data(f_clean_rare)), permutations = 999) #by="margin" would test each term against each other
test.adonis #R² is the percentage of variance explained by the rounds
adonis_adjusted = p.adjust(test.adonis$`Pr(>F)`, method="BH")
beta_forearm_batch = betaplot_batch_f + annotation_custom(grobTree(textGrob("R² = 0.0755, p = 0.187", x=0.03, y=0.98, hjust=0, gp=gpar(col="black", fontsize=8))))

betaplot_batch_h= plot_ordination(hands_clean_rare, pcoa_hands, color = "extraction_batch") +
  stat_ellipse() + 
  labs(color="Batch") +
  scale_color_brewer(palette="Dark2") + geom_point(size=2) + theme_bw()
test.adonis = adonis2(bray_hands ~ extraction_batch, data = data.frame(sample_data(hands_clean_rare)), permutations = 999) #by="margin" would test each term against each other
test.adonis #R² is the percentage of variance explained by the rounds
adonis_adjusted = p.adjust(test.adonis$`Pr(>F)`, method="BH")
beta_hands_batch = betaplot_batch_h + annotation_custom(grobTree(textGrob("R² = 0.05568, p = 0.022", x=0.03, y=0.98, hjust=0, gp=gpar(col="black", fontsize=8))))
beta_comb_batch = ggarrange(beta_hands_batch, beta_forearm_batch, labels=c("A", "B"), ncol=2, nrow=1, common.legend=TRUE, legend="right", align="hv")

#combined beta div plot
beta_comb = ggarrange(beta_hands, beta_forearms, labels=c("A", "B"), ncol=2, nrow=1, common.legend=TRUE, legend="bottom", align="hv") + geom_point(size=2)
beta_comb = beta_comb + geom_point(size=1)
ggsave(beta_comb, file="beta_div_combined_paper1.pdf", width=7.5, height=4)

#Beta dispersion: check whether equal variance between group; if not use ANOSIM instead of permanova
f_dispersion = betadisper(bray_forearms, data.frame(sample_data(f_clean_rare))$sample_round, type = "centroid")
f_disp_test = permutest(f_dispersion, pairwise=TRUE, permutations=999)
anosim_results = anosim(bray_forearms, data.frame(sample_data(f_clean_rare))$sample_round, permutations = 999)
anosim_results$statistic #R value range between -1 and 1. Positive values suggest that there is more similarity within groups than between groups.
anosim_results$signif

anova_test_f = anova(f_dispersion) #test statistical difference in beta dispersion
TukeyHSD(f_dispersion) #which contrasts are different
plot(f_dispersion)
boxplot(f_dispersion, xlab="Sample round") #distribution of distances to centroid for each round - how similar are the soldiers to each other?

beta_disp_hands = data.frame(
  site = "hands", 
  round = data.frame(sample_data(hands_clean_rare))$sample_round,
  beta_disp = h_dispersion$distances
)
beta_disp_forearms = data.frame(
  site = "forearms", 
  round = data.frame(sample_data(f_clean_rare))$sample_round,
  beta_disp = f_dispersion$distances
)
disp_combined = rbind(beta_disp_hands, beta_disp_forearms)

#beta dispersion boxplot (soldier similarity over time)
betadisp_boxes = ggplot(disp_combined, aes(x=round, y=beta_disp, fill=round)) +
  geom_boxplot() +
  facet_wrap(~reorder(site,beta_disp)) +
  labs(x="", y="Distance to centroid") + 
  theme_bw() +
  scale_fill_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name = "Sample round", labels = c("Baseline", "Post Exercise", "3W Post Exercise")) 
ggsave(betadisp_boxes, filename="beta_disp_boxplot.pdf")
saveRDS(betadisp_boxes, file="beta_disp_boxplot.rds")

#beta div hands vs forearms (Fig. S7)
baseline_all = subset_samples(ps_sol_paper1_clean_rarefied, sample_round==1)
baseline_all = prune_taxa(taxa_sums(baseline_all)>0, baseline_all)
bray_base = phyloseq::distance(baseline_all, method="bray")
pcoa_base = ordinate(baseline_all, method="PCoA", distance=bray_base, formula = ~skinsite+extraction_batch)
plot_ordination(baseline_all, pcoa_base, color="skinsite") + stat_ellipse() + theme_bw() + geom_point(size=3) + labs(color="Skin site") + scale_color_discrete(labels=c("hands", "forearms")) -> bray_base_skinsite

#permanova test for statistical difference
test.adonis = adonis2(bray_base ~ skinsite + subjectid + extraction_batch, data = data.frame(sample_data(baseline_all)), permutations = 999, by="margin") #by="margin" would test each term against each other #by="margin" would test each term against each other
test.adonis #R² is the percentage of variance explained by the rounds
adonis_adjusted = p.adjust(test.adonis$`Pr(>F)`, method="BH")
bray_base_skinsite <- bray_base_skinsite + annotation_custom(grobTree(textGrob("R² = 0.08313, p = 0.003", x=0.03, y=0.06, hjust=0, gp=gpar(col="black", fontsize=10))))

# Combined alpha and beta div for skin site comparison baseline (Fig. S7)
ggarrange(alpha_skinsite_base, bray_base_skinsite, labels=c("A", "B"), ncol=1, nrow=2) -> base_hand_v_forearm
ggsave(base_hand_v_forearm, filename = "baseline_hand_vs_forearm.pdf", height=5, width=6)
             


#beta div intra vs inter
bray = phyloseq::distance(ps_sol_paper1_clean_rarefied, method="bray")
bray_mat = as.matrix(bray) %>% as.data.frame()
rownames(bray_mat) = as.numeric(rownames(bray_mat))
meta_sol_paper1 = data.frame(sample_data(ps_sol_paper1_clean_rarefied))
meta_sol_paper1$sampleid = as.character(meta_sol_paper1$sampleid)

#within vs between subject distance
bray_within_between_samesite = bray_mat %>% 
  rownames_to_column("sample_id_1") %>% 
  pivot_longer(-sample_id_1, names_to="sample_id_2", values_to = "dist") %>% 
  left_join(meta_sol_paper1 %>% select(
    sample_id_1 = sampleid, subject_id_1 = subjectid, sample_round_1 = sample_round, skinsite_1 = skinsite)) %>% 
  left_join(meta_sol_paper1 %>% select(
    sample_id_2 = sampleid, subject_id_2 = subjectid, sample_round_2 = sample_round, skinsite_2 = skinsite)) %>% 
  mutate(same_individual = case_when(subject_id_1 == subject_id_2 ~"Within",
                                     subject_id_1 != subject_id_2 ~"Between")) %>% 
  mutate(same_site = case_when(skinsite_1 == skinsite_2 ~"Same_skinsite",
                               skinsite_1 != skinsite_2 ~"Hands vs Forearms")) %>% 
  mutate(same_time = case_when(sample_round_1 == sample_round_2 ~"Same_time",
                               sample_round_1 != sample_round_2 ~"Over time")) %>%
  filter(sample_id_1 != sample_id_2) ##remove dist to the same sample

inter_intra_distances <- bray_within_between_samesite %>% 
  mutate(pair_id = paste(pmin(sample_id_1, sample_id_2), pmax(sample_id_1, sample_id_2), sep="_")) %>% 
  distinct(pair_id, .keep_all = TRUE)

bray_within_between_samesite %>% 
  group_by(same_individual, same_site, subject_id_1, skinsite_1) %>% 
  summarize(mean_dist= mean(dist)) %>% 
  ggplot(aes(x=reorder(same_individual, mean_dist),y=mean_dist, fill=same_individual)) + 
  geom_boxplot() + labs(y="Bray-Curtis distance") + scale_fill_manual(values=c("#9e4f4a", "#ecb775")) +
  facet_wrap(~skinsite_1, labeller= labeller(skinsite_1 = function(variable) { #facet labels are capitalized
    return(tools::toTitleCase(as.character(variable)))})) -> intra_inter_samesite_plot

bray_within_between_samesite %>% 
  filter(skinsite_1 != skinsite_2) %>% 
  group_by(same_individual, same_site, subject_id_1, skinsite_1) %>% 
  summarize(mean_dist= mean(dist)) %>% 
  ggplot(aes(x=reorder(same_individual, mean_dist),y=mean_dist, fill=same_individual)) + 
  geom_boxplot() + facet_wrap(~same_site) + labs(y="Bray-Curtis distance") + scale_fill_manual(values=c("#9e4f4a", "#ecb775")) -> intra_inter_h_vs_f_plot

# Combined plot
intra_inter_beta_plot <- ggarrange(intra_inter_samesite_plot + rremove("xlab") + theme(legend.position = "none", text = element_text(size=16)), 
                                   intra_inter_h_vs_f_plot + rremove("ylab") + rremove("xlab") + 
                                   theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none", text = element_text(size=16)), widths=c(2,1))
ggsave(intra_inter_beta_plot, file="intra_inter_betaplot.pdf", width=11, height=5)


# Beta div intra-individual over time Fig. 5
#load("inter_intra_distances.rds")
inter_intra_distances <- readRDS("PATH/inter_intra_distances.rds")

#between-distance over time, for each skinsite (Fig 5C)
inter_intra_distances %>% as.data.frame() %>% filter(same_site=="Same_skinsite" & same_individual=="Between" & same_time=="Same_time") -> inter_dist
inter_dist %>% ggplot(aes(x=sample_round_1, y=dist, fill=sample_round_1)) + geom_boxplot() + scale_fill_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name="Sample round", labels=c("Baseline", "Post exercise", "3W Post exercise")) + labs(x="", y="Bray-Curtis distance (inter-individual)") + theme_bw() + facet_wrap(~skinsite_1, labeller=labeller(skinsite_1 = c(hands="hands", forearm="forearms"))) -> inter_dist_boxes
saveRDS(inter_dist_boxes, file="inter_dist_boxes.rds")

fig5_a <- inter_intra_distances %>%
  as.data.frame() %>%
  filter(same_site == "Same_skinsite" & same_time == "Over time" & same_individual == "Within") %>%
  group_by(sample_round_1) %>%
  ggplot(aes(x = reorder(skinsite_1, dist), y = dist, fill = skinsite_1)) + 
  geom_boxplot() + 
  scale_fill_discrete(name = "Skin site", labels = c("hands", "forearms")) + 
  labs(x = "", y = "Bray-Curtis distance (intra-individual)") + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  ylim(0.35, 0.97)
                        
fig5_b <- inter_intra_distances %>%
  as.data.frame() %>%
  filter(same_site == "Hands vs Forearms" & same_time == "Same_time" & same_individual == "Within") %>%
  group_by(sample_round_1) %>%
  ggplot(aes(x = sample_round_1, y = dist, fill = sample_round_1)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("azure4", "deepskyblue4", "darkgoldenrod3"), name = "Sample round", labels = c("Baseline", "Post Exercise", "3W Post Exercise")) + 
  labs(x="", y="Bray-Curtis distance (H-to-F)") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  ylim(0.35, 0.97)
          
fig5_c <- readRDS("PATH/inter_dist_boxes.rds") + ylim(0.35, 0.97) + theme(axis.text.x = element_blank())
fig5_d <- readRDS("PATH/beta_disp_boxplot.rds") + ylim(0.35, 0.97) + theme(axis.text.x = element_blank())
fig5_all <- ggarrange(fig5_a + theme(legend.position = "right", legend.text = element_text(size=11), legend.title=element_text(size=12)), 
                     fig5_b + theme(legend.position = "right", legend.text = element_text(size=11), legend.title=element_text(size=12)), 
                     fig5_c + theme(legend.position = "none", strip.text=element_text(size=11)), 
                     fig5_d + theme(legend.position = "none", strip.text=element_text(size=11)), 
                     ncol=2, nrow=2, labels=c("A", "B", "C", "D"))
ggsave(fig5_all, dpi=300, filename = "Figure5.pdf", useDingbats=FALSE, width=7, height=7, units="in")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

sessionInfo()

