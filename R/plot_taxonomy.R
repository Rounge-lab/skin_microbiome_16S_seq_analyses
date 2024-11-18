#Loading libraries
library(dplyr)
library(phyloseq) 
library(ggplot2)
library(ggpubr)

#load non-normalized, clean phyloseq object
load("/PATH/ps_sol_paper1_clean.RData") 

#Baseline, relative abundances
ps_baseline <- ps_soldiers_paper1_clean %>%
  prune_samples(sample_data(.)$sample_round == 1, .) %>% #baseline = sample round 1
  prune_taxa(taxa_sums(.) != 0, .) %>% #keep only taxa with non-zero sums
  transform_sample_counts(function(x) x / sum(x) * 100) #relative abundances

#phyloseq object: hands only
ps_hands = <- ps_baseline %>%
  subset_samples(skinsite == "hands") %>%
  prune_taxa(taxa_sums(.) != 0, .)
#phyloseq object: forearms only
ps_forearms <- ps_baseline %>%
  subset_samples(skinsite == "forearm") %>%
  prune_taxa(taxa_sums(.) != 0, .) %>%
  { 
    bad_f <- sort(sample_sums(.)) %>% head(15) %>% names() #the 15 "bad" forearm samples with lowest seq depth
    prune_samples(!(sample_data(.)$sampleid %in% bad_f), .) #keep only the 30 best samples
  }

# Function to glom by taxonomic rank
glom_data <- function(physeq_obj, taxrank) {
  tax_glom(physeq_obj, taxrank = taxrank, NArm = FALSE) %>%
    psmelt()
}

# Function to order taxonomy, show those with certain abundance
process_data <- function(melted_data, min_abundance) {
  melted_data %>%
    mutate(Genus = as.character(Genus)) %>%
    group_by(sample_round, Genus) %>%
    mutate(median = median(Abundance)) %>%
    ungroup() %>%
    mutate(Genus = ifelse(median <= min_abundance, "Other", Genus)) %>%
    mutate(subjectid = as.factor(subjectid)) %>%
    group_by(subjectid, sample_round, Genus) %>%
    summarise(Abundance = sum(Abundance), .groups = 'drop')
}

# Function to plot the processed data
plot_data <- function(summarized_data) {
  theme_set(theme_bw())
  ggplot(summarized_data, aes(x = subjectid, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", aes(fill = Genus)) +
    labs(x = NULL, y = "Relative abundance [%]") +
    theme(axis.text.x = element_blank(),  # Remove x-axis text (subject IDs)
          axis.ticks.x = element_blank()) + # Remove x-axis ticks
    facet_wrap(~sample_round, scales = "free_x", nrow = 3)
}

# Function to handle multiple Phyloseq objects
process_phyloseq_objects <- function(physeq_objs, taxrank, min_abundance) {
  plots <- lapply(physeq_objs, function(obj) {
    melted_data <- glom_data(obj, taxrank)
    summarized_data <- process_data(melted_data, min_abundance)
    plot_data(summarized_data)
  })
  
  return(plots)
}

# Run functions on chosen phyloseq objects
physeq_objs <- list(ps_hands, ps_forearms) #insert phyloseq objects
plots <- process_phyloseq_objects(physeq_objs, taxrank = "Genus", min_abundance = 2.5) #choose taxonomic rank and minimum abundance to plot

# Combine plots into one figure with labels A, B, C, etc and individual legends
n_plots <- length(plots)
combined_plot <- ggarrange(plotlist = plots, labels = LETTERS[1:n_plots], ncol = ceiling(sqrt(n_plots)), nrow = ceiling(n_plots / ncol), common.legend = FALSE, align="hv")

# Print the combined taxonomic plot
print(combined_plot)



sessionInfo()
