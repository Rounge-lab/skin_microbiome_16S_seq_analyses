#Script to create heatmaps for the top differentially abundant taxa with ComplexHeatmap (doi: https://doi.org/10.1002/imt2.43)

# Loading libraries -----------------------------------------------------------------------------------------------------------------------
library(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(vegan)
library(scales)
library(phyloseq)
library(ggplot2)
library(ComplexHeatmap)
library(tibble)
library(grid)
# ------------------------------------------------------------------------------------------------------------------------------------------

#load DESeq2 results
diff_abund_res <- read_rds("diff_abund_res.rds")  

# Function to read and transform OTU data
read_and_transform_otu <- function(file_path) {
  otu_data <- read_rds(file_path) %>%
    rownames_to_column("ASV_id") %>%
    tibble() %>%
    mutate(across(-ASV_id, .fns = function(x) log10(x + 1)))
  return(otu_data)
}

# Function to filter differential abundance results by significance and log2FC
filter_diff_abund_res <- function(diff_abund_res, contrast, sites, padj_threshold = 0.1, abs_fc_threshold = 0.2, top_n = 15) {
  filtered_data <- diff_abund_res %>%
    mutate(abs_fc = abs(log2FoldChange)) %>%
    filter(contrast == contrast, 
           site %in% sites,
           padj < padj_threshold,
           abs_fc >= abs_fc_threshold) %>%
    slice_max(order_by = abs_fc, n = top_n) %>%
    arrange(desc(log2FoldChange)) %>%
    mutate(ASV_id = factor(ASV_id, levels = ASV_id))
  return(filtered_data)
}

# Function to merge differential abundance data with ASV abundance
merge_with_otu <- function(diff_abund_data, otu_data) {
  merged_data <- diff_abund_data %>%
    select(ASV_id) %>%
    left_join(otu_data, by = "ASV_id")
  return(merged_data)
}

# Function to perform Hierarchical clustering on abundance data
perform_clustering <- function(abundance_data) {
  clustering_result <- abundance_data %>%
    column_to_rownames("ASV_id") %>%
    as.matrix() %>%
    t() %>%
    vegdist(method = "bray") %>% #Bray-Curtis dissimilarity
    hclust(method = "ward.D2") # Ward D2 clustering algorithm
  return(clustering_result) #return object of class hclust to use for dendrogram in Heatmap
}

# Function to prepare data matrix for Heatmap()
prepare_heatmap_matrix <- function(abundance_data) {
  heatmap_matrix <- abundance_data %>%
    column_to_rownames("ASV_id") %>%
    as.matrix()
    #heatmap_matrix <- heatmap_matrix[nrow(heatmap_matrix):1, ] #reverse the order of the taxa if needed for plotting
  return(heatmap_matrix)
}


# Loading pre-saved ASV tables
otu_hands <- read_and_transform_otu("/PATH/otu_hands.rds") 
otu_forearms <- read_and_transform_otu("/PATH/otu_forearms.rds") 

# Define parameters for the skin site and contrast you are interested in
skin_sites <- c("h", "f")
contrast_of_interest <- "2v1"

# Initialize lists to store results for each site
filtered_data_list <- list()
merged_data_list <- list()
clustering_result_list <- list()
heatmap_matrix_list <- list()

# Process each site separately
for (site in skin_sites) {
  
  # Determine the correct OTU data to use based on the skin site
  if (site == "h") {
    otu_data <- otu_hands
  } else if (site == "f") {
    otu_data <- otu_forearms
  }

  filtered_data <- filter_diff_abund_res(diff_abund_res, contrast_of_interest, site)
  merged_data <- merge_with_otu(filtered_data, otu_data)
  clustering_result <- perform_clustering(merged_data)
  heatmap_matrix <- prepare_heatmap_matrix(merged_data)
  
  # Store results in lists
  filtered_data_list[[site]] <- filtered_data
  merged_data_list[[site]] <- merged_data
  clustering_result_list[[site]] <- clustering_result
  heatmap_matrix_list[[site]] <- heatmap_matrix
}

#Access matrices to be used for heatmaps below
tmp_abund_matrix_h <- heatmap_matrix_list[["h"]]
tmp_abund_matrix_f <- heatmap_matrix_list[["f"]]




# Function to load rarefied phyloseq object and keep only DA ASVs
load_and_prune <- function(file_path, tmp_abund_matrix, ps_name) {
  load(file_path)
  ps_rarefied <- get(ps_name) #retrieve the phyloseq object using provided name
  pruned_ps <- prune_taxa(rownames(tmp_abund_matrix), ps_rarefied) #keep only taxa significantly diff abund
  return(pruned_ps) 
}

# Load and prune phyloseq objects
ps_DA_h <- load_and_prune("/PATH/hands_paper1_clean_rarefied.RData", tmp_abund_matrix_h, "hands_clean_rare")
ps_DA_f <- load_and_prune("/PATH/forearms_paper1_clean_rarefied.RData", tmp_abund_matrix_f, "f_clean_rare")

# Function to create metadata dataframe (to access sample round in heatmap)
create_metadata <- function(ps_DA) {
  data.frame(sample_data(ps_DA))
}

# Create metadata dataframes
meta_h <- create_metadata(ps_DA_h)
meta_f <- create_metadata(ps_DA_f)


# Function to create annotation column
create_annotation_col <- function(meta_data) {
  annotation_col <- data.frame(`Sample round` = as.factor(meta_data$sample_round), check.names = FALSE)
  annotation_col$`Sample round` <- revalue(annotation_col$`Sample round`, c("1" = "Baseline", "2" = "Post exercise", "3" = "3W Post exercise"))
  rownames(annotation_col) <- rownames(meta_data) #matching sample ids as rownames
  return(annotation_col)
}


# Create annotation columns
annotation_col_h <- create_annotation_col(meta_h)
annotation_col_f <- create_annotation_col(meta_f)

#color metadata
ann_colors = list(
  `Sample round` = c("Baseline" = "azure4", "Post exercise" = "deepskyblue4", "3W Post exercise" = "darkgoldenrod3")
  )

# Function to create heatmap annotations
create_heatmap_annotation <- function(annotation_df, annotation_colors) {
  HeatmapAnnotation(
    df = annotation_df,
    col = annotation_colors,
    simple_anno_size = unit(3.5, "mm"), # height of Sample round boxes
    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"),
    annotation_legend_param = list(
      title = "Sample Round",
      nrow = 1, # horizontal direction of rounds
      title_position = "lefttop",
      labels_gp = gpar(fontsize = 7.5),
      title_gp = gpar(fontsize = 7.5, fontface = "bold")
    ),
    show_annotation_name = FALSE # hide annotation next to dendrogram
  )
}

# Create annotation for each skin site
ha_f <- create_heatmap_annotation(annotation_col_f, ann_colors_f)
ha_h <- create_heatmap_annotation(annotation_col_h, ann_colors_h)


#Color scale bar metrics
create_heatmap_legend_param <- function() {
  list(
    title = "Abundance [log(reads+1)]", #scale bar title and breaks, at=c(0,5,10,15), 
    labels_gp = gpar(fontsize = 7.5),
    title_gp = gpar(fontsize = 7.5, fontface = "bold"), #log2(x+1) size and bold
    direction = "horizontal", 
    title_position = "lefttop",
    grid_height = unit(0.2, "cm"), #height of color scale bar
    legend_width = unit(3, "cm") #width of scale bar
  )
}

# Scale bar for hands and forearms
heatmap_legend_param_h <- create_heatmap_legend_param()
heatmap_legend_param_f <- create_heatmap_legend_param()


# Function to create and save heatmap
create_and_save_heatmap <- function(matrix, legend_param, clustering, ha, filename) {
  # Create heatmap object
  ht <- Heatmap(matrix, 
                col = colorRampPalette(colors = c("gray87", muted("red")))(250), 
                heatmap_legend_param = legend_param,
                cluster_columns = as.dendrogram(clustering), # use the clustering already performed
                cluster_rows = FALSE,
                show_row_dend = TRUE,
                show_column_dend = TRUE,
                column_dend_height = column_dend_height, 
                column_dend_gp = gpar(fontface = "bold"),
                top_annotation = annotation,
                show_column_names = FALSE,
                show_row_names = FALSE)

  # Save heatmap object as grob for saving image
  ht_grob <- grid.grabExpr(draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "top", legend_grouping = "original")) #scale bar below, column annotation on right side 
  ggsave(filename = filename, plot = ht_grob, width = 4, height = 3.2, useDingbats = FALSE, dpi = 300)
  return(ht_grob) #return heatmap grob
}

# Create heatmap for each skin site
ht_grob_h = create_and_save_heatmap(tmp_abund_matrix_h, heatmap_legend_param_h, tmp_hcl_h, ha_h, "heatmap_DA_hands.pdf", column_dend_height = unit(15, "mm")) 
ht_grob_f = create_and_save_heatmap(tmp_abund_matrix_f, heatmap_legend_param_f, tmp_hcl_f, ha_f, "heatmap_DA_forearms.pdf", column_dend_height = unit(13, "mm")) 


# Combine skin site heatmaps into one figure, Fig. 5c-d
heatmap_combined = ggarrange(ht_grob, ht_grob_f, ncol=1, nrow=2, labels=c("C","D"), font.label=list(size=18, color="black", face="bold")) #panel c-d
saveRDS(heatmap_combined, file="heatmap_combined.rds") 

# Making a combined figure 5 by adding DA plot
load("DAplot.rds") #panel a-b created in plot_DA.R
ggarrange(DAplot, heatmap_combined, common.legend = F, align="hv", legend = "right", ncol=2, widths=c(1.8,1)) -> figure5
ggsave(figure5, filename="Fig5_abcd.pdf", useDingbats=FALSE, dpi=300)


sessionInfo()