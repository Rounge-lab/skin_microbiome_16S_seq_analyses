# Loading libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)

# make sure to have a prevalence data set
read_rds(paste0("f", "_prevalence.rds")) %>% 
  as.data.frame() %>% 
  rownames_to_column("ASV_id") %>% 
  tibble()

# Function to read and process differential abundance results
diff_abund_res <- 
  lapply(c("f", "h"), function(site) { # Loop over both skin sites
    lapply(c("2v1", "3v1", "3v2"), function(contrast) { #loop over sample round contrasts
      read_rds(paste0("", site, "_", contrast, "_shrunk.rds")) %>% #read DA result files with shrunk log2FC values
        as.data.frame() %>% 
        rownames_to_column("ASV_id") %>% 
        tibble() %>% 
        mutate(site = site,
               contrast = contrast)
    }) %>% 
      bind_rows() #combine results for all contrasts
  }) %>% 
  bind_rows() #combine results for skin sites
  
# Function to read and process prevalence data
prevalence <- 
  lapply(c("f", "h"), function(site) {
    read_rds(paste0(site, "_prevalence.rds")) %>% 
      as.data.frame() %>% 
      rownames_to_column("ASV_id") %>% 
      tibble() %>% 
      dplyr::rename(prevalence = 2) %>% # Rename the second column to 'prevalence'
      mutate(site = site)
  }) %>% 
  bind_rows() #%>% 
  # pivot_wider(names_from = dataset, values_from = prevalence, values_fill = 0)

#Plotting differential abundance results
make_diff_abund_plot <- function(da_res_tab, 
                                 s, 
                                 select_by_abs = NULL,
                                 min_abs_fc = 0,
                                 contrast_order = c("2v1", "3v2", "3v1"), 
                                 order_by_contrast = "2v1", #it was 3v1
                                 min_abs_in_contrast = TRUE, ## Select based on min_abs in contrast - 
                                                              ## FALSE means feature with any abs_diff > min abs will be kept
                                 sign_in_contrast = TRUE, ## Only select those with significant difference in selected contrast
                                 sign_lvl = 0.1, #FDR p-value threshold
                                 plot_type = "three_panel", #"one_panel", ## Set this to anything else than "one_panel" to split panels by contrast
                                 taxonomy_label = "Genus") {
  
  da_res_tab <-
    da_res_tab %>% 
    mutate(abs_diff = abs(log2FoldChange))
  
  if (!is.null(select_by_abs)) {
    da_res_tab <-
      da_res_tab %>% 
      mutate(log_p = log10(padj)) %>% 
      group_by(site, contrast) %>%
      arrange(desc(abs_diff), desc(log_p)) %>% 
      dplyr::mutate(order_abs_diff = seq(n())) %>% 
      group_by(site) %>% 
      # filter(ASV_id %in% ASV_id[ order_abs_diff %in% (1:select_by_abs) & contrast %in% order_by_contrast & abs_diff > min_abs_fc]) %>%
      filter(ASV_id %in% ASV_id[ order_abs_diff %in% (1:select_by_abs) & contrast %in% order_by_contrast]) %>%
      ungroup()
  }
  
  if (min_abs_in_contrast) {
    da_res_tab <-
      da_res_tab %>% 
      group_by(site) %>% 
      filter(ASV_id %in% ASV_id[ contrast %in% order_by_contrast & abs_diff >= min_abs_fc ]) %>%
      ungroup()  
  } 
  
  if (sign_in_contrast) {
    da_res_tab <-
      da_res_tab %>% 
      group_by(site) %>% 
      filter(ASV_id %in% ASV_id[ contrast %in% order_by_contrast &  padj < sign_lvl]) %>%
      ungroup()  
  }
  
  
  tmp <-
    da_res_tab %>% 
    group_by(ASV_id, site) %>% 
    filter(any(padj < sign_lvl & abs_diff >= min_abs_fc)) %>% 
    ungroup() %>% 
    mutate(ASV_id = factor(ASV_id, 
                           levels = da_res_tab %>% 
                             filter(site %in% s, contrast %in% order_by_contrast) %>% 
                             arrange(log2FoldChange) %>% 
                             pull(ASV_id))) %>% 
    mutate(contrast = factor(contrast, levels = contrast_order)) %>% 
    filter(site %in% s) %>% 
    mutate(sign = ifelse(padj < sign_lvl, "sign", "non-sign")) %>% 
    dplyr::rename(tax_lab = !!taxonomy_label)
  
  # print(tmp[,taxonomy_label])
  # print(tmp$tax_lab)
  # print(taxonomy_label)
  
  
  if (plot_type == "one_panel") {
    tmp %>% 
      ggplot(aes(x = log2FoldChange, 
                 xmin = log2FoldChange-lfcSE, 
                 xmax = log2FoldChange+lfcSE,
                 y = ASV_id,
                 color = contrast, 
                 alpha = sign)) +
      geom_vline(xintercept = 0, linetype = 2, color = "gray") +
      geom_pointrange(position = position_dodge(width = 0.6)) +
      scale_alpha_manual(values = c(0.2,1)) +
      scale_y_discrete(labels = function(hex) enframe(hex, value = "ASV_id") %>% left_join(tmp %>% select(ASV_id, tax_lab) %>% distinct(), by = "ASV_id") %>% pull(tax_lab)) +
      theme_bw() +
      labs(y="ASV") +
      theme(axis.title = element_text(size=16), axis.text=element_text(size=12), 
            legend.text = element_text(size=12), legend.title = element_text(size=14))#+
  } else { #three panels
    tmp %>% mutate(contrast = recode(contrast, "2v1" = "Post exercise vs Baseline", "3v2" = "3W Post exercise vs Post exercise", "3v1" = "3W Post exercise vs Baseline")) %>% 
      ggplot(aes(x = log2FoldChange, 
                 xmin = log2FoldChange-lfcSE, 
                 xmax = log2FoldChange+lfcSE,
                 y = ASV_id,
                 color = Phylum, #color by phylum
                 alpha=sign)) + #significance as bold or not
      geom_vline(xintercept = 0, linetype = 2, color = "gray") +
      geom_pointrange() +
      geom_point() +
      scale_alpha_manual(values = c(0.2,1)) +
      scale_y_discrete(labels = function(hex) enframe(hex, value = "ASV_id") %>% left_join(tmp %>% select(
        ASV_id, tax_lab) %>% distinct(), by = "ASV_id") %>% pull(tax_lab)) +
      theme_bw() +
      scale_color_manual(values = c("p__Actinobacteriota" ="#60a75a", "p__Firmicutes" ="#fd9d53", 
                                    "p__Proteobacteria" = "#4b0082", "p__Bacteroidota" = "#CC79A7"), 
                         labels = c("p__Actinobacteriota" = "Actinomycetota", "p__Bacteroidota" = "Bacteroidota", 
                         "p__Firmicutes" = "Bacillota", "p__Proteobacteria" = "Pseudomonadota")) +
      labs(y="ASV")+
      theme(axis.title = element_text(size=11, face="bold"), axis.text=element_text(size=11), axis.text.y=element_text(face="italic"), 
            legend.text = element_text(size=10), legend.title = element_text(size=11, face="bold")) +
      facet_wrap(~contrast)
  }
  
}

#Generate hand DA plots (Fig. 5a)
(diff_abund_res %>% 
    make_diff_abund_plot(s = "h", 
                         order_by_contrast = "2v1", 
                         select_by_abs = 15, #top 15 by absolute log2FC
                         min_abs_fc = 0.2, 
                         sign_in_contrast = TRUE,
                         min_abs_in_contrast = TRUE) -> h_plot)
h_plot = h_plot + guides(color= guide_legend(order=1), alpha=guide_legend(order=2, title=NULL)) + theme(legend.position = "top")  #remove title of significance label

#Generate forearm DA plots (Fig. 5b)
(diff_abund_res %>% 
    make_diff_abund_plot(s = "f", 
                         order_by_contrast = "2v1", 
                         select_by_abs = NULL, 
                         min_abs_fc = 0.2, 
                         sign_in_contrast = TRUE,
                         min_abs_in_contrast = TRUE) -> f_plot)
f_plot = f_plot + guides(color=guide_legend(order=1), alpha=guide_legend(order=2, title=NULL)) + theme(legend.position = "none") #remove title of significance label

#Combine both hands and forearms DA plots
ggarrange(h_plot, f_plot, common.legend = TRUE, align="hv", labels=c("",""), 
          font.label=list(size=18, color="black", face="bold"), legend = "top", ncol=1) -> DAplot
ggsave(DAplot, filename="Fig5_a_b.pdf", useDingbats=FALSE, width=7.5, height=6.5, dpi=300)
saveRDS(DAplot, file="DAplot.rds") #save to use in heatmaps.R to create full figure 5

#Save tsv file with full differential abundance results (Table S6)
diff_abund_res %>% 
  mutate(site = case_when(site %in% "f" ~ "forearm",
                          site %in% "h" ~ "hand")) %>% 
  write_tsv("differential_abundance_results.tsv")


## Volcano plot
diff_abund_res %>% 
  mutate(sign = ifelse(padj < 0.01, "sign", "non-sign")) %>% 
  left_join(prevalence) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = cut(x = prevalence, breaks = c(0,5,10,Inf)), alpha = sign)) +
  scale_alpha_manual(values = c(0.2,1)) +
  geom_point() +
  facet_wrap(~contrast+site, ncol = 2)

## MA plot
diff_abund_res %>% 
  mutate(sign = ifelse(padj < 0.01, "sign", "non-sign")) %>% 
  left_join(prevalence) %>% 
  ggplot(aes(x = baseMean, y = log2FoldChange, color = cut(x = prevalence, breaks = c(0,5,10,Inf)), alpha = sign)) +
  scale_alpha_manual(values = c(0.2,1)) +
  geom_point() +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~contrast+site, ncol = 2)


sessionInfo()