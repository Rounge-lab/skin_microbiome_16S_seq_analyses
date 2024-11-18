#!/usr/bin/Rscript

#install.packages("tidyr")
#install.packages("stringr")
library(tidyr)
library(stringr)

#import taxonomy of ASVs for all 217 samples
taxonomy = read.table("/PATH/taxonomy_all_above_5000_reads_f270_r210.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
taxonomy = as.data.frame(taxonomy)
tax_levels=strsplit(taxonomy$Taxon, "; ") #split into rank columns
tax_matrix = do.call(rbind, lapply(tax_levels, function(x) {
  levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_levels = c(x, rep(NA, length(levels) - length(x)))
  setNames(tax_levels, levels)
}))
tax.clean = as.data.frame(tax_matrix)
rownames(tax.clean) = taxonomy$`Feature ID` #add ASV IDs as row names
tax.clean[is.na(tax.clean)] = "" #remove NAs
tax.clean[tax.clean=="__"] <- ""

#Remove prefix to taxonomic levels (optional)
# tax.clean <- data.frame(row.names = taxonomy$`Feature ID`, #add ASV IDs as row names)
#                         Kingdom = str_replace(tax_matrix[,1], "k__",""),
#                         Phylum = str_replace(tax_matrix[,2], "p__",""),
#                         Class = str_replace(tax_matrix[,3], "c__",""),
#                         Order = str_replace(tax_matrix[,4], "o__",""),
#                         Family = str_replace(tax_matrix[,5], "f__",""),
#                         Genus = str_replace(tax_matrix[,6], "g__",""),
#                         Species = str_replace(tax_matrix[,7], "s__",""),
#                         stringsAsFactors = FALSE)

# Loop through rows and fill missing classifications with "Unclassified" labels based on the last known classification level
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){  # If the Species is classified (not empty), keep it as is
    tax.clean$Species[i] <- tax.clean$Species[i]
    # For each else if statement, check if a taxonomic classification (Phylum, Class, Order, Family, Genus) is empty.
    # If empty, create an "Unclassified" label using the previous higher taxonomic level and assign it to all lower levels.
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = "_")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = "_")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = "_")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = "_")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = "_")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified",tax.clean$Genus[i], sep = "_")
  }
}

save(tax.clean, file="taxonomy_all_above_5000_f270_r210.RData")
write.table(tax.clean, "taxonomy_clean.tsv", sep="\t", quote=FALSE, row.names=FALSE)

cat("Done with the script\n")


