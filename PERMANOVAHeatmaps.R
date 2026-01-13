# This script contains code utilized to generate analyses and plots for the manuscript:
#"Millennia of stability, decades of change: tracing 5,000 years of marine biodiversity shifts in the Red Sea"
# By Rivera Rosas et al., 2026.

# Load necessary packages:
library(RColorBrewer)
library(decontam); packageVersion("decontam")
library(phyloseq)
library(speedyseq)
library(ggplot2)
library(metagMisc)
library(ape)
library(vegan)
library(picante)
library(cowplot)
library(dplyr)
library(DESeq2)
library(BiocManager)
library(microbiome)
library(plyr)
library(dplyr)
library(tidyverse)
library(plotly)
library(phyloseq)
library(microViz)
library(ggplot2)
library(fastDummies)
library(vegan)
library(zCompositions)
library(phyloseq)
library(data.table)
library(htmlwidgets)
library(iNEXT)
library(biomehorizon)
library(microViz)
library(ComplexHeatmap)
library(compositions)
library(reshape2)
library(MARSS)
library(matrixStats)
library(crqa)
library(nonlinearTseries)
library(adespatial)
library(mvabund)
library(biwavelet)
library(WaveletComp)
library(MetaCycle)
library(zoo)
library(lomb)
library(patchwork)
library(caret)
library(reshape2)


## Figure 4.- Pairwise PERMANOVA Heatmaps of community composition across time periods
# Make sure to set correct directory where .rds files are located

# First, Euka02 (18Sv7 primer)
phyeu<-readRDS("physeqeukarsdefinal.rds")
phyeu
phyclean<-subset_taxa(phyeu, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phyclean2<-subset_taxa(phyclean2, Class!="NA")
phyclean2

# Now running PERMANOVA and pairwise comparisons
meta <- as(sample_data(phyclean2), "data.frame")  # Instead of as.data.frame()

otu <- as.data.frame(otu_table(phyclean2))
otu_clr <- apply(otu + 0.01, 2, clr)
otu_t <- t(otu_clr)

# Ensure rownames match
meta <- meta[rownames(otu_t), , drop = FALSE]
meta$C14_Climate_Name <- as.factor(meta$C14_Climate_Name)

dist_matrix <- dist(otu_t, method = "euclidean")

# Now run PERMANOVA
permanova_result <- adonis2(dist_matrix ~ C14_Climate_Name, data = meta, permutations = 999)
print(permanova_result)

# Now run PERMANOVA for zone
permanova_result <- adonis2(dist_matrix ~ Zone, data = meta, permutations = 999)
print(permanova_result)

# Get all unique time period combinations
time_periods <- unique(meta$C14_Climate_Name)
pairwise_results <- data.frame()

# Loop through all unique pairs of time periods
for (i in 1:(length(time_periods) - 1)) {
  for (j in (i + 1):length(time_periods)) {
    tp1 <- time_periods[i]
    tp2 <- time_periods[j]
    
    # Subset metadata and clr table
    subset_samples <- rownames(meta)[meta$C14_Climate_Name %in% c(tp1, tp2)]
    sub_meta <- meta[subset_samples, , drop = FALSE]
    sub_otu <- otu_t[subset_samples, , drop = FALSE]
    
    # Recalculate distance matrix for subset
    dist_sub <- dist(sub_otu)
    
    # Run PERMANOVA
    result <- adonis2(dist_sub ~ C14_Climate_Name, data = sub_meta, permutations = 999)
    
    # Save results
    pairwise_results <- rbind(pairwise_results, data.frame(
      Group1 = tp1,
      Group2 = tp2,
      F_value = result$F[1],
      R2 = result$R2[1],
      p_value = result$`Pr(>F)`[1]
    ))
  }
}

# Adjust p-values (optional)
pairwise_results$p_adjusted <- p.adjust(pairwise_results$p_value, method = "BH")

# View top results
print(pairwise_results %>% arrange(p_adjusted))

# 1. Create all combinations
all_combos <- expand.grid(Group1 = time_periods, Group2 = time_periods, stringsAsFactors = FALSE)

# 2. Join p-values in both directions (symmetric)
merged_results <- merge(all_combos, pairwise_results[, c("Group1", "Group2", "p_adjusted")], 
                        by = c("Group1", "Group2"), all.x = TRUE)

# 3. Fill in the upper triangle by swapping Group1/Group2
missing_upper <- is.na(merged_results$p_adjusted)
merged_results$p_adjusted[missing_upper] <- pairwise_results$p_adjusted[
  match(paste(merged_results$Group2[missing_upper], merged_results$Group1[missing_upper]),
        paste(pairwise_results$Group1, pairwise_results$Group2))
]


# 4. Final heatmap
v7hm<-ggplot(merged_results, aes(Group1, Group2, fill = p_adjusted)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "#1098ad", 
    high = "#e3fafc", 
    limits = c(0, 1),
    name = "Adj. p-value"
  ) +
  geom_text(
    aes(label = ifelse(p_adjusted < 0.05, sprintf("%.2f", p_adjusted), "")),
    color = "black", size = 4
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  labs(
    title = "18S v7 Primer Pairwise PERMANOVA Heatmap",
    x = "Time Period", y = "Time Period"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title= element_text(size=16),
    axis.text = element_text(size = 14),
    legend.text= element_text(size=12),
    legend.title= element_text(size=14),
    plot.title = element_text(hjust = 0.5, size = 16)
  )
v7hm



## Now for COI
phycoi<-readRDS("physeqco1rsdefinal.rds")
phycoi
phyclean<-subset_taxa(phycoi, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean2<-subset_taxa(phyclean, Class!="NA")

# Now running PERMANOVA and pairwise comparisons
meta <- as(sample_data(phyclean2), "data.frame")  # Instead of as.data.frame()

otu <- as.data.frame(otu_table(phyclean2))
otu_clr <- apply(otu + 0.01, 2, clr)
otu_t <- t(otu_clr)

# Ensure rownames match
meta <- meta[rownames(otu_t), , drop = FALSE]
meta$C14_Climate_Name <- as.factor(meta$C14_Climate_Name)

dist_matrix <- dist(otu_t, method = "euclidean")

# Now run PERMANOVA for zone
permanova_result <- adonis2(dist_matrix ~ Zone, data = meta, permutations = 999)
print(permanova_result)

# Now run PERMANOVA for time
permanova_result <- adonis2(dist_matrix ~ C14_Climate_Name, data = meta, permutations = 999)
print(permanova_result)

# Get all unique time period combinations
time_periods <- unique(meta$C14_Climate_Name)
pairwise_results <- data.frame()

# Loop through all unique pairs of time periods
for (i in 1:(length(time_periods) - 1)) {
  for (j in (i + 1):length(time_periods)) {
    tp1 <- time_periods[i]
    tp2 <- time_periods[j]
    
    # Subset metadata and clr table
    subset_samples <- rownames(meta)[meta$C14_Climate_Name %in% c(tp1, tp2)]
    sub_meta <- meta[subset_samples, , drop = FALSE]
    sub_otu <- otu_t[subset_samples, , drop = FALSE]
    
    # Recalculate distance matrix for subset
    dist_sub <- dist(sub_otu)
    
    # Run PERMANOVA
    result <- adonis2(dist_sub ~ C14_Climate_Name, data = sub_meta, permutations = 999)
    
    # Save results
    pairwise_results <- rbind(pairwise_results, data.frame(
      Group1 = tp1,
      Group2 = tp2,
      F_value = result$F[1],
      R2 = result$R2[1],
      p_value = result$`Pr(>F)`[1]
    ))
  }
}

# Adjust p-values (optional)
pairwise_results$p_adjusted <- p.adjust(pairwise_results$p_value, method = "BH")

# View top results
print(pairwise_results %>% arrange(p_adjusted))

library(reshape2)
library(ggplot2)

# 1. Create all combinations
all_combos <- expand.grid(Group1 = time_periods, Group2 = time_periods, stringsAsFactors = FALSE)

# 2. Join p-values in both directions (symmetric)
merged_results <- merge(all_combos, pairwise_results[, c("Group1", "Group2", "p_adjusted")], 
                        by = c("Group1", "Group2"), all.x = TRUE)

# 3. Fill in the upper triangle by swapping Group1/Group2
missing_upper <- is.na(merged_results$p_adjusted)
merged_results$p_adjusted[missing_upper] <- pairwise_results$p_adjusted[
  match(paste(merged_results$Group2[missing_upper], merged_results$Group1[missing_upper]),
        paste(pairwise_results$Group1, pairwise_results$Group2))
]


# 4. Final heatmap
co1hm<-ggplot(merged_results, aes(Group1, Group2, fill = p_adjusted)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "#0ca678", 
    high = "#c3fae8", 
    limits = c(0, 1),
    name = "Adj. p-value"
  ) +
  geom_text(
    aes(label = ifelse(p_adjusted < 0.05, sprintf("%.2f", p_adjusted), "")),
    color = "black", size = 4
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  labs(
    title = "COI Primer Pairwise PERMANOVA Heatmap",
    x = "Time Period", y = "Time Period"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title= element_text(size=16),
    axis.text = element_text(size = 14),
    legend.text= element_text(size=12),
    legend.title= element_text(size=14),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

co1hm

combined_heat<-v7hm|co1hm
