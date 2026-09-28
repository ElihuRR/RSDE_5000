# This script contains code utilized to generate analyses and plots for the manuscript:
#"Millennia of stability, decades of change: tracing 5,000 years of marine biodiversity shifts in the Red Sea"
# By Rivera Rosas et al., 2026.

# Load necessary packages:
library(RColorBrewer)
library(phyloseq)
library(speedyseq)
library(ggplot2)
library(metagMisc)
library(ape)
library(vegan)
library(picante)
library(cowplot)
library(DESeq2)
library(patchwork)
library(BiocManager)
library(microbiome)
library(devtools)
library(plyr)
library(dplyr)
library(tidyverse)
library(plotly)
library(vegan)
library(zCompositions)
library(phyloseq)
library(data.table)
library(htmlwidgets)
library(iNEXT)
library(microViz)
library(ComplexHeatmap)
library(mvabund)
library(pairwiseAdonis)
library(compositions)
library(Hmisc)
library(mice)
library(ggcorrplot)
library(caret)
library(ALDEx2)
library(fantaxtic)
library(ggnested)
library(phylosmith)
library(igraph)
library(plotly)
library(ggrepel)

## Colors to be utilized throughout:

fancycol<-c("#fddba4","#50d5fe","#f3988a","#b7e899","#f999b3","#7de8c7","#f2b6f6","#b3fec6","#faa73d","#a2b870",
            "#c4befe","#edf55e","#88c7fe","#dfcb81","#24fefd","#fdfbb4","#d460c3","#febb90","#23c6c8","#fea6ab",
            "#5bccaf","#fed3f2","#a2e8ad","#c8a1bf","#f8cd56","#8cb1ca","#c8fe7a","#c9a590","#b9fefb","#fe886b",
            "#fe8cb3","#f4e0c6","#7eb6b6","#9fb47c","#8c7acb","#d8c3f2","#b9fee1","#8ab6a1","#c0e1cf","#28fd59")

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

zone_colors <- c("North"   = "#2ca25f",
                 "Central" = "#8856a7",
                 "South"   = "#e34a33")


## Loading and setting up the objects needed:
phyeu<-readRDS("physeqeukarsdefinal.rds")
phyeu
phyclean<-subset_taxa(phyeu, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean
tax_table(phyclean)
# A total of 923 taxa were eliminated, corresponds to around 50% of reads 
# Of note is that with Euka02 specifically, there are still many Phyla that contain NA even after this cleaning, therefore:
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phyclean2<-subset_taxa(phyclean2, Class!="NA")
phyclean2


phyno<-subset_samples(phyclean, Zone=="North")
phyno<-prune_taxa(taxa_sums(phyno)> 1, phyno)

phyce<-subset_samples(phyclean, Zone=="Central")
phyce<-prune_taxa(taxa_sums(phyce)> 1, phyce)

physu<-subset_samples(phyclean, Zone=="South")
physu<-prune_taxa(taxa_sums(physu)> 1, physu)


##### 1.- Hill-number alpha diversity

## Setting up Hill number diversity for regression

## Fast bootstrapped Hill numbers
## Will directly resample read counts per sample

## Core Hill number functions 
hill_q0 <- function(x) sum(x > 0)                          # Richness
hill_q1 <- function(x) {                                    # exp(Shannon)
  p <- x[x > 0] / sum(x)
  exp(-sum(p * log(p)))
}
hill_q2 <- function(x) {                                    # Inverse Simpson
  p <- x[x > 0] / sum(x)
  1 / sum(p^2)
}

## Bootstrap wrapper for one sample: 
bootstrap_hills <- function(counts, nboot = 500, ci = 0.95) {
  
  counts  <- counts[counts > 0]   # drop absent taxa
  n_reads <- sum(counts)
  probs   <- counts / n_reads
  alpha   <- (1 - ci) / 2
  
  # Resample reads with replacement, recompute Hill numbers
  boot_mat <- replicate(nboot, {
    boot_counts <- as.vector(rmultinom(1, size = n_reads, prob = probs))
    c(q0 = hill_q0(boot_counts),
      q1 = hill_q1(boot_counts),
      q2 = hill_q2(boot_counts))
  })
  
  # Observed values (no resampling)
  obs <- c(q0 = hill_q0(counts),
           q1 = hill_q1(counts),
           q2 = hill_q2(counts))
  
  # CI from bootstrap distribution
  tibble(
    q        = c(0L, 1L, 2L),
    Hill     = obs,
    Hill_LCL = apply(boot_mat, 1, quantile, probs = alpha),
    Hill_UCL = apply(boot_mat, 1, quantile, probs = 1 - alpha)
  )
}

## Apply across all samples in a phyloseq object 

compute_hills_fast <- function(physeq, zone_name, nboot = 500) {
  
  otu <- as.data.frame(otu_table(physeq))
  if (!taxa_are_rows(physeq)) otu <- t(otu)
  
  # Fix: as.data.frame() wrapping BEFORE rownames_to_column and select
  meta <- as.data.frame(as.matrix(sample_data(physeq))) %>%
    rownames_to_column("SampleID") %>%
    dplyr::select(SampleID, New_YearC14) %>%
    mutate(New_YearC14 = as.numeric(New_YearC14))  # coerce – stored as char in some phyloseq builds
  
  message("Bootstrapping Hill numbers for ", zone_name,
          " (", ncol(otu), " samples x ", nboot, " boots)...")
  
  results <- lapply(colnames(otu), function(s) {
    bootstrap_hills(otu[, s], nboot = nboot) %>%
      mutate(SampleID = s)
  }) %>%
    bind_rows()
  
  left_join(results, meta, by = "SampleID") %>%
    mutate(Zone = zone_name)
}

## Run – with nboot = 500
hills_no <- compute_hills_fast(phyno, "North",   nboot = 500)
hills_ce <- compute_hills_fast(phyce, "Central", nboot = 500)
hills_su <- compute_hills_fast(physu, "South",   nboot = 500)

hills_all <- bind_rows(hills_no, hills_ce, hills_su) %>%
  mutate(
    Zone = factor(Zone, levels = c("North", "Central", "South")),
    q_label = factor(q,
                     levels = c(0L, 1L, 2L),
                     labels = c("q = 0  (Richness)",
                                "q = 1  (Shannon-based)",
                                "q = 2  (Simpson-based)"))
  )

## Plot function 
build_hill_plot <- function(q_val, y_label) {
  
  df <- hills_all %>% filter(q == q_val)
  
  ggplot(df, aes(x = New_YearC14, y = Hill,
                 color = Zone, fill = Zone)) +
    
    geom_ribbon(aes(ymin = Hill_LCL, ymax = Hill_UCL),
                alpha = 0.18, color = NA) +
    
    geom_point(size = 2.8, shape = 21,
               stroke = 0.4, alpha = 0.75) +
    
    geom_smooth(method      = "gam",
                formula     = y ~ s(x, bs = "cs"),
                method.args = list(method = "REML"),
                se          = FALSE,
                linewidth   = 1.0) +
    
    scale_color_manual(values = zone_colors) +
    scale_fill_manual(values  = zone_colors) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 10),  # more breaks per panel
      expand = c(0.03, 0.03)                      
    ) +
    
    facet_wrap(~ Zone, ncol = 1, scales = "free_x") +
    
    labs(x = "Year (BCE/CE)", y = y_label) +
    
    theme_bw(base_size = 11) +
    theme( 
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text         = element_text(face = "bold", size = 18),
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey90"),
      axis.text.x        = element_text(angle = 30, hjust = 1, size = 14),
      axis.text.y        = element_text (size=14),
      axis.title         = element_text (size=18)
    )
}

p_D0Eu <- build_hill_plot(0L, "Taxon Richness (q = 0)")
p_D1Eu <- build_hill_plot(1L, "Effective No. Common Taxa (q = 1)")
p_D2Eu <- build_hill_plot(2L, "Effective No. Dominant Taxa (q = 2)")



#### 2.- Relative abundance per zone:

## Fix: average abundance across samples at the same year × zone × phylum ----
## This prevents multiple cores at similar years from stacking additively

phylum_glom <- tax_glom(phyclean2, taxrank = "Phylum")

## Transform to relative abundance PER SAMPLE
## Each sample sums to 1 — this is done per sample, not pooled
phylum_rel <- transform_sample_counts(phylum_glom, function(x) x / sum(x))

## Melt to long format to keep individual samples intact
## Step 1: Extract metadata BEFORE psmelt corrupts it
## Pull directly from sample_data as a clean data frame
meta_clean <- as.data.frame(as.matrix(sample_data(phylum_rel))) %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(Sample, New_YearC14, Zone, Dive_Name) %>%
  mutate(New_YearC14 = as.numeric(New_YearC14),
         Zone        = factor(Zone, levels = c("North", "Central", "South")))


## Step 2: Melt WITHOUT relying on psmelt for metadata
## Extract OTU table and taxonomy separately
otu_mat <- as.data.frame(otu_table(phylum_rel))
if (!taxa_are_rows(phylum_rel)) otu_mat <- t(otu_mat)

tax_mat <- as.data.frame(tax_table(phylum_rel)) %>%
  rownames_to_column("OTU") %>%
  dplyr::select(OTU, Phylum)

## Pivot OTU table to long format manually
phylum_melt <- psmelt(phylum_rel) %>%
  mutate(
    Phylum      = as.character(Phylum),
    New_YearC14 = as.numeric(as.character(New_YearC14)),
    Zone        = factor(as.character(Zone),
                         levels = c("North", "Central", "South"))
  ) %>%
  ## speedyseq renames Sample column to sample_Sample so we standarize it
  dplyr::rename(SampleID = sample_Sample)


## Verify years are correct now
phylum_melt %>%
  group_by(Zone) %>%
  summarise(
    min_yr = min(New_YearC14, na.rm = TRUE),
    max_yr = max(New_YearC14, na.rm = TRUE),
    n_samples = n_distinct(Sample)
  ) %>%
  print()

bar_df <- phylum_melt %>%
  mutate(
    Phylum_grouped = ifelse(Phylum %in% top_phyla,
                            Phylum,
                            "Other (<10% total)"),
    Phylum_grouped = factor(
      Phylum_grouped,
      levels = c(sort(top_phyla), "Other (<10% total)")
    )
  ) %>%
  group_by(SampleID, Phylum_grouped, Zone, New_YearC14, Dive_Name) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Zone, New_YearC14, Phylum_grouped) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Zone = factor(Zone, levels = c("North", "Central", "South")))


phylum_levels <- c(sort(top_phyla), "Other (<10% total)")
phylum_colors <- setNames(fancycol[seq_along(phylum_levels)], phylum_levels)

bar_df_plot <- bar_df %>%
  mutate(
    Year_factor = factor(
      New_YearC14,
      levels = sort(unique(as.numeric(New_YearC14)))
    )
  )

p_barEu <- ggplot(bar_df_plot,
                  aes(x    = Year_factor,
                      y    = Abundance,
                      fill = Phylum_grouped)) +
  
  geom_bar(stat     = "identity",
           position = "stack") +    ## no width needed — discrete bars
  
  scale_fill_manual(values = phylum_colors, name = "Phylum") +
  
  ## Label x-axis with actual years, rotating for readability
  scale_x_discrete(
    labels = function(x) {
      x <- as.numeric(as.character(x))
      ifelse(x < 0,
             paste0(abs(round(x)), " BCE"),
             paste0(round(x), " CE"))
    }
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = 0),
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.01)
  ) +
  
  facet_wrap(~ Zone, ncol = 1, scales = "free_x") +
  
  labs(x = "Year (BCE/CE)", y = "Relative abundance") +
  
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 15),
    axis.text.x      = element_text(angle = 90, hjust = 1,
                                    vjust = 0.5, size = 12),
    axis.text.y      = element_text(size = 9),
    axis.line        = element_line(colour = "black"),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 16),
    legend.text      = element_text(size = 14)
  ) +
  
  guides(fill = guide_legend(nrow = 3))

p_barEu


ggsave("Fig_Barplot_Phylum_18S.png",
       plot   = p_barEu,
       width  = 42,    ## wider to accommodate many discrete bars
       height = 40,
       units  = "cm",
       dpi = 300)
