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
library(indicspecies)
library(tibble)


## Figures S5 and S6.- IndVal Analyses
# Make sure to set correct directory where .rds files are located

# Euka02
setwd("/Users/riverade/Desktop/RSDE eDNA/New Euka02 RSDE") 
phyeu<-readRDS("physeqeukarsdefinal.rds")
phyeu
# For graphics without NAs:
phyclean<-subset_taxa(phyeu, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phyclean2<-subset_taxa(phyclean2, Class!="NA")

# Filtering so rare taxa don't influence the indVal as much. 
prevalence_threshold <- 0.03 * nsamples(phyclean2)
physpe <- prune_taxa(taxa_sums(otu_table(phyclean2)) >= 50 & 
                       rowSums(otu_table(phyclean2) > 0) >= prevalence_threshold, phyclean2)
physpe

# Get abundance table
otu_mat <- as.data.frame(otu_table(physpe))
if (taxa_are_rows(physpe)) {
  otu_mat <- t(otu_mat)
}

# Get sample groupings 
groups <- sample_data(physpe)$C14_Climate_Name

# Transform, trying with relative abundance and presence absence
otu_rel <- sweep(otu_mat, 1, rowSums(otu_mat), FUN = "/")
otu_pa <- ifelse(otu_mat > 0, 1, 0)

# Running indval
indval_resulteusp <- multipatt(otu_pa, groups, func = "IndVal.g", control = how(nperm = 999))

# To view summary
indval_resulteusp
summary(indval_resulteusp)
saveRDS(indval_resulteusp, "Euka02SpeciesIndValResult.rds")

# Assigning indicator taxa for Presence/Absence according to the summary of results

# To know what each "spXXXX" is:
otu_to_taxa <- data.frame(
  OTU = taxa_names(physpe),
  Kingdom = tax_table(physpe)[, "Kingdom"],
  Phylum = tax_table(physpe)[, "Phylum"],
  Class = tax_table(physpe)[, "Class"],
  Order = tax_table(physpe)[, "Order"],
  Family = tax_table(physpe)[, "Family"],
  Genus = tax_table(physpe)[, "Genus"],
  Species = tax_table(physpe)[, "Species"]
)
# Filter to only these and view their Class
indicator_taxa_eusp <- otu_to_taxa %>%
  filter(OTU %in% indicator_taxaeusp)
indicator_taxa_eusp

# Melting the phyloseq object
phy_melted <- psmelt(physpe)

# Subset to only indicator taxa
phy_melted_indeusp <- phy_melted[phy_melted$OTU %in% indicator_taxaeusp, ]

time_order <- c("Late Holocene Arid Transition", 
                "Bronze Age Recovery", 
                "4.2 ka Drought", 
                "Bronze Age Warm Peak", 
                "3.2 ka Late Bronze Age Drought", 
                "Iron Age Cold Epoch", 
                "Roman Warm Period", 
                "Late Antique Little Ice Age", 
                "Medieval Warm Period", 
                "Little Ice Age", 
                "Early Industrial Revolution", 
                "Recent Anthropogenic Warming")

phy_melted_indeusp$C14_Climate_Name <- factor(phy_melted_indeusp$C14_Climate_Name, levels = time_order)
phy_melted_indeusp$FacetLabel <- phy_melted_indeusp$Class

# To get the full table of Indval
# Get stats and p-values
indval_df <- as.data.frame(indval_resulteusp$sign) %>%
  rownames_to_column("OTU")

indval_filtered <- indval_df %>%
  filter(stat > 0.7, !is.na(p.value), p.value <= 0.05)

# Getting the taxonomic table from the phyloseq file and joining it to this new IndVal table
# To know what each "spXXXX" is:
otu_to_taxa <- as.data.frame(as.matrix(tax_table(physpe)))
otu_to_taxa$OTU <- rownames(otu_to_taxa)

indval_annotated <- indval_filtered %>%
  dplyr::left_join(otu_to_taxa, by = "OTU")

# Pivot to long format
indval_long <- indval_annotated %>%
  pivot_longer(cols = starts_with("s."), names_to = "TimePeriod", values_to = "Indicator") %>%
  filter(Indicator == 1)

# Fix time period order here 
time_order <- c("Late Holocene Arid Transition", 
                "Bronze Age Recovery", 
                "4.2 ka Drought", 
                "Bronze Age Warm Peak", 
                "3.2 ka Late Bronze Age Drought", 
                "Iron Age Cold Epoch", 
                "Roman Warm Period", 
                "Late Antique Little Ice Age", 
                "Medieval Warm Period", 
                "Little Ice Age", 
                "Early Industrial Revolution", 
                "Recent Anthropogenic Warming")

# Make sure TimePeriod is treated as factor with the right order
indval_long$TimePeriod <- gsub("^s\\.", "", indval_long$TimePeriod)
indval_long$TimePeriod <- factor(indval_long$TimePeriod, levels = time_order)

# Plot
euspeplot<-ggplot(indval_long, aes(x = TimePeriod, y = OTU)) +
  geom_point(aes(color = Phylum), size = 4) +
  scale_color_d3(palette = "category20") +
  theme_minimal() +
  theme(axis.title= element_text(size=14)) +
  theme(axis.text.y= element_text(size=12)) +
  theme(legend.title=element_text(size=16), legend.text=element_text(size=14)) +
  theme(axis.text.x = element_text(size= 12, angle = 45, hjust = 1)) +
  labs(title = "18S v7 Indicator OTUs Across Time Periods", y = "OTUs", x = "Time Period")

euspeplot



## Now for COI 
phyco<-readRDS("physeqco1rsdefinal.rds")
phyco
# For graphics without NAs:
phyclean<-subset_taxa(phyco, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phyclean2<-subset_taxa(phyclean, Class!="NA")
phyclean2

# Filtering so rare taxa don't influence the indVal as much. 
prevalence_threshold <- 0.03 * nsamples(phyclean2)
physpe<- prune_taxa(taxa_sums(otu_table(phyclean2)) >= 50 & 
                      rowSums(otu_table(phyclean2) > 0) >= prevalence_threshold, phyclean2)
physpe <- subset_taxa(physpe, Class != "Collembola")
physpe

# Get abundance table
otu_mat <- as.data.frame(otu_table(physpe))
if (taxa_are_rows(physpe)) {
  otu_mat <- t(otu_mat)
}

# Get sample groupings 
groups <- sample_data(physpe)$C14_Climate_Name

# Transform, trying with relative abundance and presence absence
otu_rel <- sweep(otu_mat, 1, rowSums(otu_mat), FUN = "/")
otu_pa <- ifelse(otu_mat > 0, 1, 0)

# Running indval
indval_resultcoisp <- multipatt(otu_pa, groups, func = "IndVal.g", control = how(nperm = 999))

# To view summary
indval_resultcoisp
summary(indval_resultcoisp)
saveRDS(indval_resultcoisp, "COISpeciesIndValResult.rds")

# To know what each "spXXXX" is:
otu_to_taxa <- data.frame(
  OTU = taxa_names(phy_coisp),
  Kingdom = tax_table(phy_coisp)[, "Kingdom"],
  Phylum = tax_table(phy_coisp)[, "Phylum"],
  Class = tax_table(phy_coip)[, "Class"],
  Order = tax_table(phy_coisp)[, "Order"],
  Family = tax_table(phy_coisp)[, "Family"],
  Genus = tax_table(phy_coisp)[, "Genus"],
  Species = tax_table(phy_coisp)[, "Species"]
)
# Filter to only these and view their Class
indicator_taxa_coisp <- otu_to_taxa %>%
  filter(OTU %in% indicator_taxacoisp)
indicator_taxa_coisp

# Melting the phyloseq object
phy_melted <- psmelt(physpe)

# Subset to only indicator taxa
phy_melted_indcoisp <- phy_melted[phy_melted$OTU %in% indicator_taxacoisp, ]

time_order <- c("Late Holocene Arid Transition", 
                "Bronze Age Recovery", 
                "4.2 ka Drought", 
                "Bronze Age Warm Peak", 
                "3.2 ka Late Bronze Age Drought", 
                "Iron Age Cold Epoch", 
                "Roman Warm Period", 
                "Late Antique Little Ice Age", 
                "Medieval Warm Period", 
                "Little Ice Age", 
                "Early Industrial Revolution", 
                "Recent Anthropogenic Warming")

phy_melted_indcoisp$C14_Climate_Name <- factor(phy_melted_indcoisp$C14_Climate_Name, levels = time_order)
phy_melted_indcoisp$FacetLabel <- phy_melted_indcoisp$Class


# To get the full table of Indval
# Get stats and p-values
indval_df <- as.data.frame(indval_resultcoisp$sign) %>%
  rownames_to_column("OTU")

indval_filtered <- indval_df %>%
  filter(stat > 0.7, !is.na(p.value), p.value <= 0.05)

# Getting the taxonomic table from the phyloseq file and joining it to this new IndVal table
# To know what each "spXXXX" is:
otu_to_taxa <- as.data.frame(as.matrix(tax_table(physpe)))
otu_to_taxa$OTU <- rownames(otu_to_taxa)

indval_annotated <- indval_filtered %>%
  dplyr::left_join(otu_to_taxa, by = "OTU")

write.csv(indval_annotated, "IndValResultsCOISp.csv")

# Pivot to long format
indval_long <- indval_annotated %>%
  pivot_longer(cols = starts_with("s."), names_to = "TimePeriod", values_to = "Indicator") %>%
  filter(Indicator == 1)

# Fix time period order here 
time_order <- c("Late Holocene Arid Transition", 
                "Bronze Age Recovery", 
                "4.2 ka Drought", 
                "Bronze Age Warm Peak", 
                "3.2 ka Late Bronze Age Drought", 
                "Iron Age Cold Epoch", 
                "Roman Warm Period", 
                "Late Antique Little Ice Age", 
                "Medieval Warm Period", 
                "Little Ice Age", 
                "Early Industrial Revolution", 
                "Recent Anthropogenic Warming")

# Make sure TimePeriod is treated as factor with the right order
indval_long$TimePeriod <- gsub("^s\\.", "", indval_long$TimePeriod)
indval_long$TimePeriod <- factor(indval_long$TimePeriod, levels = time_order)


# Plot
coispeplot<-ggplot(indval_long, aes(x = TimePeriod, y = OTU)) +
  geom_point(aes(color = Phylum), size = 4) +
  scale_color_d3(palette = "category20") +
  theme_minimal() +
  theme(axis.title= element_text(size=14)) +
  theme(axis.text.y= element_text(size=12)) +
  theme(legend.title=element_text(size=16), legend.text=element_text(size=14)) +
  theme(axis.text.x = element_text(size= 12, angle = 45, hjust = 1)) +
  labs(title = "COI Indicator OTUs Across Time Periods", y = "OTUs", x = "Time Period")


coispeplot 