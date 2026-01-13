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
library(mia)
library(mice)
library(MetBrewer)


## Figure 7 and S4.- RDAs by time periods and by zone
# Make sure to set correct directory where .rds files are located

# First, Euka02 (18Sv7 primer)

phyeu<-readRDS("physeqeukarsdefinal.rds")
phyeu

# Remove taxa with unreliable taxonomy, as well as those that contain NA in Phylum
phyclean<-subset_taxa(phyeu, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phyclean2<-subset_taxa(phyclean2, Class!="NA")
phyclean2

# Adding the sample data that will be used (this one has been trimmed to only have the relevant numeric variables
metadata<-as.data.frame(read.csv("RSDEOnlyNumericMetadataSubsetEuka02.csv"))

# Imputing data to fill the NAs
sample_names<-metadata$SampleName
metadatanum<-metadata[,-1]
imputed_metadata <- mice(metadatanum, method = "pmm", m = 1, maxit = 5)
metadata_clean <- complete(imputed_metadata)
metadata_clean
metadata_clean <- cbind(SampleName = sample_names, metadata_clean)

meta_dummies <- dummy_cols(
  metadata_clean,
  select_columns = c("Zone", "C14_Climate_Name"),
  remove_first_dummy = TRUE,
)
rownames(meta_dummies) <- meta_dummies$SampleName
meta_dummies$SampleName <- NULL

# Fix characters in the column names, otherwise the stats methods give errors
colnames(meta_dummies) <- gsub(" ", "_", colnames(meta_dummies))
colnames(meta_dummies) <- gsub("'", "", colnames(meta_dummies))
colnames(meta_dummies) <- gsub("-", "_", colnames(meta_dummies))

# Overwriting the metadata in phyloseq object
sample_data(phyclean2) <- meta_dummies

# Extract OTU table from phyloseq. Important: transform as you see fit (I aggregate on Phylum). Needed as part of CLR transform.
phylum_fixed <- phyclean2 %>%
  tax_fix() %>%
  tax_glom("Phylum", NArm = FALSE) # Keeps all even if upstream taxonomy differs
comm <- tax_transform(phylum_fixed, trans = "clr") %>%
  otu_table() %>%
  t()

# Do RDA with all variables to then check which to use
rda <- rda(comm ~ Zone+C14_Climate_Name+Longitude+Latitude+New_YearC14+Depth+New_YearC14+RelSeaLvl+TN_Percentage+TOC_Percentage+CN+
             TC_Percentage+CaCO3_Percentage+Opal_Percentage+As+Cd+Cr+Cu+Mn+Mo+Ni+P+NP+Pb+Sb+Se+Sn+Sr+V+Zn+Al_Percentage+
             Ca_Percentage+Fe_Percentage+K_Percentage+Mg_Percentage+Na_Percentage+S_Percentage+d0.5μm+Si_Percentage+Cl_Percentage+
             Ti_Percentage, data = meta_dummies)

# 1. Check collinearity using VIF (Variance Inflation Factor):
vif.cca(rda)
#A VIF > 10 typically suggests problematic multicollinearity.
# 2. Do a global permutation test for the full model. Only continue if significant
anova(rda, permutations = 999)
# 3. Stepwise selection of significant variables
null_model <- rda(comm ~ 1, data = meta_dummies)
step_model <- ordistep(null_model, scope = formula(rda), direction = "both", permutations = 999)
# 4. Check significance of individual variables:
anov<-anova(step_model, by = "term", permutations = 999)


#Save the step model object as an RDS object to avoid having to do it again (lots of processing time)
saveRDS(step_model, "rda_step_modeleuka02.rds")
step_model<-readRDS("rda_step_modeleuka02.rds")

# Now with this, we can proceed with the ellaboration of graphics with only the relevant variables
anov
# Based on this, the significant variables will be:
# Depth, Zone, K_Percentage, Ni, TOC_Percentage, Al_Percentage, Mn, C14_Climate_Name?, SR, NP, V, 
significant_vars <- c(
  "Depth", "K_Percentage", "Ni", 
  "Na_Percentage", "Sr",
  "P"
)

# Zone and C14_Climate_Name will be used not for coloring the samples in the RDA
# Transform and calculate constrained ordination (RDA)
# Convert grouping variables to factors
sample_data(phyclean2)$Zone <- factor(sample_data(phyclean2)$Zone)
sample_data(phyclean2)$C14_Climate_Name <- factor(sample_data(phyclean2)$C14_Climate_Name)


# Run ordination calculation (CLR will be done automatically)
# Fix convergent taxonomy issues manually by specifying unknowns
phyclean2_fixed <- phyclean2 %>%
  tax_fix(unknowns = c(
    "Amoebozoa", "Annelida", "Apicomplexa", "Arthropoda", "Ascomycota", "Basidiomycota",
    "Bigyra", "Cercozoa", "Charophyta", "Chlorophyta", "Choanozoa", "Chordata",
    "Chytridiomycota", "Ciliophora", "Cnidaria", "Cryptophyta", "Dinoflagellata", "Gyrista",
    "Haptophyta", "Hyphochytriomycota", "Mollusca", "Myzozoa", "Nematoda", "Nibbleridea",
    "Platyhelminthes", "Porifera", "Prasinodermophyta", "Radiozoa", "Retaria", "Rhodophyta",
    "Sulcozoa", "Zygomycota"
  ))

# Now aggregate at Phylum
phylum_fixed <- phyclean2_fixed %>%
  tax_glom("Phylum", NArm = FALSE)

sample_data(phylum_fixed) <- meta_dummies
colnames(sample_data(phylum_fixed))
sample_data(phylum_fixed)$Zone <- factor(
  sample_data(phylum_fixed)$Zone,
  levels = zone_order
)

rda_zoneeuk <- phylum_fixed %>%
  tax_transform(trans = "clr") %>%
  ord_calc(
    constraints = significant_vars,
    method="RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Zone", fill = "Zone",
    shape = "Zone", alpha = 0.5,
    size = 3,
    constraint_vec_style = vec_constraint(1.5, alpha = 0.8), 
    constraint_lab_style = constraint_lab_style(size = 5, alpha=0.5)
  )  + 
  scale_color_manual(values = newzonecol) +
  scale_fill_manual(values = newzonecol)+
  scale_shape_girafe_filled() +
  ggside::geom_xsideboxplot(aes(fill = Zone, y = Zone), orientation = "y", varwidth=FALSE, show.legend=FALSE) +
  ggside::geom_ysideboxplot(aes(fill = Zone, x = Zone), orientation = "x", varwidth=FALSE, show.legend=FALSE) +
  ggside::theme_ggside_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size=18),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.title = element_text(face = "bold", size=16),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 14)
  ) +
  labs(
    x = "RDA1 [6.4%]", 
    y = "RDA2 [3.5%]", 
    title = "Redundancy Analysis (RDA) of Eukaryotic Communities according to the 18S v7 dataset",
    caption = NULL
  )

rda_zoneeuk

# Now running it for time period.

# Need to add many different shapes, thus:
# Generate distinct shapes (R allows up to 25 default shapes: 0-25)
n_climates <- length(unique(sample_data(phylum_fixed)$C14_Climate_Name))
shapes <- rep(0:25, length.out = n_climates)
names(c25) <- levels(factor(sample_data(phylum_fixed)$C14_Climate_Name))
names(shapes) <- levels(factor(sample_data(phylum_fixed)$C14_Climate_Name))
sample_data(phylum_fixed)$C14_Climate_Name <- factor(
  sample_data(phylum_fixed)$C14_Climate_Name,
  levels = time_order
)

rda_timeeuk <- phylum_fixed %>%
  tax_transform(trans = "clr") %>%
  ord_calc(
    constraints = significant_vars,
    method = "RDA"
  ) %>%
  ord_plot(
    axes = c(1, 2),
    colour = "C14_Climate_Name", fill = "C14_Climate_Name",
    shape = "C14_Climate_Name", alpha = 0.5,
    size = 3,
    constraint_vec_style = vec_constraint(1.5, alpha = 0.8),
    constraint_lab_style = constraint_lab_style(size = 5, alpha = 0.5)
  ) +
  scale_color_manual(values = c25) +
  scale_fill_manual(values = c25) +
  scale_shape_manual(values = shapes) +
  ggside::geom_xsideboxplot(aes(fill = C14_Climate_Name, y = C14_Climate_Name), orientation = "y", varwidth = FALSE, show.legend = FALSE) +
  ggside::geom_ysideboxplot(aes(fill = C14_Climate_Name, x = C14_Climate_Name), orientation = "x", varwidth = FALSE, show.legend = FALSE) +
  ggside::theme_ggside_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.title = element_text(face = "bold", size=16),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 14)
  ) +
  labs(
    x = "RDA1 [6.4%]",
    y = "RDA2 [3.5%]",
    title = "Redundancy Analysis (RDA) of Eukaryotic Communities by C14 Climate Periods",
    caption = "168 samples & 1804 taxa (Phylum); tax_transform = clr"
  )

# Display
rda_timeeuk



## Now for COI
phyco<-readRDS("physeqco1rsdefinal.rds")
phyco
# For graphics without NAs:
phyclean<-subset_taxa(phyco, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phyclean2<-subset_taxa(phyclean, Class!="NA")
phyclean2

# Adding the sample data that will be used (this one has been trimmed to only have the relevant variables)
metadata<-as.data.frame(read.csv("RSDEWorkingOnlyNumericCO1SamData.csv"))

# Imputing data to fill the NAs
sample_names<-metadata$SampleName
metadatanum<-metadata[,-1]
imputed_metadata <- mice(metadatanum, method = "pmm", m = 1, maxit = 5)
metadata_clean <- complete(imputed_metadata)
metadata_clean
metadata_clean <- cbind(SampleName = sample_names, metadata_clean)

meta_dummies <- dummy_cols(
  metadata_clean,
  select_columns = c("Zone", "C14_Climate_Name"),
  remove_first_dummy = TRUE,
)
rownames(meta_dummies) <- meta_dummies$SampleName
meta_dummies$SampleName <- NULL

# I have to fix characters in the column names, otehrwise the stats methods give errors
colnames(meta_dummies) <- gsub(" ", "_", colnames(meta_dummies))
colnames(meta_dummies) <- gsub("'", "", colnames(meta_dummies))
colnames(meta_dummies) <- gsub("-", "_", colnames(meta_dummies))

# Overwriting the metadata in phyloseq object
sample_data(phyclean2) <- meta_dummies

# Extract OTU table from phyloseq. Important: transform as you see fit (I aggregate on Phylum). Needed as part of CLR transform.
phylum_fixed <- phyclean2 %>%
  tax_fix() %>%
  tax_glom("Phylum", NArm = FALSE) # Keeps all even if upstream taxonomy differs
comm <- tax_transform(phylum_fixed, trans = "clr") %>%
  otu_table() %>%
  t()

# Do RDA with all variables to then check which to use
rda <- rda(comm ~ Zone+C14_Climate_Name+Longitude+Latitude+New_YearC14+Depth+New_YearC14+RelSeaLvl+TN_Percentage+TOC_Percentage+CN+
             TC_Percentage+CaCO3_Percentage+Opal_Percentage+As+Cd+Cr+Cu+Mn+Mo+Ni+P+NP+Pb+Sb+Se+Sn+Sr+V+Zn+Al_Percentage+
             Ca_Percentage+Fe_Percentage+K_Percentage+Mg_Percentage+Na_Percentage+S_Percentage+Si_Percentage+Cl_Percentage+
             Ti_Percentage, data = meta_dummies)

# 1. Check collinearity using VIF (Variance Inflation Factor):
vif.cca(rda)
#A VIF > 10 typically suggests problematic multicollinearity.
# 2. Do a global permutation test for the full model. Only continue if significant
anova(rda, permutations = 999)
# 3. Stepwise selection of significant variables
null_model <- rda(comm ~ 1, data = meta_dummies)
step_model <- ordistep(null_model, scope = formula(rda), direction = "both", permutations = 999)
# 4. Check significance of individual variables:
anov<-anova(step_model, by = "term", permutations = 999)

#Save the step model object as an RDS object to avoid having to do it again (lots of processing time)
saveRDS(step_model, "rda_step_modelCOI.rds")
step_model<-readRDS("rda_step_modelCOI.rds")

# Now with this, we can proceed with the ellaboration of graphics with only the relevant variables
anov
# Based on this, the significant variables will be:
# 
significant_vars <- c("Depth", "CN", "P", "Opal_Percentage", 
                      "Ni", "Mg_Percentage")

# Zone and C14_Climate_Name will be used not for coloring the samples in the RDA
# Transform and calculate constrained ordination (RDA)
# Convert grouping variables to factors
sample_data(phyclean2)$Zone <- factor(sample_data(phyclean2)$Zone)
sample_data(phyclean2)$C14_Climate_Name <- factor(sample_data(phyclean2)$C14_Climate_Name)


# Run ordination calculation (CLR will be done automatically)
# Fix convergent taxonomy issues manually by specifying unknowns
phyclean2_fixed <- phyclean2 %>%
  tax_fix(unknowns = c(
    "Amoebozoa", "Annelida", "Apicomplexa", "Arthropoda", "Ascomycota", "Basidiomycota",
    "Bigyra", "Cercozoa", "Charophyta", "Chlorophyta", "Choanozoa", "Chordata",
    "Chytridiomycota", "Ciliophora", "Cnidaria", "Cryptophyta", "Dinoflagellata", "Gyrista",
    "Haptophyta", "Hyphochytriomycota", "Mollusca", "Myzozoa", "Nematoda", "Nibbleridea",
    "Platyhelminthes", "Porifera", "Prasinodermophyta", "Radiozoa", "Retaria", "Rhodophyta",
    "Sulcozoa", "Zygomycota"
  ))

# Now aggregate at Phylum
phylum_fixed <- phyclean2_fixed %>%
  tax_glom("Phylum", NArm = FALSE)

sample_data(phylum_fixed) <- meta_dummies
colnames(sample_data(phylum_fixed))
sample_data(phylum_fixed)$Zone <- factor(
  sample_data(phylum_fixed)$Zone,
  levels = zone_order
)

rda_zonecoi <- phylum_fixed %>%
  tax_transform(trans = "clr") %>%
  ord_calc(
    constraints = significant_vars,
    method="RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Zone", fill = "Zone",
    shape = "Zone", alpha = 0.5,
    size = 2,
    constraint_vec_style = vec_constraint(1.5, alpha = 0.8), 
    constraint_lab_style = constraint_lab_style(size = 5, alpha=0.5)
  )  + 
  scale_color_manual(values = newzonecol) +
  scale_fill_manual(values = newzonecol)+
  scale_shape_girafe_filled() +
  ggside::geom_xsideboxplot(aes(fill = Zone, y = Zone), orientation = "y", varwidth=FALSE, show.legend=FALSE) +
  ggside::geom_ysideboxplot(aes(fill = Zone, x = Zone), orientation = "x", varwidth=FALSE, show.legend=FALSE) +
  ggside::theme_ggside_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size=16),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.title = element_text(face = "bold", size=16),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 14)
  ) +
  labs(
    x = "RDA1 [10.2%]",
    y = "RDA2 [3.5%]",
    title = "Redundancy Analysis (RDA) of Eukaryotic Communities by C14 Climate Periods",
    caption = NULL
  )

rda_zonecoi


# Now running it for time period.

# Need to add many different shapes, thus:
# Generate distinct shapes (R allows up to 25 default shapes: 0-25)
n_climates <- length(unique(sample_data(phylum_fixed)$C14_Climate_Name))
shapes <- rep(0:25, length.out = n_climates)
names(c25) <- levels(factor(sample_data(phylum_fixed)$C14_Climate_Name))
names(shapes) <- levels(factor(sample_data(phylum_fixed)$C14_Climate_Name))
sample_data(phylum_fixed)$C14_Climate_Name <- factor(
  sample_data(phylum_fixed)$C14_Climate_Name,
  levels = time_order
)

rda_timecoi <- phylum_fixed %>%
  tax_transform(trans = "clr") %>%
  ord_calc(
    constraints = significant_vars,
    method = "RDA"
  ) %>%
  ord_plot(
    axes = c(1, 2),
    colour = "C14_Climate_Name", fill = "C14_Climate_Name",
    shape = "C14_Climate_Name", alpha = 0.5,
    size = 2,
    constraint_vec_style = vec_constraint(1.5, alpha = 0.8),
    constraint_lab_style = constraint_lab_style(size = 5, alpha = 0.5)
  ) +
  scale_color_manual(values = c25) +
  scale_fill_manual(values = c25) +
  scale_shape_manual(values = shapes) +
  ggside::geom_xsideboxplot(aes(fill = C14_Climate_Name, y = C14_Climate_Name), orientation = "y", varwidth = FALSE, show.legend = FALSE) +
  ggside::geom_ysideboxplot(aes(fill = C14_Climate_Name, x = C14_Climate_Name), orientation = "x", varwidth = FALSE, show.legend = FALSE) +
  ggside::theme_ggside_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(face = "bold", size=16),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 14)
  ) +
  labs(
    x = "RDA1 [10.2%]",
    y = "RDA2 [3.5%]",
    title = "Redundancy Analysis (RDA) of Eukaryotic Communities by C14 Climate Periods",
    caption = "159 samples & 1263 taxa (Phylum); tax_transform = clr"
  )

# Display
rda_timecoi


## Putting plots together
rda_zonecoi|rda_timecoi
