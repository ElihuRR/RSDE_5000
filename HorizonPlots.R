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

## Figure 5 and 6.- Horizon plots of metazoan taxa in time
# Make sure to set correct directory where .rds files are located

## First Euka02
# Data has been edited in Excel due to constraints in R and is now merged OTU with taxa
# Data was first subset to only metazoans.

# Now running the new prepanel with this version edited for phyla and class. However, taxa table is not required anymore  
phylum_otu<-read.csv("PhylumOTUMerged.csv")

# Now adding the metadata
meta<-read.csv("Meta_Trial5v2.csv", header=T)
meta

prepanel_dataall <- prepanel(
  otudata = phylum_otu,
  metadata = meta,
  subj = "All",
  thresh_prevalence = 30,
  thresh_abundance = 0.2)

# Generate the plot
all <- horizonplot(prepanel_dataall, aesthetics = horizonaes(
  title = "Red Sea Horizon Plot",
  xlabel = "Time",
  ylabel = "Metazoan Taxa",
  legendTitle = "Quartiles Relative to Taxon Median",
  legendPosition = "bottom"
))

all


## Now for COI
phylum_otu<-read.csv("PhylumOTUCO1RSDE.csv")
metad<-read.csv("MetaBiomeHorizon1.csv")

prepanel_dataall <- prepanel(
  otudata = phylum_otu,
  metadata = metad,
  subj = "All",
  thresh_prevalence = 30,
  thresh_abundance = 0.2)

str(prepanel_dataall)

# Generate the plot
all <- horizonplot(prepanel_dataall, aesthetics = horizonaes(
  title = "Red Sea Horizon Plot",
  xlabel = "Time",
  ylabel = "Animal Taxa",
  legendTitle = "Quartiles Relative to Taxon Median",
  legendPosition = "bottom"
))

all
