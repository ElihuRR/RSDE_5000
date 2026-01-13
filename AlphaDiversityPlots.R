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


## Figure 3.- Hill number alpha-diversity in the Eastern Red Sea aggregated by time periods
# Make sure to set correct directory where .rds files are located

# Color palette for plots:
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


## First Euka02
phyeu<-readRDS("physeqeukarsdefinal.rds")
phyeu

# For graphics without NAs:
phyclean<-subset_taxa(phyeu, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phyclean2<-subset_taxa(phyclean2, Class!="NA")
phyclean2

# Defining the order of the periods as factors in the sample data of the physeq
sample_data(phyclean2)$C14_Climate_Name <- factor(sample_data(phyclean2)$C14_Climate_Name, 
                                                  levels = c("Late Holocene Arid Transition", 
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
                                                             "Recent Anthropogenic Warming"))

## STEP 1 – get the alpha metrics
alpha_raw <- estimate_richness(
  phyclean2,
  measures = c("Observed", "Shannon", "Simpson")
)

## (optional sanity check: should be same order as sample_names)
# all(rownames(alpha_raw) == sample_names(phyclean2))

## STEP 2 – write them into sample_data(phyclean2)
sd <- as.data.frame(sample_data(phyclean2))

sd$Observed <- alpha_raw$Observed
sd$Shannon  <- alpha_raw$Shannon
sd$Simpson  <- alpha_raw$Simpson

sample_data(phyclean2) <- sd  # put back into the phyloseq object

## STEP 3 – pull out sample_data as a regular data frame
## Pull sample_data and compute Hill numbers
meta2 <- as(sample_data(phyclean2), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(
    D0 = Observed,                        # Hill q = 0 (richness)
    D1 = exp(Shannon),                    # Hill q = 1 (Shannon-based)
    D2 = ifelse(Simpson == 1,
                NA_real_,
                1 / (1 - Simpson))        # Hill q = 2 (Simpson-based)
  )

## STEP 4 – long format for plotting
alpha_long <- meta2 %>%
  dplyr::select(SampleID, C14_Climate_Name, D0, D1, D2) %>%
  pivot_longer(
    cols      = c(D0, D1, D2),
    names_to  = "Metric",
    values_to = "Hill_diversity"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("D0", "D1", "D2"),
      labels = c("q = 0 (richness)",
                 "q = 1 (Shannon-based)",
                 "q = 2 (Simpson-based)")
    )
  )

## Function for creating the plot ##

panel_base <- function(df, panel_title) {
  ggplot(
    df,
    aes(x = C14_Climate_Name,
        y = Hill_diversity)
  ) +
    geom_boxplot(
      aes(fill   = C14_Climate_Name,
          colour = C14_Climate_Name),
      alpha        = 0.2,
      outlier.alpha = 0.6
    ) +
    geom_jitter(
      aes(colour = C14_Climate_Name),
      width       = 0.15,
      size        = 1.5,
      alpha       = 0.7,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c25, name = NULL) +
    scale_color_manual(values = c25, guide = "none") +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(
      title = panel_title,
      x = "Time Period",
      y = "",
    ) +
    theme_bw() +
    theme(
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.text.x      = element_blank(),
      axis.text.y      = element_text(size = 12),
      axis.title.x     = element_text(size = 16, face = "bold"),
      axis.title.y     = element_text(size = 16, face = "bold"),
      plot.title       = element_text(size = 14, face = "bold"),
      legend.text      = element_text(size = 18)
      # NOTE: no legend.position here – patchwork will handle it
    )
}

p_q0 <- panel_base(
  filter(alpha_long, Metric == "q = 0 (richness)"),
  "q = 0 (richness)"
)

p_q1 <- panel_base(
  filter(alpha_long, Metric == "q = 1 (Shannon-based)"),
  "q = 1 (Shannon-based)"
)

p_q2 <- panel_base(
  filter(alpha_long, Metric == "q = 2 (Simpson-based)"),
  "q = 2 (Simpson-based)"
)


combinedeuka <- (p_q0 | p_q1 | p_q2) +
  plot_layout(guides = "collect") +          # <- single shared legend
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  theme(legend.position = "bottom")          # put shared legend at bottom

combinedeuka


## Now we do COI
phyco<-readRDS("physeqco1rsdefinal.rds")
phyco
sample_data(phyco)

phyclean<-subset_taxa(phyco, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean2<-subset_taxa(phyclean, Class!="NA")
# Sanity check
sample_data(phyclean2)
tax_table(phyclean2)


## STEP 1 – get the alpha metrics
alpha_raw <- estimate_richness(
  phyclean2,
  measures = c("Observed", "Shannon", "Simpson")
)

## (optional sanity check: should be same order as sample_names)
# all(rownames(alpha_raw) == sample_names(phyclean2))

## STEP 2 – write them into sample_data(phyclean2)
sd <- as.data.frame(sample_data(phyclean2))

sd$Observed <- alpha_raw$Observed
sd$Shannon  <- alpha_raw$Shannon
sd$Simpson  <- alpha_raw$Simpson

sample_data(phyclean2) <- sd  # put back into the phyloseq object

## STEP 3 – pull out sample_data as a regular data frame
## Pull sample_data and compute Hill numbers
meta2 <- as(sample_data(phyclean2), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  mutate(
    D0 = Observed,                        # Hill q = 0 (richness)
    D1 = exp(Shannon),                    # Hill q = 1 (Shannon-based)
    D2 = ifelse(Simpson == 1,
                NA_real_,
                1 / (1 - Simpson))        # Hill q = 2 (Simpson-based)
  )

## STEP 4 – long format for plotting
alpha_long <- meta2 %>%
  dplyr::select(SampleID, C14_Climate_Name, D0, D1, D2) %>%
  pivot_longer(
    cols      = c(D0, D1, D2),
    names_to  = "Metric",
    values_to = "Hill_diversity"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("D0", "D1", "D2"),
      labels = c("q = 0 (richness)",
                 "q = 1 (Shannon-based)",
                 "q = 2 (Simpson-based)")
    )
  )


## Creating the plots for both Euka02 + COI

p_q0c <- panel_base(
  filter(alpha_long, Metric == "q = 0 (richness)"),
  "q = 0 (richness)"
)

p_q1c <- panel_base(
  filter(alpha_long, Metric == "q = 1 (Shannon-based)"),
  "q = 1 (Shannon-based)"
)

p_q2c <- panel_base(
  filter(alpha_long, Metric == "q = 2 (Simpson-based)"),
  "q = 2 (Simpson-based)"
)

## Now putting the plots together

row_euka <- p_q0 | p_q1 | p_q2   # top row
row_coi <- p_q0c | p_q1c | p_q2c

row_euka <- row_euka & theme(legend.position = "none")
row_coi <- row_coi & theme(legend.position = "none")

combined_6 <- (p_q0 | p_q1 | p_q2 |
                 p_q0c | p_q1c | p_q2c) +
  plot_layout(nrow = 2, guides = "collect") +     # collect -> 1 legend
  plot_annotation(
    tag_levels = "a",                             # a), b), c), d), e), f)
    tag_suffix = ")"
  ) &
  theme(legend.position = "bottom")

combined_6
