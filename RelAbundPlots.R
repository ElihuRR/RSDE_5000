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


## Figure 2.- Relative abundances of phyla in the Eastern Red Sea aggregated by time periods
# Make sure to set correct directory where .rds files are located

# First, Euka02 (18Sv7 primer)
phyeu<-readRDS("physeqeukarsdefinal.rds")
phyeu

# Remove taxa with unreliable taxonomy, as well as those that contain NA in Phylum
phyclean<-subset_taxa(phyeu, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyeu<-subset_taxa(phyclean, Phylum!="NA")

# Defining the order for the time periods:
periodnames<-c("Late Holocene Arid Transition","Bronze Age Recovery","4.2 ka Drought","Bronze Age Warm Peak","3.2 ka Late Bronze Age Drought","Iron Age Cold Epoch","Roman Warm Period","Late Antique Little Ice Age","Medieval Warm Period","Little Ice Age","Early Industrial Revolution","Recent Anthropogenic Warming")

# Step 1: Agglomerate to Phylum
y1 <- tax_glom(phyeu, taxrank = 'Phylum')
# Step 2: Merge samples by C14_Climate_Name
y2 <- merge_samples(y1, group = "C14_Climate_Name")

# Step 3: Rebuild sample_data to retain C14_Climate_Name as a variable
# (since merge_samples drops all sample_data except sample names)
# Create a new sample_data with Time Periods as rownames
meta <- sample_data(phyeu)
times <- unique(meta$C14_Climate_Name)
new_meta <- data.frame(C14_Climate_Name = times)
rownames(new_meta) <- as.character(times)
sample_data(y2) <- sample_data(new_meta)

# Step 4: Obtain relative abundance of each phylum
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4eu <- psmelt(y3) # create dataframe from phyloseq object
y4eu$Phylum <- as.character(y4eu$Phylum) #convert to character
y4eu$C14_Climate_Name <- factor(as.character(y4eu$C14_Climate_Name), levels = periodnames)

# Step 5: # Identify low-abundance phyla (≤1% across the whole dataset)
total_abund <- y4eu %>%
  group_by(Phylum) %>%
  summarise(total_abundance = sum(Abundance, na.rm = TRUE)) %>%
  mutate(relative_abundance = total_abundance / sum(total_abundance))

low_abund_phyla <- total_abund %>%
  filter(relative_abundance <= 0.01) %>%
  pull(Phylum)

y4eu$Phylum_grouped <- ifelse(y4eu$Phylum %in% low_abund_phyla, "Other(<1% total)", y4eu$Phylum)
y4eu$Phylum_grouped <- factor(y4eu$Phylum_grouped, levels = sort(unique(y4eu$Phylum_grouped)))



# Now, same process with COI primer
phycoi<-readRDS("physeqco1rsdefinal.rds")
phycoi

# Step 1: Agglomerate to Phylum
y1 <- tax_glom(phycoi, taxrank = 'Phylum')
# Step 2: Merge samples by C14_Climate_Name
y2 <- merge_samples(y1, group = "C14_Climate_Name")

# Step 3: Rebuild sample_data to retain C14_Climate_Name as a variable
# (since merge_samples drops all sample_data except sample names)
# Create a new sample_data with Time Period as rownames
meta <- sample_data(phyeu)
times <- unique(meta$C14_Climate_Name)
new_meta <- data.frame(C14_Climate_Name = times)
rownames(new_meta) <- as.character(times)
sample_data(y2) <- sample_data(new_meta)

# Step 4: Obtain relative abundance of each phylum
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4co <- psmelt(y3) # create dataframe from phyloseq object
y4co$Phylum <- as.character(y4co$Phylum) #convert to character
y4co$C14_Climate_Name <- factor(as.character(y4co$C14_Climate_Name), levels = periodnames)

# Step 5: Identify low-abundance phyla (≤1% across the whole dataset)
total_abund <- y4co %>%
  group_by(Phylum) %>%
  summarise(total_abundance = sum(Abundance, na.rm = TRUE)) %>%
  mutate(relative_abundance = total_abundance / sum(total_abundance))

low_abund_phyla <- total_abund %>%
  filter(relative_abundance <= 0.01) %>%
  pull(Phylum)

y4co$Phylum_grouped <- ifelse(y4co$Phylum %in% low_abund_phyla, "Other(<1% total)", y4co$Phylum)
y4co$Phylum_grouped <- factor(y4co$Phylum_grouped, levels = sort(unique(y4co$Phylum_grouped)))



## Now, before plotting, we generate a color palette that is consistent across the two primers:

# Color palette for the different phyla, has about 40 different colors #
fancycol<-c("#fddba4","#50d5fe","#f3988a","#b7e899","#f999b3","#7de8c7","#f2b6f6","#b3fec6","#faa73d","#a2b870",
            "#c4befe","#edf55e","#88c7fe","#dfcb81","#24fefd","#fdfbb4","#d460c3","#febb90","#23c6c8","#fea6ab",
            "#5bccaf","#fed3f2","#a2e8ad","#c8a1bf","#f8cd56","#8cb1ca","#c8fe7a","#c9a590","#b9fefb","#fe886b",
            "#fe8cb3","#f4e0c6","#7eb6b6","#9fb47c","#8c7acb","#d8c3f2","#b9fee1","#8ab6a1","#c0e1cf","#28fd59")

# New colors, assigned to each phylum:
# Step 1: Combine phylum lists from both Euka and CO1 datasets
all_phyla <- unique(c(unique(y4eu$Phylum_grouped), unique(y4co$Phylum_grouped)))
# Step 2: Sort phyla alphabetically (for reproducibility)
all_phyla <- sort(all_phyla)
# Step 3: Create named color vector
# Use length(all_phyla) to ensure you’re not exceeding available colors
fancycol_vector <- setNames(fancycol[1:length(all_phyla)], all_phyla)


# Now we plot
# 18S v7
eukrel <- ggplot(data = y4eu, aes(x = C14_Climate_Name, y = Abundance, fill = Phylum_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = fancycol_vector) +  # Ensure "1% or less phyla" has a color
  labs(x = "Decade", y = "Relative abundance", fill="Phylum") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 32, face = "bold"),
    strip.text = element_text(size = 14),
    legend.position = "right",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  guides(fill = guide_legend(nrow = 18)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
eukrel<-eukrel+labs(x=NULL) + ggtitle("a)") + theme(plot.title = element_text(size = 22))
eukrel

# COI
coirel <- ggplot(data = y4co, aes(x = C14_Climate_Name, y = Abundance, fill = Phylum_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = fancycol_vector) +  # Ensure "1% or less phyla" has a color
  labs(x = "Time Period", y = "Relative abundance", fill="Phylum") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 32, face = "bold"),
    strip.text = element_text(size = 14),
    legend.position = "right",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  guides(fill = guide_legend(nrow = 18)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
coirel<-coirel+ ggtitle("b)") + theme(plot.title = element_text(size = 22))
coirel

# Now putting the plots together
phyrelplots<-eukrel/coirel
phyrelplots