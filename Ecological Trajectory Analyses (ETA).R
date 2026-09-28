# This script contains code utilized to generate analyses and plots for the manuscript:
#"Millennia of stability, decades of change: tracing 5,000 years of marine biodiversity shifts in the Red Sea"
# By Rivera Rosas et al., 2026.

# Load necessary packages:
library(ecotraj)
library(RColorBrewer)
library(phyloseq)
library(speedyseq)
library(ggplot2)
library(ape)
library(vegan)
library(picante)
library(cowplot)
library(DESeq2)
library(patchwork)
library(BiocManager)
library(microbiome)

## 1.- Building the distance matrix 
## This is the single most important input to ALL ecotraj analyses.
## Everything ecotraj does is geometry computed ON this matrix.
## 
## We use the same CLR-Euclidean approach

# Load the data and clean:
phyeu<-readRDS("physeqeukarsdefinal.rds")
phyclean<-subset_taxa(phyeu, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean
tax_table(phyclean)
# A total of 923 taxa were eliminated, corresponds to around 50% of reads 
# Of note is that with Euka02 specifically, there are still many Phyla that contain NA even after this cleaning, therefore:
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phyclean2<-subset_taxa(phyclean2, Class!="NA")
phyclean2


## Step 1a: Extract OTU table from the cleaned phyloseq object

otu_full <- as.data.frame(otu_table(phyclean2))
if (!taxa_are_rows(phyclean2)) otu_full <- t(otu_full)
## otu_full is now: taxa as rows, samples as columns

## Step 1b: CLR transform with pseudocount
otu_clr <- apply(otu_full + 0.01, 2,
                 function(x) log(x) - mean(log(x)))
## Each column is now a CLR-transformed community vector for one sample

## Step 1c: Compute Euclidean distance matrix
## dist() works on ROWS, so we transpose first
## The result is a sample × sample distance matrix
D <- dist(t(otu_clr), method = "euclidean")

## Quick sanity check — should be n_samples × n_samples
cat("Distance matrix dimensions:", attr(D, "Size"), "×", attr(D, "Size"), "\n")
## Expected: 175 × 175

## Step 1d: Extract and prepare metadata
meta_full <- as.data.frame(as.matrix(sample_data(phyclean2))) %>%
  rownames_to_column("SampleID") %>%
  mutate(
    New_YearC14 = as.numeric(New_YearC14),
    Zone        = factor(Zone, levels = c("North", "Central", "South"))
  )


## CHUNK 2 — Build the trajectory object 

core_col <- "Dive_Name"   

## Step 2b: ensure metadata is sorted consistently with the distance matrix
## The row/column order of D must match the order of samples in meta_full
## dist() preserves the row order of the input matrix, which came from
## t(otu_clr) — so we align meta_full to that same order

sample_order <- rownames(as.matrix(D))   # extract sample order from D

meta_ordered <- meta_full %>%
  filter(SampleID %in% sample_order) %>%
  slice(match(sample_order, SampleID))   # enforce identical order

## Sanity check — these must match exactly before proceeding
stopifnot(all(meta_ordered$SampleID == sample_order))
cat("Sample order confirmed: OK\n")


## Step 2c: build the sites and surveys vectors
## Sites = core ID repeated for each sample in that core
## Surveys = integer rank within each core, ordered oldest -> most recent
## We sort by year WITHIN each core so survey 1 = oldest layer

meta_ordered <- meta_ordered %>%
  group_by(across(all_of(core_col))) %>%
  arrange(New_YearC14, .by_group = TRUE) %>%  # oldest layer = survey 1
  mutate(survey_order = row_number()) %>%      # 1, 2, 3... per core
  ungroup()

## Extract the two vectors ecotraj needs
sites_vec   <- meta_ordered[[core_col]]        # core ID per sample
surveys_vec <- meta_ordered$survey_order        # integer order per sample


## Quick check — each core should have surveys 1:n with no gaps
meta_ordered %>%
  group_by(across(all_of(core_col))) %>%
  summarise(
    n_samples  = n(),
    min_survey = min(survey_order),
    max_survey = max(survey_order),
    year_start = min(New_YearC14),
    year_end   = max(New_YearC14),
    zone       = first(Zone)
  ) %>%
  print()


## Step 2d: re-sort the distance matrix to match meta_ordered
## After the arrange() above, meta_ordered rows may have shifted
## D must reflect the same new order

new_order <- match(meta_ordered$SampleID, sample_order)
D_ordered <- as.dist(as.matrix(D)[new_order, new_order])


## Step 2e: define the trajectory object
## This is the core ecotraj function, it stores the geometric information for all downstream analyses

## Build a times vector aligned to meta_ordered
times_vec <- meta_ordered$New_YearC14

traj <- defineTrajectories(
  d       = D_ordered,
  sites   = sites_vec,
  surveys = surveys_vec,
  times   = times_vec   
)

print(traj)


## Step 2f: add metadata back as trajectory attributes
## ecotraj doesn't store years internally, we attach them separately
## for use in plotting and rate-of-change calculations later
traj_meta <- meta_ordered %>%
  dplyr::select(SampleID, 
                all_of(core_col), 
                Zone, 
                New_YearC14, 
                survey_order)

cat("\nTrajectory object built successfully.\n")
cat("Cores:", length(unique(sites_vec)), "\n")
cat("Total samples:", length(sites_vec), "\n")



## CHUNK 3 — Segment lengths over time 
##
## trajectoryLengths() computes:
##   - "segments": a matrix of distances between each consecutive pair
##                 of surveys within every trajectory (core)
##                 rows = trajectories, columns = segments (1→2, 2→3, etc.)
##   - "total":    the sum of all segment lengths per trajectory —
##                 the total distance the community has travelled
##
## This is purely geometric, it reads distances from D_ordered through
## the trajectory structure defined in traj.

traj
seg_lengths <- trajectoryLengths(traj)
seg_speeds <- trajectorySpeeds(traj)

print(seg_speeds)
print(seg_lengths)

## Step 3a: reshape into long format for plotting 
## The segment matrix has cores as rows and segment positions as columns
## so we need to attach the actual calendar years to each segment midpoint

seg_df <- seg_lengths %>%
  rownames_to_column("Core") %>%
  pivot_longer(
    cols      = starts_with("S"),   # selects S1-S24, excludes Path
    names_to  = "segment_index",
    values_to = "segment_length"
  ) %>%
  filter(!is.na(segment_length)) %>%
  mutate(segment_index = as.integer(str_extract(segment_index, "\\d+")))

## Reshape speeds the same way as seg_lengths
## The structure is identical, same columns S1-S24 plus Path (total speed)
speed_df <- seg_speeds %>%
  rownames_to_column("Core") %>%
  pivot_longer(
    cols      = starts_with("S"),
    names_to  = "segment_index",
    values_to = "segment_speed"
  ) %>%
  filter(!is.na(segment_speed)) %>%
  mutate(segment_index = as.integer(str_extract(segment_index, "\\d+")))


## Step 3b: attach year midpoints 
## Each segment connects survey i to survey i+1
## The midpoint year = mean of those two survey years
## We pull this from traj_meta built in Chunk 2

year_midpoints <- traj_meta %>%
  group_by(Dive_Name) %>%
  arrange(survey_order, .by_group = TRUE) %>%
  mutate(
    ## midpoint between this survey and the next
    year_mid = (New_YearC14 + lead(New_YearC14)) / 2,
    ## segment index = the "from" survey number
    segment_index = survey_order
  ) %>%
  filter(!is.na(year_mid)) %>%   # last survey has no "next", drop it
  dplyr::select(Core = Dive_Name, Zone, segment_index, year_mid) %>%
  ungroup()

seg_plot_df <- left_join(seg_df, year_midpoints,
                         by = c("Core", "segment_index")) %>%
  mutate(Zone = factor(Zone, levels = c("North", "Central", "South")))

## Join year midpoints (reuse year_midpoints from Chunk 3)
speed_plot_df <- left_join(speed_df, year_midpoints,
                           by = c("Core", "segment_index")) %>%
  mutate(Zone = factor(Zone, levels = c("North", "Central", "South")))

## Quick check: no NAs in year_mid after join?
stopifnot(all(!is.na(seg_plot_df$year_mid)))

## Step 3c: total path lengths 
total_lengths <- seg_lengths %>%
  rownames_to_column("Core") %>%
  dplyr::select(Core, total_length = Path) %>%
  left_join(
    traj_meta %>%
      distinct(Dive_Name, Zone) %>%
      rename(Core = Dive_Name),
    by = "Core"
  ) %>%
  arrange(Zone, Core)

print(total_lengths)


## Step 3d: plot segment lengths over time for a quick check
zone_colors <- c("North"   = "#2ca25f",
                 "Central" = "#8856a7",
                 "South"   = "#e34a33")

p_seg <- ggplot(seg_plot_df,
                aes(x = year_mid, y = segment_length,
                    color = Zone, group = Core)) +
  
  ## Individual core trajectories as semi-transparent lines + points
  ## group = Core ensures each core gets its own connected line
  ## This is important — you don't want segments from different cores joined
  geom_line(alpha  = 0.35, linewidth = 0.5) +
  geom_point(aes(shape = Core),
             size  = 1.8,
             alpha = 0.55) +
  
  ## Zone-level GAM smooth across all cores in that zone
  ## This is fit across cores, treating year_mid as continuous predictor
  ## se = TRUE shows the cross-core uncertainty in the zone trend
  geom_smooth(aes(group = Zone),   # override group = Core for the smooth
              method      = "gam",
              formula     = y ~ s(x, bs = "cs"),
              method.args = list(method = "REML"),
              se          = TRUE,
              alpha       = 0.15,
              linewidth   = 0.9) +
  
  scale_color_manual(values = zone_colors) +
  
  scale_x_continuous(
    breaks = function(lims) seq(lims[1], lims[2], length.out = 10),
    labels = function(x) round(x, 0),
    expand = c(0.01, 0)
  ) +
  
  facet_wrap(~ Zone, ncol = 1, scales = "free_x") +
  
  labs(
    x     = "Year (BCE/CE)",
    y     = "Segment length\n(CLR–Euclidean dissimilarity)",
    color = "Region",
    shape = "Core"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold", size = 11),
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 8),
    axis.title.y       = element_text(size = 9)
  )

p_seg



## CHUNK 4

## Step 4a: compute NC natively
nc_native <- trajectoryLengths(traj, relativeToInitial = TRUE)

## Step 4b: reshape the data
nc_df <- nc_native %>%
  rownames_to_column("Core") %>%
  pivot_longer(
    cols      = starts_with("Lt"),
    names_to  = "survey_label",
    values_to = "NC"
  ) %>%
  filter(!is.na(NC)) %>%
  mutate(survey_index = as.integer(str_extract(survey_label, "(?<=_t)\\d+")))

## Step 4c: join calendar years and zone from traj_meta
nc_df <- nc_df %>%
  left_join(
    traj_meta %>%
      dplyr::rename(Core = Dive_Name) %>%
      dplyr::select(Core, survey_order, New_YearC14, Zone),
    by = c("Core", "survey_index" = "survey_order")
  ) %>%
  mutate(Zone = factor(Zone, levels = c("North", "Central", "South")))


## Step 4d: compute RDT, there is no native function, so we do manual derivation
## RDT(n) = NC(n-1) - NC(n)
## Negative = departing further from baseline
## Positive = recovering back toward baseline
nc_df <- nc_df %>%
  group_by(Core) %>%
  arrange(survey_index, .by_group = TRUE) %>%
  mutate(RDT = lag(NC) - NC) %>%
  ungroup()

## Summary by zone
nc_df %>%
  group_by(Zone) %>%
  summarise(
    mean_NC       = mean(NC,  na.rm = TRUE),
    mean_RDT      = mean(RDT, na.rm = TRUE),
    pct_departing = mean(RDT < 0, na.rm = TRUE) * 100
  ) %>% print()


## CHUNK 5 - Combined figure: Speed + NC + RDT per zone 
## Three panels per zone: segment speed (top), NC (middle), RDT (bottom)
## Faceted by zone across columns

library(patchwork)   # for combining ggplot panels
library(ggh4x)       # for nested facets 

## Step 5a: prepare a clean combined data frame 
## Speeds are at segment midpoints (year_mid)
## NC and RDT are at survey points (New_YearC14)
## We keep them separate and bind with a metric label for faceting

speed_long <- speed_plot_df %>%
  dplyr::select(Core, Zone, year_mid, segment_speed) %>%
  dplyr::rename(Year = year_mid, Value = segment_speed) %>%
  mutate(Metric = "Speed\n(dissimilarity yr⁻¹)")

nc_long <- nc_df %>%
  dplyr::select(Core, Zone, New_YearC14, NC) %>%
  dplyr::rename(Year = New_YearC14, Value = NC) %>%
  mutate(Metric = "Net Change\n(distance from baseline)")

rdt_long <- nc_df %>%
  dplyr::select(Core, Zone, New_YearC14, RDT) %>%
  filter(!is.na(RDT)) %>%
  dplyr::rename(Year = New_YearC14, Value = RDT) %>%
  mutate(Metric = "RDT\n(recovering ↑  /  departing ↓)")

all_metrics <- bind_rows(speed_long, nc_long, rdt_long) %>%
  mutate(
    Zone = factor(Zone, levels = c("North", "Central", "South")),
    ## Control panel order top to bottom
    Metric = factor(Metric, levels = c(
      "Speed\n(dissimilarity yr⁻¹)",
      "Net Change\n(distance from baseline)",
      "RDT\n(recovering ↑  /  departing ↓)"
    ))
  )


## Step 5b: build final plot 
zone_colors <- c("North"   = "#2ca25f",
                 "Central" = "#8856a7",
                 "South"   = "#e34a33")

## RDT needs a horizontal reference line at 0 — departing below, recovering above
## We add this only to the RDT panel using a conditional data argument
rdt_ref <- data.frame(
  Metric = factor("RDT\n(recovering ↑  /  departing ↓)",
                  levels = levels(all_metrics$Metric))
)


## Panel label data frame: one row per Metric × Zone combination
## Labels placed in top-left of each panel

row_labels <- c(
  "Speed\n(dissimilarity yr\u207b\u00b9)"         = "a) Speed (dissimilarity yr\u207b\u00b9)",
  "Net Change\n(distance from baseline)"           = "b) Net Change (distance from t\u2080)",
  "RDT\n(recovering \u2191  /  departing \u2193)" = "c) Rate of Directional Turnover"
)
max_display_year <- 2024

p_combined <- ggplot(all_metrics,
                     aes(x = Year, y = Value,
                         color = Zone, fill = Zone,
                         group = Core)) +
  
  geom_line(alpha = 0.30, linewidth = 0.45) +
  
  geom_point(aes(shape = Core), size = 2.5, alpha = 0.70) +
  
  geom_hline(
    data        = rdt_ref,
    aes(yintercept = 0),
    linetype    = "dashed",
    linewidth   = 0.45,
    color       = "grey40",
    inherit.aes = FALSE
  ) +
  
  geom_smooth(
    aes(group = Zone),
    method      = "gam",
    formula     = y ~ s(x, bs = "cs"),
    method.args = list(method = "REML"),
    se          = TRUE,
    alpha       = 0.15,
    linewidth   = 1.0
  ) +
  
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values  = zone_colors) +
  
  scale_shape_manual(
    values = c(
      "CHR0183" = 16, "CHR0191" = 17, "CHR0256" = 15,
      "CHR0264" = 18, "CHR0280" = 1,  "CHR0287" = 2,
      "CHR0292" = 0
    )
  ) +
  
  ## Fix 2: pretty breaks + hard cap at 2024
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 6),
    labels = function(x) round(x, 0),
    limits = function(x) c(x[1], min(x[2], max_display_year)),
    expand = expansion(mult = c(0.03, 0.01))
  ) +
  
  ## switch = "y" moves row strips to the LEFT so they act as y-axis labels
  ## labeller renders the full descriptive titles defined in row_labels
  facet_grid(
    Metric ~ Zone,
    scales   = "free",
    switch   = "y",
    labeller = labeller(Metric = row_labels)
  ) +
  
  labs(
    x     = "Year (BCE/CE)",
    y     = NULL,        
    color = "Region",
    fill  = "Region",
    shape = "Core"
  ) +
  
  guides(color = "none", fill = "none") +
  
  theme_bw(base_size = 13) +
  theme(
    strip.background    = element_rect(fill = "grey92", color = NA),
    strip.text.x        = element_text(face = "bold", size = 16),
    
    strip.text.y.left = element_text(face  = "bold",
                                     size  = 16,
                                     angle = 90,   
                                     hjust = 0.5,
                                     vjust = 0.5),
    
    strip.placement     = "outside",
    strip.background.y  = element_blank(),  
    
    legend.position     = "bottom",
    legend.box          = "horizontal",
    legend.text         = element_text(size = 12),
    legend.title        = element_text(size = 14, face = "bold"),
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_line(color = "grey92"),
    axis.text.x         = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y         = element_text(size = 16),
    axis.title.x        = element_text(size = 18, face= 'bold'),
    axis.title.y        = element_text(size = 18, face= 'bold'),
    
    ## Extra left margin so the two-line strip labels don't get clipped
    plot.margin         = margin(t = 5, r = 10, b = 5, l = 15, unit = "pt")
  )

p_combined

ggsave("Fig_ETA_combined_Euka02.pdf",
       plot   = p_combined,
       width  = 14,
       height = 14,
       device = "pdf")


ggsave("Fig_ETA_combined_Euka02.png",
       plot   = p_combined,
       width  = 14,
       height = 14,
       dpi=300)


## CHUNK 6-Trajectory PCoA

## Step 1: run trajectoryPCoA() to extract coordinates.
## We redirect the base R plot to a null device since we are
## rebuilding it in ggplot2. The return value is the cmdscale
## coordinate matrix: one row per sample in distance matrix order.

## eig = TRUE returns eigenvalues so we can compute % variance explained.
## k = ncol - 1 gets all axes; we use axes 1 and 2 for the figure.

## Suppress the base R plot, capture coordinates
invisible(pdf(NULL))
pcoa_result <- trajectoryPCoA(traj, axes = c(1, 2))
invisible(dev.off())


## Run cmdscale directly on the trajectory distance matrix
## eig = TRUE returns eigenvalues for variance explained calculation
## k = 2 gives us the two axes we need for plotting
pcoa_eig <- cmdscale(traj$d, k = 2, eig = TRUE)

## Extract coordinates — rows match traj$metadata order exactly
pcoa_coords <- as.data.frame(pcoa_eig$points) %>%
  setNames(c("PC1", "PC2"))

## Variance explained
total_var    <- sum(pcoa_eig$eig[pcoa_eig$eig > 0])
pc1_var      <- round(pcoa_eig$eig[1] / total_var * 100, 1)
pc2_var      <- round(pcoa_eig$eig[2] / total_var * 100, 1)

cat("PC1:", pc1_var, "%\nPC2:", pc2_var, "%\n")
cat("Total (PC1+PC2):", round(pc1_var + pc2_var, 1), "%\n")

## Attach metadata, row order of cmdscale output matches traj$metadata
pcoa_df <- pcoa_coords %>%
  bind_cols(traj$metadata) %>%
  rename(Core = sites, Survey = surveys, Year = times) %>%
  left_join(
    traj_meta %>% distinct(Dive_Name, Zone) %>% rename(Core = Dive_Name),
    by = "Core"
  ) %>%
  mutate(Zone = factor(Zone, levels = c("North", "Central", "South")))

## Quick sanity check
cat("\nRows in pcoa_df:", nrow(pcoa_df), "\n")
cat("Cores present:", paste(unique(pcoa_df$Core), collapse = ",
                            "), "\n")
head(pcoa_df)

## Step 2: Building the data layers

## Arrow segments with year of the originating survey
arrows_df <- pcoa_df %>%
  arrange(Core, Survey) %>%
  group_by(Core, Zone) %>%
  mutate(
    PC1_end = lead(PC1),
    PC2_end = lead(PC2)
  ) %>%
  filter(!is.na(PC1_end)) %>%
  ungroup()

## Start points (open circle) and end points (filled diamond) per core
start_df <- pcoa_df %>%
  group_by(Core, Zone) %>%
  slice(1) %>%
  ungroup()

end_df <- pcoa_df %>%
  group_by(Core, Zone) %>%
  slice(n()) %>%
  ungroup()

## Convex hull per zone — shows spatial occupation and overlap
hull_df <- pcoa_df %>%
  group_by(Zone) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()


## Step 3: Creating the figure

pcoa_plotEu<-ggplot() +
  
  ## Zone convex hulls — border color derived from fill, no separate color aes
  geom_polygon(
    data      = hull_df,
    aes(x = PC1, y = PC2,
        fill  = Zone,
        color = after_scale(fill)),   ## border matches fill, no color scale conflict
    alpha     = 0.07,
    linewidth = 0.4,
    linetype  = "dashed"
  ) +
  
  geom_segment(
    data = arrows_df,
    aes(x = PC1, y = PC2,
        xend  = PC1_end, yend = PC2_end,
        color = Year),
    arrow     = arrow(length = unit(0.12, "cm"), type = "closed"),
    alpha     = 0.70,
    linewidth = 0.40
  ) +
  
  ## Start point per core — open circle
  geom_point(
    data  = start_df,
    aes(x = PC1, y = PC2, fill = Zone),
    shape  = 21,
    size   = 3.0,
    stroke = 1.0,
    color  = "grey20"
  ) +
  
  ## End point per core — filled diamond
  geom_point(
    data  = end_df,
    aes(x = PC1, y = PC2, fill = Zone),
    shape  = 23,
    size   = 3.0,
    stroke = 0.6,
    color  = "grey20"
  ) +
  
  scale_color_viridis_c(
    option = "turbo",
    name   = "Year (BCE/CE)",
    breaks = c(-3000, -1500, 0, 1000, 2000),
    labels = c("-3000", "-1500", "0", "1000", "2000")
  ) +
  
  scale_fill_manual(values = zone_colors, name = "Zone") +
  
  labs(
    x       = paste0("PC1 (", pc1_var, "%)"),
    y       = paste0("PC2 (", pc2_var, "%)")
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.caption     = element_text(size = 10, hjust = 0,
                                    color = "grey40",
                                    margin = margin(t = 6))
  )

pcoa_plotEu

ggsave("Fig_PCoA_ETA_18S.pdf",
       plot   = pcoa_plot,
       width  = 22,    
       height = 14,
       units  = "cm",
       device = "pdf")
