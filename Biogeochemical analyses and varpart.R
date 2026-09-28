# This script contains code utilized to generate analyses and plots for the manuscript:
#"Millennia of stability, decades of change: tracing 5,000 years of marine biodiversity shifts in the Red Sea"
# By Rivera Rosas et al., 2026.

# Load necessary packages:

library(phyloseq); library(speedyseq)
library(vegan)                        
library(tidyverse)                     
library(patchwork)
library(scales)

set.seed(1234)                         # permutation reproducibility

zone_colors <- c("North" = "#2ca25f", "Central" = "#8856a7", "South" = "#e34a33")

## Global per-core palette — one unique colour per core across ALL zones
core_colors <- c(
  CHR0280 = "#1b9e77", CHR0287 = "#d95f02", CHR0292 = "#7570b3",   # North
  CHR0183 = "#e7298a", CHR0256 = "#1f78b4", CHR0264 = "#e6ab02",   # Central
  CHR0191 = "#a6761d"                                              # South
)


## phyloseq object: Euka02 
phyeu <- readRDS("physeqeukarsdefinal.rds")
phyclean  <- subset_taxa(phyeu,   lowest_rank != "No match in database")
phyclean  <- subset_taxa(phyclean, lowest_rank != "Taxonomy unreliable")
phyclean2 <- subset_taxa(phyclean,  Phylum != "NA")
phyclean2 <- subset_taxa(phyclean2, Class  != "NA")
phyeufinal <- phyclean2
phyeufinal


##  phyloseq object: COI 
phyco<-readRDS("physeqco1rsdefinal.rds")
phyco
phyclean<-subset_taxa(phyco, lowest_rank!="No match in database")
phyclean<-subset_taxa(phyclean, lowest_rank!="Taxonomy unreliable")
phyclean2<-subset_taxa(phyclean, Phylum!="NA")
phycoifinal<-subset_taxa(phyclean2, Class!="NA")


## Full sample data:
samdata <- read.csv("RSDEFinalFinalPb210SampleData.csv")

## The 7 dated cores actually used (CHR0219 / CHR0229 excluded: no chronology)
dated_cores <- c("CHR0183","CHR0191","CHR0256","CHR0264","CHR0280","CHR0287","CHR0292")


## A.  PROXY TABLE + DERIVED PROXIES  

## We do NOT average a band into one value. We use established RATIOS
## (Ti/Ca, Fe/Al, Si/Al, Sr/Ca) and Al-normalised enrichment factors (EF) for
## the metals, each referenced to that core's pre-industrial baseline.


## A.0  robust column getter because read.csv is buggy
getcol <- function(df, want) {
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  idx <- match(norm(want), norm(names(df)))
  if (is.na(idx))
    stop("getcol(): column not found: '", want, "'\n  Available:\n  ",
         paste(names(df), collapse = ", "))
  suppressWarnings(as.numeric(df[[idx]]))
}

sd_use <- samdata %>% dplyr::filter(Dive_Name %in% dated_cores)

proxy <- data.frame(
  SampleName = sd_use$SampleName,
  Core       = sd_use$Dive_Name,
  Zone       = factor(sd_use$Zone, levels = c("North","Central","South")),
  Depth      = getcol(sd_use, "Depth"),
  Year       = getcol(sd_use, "New_YearC14"),          # negative = BCE
  
  ## Band 1 — productivity / organic matter
  TOC   = getcol(sd_use, "TOC_Percentage"),
  Opal  = getcol(sd_use, "Opal_Percentage"),
  CaCO3 = getcol(sd_use, "CaCO3_Percentage"),
  CN    = getcol(sd_use, "C/N"),
  Sr_Ca = getcol(sd_use, "Sr(mg/kg)") / (getcol(sd_use, "Ca_Percentage") * 10000),
  
  ## Band 2 — terrigenous / climate (dust–aridity)
  Ti_Ca = getcol(sd_use, "Ti_Percentage") / getcol(sd_use, "Ca_Percentage"),
  Si_Al = getcol(sd_use, "Si_Percentage") / getcol(sd_use, "Al_Percentage"),
  grain = getcol(sd_use, "d0.5(μm)"),
  
  ## Band 3 — redox / oxygenation  (Fe_Al belongs here, not terrigenous)
  Mo    = getcol(sd_use, "Mo(mg/kg)"),
  V     = getcol(sd_use, "V(mg/kg)"),
  Mn    = getcol(sd_use, "Mn(mg/kg)"),
  S_pct = getcol(sd_use, "S_Percentage"),
  Fe_Al = getcol(sd_use, "Fe_Percentage") / getcol(sd_use, "Al_Percentage"),
  
  ## raw inputs retained for enrichment factors
  Al = getcol(sd_use, "Al_Percentage"),
  Pb = getcol(sd_use, "Pb(mg/kg)"),
  Zn = getcol(sd_use, "Zn(mg/kg)"),
  Cu = getcol(sd_use, "Cu(mg/kg)"),
  Cd = getcol(sd_use, "Cd(mg/kg)"),
  Ni = getcol(sd_use, "Ni(mg/kg)"),
  
  ## extra context (SI / robustness only)
  P  = getcol(sd_use, "P(mg/kg)"),
  Ca = getcol(sd_use, "Ca_Percentage"),
  Ti = getcol(sd_use, "Ti_Percentage"),
  
  stringsAsFactors = FALSE
)


## A.1  enrichment factors:
## EF = (Me/Al)_sample / (Me/Al)_baseline;  baseline = pre-1800 CE median of Me/Al WITHIN EACH CORE (fallback: 3 oldest samples if no pre-1800 material).
add_EF <- function(df, metals = c("Pb","Zn","Cu","Cd","Ni"), cutoff = 1800) {
  stale <- grep("_(AlRatio|base|EF)$", names(df), value = TRUE)
  if (length(stale)) df <- df[, setdiff(names(df), stale), drop = FALSE]
  
  need <- c(metals, "Al", "Core", "Year")
  miss <- setdiff(need, names(df))
  if (length(miss))
    stop("add_EF(): missing column(s): ", paste(miss, collapse = ", "))
  
  for (m in metals) {
    ratio <- df[[m]] / df$Al
    base <- tapply(seq_len(nrow(df)), df$Core, function(i) {
      pre <- ratio[i][!is.na(df$Year[i]) & df$Year[i] < cutoff]
      if (sum(is.finite(pre)) >= 1) return(stats::median(pre, na.rm = TRUE))
      o <- i[order(df$Year[i])]
      stats::median(ratio[o][seq_len(min(3, length(o)))], na.rm = TRUE)
    })
    df[[paste0(m, "_EF")]] <- ratio / base[as.character(df$Core)]
  }
  df
}



## B.  PER-ZONE STRATIGRAPHIC FIGURE (vertical: age on Y, proxies as columns)
## Each proxy is its own curve, grouped visually into bands, shared age axis,
## cores coloured, recent (>=1850 CE) departure shaded. 

proxy_levels <- c("TOC","Opal","CaCO3","CN","Ti_Ca","Si_Al","grain","Mo","V","Fe_Al","Pb_EF","Zn_EF","Ni_EF")
proxy_labs   <- c(TOC="TOC (%)", Opal="Opal (%)", CaCO3="CaCO3 (%)", CN="C/N",
                  Ti_Ca="Ti/Ca", Si_Al="Si/Al", grain="Grain (um)",
                  Mo="Mo (mg/kg)", V="V (mg/kg)", Fe_Al="Fe/Al",
                  Pb_EF="Pb EF", Zn_EF="Zn EF", Ni_EF="Ni_EF")

proxy_long <- proxy %>%
  dplyr::select(SampleName, Core, Zone, Year, dplyr::all_of(proxy_levels)) %>%
  tidyr::pivot_longer(dplyr::all_of(proxy_levels), names_to = "proxy", values_to = "value") %>%
  dplyr::mutate(proxy = factor(proxy, levels = proxy_levels, labels = proxy_labs[proxy_levels]))

plot_zone_vertical <- function(z) {
  d <- proxy_long %>% dplyr::filter(Zone == z, !is.na(value)) %>% dplyr::arrange(Core, Year)
  ggplot(d, aes(x = value, y = Year, colour = Core)) +
    geom_rect(inherit.aes = FALSE, ymin = 1850, ymax = 2010,
              xmin = -Inf, xmax = Inf, fill = "grey88", alpha = 0.04) +
    geom_path(aes(group = Core), linewidth = 0.4, alpha = 0.6) +
    geom_point(size = 1.2) +
    facet_grid(. ~ proxy, scales = "free_x") +
    scale_colour_manual(values = core_colors) +
    scale_x_continuous(n.breaks = 5, guide = guide_axis(check.overlap = TRUE)) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 6),
      labels = function(y) ifelse(y < 0, paste0(abs(y), " BCE"),
                                  ifelse(y == 0, "0", paste0(y, " CE")))) +
    labs(title = paste0(z, " Red Sea"), x = NULL, y = "Year", colour = "Core") +
    theme_bw(base_size = 10) +
    theme(strip.text.x     = element_text(size = 12),
          strip.background = element_rect(fill = "grey95", colour = NA),
          panel.grid.minor = element_blank(),
          axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y      = element_text(size = 12),
          axis.title       = element_text(size = 14),
          legend.title     = element_text(size = 14),
          legend.text      = element_text(size = 12),
          legend.position  = "bottom",
          plot.title       = element_text(size = 16, face = "bold"))
}

strat_north_v   <- plot_zone_vertical("North")
strat_central_v <- plot_zone_vertical("Central")
strat_south_v   <- plot_zone_vertical("South")

ggsave(paste0(out_dir, "Fig_strat_North.png"),   strat_north_v,   width = 16, height = 10, dpi = 300)
ggsave(paste0(out_dir, "Fig_strat_Central.png"), strat_central_v, width = 16, height = 10, dpi = 300)
ggsave(paste0(out_dir, "Fig_strat_South.png"),   strat_south_v,   width = 16, height = 10, dpi = 300)




## C.  PROXY DIAGNOSTICS to verify band assignment
## Performs the test Tribovillard et al. (2006) prescribe and Algeo & Liu (2020)
## reaffirm: is an element's variation detrital (=> not interpretable raw) or
## hydrogenous/authigenic (=> interpretable) Reassignment is allowed ONLY where
## a dual-source mechanism is documented. 

band_def <- list(
  Productivity  = c("TOC","Opal","CaCO3","CN","Sr_Ca"),
  Terrigenous   = c("Ti_Ca","Si_Al","grain"),
  Redox         = c("Mo","V","Mn","S_pct","Fe_Al"),
  Anthropogenic = c("Pb_EF","Zn_EF","Cu_EF","Cd_EF", "Ni_EF")
)

## TEST 1: detrital vs hydrogenous: r(raw element, Al) 
raw_elements <- c("Pb","Zn","Cu","Cd","Ni","Mo","V","Mn")
t1 <- data.frame(Element = raw_elements,
                 r_with_Al = sapply(raw_elements, function(e)
                   round(cor(proxy[[e]], proxy$Al, use = "pairwise.complete.obs"), 3)))
t1$Interpretation <- ifelse(abs(t1$r_with_Al) > 0.7, "detrital-dominated",
                            ifelse(abs(t1$r_with_Al) > 0.4, "mixed", "hydrogenous/authigenic"))
cat("\n=== TEST 1: r(element, Al)  (high |r| => detrital => normalise) ===\n")
print(t1, row.names = FALSE)

## TEST 2: does each proxy behave as its band claims?
all_vars <- unique(unlist(band_def)); all_vars <- all_vars[all_vars %in% names(proxy)]
cmat <- cor(proxy[, all_vars], use = "pairwise.complete.obs")
cat("\n=== TEST 2: diagnostic pairs ===\n")
chk <- function(a,b,note) if (all(c(a,b) %in% names(proxy)))
  cat(sprintf("  r(%-6s,%-6s)=%+.3f  %s\n", a, b,
              cor(proxy[[a]], proxy[[b]], use="pairwise.complete.obs"), note))
chk("Fe_Al","Mo","Fe/Al redox? (high => authigenic, put in Redox)")
chk("S_pct","Mo","S redox? (high => redox)")
chk("S_pct","CaCO3","S evaporitic? (high => gypsum)")
chk("Si_Al","Opal","Si/Al biogenic? (low => detrital quartz)")
chk("Ni","TOC","Ni OM-delivered? (Tribovillard)") 
chk("Ni_EF","Pb_EF","Ni excess tracks Pb excess? (=> anthropogenic)")
chk("Ti_Ca","CaCO3","Ti/Ca vs CaCO3 (shared-Ca anticorrelation)")


## TEST 3: within-band coherence
cat("\n=== TEST 3: within-band mean |r| ===\n")
for (b in names(band_def)) {
  v <- band_def[[b]][band_def[[b]] %in% names(proxy)]
  if (length(v) < 2) next
  cm <- cor(proxy[, v], use = "pairwise.complete.obs")
  cat(sprintf("  %-14s (%d vars) mean|r|=%.3f\n", b, length(v),
              mean(abs(cm[upper.tri(cm)]), na.rm = TRUE)))
}

## TEST 4: cross-band leakage (|r|>0.6). Algeo&Liu: Pb/Zn EF can read redox
cat("\n=== TEST 4: cross-band pairs |r|>0.6 ===\n")
lookup <- utils::stack(band_def); names(lookup) <- c("var","band")
leak <- data.frame()
for (i in seq_along(all_vars)) for (j in seq_len(i-1)) {
  a <- all_vars[i]; b <- all_vars[j]
  if (identical(as.character(lookup$band[match(a,lookup$var)]),
                as.character(lookup$band[ match(b,lookup$var)]))) next
  r <- cmat[a,b]
  if (!is.na(r) && abs(r) > 0.6)
    leak <- rbind(leak, data.frame(var1=a, var2=b, r=round(r,3)))
}
if (nrow(leak)) print(leak[order(-abs(leak$r)),], row.names=FALSE) else
  cat("  none — bands cleanly separated.\n")



## D.  VARIATION PARTITIONING BY BAND  

NPERM <- 9999            # final-run permutations (use 999 while drafting)

## D.1  prep: clr community matrix + aligned proxy function to use for both primers
prep_primer <- function(physeq, proxy) {
  otu <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu <- t(otu)          # -> samples x taxa
  otu_ps <- otu; otu_ps[otu_ps == 0] <- 0.01        # pseudocount
  clr <- t(apply(otu_ps, 1, function(r) { lr <- log(r); lr - mean(lr) }))
  
  common <- intersect(rownames(clr), proxy$SampleName)
  if (length(common) < 25)
    stop("prep_primer(): only ", length(common),
         " samples shared between this phyloseq object and the proxy table.\n",
         "  Check that sample_names() match proxy$SampleName for this primer.")
  
  clr    <- clr[common, , drop = FALSE]
  otu_ps <- otu_ps[common, , drop = FALSE]
  px     <- proxy[match(common, proxy$SampleName), , drop = FALSE]
  stopifnot(identical(rownames(clr), px$SampleName))
  cat(sprintf("  prep: %d samples shared (of %d proxy rows)\n", length(common), nrow(proxy)))
  list(clr = clr, otu_ps = otu_ps, px = px)
}


## D.2  band builder function
make_bands <- function(px, rows) {
  sub <- px[rows, , drop = FALSE]
  grab <- function(cols) {
    m <- as.matrix(sub[, cols, drop = FALSE])
    keep <- apply(m, 2, function(x) all(is.finite(x)) && stats::sd(x) > 1e-8)
    if (!any(keep)) return(NULL)
    scale(m[, keep, drop = FALSE])
  }
  list(prod  = grab(c("TOC","Opal","CaCO3","CN","Sr_Ca")),
       terr  = grab(c("Ti_Ca","Si_Al","grain")),
       redox = grab(c("Mo","V","Mn","S_pct","Fe_Al")),
       anth  = grab(c("Pb_EF","Zn_EF","Cu_EF","Cd_EF","Ni_EF")))
}


## D.3  helpers: unique-fraction test, and clean 2-table shared fraction
## unique fraction significance via partial RDA (adjR2 point estimate is taken
## from varpart itself, below, so the table matches the varpart object exactly)
test_unique <- function(Y, band, others, nperm = NPERM) {
  if (is.null(band)) return(c(F = NA, p = NA))
  W <- do.call(cbind, others)
  m <- vegan::rda(Y, band, W)
  a <- anova(m, permutations = nperm)
  c(F = a$F[1], p = a$`Pr(>F)`[1])
}

## D.4  one partition -> tidy results data.frame

run_partition <- function(dat, rows, label, primer,
                          do_plot = FALSE, file = NULL, nperm = NPERM) {
  n <- sum(rows)
  cat(sprintf("\n---- %s | %s  (n = %d) ----\n", primer, label, n))
  if (n < 25) { cat("  n < 25 — skipped.\n"); return(NULL) }
  if (n < 55) cat("  NOTE: low n/p; interpret DIRECTION, not exact values.\n")
  
  clr <- dat$clr; otu_ps <- dat$otu_ps; px <- dat$px
  present <- colSums(otu_ps[rows, , drop = FALSE] > 0.01) > 0   # raw-presence filter
  Y <- clr[rows, present, drop = FALSE]
  b <- make_bands(px, rows)
  if (any(vapply(b, is.null, logical(1))))
    cat("  WARNING: a band was constant in this subset.\n")
  
  vp <- vegan::varpart(Y, b$prod, b$terr, b$redox, b$anth)
  
  ## point estimates: the 4 unique fractions straight from varpart
  uniq <- vp$part$indfract$Adj.R.square[1:4]
  dfb  <- vp$part$indfract$Df[1:4]
  tot  <- vp$part$fract$Adj.R.square[nrow(vp$part$fract)]  # all-bands adj R2
  
  ## significance of each unique fraction
  tp <- rbind(
    test_unique(Y, b$prod,  list(b$terr, b$redox, b$anth), nperm),
    test_unique(Y, b$terr,  list(b$prod, b$redox, b$anth), nperm),
    test_unique(Y, b$redox, list(b$prod, b$terr,  b$anth), nperm),
    test_unique(Y, b$anth,  list(b$prod, b$terr,  b$redox), nperm))
  
  ## whole-model global test
  gm  <- anova(vegan::rda(Y, cbind(b$prod, b$terr, b$redox, b$anth)), permutations = nperm)
  
  res <- data.frame(
    Primer   = primer, Window = label, n = n,
    Fraction = c("Productivity","Terrigenous","Redox","Anthropogenic","Whole model"),
    adjR2    = round(c(uniq, sh_ra, tot), 4),
    Df       = c(dfb, NA, gm$Df[1]),
    F        = round(c(tp[,"F"], NA, gm$F[1]), 3),
    p        = c(tp[,"p"], NA, gm$`Pr(>F)`[1]),
    stringsAsFactors = FALSE)
  print(res, row.names = FALSE)
  
  if (do_plot && !is.null(file)) {
    png(file, width = 1600, height = 1600, res = 220)
    plot(vp, bg = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"),
         Xnames = c("Prod","Terr","Redox","Anth")); title(paste(primer, label)); dev.off()
  }
  res
}

## D.5  run both primers x {whole record, pre-1800, post-1800}

primers <- list("18S v7" = phyeufinal, "COI" = phycoifinal)

all_results <- list()
for (pr in names(primers)) {
  dat <- prep_primer(primers[[pr]], proxy)
  wr  <- rep(TRUE, nrow(dat$px))
  pre <- dat$px$Year <  1800
  pos <- dat$px$Year >= 1800
  all_results[[paste(pr,"whole")]] <- run_partition(dat, wr,  "Whole record",  pr)
  all_results[[paste(pr,"pre")]]   <- run_partition(dat, pre, "Pre-1800 CE",   pr)
  all_results[[paste(pr,"post")]]  <- run_partition(dat, pos, "Post-1800 CE",  pr)
}


## D.6  assemble the publication table (CSV; drop into manuscript / SI)
varpart_table <- do.call(rbind, all_results)
rownames(varpart_table) <- NULL
cat("\n\n#### COMBINED VARPART TABLE ####\n")
print(varpart_table, row.names = FALSE)
write.csv(varpart_table, "varpart_results_both_primers.csv", row.names = FALSE)

