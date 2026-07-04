CASE — Reproducible Analysis (Coastal Acidification & Sewage Effluent)
================
Jonathan Puritz
2026-07-04

- [Setup](#setup)
  - [Software environment](#software-environment)
- [Fig 1 — Experimental design and replicated
  survival](#fig-1--experimental-design-and-replicated-survival)
  - [Survival statistics](#survival-statistics)
  - [Fig 1 — four-spawn survival
    strip](#fig-1--four-spawn-survival-strip)
  - [Fig 1 — assemble with PowerPoint
    schematic](#fig-1--assemble-with-powerpoint-schematic)
    - [Fig S1 — Final survival across T0
      denominators](#fig-s1--final-survival-across-t0-denominators)
    - [Fig S2 — Larval morphology and development
      (4-panel)](#fig-s2--larval-morphology-and-development-4-panel)
  - [Filtering](#filtering)
    - [Create Unified Sync file](#create-unified-sync-file)
    - [Add coverage stats to sync file and filter by minimum
      coverage](#add-coverage-stats-to-sync-file-and-filter-by-minimum-coverage)
  - [Background allele frequencies (10,000 random
    loci)](#background-allele-frequencies-10000-random-loci)
    - [Venn Diagram for called loci](#venn-diagram-for-called-loci)
  - [Differentiate Singleton Loci](#differentiate-singleton-loci)
- [Individual Blocks for outlier
  testing](#individual-blocks-for-outlier-testing)
  - [Calculate P-values](#calculate-p-values)
  - [Across all blocks](#across-all-blocks)
  - [Determine Significant Loci](#determine-significant-loci)
    - [Proximity Control Loci filter](#proximity-control-loci-filter)
    - [Tier-specific CON filtering and treatment
      splits](#tier-specific-con-filtering-and-treatment-splits)
    - [Tier-specific gene annotation](#tier-specific-gene-annotation)
    - [Outlier-overlap LOC lists (for the Fig 4
      Venn)](#outlier-overlap-loc-lists-for-the-fig-4-venn)
  - [Gene Ontology Enrichment Analysis with
    topGO](#gene-ontology-enrichment-analysis-with-topgo)
    - [Download and prepare GO
      annotation](#download-and-prepare-go-annotation)
    - [Create gene-to-GO mapping file](#create-gene-to-go-mapping-file)
    - [Extract GO terms for candidate genes (by stressor
      category)](#extract-go-terms-for-candidate-genes-by-stressor-category)
    - [Extract GO terms for tier-specific candidate
      genes](#extract-go-terms-for-tier-specific-candidate-genes)
    - [Load libraries and gene-to-GO
      data](#load-libraries-and-gene-to-go-data)
    - [Define GO enrichment function](#define-go-enrichment-function)
    - [Run GO enrichment for all four
      categories](#run-go-enrichment-for-all-four-categories)
    - [Run GO enrichment for Core
      loci](#run-go-enrichment-for-core-loci)
    - [Run GO enrichment for Private
      loci](#run-go-enrichment-for-private-loci)
    - [Run GO enrichment for Convergent
      loci](#run-go-enrichment-for-convergent-loci)
    - [Save results](#save-results)
    - [Save tier-specific topGO
      results](#save-tier-specific-topgo-results)
  - [Data-driven GO clusters (mechanism themes for Fig 4D + Table
    S)](#data-driven-go-clusters-mechanism-themes-for-fig-4d--table-s)
  - [Outlier allele-frequency matrices (for the Fig 2
    PCAs)](#outlier-allele-frequency-matrices-for-the-fig-2-pcas)
- [Fig 2 — Outlier loci & cross-spawn
  repeatability](#fig-2--outlier-loci--cross-spawn-repeatability)
  - [Figure parameters and guard](#figure-parameters-and-guard)
  - [2.1 PCAs — ΔAF-from-initial (all blocks) + three per-block outlier
    PCAs](#21-pcas--δaf-from-initial-all-blocks--three-per-block-outlier-pcas)
  - [2.2 Reproducibility bars (3b) + lethality coupling
    (3c)](#22-reproducibility-bars-3b--lethality-coupling-3c)
  - [2.3 Within-gene allele turnover (3E-A) — same genes, different
    SNPs](#23-within-gene-allele-turnover-3e-a--same-genes-different-snps)
  - [Fig S3 — Genome-wide PCA (10,000 random
    loci)](#fig-s3--genome-wide-pca-10000-random-loci)
  - [Fig S4 — Tier-specific PCA](#fig-s4--tier-specific-pca)
  - [Tier-specific PCA](#tier-specific-pca)
    - [Extract tier-specific sync
      files](#extract-tier-specific-sync-files)
    - [All-samples PCA helper
      function](#all-samples-pca-helper-function)
    - [Run all-samples PCA for each
      tier](#run-all-samples-pca-for-each-tier)
- [Fig 3 — Magnitude synergy: control-anchored
  additivity](#fig-3--magnitude-synergy-control-anchored-additivity)
  - [3.1 Additivity scatter and per-panel orthogonal
    fits](#31-additivity-scatter-and-per-panel-orthogonal-fits)
  - [Fig S5 — Additivity scatter faceted by replication
    tier](#fig-s5--additivity-scatter-faceted-by-replication-tier)
- [Fig 4 — Compositional synergy (novel
  targets)](#fig-4--compositional-synergy-novel-targets)
  - [4.1 Block-12 Manhattan, stacked across CA / SE /
    CASE](#41-block-12-manhattan-stacked-across-ca--se--case)
  - [4.1b Per-spawn Manhattan stacks (B10, B11) —
    supplement](#41b-per-spawn-manhattan-stacks-b10-b11--supplement)
  - [4.2 GO enrichment — custom plot from the TopGO
    CSVs](#42-go-enrichment--custom-plot-from-the-topgo-csvs)
  - [4.3 One Venn — SNP and gene counts together (custom
    ggplot)](#43-one-venn--snp-and-gene-counts-together-custom-ggplot)
  - [4.4b Mechanism panel D — data-driven stress GO
    clusters](#44b-mechanism-panel-d--data-driven-stress-go-clusters)
  - [4.5 Assemble Fig 4 (patchwork)](#45-assemble-fig-4-patchwork)
  - [Fig S6 — Per-spawn outlier overlap
    (Venns)](#fig-s6--per-spawn-outlier-overlap-venns)
  - [Fig S9 — Full GO enrichment by stressor
    (topGO)](#fig-s9--full-go-enrichment-by-stressor-topgo)
  - [Fig S13 — GO-term overlap and biological-theme
    enrichment](#fig-s13--go-term-overlap-and-biological-theme-enrichment)
- [Fig 5 — Diversity, architecture, and convergent function (CA / SE /
  CASE)](#fig-5--diversity-architecture-and-convergent-function-ca--se--case)
  - [Fig S15 — Top enriched GO terms per tier
    (topGO)](#fig-s15--top-enriched-go-terms-per-tier-topgo)
  - [Fig S16 — Tier architecture by spawn
    (3x3)](#fig-s16--tier-architecture-by-spawn-3x3)
- [Supplementary tables](#supplementary-tables)
  - [Table S3 — GO enrichment terms by stressor (one row per
    term)](#table-s3--go-enrichment-terms-by-stressor-one-row-per-term)
    - [Table S3a — All Outlier Loci](#table-s3a--all-outlier-loci)
    - [Table S3b — CA Stressor](#table-s3b--ca-stressor)
    - [Table S3c — SE Stressor](#table-s3c--se-stressor)
    - [Table S3d — CASE Stressor
      (CA+SE)](#table-s3d--case-stressor-case)
  - [Table 1 — All data-driven GO clusters (main
    text)](#table-1--all-data-driven-go-clusters-main-text)

> Single reproducible pipeline for the CASE study. The document is
> organized around the five main figures (Fig 1 survival/phenotype → Fig
> 2 outliers & declining reproducibility → Fig 3 magnitude synergy:
> control-anchored additivity → Fig 4 compositional synergy: novel
> targets → Fig 5 diversity, architecture & convergent function); a
> shared genomic-substrate section precedes the figures because Figs 2–5
> all depend on the same allele-frequency, CMH and locus-tier objects.
> Supplementary figures (S1–S15) and tables follow each result. Heavy
> preprocessing steps that require external tools or large intermediate
> files are kept as `eval=FALSE` documented commands.
>
> Outputs: main figures → `figures/main/`, supplements →
> `figures/supp/`, tables → `results/tables/`.

# Setup

``` r
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(qvalue)
library(stringr)
library(pcadapt)
library(poolSeq)
library(ACER)
library(scales)
library(data.table)
library(patchwork)
library(ggridges)
library(ggrain)
library(grid)
library(ggman)
library(VennDiagram)
library(topGO)
library(kableExtra)

# ---- output directories -----------------------------------------------------
out_dir  <- "figures/main"
supp_dir <- "figures/supp"
tab_dir  <- "results/tables"
for (d in c(out_dir, supp_dir, tab_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

save_case <- function(plot, file, width = 10, height = 6, res = 300) {
  png(file.path(out_dir, file), type = "cairo", units = "in",
      width = width, height = height, res = res, bg = "transparent",
      family = "DejaVu Sans")
  print(plot); invisible(dev.off())
}
save_supp <- function(plot, file, width = 10, height = 6, res = 300) {
  png(file.path(supp_dir, file), type = "cairo", units = "in",
      width = width, height = height, res = res, bg = "transparent",
      family = "DejaVu Sans")
  print(plot); invisible(dev.off())
}

# ---- palettes (Okabe-Ito; CON/CA/SE/CASE order) -----------------------------
treat_levels <- c("CON","CA","SE","CASE")
treat_cols   <- c(CON="#999999", CA="#E69F00", SE="#0072B2", CASE="#009E73", IS="sky blue")
cbPalette        <- c("#0072B2","#56B4E9","#D55E00","#E69F00","#009E73","#999999","#F0E442","#CC79A7")
cbPaletteSmall   <- c("#E69F00","#009E73","#999999","sky blue","#0072B2")
cbPaletteSmall3  <- c("#E69F00","#0072B2","#009E73")        # CA, SE, CASE (Venn fills)
cbPaletteSmall4  <- c("#999999","#E69F00","#0072B2","#009E73")
cbPaletteSmall5  <- c("#999999","#E69F00","#0072B2","#009E73")
spawn_levels <- c("B10","B11","B12")
# Spawn palette now reuses the Tier palette (Core/Convergent/Private colours),
# repurposed 1:1 for B10/B11/B12.
spawn_cols   <- c(B10="#CC79A7", B11="#D55E00", B12="#F0E442")
spawn_shapes <- c(B10=21, B11=22, B12=23)
stat_cols    <- c(Observed="#D7191C", `CON baseline`="#999999", `Random null`="#4D4D4D")
tier_levels  <- c("Core","Convergent","Private","Background")
tier_cols    <- c(Core="#0072B2", Convergent="#E69F00", `Private`="#CC79A7", Background="grey80")
tier_shape   <- c(Core=21, Convergent=23, `Private`=25)

theme_case <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right",
          legend.title = element_text(face = "bold", size = base_size + 2),
          legend.text  = element_text(size = base_size * 0.8 + 2))
}

# ---- shared helper functions ------------------------------------------------
# Convert per-treatment CMH p-values to q-values (pi0 = 1 => Benjamini-Hochberg)
# and relabel chromosomes 1-10 from their RefSeq accessions.
Qvalue_convert <- function(table) {
  table <- spread(table, GROUP, PVAL)
  table$CHR <- table$CHROM
  table %>%
    mutate(CHR = str_replace(CHR, "NC_035780.1", "1"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035781.1", "2"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035782.1", "3"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035783.1", "4"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035784.1", "5"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035785.1", "6"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035786.1", "7"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035787.1", "8"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035788.1", "9"))  %>%
    mutate(CHR = str_replace(CHR, "NC_035789.1", "10")) -> pp1
  pp1 <- unite(pp1, "SNP", c("CHROM","BP"), remove = FALSE)
  pp1$CHR   <- as.numeric(pp1$CHR)
  pp1$QCA   <- qvalue(pp1$PCA,   pi0 = 1)$qvalues
  pp1$QCON  <- qvalue(pp1$PCON,  pi0 = 1)$qvalues
  pp1$QSE   <- qvalue(pp1$PSE,   pi0 = 1)$qvalues
  pp1$QCASE <- qvalue(pp1$PCASE, pi0 = 1)$qvalues
  return(pp1)
}

# Flag and subset treatment-significant loci (excluding CON-significant ones).
Significat_subset <- function(pv, alpha, alpha2) {
  pv$Sig.CASE <- pv$QCASE < alpha
  pv$Sig.CA   <- pv$QCA   < alpha
  pv$Sig.SE   <- pv$QSE   < alpha
  ppsig <- subset(pv, QCA < alpha | QCASE < alpha | QSE < alpha)
  ppsig <- subset(ppsig, QCON > alpha2)
  print(nrow(pv)); print(nrow(ppsig)); print(nrow(ppsig) / nrow(pv))
  return(ppsig)
}

# Collapse a per-block table to one row per SNP (minimum q across blocks).
group_and_average <- function(df) {
  df %>%
    dplyr::group_by(SNP) %>%
    dplyr::summarize(across(c(PCA, PCASE, PCON, PSE, CHR, QCA, QCON, QSE, QCASE),
                            ~ if (all(is.na(.))) NA else min(., na.rm = TRUE)),
                     Sig.CASE = any(Sig.CASE, na.rm = TRUE),
                     Sig.CA = any(Sig.CA, na.rm = TRUE), Sig.SE = any(Sig.SE, na.rm = TRUE),
                     CHROM = dplyr::first(CHROM), BP = dplyr::first(BP), .groups = "drop")
}
```

## Software environment

``` bash
# Popoolation2 helper scripts and conda environments (run once).
cd ../scripts && git clone https://github.com/ToBoDev/assessPool.git
mamba env create --file ../CASE_environment.yaml
mamba env create --file ../random_draw_environment.yaml
```

# Fig 1 — Experimental design and replicated survival

## Survival statistics

``` r
mort_file <- "../data/CASE_FINAL_Mortality.txt"
mort_all  <- fread(mort_file)                                  # keeps original (spaced/%) col names
# The file gives survival relative to four T0 denominators; rename for convenience.
setnames(mort_all, "%_Initial Calculated", "surv_Calculated")
setnames(mort_all, "%_Initial_Composite",  "surv_Composite")
setnames(mort_all, "%_Initial_Average",    "surv_Average")
setnames(mort_all, "%_Initial_Seeded",     "surv_Seeded")
mort_all[, Treatment := factor(Treatment, levels = c("CON","CA","SE","CASE"))]
final_data_end <- mort_all[Time > 12]                          # end timepoint (Time 24; all four spawns incl. B6)

# PRIMARY survival for ALL statistics = the COMPOSITE denominator. Initial density is estimated
# with sampling error (under-homogenization of near-full containers), so the composite T0 = mean of
# the stocking (seeded) density, the jar's own calculated T0 density, and the block-mean calculated
# T0 density — the best accuracy/error trade-off. The Calculated/Average/Seeded variants are kept
# only for the denominator-robustness display in Fig 1.
final_data_end[, Survival := surv_Composite]
final_data_end[, Survival_centered := Survival - mean(Survival, na.rm = TRUE), by = Block]

print(final_data_end[, .(mean_centered = mean(Survival_centered, na.rm = TRUE)), by = Treatment])
```

    ##    Treatment mean_centered
    ##       <fctr>         <num>
    ## 1:      CASE     -7.156910
    ## 2:       CON      8.164350
    ## 3:        SE     -6.565336
    ## 4:        CA      5.557896

``` r
fit <- aov(Survival_centered ~ Treatment, data = final_data_end)
print(summary(fit)); print(TukeyHSD(fit, "Treatment"))
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Treatment    3   2302   767.5   7.392 0.000409 ***
    ## Residuals   44   4568   103.8                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = Survival_centered ~ Treatment, data = final_data_end)
    ## 
    ## $Treatment
    ##                 diff       lwr       upr     p adj
    ## CA-CON    -2.6064539 -13.71311  8.500203 0.9229719
    ## SE-CON   -14.7296860 -25.83634 -3.623029 0.0051062
    ## CASE-CON -15.3212595 -26.42792 -4.214602 0.0033899
    ## SE-CA    -12.1232322 -23.22989 -1.016575 0.0275706
    ## CASE-CA  -12.7148056 -23.82146 -1.608148 0.0191620
    ## CASE-SE   -0.5915734 -11.69823 10.515084 0.9989579

``` r
# --- Factorial CA x SE GLMM on the composite survival (magnitude synergy at the phenotype) ---
# Betabinomial with the composite T0 count as the binomial denominator: survivors_i = round(p_i*N0_i),
# N0_i = composite T0 larvae recovered per jar as 100 * Calculated(end) / composite-survival%. The
# composite T0 count is a density estimate, so this weights jars by their estimated initial size; an
# ordered-beta cross-check on the proportion (correct family for continuous %) is printed alongside.
suppressMessages(library(glmmTMB))
sd_glmm <- copy(final_data_end)
sd_glmm[, CA := factor(fifelse(Treatment %in% c("CA","CASE"), "1", "0"))]
sd_glmm[, SE := factor(fifelse(Treatment %in% c("SE","CASE"), "1", "0"))]
sd_glmm[, Block := factor(Block)]
sd_glmm[, p := pmin(pmax(Survival / 100, 1e-4), 1 - 1e-4)]
sd_glmm[, N0 := pmax(round(`Calculated Larvae` * 100 / Survival), 1L)]   # composite T0 count
sd_glmm[, succ := pmin(round(p * N0), N0)][, fail := N0 - succ]
m_bb <- glmmTMB(cbind(succ, fail) ~ CA * SE + (1 | Block), family = betabinomial, data = sd_glmm)
cat("\n--- Betabinomial GLMM (composite survival): CA x SE interaction ---\n")
```

    ## 
    ## --- Betabinomial GLMM (composite survival): CA x SE interaction ---

``` r
print(summary(m_bb)$coefficients$cond)
```

    ##                Estimate Std. Error    z value     Pr(>|z|)
    ## (Intercept) -0.57472411  0.3923442 -1.4648468 0.1429627166
    ## CA1         -0.07035818  0.1812099 -0.3882690 0.6978169887
    ## SE1         -0.77943577  0.2078715 -3.7496031 0.0001771147
    ## CA1:SE1      0.04243442  0.2909081  0.1458688 0.8840249552

``` r
surv_glmm_tab <- as.data.table(summary(m_bb)$coefficients$cond, keep.rownames = "term")
fwrite(surv_glmm_tab, file.path(tab_dir, "Fig1_survival_GLMM.csv"))
cat(sprintf("\nInteraction CA1:SE1: beta = %.3f, p = %.3g  (n.s. => phenotype is not synergistic)\n",
            surv_glmm_tab[term == "CA1:SE1", Estimate], surv_glmm_tab[term == "CA1:SE1", `Pr(>|z|)`]))
```

    ## 
    ## Interaction CA1:SE1: beta = 0.042, p = 0.884  (n.s. => phenotype is not synergistic)

``` r
m_ob <- tryCatch(glmmTMB(p ~ CA * SE + (1 | Block), family = ordbeta, data = sd_glmm),
                 error = function(e) NULL)
if (!is.null(m_ob)) { cat("\n--- ordered-beta cross-check (proportion) ---\n"); print(summary(m_ob)$coefficients$cond) }
```

## Fig 1 — four-spawn survival strip

The PowerPoint experiment schematic is placed above this panel in the
final figure.

``` r
# Composite-denominator survival at the end timepoint (the value used for all statistics). The
# sensitivity to the T0 denominator choice is shown separately in Fig S3.
mort_end2 <- final_data_end[, .(Block, Treatment,
  Spawn = factor(paste0("B", Block), levels = paste0("B", sort(unique(Block)))),
  Survival = surv_Composite)]
fig1b <- ggplot(mort_end2, aes(Treatment, Survival, fill = Treatment)) +
  geom_boxplot(notch = FALSE, outlier.colour = NA, alpha = 0.75, colour = "black", width = 0.25) +
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, colour = "black") +
  geom_jitter(aes(colour = Treatment), width = 0.08, height = 0, size = 2, alpha = 0.8) +
  scale_fill_manual(values = treat_cols) + scale_colour_manual(values = treat_cols) +
  facet_wrap(~ Spawn, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Percent survival (composite T0)") +
  theme_case() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "none")
fig1b
```

![](Final_reproducible_analysis_files/figure-gfm/fig1-survival-1.png)<!-- -->

## Fig 1 — assemble with PowerPoint schematic

``` r
fig1_top_png <- png::readPNG("../data/Fig1Panel1_top.png")
fig1_top <- patchwork::wrap_elements(full = grid::rasterGrob(fig1_top_png, interpolate = TRUE))

fig1 <- fig1_top / fig1b +
  plot_layout(heights = c(1, 0.6)) 
  #plot_annotation(tag_levels = "A")
fig1
```

![](Final_reproducible_analysis_files/figure-gfm/fig1-assemble-1.png)<!-- -->

``` r
save_case(fig1, "Fig1_experiment_survival.png", width = 14, height = 9)
```

### Fig S1 — Final survival across T0 denominators

Sensitivity of the mortality estimate to the initial-density (T0)
reference. The main figure (Fig 1) uses the composite denominator; here
the end-timepoint survival is shown under all four denominators, faceted
by spawn (rows) and denominator (columns), so the treatment order can be
seen to hold regardless of the reference.

``` r
mort_dn <- melt(final_data_end[, .(Block,
    Spawn = factor(paste0("B", Block), levels = paste0("B", sort(unique(Block)))), Treatment,
    Composite = surv_Composite, Calculated = surv_Calculated, Average = surv_Average, Seeded = surv_Seeded)],
  id.vars = c("Block","Spawn","Treatment"), variable.name = "Denominator", value.name = "Survival")
mort_dn[, Denominator := factor(Denominator, levels = c("Composite","Calculated","Average","Seeded"))]
figS1mort <- ggplot(mort_dn, aes(Treatment, Survival, fill = Treatment)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.7, colour = "black", width = 0.55) +
  geom_jitter(aes(colour = Treatment), width = 0.12, height = 0, size = 1.1, alpha = 0.7) +
  scale_fill_manual(values = treat_cols) + scale_colour_manual(values = treat_cols, guide = "none") +
  facet_grid(Spawn ~ Denominator, scales = "free_y") +
  labs(x = NULL, y = "Percent survival",
       title = "Fig S1 — final survival by treatment across T0 denominators and spawns") +
  theme_case(11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "none")
save_supp(figS1mort, "FigS1_mortality_denominators.png", width = 12, height = 11)
figS1mort
```

![](Final_reproducible_analysis_files/figure-gfm/figS1-mortality-denominators-1.png)<!-- -->
\## Larval development and size

``` r
dev_size <- read.table("../data/CASE_LNDA.txt", sep = "\t", header = TRUE)
dev_size <- filter(dev_size, Time != "T04")
dev_size <- filter(dev_size, N.D.A != "NA")
dev_size$Treatment <- factor(dev_size$Treatment, levels = c("CON","CA","SE","CASE"))
dev_size$N.D.A     <- factor(dev_size$N.D.A, levels = c("Dead","Deformed","Normal"))
NDA.plot <- ggplot(dev_size, aes(x = Treatment, fill = N.D.A)) +
  geom_bar(position = "fill", na.rm = TRUE) +
  facet_grid(rows = vars(Block), cols = vars(Time), scales = "free") +
  labs(y = "Proportion", title = "Larval development (normal / deformed / dead)") +
  theme_classic()
```

``` r
final_size <- read.table("../data/CASE_LARVAL_SIZE.txt", sep = "\t", header = TRUE)
final_size <- filter(final_size, Time != "T04")
final_size <- filter(final_size, Larval.Stage != "Embryo")
final_size <- filter(final_size, Normal.Deformed == "Normal")
final_size <- filter(final_size, Alive.Dead == "Alive")
final_size$Treatment <- factor(final_size$Treatment, levels = c("CON","CA","SE","CASE"))

ShellHeight <- ggplot(final_size, aes(x = Time, y = Height, fill = Treatment)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(rows = vars(Larval.Stage), cols = vars(Block), scales = "free") +
  scale_fill_manual(values = cbPaletteSmall5) +
  labs(title = "Shell height") + theme_classic()

HingeLength <- ggplot(final_size, aes(x = Time, y = Length, fill = Treatment)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(rows = vars(Larval.Stage), cols = vars(Block), scales = "free") +
  scale_fill_manual(values = cbPaletteSmall5) +
  labs(title = "Hinge length") + theme_classic()

Circumference <- ggplot(final_size, aes(x = Time, y = Circumference, fill = Treatment)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(rows = vars(Larval.Stage), cols = vars(Block), scales = "free") +
  scale_fill_manual(values = cbPaletteSmall5) +
  labs(title = "Circumference") + theme_classic()
```

### Fig S2 — Larval morphology and development (4-panel)

``` r
figS2 <- (ShellHeight + HingeLength) / (Circumference + NDA.plot) +
  plot_annotation(tag_levels = "A",
                  title = "Fig S1 — Larval morphology and development")
save_supp(figS2, "FigS2_morphology.png", width = 16, height = 13)
figS2
```

![](Final_reproducible_analysis_files/figure-gfm/figS2-morphology-1.png)<!-- -->
\# Genomic data — shared substrate

These steps build the allele-frequency, CMH q-value, and locus-tier
objects that all of Figs 2–4 depend on. Steps that need large
intermediate files or external tools (bcftools, popoolation2, bedtools,
dDocent) are kept as `eval=FALSE`.

## Filtering

``` bash
source activate CASE
cd ../raw.vcf
bash ../scripts/CASE_PVCF_filter2.sh CASE 64 20
```

``` bash
source activate CASE
cd ../raw.vcf 
bcftools index CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz
bcftools query -l CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz > samples
```

### Create Unified Sync file

Each sequencing run had some extra samples that need to be removed
before starting the analysis Each block was filtered for missing data
less than 10% and MAF of \> 0.015 based on read counts using VCF file
first.

``` bash
source activate CASE

cd ../raw.vcf
bcftools view -S <(grep B11 samples | grep -v 1.G) CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz | bcftools view -i 'F_MISSING<0.1'  | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' | bcftools query -f '%CHROM\t%POS\n' > B11.pos

bcftools view -S <(grep B10 samples | grep -v J17 | grep -v J05) CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz | bcftools view -i 'F_MISSING<0.1' | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' | bcftools query -f '%CHROM\t%POS\n' > B10.pos

bcftools view -S <(grep B12 samples) CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz | bcftools view -i 'F_MISSING<0.1' | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' | bcftools query -f '%CHROM\t%POS\n'  > B12.pos

cat B*.pos | sort | uniq > fil.pos

bcftools view -R fil.pos -S <(grep -v 1.GB11 samples | grep -v J17B10| grep -v J05B10) --threads 40 -m2 -M 4 CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz -O z -o SNP.CASE.TRSdp.20.B90.2a.perp.vcf.gz
```

#### Create Sync files

``` bash
source activate CASE

bcftools view --threads 40 ../raw.vcf/SNP.CASE.TRSdp.20.B90.2a.perp.vcf.gz  | mawk -F'\t' -v OFS='\t' '{ for(i=1;i<=NF;i++) if($i==".:.:.:.:.:.:.:.") $i="."; print }' > temp.vcf
python2 ../scripts/VCFtoPopPool.py temp.vcf CASE.All.Blocks.sync 
rm temp.vcf
```

#### This sort is to maintain reproducibility across older analysis

``` bash
mawk 'BEGIN {OFS="\t"
    order["CHROM"] = 0
    order["NC_035789.1"] = 1
    order["NC_035780.1"] = 2
    order["NC_035781.1"] = 3
    order["NC_035782.1"] = 4
    order["NC_035783.1"] = 5
    order["NC_035784.1"] = 6
    order["NC_035785.1"] = 7
    order["NC_035786.1"] = 8
    order["NC_035787.1"] = 9
    order["NC_035788.1"] = 10
}
{ print order[$1], $2, NR, $0 }' CASE.All.Blocks.sync | sort -k1,1n -k2,2n -k3,3n | cut -f4- > test.sync

mv test.sync CASE.All.Blocks.sync
```

``` bash
cat <(echo -e "CHROM\tPOS") ../raw.vcf/B10.pos > B10.pos
cat <(echo -e "CHROM\tPOS") ../raw.vcf/B11.pos > B11.pos
cat <(echo -e "CHROM\tPOS") ../raw.vcf/B12.pos > B12.pos

mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' B10.pos CASE.All.Blocks.sync | mawk 'NR==1{for(i=1;i<=NF;i++)if(i<=3||$i~/B10$/)k[i]}{o="";for(i=1;i<=NF;i++)if(i in k)o=o?o OFS $i:$i;print o}' OFS='\t' | mawk '!/\t\.\t/' > CASE.Block10.sync
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' B11.pos CASE.All.Blocks.sync | mawk 'NR==1{for(i=1;i<=NF;i++)if(i<=3||$i~/B11$/)k[i]}{o="";for(i=1;i<=NF;i++)if(i in k)o=o?o OFS $i:$i;print o}' OFS='\t' | mawk '!/\t\.\t/' > CASE.Block11.sync
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' B12.pos CASE.All.Blocks.sync | mawk 'NR==1{for(i=1;i<=NF;i++)if(i<=3||$i~/B12$/)k[i]}{o="";for(i=1;i<=NF;i++)if(i in k)o=o?o OFS $i:$i;print o}' OFS='\t' | mawk '!/\t\.\t/' > CASE.Block12.sync
```

### Add coverage stats to sync file and filter by minimum coverage

``` bash
source activate CASE
mawk -f ../scripts/add_cov_sync CASE.Block10.sync | mawk '$20 > 9 && $22 > 24' > CASE.dp20.Block10.cov.sync &
mawk -f ../scripts/add_cov_sync CASE.Block11.sync | mawk '$24 > 9 && $26 > 24' > CASE.dp20.Block11.cov.sync &
mawk -f ../scripts/add_cov_sync CASE.Block12.sync | mawk '$24 > 9 && $26 > 24' > CASE.dp20.Block12.cov.sync &
#mawk -f ../scripts/add_cov_sync CASE.All.Blocks.sync | mawk '$66 > 9 && $68 > 24' > CASE.All.Blocks.cov.sync
mawk -f ../scripts/add_cov_sync CASE.All.Blocks.sync | mawk '$60 > 9 && $62 > 24' | mawk '!/\t\.\t/'  > CASE.All.Blocks.cov.sync

cut -f1,2 CASE.dp20.Block1*.cov.sync | sort | uniq -c | mawk '$1 > 2' | mawk '{print $2 "\t" $3}' > filtered.loci.allthreeblocks

cat <(echo "SNP") <(cat Block1*.pos | cut -f1,2 | sort | uniq -c | mawk '$1 <2' | mawk '{print $2"_"$3}') > singleton.loci
```

## Background allele frequencies (10,000 random loci)

Used for the genome-wide PCA (Fig S3).

``` bash
source activate CASE

mawk '!/CHR/' filtered.loci.allthreeblocks | shuf -n 50000  | shuf -n 10000 > rand.loci

mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' rand.loci CASE.All.Blocks.cov.sync | cut --complement -f60- > input.all.rand.sync
```

``` r
input.all.rand <- read.sync(file = "input.all.rand.sync", gen = seq(1,56), repl = seq(1,56), polarization = "rising")
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
input.all.rand.af <- af(input.all.rand, gen = seq(1,56), repl = seq(1,56))
treatment.labels <- c(rep("CA",11), rep("CASE",11), rep("CON",11), rep("IS",12), rep("SE",11))
spawn.names <- c("B12","B11","B12","B10","B11","B10","B11","B12","B10","B12","B11","B10","B11","B12","B12","B10","B12","B11","B10","B11","B11","B12","B10","B11","B12","B11","B10","B12","B10","B11","B12","B11","B12","B10","B11","B10","B11","B10","B11","B10","B11","B12","B12","B12","B12","B11","B12","B10","B10","B11","B11","B10","B12","B12","B11","B12")
colnames(input.all.rand.af) <- paste(treatment.labels, spawn.names, sep = "_")
```

### Venn Diagram for called loci

``` r
# Read and preprocess data
B10.loci <- read.table("Block10.pos", header = TRUE)
colnames(B10.loci) <- c("CHROM", "BP")
B10.loci <- B10.loci %>% mutate(BLOCK_10 = TRUE)

B11.loci <- read.table("Block11.pos", header = TRUE)
colnames(B11.loci) <- c("CHROM", "BP")
B11.loci <- B11.loci %>% mutate(BLOCK_11 = TRUE)

B12.loci <- read.table("Block12.pos", header = TRUE)
colnames(B12.loci) <- c("CHROM", "BP")
B12.loci <- B12.loci %>% mutate(BLOCK_12 = TRUE)

# Full join all data frames to preserve all CHROM and BP combinations
pos_df <- B10.loci %>%
  full_join(B11.loci %>% dplyr::select(CHROM, BP, BLOCK_11), by = c("CHROM", "BP")) %>%
  full_join(B12.loci %>% dplyr::select(CHROM, BP, BLOCK_12), by = c("CHROM", "BP"))

# Replace NA values in the BLOCK_* columns with FALSE
pos_df <- pos_df %>%
  mutate(
    BLOCK_10 = ifelse(is.na(BLOCK_10), FALSE, BLOCK_10),
    BLOCK_11 = ifelse(is.na(BLOCK_11), FALSE, BLOCK_11),
    BLOCK_12 = ifelse(is.na(BLOCK_12), FALSE, BLOCK_12)
  )


loci.b10 <- nrow(subset(pos_df,BLOCK_10 == TRUE))
loci.b11 <- nrow(subset(pos_df,BLOCK_11 == TRUE))
loci.b12 <- nrow(subset(pos_df,BLOCK_12 == TRUE))

loci.b10.b11 <- nrow(subset(pos_df,BLOCK_10 == TRUE & BLOCK_11 ==TRUE))
loci.b10.b12 <- nrow(subset(pos_df,BLOCK_10 == TRUE & BLOCK_12 ==TRUE))
loci.b11.b12 <- nrow(subset(pos_df,BLOCK_11 == TRUE & BLOCK_12 ==TRUE))
loci.b10.b11.b12 <- nrow(subset(pos_df,BLOCK_10 == TRUE & BLOCK_11 ==TRUE & BLOCK_12 ==TRUE))

png(filename="VenSNPs.png", type="cairo",units="px", width=4000, height=4000, res=500, bg="transparent")

venn.plot <- draw.triple.venn(
  area1=loci.b10, area2=loci.b11, area3=loci.b12, 
  n12= loci.b10.b11, n13=loci.b10.b12, 
  n23=loci.b11.b12,scaled=TRUE,
  n123=loci.b10.b11.b12, cex=rep(2, 7),
  category = c("Block 10", "Block 11", "Block 12"),
  fill = cbPaletteSmall3,cat.cex = rep(2, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

## Differentiate Singleton Loci

``` bash
source activate CASE

tail -n +2 singleton.loci | awk '{
  match($1, /_[0-9]+$/)
  print substr($1, 1, RSTART-1) "\t" substr($1, RSTART+1)
}' > sing.pos.txt

mawk 'BEGIN {
  FS = "\t"
  maf_min   = 0.015
  depth_min = 20
  while ((getline line < "sing.pos.txt") > 0) {
    split(line, p, "\t")
    sing[p[1] FS p[2]] = 1
  }
  close("sing.pos.txt")
  print "SNP\tdepth_b10\tmaf_b10\tdepth_b11\tmaf_b11\tdepth_b12\tmaf_b12\tdropout_type"
}
NR == 1 {
  nb10 = nb11 = nb12 = 0
  for (i = 1; i <= NF; i++) {
    if ($i ~ /B10/) c10[++nb10] = i
    if ($i ~ /B11/) c11[++nb11] = i
    if ($i ~ /B12/) c12[++nb12] = i
  }
  next
}
!(($1 FS $2) in sing) { next }
function block_maf(cols, ncols,    i, ct, d, d_total, minor_total, max_ct) {
  d_total = 0; minor_total = 0
  for (i = 1; i <= ncols; i++) {
    split($(cols[i]), ct, ":")
    d = ct[1]+ct[2]+ct[3]+ct[4]
    d_total += d
    if (d > 0) {
      max_ct = ct[1]
      if (ct[2] > max_ct) max_ct = ct[2]
      if (ct[3] > max_ct) max_ct = ct[3]
      if (ct[4] > max_ct) max_ct = ct[4]
      minor_total += d - max_ct
    }
  }
  bd_out   = d_total
  bmaf_out = (d_total > 0) ? minor_total / d_total : -1
}
{
  block_maf(c10, nb10); d10=bd_out; m10=bmaf_out
  block_maf(c11, nb11); d11=bd_out; m11=bmaf_out
  block_maf(c12, nb12); d12=bd_out; m12=bmaf_out
  snp = $1"_"$2

  # Near-fixed in an absent block = adequate depth but MAF below threshold
  lv10 = (d10>=depth_min && m10>=0 && m10<maf_min) ? 1 : 0
  lv11 = (d11>=depth_min && m11>=0 && m11<maf_min) ? 1 : 0
  lv12 = (d12>=depth_min && m12>=0 && m12<maf_min) ? 1 : 0

  type = (lv10||lv11||lv12) ? "low_variation" : "low_coverage"
  printf "%s\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%s\n", \
    snp, d10, m10, d11, m11, d12, m12, type
}' CASE.All.Blocks.sync > singleton.loci.classified.txt

echo "Breakdown by dropout type:"
cut -f8 singleton.loci.classified.txt | sort | uniq -c

awk 'NR==1 || $8=="low_variation" {print $1}' singleton.loci.classified.txt > singleton.loci.lowvar
awk 'NR==1 || $8=="low_coverage"  {print $1}' singleton.loci.classified.txt > singleton.loci.lowcov
```

    ## Breakdown by dropout type:
    ##       1 dropout_type
    ##  185068 low_coverage
    ##  290089 low_variation

# Individual Blocks for outlier testing

#### B10

``` bash
source activate CASE

bash ../scripts/sum.sh <(cut -f1,2,3,13-16 CASE.dp20.Block10.cov.sync) > CASE.B10.IS.sync
paste CASE.B10.IS.sync <(cut -f4 CASE.B10.IS.sync) <(cut -f4 CASE.B10.IS.sync)  > CASE.B10.ISsum3.sync

AV_COV=$(mawk -f ../scripts/add_cov_sync <(cut -f13-16 --complement CASE.dp20.Block10.cov.sync) | mawk '{sum=sum+$21} END {print sum/NR}')

echo $AV_COV

mawk -f ../scripts/add_cov_sync CASE.B10.ISsum3.sync | mawk -v x=$AV_COV '$7 >= x' | mawk '!/CHR/' | cut -f1-6 > CASE.B10.ISsum3nh.sync

mawk -f ../scripts/add_cov_sync_IS CASE.B10.ISsum3.sync | mawk -v x=$AV_COV '$7 >= x'> CASE.B10.ISsum3nh.cov.sync


paste CASE.dp20.Block10.cov.sync <(mawk -f ../scripts/add_cov_sync CASE.B10.ISsum3.sync ) | mawk -v x=$AV_COV '$29 >= x' | cut -f1-22 > temp.sync

mv temp.sync CASE.dp20.Block10.COV.sync

mawk -f ../scripts/add_cov_sync <(cut -f13-16,20- --complement CASE.dp20.Block10.COV.sync) > temp.cov


paste temp.cov <(cut -f7 CASE.B10.ISsum3nh.cov.sync) | mawk -v OFS='\t' '{if ($18 > $19) {$18=$19; print $0} else {print $0}}' | cut -f1-18> CASE.dp20.Block10.COV
```

``` bash
source activate random_draw


python ../scripts/sub_sample.py CASE.B10.ISsum3nh.cov.sync CASE.B10.ISsum3nh.sync CASE.B10.ISs.sync

cat <(echo -e "CHROM\tPOS\tREF\tIS_RS1\tIS_RS2\tIS_RS3") CASE.B10.ISs.sync > CASE.B10.ISsum.sync

paste <(cut -f1,2,3,4 CASE.B10.ISsum.sync) <(cut -f4 CASE.dp20.Block10.COV.sync) <(cut -f5 CASE.B10.ISsum.sync) <(cut -f5 CASE.dp20.Block10.COV.sync) <(cut -f6 CASE.B10.ISsum.sync) <(cut -f6 CASE.dp20.Block10.COV.sync)> CA.B10.sync

head -1 CA.B10.sync

mawk '!/CHR/' CA.B10.sync > CA.B10.input

paste <(cut -f1,2,3,4 CASE.B10.ISsum.sync)  <(cut -f7 CASE.dp20.Block10.COV.sync) <(cut -f5 CASE.B10.ISsum.sync) <(cut -f8 CASE.dp20.Block10.COV.sync) <(cut -f6 CASE.B10.ISsum.sync) <(cut -f9 CASE.dp20.Block10.COV.sync)> CASE.B10.sync

head -1 CASE.B10.sync
tail -n +2 CASE.B10.sync > CASE.B10.input


paste <(cut -f1,2,3,4 CASE.B10.ISsum.sync) <(cut -f17 CASE.dp20.Block10.COV.sync) <(cut -f5 CASE.B10.ISsum.sync) <(cut -f19 CASE.dp20.Block10.COV.sync) <(cut -f6 CASE.B10.ISsum.sync) <(cut -f18 CASE.dp20.Block10.COV.sync) > SE.B10.sync

head -1 SE.B10.sync
tail -n +2 SE.B10.sync > SE.B10.input


paste <(cut -f1,2,3,4 CASE.B10.ISsum.sync) <(cut -f10 CASE.dp20.Block10.COV.sync) <(cut -f5 CASE.B10.ISsum.sync) <(cut -f11 CASE.dp20.Block10.COV.sync) <(cut -f6 CASE.B10.ISsum.sync) <(cut -f12 CASE.dp20.Block10.COV.sync)> CON.B10.sync

head -1 CON.B10.sync
tail -n +2 CON.B10.sync | mawk '!/0:0:0:0:0:0/' > CON.B10.input
```

#### B11

``` bash
AV_COV=$(mawk -f ../scripts/add_cov_sync <(cut -f16-19 --complement CASE.dp20.Block11.cov.sync) | mawk '{sum=sum+$22} END {print sum/NR}')

echo $AV_COV

bash ../scripts/sum.sh <(cut -f1,2,3,16-19 CASE.dp20.Block11.cov.sync) > CASE.B11.IS.sync
paste CASE.B11.IS.sync <(cut -f4 CASE.B11.IS.sync) <(cut -f4 CASE.B11.IS.sync) <(cut -f4 CASE.B11.IS.sync) > CASE.B11.ISsum4.sync

mawk -f ../scripts/add_cov_sync CASE.B11.ISsum4.sync | mawk -v x=$AV_COV '$8 >= x' | mawk '!/CHR/' | cut -f1-7 > CASE.B11.ISsum4nh.sync

paste CASE.dp20.Block11.cov.sync <(mawk -f ../scripts/add_cov_sync CASE.B11.ISsum4.sync ) | mawk -v x=$AV_COV '$34 >= x' | cut -f1-23 > temp.sync

mv temp.sync CASE.dp20.Block11.COV.sync

mawk -f ../scripts/add_cov_sync_IS CASE.B11.ISsum4.sync | mawk -v x=$AV_COV '$8 >= x'  > CASE.B11.ISsum4.cov.sync
mawk -f ../scripts/add_cov_sync <(cut -f16-19 --complement CASE.dp20.Block11.COV.sync) > CASE.dp20.Block11.COV

source activate random_draw


python ../scripts/sub_sample.py CASE.B11.ISsum4.cov.sync CASE.B11.ISsum4nh.sync CASE.B11.ISs.sync

cat <(echo -e "CHROM\tPOS\tREF\tIS_RS1\tIS_RS2\tIS_RS3\tIS_R4") CASE.B11.ISs.sync > CASE.B11.ISsum.sync

paste <(cut -f1,2,3,4 CASE.B11.ISsum.sync) <(cut -f4 CASE.dp20.Block11.COV.sync) <(cut -f5 CASE.B11.ISsum.sync) <(cut -f 5 CASE.dp20.Block11.COV.sync) <(cut -f6 CASE.B11.ISsum.sync) <(cut -f6 CASE.dp20.Block11.COV.sync) <(cut -f7 CASE.B11.ISsum.sync) <(cut -f7 CASE.dp20.Block11.COV.sync) > CA.B11.sync

head -1 CA.B11.sync
tail -n +2 CA.B11.sync > CA.B11.input

paste <(cut -f1,2,3,4 CASE.B11.ISsum.sync) <(cut -f8 CASE.dp20.Block11.COV.sync) <(cut -f5 CASE.B11.ISsum.sync) <(cut -f 9 CASE.dp20.Block11.COV.sync) <(cut -f6 CASE.B11.ISsum.sync) <(cut -f10 CASE.dp20.Block11.COV.sync) <(cut -f7 CASE.B11.ISsum.sync) <(cut -f11 CASE.dp20.Block11.COV.sync) > CASE.B11.sync

head -1 CASE.B11.sync
tail -n +2 CASE.B11.sync > CASE.B11.input

paste <(cut -f1,2,3,4 CASE.B11.ISsum.sync) <(cut -f20 CASE.dp20.Block11.COV.sync) <(cut -f5 CASE.B11.ISsum.sync) <(cut -f 21 CASE.dp20.Block11.COV.sync) <(cut -f6 CASE.B11.ISsum.sync) <(cut -f22 CASE.dp20.Block11.COV.sync) <(cut -f7 CASE.B11.ISsum.sync) <(cut -f23 CASE.dp20.Block11.COV.sync) > SE.B11.sync

head -1 SE.B11.sync
tail -n +2 SE.B11.sync > SE.B11.input

paste <(cut -f1,2,3,4 CASE.B11.ISsum.sync) <(cut -f12 CASE.dp20.Block11.COV.sync) <(cut -f5 CASE.B11.ISsum.sync) <(cut -f 13 CASE.dp20.Block11.COV.sync) <(cut -f6 CASE.B11.ISsum.sync) <(cut -f14 CASE.dp20.Block11.COV.sync) <(cut -f7 CASE.B11.ISsum.sync) <(cut -f15 CASE.dp20.Block11.COV.sync) > CON.B11.sync

head -1 CON.B11.sync
tail -n +2 CON.B11.sync > CON.B11.input
```

#### B12

``` bash
source activate CASE

bash ../scripts/sum.sh <(cut -f1,2,3,16-19 CASE.dp20.Block12.cov.sync) > CASE.B12.IS.sync
paste CASE.B12.IS.sync <(cut -f4 CASE.B12.IS.sync) <(cut -f4 CASE.B12.IS.sync) <(cut -f4 CASE.B12.IS.sync) > CASE.B12.ISsum4.sync

AV_COV=$(mawk -f ../scripts/add_cov_sync <(cut -f16-19 --complement CASE.dp20.Block12.cov.sync) | mawk '{sum=sum+$22} END {print sum/NR}')

echo $AV_COV

mawk -f ../scripts/add_cov_sync CASE.B12.ISsum4.sync | mawk -v x=$AV_COV '$8 >= x' | mawk '!/CHR/' | cut -f1-7 > CASE.B12.ISsum4nh.sync

paste CASE.dp20.Block12.cov.sync <(mawk -f ../scripts/add_cov_sync CASE.B12.ISsum4.sync ) | mawk -v x=$AV_COV '$34 >= x' | cut -f1-26 > temp.sync

mv -f temp.sync CASE.dp20.Block12.COV.sync

mawk -f ../scripts/add_cov_sync_IS CASE.B12.ISsum4.sync | mawk -v x=$AV_COV '$8 >= x'  > CASE.B12.ISsum4nh.cov.sync

mawk -f ../scripts/add_cov_sync <(cut -f16-19,24- --complement CASE.dp20.Block12.COV.sync) > CASE.dp20.Block12.COV


paste CASE.dp20.Block12.COV <(cut -f8 CASE.B12.ISsum4nh.cov.sync) | mawk -v OFS='\t' '{if ($22 > $23) {$22=$23; print $0} else {print $0}}' > temp.cov.sync

cut -f1-22 temp.cov.sync > CASE.dp20.Block12.COV
rm temp.cov.sync
```

``` bash
source activate random_draw

python ../scripts/sub_sample.py CASE.B12.ISsum4nh.cov.sync CASE.B12.ISsum4nh.sync CASE.B12.ISs.sync

cat <(echo -e "CHROM\tPOS\tREF\tIS_RS1\tIS_RS2\tIS_RS3\tIS_R4") CASE.B12.ISs.sync > CASE.B12.ISsum.sync

paste <(cut -f1,2,3,4 CASE.B12.ISsum.sync) <(cut -f4 CASE.dp20.Block12.COV.sync) <(cut -f5 CASE.B12.ISsum.sync) <(cut -f 5 CASE.dp20.Block12.COV.sync) <(cut -f6 CASE.B12.ISsum.sync) <(cut -f6 CASE.dp20.Block12.COV.sync) <(cut -f7 CASE.B12.ISsum.sync) <(cut -f7 CASE.dp20.Block12.COV.sync) > CA.B12.sync

head -1 CA.B12.sync
tail -n +2 CA.B12.sync > CA.B12.input

paste <(cut -f1,2,3,4 CASE.B12.ISsum.sync) <(cut -f8 CASE.dp20.Block12.COV.sync) <(cut -f5 CASE.B12.ISsum.sync) <(cut -f 9 CASE.dp20.Block12.COV.sync) <(cut -f6 CASE.B12.ISsum.sync) <(cut -f10 CASE.dp20.Block12.COV.sync) <(cut -f7 CASE.B12.ISsum.sync) <(cut -f11 CASE.dp20.Block12.COV.sync) > CASE.B12.sync

head -1 CASE.B12.sync
tail -n +2 CASE.B12.sync > CASE.B12.input

paste <(cut -f1,2,3,4 CASE.B12.ISsum.sync) <(cut -f20 CASE.dp20.Block12.COV.sync) <(cut -f5 CASE.B12.ISsum.sync) <(cut -f 21 CASE.dp20.Block12.COV.sync) <(cut -f6 CASE.B12.ISsum.sync) <(cut -f22 CASE.dp20.Block12.COV.sync) <(cut -f7 CASE.B12.ISsum.sync) <(cut -f23 CASE.dp20.Block12.COV.sync) > SE.B12.sync

head -1 SE.B12.sync
tail -n +2 SE.B12.sync > SE.B12.input

paste <(cut -f1,2,3,4 CASE.B12.ISsum.sync) <(cut -f12 CASE.dp20.Block12.COV.sync) <(cut -f5 CASE.B12.ISsum.sync) <(cut -f 13 CASE.dp20.Block12.COV.sync) <(cut -f6 CASE.B12.ISsum.sync) <(cut -f14 CASE.dp20.Block12.COV.sync) <(cut -f7 CASE.B12.ISsum.sync) <(cut -f15 CASE.dp20.Block12.COV.sync) > CON.B12.sync

head -1 CON.B12.sync
tail -n +2 CON.B12.sync > CON.B12.input
```

### Calculate P-values

#### Modified CMH

``` r
ca.b10.sync <- read.sync(file="CA.B10.input", gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
ca.b10.cov <- poolSeq::coverage(ca.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
ca.b10.af <- af(ca.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))

case.b10.sync <- read.sync(file="CASE.B10.input", gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
case.b10.cov <- poolSeq::coverage(case.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
case.b10.af <- af(case.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))

se.b10.sync <- read.sync(file="SE.B10.input", gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
se.b10.cov <- poolSeq::coverage(se.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
se.b10.af <- af(se.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))

con.b10.sync <- read.sync(file="CON.B10.input", gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
con.b10.cov <- poolSeq::coverage(con.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
con.b10.af <- af(con.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

``` r
ca.b10.pvals <- adapted.cmh.test(freq=ca.b10.af, coverage=ca.b10.cov, Ne=rep(250, 3), gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3), poolSize=rep(c(40000,10000), ncol(ca.b10.af)/2), MeanStart = TRUE, mincov =20)

pcadf <- as.data.table(cbind(row.names(ca.b10.af),ca.b10.pvals))
colnames(pcadf) <- c("Locus","PVAL")
pcadf$GROUP <- "PCA"


case.b10.pvals <- adapted.cmh.test(freq=case.b10.af, coverage=case.b10.cov, Ne=rep(250, 3), gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3), poolSize=rep(c(40000,10000), ncol(case.b10.af)/2), MeanStart = TRUE, mincov =20)

pcasedf <- as.data.table(cbind(row.names(case.b10.af),case.b10.pvals))
colnames(pcasedf) <- c("Locus","PVAL")
pcasedf$GROUP <- "PCASE"


se.b10.pvals <- adapted.cmh.test(freq=se.b10.af, coverage=se.b10.cov, Ne=rep(250, 3), gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3), poolSize=rep(c(40000,10000), ncol(se.b10.af)/2), MeanStart = TRUE, mincov =20)

psedf <- as.data.table(cbind(row.names(se.b10.af),se.b10.pvals))
colnames(psedf) <- c("Locus","PVAL")
psedf$GROUP <- "PSE"


con.b10.pvals <- adapted.cmh.test(freq=con.b10.af, coverage=con.b10.cov, Ne=rep(250, 3), gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3), poolSize=rep(c(40000,10000), ncol(con.b10.af)/2), MeanStart = TRUE, mincov =20)

pcondf <- as.data.table(cbind(row.names(con.b10.af),con.b10.pvals))
colnames(pcondf) <- c("Locus","PVAL")
pcondf$GROUP <- "PCON"

B10.pval.table <- bind_rows(pcadf,pcasedf,pcondf,psedf)
B10.pval.table <- B10.pval.table %>%
  separate(Locus, into = c("CHROM", "BP"), sep = "\\.(?=[^.]+$)")
B10.pval.table$PVAL <- as.numeric(B10.pval.table$PVAL)
```

``` r
ca.b11.sync <- read.sync(file="CA.B11.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
ca.b11.cov <- poolSeq::coverage(ca.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
ca.b11.af <- af(ca.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

case.b11.sync <- read.sync(file="CASE.B11.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
case.b11.cov <- poolSeq::coverage(case.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
case.b11.af <- af(case.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

se.b11.sync <- read.sync(file="SE.B11.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
se.b11.cov <- poolSeq::coverage(se.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
se.b11.af <- af(se.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

con.b11.sync <- read.sync(file="CON.B11.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
con.b11.cov <- poolSeq::coverage(con.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
con.b11.af <- af(con.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

``` r
ca.b11.pvals <- adapted.cmh.test(freq=ca.b11.af, coverage=ca.b11.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(ca.b11.af)/2), MeanStart = TRUE, mincov =14)

pcadf <- as.data.table(cbind(row.names(ca.b11.af),ca.b11.pvals))
colnames(pcadf) <- c("Locus","PVAL")
pcadf$GROUP <- "PCA"

case.b11.pvals <- adapted.cmh.test(freq=case.b11.af, coverage=case.b11.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(case.b11.af)/2), MeanStart = TRUE, mincov =14)

pcasedf <- as.data.table(cbind(row.names(case.b11.af),case.b11.pvals))
colnames(pcasedf) <- c("Locus","PVAL")
pcasedf$GROUP <- "PCASE"

se.b11.pvals <- adapted.cmh.test(freq=se.b11.af, coverage=se.b11.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(se.b11.af)/2), MeanStart = TRUE, mincov =14)

psedf <- as.data.table(cbind(row.names(se.b11.af),se.b11.pvals))
colnames(psedf) <- c("Locus","PVAL")
psedf$GROUP <- "PSE"

con.b11.pvals <- adapted.cmh.test(freq=con.b11.af, coverage=con.b11.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(con.b11.af)/2), MeanStart = TRUE, mincov =14)

pcondf <- as.data.table(cbind(row.names(con.b11.af),con.b11.pvals))
colnames(pcondf) <- c("Locus","PVAL")
pcondf$GROUP <- "PCON"


B11.pval.table <- bind_rows(pcadf,pcasedf,pcondf,psedf)
B11.pval.table <- B11.pval.table %>%
  separate(Locus, into = c("CHROM", "BP"), sep = "\\.(?=[^.]+$)")
B11.pval.table$PVAL <- as.numeric(B11.pval.table$PVAL)
```

``` r
ca.b12.sync <- read.sync(file="CA.B12.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
ca.b12.cov <- poolSeq::coverage(ca.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
ca.b12.af <- af(ca.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

case.b12.sync <- read.sync(file="CASE.B12.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
case.b12.cov <- poolSeq::coverage(case.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
case.b12.af <- af(case.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

se.b12.sync <- read.sync(file="SE.B12.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
se.b12.cov <- poolSeq::coverage(se.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
se.b12.af <- af(se.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

con.b12.sync <- read.sync(file="CON.B12.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
con.b12.cov <- poolSeq::coverage(con.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
con.b12.af <- af(con.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

``` r
ca.b12.pvals <- adapted.cmh.test(freq=ca.b12.af, coverage=ca.b12.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(ca.b12.af)/2), MeanStart = TRUE, mincov =20)

pcadf <- as.data.table(cbind(row.names(ca.b12.af),ca.b12.pvals))
colnames(pcadf) <- c("Locus","PVAL")
pcadf$GROUP <- "PCA"

case.b12.pvals <- adapted.cmh.test(freq=case.b12.af, coverage=case.b12.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(case.b12.af)/2), MeanStart = TRUE, mincov =20)

pcasedf <- as.data.table(cbind(row.names(case.b12.af),case.b12.pvals))
colnames(pcasedf) <- c("Locus","PVAL")
pcasedf$GROUP <- "PCASE"

se.b12.pvals <- adapted.cmh.test(freq=se.b12.af, coverage=se.b12.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(se.b12.af)/2), MeanStart = TRUE, mincov =20)

psedf <- as.data.table(cbind(row.names(se.b12.af),se.b12.pvals))
colnames(psedf) <- c("Locus","PVAL")
psedf$GROUP <- "PSE"


con.b12.pvals <- adapted.cmh.test(freq=con.b12.af, coverage=con.b12.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(con.b12.af)/2), MeanStart = TRUE, mincov =20)

pcondf <- as.data.table(cbind(row.names(con.b12.af),con.b12.pvals))
colnames(pcondf) <- c("Locus","PVAL")
pcondf$GROUP <- "PCON"

B12.pval.table <- bind_rows(pcadf,pcasedf,pcondf,psedf)
B12.pval.table <- B12.pval.table %>%
  separate(Locus, into = c("CHROM", "BP"), sep = "\\.(?=[^.]+$)")
B12.pval.table$PVAL <- as.numeric(B12.pval.table$PVAL)
```

### Across all blocks

#### Create sync file

``` bash
source activate CASE

head -1 CASE.All.Blocks.cov.sync > h.temp
tail -n +2 CASE.All.Blocks.cov.sync | sort -k1,1 -k2,2n > b.temp

cat h.temp b.temp > CASE.All.Blocks.cov.s.sync
rm h.temp b.temp
cp CASE.All.Blocks.cov.sync CASE.All.Blocks.cov.s.sync

bash ../scripts/sum.sh <(cut -f1,2,3,37,39,41,43 CASE.All.Blocks.cov.s.sync) > CASE.B10.AB.ISsum.sync
bash ../scripts/sum.sh <(cut -f1,2,3,38,40,42,44 CASE.All.Blocks.cov.s.sync) > CASE.B11.AB.ISsum.sync
bash ../scripts/sum.sh <(cut -f1,2,3,45-48 CASE.All.Blocks.cov.s.sync) > CASE.B12.AB.ISsum.sync


paste CASE.All.Blocks.cov.s.sync <(mawk -f ../scripts/add_cov_sync CASE.B10.AB.ISsum.sync | cut -f5 ) <(mawk -f ../scripts/add_cov_sync CASE.B11.AB.ISsum.sync | cut -f 5) <(mawk -f ../scripts/add_cov_sync CASE.B12.AB.ISsum.sync | cut -f5) | mawk '$63 >72 && $64 > 72 && $65 > 72 && $61 > 1' | cut -f1-62 > CASE.All.Blocks.cov.fil.sync

bash ../scripts/sum.sh <(cut -f1,2,3,37,39,41,43 CASE.All.Blocks.cov.fil.sync) > CASE.B10.AB.ISsum.sync
bash ../scripts/sum.sh <(cut -f1,2,3,38,40,42,44 CASE.All.Blocks.cov.fil.sync) > CASE.B11.AB.ISsum.sync
bash ../scripts/sum.sh <(cut -f1,2,3,45-48 CASE.All.Blocks.cov.fil.sync) > CASE.B12.AB.ISsum.sync

mawk -f ../scripts/add_cov_sync <(cut -f1-3,7,9,12,15,19,22,26,30,32,51,52,55  CASE.All.Blocks.cov.fil.sync ) > B10.AB.cov.sync
mawk -f ../scripts/add_cov_sync <(cut -f1-3,5,10,14,16,23,24,27,29,35,49,54,58 CASE.All.Blocks.cov.fil.sync) > B11.AB.cov.sync
#mawk -f ../scripts/add_cov_sync <(cut -f1-3,5,8,10,14,16,21,23,24,27,29,33,35,49,53,54,58 CASE.All.Blocks.cov.fil.sync) > B11.AB.cov.sync

mawk -f ../scripts/add_cov_sync <(cut -f1-3,4,6,13,17,18,20,28,31,34,56,57,59 CASE.All.Blocks.cov.fil.sync) > B12.AB.cov.sync
#mawk -f ../scripts/add_cov_sync <(cut -f1-3,4,6,11,13,17,18,20,25,28,31,34,36,50,56,57,59 CASE.All.Blocks.cov.fil.sync) > B12.AB.cov.sync

#paste CASE.B10.AB.ISsum.sync <(cut -f4 CASE.B10.AB.ISsum.sync) <(cut -f4 CASE.B10.AB.ISsum.sync)   > CASE.B10.AB.ISsum3.sync

#paste CASE.B11.AB.ISsum.sync <(cut -f4 CASE.B11.AB.ISsum.sync) <(cut -f4 CASE.B11.AB.ISsum.sync)  > CASE.B11.AB.ISsum4.sync

#paste CASE.B12.AB.ISsum.sync <(cut -f4 CASE.B12.AB.ISsum.sync) <(cut -f4 CASE.B12.AB.ISsum.sync)  > CASE.B12.AB.ISsum4.sync

paste CASE.B10.AB.ISsum.sync <(cut -f4 CASE.B10.AB.ISsum.sync) <(cut -f4 CASE.B10.AB.ISsum.sync) <(cut -f4 CASE.B10.AB.ISsum.sync)  > CASE.B10.AB.ISsum3.sync
paste CASE.B11.AB.ISsum.sync <(cut -f4 CASE.B11.AB.ISsum.sync) <(cut -f4 CASE.B11.AB.ISsum.sync) <(cut -f4 CASE.B11.AB.ISsum.sync)  > CASE.B11.AB.ISsum4.sync
paste CASE.B12.AB.ISsum.sync <(cut -f4 CASE.B12.AB.ISsum.sync) <(cut -f4 CASE.B12.AB.ISsum.sync) <(cut -f4 CASE.B12.AB.ISsum.sync)  > CASE.B12.AB.ISsum4.sync


mawk -f ../scripts/add_cov_sync_IS CASE.B10.AB.ISsum3.sync  > CASE.B10.AB.ISsum3.cov.sync
mawk -f ../scripts/add_cov_sync_IS CASE.B11.AB.ISsum4.sync  > CASE.B11.AB.ISsum4.cov.sync
mawk -f ../scripts/add_cov_sync_IS CASE.B12.AB.ISsum4.sync  > CASE.B12.AB.ISsum4.cov.sync

paste B10.AB.cov.sync <(cut -f7 CASE.B10.AB.ISsum3.cov.sync) | mawk -v OFS='\t' '{if (NR > 1 && $19 > $20) {$19=$20; print $0} else {print $0}}' > temp.cov.sync
cut -f1-18 temp.cov.sync > B10.AB.cov.sync

paste B11.AB.cov.sync <(cut -f7 CASE.B11.AB.ISsum4.cov.sync) | mawk -v OFS='\t' '{if (NR > 1 && $19 > $20) {$19=$20; print $0} else {print $0}}' > temp.cov.sync
cut -f1-18 temp.cov.sync > B11.AB.cov.sync

paste B12.AB.cov.sync <(cut -f7 CASE.B12.AB.ISsum4.cov.sync) | mawk -v OFS='\t' '{if (NR > 1 && $19 > $20) {$19=$20; print $0} else {print $0}}' > temp.cov.sync
cut -f1-18 temp.cov.sync > B12.AB.cov.sync
rm temp.cov.sync
```

``` bash
source activate random_draw
python ../scripts/sub_sample.py CASE.B10.AB.ISsum3.cov.sync <( mawk '!/CHR/' CASE.B10.AB.ISsum3.sync) CASE.B10.ISs.AB.sync
python ../scripts/sub_sample.py CASE.B11.AB.ISsum4.cov.sync <( mawk '!/CHR/' CASE.B11.AB.ISsum4.sync) CASE.B11.ISs.AB.sync
python ../scripts/sub_sample.py CASE.B12.AB.ISsum4.cov.sync <(mawk '!/CHR/' CASE.B12.AB.ISsum4.sync) CASE.B12.ISs.AB.sync
```

``` bash
cat <(echo -e "CHROM\tPOS\tREF\tISB11_RS1\tISB11_RS2\tISB11_RS3\tISB11_RS4") CASE.B11.ISs.AB.sync > CASE.B11.ISsum.AB.sync
cat <(echo -e "CHROM\tPOS\tREF\tISB10_RS1\tISB10_RS2\tISB10_RS3\tISB10_RS4") CASE.B10.ISs.AB.sync > CASE.B10.ISsum.AB.sync
cat <(echo -e "CHROM\tPOS\tREF\tISB12_RS1\tISB12_RS2\tISB12_RS3\tISB12_RS4") CASE.B12.ISs.AB.sync > CASE.B12.ISsum.AB.sync

paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f15 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f19 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f22 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f16 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f21 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f23 CASE.All.Blocks.cov.fil.sync) <(cut -f7 CASE.B11.ISsum.AB.sync) <(cut -f24 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f17 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f18 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f20 CASE.All.Blocks.cov.fil.sync) <(cut -f7 CASE.B12.ISsum.AB.sync) <(cut -f25 CASE.All.Blocks.cov.fil.sync) > CASE.AB.sync


tail -n +2 CASE.AB.sync | mawk '$2 != 74609903 && $2 != 93711267 && $2 != 36030885' > CASE.AB.input 

paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f7 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f9 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f12 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f5 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f8 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f10 CASE.All.Blocks.cov.fil.sync) <(cut -f7 CASE.B11.ISsum.AB.sync) <(cut -f14 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f4 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f6 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f11 CASE.All.Blocks.cov.fil.sync) <(cut -f7 CASE.B12.ISsum.AB.sync) <(cut -f13 CASE.All.Blocks.cov.fil.sync) > CA.AB.sync

#paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f7 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f9 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f12 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f5 CASE.All.Blocks.cov.fil.sync)  <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f10 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f14 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f4 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f6 CASE.All.Blocks.cov.fil.sync)  <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f13 CASE.All.Blocks.cov.fil.sync) > CA.AB.sync

tail -n +2 CA.AB.sync | mawk '$2 != 93711267 && $2 != 59024150 && $2 != 2532550'> CA.AB.input 

paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f26 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f30 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f32 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f27 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f29 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f33 CASE.All.Blocks.cov.fil.sync) <(cut -f7 CASE.B11.ISsum.AB.sync) <(cut -f35 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f28 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f31 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f34 CASE.All.Blocks.cov.fil.sync) <(cut -f7 CASE.B12.ISsum.AB.sync) <(cut -f36 CASE.All.Blocks.cov.fil.sync) > CON.AB.sync


tail -n +2 CON.AB.sync | mawk '$2 != 34804390 && $2 != 54895214 && $2 !=  49000140 && $2 != 15140923' > CON.AB.input 

paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f51 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f52 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f55 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f49 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f53 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f54 CASE.All.Blocks.cov.fil.sync) <(cut -f7 CASE.B11.ISsum.AB.sync) <(cut -f58 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f50 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f56 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f57 CASE.All.Blocks.cov.fil.sync) <(cut -f7 CASE.B12.ISsum.AB.sync) <(cut -f59 CASE.All.Blocks.cov.fil.sync) > SE.AB.sync


tail -n +2 SE.AB.sync | mawk '$2 != 67902213 && $2 != 19828063 && $2 != 93711267 ' > SE.AB.input 
```

#### Calculate p-values

``` r
ca.ab.sync <- read.sync(file="CA.AB.input",gen=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),repl=c(1,1,2,2,3,3,4 ,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
ca.ab.cov <- poolSeq::coverage(ca.ab.sync,gen=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),repl=c(1,1,2,2,3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
ca.ab.af <- af(ca.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))

case.ab.sync <- read.sync(file="CASE.AB.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
case.ab.cov <- poolSeq::coverage(case.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
case.ab.af <- af(case.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))

se.ab.sync <- read.sync(file="SE.AB.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
se.ab.cov <- poolSeq::coverage(se.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
se.ab.af <- af(se.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))

con.ab.sync <- read.sync(file="CON.AB.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
con.ab.cov <- poolSeq::coverage(con.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
con.ab.af <- af(con.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ))
```

``` r
ca.ab.pvals <- adapted.cmh.test(freq=ca.ab.af, coverage=ca.ab.cov, Ne=c(250,250,250,1000,1000,1000,1000,1000,1000,1000,1000), gen=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),repl=c(1,1,2,2,3,3,4 ,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ), poolSize=rep(c(40000,10000), ncol(ca.ab.af)/2), MeanStart = TRUE, mincov =14)

pcadf <- as.data.table(cbind(row.names(ca.ab.af),ca.ab.pvals))
colnames(pcadf) <- c("Locus","PVAL")
pcadf$GROUP <- "PCA"

case.ab.pvals <- adapted.cmh.test(freq=case.ab.af, coverage=case.ab.cov, Ne=c(250,250,250,1000,1000,1000,1000,1000,1000,1000,1000), gen=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),repl=c(1,1,2,2,3,3,4 ,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ), poolSize=rep(c(40000,10000), ncol(case.ab.af)/2), MeanStart = TRUE, mincov =14)

pcasedf <- as.data.table(cbind(row.names(case.ab.af),case.ab.pvals))
colnames(pcasedf) <- c("Locus","PVAL")
pcasedf$GROUP <- "PCASE"

se.ab.pvals <- adapted.cmh.test(freq=se.ab.af, coverage=se.ab.cov, Ne=c(250,250,250,1000,1000,1000,1000,1000,1000,1000,1000), gen=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),repl=c(1,1,2,2,3,3,4 ,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ), poolSize=rep(c(40000,10000), ncol(se.ab.af)/2), MeanStart = TRUE, mincov =14)

psedf <- as.data.table(cbind(row.names(se.ab.af),se.ab.pvals))
colnames(psedf) <- c("Locus","PVAL")
psedf$GROUP <- "PSE"

con.ab.pvals <- adapted.cmh.test(freq=con.ab.af, coverage=con.ab.cov, Ne=c(250,250,250,1000,1000,1000,1000,1000,1000,1000,1000), gen=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1),repl=c(1,1,2,2,3,3,4 ,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11 ), poolSize=rep(c(40000,10000), ncol(con.ab.af)/2), MeanStart = TRUE, mincov =14)

pcondf <- as.data.table(cbind(row.names(con.ab.af),con.ab.pvals)) 
colnames(pcondf) <- c("Locus","PVAL")
pcondf$GROUP <- "PCON"

AB.pval.table <- bind_rows(pcadf,pcasedf,pcondf,psedf)
AB.pval.table <- AB.pval.table %>%
  separate(Locus, into = c("CHROM", "BP"), sep = "\\.(?=[^.]+$)")
AB.pval.table$PVAL <- as.numeric(AB.pval.table$PVAL)
```

## Determine Significant Loci

``` r
# ── Step 1: Convert raw p-values to q-values for each block ──────────────────
# Qvalue_convert applies the qvalue package (pi0 = 1) to CMH p-value tables,
# returning a data frame with q-values for each treatment (QCA, QCASE, QSE, QCON)

Qvalue_convert(B10.pval.table) -> B10.pv
Qvalue_convert(B11.pval.table) -> B11.pv
Qvalue_convert(B12.pval.table) -> B12.pv
Qvalue_convert(AB.pval.table) -> AB.pv
```

``` r
# ── Step 2: Define FDR thresholds ────────────────────────────────────────────
# alpha  = primary FDR threshold applied to treatment q-values (QCA, QCASE, QSE)
# alpha2 = CON q-value floor: SNPs must have QCON > alpha2 to exclude loci
#          that are also differentiating in the control (likely artifacts or
#          family effects rather than treatment-driven selection)
#
# Per-block alpha corrections (applied to Convergent and Private tiers):
#   B10.alpha = 0.001 — 10x stricter than base; B10 p-values are inflated by
#               genetic drift from fewer spawning broodstock, so a tighter
#               threshold is needed to prevent B10 from dominating non-core tiers.
#   B11.alpha = 0.1   — 10x looser than base; B11 has lower coverage and fewer
#               loci, so statistical power is reduced and the threshold is relaxed.
#   Core tiers (Tier A/B) do not apply these corrections because requiring signal
#   across multiple independent blocks already controls for private drift.

# ── Helper: apply_treatment_filter ───────────────────────────────────────────
# After subsetting to SNPs significant for one treatment, zero out the Sig.*
# flags for the other two treatments. This ensures that when per-treatment
# subsets are combined with bind_rows(), each row is attributed to a single
# stressor and doesn't carry spurious significance flags from the subsetting
# source data frame.
apply_treatment_filter <- function(df, keep_treatment) {
  treatments <- c("Sig.CA", "Sig.CASE", "Sig.SE")
  other <- setdiff(treatments, keep_treatment)
  df[other] <- FALSE
  df
}

# top1pct_filter: within the already-significant loci (df), retain only those
# whose treatment q-value falls in the top 1% of that significant set.
# Quantiles are computed only on loci with QCON > 0.1 so that CON-significant
# loci with inflated treatment signal do not pull the threshold up against
# genuinely treatment-specific outliers.
top1pct_filter <- function(df) {
  df_con <- subset(df, QCON > 0.1)
  subset(df,
    (Sig.CASE == TRUE & QCASE < quantile(df_con$QCASE, na.rm = TRUE, probs = 0.01)) |
    (Sig.CA   == TRUE & QCA   < quantile(df_con$QCA,   na.rm = TRUE, probs = 0.01)) |
    (Sig.SE   == TRUE & QSE   < quantile(df_con$QSE,   na.rm = TRUE, probs = 0.01))
  )
}

# ── Step 3: Pool all three blocks (unfiltered) ───────────────────────────────
# pp1 is used later to identify CON-significant SNPs across the per-block data
pp1 <- rbind(B10.pv, B11.pv, B12.pv)

# ══════════════════════════════════════════════════════════════════════════════
# TIER A: q < 0.10, significant in ALL THREE BLOCKS (Sig.Loci.3)
# ══════════════════════════════════════════════════════════════════════════════
alpha  <- 0.10
alpha2 <- 0.10


ppsig.B10.1 <- Significat_subset(B10.pv, alpha,     alpha2)
```

    ## [1] 495239
    ## [1] 85765
    ## [1] 0.173179

``` r
ppsig.B11.1 <- Significat_subset(B11.pv, alpha, alpha2)
```

    ## [1] 175041
    ## [1] 14925
    ## [1] 0.08526574

``` r
ppsig.B12.1 <- Significat_subset(B12.pv, alpha,     alpha2)
```

    ## [1] 580904
    ## [1] 77241
    ## [1] 0.1329669

``` r
ppsig.B10.1$BLOCK <- 10
ppsig.B11.1$BLOCK <- 11
ppsig.B12.1$BLOCK <- 12

ppsig3.FDR10 <- rbind(ppsig.B10.1, ppsig.B11.1, ppsig.B12.1)

# Retain only SNPs significant for a given treatment in all 3 blocks (n() > 2),
# then zero out the other treatment flags before combining
all_sig.CA   <- subset(ppsig3.FDR10, Sig.CA   == TRUE) %>%
                  group_by(SNP) %>% filter(n() > 2) %>%
                  apply_treatment_filter("Sig.CA")

all_sig.CASE <- subset(ppsig3.FDR10, Sig.CASE == TRUE) %>%
                  group_by(SNP) %>% filter(n() > 2) %>%
                  apply_treatment_filter("Sig.CASE")

all_sig.SE   <- subset(ppsig3.FDR10, Sig.SE   == TRUE) %>%
                  group_by(SNP) %>% filter(n() > 2) %>%
                  apply_treatment_filter("Sig.SE")

all_sig <- bind_rows(all_sig.CA, all_sig.CASE, all_sig.SE)
write.table(all_sig, "Sig.Loci.3", sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TIER B: q < 0.05, significant in AT LEAST TWO BLOCKS (Sig.Loci.2)
# ══════════════════════════════════════════════════════════════════════════════
alpha  <- 0.05
alpha2 <- 0.10
B11.alpha <- 0.1 # Block 11 correction: 0.10

ppsig.B10.05 <- Significat_subset(B10.pv, alpha,     alpha2)
```

    ## [1] 495239
    ## [1] 59545
    ## [1] 0.1202349

``` r
ppsig.B11.05 <- Significat_subset(B11.pv, B11.alpha, alpha2)
```

    ## [1] 175041
    ## [1] 14925
    ## [1] 0.08526574

``` r
ppsig.B12.05 <- Significat_subset(B12.pv, alpha,     alpha2)
```

    ## [1] 580904
    ## [1] 42661
    ## [1] 0.07343898

``` r
ppsig.B10.05$BLOCK <- 10
ppsig.B11.05$BLOCK <- 11
ppsig.B12.05$BLOCK <- 12

ppsig3.FDR05 <- rbind(ppsig.B10.05, ppsig.B11.05, ppsig.B12.05)

# Retain SNPs significant for a given treatment in at least 2 blocks (n() > 1)
multi_sig.CA   <- subset(ppsig3.FDR05, Sig.CA   == TRUE) %>%
                    group_by(SNP) %>% filter(n() > 1) %>%
                    apply_treatment_filter("Sig.CA")

multi_sig.CASE <- subset(ppsig3.FDR05, Sig.CASE == TRUE) %>%
                    group_by(SNP) %>% filter(n() > 1) %>%
                    apply_treatment_filter("Sig.CASE")

multi_sig.SE   <- subset(ppsig3.FDR05, Sig.SE   == TRUE) %>%
                    group_by(SNP) %>% filter(n() > 1) %>%
                    apply_treatment_filter("Sig.SE")

multi_sig <- rbind(multi_sig.CA, multi_sig.CASE, multi_sig.SE)
write.table(multi_sig, "Sig.Loci.2", sep = "\t", row.names = FALSE, quote = FALSE)

# ── Per-block q < 0.01 intermediates ─────────────────────────────────────────
# ppsig3.FDR01 feeds PRIVATE TIER and CONVERGENT TIER:
alpha  <- 0.01
alpha2 <- 0.10
B11.alpha <- 0.1  # Block 11 correction: 0.1
B10.alpha <- 0.001 # Block 10 correction: 0.001

ppsig.B10.01 <- Significat_subset(B10.pv, B10.alpha,     alpha2)
```

    ## [1] 495239
    ## [1] 13817
    ## [1] 0.02789966

``` r
ppsig.B11.01 <- Significat_subset(B11.pv, B11.alpha, alpha2)
```

    ## [1] 175041
    ## [1] 14925
    ## [1] 0.08526574

``` r
ppsig.B12.01 <- Significat_subset(B12.pv, alpha,     alpha2)
```

    ## [1] 580904
    ## [1] 12544
    ## [1] 0.02159393

``` r
ppsig.B10.01$BLOCK <- 10
ppsig.B11.01$BLOCK <- 11
ppsig.B12.01$BLOCK <- 12

ppsig3.FDR01 <- rbind(ppsig.B10.01, ppsig.B11.01, ppsig.B12.01)
write.table(ppsig3.FDR01,    "Sig.Loci.FDR01.1",       sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ppsig.B10.01,    "Sig.Loci.FDR01.Block10",  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ppsig.B11.01,    "Sig.Loci.FDR01.Block11",  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ppsig.B12.01,    "Sig.Loci.FDR01.Block12",  sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# PRIVATE TIER: private alpha thresholds AND top 1% of q-values
# ══════════════════════════════════════════════════════════════════════════════
# top1pct_filter selects the top 1% strongest signals among already-significant
# loci (quantile computed on CON-clean subset; QCON > 0.1).
# Block alphas: B10 = 0.001 (drift correction), B11 = 0.1 (power correction),
# B12 = 0.01 (base). Restricted to singletons via semi_join below.

ppsig.B10 <- Significat_subset(B10.pv, B10.alpha,     alpha2) %>% top1pct_filter()
```

    ## [1] 495239
    ## [1] 13817
    ## [1] 0.02789966

``` r
ppsig.B11 <- Significat_subset(B11.pv, B11.alpha, alpha2) %>% top1pct_filter()
```

    ## [1] 175041
    ## [1] 14925
    ## [1] 0.08526574

``` r
ppsig.B12 <- Significat_subset(B12.pv, alpha,     alpha2) %>% top1pct_filter()
```

    ## [1] 580904
    ## [1] 12544
    ## [1] 0.02159393

``` r
ppsig.B10$BLOCK <- 10
ppsig.B11$BLOCK <- 11
ppsig.B12$BLOCK <- 12

singleton   <- read.table("singleton.loci",         header = TRUE)
sing.lowvar <- read.table("singleton.loci.lowvar",  header = TRUE)
sing.lowcov <- read.table("singleton.loci.lowcov",  header = TRUE)

Priv_all       <- rbind(ppsig.B10, ppsig.B11, ppsig.B12)
ppsig1.s        <- semi_join(Priv_all, singleton,   by = "SNP")
ppsig1.s.lowvar <- semi_join(Priv_all, sing.lowvar, by = "SNP")
ppsig1.s.lowcov <- semi_join(Priv_all, sing.lowcov, by = "SNP")

# ══════════════════════════════════════════════════════════════════════════════
# CONVERGENT TIER: per-block q < 0.01 in ≥1 block AND q < 0.0001 in all-blocks (AB)
# ══════════════════════════════════════════════════════════════════════════════
# These loci have modest per-block signal that cross-validates against the
# pooled all-blocks analysis. Detected only in the pooled across-spawn analysis rather than any single
# spawning block — indicating consistent but weak selection below individual-block
# detection power. BLOCK 33 is a dummy identifier for the AB dataset.
ppsig1.ab       <- Significat_subset(AB.pv, 0.0001, alpha2)
```

    ## [1] 396601
    ## [1] 1917
    ## [1] 0.004833573

``` r
ppsig1.ab$BLOCK <- 33

conv.SE <- bind_rows(subset(ppsig3.FDR01, Sig.SE   == TRUE),
                    subset(ppsig1.ab,    Sig.SE   == TRUE)) %>%
          group_by(SNP) %>% filter(n() > 1) %>%
          apply_treatment_filter("Sig.SE")

conv.CASE <- bind_rows(subset(ppsig3.FDR01, Sig.CASE == TRUE),
                      subset(ppsig1.ab,    Sig.CASE == TRUE)) %>%
            group_by(SNP) %>% filter(n() > 1) %>%
            apply_treatment_filter("Sig.CASE")

conv.CA <- bind_rows(subset(ppsig3.FDR01, Sig.CA   == TRUE),
                    subset(ppsig1.ab,    Sig.CA   == TRUE)) %>%
          group_by(SNP) %>% filter(n() > 1) %>%
          apply_treatment_filter("Sig.CA")

convergent.sig <- bind_rows(conv.SE, conv.CASE, conv.CA)

# ══════════════════════════════════════════════════════════════════════════════
# CON FILTER
# ══════════════════════════════════════════════════════════════════════════════
conpp1   <- subset(pp1,   QCON < 0.01)
conpp.AB <- subset(AB.pv, QCON < 0.001)
all_con  <- bind_rows(conpp1, conpp.AB) %>% group_by(SNP) %>% filter(n() > 1)
all_con  <- subset(all_con, QCON / QCASE < 10 & QCON / QCA < 10 & QCON / QSE < 10)
write.table(all_con, "ALL.CON.SNPS", sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# FINAL UNION: all three tiers combined, CON-filtered
# ══════════════════════════════════════════════════════════════════════════════
total.sig <- anti_join(
  bind_rows(all_sig, multi_sig, convergent.sig, ppsig1.s),
  all_con,
  by = "SNP"
)
write.table(total.sig, "Total.Significant.uf.Loci", sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# THREE-WAY SPLIT of significant loci
# ══════════════════════════════════════════════════════════════════════════════
# Tiers are assigned hierarchically (Core > Convergent > Private) so that
# each locus appears in exactly one tier.
#
# Exception: low-variation singletons that appear in the Convergent candidates
# are reclassified as Private. For these loci the Convergent all-blocks
# signal is driven by a single polymorphic spawn — the other two blocks are
# near-fixed (MAF < 0.015), so the AB result is not cross-block validation.
# priv_source records whether a Private locus originated from the
# Convergent candidates ("Convergent") or directly from per-block outlier filtering ("Single").

# Core (Tiers A and B): strong replicated signal across multiple independent
# spawning blocks. Consistent response regardless of private variation
# in effluent composition or larval pool standing genetic variation.
core.sig <- anti_join(
  bind_rows(all_sig, multi_sig),
  all_con,
  by = "SNP"
)
write.table(core.sig, "Sig.Loci.Core", sep = "\t", row.names = FALSE, quote = FALSE)

# Convergent: Core excluded first; then lowvar singletons split off to Private.
# Remaining loci cross-validate a per-block signal against the all-blocks analysis
# across spawns with adequate MAF — genuine weak consistent selection.
conv_candidates <- anti_join(convergent.sig, core.sig, by = "SNP")
conv_to_priv      <- semi_join(conv_candidates, sing.lowvar, by = "SNP") %>%
                  anti_join(all_con, by = "SNP") %>%
                  mutate(priv_source = "Convergent")
convergent.sig  <- anti_join(conv_candidates, sing.lowvar, by = "SNP") %>%
                  anti_join(all_con, by = "SNP")
write.table(convergent.sig, "Sig.Loci.Convergent", sep = "\t", row.names = FALSE, quote = FALSE)

# Private: two sources combined, annotated by priv_source.
# "Single"    — top-5% per-block outliers among singleton loci (testable in 1–2 blocks
#               only due to low coverage or near-fixation in absent spawns).
# "Convergent" — lowvar singletons reclassified from Convergent: near-fixed in absent
#               blocks, so AB signal was single-spawn-driven rather than cross-validated.
priv_single <- anti_join(
  anti_join(
    anti_join(ppsig1.s, core.sig, by = "SNP"),
    convergent.sig, by = "SNP"
  ),
  all_con,
  by = "SNP"
) %>% mutate(priv_source = "Single")

private.sig <- bind_rows(priv_single, conv_to_priv)
write.table(private.sig, "Sig.Loci.Private", sep = "\t", row.names = FALSE, quote = FALSE)
```

### Proximity Control Loci filter

``` bash
source activate CASE

cut -f2,3 Total.Significant.uf.Loci | tail -n +2 | sort -k1,2 | uniq | mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > CASE.Total.Significant.uf.loci.bed

bedtools intersect -wb -a  <(mawk '{print $2"\t"$3-1"\t"$3}' ALL.CON.SNPS ) -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > CON.LOCI

bedtools intersect -wb -a  <(mawk '{print $2"\t"$3-1"\t"$3}' ALL.CON.SNPS ) -b ~/CASE/analysis/sorted.ref3.0.gene.bed | cut -f4- | bedtools sort | bedtools merge > con.loci.filter.bed

cat <(echo -e "SNP\tCHROM\tBP") <(bedtools intersect -a CASE.Total.Significant.uf.loci.bed -b con.loci.filter.bed | mawk '{print $1"_"$3"\t"$1"\t"$3}') > ToFilter.Con.Loci.txt
```

``` r
con_filter <- read.table("ToFilter.Con.Loci.txt",header=TRUE)

total.sig <- unique(anti_join(total.sig,con_filter, by = "SNP"))

write.table(total.sig, "Total.Significant.Loci", sep="\t", row.names = FALSE, quote = FALSE)


grouped_total.sig <- group_and_average(total.sig)

write.table(grouped_total.sig, "Total.Significant.Grouped.Loci", sep="\t", row.names = FALSE, quote = FALSE)

totalsig.CA <- subset(total.sig, Sig.CA == TRUE & Sig.SE == FALSE )
totalsig.SE <- subset(total.sig, Sig.CA == FALSE & Sig.SE == TRUE )
totalsig.both <- subset(total.sig, Sig.CA == TRUE & Sig.SE == TRUE )
totalsig.CASE <- subset(total.sig, Sig.CASE == TRUE )

write.table(totalsig.CA, "sig.ca", sep="\t", row.names = FALSE, quote = FALSE)
write.table(totalsig.SE, "sig.se", sep="\t", row.names = FALSE, quote = FALSE)
write.table(totalsig.both, "sig.both", sep="\t", row.names = FALSE, quote = FALSE)
write.table(totalsig.CASE, "sig.case", sep="\t", row.names = FALSE, quote = FALSE)
```

``` bash
cat <(echo -e "SNP\tGroup") <(mawk '{if (NR > 1) print $1"\tSE"}' sig.se ) <(mawk '{if (NR > 1) print $1"\tCA"}' sig.ca ) <(mawk '{if (NR > 1) print $1"\tBoth"}' sig.both ) > sig.table
```

### Tier-specific CON filtering and treatment splits

``` bash
source activate CASE

# ─── bedtools CON-filtering for each locus tier ──────────────────────────────
# con.loci.filter.bed (built above from ALL.CON.SNPS x gene model) is
# reference-derived and applies equally to all tiers. core.sig, convergent.sig,
# and private.sig were already R-level filtered with anti_join(... , all_con);
# this step removes any remaining loci overlapping CON-associated gene regions.

# Core (Tiers A and B: strong replicated signal across multiple blocks)
cut -f2,3 Sig.Loci.Core | tail -n +2 | sort -k1,2 | uniq | \
  mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > Core.loci.bed
cat <(echo -e "SNP\tCHROM\tBP") \
  <(bedtools intersect -a Core.loci.bed -b con.loci.filter.bed | \
    mawk '{print $1"_"$3"\t"$1"\t"$3}') > ToFilter.Con.Core.Loci.txt

# Convergent (q<0.01 per-block AND q<0.0001 all-blocks; consistent weak signal)
cut -f2,3 Sig.Loci.Convergent | tail -n +2 | sort -k1,2 | uniq | \
  mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > Conv.loci.bed
cat <(echo -e "SNP\tCHROM\tBP") \
  <(bedtools intersect -a Conv.loci.bed -b con.loci.filter.bed | \
    mawk '{print $1"_"$3"\t"$1"\t"$3}') > ToFilter.Con.Conv.Loci.txt

# Private (top 5% outliers among singletons; filename quoted due to hyphen)
cut -f2,3 "Sig.Loci.Private" | tail -n +2 | sort -k1,2 | uniq | \
  mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > Priv.loci.bed
cat <(echo -e "SNP\tCHROM\tBP") \
  <(bedtools intersect -a Priv.loci.bed -b con.loci.filter.bed | \
    mawk '{print $1"_"$3"\t"$1"\t"$3}') > ToFilter.Con.Priv.Loci.txt
```

``` r
# ─── Secondary CON filter (bedtools-derived) for each tier ───────────────────
# Applies the gene-region-based CON filter (bedtools intersect result) on top
# of the R-level anti_join filter already applied during locus definition.

# Core
con_filter.core <- read.table("ToFilter.Con.Core.Loci.txt", header = TRUE)
core.sig <- unique(anti_join(core.sig, con_filter.core, by = "SNP"))
write.table(core.sig, "Core.Significant.Loci",
            sep = "\t", row.names = FALSE, quote = FALSE)
grouped_core.sig <- group_and_average(core.sig)
write.table(grouped_core.sig, "Core.Significant.Grouped.Loci",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Convergent
con_filter.conv <- read.table("ToFilter.Con.Conv.Loci.txt", header = TRUE)
convergent.sig <- unique(anti_join(convergent.sig, con_filter.conv, by = "SNP"))
write.table(convergent.sig, "Conv.Significant.Loci",
            sep = "\t", row.names = FALSE, quote = FALSE)
grouped_convergent.sig <- group_and_average(convergent.sig)
write.table(grouped_convergent.sig, "Conv.Significant.Grouped.Loci",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Private
con_filter.priv <- read.table("ToFilter.Con.Priv.Loci.txt", header = TRUE)
private.sig <- unique(anti_join(private.sig, con_filter.priv, by = "SNP"))
write.table(private.sig, "Priv.Significant.Loci",
            sep = "\t", row.names = FALSE, quote = FALSE)
grouped_private.sig <- group_and_average(private.sig)
write.table(grouped_private.sig, "Priv.Significant.Grouped.Loci",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ─── Per-treatment splits for each tier ──────────────────────────────────────
# Mirrors the total.sig treatment splits above.
# Sig.CA, Sig.SE, Sig.CASE are logical columns; "both" = loci sig in both CA and SE.
# Empty subsets are possible for sparse tiers in some treatment categories.

# Core
coresig.CA   <- subset(core.sig, Sig.CA == TRUE  & Sig.SE == FALSE)
coresig.SE   <- subset(core.sig, Sig.CA == FALSE & Sig.SE == TRUE)
coresig.both <- subset(core.sig, Sig.CA == TRUE  & Sig.SE == TRUE)
coresig.CASE <- subset(core.sig, Sig.CASE == TRUE)
write.table(coresig.CA,   "sig.core.ca",   sep = "\t", row.names = FALSE, quote = FALSE)
write.table(coresig.SE,   "sig.core.se",   sep = "\t", row.names = FALSE, quote = FALSE)
write.table(coresig.both, "sig.core.both", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(coresig.CASE, "sig.core.case", sep = "\t", row.names = FALSE, quote = FALSE)

# Convergent
convsig.CA   <- subset(convergent.sig, Sig.CA == TRUE  & Sig.SE == FALSE)
convsig.SE   <- subset(convergent.sig, Sig.CA == FALSE & Sig.SE == TRUE)
convsig.both <- subset(convergent.sig, Sig.CA == TRUE  & Sig.SE == TRUE)
convsig.CASE <- subset(convergent.sig, Sig.CASE == TRUE)
write.table(convsig.CA,   "sig.conv.ca",   sep = "\t", row.names = FALSE, quote = FALSE)
write.table(convsig.SE,   "sig.conv.se",   sep = "\t", row.names = FALSE, quote = FALSE)
write.table(convsig.both, "sig.conv.both", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(convsig.CASE, "sig.conv.case", sep = "\t", row.names = FALSE, quote = FALSE)

# Private
bssig.CA   <- subset(private.sig, Sig.CA == TRUE  & Sig.SE == FALSE)
bssig.SE   <- subset(private.sig, Sig.CA == FALSE & Sig.SE == TRUE)
bssig.both <- subset(private.sig, Sig.CA == TRUE  & Sig.SE == TRUE)
bssig.CASE <- subset(private.sig, Sig.CASE == TRUE)
write.table(bssig.CA,   "sig.priv.ca",   sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bssig.SE,   "sig.priv.se",   sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bssig.both, "sig.priv.both", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bssig.CASE, "sig.priv.case", sep = "\t", row.names = FALSE, quote = FALSE)
```

``` bash
# ─── Per-tier sig.table files (mirrors sig.table for total) ──────────────────
for TIER in Core Conv Priv; do
  TIER_LC="${TIER,,}"
  cat <(echo -e "SNP\tGroup") \
    <(mawk '{if (NR > 1) print $1"\tSE"}'   sig.${TIER_LC}.se  ) \
    <(mawk '{if (NR > 1) print $1"\tCA"}'   sig.${TIER_LC}.ca  ) \
    <(mawk '{if (NR > 1) print $1"\tBoth"}' sig.${TIER_LC}.both) > sig.${TIER_LC}.table
done
```

``` bash
source activate CASE

cut -f2,3 Total.Significant.Loci | tail -n +2 | sort -k1,2 | uniq | mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > CASE.Total.Significant.loci.bed
bedtools slop -b 100 -i CASE.Total.Significant.loci.bed  -g ../reference.fasta.fai  | bedtools merge -i - >  CASE.Total.Significant.loci.intervals.bed


bedtools intersect -wb -a CASE.Total.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.TOTAL.CASE.LOC
```

``` bash
source activate CASE

mawk '$14 == "TRUE" ' Total.Significant.Loci | sort -k1,2 | mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CA.Significant.loci.bed
bedtools intersect -wb -a CA.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.1.CA.LOC

mawk '$15 =="TRUE"' Total.Significant.Loci | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}'| uniq  > SE.Significant.loci.bed
bedtools intersect -wb -a SE.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.1.SE.LOC

mawk '$13 == "TRUE"' Total.Significant.Loci | sort -k1,2 |   mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CASE.Significant.loci.bed
bedtools intersect -wb -a CASE.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.1.CASE.LOC

mawk '$11 == "TRUE" && $13 == "TRUE"  && $12 == "FALSE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > SE.CASE.Significant.loci.bed

bedtools intersect -wb -a SE.CASE.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.SE.CASE.LOC

mawk '$11 == "TRUE" && $13 == "TRUE"  && $12 == "TRUE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > SE.CASE.CA.Significant.loci.bed

bedtools intersect -wb -a SE.CASE.CA.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.SE.CASE.CA.LOC

mawk '$11 == "FALSE" && $13 == "TRUE"  && $12 == "TRUE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > SE.CA.Significant.loci.bed

mawk '$11 == "TRUE" && $13 == "FALSE"  && $12 == "TRUE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > CA.CASE.Significant.loci.bed

mawk '$11 == "FALSE" && $13 == "TRUE"  && $12 == "FALSE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > SE.ONLY.Significant.loci.bed

mawk '$11 == "TRUE" && $13 == "FALSE"  && $12 == "FALSE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > CASE.ONLY.Significant.loci.bed

mawk '$11 == "FALSE" && $13 == "FALSE"  && $12 == "TRUE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > CA.ONLY.Significant.loci.bed

cat CA.ONLY.Significant.loci.bed CA.CASE.Significant.loci.bed| sort|  uniq  > CA.factor.Significant.loci.bed

cat SE.ONLY.Significant.loci.bed SE.CASE.Significant.loci.bed| sort | uniq  > SE.factor.Significant.loci.bed
```

``` bash
source activate CASE


#B10
mawk '$14 == "TRUE"' Sig.Loci.FDR01.Block10 | sort -k1,2  | mawk '{print $2 "\t" $3-1 "\t" $3}'| uniq  > CA.Significant.B10.loci.bed
bedtools intersect -wb -a CA.Significant.B10.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B10.CA.LOC

mawk '$15 == "TRUE" ' Sig.Loci.FDR01.Block10 | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > SE.Significant.B10.loci.bed
bedtools intersect -wb -a SE.Significant.B10.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B10.SE.LOC

mawk '$13 == "TRUE" ' Sig.Loci.FDR01.Block10 | sort -k1,2 |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CASE.Significant.B10.loci.bed
bedtools intersect -wb -a CASE.Significant.B10.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B10.CASE.LOC

#B11
mawk '$14 == "TRUE"' Sig.Loci.FDR01.Block11 | sort -k1,2  | mawk '{print $2 "\t" $3-1 "\t" $3}'  | uniq > CA.Significant.B11.loci.bed
bedtools intersect -wb -a CA.Significant.B11.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B11.CA.LOC

mawk '$15 == "TRUE" ' Sig.Loci.FDR01.Block11 | sort -k1,2 |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq > SE.Significant.B11.loci.bed
bedtools intersect -wb -a SE.Significant.B11.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B11.SE.LOC

mawk '$13 == "TRUE" ' Sig.Loci.FDR01.Block11 | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CASE.Significant.B11.loci.bed
bedtools intersect -wb -a CASE.Significant.B11.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B11.CASE.LOC

#B12
mawk '$14 == "TRUE" ' Sig.Loci.FDR01.Block12 | sort -k1,2 | mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CA.Significant.B12.loci.bed
bedtools intersect -wb -a CA.Significant.B12.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B12.CA.LOC

mawk '$15 == "TRUE" ' Sig.Loci.FDR01.Block12 | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > SE.Significant.B12.loci.bed
bedtools intersect -wb -a SE.Significant.B12.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B12.SE.LOC

mawk '$13 == "TRUE" ' Sig.Loci.FDR01.Block12 | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CASE.Significant.B12.loci.bed
bedtools intersect -wb -a CASE.Significant.B12.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B12.CASE.LOC
```

### Tier-specific gene annotation

``` r
# ─── Column index verification ────────────────────────────────────────────────
# The bash chunks below use fixed column indices. This block verifies those
# indices against the actual file headers before any gene annotation runs.
# stop() halts the pipeline with a diagnostic message on mismatch (e.g., if
# group_and_average column order changes upstream).
#
# Expected schema for *.Significant.Loci (same as Total.Significant.Loci):
#   col 13 = Sig.CASE, col 14 = Sig.CA, col 15 = Sig.SE
# Expected schema for *.Significant.Grouped.Loci (group_and_average output):
#   col 11 = Sig.CASE, col 12 = Sig.CA, col 13 = Sig.SE, col 14 = CHROM, col 15 = BP

verify_col_indices <- function(filepath, expected) {
  hdr <- colnames(read.table(filepath, header = TRUE, nrows = 1, sep = "\t"))
  for (nm in names(expected)) {
    idx <- expected[[nm]]
    actual <- if (length(hdr) >= idx) hdr[idx] else "<missing>"
    if (actual != nm) {
      stop(sprintf(
        "Column index mismatch in '%s': expected col %d = '%s', found '%s'.\nFull header: %s",
        filepath, idx, nm, actual,
        paste(sprintf("[%d]%s", seq_along(hdr), hdr), collapse = " ")
      ))
    }
  }
  invisible(TRUE)
}

loci_expected    <- list(Sig.CASE = 13L, Sig.CA = 14L, Sig.SE = 15L)
grouped_expected <- list(Sig.CASE = 11L, Sig.CA = 12L, Sig.SE = 13L, CHROM = 14L, BP = 15L)

# Baseline check on total files to confirm column schema before checking tiers
verify_col_indices("Total.Significant.Loci",         loci_expected)
verify_col_indices("Total.Significant.Grouped.Loci", grouped_expected)
message("Total: column indices verified")
```

    ## Total: column indices verified

``` r
for (tier in c("Core", "Conv", "Priv")) {
  verify_col_indices(paste0(tier, ".Significant.Loci"),         loci_expected)
  verify_col_indices(paste0(tier, ".Significant.Grouped.Loci"), grouped_expected)
  message(tier, ": column indices verified")
}
```

    ## Core: column indices verified

    ## Conv: column indices verified

    ## Priv: column indices verified

``` r
message("All column index checks passed. Proceeding to gene annotation.")
```

    ## All column index checks passed. Proceeding to gene annotation.

``` bash
source activate CASE

GENE_BED=~/CASE/analysis/sorted.ref3.0.gene.bed

# ─── Total LOC gene lists per tier ───────────────────────────────────────────

cut -f2,3 Core.Significant.Loci | tail -n +2 | sort -k1,2 | uniq | \
  mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > Core.Significant.loci.bed
bedtools intersect -wb -a Core.Significant.loci.bed -b $GENE_BED | \
  grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.Core.TOTAL.LOC

cut -f2,3 Priv.Significant.Loci | tail -n +2 | sort -k1,2 | uniq | \
  mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > Priv.Significant.loci.bed
bedtools intersect -wb -a Priv.Significant.loci.bed -b $GENE_BED | \
  grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.Priv.TOTAL.LOC

cut -f2,3 Conv.Significant.Loci | tail -n +2 | sort -k1,2 | uniq | \
  mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > Conv.Significant.loci.bed
bedtools intersect -wb -a Conv.Significant.loci.bed -b $GENE_BED | \
  grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.Conv.TOTAL.LOC
```

``` bash
source activate CASE

GENE_BED=~/CASE/analysis/sorted.ref3.0.gene.bed

# ─── Per-treatment LOC lists per tier ────────────────────────────────────────
# Column indices (verified by the R chunk above):
#   *.Significant.Loci:         col 13 = Sig.CASE, col 14 = Sig.CA, col 15 = Sig.SE
#   *.Significant.Grouped.Loci: col 11 = Sig.CASE, col 12 = Sig.CA, col 13 = Sig.SE
#                               col 14 = CHROM, col 15 = BP
# Empty output files are expected for sparse tiers in some treatment categories.

for TIER in Core Conv Priv; do
  LOCI="${TIER}.Significant.Loci"

  # CA (col 14)
  mawk '$14 == "TRUE"' $LOCI | sort -k1,2 | mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq \
    > CA.${TIER}.loci.bed
  bedtools intersect -wb -a CA.${TIER}.loci.bed -b $GENE_BED | \
    grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq \
    > Sig.loci.${TIER}.CA.LOC

  # SE (col 15)
  mawk '$15 == "TRUE"' $LOCI | sort -k1,2 | mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq \
    > SE.${TIER}.loci.bed
  bedtools intersect -wb -a SE.${TIER}.loci.bed -b $GENE_BED | \
    grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq \
    > Sig.loci.${TIER}.SE.LOC

  # CASE (col 13)
  mawk '$13 == "TRUE"' $LOCI | sort -k1,2 | mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq \
    > CASE.${TIER}.loci.bed
  bedtools intersect -wb -a CASE.${TIER}.loci.bed -b $GENE_BED | \
    grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq \
    > Sig.loci.${TIER}.CASE.LOC
done
```

``` bash
source activate CASE

GENE_BED=~/CASE/analysis/sorted.ref3.0.gene.bed

# ─── Pairwise and triple-treatment LOC lists per tier ────────────────────────
# Grouped.Loci column indices: col 11 = Sig.CASE, col 12 = Sig.CA, col 13 = Sig.SE
#                              col 14 = CHROM, col 15 = BP
# Core tier will have the most populated combination sets; Convergent and Private
# output may be empty in some categories, which is expected and interpretable.

for TIER in Core Conv Priv; do
  GLOCI="${TIER}.Significant.Grouped.Loci"

  # SE+CASE but not CA
  mawk '$11 == "TRUE" && $13 == "TRUE" && $12 == "FALSE"' $GLOCI | \
    mawk '{print $14 "\t" $15-1 "\t" $15}' | uniq > SE.CASE.${TIER}.loci.bed
  bedtools intersect -wb -a SE.CASE.${TIER}.loci.bed -b $GENE_BED | \
    grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq \
    > Sig.loci.${TIER}.SE.CASE.LOC

  # SE+CASE+CA
  mawk '$11 == "TRUE" && $13 == "TRUE" && $12 == "TRUE"' $GLOCI | \
    mawk '{print $14 "\t" $15-1 "\t" $15}' | uniq > SE.CASE.CA.${TIER}.loci.bed
  bedtools intersect -wb -a SE.CASE.CA.${TIER}.loci.bed -b $GENE_BED | \
    grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq \
    > Sig.loci.${TIER}.SE.CASE.CA.LOC

  # CA factor: CA-only + CA+CASE loci (all loci responding to acidification)
  mawk '$12 == "TRUE" && $11 == "FALSE" && $13 == "FALSE"' $GLOCI | \
    mawk '{print $14 "\t" $15-1 "\t" $15}' | uniq > CA.ONLY.${TIER}.loci.bed
  mawk '$12 == "TRUE" && $11 == "TRUE"  && $13 == "FALSE"' $GLOCI | \
    mawk '{print $14 "\t" $15-1 "\t" $15}' | uniq > CA.CASE.${TIER}.loci.bed
  cat CA.ONLY.${TIER}.loci.bed CA.CASE.${TIER}.loci.bed | sort | uniq > CA.factor.${TIER}.loci.bed
  bedtools intersect -wb -a CA.factor.${TIER}.loci.bed -b $GENE_BED | \
    grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq \
    > Sig.loci.${TIER}.CA.factor.LOC

  # SE factor: SE-only + SE+CASE loci (all loci responding to sewage effluent)
  mawk '$13 == "TRUE" && $11 == "FALSE" && $12 == "FALSE"' $GLOCI | \
    mawk '{print $14 "\t" $15-1 "\t" $15}' | uniq > SE.ONLY.${TIER}.loci.bed
  mawk '$13 == "TRUE" && $11 == "TRUE"  && $12 == "FALSE"' $GLOCI | \
    mawk '{print $14 "\t" $15-1 "\t" $15}' | uniq > SE.CASE.ONLY.${TIER}.loci.bed
  cat SE.ONLY.${TIER}.loci.bed SE.CASE.ONLY.${TIER}.loci.bed | sort | uniq > SE.factor.${TIER}.loci.bed
  bedtools intersect -wb -a SE.factor.${TIER}.loci.bed -b $GENE_BED | \
    grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq \
    > Sig.loci.${TIER}.SE.factor.LOC
done
```

### Outlier-overlap LOC lists (for the Fig 4 Venn)

Pairwise and triple gene-list intersections used by the combined
SNP/gene Venn.

``` bash
source activate CASE
cat Sig.loci.1.CA.LOC Sig.loci.1.SE.LOC   | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.1.CA.SE.LOC
cat Sig.loci.1.CA.LOC Sig.loci.1.CASE.LOC | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.1.CA.CASE.LOC
cat Sig.loci.1.SE.LOC Sig.loci.1.CASE.LOC | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.1.SE.CASE.LOC
cat Sig.loci.1.SE.LOC Sig.loci.1.CASE.LOC Sig.loci.1.CA.LOC | sort | uniq -c | mawk '$1 > 2' | mawk '{print $2}' > Sig.loci.1.SE.CASE.CA.LOC
```

## Gene Ontology Enrichment Analysis with topGO

### Download and prepare GO annotation

``` bash
source activate CASE
mkdir -p ./GO
cd ./GO

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_gene_ontology.gaf.gz

gunzip GCF_002022765.2_C_virginica-3.0_gene_ontology.gaf.gz

mv GCF_002022765.2_C_virginica-3.0_gene_ontology.gaf gene_ontology.gaf
```

### Create gene-to-GO mapping file

Making a tab file with gene ID and GO IDs. Each gene has one row with
all associated GO terms separated by semicolons.

``` bash
cd ./GO

awk 'NR>9 {print $3 "\t" $5}' gene_ontology.gaf | 
sed '1i gene_name\tGeneOntologyIDs' | 
sort | awk '
NR==1 {print; next}
$1!=prev {
if (NR>2) print prev "\t" goids;
prev=$1;
goids=$2;
next
}
{
goids=goids "; " $2
}
END {print prev "\t" goids}
' > gene2go_full.tab
```

### Extract GO terms for candidate genes (by stressor category)

``` bash
cd ./GO

# All outlier loci
grep -F -f ../Sig.loci.TOTAL.CASE.LOC gene2go_full.tab > ALL_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' ALL_cand_genes_GO.tab > ALL_cand_genes_GO_head.tab

# CA stressor loci
grep -F -f ../Sig.loci.1.CA.LOC gene2go_full.tab > CA_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' CA_cand_genes_GO.tab > CA_cand_genes_GO_head.tab

# SE stressor loci
grep -F -f ../Sig.loci.1.SE.LOC gene2go_full.tab > SE_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' SE_cand_genes_GO.tab > SE_cand_genes_GO_head.tab

# CASE (CA+SE combined) stressor loci
grep -F -f ../Sig.loci.1.CASE.LOC gene2go_full.tab > CASE_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' CASE_cand_genes_GO.tab > CASE_cand_genes_GO_head.tab
```

### Extract GO terms for tier-specific candidate genes

``` bash
cd ./GO

# --- Core loci (Tiers A–C, multi-block) ---
grep -F -f ../Sig.loci.Core.TOTAL.LOC  gene2go_full.tab > Core_ALL_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Core_ALL_cand_genes_GO.tab  > Core_ALL_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Core.CA.LOC     gene2go_full.tab > Core_CA_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Core_CA_cand_genes_GO.tab   > Core_CA_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Core.SE.LOC     gene2go_full.tab > Core_SE_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Core_SE_cand_genes_GO.tab   > Core_SE_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Core.CASE.LOC   gene2go_full.tab > Core_CASE_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Core_CASE_cand_genes_GO.tab > Core_CASE_cand_genes_GO_head.tab

# --- Private loci ---
grep -F -f ../Sig.loci.Priv.TOTAL.LOC    gene2go_full.tab > Priv_ALL_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Priv_ALL_cand_genes_GO.tab    > Priv_ALL_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Priv.CA.LOC       gene2go_full.tab > Priv_CA_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Priv_CA_cand_genes_GO.tab     > Priv_CA_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Priv.SE.LOC       gene2go_full.tab > Priv_SE_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Priv_SE_cand_genes_GO.tab     > Priv_SE_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Priv.CASE.LOC     gene2go_full.tab > Priv_CASE_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Priv_CASE_cand_genes_GO.tab   > Priv_CASE_cand_genes_GO_head.tab

# --- Convergent loci (per-block + all-blocks cross-validated) ---
# Note: these files may be empty if convergent loci map to no annotated genes;
# downstream R code handles empty input gracefully via run_topGO_analysis().
grep -F -f ../Sig.loci.Conv.TOTAL.LOC  gene2go_full.tab > Conv_ALL_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Conv_ALL_cand_genes_GO.tab  > Conv_ALL_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Conv.CA.LOC     gene2go_full.tab > Conv_CA_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Conv_CA_cand_genes_GO.tab   > Conv_CA_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Conv.SE.LOC     gene2go_full.tab > Conv_SE_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Conv_SE_cand_genes_GO.tab   > Conv_SE_cand_genes_GO_head.tab

grep -F -f ../Sig.loci.Conv.CASE.LOC   gene2go_full.tab > Conv_CASE_cand_genes_GO.tab
awk 'BEGIN {print "gene_name\tGeneOntologyIDs"} {print}' Conv_CASE_cand_genes_GO.tab > Conv_CASE_cand_genes_GO_head.tab
```

### Load libraries and gene-to-GO data

``` r
all_gene2go <- read.delim("./GO/gene2go_full.tab", sep="\t")
gene2go_topgo <- readMappings("./GO/gene2go_full.tab", IDsep=";", sep="\t")
all_genes <- as.character(unique(all_gene2go$gene_name))

topDiffGenes <- function(allScore) {
  return(allScore < 0.05)
}
```

### Define GO enrichment function

This function runs topGO for all three ontologies (BP, CC, MF), extracts
significant GO terms, maps them back to candidate genes, and returns the
combined results.

``` r
run_topGO_analysis <- function(cand_gene_file, all_genes, gene2go_topgo, label) {
  
  # Load candidate gene-to-GO mapping
  cand_gene2go <- read.delim(cand_gene_file, sep="\t")
  cand_genes <- as.character(unique(cand_gene2go$gene_name))
  
  cat(paste0("\n--- ", label, " ---\n"))
  cat(paste0("Number of candidate genes with GO annotations: ", length(cand_genes), "\n"))
  
  # Create gene list
  GeneList <- factor(as.integer(all_genes %in% cand_genes))
  names(GeneList) <- all_genes

  # topGO's factor path requires both classes present (a 0/1 split). If no
  # candidate gene falls in the annotated universe (or, degenerately, all do),
  # the factor collapses to one level and new("topGOdata") errors with
  # "allGenes must be a factor with 2 levels". Skip gracefully and return NULL
  # (handled downstream), printing the count so the cause is visible on knit.
  n_in <- sum(all_genes %in% cand_genes)
  cat(paste0("Candidate genes in annotated universe: ", n_in, " of ", length(all_genes), "\n"))
  if (nlevels(GeneList) < 2) {
    cat(paste0("Skipping ", label, ": gene selection has <2 levels (",
               n_in, " candidates in universe) -- not testable.\n"))
    return(NULL)
  }

  # Separate GO terms for gene-level mapping later
  cand_gene2go_sep <- cand_gene2go %>%
    separate_rows(GeneOntologyIDs, sep = ";")
  cand_gene2go_sep$GeneOntologyIDs <- trimws(cand_gene2go_sep$GeneOntologyIDs)
  
  results_list <- list()
  
  # --- Biological Processes ---
  GO_BP <- new("topGOdata", ontology="BP", gene2GO=gene2go_topgo, 
               allGenes=GeneList, annot = annFUN.gene2GO, geneSel=topDiffGenes)
  GO_BP_FE <- runTest(GO_BP, algorithm="weight01", statistic="fisher")
  GO_BP_table <- GenTable(GO_BP, Fisher = GO_BP_FE, orderBy = "Fisher", topNodes = 100, numChar = 100)
  GO_BP_table$Fisher <- as.numeric(GO_BP_table$Fisher)
  GO_BP_sig <- GO_BP_table[GO_BP_table$Fisher < 0.05,]
  
  if(nrow(GO_BP_sig) > 0) {
    GO_BP_sig$GO.ID <- trimws(GO_BP_sig$GO.ID)
    GO_BP_sig_gene <- cand_gene2go_sep %>%
      left_join(GO_BP_sig, by = c("GeneOntologyIDs" = "GO.ID")) %>%
      na.omit()
    GO_BP_sig_gene$ontology <- "Biological Processes"
    results_list[["BP"]] <- GO_BP_sig_gene
  }
  
  # --- Cellular Components ---
  GO_CC <- new("topGOdata", ontology="CC", gene2GO=gene2go_topgo, 
               allGenes=GeneList, annot = annFUN.gene2GO, geneSel=topDiffGenes)
  GO_CC_FE <- runTest(GO_CC, algorithm="weight01", statistic="fisher")
  GO_CC_table <- GenTable(GO_CC, Fisher = GO_CC_FE, orderBy = "Fisher", topNodes = 100, numChar = 100)
  GO_CC_table$Fisher <- as.numeric(GO_CC_table$Fisher)
  GO_CC_sig <- GO_CC_table[GO_CC_table$Fisher < 0.05,]
  
  if(nrow(GO_CC_sig) > 0) {
    GO_CC_sig$GO.ID <- trimws(GO_CC_sig$GO.ID)
    GO_CC_sig_gene <- cand_gene2go_sep %>%
      left_join(GO_CC_sig, by = c("GeneOntologyIDs" = "GO.ID")) %>%
      na.omit()
    GO_CC_sig_gene$ontology <- "Cellular Components"
    results_list[["CC"]] <- GO_CC_sig_gene
  }
  
  # --- Molecular Functions ---
  GO_MF <- new("topGOdata", ontology="MF", gene2GO=gene2go_topgo, 
               allGenes=GeneList, annot = annFUN.gene2GO, geneSel=topDiffGenes)
  GO_MF_FE <- runTest(GO_MF, algorithm="weight01", statistic="fisher")
  GO_MF_table <- GenTable(GO_MF, Fisher = GO_MF_FE, orderBy = "Fisher", topNodes = 100, numChar = 100)
  GO_MF_table$Fisher <- as.numeric(GO_MF_table$Fisher)
  GO_MF_sig <- GO_MF_table[GO_MF_table$Fisher < 0.05,]
  
  if(nrow(GO_MF_sig) > 0) {
    GO_MF_sig$GO.ID <- trimws(GO_MF_sig$GO.ID)
    GO_MF_sig_gene <- cand_gene2go_sep %>%
      left_join(GO_MF_sig, by = c("GeneOntologyIDs" = "GO.ID")) %>%
      na.omit()
    GO_MF_sig_gene$ontology <- "Molecular Functions"
    results_list[["MF"]] <- GO_MF_sig_gene
  }
  
  # Combine all ontologies
  if(length(results_list) > 0) {
    combined <- do.call(rbind, results_list)
    combined <- combined %>% mutate(prop.sig.genes = Significant/Annotated)
    combined <- unique(combined)
    combined$stressor <- label
    return(combined)
  } else {
    cat(paste0("No significant GO terms found for ", label, "\n"))
    return(NULL)
  }
}
```

### Run GO enrichment for all four categories

``` r
# All outlier loci
GO_ALL <- run_topGO_analysis("./GO/ALL_cand_genes_GO_head.tab", all_genes, gene2go_topgo, "All Outliers")
```

    ## 
    ## --- All Outliers ---
    ## Number of candidate genes with GO annotations: 1289
    ## Candidate genes in annotated universe: 1289 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 2376 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  18 nodes to be scored   (0 eliminated genes)

    ## 
    ##   Level 14:  34 nodes to be scored   (7 eliminated genes)

    ## 
    ##   Level 13:  58 nodes to be scored   (81 eliminated genes)

    ## 
    ##   Level 12:  107 nodes to be scored  (304 eliminated genes)

    ## 
    ##   Level 11:  156 nodes to be scored  (827 eliminated genes)

    ## 
    ##   Level 10:  230 nodes to be scored  (2721 eliminated genes)

    ## 
    ##   Level 9:   289 nodes to be scored  (4292 eliminated genes)

    ## 
    ##   Level 8:   312 nodes to be scored  (6029 eliminated genes)

    ## 
    ##   Level 7:   354 nodes to be scored  (7564 eliminated genes)

    ## 
    ##   Level 6:   336 nodes to be scored  (9944 eliminated genes)

    ## 
    ##   Level 5:   262 nodes to be scored  (11450 eliminated genes)

    ## 
    ##   Level 4:   138 nodes to be scored  (13167 eliminated genes)

    ## 
    ##   Level 3:   61 nodes to be scored   (14140 eliminated genes)

    ## 
    ##   Level 2:   16 nodes to be scored   (14588 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14781 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 603 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  23 nodes to be scored   (43 eliminated genes)

    ## 
    ##   Level 13:  33 nodes to be scored   (74 eliminated genes)

    ## 
    ##   Level 12:  55 nodes to be scored   (247 eliminated genes)

    ## 
    ##   Level 11:  75 nodes to be scored   (544 eliminated genes)

    ## 
    ##   Level 10:  78 nodes to be scored   (1217 eliminated genes)

    ## 
    ##   Level 9:   58 nodes to be scored   (3124 eliminated genes)

    ## 
    ##   Level 8:   65 nodes to be scored   (3726 eliminated genes)

    ## 
    ##   Level 7:   37 nodes to be scored   (4282 eliminated genes)

    ## 
    ##   Level 6:   44 nodes to be scored   (7516 eliminated genes)

    ## 
    ##   Level 5:   35 nodes to be scored   (7809 eliminated genes)

    ## 
    ##   Level 4:   47 nodes to be scored   (9200 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (10237 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (11050 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11114 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 935 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (29 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  24 nodes to be scored   (45 eliminated genes)

    ## 
    ##   Level 9:   54 nodes to be scored   (105 eliminated genes)

    ## 
    ##   Level 8:   89 nodes to be scored   (904 eliminated genes)

    ## 
    ##   Level 7:   162 nodes to be scored  (3325 eliminated genes)

    ## 
    ##   Level 6:   209 nodes to be scored  (4136 eliminated genes)

    ## 
    ##   Level 5:   172 nodes to be scored  (7989 eliminated genes)

    ## 
    ##   Level 4:   146 nodes to be scored  (11102 eliminated genes)

    ## 
    ##   Level 3:   48 nodes to be scored   (13887 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (15364 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19529 eliminated genes)

``` r
# CA stressor
GO_CA <- run_topGO_analysis("./GO/CA_cand_genes_GO_head.tab", all_genes, gene2go_topgo, "CA")
```

    ## 
    ## --- CA ---
    ## Number of candidate genes with GO annotations: 282
    ## Candidate genes in annotated universe: 282 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 1060 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  9 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  16 nodes to be scored   (43 eliminated genes)

    ## 
    ##   Level 12:  32 nodes to be scored   (198 eliminated genes)

    ## 
    ##   Level 11:  48 nodes to be scored   (452 eliminated genes)

    ## 
    ##   Level 10:  78 nodes to be scored   (2118 eliminated genes)

    ## 
    ##   Level 9:   123 nodes to be scored  (3453 eliminated genes)

    ## 
    ##   Level 8:   127 nodes to be scored  (4485 eliminated genes)

    ## 
    ##   Level 7:   152 nodes to be scored  (6461 eliminated genes)

    ## 
    ##   Level 6:   164 nodes to be scored  (9082 eliminated genes)

    ## 
    ##   Level 5:   158 nodes to be scored  (10468 eliminated genes)

    ## 
    ##   Level 4:   90 nodes to be scored   (12623 eliminated genes)

    ## 
    ##   Level 3:   44 nodes to be scored   (13798 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (14431 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14712 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 301 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  6 nodes to be scored    (8 eliminated genes)

    ## 
    ##   Level 13:  10 nodes to be scored   (18 eliminated genes)

    ## 
    ##   Level 12:  22 nodes to be scored   (104 eliminated genes)

    ## 
    ##   Level 11:  32 nodes to be scored   (252 eliminated genes)

    ## 
    ##   Level 10:  36 nodes to be scored   (945 eliminated genes)

    ## 
    ##   Level 9:   33 nodes to be scored   (2733 eliminated genes)

    ## 
    ##   Level 8:   32 nodes to be scored   (3434 eliminated genes)

    ## 
    ##   Level 7:   25 nodes to be scored   (4071 eliminated genes)

    ## 
    ##   Level 6:   21 nodes to be scored   (7219 eliminated genes)

    ## 
    ##   Level 5:   19 nodes to be scored   (7779 eliminated genes)

    ## 
    ##   Level 4:   30 nodes to be scored   (9153 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (10072 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (10977 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11113 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 400 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  8 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 9:   24 nodes to be scored   (78 eliminated genes)

    ## 
    ##   Level 8:   32 nodes to be scored   (714 eliminated genes)

    ## 
    ##   Level 7:   55 nodes to be scored   (2813 eliminated genes)

    ## 
    ##   Level 6:   75 nodes to be scored   (3478 eliminated genes)

    ## 
    ##   Level 5:   79 nodes to be scored   (6669 eliminated genes)

    ## 
    ##   Level 4:   72 nodes to be scored   (9468 eliminated genes)

    ## 
    ##   Level 3:   32 nodes to be scored   (12788 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (14491 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19323 eliminated genes)

``` r
# SE stressor
GO_SE <- run_topGO_analysis("./GO/SE_cand_genes_GO_head.tab", all_genes, gene2go_topgo, "SE")
```

    ## 
    ## --- SE ---
    ## Number of candidate genes with GO annotations: 559
    ## Candidate genes in annotated universe: 559 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 1617 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  10 nodes to be scored   (0 eliminated genes)

    ## 
    ##   Level 14:  21 nodes to be scored   (5 eliminated genes)

    ## 
    ##   Level 13:  32 nodes to be scored   (36 eliminated genes)

    ## 
    ##   Level 12:  60 nodes to be scored   (244 eliminated genes)

    ## 
    ##   Level 11:  90 nodes to be scored   (607 eliminated genes)

    ## 
    ##   Level 10:  146 nodes to be scored  (2331 eliminated genes)

    ## 
    ##   Level 9:   195 nodes to be scored  (3667 eliminated genes)

    ## 
    ##   Level 8:   198 nodes to be scored  (5352 eliminated genes)

    ## 
    ##   Level 7:   243 nodes to be scored  (7062 eliminated genes)

    ## 
    ##   Level 6:   236 nodes to be scored  (9479 eliminated genes)

    ## 
    ##   Level 5:   205 nodes to be scored  (11110 eliminated genes)

    ## 
    ##   Level 4:   109 nodes to be scored  (12908 eliminated genes)

    ## 
    ##   Level 3:   52 nodes to be scored   (14038 eliminated genes)

    ## 
    ##   Level 2:   16 nodes to be scored   (14547 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14768 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 452 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  9 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 13:  20 nodes to be scored   (61 eliminated genes)

    ## 
    ##   Level 12:  32 nodes to be scored   (153 eliminated genes)

    ## 
    ##   Level 11:  56 nodes to be scored   (406 eliminated genes)

    ## 
    ##   Level 10:  61 nodes to be scored   (998 eliminated genes)

    ## 
    ##   Level 9:   52 nodes to be scored   (2893 eliminated genes)

    ## 
    ##   Level 8:   53 nodes to be scored   (3639 eliminated genes)

    ## 
    ##   Level 7:   29 nodes to be scored   (4240 eliminated genes)

    ## 
    ##   Level 6:   25 nodes to be scored   (7460 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (7804 eliminated genes)

    ## 
    ##   Level 4:   40 nodes to be scored   (9178 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (10233 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (11046 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11113 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 567 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (29 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (43 eliminated genes)

    ## 
    ##   Level 9:   28 nodes to be scored   (86 eliminated genes)

    ## 
    ##   Level 8:   44 nodes to be scored   (687 eliminated genes)

    ## 
    ##   Level 7:   76 nodes to be scored   (2948 eliminated genes)

    ## 
    ##   Level 6:   125 nodes to be scored  (3846 eliminated genes)

    ## 
    ##   Level 5:   113 nodes to be scored  (7114 eliminated genes)

    ## 
    ##   Level 4:   98 nodes to be scored   (10560 eliminated genes)

    ## 
    ##   Level 3:   45 nodes to be scored   (13568 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (15048 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19520 eliminated genes)

``` r
# CASE (CA+SE combined) stressor
GO_CASE <- run_topGO_analysis("./GO/CASE_cand_genes_GO_head.tab", all_genes, gene2go_topgo, "CASE")
```

    ## 
    ## --- CASE ---
    ## Number of candidate genes with GO annotations: 942
    ## Candidate genes in annotated universe: 942 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 2087 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  22 nodes to be scored   (2 eliminated genes)

    ## 
    ##   Level 13:  46 nodes to be scored   (53 eliminated genes)

    ## 
    ##   Level 12:  90 nodes to be scored   (284 eliminated genes)

    ## 
    ##   Level 11:  130 nodes to be scored  (801 eliminated genes)

    ## 
    ##   Level 10:  191 nodes to be scored  (2611 eliminated genes)

    ## 
    ##   Level 9:   249 nodes to be scored  (4175 eliminated genes)

    ## 
    ##   Level 8:   275 nodes to be scored  (5863 eliminated genes)

    ## 
    ##   Level 7:   315 nodes to be scored  (7438 eliminated genes)

    ## 
    ##   Level 6:   310 nodes to be scored  (9887 eliminated genes)

    ## 
    ##   Level 5:   246 nodes to be scored  (11382 eliminated genes)

    ## 
    ##   Level 4:   130 nodes to be scored  (13155 eliminated genes)

    ## 
    ##   Level 3:   57 nodes to be scored   (14086 eliminated genes)

    ## 
    ##   Level 2:   16 nodes to be scored   (14538 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14746 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 525 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  22 nodes to be scored   (37 eliminated genes)

    ## 
    ##   Level 13:  28 nodes to be scored   (64 eliminated genes)

    ## 
    ##   Level 12:  48 nodes to be scored   (246 eliminated genes)

    ## 
    ##   Level 11:  68 nodes to be scored   (501 eliminated genes)

    ## 
    ##   Level 10:  67 nodes to be scored   (1179 eliminated genes)

    ## 
    ##   Level 9:   49 nodes to be scored   (3092 eliminated genes)

    ## 
    ##   Level 8:   51 nodes to be scored   (3630 eliminated genes)

    ## 
    ##   Level 7:   31 nodes to be scored   (4183 eliminated genes)

    ## 
    ##   Level 6:   40 nodes to be scored   (7372 eliminated genes)

    ## 
    ##   Level 5:   31 nodes to be scored   (7736 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (9189 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (10235 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (11033 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11114 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 809 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (29 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  18 nodes to be scored   (40 eliminated genes)

    ## 
    ##   Level 9:   46 nodes to be scored   (93 eliminated genes)

    ## 
    ##   Level 8:   74 nodes to be scored   (689 eliminated genes)

    ## 
    ##   Level 7:   135 nodes to be scored  (3288 eliminated genes)

    ## 
    ##   Level 6:   179 nodes to be scored  (4030 eliminated genes)

    ## 
    ##   Level 5:   151 nodes to be scored  (7759 eliminated genes)

    ## 
    ##   Level 4:   135 nodes to be scored  (10944 eliminated genes)

    ## 
    ##   Level 3:   43 nodes to be scored   (13742 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (15336 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19509 eliminated genes)

### Run GO enrichment for Core loci

``` r
# Core loci: multi-block, Tiers A–C; highest-confidence locus set
GO_Core_ALL  <- run_topGO_analysis("./GO/Core_ALL_cand_genes_GO_head.tab",  all_genes, gene2go_topgo, "Core: All Outliers")
```

    ## 
    ## --- Core: All Outliers ---
    ## Number of candidate genes with GO annotations: 1009
    ## Candidate genes in annotated universe: 1009 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 2093 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  18 nodes to be scored   (0 eliminated genes)

    ## 
    ##   Level 14:  30 nodes to be scored   (7 eliminated genes)

    ## 
    ##   Level 13:  52 nodes to be scored   (81 eliminated genes)

    ## 
    ##   Level 12:  96 nodes to be scored   (285 eliminated genes)

    ## 
    ##   Level 11:  128 nodes to be scored  (797 eliminated genes)

    ## 
    ##   Level 10:  198 nodes to be scored  (2649 eliminated genes)

    ## 
    ##   Level 9:   250 nodes to be scored  (4188 eliminated genes)

    ## 
    ##   Level 8:   265 nodes to be scored  (5868 eliminated genes)

    ## 
    ##   Level 7:   303 nodes to be scored  (7429 eliminated genes)

    ## 
    ##   Level 6:   296 nodes to be scored  (9742 eliminated genes)

    ## 
    ##   Level 5:   245 nodes to be scored  (11298 eliminated genes)

    ## 
    ##   Level 4:   133 nodes to be scored  (13109 eliminated genes)

    ## 
    ##   Level 3:   59 nodes to be scored   (14136 eliminated genes)

    ## 
    ##   Level 2:   15 nodes to be scored   (14588 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14781 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 541 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  17 nodes to be scored   (43 eliminated genes)

    ## 
    ##   Level 13:  29 nodes to be scored   (70 eliminated genes)

    ## 
    ##   Level 12:  45 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 11:  70 nodes to be scored   (491 eliminated genes)

    ## 
    ##   Level 10:  69 nodes to be scored   (1157 eliminated genes)

    ## 
    ##   Level 9:   55 nodes to be scored   (3105 eliminated genes)

    ## 
    ##   Level 8:   56 nodes to be scored   (3650 eliminated genes)

    ## 
    ##   Level 7:   33 nodes to be scored   (4260 eliminated genes)

    ## 
    ##   Level 6:   40 nodes to be scored   (7486 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (7809 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (9200 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (10236 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (11049 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11114 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 814 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (29 eliminated genes)

    ## 
    ##   Level 11:  6 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  20 nodes to be scored   (45 eliminated genes)

    ## 
    ##   Level 9:   48 nodes to be scored   (104 eliminated genes)

    ## 
    ##   Level 8:   73 nodes to be scored   (772 eliminated genes)

    ## 
    ##   Level 7:   131 nodes to be scored  (3245 eliminated genes)

    ## 
    ##   Level 6:   178 nodes to be scored  (4025 eliminated genes)

    ## 
    ##   Level 5:   155 nodes to be scored  (7727 eliminated genes)

    ## 
    ##   Level 4:   133 nodes to be scored  (10840 eliminated genes)

    ## 
    ##   Level 3:   46 nodes to be scored   (13829 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (15297 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19527 eliminated genes)

``` r
GO_Core_CA   <- run_topGO_analysis("./GO/Core_CA_cand_genes_GO_head.tab",   all_genes, gene2go_topgo, "Core: CA")
```

    ## 
    ## --- Core: CA ---
    ## Number of candidate genes with GO annotations: 183
    ## Candidate genes in annotated universe: 183 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 743 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  9 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  13 nodes to be scored   (43 eliminated genes)

    ## 
    ##   Level 12:  21 nodes to be scored   (198 eliminated genes)

    ## 
    ##   Level 11:  24 nodes to be scored   (384 eliminated genes)

    ## 
    ##   Level 10:  45 nodes to be scored   (1980 eliminated genes)

    ## 
    ##   Level 9:   83 nodes to be scored   (3251 eliminated genes)

    ## 
    ##   Level 8:   80 nodes to be scored   (4103 eliminated genes)

    ## 
    ##   Level 7:   101 nodes to be scored  (5770 eliminated genes)

    ## 
    ##   Level 6:   113 nodes to be scored  (7981 eliminated genes)

    ## 
    ##   Level 5:   117 nodes to be scored  (9819 eliminated genes)

    ## 
    ##   Level 4:   79 nodes to be scored   (12098 eliminated genes)

    ## 
    ##   Level 3:   40 nodes to be scored   (13738 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (14383 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14690 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 246 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  3 nodes to be scored    (8 eliminated genes)

    ## 
    ##   Level 13:  6 nodes to be scored    (18 eliminated genes)

    ## 
    ##   Level 12:  17 nodes to be scored   (84 eliminated genes)

    ## 
    ##   Level 11:  27 nodes to be scored   (160 eliminated genes)

    ## 
    ##   Level 10:  29 nodes to be scored   (666 eliminated genes)

    ## 
    ##   Level 9:   29 nodes to be scored   (2548 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (3209 eliminated genes)

    ## 
    ##   Level 7:   19 nodes to be scored   (4011 eliminated genes)

    ## 
    ##   Level 6:   19 nodes to be scored   (7094 eliminated genes)

    ## 
    ##   Level 5:   14 nodes to be scored   (7722 eliminated genes)

    ## 
    ##   Level 4:   25 nodes to be scored   (9150 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (10065 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (10960 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11113 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 337 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (77 eliminated genes)

    ## 
    ##   Level 8:   27 nodes to be scored   (644 eliminated genes)

    ## 
    ##   Level 7:   48 nodes to be scored   (2691 eliminated genes)

    ## 
    ##   Level 6:   64 nodes to be scored   (3088 eliminated genes)

    ## 
    ##   Level 5:   65 nodes to be scored   (6267 eliminated genes)

    ## 
    ##   Level 4:   59 nodes to be scored   (9370 eliminated genes)

    ## 
    ##   Level 3:   30 nodes to be scored   (12645 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (14249 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19282 eliminated genes)

``` r
GO_Core_SE   <- run_topGO_analysis("./GO/Core_SE_cand_genes_GO_head.tab",   all_genes, gene2go_topgo, "Core: SE")
```

    ## 
    ## --- Core: SE ---
    ## Number of candidate genes with GO annotations: 437
    ## Candidate genes in annotated universe: 437 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 1372 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  10 nodes to be scored   (0 eliminated genes)

    ## 
    ##   Level 14:  18 nodes to be scored   (5 eliminated genes)

    ## 
    ##   Level 13:  27 nodes to be scored   (36 eliminated genes)

    ## 
    ##   Level 12:  52 nodes to be scored   (195 eliminated genes)

    ## 
    ##   Level 11:  79 nodes to be scored   (524 eliminated genes)

    ## 
    ##   Level 10:  127 nodes to be scored  (2267 eliminated genes)

    ## 
    ##   Level 9:   164 nodes to be scored  (3621 eliminated genes)

    ## 
    ##   Level 8:   161 nodes to be scored  (5195 eliminated genes)

    ## 
    ##   Level 7:   197 nodes to be scored  (6849 eliminated genes)

    ## 
    ##   Level 6:   196 nodes to be scored  (9297 eliminated genes)

    ## 
    ##   Level 5:   176 nodes to be scored  (10977 eliminated genes)

    ## 
    ##   Level 4:   99 nodes to be scored   (12810 eliminated genes)

    ## 
    ##   Level 3:   48 nodes to be scored   (13981 eliminated genes)

    ## 
    ##   Level 2:   14 nodes to be scored   (14547 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14753 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 398 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  7 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 13:  14 nodes to be scored   (39 eliminated genes)

    ## 
    ##   Level 12:  25 nodes to be scored   (126 eliminated genes)

    ## 
    ##   Level 11:  50 nodes to be scored   (266 eliminated genes)

    ## 
    ##   Level 10:  55 nodes to be scored   (925 eliminated genes)

    ## 
    ##   Level 9:   46 nodes to be scored   (2843 eliminated genes)

    ## 
    ##   Level 8:   44 nodes to be scored   (3547 eliminated genes)

    ## 
    ##   Level 7:   26 nodes to be scored   (4159 eliminated genes)

    ## 
    ##   Level 6:   25 nodes to be scored   (7251 eliminated genes)

    ## 
    ##   Level 5:   27 nodes to be scored   (7783 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (9178 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (10233 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (11045 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11112 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 516 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (29 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (43 eliminated genes)

    ## 
    ##   Level 9:   28 nodes to be scored   (85 eliminated genes)

    ## 
    ##   Level 8:   42 nodes to be scored   (613 eliminated genes)

    ## 
    ##   Level 7:   68 nodes to be scored   (2948 eliminated genes)

    ## 
    ##   Level 6:   107 nodes to be scored  (3845 eliminated genes)

    ## 
    ##   Level 5:   102 nodes to be scored  (7077 eliminated genes)

    ## 
    ##   Level 4:   90 nodes to be scored   (10384 eliminated genes)

    ## 
    ##   Level 3:   44 nodes to be scored   (13471 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (14932 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19520 eliminated genes)

``` r
GO_Core_CASE <- run_topGO_analysis("./GO/Core_CASE_cand_genes_GO_head.tab", all_genes, gene2go_topgo, "Core: CASE")
```

    ## 
    ## --- Core: CASE ---
    ## Number of candidate genes with GO annotations: 688
    ## Candidate genes in annotated universe: 688 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 1775 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  18 nodes to be scored   (2 eliminated genes)

    ## 
    ##   Level 13:  40 nodes to be scored   (42 eliminated genes)

    ## 
    ##   Level 12:  80 nodes to be scored   (265 eliminated genes)

    ## 
    ##   Level 11:  103 nodes to be scored  (752 eliminated genes)

    ## 
    ##   Level 10:  157 nodes to be scored  (2540 eliminated genes)

    ## 
    ##   Level 9:   205 nodes to be scored  (4063 eliminated genes)

    ## 
    ##   Level 8:   218 nodes to be scored  (5635 eliminated genes)

    ## 
    ##   Level 7:   258 nodes to be scored  (7158 eliminated genes)

    ## 
    ##   Level 6:   265 nodes to be scored  (9633 eliminated genes)

    ## 
    ##   Level 5:   228 nodes to be scored  (11190 eliminated genes)

    ## 
    ##   Level 4:   124 nodes to be scored  (13080 eliminated genes)

    ## 
    ##   Level 3:   55 nodes to be scored   (14062 eliminated genes)

    ## 
    ##   Level 2:   15 nodes to be scored   (14520 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14746 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 467 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  17 nodes to be scored   (37 eliminated genes)

    ## 
    ##   Level 13:  25 nodes to be scored   (64 eliminated genes)

    ## 
    ##   Level 12:  40 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 11:  61 nodes to be scored   (464 eliminated genes)

    ## 
    ##   Level 10:  55 nodes to be scored   (1139 eliminated genes)

    ## 
    ##   Level 9:   48 nodes to be scored   (3060 eliminated genes)

    ## 
    ##   Level 8:   44 nodes to be scored   (3558 eliminated genes)

    ## 
    ##   Level 7:   28 nodes to be scored   (4165 eliminated genes)

    ## 
    ##   Level 6:   35 nodes to be scored   (7340 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (7736 eliminated genes)

    ## 
    ##   Level 4:   38 nodes to be scored   (9189 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (10234 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (11033 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11114 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 680 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (29 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  16 nodes to be scored   (40 eliminated genes)

    ## 
    ##   Level 9:   42 nodes to be scored   (93 eliminated genes)

    ## 
    ##   Level 8:   58 nodes to be scored   (627 eliminated genes)

    ## 
    ##   Level 7:   102 nodes to be scored  (3220 eliminated genes)

    ## 
    ##   Level 6:   141 nodes to be scored  (3916 eliminated genes)

    ## 
    ##   Level 5:   130 nodes to be scored  (7332 eliminated genes)

    ## 
    ##   Level 4:   122 nodes to be scored  (10593 eliminated genes)

    ## 
    ##   Level 3:   41 nodes to be scored   (13647 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (15261 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19507 eliminated genes)

### Run GO enrichment for Private loci

``` r
# Private loci: per-block outlier singletons and reclassified Convergent loci
GO_Priv_ALL  <- run_topGO_analysis("./GO/Priv_ALL_cand_genes_GO_head.tab",  all_genes, gene2go_topgo, "Private: All Outliers")
```

    ## 
    ## --- Private: All Outliers ---
    ## Number of candidate genes with GO annotations: 295
    ## Candidate genes in annotated universe: 295 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 1196 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  21 nodes to be scored   (24 eliminated genes)

    ## 
    ##   Level 12:  36 nodes to be scored   (178 eliminated genes)

    ## 
    ##   Level 11:  61 nodes to be scored   (511 eliminated genes)

    ## 
    ##   Level 10:  91 nodes to be scored   (2123 eliminated genes)

    ## 
    ##   Level 9:   138 nodes to be scored  (3403 eliminated genes)

    ## 
    ##   Level 8:   155 nodes to be scored  (4549 eliminated genes)

    ## 
    ##   Level 7:   178 nodes to be scored  (6468 eliminated genes)

    ## 
    ##   Level 6:   184 nodes to be scored  (9176 eliminated genes)

    ## 
    ##   Level 5:   173 nodes to be scored  (10750 eliminated genes)

    ## 
    ##   Level 4:   92 nodes to be scored   (12787 eliminated genes)

    ## 
    ##   Level 3:   43 nodes to be scored   (13912 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (14508 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14713 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 324 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  13 nodes to be scored   (5 eliminated genes)

    ## 
    ##   Level 12:  26 nodes to be scored   (73 eliminated genes)

    ## 
    ##   Level 11:  40 nodes to be scored   (308 eliminated genes)

    ## 
    ##   Level 10:  41 nodes to be scored   (773 eliminated genes)

    ## 
    ##   Level 9:   28 nodes to be scored   (2694 eliminated genes)

    ## 
    ##   Level 8:   35 nodes to be scored   (3466 eliminated genes)

    ## 
    ##   Level 7:   25 nodes to be scored   (3961 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (7213 eliminated genes)

    ## 
    ##   Level 5:   21 nodes to be scored   (7738 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (9142 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (10068 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (11032 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11114 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 452 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (28 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (44 eliminated genes)

    ## 
    ##   Level 8:   34 nodes to be scored   (577 eliminated genes)

    ## 
    ##   Level 7:   65 nodes to be scored   (2743 eliminated genes)

    ## 
    ##   Level 6:   91 nodes to be scored   (3562 eliminated genes)

    ## 
    ##   Level 5:   90 nodes to be scored   (6951 eliminated genes)

    ## 
    ##   Level 4:   94 nodes to be scored   (9949 eliminated genes)

    ## 
    ##   Level 3:   35 nodes to be scored   (13205 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (14852 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19456 eliminated genes)

``` r
GO_Priv_CA   <- run_topGO_analysis("./GO/Priv_CA_cand_genes_GO_head.tab",   all_genes, gene2go_topgo, "Private: CA")
```

    ## 
    ## --- Private: CA ---
    ## Number of candidate genes with GO annotations: 90
    ## Candidate genes in annotated universe: 90 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 598 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  14 nodes to be scored   (111 eliminated genes)

    ## 
    ##   Level 11:  28 nodes to be scored   (244 eliminated genes)

    ## 
    ##   Level 10:  37 nodes to be scored   (1867 eliminated genes)

    ## 
    ##   Level 9:   52 nodes to be scored   (2998 eliminated genes)

    ## 
    ##   Level 8:   63 nodes to be scored   (3790 eliminated genes)

    ## 
    ##   Level 7:   86 nodes to be scored   (5148 eliminated genes)

    ## 
    ##   Level 6:   97 nodes to be scored   (7680 eliminated genes)

    ## 
    ##   Level 5:   103 nodes to be scored  (9544 eliminated genes)

    ## 
    ##   Level 4:   60 nodes to be scored   (11983 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (13431 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (14231 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14671 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 177 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  9 nodes to be scored    (20 eliminated genes)

    ## 
    ##   Level 11:  18 nodes to be scored   (92 eliminated genes)

    ## 
    ##   Level 10:  20 nodes to be scored   (426 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (1975 eliminated genes)

    ## 
    ##   Level 8:   19 nodes to be scored   (2514 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (3397 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (6488 eliminated genes)

    ## 
    ##   Level 5:   13 nodes to be scored   (7441 eliminated genes)

    ## 
    ##   Level 4:   20 nodes to be scored   (9053 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (10058 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (10948 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11106 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 229 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (23 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (23 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (44 eliminated genes)

    ## 
    ##   Level 8:   15 nodes to be scored   (577 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (2640 eliminated genes)

    ## 
    ##   Level 6:   39 nodes to be scored   (3233 eliminated genes)

    ## 
    ##   Level 5:   45 nodes to be scored   (6110 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (8625 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (11890 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (13931 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19319 eliminated genes)

``` r
GO_Priv_SE   <- run_topGO_analysis("./GO/Priv_SE_cand_genes_GO_head.tab",   all_genes, gene2go_topgo, "Private: SE")
```

    ## 
    ## --- Private: SE ---
    ## Number of candidate genes with GO annotations: 101
    ## Candidate genes in annotated universe: 101 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 682 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  9 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  14 nodes to be scored   (129 eliminated genes)

    ## 
    ##   Level 11:  23 nodes to be scored   (291 eliminated genes)

    ## 
    ##   Level 10:  33 nodes to be scored   (1866 eliminated genes)

    ## 
    ##   Level 9:   59 nodes to be scored   (2910 eliminated genes)

    ## 
    ##   Level 8:   73 nodes to be scored   (3530 eliminated genes)

    ## 
    ##   Level 7:   102 nodes to be scored  (4674 eliminated genes)

    ## 
    ##   Level 6:   125 nodes to be scored  (6675 eliminated genes)

    ## 
    ##   Level 5:   124 nodes to be scored  (8606 eliminated genes)

    ## 
    ##   Level 4:   69 nodes to be scored   (12241 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (13583 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (14196 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14672 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 186 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  13 nodes to be scored   (13 eliminated genes)

    ## 
    ##   Level 11:  21 nodes to be scored   (118 eliminated genes)

    ## 
    ##   Level 10:  24 nodes to be scored   (516 eliminated genes)

    ## 
    ##   Level 9:   21 nodes to be scored   (1969 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (2750 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (3649 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (7093 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (7543 eliminated genes)

    ## 
    ##   Level 4:   20 nodes to be scored   (9137 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (10056 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (10971 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11107 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 221 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (23 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (23 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (44 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (577 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (2518 eliminated genes)

    ## 
    ##   Level 6:   43 nodes to be scored   (2825 eliminated genes)

    ## 
    ##   Level 5:   43 nodes to be scored   (5687 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (8090 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (11239 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (13676 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (18819 eliminated genes)

``` r
GO_Priv_CASE <- run_topGO_analysis("./GO/Priv_CASE_cand_genes_GO_head.tab", all_genes, gene2go_topgo, "Private: CASE")
```

    ## 
    ## --- Private: CASE ---
    ## Number of candidate genes with GO annotations: 210
    ## Candidate genes in annotated universe: 210 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 1058 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  17 nodes to be scored   (24 eliminated genes)

    ## 
    ##   Level 12:  30 nodes to be scored   (160 eliminated genes)

    ## 
    ##   Level 11:  51 nodes to be scored   (479 eliminated genes)

    ## 
    ##   Level 10:  82 nodes to be scored   (2061 eliminated genes)

    ## 
    ##   Level 9:   122 nodes to be scored  (3316 eliminated genes)

    ## 
    ##   Level 8:   130 nodes to be scored  (4455 eliminated genes)

    ## 
    ##   Level 7:   155 nodes to be scored  (6017 eliminated genes)

    ## 
    ##   Level 6:   165 nodes to be scored  (8821 eliminated genes)

    ## 
    ##   Level 5:   157 nodes to be scored  (10491 eliminated genes)

    ## 
    ##   Level 4:   85 nodes to be scored   (12540 eliminated genes)

    ## 
    ##   Level 3:   41 nodes to be scored   (13846 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (14448 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14691 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 281 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  10 nodes to be scored   (5 eliminated genes)

    ## 
    ##   Level 12:  22 nodes to be scored   (72 eliminated genes)

    ## 
    ##   Level 11:  37 nodes to be scored   (251 eliminated genes)

    ## 
    ##   Level 10:  37 nodes to be scored   (738 eliminated genes)

    ## 
    ##   Level 9:   24 nodes to be scored   (2534 eliminated genes)

    ## 
    ##   Level 8:   28 nodes to be scored   (3395 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (3895 eliminated genes)

    ## 
    ##   Level 6:   19 nodes to be scored   (7100 eliminated genes)

    ## 
    ##   Level 5:   18 nodes to be scored   (7700 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (9142 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (10066 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (10910 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11058 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 397 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (28 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (43 eliminated genes)

    ## 
    ##   Level 8:   30 nodes to be scored   (503 eliminated genes)

    ## 
    ##   Level 7:   56 nodes to be scored   (2506 eliminated genes)

    ## 
    ##   Level 6:   79 nodes to be scored   (3477 eliminated genes)

    ## 
    ##   Level 5:   79 nodes to be scored   (6858 eliminated genes)

    ## 
    ##   Level 4:   86 nodes to be scored   (9736 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (12846 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (14524 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (18903 eliminated genes)

### Run GO enrichment for Convergent loci

``` r
# Convergent loci: per-block + all-blocks cross-validated; may return NULL if too few annotated genes
GO_Conv_ALL  <- run_topGO_analysis("./GO/Conv_ALL_cand_genes_GO_head.tab",  all_genes, gene2go_topgo, "Convergent: All Outliers")
```

    ## 
    ## --- Convergent: All Outliers ---
    ## Number of candidate genes with GO annotations: 265
    ## Candidate genes in annotated universe: 265 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 1227 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  15 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 12:  36 nodes to be scored   (121 eliminated genes)

    ## 
    ##   Level 11:  64 nodes to be scored   (422 eliminated genes)

    ## 
    ##   Level 10:  100 nodes to be scored  (2079 eliminated genes)

    ## 
    ##   Level 9:   135 nodes to be scored  (3494 eliminated genes)

    ## 
    ##   Level 8:   145 nodes to be scored  (4911 eliminated genes)

    ## 
    ##   Level 7:   185 nodes to be scored  (6187 eliminated genes)

    ## 
    ##   Level 6:   196 nodes to be scored  (8787 eliminated genes)

    ## 
    ##   Level 5:   178 nodes to be scored  (10509 eliminated genes)

    ## 
    ##   Level 4:   102 nodes to be scored  (12793 eliminated genes)

    ## 
    ##   Level 3:   49 nodes to be scored   (13949 eliminated genes)

    ## 
    ##   Level 2:   15 nodes to be scored   (14477 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14746 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 309 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  9 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 13:  12 nodes to be scored   (49 eliminated genes)

    ## 
    ##   Level 12:  18 nodes to be scored   (135 eliminated genes)

    ## 
    ##   Level 11:  36 nodes to be scored   (394 eliminated genes)

    ## 
    ##   Level 10:  41 nodes to be scored   (840 eliminated genes)

    ## 
    ##   Level 9:   35 nodes to be scored   (2749 eliminated genes)

    ## 
    ##   Level 8:   32 nodes to be scored   (3417 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (4043 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (7163 eliminated genes)

    ## 
    ##   Level 5:   19 nodes to be scored   (7638 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (9064 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (10049 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (10904 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11107 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 399 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (0 eliminated genes)

    ## 
    ##   Level 8:   31 nodes to be scored   (522 eliminated genes)

    ## 
    ##   Level 7:   58 nodes to be scored   (2789 eliminated genes)

    ## 
    ##   Level 6:   92 nodes to be scored   (3559 eliminated genes)

    ## 
    ##   Level 5:   79 nodes to be scored   (6814 eliminated genes)

    ## 
    ##   Level 4:   75 nodes to be scored   (9809 eliminated genes)

    ## 
    ##   Level 3:   35 nodes to be scored   (13189 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (14715 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19425 eliminated genes)

``` r
GO_Conv_CA   <- run_topGO_analysis("./GO/Conv_CA_cand_genes_GO_head.tab",   all_genes, gene2go_topgo, "Convergent: CA")
```

    ## 
    ## --- Convergent: CA ---
    ## Number of candidate genes with GO annotations: 18
    ## Candidate genes in annotated universe: 18 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 289 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (86 eliminated genes)

    ## 
    ##   Level 11:  10 nodes to be scored   (105 eliminated genes)

    ## 
    ##   Level 10:  14 nodes to be scored   (1022 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (2546 eliminated genes)

    ## 
    ##   Level 8:   20 nodes to be scored   (2955 eliminated genes)

    ## 
    ##   Level 7:   37 nodes to be scored   (3189 eliminated genes)

    ## 
    ##   Level 6:   52 nodes to be scored   (4347 eliminated genes)

    ## 
    ##   Level 5:   55 nodes to be scored   (6522 eliminated genes)

    ## 
    ##   Level 4:   41 nodes to be scored   (9103 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (11454 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (13965 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14514 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 54 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (354 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (456 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (590 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (5179 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (7210 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (8684 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (10022 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (10845 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11048 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 82 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (460 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (989 eliminated genes)

    ## 
    ##   Level 6:   14 nodes to be scored   (1050 eliminated genes)

    ## 
    ##   Level 5:   15 nodes to be scored   (4728 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (6246 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (7941 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (8618 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (16625 eliminated genes)

``` r
GO_Conv_SE   <- run_topGO_analysis("./GO/Conv_SE_cand_genes_GO_head.tab",   all_genes, gene2go_topgo, "Convergent: SE")
```

    ## 
    ## --- Convergent: SE ---
    ## Number of candidate genes with GO annotations: 54
    ## Candidate genes in annotated universe: 54 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 454 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  7 nodes to be scored    (19 eliminated genes)

    ## 
    ##   Level 11:  15 nodes to be scored   (60 eliminated genes)

    ## 
    ##   Level 10:  29 nodes to be scored   (1147 eliminated genes)

    ## 
    ##   Level 9:   44 nodes to be scored   (2838 eliminated genes)

    ## 
    ##   Level 8:   46 nodes to be scored   (3614 eliminated genes)

    ## 
    ##   Level 7:   64 nodes to be scored   (4647 eliminated genes)

    ## 
    ##   Level 6:   77 nodes to be scored   (6966 eliminated genes)

    ## 
    ##   Level 5:   76 nodes to be scored   (8018 eliminated genes)

    ## 
    ##   Level 4:   52 nodes to be scored   (11645 eliminated genes)

    ## 
    ##   Level 3:   28 nodes to be scored   (13235 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (13929 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14149 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 137 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (22 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (84 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (149 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (469 eliminated genes)

    ## 
    ##   Level 9:   16 nodes to be scored   (1855 eliminated genes)

    ## 
    ##   Level 8:   17 nodes to be scored   (2425 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (3507 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (6761 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (7520 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (9016 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (10032 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (10137 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (10239 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 169 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   11 nodes to be scored   (460 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (2370 eliminated genes)

    ## 
    ##   Level 6:   36 nodes to be scored   (2885 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (6039 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (8221 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (10598 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (12531 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (17640 eliminated genes)

``` r
GO_Conv_CASE <- run_topGO_analysis("./GO/Conv_CASE_cand_genes_GO_head.tab", all_genes, gene2go_topgo, "Convergent: CASE")
```

    ## 
    ## --- Convergent: CASE ---
    ## Number of candidate genes with GO annotations: 218
    ## Candidate genes in annotated universe: 218 of 22614

    ## 
    ## Building most specific GOs .....

    ##  ( 2565 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 5240 GO terms and 11520 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 14789 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 1085 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  13 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 12:  32 nodes to be scored   (121 eliminated genes)

    ## 
    ##   Level 11:  54 nodes to be scored   (388 eliminated genes)

    ## 
    ##   Level 10:  79 nodes to be scored   (2049 eliminated genes)

    ## 
    ##   Level 9:   112 nodes to be scored  (3407 eliminated genes)

    ## 
    ##   Level 8:   126 nodes to be scored  (4674 eliminated genes)

    ## 
    ##   Level 7:   160 nodes to be scored  (6051 eliminated genes)

    ## 
    ##   Level 6:   176 nodes to be scored  (8597 eliminated genes)

    ## 
    ##   Level 5:   165 nodes to be scored  (10347 eliminated genes)

    ## 
    ##   Level 4:   97 nodes to be scored   (12770 eliminated genes)

    ## 
    ##   Level 3:   49 nodes to be scored   (13926 eliminated genes)

    ## 
    ##   Level 2:   15 nodes to be scored   (14475 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (14746 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 834 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1151 GO terms and 2236 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 13587 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 292 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 16:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 15:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  8 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 13:  12 nodes to be scored   (45 eliminated genes)

    ## 
    ##   Level 12:  18 nodes to be scored   (124 eliminated genes)

    ## 
    ##   Level 11:  33 nodes to be scored   (394 eliminated genes)

    ## 
    ##   Level 10:  38 nodes to be scored   (840 eliminated genes)

    ## 
    ##   Level 9:   31 nodes to be scored   (2629 eliminated genes)

    ## 
    ##   Level 8:   29 nodes to be scored   (3386 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (4021 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (7160 eliminated genes)

    ## 
    ##   Level 5:   19 nodes to be scored   (7638 eliminated genes)

    ## 
    ##   Level 4:   30 nodes to be scored   (9064 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (10049 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (10904 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (11107 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 1923 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 2502 GO terms and 3298 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 19601 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 374 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (0 eliminated genes)

    ## 
    ##   Level 8:   30 nodes to be scored   (522 eliminated genes)

    ## 
    ##   Level 7:   52 nodes to be scored   (2789 eliminated genes)

    ## 
    ##   Level 6:   81 nodes to be scored   (3558 eliminated genes)

    ## 
    ##   Level 5:   76 nodes to be scored   (6806 eliminated genes)

    ## 
    ##   Level 4:   72 nodes to be scored   (9724 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (12956 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (14550 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (19425 eliminated genes)

### Save results

``` r
if(!is.null(GO_ALL))  write.csv(GO_ALL,  "./GO/ALL_GO_en_sig_gene.csv", row.names = FALSE)
if(!is.null(GO_CA))   write.csv(GO_CA,   "./GO/CA_GO_en_sig_gene.csv", row.names = FALSE)
if(!is.null(GO_SE))   write.csv(GO_SE,   "./GO/SE_GO_en_sig_gene.csv", row.names = FALSE)
if(!is.null(GO_CASE)) write.csv(GO_CASE, "./GO/CASE_GO_en_sig_gene.csv", row.names = FALSE)
```

### Save tier-specific topGO results

``` r
# Core
if(!is.null(GO_Core_ALL))  write.csv(GO_Core_ALL,  "./GO/Core_ALL_GO_en_sig_gene.csv",  row.names = FALSE)
if(!is.null(GO_Core_CA))   write.csv(GO_Core_CA,   "./GO/Core_CA_GO_en_sig_gene.csv",   row.names = FALSE)
if(!is.null(GO_Core_SE))   write.csv(GO_Core_SE,   "./GO/Core_SE_GO_en_sig_gene.csv",   row.names = FALSE)
if(!is.null(GO_Core_CASE)) write.csv(GO_Core_CASE, "./GO/Core_CASE_GO_en_sig_gene.csv", row.names = FALSE)

# Private
if(!is.null(GO_Priv_ALL))  write.csv(GO_Priv_ALL,  "./GO/Priv_ALL_GO_en_sig_gene.csv",  row.names = FALSE)
if(!is.null(GO_Priv_CA))   write.csv(GO_Priv_CA,   "./GO/Priv_CA_GO_en_sig_gene.csv",   row.names = FALSE)
if(!is.null(GO_Priv_SE))   write.csv(GO_Priv_SE,   "./GO/Priv_SE_GO_en_sig_gene.csv",   row.names = FALSE)
if(!is.null(GO_Priv_CASE)) write.csv(GO_Priv_CASE, "./GO/Priv_CASE_GO_en_sig_gene.csv", row.names = FALSE)

# Convergent
if(!is.null(GO_Conv_ALL))  write.csv(GO_Conv_ALL,  "./GO/Conv_ALL_GO_en_sig_gene.csv",  row.names = FALSE)
if(!is.null(GO_Conv_CA))   write.csv(GO_Conv_CA,   "./GO/Conv_CA_GO_en_sig_gene.csv",   row.names = FALSE)
if(!is.null(GO_Conv_SE))   write.csv(GO_Conv_SE,   "./GO/Conv_SE_GO_en_sig_gene.csv",   row.names = FALSE)
if(!is.null(GO_Conv_CASE)) write.csv(GO_Conv_CASE, "./GO/Conv_CASE_GO_en_sig_gene.csv", row.names = FALSE)
```

## Data-driven GO clusters (mechanism themes for Fig 4D + Table S)

Groups the enriched GO terms (pooled across CA/SE/CASE) by overlap of
the candidate genes annotated to them (Jaccard \>= 0.10), finds
communities with `igraph` greedy modularity, and tags an interpretive
set of stress-relevant clusters. The clustering is data-driven and
reproducible; the stress selection/naming is a curated overlay anchored
on stable GO IDs. Builds `go_clusters_stress` (Fig 4D) and
`go_clusters_full` (Table S).

``` r
library(igraph)
go_list <- Filter(Negate(is.null), list(GO_CA, GO_SE, GO_CASE))
go_pool <- as.data.table(rbindlist(go_list, fill = TRUE))
setnames(go_pool, "GeneOntologyIDs", "GO")
go_pool[, Fisher := suppressWarnings(as.numeric(Fisher))]
go_pool <- go_pool[is.finite(Fisher)]

# one row per GO term: representative text, min Fisher, per-treatment enrichment, gene set
term_info  <- go_pool[, .(Term = Term[1L], Fisher = min(Fisher),
                          CA = as.integer(any(stressor=="CA")),
                          SE = as.integer(any(stressor=="SE")),
                          CASE = as.integer(any(stressor=="CASE"))), by = GO]
gene_by_go <- go_pool[, .(g = list(unique(gene_name))), by = GO]
term_info  <- merge(term_info, gene_by_go, by = "GO")
term_info[, ngenes := lengths(g)]
clG <- term_info[ngenes >= 2]                          # need >=2 genes to share

# gene-overlap Jaccard edge list (n ~ 180 terms; nested loop is fine)
THR <- 0.10
gl <- clG$g; ids <- clG$GO
from <- character(0); to <- character(0); w <- numeric(0)
if (length(ids) > 1L) for (i in 1:(length(ids)-1L)) {
  gi <- gl[[i]]
  for (j in (i+1L):length(ids)) {
    inter <- length(intersect(gi, gl[[j]]))
    if (inter == 0L) next
    jac <- inter / length(union(gi, gl[[j]]))
    if (jac >= THR) { from <- c(from, ids[i]); to <- c(to, ids[j]); w <- c(w, jac) }
  }
}
g <- graph_from_data_frame(data.frame(from = from, to = to, weight = w),
                           directed = FALSE, vertices = data.frame(name = ids))
comm <- cluster_fast_greedy(g, weights = E(g)$weight)
clG[, cluster_id := as.integer(membership(comm)[GO])]
clG <- clG[cluster_id %in% clG[, .N, by = cluster_id][N >= 3, cluster_id]]   # >=3 terms

# curated stress overlay: category -> stable anchor GO IDs (a cluster is that
# category if it contains ANY anchor). Robust to representative-term tie-breaks.
stress_cats <- list(
  "Protein folding / chaperone"                     = c("GO:0051082","GO:0006457","GO:0005832"),
  "Ubiquitin-proteasome / translational repression" = c("GO:0030371","GO:0043161","GO:0000209","GO:0061630"),
  "Ubiquitin ligase / proteolysis"                  = c("GO:0031625","GO:0019941","GO:0030162"),
  "Acid-base / V-type ATPase"                        = c("GO:0046961","GO:1902600","GO:0016471","GO:0033179"),
  "Na/K ionoregulation"                              = c("GO:0006883","GO:0036376","GO:0030007","GO:1990573"),
  "Mitochondrial calcium"                            = c("GO:0036444","GO:1990246","GO:0051560"),
  "Cilium / ciliary motility"                        = c("GO:0001669","GO:0036064","GO:0060294","GO:0031514","GO:0003341"),
  "Stress granule / mRNA processing"                 = c("GO:0006376","GO:0010494","GO:0008143","GO:0003730"),
  "Translation initiation"                           = c("GO:0003729","GO:0003743","GO:0001732"),
  "Steroid / progesterone metabolism"               = c("GO:0004508","GO:0047442","GO:0042448"))
cat_of <- function(gos) { for (nm in names(stress_cats)) if (any(gos %in% stress_cats[[nm]])) return(nm); NA_character_ }
cl_sum <- clG[, .(category = cat_of(GO), rep_term = Term[which.min(Fisher)], minF = min(Fisher),
                  CA = sum(CA), SE = sum(SE), CASE = sum(CASE), n = .N), by = cluster_id]
cl_sum[, cluster := data.table::fifelse(is.na(category), paste0("Other: ", rep_term), category)]
cl_sum[, stress  := !is.na(category)]
clG <- merge(clG, cl_sum[, .(cluster_id, cluster, stress, minF)], by = "cluster_id")

go_clusters_full   <- clG[order(minF, Fisher), .(cluster, stress, GO, Term, Fisher, CA, SE, CASE)]
go_clusters_stress <- cl_sum[stress == TRUE][order(minF), .(cluster, CA, SE, CASE, n, minF)]

# all clusters x treatment and x tier (for Fig S12)
go_clusters_trt <- clG[, .(CA = sum(CA), SE = sum(SE), CASE = sum(CASE),
                           minF = minF[1L], stress = stress[1L]), by = cluster]
tier_core <- unique(GO_Core_ALL$GeneOntologyIDs)
tier_conv <- unique(GO_Conv_ALL$GeneOntologyIDs)
tier_priv <- unique(GO_Priv_ALL$GeneOntologyIDs)
go_clusters_tier <- clG[, .(Core = sum(GO %in% tier_core), Convergent = sum(GO %in% tier_conv),
                            Private = sum(GO %in% tier_priv), minF = minF[1L], stress = stress[1L]), by = cluster]

fwrite(go_clusters_full,   file.path(tab_dir, "TableS_GO_clusters_full.csv"))
fwrite(go_clusters_stress, file.path(tab_dir, "TableS_GO_clusters_stress.csv"))
cat("data-driven clusters:", uniqueN(clG$cluster_id),
    "| stress-highlighted:", nrow(go_clusters_stress), "\n")
```

    ## data-driven clusters: 18 | stress-highlighted: 10

## Outlier allele-frequency matrices (for the Fig 2 PCAs)

``` bash
source activate CASE
# all-sample and per-block sync files restricted to the Total.Significant.Loci set
mawk '!/CHR/' Total.Significant.Loci | cut -f2,3 | sort | uniq > TSL.pca.filter
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' TSL.pca.filter CASE.All.Blocks.cov.sync | cut --complement -f60- > outlier.pca.sync
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' TSL.pca.filter CASE.All.Blocks.sync     | cut --complement -f60- > outlier.pca.block.sync
```

``` r
treatment.labels <- c(rep("CA",11), rep("CASE",11), rep("CON",11), rep("IS",12), rep("SE",11))
spawn.names <- c("B12","B11","B12","B10","B11","B10","B11","B12","B10","B12","B11","B10","B11","B12","B12","B10","B12","B11","B10","B11","B11","B12","B10","B11","B12","B11","B10","B12","B10","B11","B12","B11","B12","B10","B11","B10","B11","B10","B11","B10","B11","B12","B12","B12","B12","B11","B12","B10","B10","B11","B11","B10","B12","B12","B11","B12")
# all 56 samples over outlier loci
input.all.outlier <- read.sync(file = "outlier.pca.sync", gen = seq(1,56), repl = seq(1,56), polarization = "rising")
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
input.all.outlier.af <- af(input.all.outlier, gen = seq(1,56), repl = seq(1,56))
colnames(input.all.outlier.af) <- paste(treatment.labels, spawn.names, sep = "_")
# per-block (raw, uncovered) sync over the same loci
input.all.outlier.block <- read.sync(file = "outlier.pca.block.sync", gen = seq(1,56), repl = seq(1,56), polarization = "rising")
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
input.all.outlier.block.af <- af(input.all.outlier.block, gen = seq(1,56), repl = seq(1,56))
colnames(input.all.outlier.block.af) <- paste(treatment.labels, spawn.names, sep = "_")
```

# Fig 2 — Outlier loci & cross-spawn repeatability

## Figure parameters and guard

`alpha`/`alpha2` are re-asserted here (the tier-definition section above
sets them internally), then the in-memory objects from the pipeline are
validated.

``` r
gene_bed <- path.expand("~/CASE/analysis/sorted.ref3.0.gene.bed")
flank    <- 1000
alpha    <- 0.05
alpha2   <- 0.10
P        <- 10000
set.seed(1)
qcol   <- c(CA = "QCA", SE = "QSE", CASE = "QCASE", CON = "QCON")
blocks <- c("B10","B11","B12")
asL    <- function(x) as.logical(x)   # tolerate logical OR "TRUE"/"FALSE" character cols
```

``` r
need <- c(as.vector(outer(c("ca","case","se","con"), c("b10","b11","b12"),
                          function(t,b) paste0(t,".",b,".af"))),
          "core.sig","convergent.sig","private.sig","grouped_total.sig",
          "B10.pv","B11.pv","B12.pv")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Knit the pipeline first (same session). Missing: ",
                        paste(miss, collapse = ", "))

pv_list <- list(B10 = as.data.table(B10.pv),
                B11 = as.data.table(B11.pv),
                B12 = as.data.table(B12.pv))
for (b in names(pv_list)) {
  d <- pv_list[[b]]
  d[, BP := as.numeric(BP)]
  if ("CHR" %in% names(d))  d[, CHR := suppressWarnings(as.numeric(CHR))]
  if (!"SNP" %in% names(d)) d[, SNP := paste0(CHROM, "_", BP)]
  pv_list[[b]] <- d[!is.na(BP)]
}
```

``` r
# SNP -> gene (LOC) map; used by the gene Venn and the within-gene turnover panel.
raw <- fread(gene_bed, header = FALSE, sep = "\t", quote = "")
attr_col <- do.call(paste, c(as.list(raw), sep = ";"))
g <- sub('.*gene=(LOC[0-9]+).*', '\\1', attr_col)
genes <- data.table(CHROM = as.character(raw$V1), start = as.numeric(raw$V2),
                    end = as.numeric(raw$V3),
                    gene = ifelse(grepl('^LOC[0-9]+$', g), g, NA_character_))
genes <- genes[!is.na(gene) & !is.na(start) & !is.na(end)]
fb <- as.numeric(flank); genes[, start := pmax(start - fb, 0)][, end := end + fb]
setkey(genes, CHROM, start, end)
all_snps <- unique(rbindlist(lapply(pv_list, function(d) d[, .(SNP, CHROM, BP)])))
all_snps[, `:=`(start = BP, end = BP)]; setkey(all_snps, CHROM, start, end)
ov <- foverlaps(all_snps, genes, type = "within", nomatch = NULL)
s1 <- unique(ov[, .(SNP, gene)])[!duplicated(SNP)]
gene_of <- setNames(s1$gene, s1$SNP)
cat("loci mapped to a gene:", length(gene_of), "of", nrow(all_snps), "\n")
```

    ## loci mapped to a gene: 690961 of 813224

``` r
# LOC -> display symbol from the bed Name= attribute where present, else the LOC id.
nm <- sub('.*Name=([^;]+).*', '\\1', attr_col)
genes_raw_gene <- sub('.*gene=(LOC[0-9]+).*', '\\1', attr_col)
sym_dt <- data.table(gene = genes_raw_gene, sym = ifelse(grepl('Name=', attr_col), nm, NA_character_))
sym_dt <- unique(sym_dt[grepl('^LOC[0-9]+$', gene)])[!duplicated(gene)]
gene_sym <- setNames(ifelse(is.na(sym_dt$sym), sym_dt$gene, sym_dt$sym), sym_dt$gene)
disp_gene <- function(loc) { out <- gene_sym[loc]; out[is.na(out)] <- loc[is.na(out)]; unname(out) }
```

## 2.1 PCAs — ΔAF-from-initial (all blocks) + three per-block outlier PCAs

``` r
# ΔAF-from-initial PCA (ported from CASE_FULL_FINAL.Rmd L3655). Treatment fill, spawn shape.
M <- input.all.outlier.af
  is_b10 <- grep("^IS_B10", colnames(M)); is_b11 <- grep("^IS_B11", colnames(M))
  is_b12 <- grep("^IS_B12", colnames(M))
  if (!all(c("ISM_B10","ISM_B11","ISM_B12") %in% colnames(M)))
    M <- cbind(M, ISM_B10 = rowMeans(M[, is_b10], na.rm = TRUE),
                  ISM_B11 = rowMeans(M[, is_b11], na.rm = TRUE),
                  ISM_B12 = rowMeans(M[, is_b12], na.rm = TRUE))
  pick <- function(suf) { i <- grep(paste0(suf, "$"), colnames(M)); i[!grepl("^IS", colnames(M)[i])] }
  b10c <- pick("B10"); b11c <- pick("B11"); b12c <- pick("B12")
  dv <- list(); dn <- character()
  for (cc in b10c) { dv[[length(dv)+1]] <- M[,"ISM_B10"] - M[,cc]; dn <- c(dn, paste0(colnames(M)[cc],"_delta")) }
  for (cc in b11c) { dv[[length(dv)+1]] <- M[,"ISM_B11"] - M[,cc]; dn <- c(dn, paste0(colnames(M)[cc],"_delta")) }
  for (cc in b12c) { dv[[length(dv)+1]] <- M[,"ISM_B12"] - M[,cc]; dn <- c(dn, paste0(colnames(M)[cc],"_delta")) }
  newM <- do.call(cbind, dv); colnames(newM) <- dn
  pr <- prcomp(t(newM), center = TRUE, scale. = FALSE)
  ve <- pr$sdev^2 / sum(pr$sdev^2) * 100
  d  <- as.data.frame(pr$x[, 1:2]); names(d) <- c("PC1","PC2")
  d$Treatment <- factor(c(rep("CA",3),rep("CASE",3),rep("CON",3),rep("SE",3),
                          rep("CA",4),rep("CASE",4),rep("CON",4),rep("SE",4),
                          rep("CA",4),rep("CASE",4),rep("CON",4),rep("SE",4)),
                        levels = treat_levels)
  d$Spawn <- factor(c(rep("B10",12), rep("B11",16), rep("B12",16)), levels = spawn_levels)
  pca_dAF <- ggplot(d, aes(PC1, PC2, fill = Treatment, shape = Spawn)) +
    geom_point(size = 6, colour = "black", alpha = 0.6) +
    scale_fill_manual(values = treat_cols, name = "Treatment",
                      limits = c("CON","CA","SE","CASE","IS"), drop = FALSE) +
    scale_shape_manual(values = spawn_shapes, name = "Spawn") +
    guides(fill = "none",  # Panel A: Spawn legend only; colour keys to spawn so it merges with Panel G's
           shape = guide_legend(override.aes = list(fill = unname(spawn_cols)))) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
         title = "Outlier loci — ΔAF from initial, all blocks") +
    theme_case(13)
pca_dAF
```

![](Final_reproducible_analysis_files/figure-gfm/pca-dAF-1.png)<!-- -->

``` r
# per-block outlier PCAs, treatment-coloured (IS shown in sky blue)
run_pca <- function(af_mat, title = "", point_size = 5, shape = 21) {
  mt <- t(as.matrix(af_mat))
  ok <- apply(mt, 2, function(col) all(is.finite(col)) && stats::var(col) > 0)
  mt <- mt[, ok, drop = FALSE]
  pr <- prcomp(mt, center = TRUE, scale. = FALSE)   # covariance PCA: let strongly-shifting outlier loci drive treatment separation
  ve <- pr$sdev^2 / sum(pr$sdev^2) * 100
  d  <- as.data.frame(pr$x[, 1:2]); names(d) <- c("PC1","PC2")
  nm <- rownames(mt)
  d$Treatment <- factor(sub("_.*", "", nm), levels = c(treat_levels, "IS"))
  ggplot(d, aes(PC1, PC2, fill = Treatment)) +
    geom_point(size = point_size, shape = shape, colour = "black", alpha = 0.8) +
    scale_fill_manual(values = treat_cols, name = "Treatment",
                      limits = c("CON","CA","SE","CASE","IS"), drop = FALSE) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), title = title) +
    theme_case(13)
}
block_pcas <- list()
pos_files <- c(B10 = "B10.pos", B11 = "B11.pos", B12 = "B12.pos")
src <- input.all.outlier.block.af; colnames(src) <- trimws(colnames(src))
for (b in blocks) {
  idx <- grep(paste0("_", b, "$"), colnames(src))
  # Restrict to loci that passed THIS block's coverage/MAF filter (its .pos file);
  # locus IDs are CHROM.POS, as in the af() rownames (matches CASE_FULL_FINAL.Rmd).
  pos  <- read.table(pos_files[[b]], header = FALSE, stringsAsFactors = FALSE)
  rows <- rownames(src) %in% paste(pos[, 1], pos[, 2], sep = ".")
  cat(sprintf("  %s: %d of %d outlier loci retained by .pos filter\n", b, sum(rows), nrow(src)))
  p <- run_pca(src[rows, idx, drop = FALSE], paste0("Outlier loci — ", b), shape = spawn_shapes[[b]])
  # guides(fill = "none") (layer-level) survives the figure-wide `& theme(...)` in fig2-assemble.
  if (b != "B10") p <- p + guides(fill = "none")   # only B10 keeps the Treatment legend
  block_pcas[[b]] <- p
}
```

    ##   B10: 3093 of 3472 outlier loci retained by .pos filter
    ##   B11: 2044 of 3472 outlier loci retained by .pos filter
    ##   B12: 2979 of 3472 outlier loci retained by .pos filter

## 2.2 Reproducibility bars (3b) + lethality coupling (3c)

``` r
# True-outlier universe per treatment: SNPs flagged Sig.<trt> TRUE in any tier
# (Core/Convergent/Private) — same definition feeding the Venn/Sig.loci
# files and the turnover panel's true_case set (cf. lines 372-375).
true_outliers <- function(trt) {
  sigcol <- paste0("Sig.", trt)
  unique(rbindlist(list(
    data.table(SNP = core.sig$SNP,       s = asL(core.sig[[sigcol]])),
    data.table(SNP = convergent.sig$SNP,  s = asL(convergent.sig[[sigcol]])),
    data.table(SNP = private.sig$SNP, s = asL(private.sig[[sigcol]]))
  ))[s == TRUE, SNP])
}
# bed-based CON-confound exclusion (same sets used for the "CON Filtered" tier
# in the Manhattan panels, lines 493-498)
con_excl <- unique(c(all_con$SNP, con_filter$SNP, singleton$SNP))
block_out_snps <- function(b, trt) {
  pv <- pv_list[[b]]
  if (trt == "CON") return(pv[get(qcol["CON"]) < alpha, SNP])
  to <- true_outliers(trt)
  s  <- pv[SNP %in% to & get(qcol[trt]) < alpha & QCON > alpha2, SNP]
  setdiff(s, con_excl)
}
shared_set <- function(sets) {
  a <- unique(sets[[1]]); b <- unique(sets[[2]]); c <- unique(sets[[3]])
  unique(c(intersect(a,b), intersect(a,c), intersect(b,c)))
}
shared_counts <- function(sets) length(shared_set(sets))

drop_na  <- function(x) x[!is.na(x)]
gmap     <- function(s) unname(gene_of[s])

# (a) SNP-level shareable universe: shared_counts() can only register a SNP
# as "reproduced" if it's present in >=2 of the 3 blocks' pv tables, so the
# random/CON null universes should be restricted to that shareable set —
# otherwise block-private SNPs dilute con_mean/random_mean and inflate
# enrichment factors.
all_snp        <- lapply(blocks, function(b) pv_list[[b]]$SNP)
shareable_snps <- shared_set(all_snp)

univ_snp <- lapply(blocks, function(b) intersect(pv_list[[b]]$SNP, shareable_snps))
con_snp  <- lapply(blocks, function(b) intersect(pv_list[[b]][QCON < alpha, SNP], shareable_snps))

# (b) Gene-level shareability: a gene can only register in the gene-level
# shared_counts() if it's reachable (via any SNP) from >=2 of the 3 blocks,
# independent of which specific SNP maps to it. Mark non-shareable genes NA
# so drop_na() excludes them from the gene-level draws.
all_gene        <- lapply(all_snp, function(s) unique(drop_na(gmap(s))))
shareable_genes <- shared_set(all_gene)

mask_unshareable_genes <- function(s) {
  g <- gmap(s)
  g[!(g %in% shareable_genes)] <- NA
  g
}
univ_gene <- lapply(univ_snp, mask_unshareable_genes)
con_gene  <- lapply(con_snp,  mask_unshareable_genes)

# Sanity check: gene-level fold-enrichment is much smaller than locus-level
# because gene_of collapses many SNPs onto few genes, so random SNP sets
# overlap at the gene level far more often by chance.
cat("Shareable SNP universe size:", length(unique(unlist(univ_snp))),
    " | Shareable gene universe size:", length(unique(unlist(drop_na(unlist(univ_gene))))), "\n")
```

    ## Shareable SNP universe size: 338070  | Shareable gene universe size: 8997

``` r
tx_test <- function(trt) {
  obs_snp <- lapply(blocks, block_out_snps, trt = trt)
  n_b <- vapply(obs_snp, length, integer(1))
  o_loc <- shared_counts(obs_snp)
  o_gen <- shared_counts(lapply(obs_snp, function(s) drop_na(gmap(s))))
  rl <- cl <- rg <- cg <- numeric(P)
  for (p in seq_len(P)) {
    ri <- lapply(seq_along(blocks), function(i) sample.int(length(univ_snp[[i]]), min(n_b[i], length(univ_snp[[i]]))))
    ci <- lapply(seq_along(blocks), function(i){ np <- length(con_snp[[i]]); sample.int(np, min(n_b[i], np), replace = n_b[i] > np) })
    rl[p] <- shared_counts(lapply(seq_along(blocks), function(i) univ_snp[[i]][ri[[i]]]))
    cl[p] <- shared_counts(lapply(seq_along(blocks), function(i) con_snp[[i]][ci[[i]]]))
    rg[p] <- shared_counts(lapply(seq_along(blocks), function(i) drop_na(univ_gene[[i]][ri[[i]]])))
    cg[p] <- shared_counts(lapply(seq_along(blocks), function(i) drop_na(con_gene[[i]][ci[[i]]])))
  }
  mk <- function(lvl, obs, rnd, conb) data.table(treatment = trt, level = lvl, obs_2of3 = obs,
    random_mean = mean(rnd), enrich_vs_random = obs/mean(rnd), p_random = (1 + sum(rnd >= obs))/(P + 1),
    con_mean = mean(conb), enrich_vs_CON = obs/mean(conb), p_CON = (1 + sum(conb >= obs))/(P + 1))
  out <- mk("locus", o_loc, rl, cl)
  out <- rbind(out, mk("gene", o_gen, rg, cg))
  out
}
core_res <- rbindlist(lapply(c("CA","SE","CASE"), tx_test))
core_res[, treatment := factor(treatment, levels = c("CA","SE","CASE"))]
fwrite(core_res, file.path(tab_dir, "Fig2_reproducibility_summary.csv"))
```

``` r
plt <- melt(core_res[, .(treatment, level, Observed = obs_2of3,
                         `CON baseline` = con_mean, `Random null` = random_mean)],
            id.vars = c("treatment","level"), variable.name = "src", value.name = "n")
plt[, src := factor(src, levels = c("Observed","CON baseline","Random null"))]
plt[, level := factor(level, levels = c("locus","gene"))]
# anchor every annotation to a common line above each facet's max (so the tall CASE/gene
# label can't ride into the panel top)
ymax_by <- plt[, .(ymax = max(n)), by = level]
lab <- merge(core_res[, .(treatment, level = factor(level, levels = c("locus","gene")),
                          enrich_vs_CON, p_CON)], ymax_by, by = "level")
# stagger annotation height by treatment so labels don't collide: CA lower, CASE higher
ymax_mult <- c(CA = 1.05, SE = 1.16, CASE = 1.30)
lab[, `:=`(y = ymax * ymax_mult[as.character(treatment)],
           txt = sprintf("%.1fx vs CON\n(p=%.4f)", enrich_vs_CON, p_CON))]
# Lollipop instead of grouped bars: a stem to each value + a point, dodged by series.
dodge <- position_dodge(width = 0.6)
fig3b <- ggplot(plt, aes(treatment, n, colour = src, group = src)) +
  geom_linerange(aes(ymin = 0, ymax = n), position = dodge, linewidth = 0.7) +
  geom_point(position = dodge, size = 3.4) +
  geom_text(data = lab, inherit.aes = FALSE, aes(treatment, y, label = txt),
            vjust = 0, size = 2.2, lineheight = 0.9) +
  facet_wrap(~ level, scales = "free_y") +
  scale_colour_manual(values = stat_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +    # headroom for the common annotation line
  labs(x = NULL, y = "Shared by >=2 of 3 spawns",
       title = "Cross-block reproducibility vs matched baselines") +
  theme_case(13) + theme(strip.background = element_blank(),
                         strip.text = element_text(face = "bold"),
                         legend.position = "right",                 # vertical (inside-panel is unreliable for facets)
                         legend.key.size = unit(0.35, "cm"), plot.title.position = "plot")

seq_blocks <- c(10, 11, 12)
# composite-denominator survival per block x treatment (final_data_end$Survival = composite)
mc <- final_data_end[Block %in% seq_blocks, .(surv = mean(Survival)), by = .(Block, Treatment)]
con_surv <- mc[Treatment == "CON", .(Block, con = surv)]
mc <- merge(mc[Treatment != "CON"], con_surv, by = "Block")
mc[, excess_mort := con - surv]
leth <- mc[, .(excess_mort = mean(excess_mort)), by = Treatment]
coup <- merge(leth, core_res[level == "locus", .(treatment, enrich_vs_CON)],
              by.x = "Treatment", by.y = "treatment")
coup[, Treatment := factor(Treatment, levels = c("CA","SE","CASE"))]
setorder(coup, excess_mort)
fig3c <- ggplot(coup, aes(excess_mort, enrich_vs_CON)) +
  geom_line(linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_point(aes(fill = Treatment), shape = 21, size = 8, colour = "black", alpha = 0.9) +
  geom_text(aes(label = Treatment), vjust = -1.6, size = 4.5) +
  scale_fill_manual(values = treat_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = "Excess larval mortality vs CON\n(mean of B10–B12, % pts)",
       y = "Reproducibility (fold vs CON)",
       title = "Complexity Increases Response, Lowers Reproducibility") +
  theme_case(13) + theme(legend.position = "none", plot.title.position = "plot")
```

## 2.3 Within-gene allele turnover (3E-A) — same genes, different SNPs

``` r
  # TRUE CASE outliers only = SNPs in the tier sets (Core/Convergent/Private) with
  # Sig.CASE, that are genome-wide significant in a given block. Far fewer points than the
  # raw per-block FDR set, so the dots no longer pile up.
  true_case <- unique(rbindlist(list(
    data.table(SNP = core.sig$SNP,       s = asL(core.sig$Sig.CASE)),
    data.table(SNP = convergent.sig$SNP,  s = asL(convergent.sig$Sig.CASE)),
    data.table(SNP = private.sig$SNP, s = asL(private.sig$Sig.CASE))))[s == TRUE, .(SNP)])
  out_long <- rbindlist(lapply(blocks, function(b) {
    pv <- pv_list[[b]]
    s  <- pv[SNP %in% true_case$SNP & QCASE < alpha & QCON > alpha2, SNP]
    data.table(SNP = s, gene = unname(gene_of[s]), block = b)[!is.na(gene)]
  }))
  # within-gene position (kb), then FOCUS on SNPs within the first 25 kb of each gene
  bp_lk <- unique(rbindlist(lapply(pv_list, function(d) d[, .(SNP, BP)])))
  out_long <- merge(out_long, bp_lk, by = "SNP", all.x = TRUE)
  out_long[, rel := (BP - min(BP)) / 1000, by = gene]
  near <- out_long[rel <= 10]
  # recurrent (>=2 spawns), readable SNP count, preferring all-three-spawn genes
  gstat <- near[, .(nblk = uniqueN(block), nsnp = uniqueN(SNP)), by = gene]
  cand  <- gstat[nblk >= 2 & nsnp %between% c(3, 12)]
  setorder(cand, -nblk, -nsnp)
  top_genes <- head(cand$gene, 8)
  if (length(top_genes) < 8) {
    extra <- setdiff(gstat[nblk >= 2][order(-nblk, -nsnp), gene], top_genes)
    top_genes <- head(c(top_genes, extra), 8)
  }
  ex <- near[gene %in% top_genes]
  disp_top <- disp_gene(top_genes)
  ex[, glab  := factor(disp_gene(gene), levels = rev(unique(disp_top)))]
  ex[, block := factor(block, levels = spawn_levels)]
  ex[, gnum  := as.integer(glab)]
  # three horizontal spawn tracks per gene (B10 top, B11 mid, B12 bottom) to declutter
  off <- c(B10 = 0.26, B11 = 0, B12 = -0.26)
  ex[, y2 := gnum + off[as.character(block)]]
  tracks <- unique(ex[, .(gnum, block, y2)])
  turnover <- ggplot() +
    geom_segment(data = tracks, aes(0, y2, xend = 10, yend = y2, colour = block),
                 linewidth = 0.3, alpha = 0.35) +
    geom_point(data = ex, aes(rel, y2, fill = block, shape = block),
               size = 2.4, colour = "black", alpha = 0.75) +
    scale_fill_manual(values = spawn_cols, guide = "none") +
    scale_colour_manual(values = spawn_cols, guide = "none") +   # track-line colour, no legend
    scale_shape_manual(values = spawn_shapes, name = "Spawn") +
    # Panel A already carries the (identical) Spawn legend in the collected guides;
    # suppress G's copy via guides() so it doesn't duplicate in the guide_area.
    guides(shape = "none") +
    scale_y_continuous(breaks = sort(unique(ex$gnum)),
                       labels = levels(ex$glab)[sort(unique(ex$gnum))]) +
    coord_cartesian(xlim = c(0, 6)) +
    labs(x = "Position within gene region (kb)", y = NULL,
         title = "Same genes, different loci across blocks", plot.title.position = "plot") +
    theme_case(13)
turnover
```

![](Final_reproducible_analysis_files/figure-gfm/fig2-turnover-1.png)<!-- -->
\## 2.4 Assemble Fig 2 (patchwork)

``` r
design2 <- "
AABBCC
DDEEFF
GGGGGH
"
fig2 <-  block_pcas$B10 + block_pcas$B11 + block_pcas$B12 + pca_dAF +
          fig3b + fig3c + free(turnover) + guide_area() +
          plot_layout(design = design2, guides = "collect") +   # collect legends into the H guide area
          plot_annotation(tag_levels = list(c("A","B","C","D","E","F","G",""))) &
          # legend.position = "right" -> vertical direction by default, suits the narrow
          # guide_area column (H) better than "bottom". Per-plot guides(fill/shape = "none")
          # (set on B11/B12 PCAs and fig3c) survive this `&` theme since guides() is
          # layer-level, not a theme element, so it can't be overridden here.
          theme(legend.position = "right",
                legend.direction = "vertical",
                legend.box = "vertical",
                plot.title = element_text(size = 12),
                legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 10),
                legend.title = element_text(face = "bold", size = 11), plot.margin = margin(t =1, r = 20, b = 5, l = 1, unit = "pt"))
save_case(fig2, "Fig2_outliers_repeatability.png", width = 15, height = 12)
fig2
```

![](Final_reproducible_analysis_files/figure-gfm/fig2-assemble-1.png)<!-- -->

## Fig S3 — Genome-wide PCA (10,000 random loci)

``` r
afr_t <- t(input.all.rand.af)
pr_r  <- prcomp(afr_t, center = TRUE, scale. = TRUE)
ve_r  <- pr_r$sdev^2 / sum(pr_r$sdev^2) * 100
dr    <- as.data.frame(pr_r$x[, 1:2]); names(dr) <- c("PC1","PC2")
dr$Treatment <- factor(sub("_.*", "", rownames(afr_t)), levels = c(treat_levels, "IS"))
dr$Spawn     <- factor(sub(".*_", "", rownames(afr_t)), levels = spawn_levels)
s2_lab <- labs(x = sprintf("PC1 (%.1f%%)", ve_r[1]), y = sprintf("PC2 (%.1f%%)", ve_r[2]))
p_s2_treat <- ggplot(dr, aes(PC1, PC2, fill = Treatment)) +
  geom_point(size = 4, shape = 21, colour = "black", alpha = 0.8) +
  scale_fill_manual(values = treat_cols, drop = FALSE) +
  s2_lab + labs(title = "By treatment") + theme_case(13)
p_s2_spawn <- ggplot(dr, aes(PC1, PC2, fill = Spawn)) +
  geom_point(size = 4, shape = 21, colour = "black", alpha = 0.8) +
  scale_fill_manual(values = spawn_cols) +
  s2_lab + labs(title = "By spawn") + theme_case(13)
# Venn of SNPs called per spawn (built above into VenSNPs.png); embed as a raster panel so the
# VennDiagram grob drops cleanly into the patchwork, the same approach used for the Fig 1 schematic.
venn_panel <- patchwork::wrap_elements(full = grid::rasterGrob(png::readPNG("VenSNPs.png"),
                                                               interpolate = TRUE))
figS2 <- venn_panel + p_s2_treat + p_s2_spawn +
  plot_annotation(tag_levels = "A",
                  title = "Fig S3 — Called-SNP overlap and genome-wide PCA (10,000 random loci)")
save_supp(figS2, "FigS3_genomewide_PCA.png", width = 18, height = 5.5)
figS2
```

![](Final_reproducible_analysis_files/figure-gfm/figS2-genomewide-pca-1.png)<!-- -->

## Fig S4 — Tier-specific PCA

## Tier-specific PCA

### Extract tier-specific sync files

``` bash
source activate CASE

# Build per-tier position filters from the Significant.Loci files created in the
# CON-filtering section. Each filter is a CHR/BP table used to subset the
# full sync file.
# NOTE: eval=TRUE (lightweight mawk subsetting of existing sync files). Must run
# so the renamed Core/Conv/Priv .outlier.pca.sync files exist for the tier PCAs.

# --- Core loci ---
mawk '!/CHR/' Core.Significant.Loci | cut -f2,3 | sort | uniq > Core.pca.filter
# All-blocks sync (cov-filtered, columns 1–56 only)
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' Core.pca.filter CASE.All.Blocks.cov.sync \
  | cut --complement -f60- > Core.outlier.pca.sync
# Raw sync (for per-block PCA loop below)
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' Core.pca.filter CASE.All.Blocks.sync \
  | cut --complement -f60- > Core.outlier.pca.block.sync

# --- Private loci ---
mawk '!/CHR/' Priv.Significant.Loci | cut -f2,3 | sort | uniq > Priv.pca.filter
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' Priv.pca.filter CASE.All.Blocks.cov.sync \
  | cut --complement -f60- > Priv.outlier.pca.sync
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' Priv.pca.filter CASE.All.Blocks.sync \
  | cut --complement -f60- > Priv.outlier.pca.block.sync

# --- Convergent loci ---
mawk '!/CHR/' Conv.Significant.Loci | cut -f2,3 | sort | uniq > Conv.pca.filter
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' Conv.pca.filter CASE.All.Blocks.cov.sync \
  | cut --complement -f60- > Conv.outlier.pca.sync
mawk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' Conv.pca.filter CASE.All.Blocks.sync \
  | cut --complement -f60- > Conv.outlier.pca.block.sync
```

### All-samples PCA helper function

``` r
# run_tier_pca(): reads a sync file for one locus tier, computes allele
# frequencies, runs PCA (all 56 samples), and saves two PNGs — one colored
# by Treatment, one by Spawn (block).
#
# Arguments:
#   sync_file   — path to the tier-filtered .sync file (cov-filtered, 56 cols)
#   tier_label  — human-readable label used in plot titles
#   out_prefix  — file-name prefix for PNG output (e.g. "PC_Core")
#
# Returns invisibly: list(pca_df, var_explained, p_treat, p_spawn)
run_tier_pca <- function(sync_file, tier_label, out_prefix) {
  cat(paste0("\n=== PCA: ", tier_label, " ===\n"))

  # Sample metadata — same order as in all existing PCA blocks
  treatment.labels <- c(rep("CA", 11), rep("CASE", 11), rep("CON", 11),
                        rep("IS", 12), rep("SE", 11))
  spawn.names <- c("B12","B11","B12","B10","B11","B10","B11","B12","B10","B12",
                   "B11","B10","B11","B12","B12","B10","B12","B11","B10","B11",
                   "B11","B12","B10","B11","B12","B11","B10","B12","B10","B11",
                   "B12","B11","B12","B10","B11","B10","B11","B10","B11","B10",
                   "B11","B12","B12","B12","B12","B11","B12","B10","B10","B11",
                   "B11","B10","B12","B12","B11","B12")

  inp    <- read.sync(file = sync_file, gen = seq(1, 56), repl = seq(1, 56),
                      polarization = "rising")
  inp.af <- af(inp, gen = seq(1, 56), repl = seq(1, 56))
  colnames(inp.af) <- paste(treatment.labels, spawn.names, sep = "_")

  af_t    <- t(inp.af)
  pca_res <- prcomp(af_t, center = TRUE, scale. = TRUE)
  pca_df  <- as.data.frame(pca_res$x)
  var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

  pca_df$Population <- spawn.names
  pca_df$Treatments <- treatment.labels

  cat(sprintf("  Loci: %d   PC1: %.1f%%   PC2: %.1f%%\n",
              nrow(inp.af), var_exp[1], var_exp[2]))

  base_theme <- list(
    theme_classic(),
    guides(alpha = "none"),
    theme(axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24)),
    labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
         title = paste0("PCA — ", tier_label))
  )

  p_treat <- ggplot(pca_df, aes(x = PC1, y = PC2,
                                fill = Treatments, shape = as.factor(Population))) +
    geom_point(size = 10, color = "black", alpha = 0.5) +
    scale_fill_manual(values = cbPaletteSmall, name = "Treatment") +
    scale_shape_manual(values = c(21, 22, 23), name = "Spawn") +
    guides(fill  = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = cbPaletteSmall[4]))) +
    base_theme

  p_spawn <- ggplot(pca_df, aes(x = PC1, y = PC2,
                                fill = as.factor(Population), shape = as.factor(Population))) +
    geom_point(size = 10, color = "black", alpha = 0.5) +
    scale_fill_manual(values = cbPaletteSmall, name = "Spawn") +
    scale_shape_manual(values = c(21, 22, 23), name = "Spawn") +
    guides(fill  = guide_legend(override.aes = list(shape = 21, alpha = 1)),
           shape = "none") +
    base_theme

  png(paste0(out_prefix, "_treat.png"), type = "cairo", units = "px",
      width = 5400, height = 3000, res = 300, bg = "transparent")
  print(p_treat); dev.off()

  png(paste0(out_prefix, "_spawn.png"), type = "cairo", units = "px",
      width = 5400, height = 3000, res = 300, bg = "transparent")
  print(p_spawn); dev.off()

  invisible(list(pca_df = pca_df, var_explained = var_exp,
                 p_treat = p_treat, p_spawn = p_spawn))
}
```

### Run all-samples PCA for each tier

``` r
# All three tiers run across all 56 samples; block-restrict Convergent/Private later
# if the all-block PCA is uninformative.
pca.core <- run_tier_pca("Core.outlier.pca.sync",  "Core Loci",           "PC_Core")
```

    ## 
    ## === PCA: Core Loci ===
    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...
    ##   Loci: 2111   PC1: 27.8%   PC2: 15.0%

``` r
pca.conv  <- run_tier_pca("Conv.outlier.pca.sync",   "Convergent Loci",     "PC_Conv")
```

    ## 
    ## === PCA: Convergent Loci ===
    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...
    ##   Loci: 493   PC1: 31.2%   PC2: 15.6%

``` r
pca.priv   <- run_tier_pca("Priv.outlier.pca.sync",    "Private Loci",        "PC_Priv")
```

    ## 
    ## === PCA: Private Loci ===
    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...
    ##   Loci: 439   PC1: 46.7%   PC2: 14.3%

``` r
# Print in-line
print(pca.core$p_treat)
```

![](Final_reproducible_analysis_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

``` r
print(pca.conv$p_treat)
```

![](Final_reproducible_analysis_files/figure-gfm/unnamed-chunk-59-2.png)<!-- -->

``` r
print(pca.priv$p_treat)
```

![](Final_reproducible_analysis_files/figure-gfm/unnamed-chunk-59-3.png)<!-- -->

``` r
figS3 <- (pca.core$p_treat + pca.conv$p_treat + pca.priv$p_treat) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
                  title = "Fig S4 — Tier-specific PCA (treatment-coloured): Core | Convergent | Private")
save_supp(figS3, "FigS4_tier_PCA.png", width = 16, height = 5.5)
figS3
```

![](Final_reproducible_analysis_files/figure-gfm/figS3-tier-pca-1.png)<!-- -->

# Fig 3 — Magnitude synergy: control-anchored additivity

Magnitude synergy asks whether the combined-stressor allele-frequency
change exceeds the *sum* of the single-stressor changes. Both axes are
anchored on the control to remove the shared lab/handling response — the
correction the pool-seq multiple-stressor precedent (Brennan et al.
2022; Burny & Schlötterer 2022) treats as essential:
`y = ΔAF_CASE − ΔAF_CON` against
`x = (ΔAF_CA − ΔAF_CON) + (ΔAF_SE − ΔAF_CON)`, coverage floor 20. Fits
use one point per locus (averaged across the blocks it appears in); the
scatter shows the per-block points coloured by spawn on a background of
a 10,000-locus subsample of non-outlier loci (the genome-wide additive
expectation). Loci are split into significance-overlap panels — four
mutually exclusive (CASE only; CASE & SE; CASE & CA; CA & SE with no
CASE) plus two reference supersets (CASE significant overall; CASE & SE
& CA). The dashed line is the additive 1:1 null; the solid black line is
the per-panel orthogonal (Deming) fit. Because outliers are ascertained
on their own significance, points sit off the additive line partly by
construction — **sub-additivity (where CA contributes) is the robust
reading; a CASE-specific supra-additive tendency is the hedged one.**

``` r
add_floor  <- 20
add_blocks <- c("b10","b11","b12")
add_design <- list(b10 = list(gen = c(0,1,0,1,0,1),       repl = c(1,1,2,2,3,3)),
                   b11 = list(gen = c(0,1,0,1,0,1,0,1),   repl = c(1,1,2,2,3,3,4,4)),
                   b12 = list(gen = c(0,1,0,1,0,1,0,1),   repl = c(1,1,2,2,3,3,4,4)))
add_us <- function(x) sub("\\.(?=[^.]+$)", "_", x, perl = TRUE)

# treatment-level ΔAF per locus per block: mean(post) - mean(initial), coverage-filtered per pool.
add_delta <- function(tr, blk) {
  af  <- as.matrix(get(paste0(tolower(tr), ".", blk, ".af")))
  cov <- as.matrix(get(paste0(tolower(tr), ".", blk, ".cov")))
  des <- add_design[[blk]]; gen <- des$gen
  i0 <- which(gen == 0); i1 <- which(gen == 1)
  mm <- function(cols) { a <- af[, cols, drop = FALSE]; c <- cov[, cols, drop = FALSE]
    a[!(is.finite(a) & is.finite(c) & c >= add_floor)] <- NA_real_; rowMeans(a, na.rm = TRUE) }
  data.table(SNP = add_us(rownames(af)), d = mm(i1) - mm(i0))
}

# per-block anchored changes (for the coloured scatter)
add_perblk <- rbindlist(lapply(add_blocks, function(blk) {
  d <- Reduce(function(a, b) merge(a, b, by = "SNP"), list(
    add_delta("CA",   blk)[, .(SNP, dCA   = d)], add_delta("SE",   blk)[, .(SNP, dSE   = d)],
    add_delta("CASE", blk)[, .(SNP, dCASE = d)], add_delta("CON",  blk)[, .(SNP, dCON  = d)]))
  d <- d[is.finite(dCA) & is.finite(dSE) & is.finite(dCASE) & is.finite(dCON)]
  d[, .(SNP, block = toupper(blk), y = dCASE - dCON, x = (dCA - dCON) + (dSE - dCON))]
}))
add_perblk[, block := factor(block, levels = toupper(add_blocks))]

# one point per locus (averaged across blocks) — used for panel membership and the fits
add_pl <- add_perblk[, .(x = mean(x), y = mean(y), n_blk = .N), by = SNP]
add_gts  <- as.data.table(grouped_total.sig)
add_flag <- function(col) unique(add_gts[as.logical(get(col)), SNP])
sc_CASE <- add_flag("Sig.CASE"); sc_CA <- add_flag("Sig.CA"); sc_SE <- add_flag("Sig.SE")
add_pl[, `:=`(sCASE = SNP %in% sc_CASE, sCA = SNP %in% sc_CA, sSE = SNP %in% sc_SE)]

# Panels: four mutually exclusive + two reference supersets (a locus can appear in more than one
# of the supersets, so these are built as explicit membership lists rather than a single factor).
add_lev <- c("CASE significant (all)","CASE only", "CASE & SE", "CASE & CA", "CASE & SE & CA","CA & SE (no CASE)")
add_panels <- list(
  `CASE only`              = add_pl[ sCASE & !sCA & !sSE, SNP],
  `CASE & SE`               = add_pl[ sCASE &  sSE & !sCA, SNP],
  `CASE & CA`               = add_pl[ sCASE &  sCA & !sSE, SNP],
  `CA & SE (no CASE)`       = add_pl[ sCA   &  sSE & !sCASE, SNP],
  `CASE significant (all)`  = add_pl[ (sCASE), SNP],
  `CASE & SE & CA`           = add_pl[ sCASE & sCA & sSE, SNP])
cat("Fig 3 loci per panel:\n")
```

    ## Fig 3 loci per panel:

``` r
print(rbindlist(lapply(add_lev, function(nm) data.table(panel = nm, n = length(add_panels[[nm]])))))
```

    ##                     panel     n
    ##                    <char> <int>
    ## 1: CASE significant (all)  2272
    ## 2:              CASE only  1941
    ## 3:              CASE & SE   232
    ## 4:              CASE & CA    70
    ## 5:         CASE & SE & CA    29
    ## 6:      CA & SE (no CASE)    64

``` r
add_plpanel <- rbindlist(lapply(add_lev, function(nm) { s <- add_panels[[nm]]
  if (!length(s)) return(NULL); add_pl[SNP %in% s, .(SNP, x, y, panel = nm)] }))
add_plpanel[, panel := factor(panel, levels = add_lev)]
add_perblk_panel <- rbindlist(lapply(add_lev, function(nm) { s <- add_panels[[nm]]
  if (!length(s)) return(NULL); add_perblk[SNP %in% s, .(SNP, block, x, y, panel = nm)] }))
add_perblk_panel[, panel := factor(panel, levels = add_lev)]
add_bgpool <- add_pl[!sCASE & !sCA & !sSE]      # non-outlier background pool
```

## 3.1 Additivity scatter and per-panel orthogonal fits

``` r
# orthogonal (Deming, error-variance ratio 1) fit: does not treat x as error-free
deming <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  sxx <- var(x); syy <- var(y); sxy <- cov(x, y)
  b <- (syy - sxx + sqrt((syy - sxx)^2 + 4 * sxy^2)) / (2 * sxy)
  c(intercept = mean(y) - b * mean(x), slope = b)
}
# one row per locus per panel for the fits (no block-level pseudoreplication)
add_fitdt  <- add_plpanel[, as.list(deming(x, y)), by = panel]        # panel, intercept, slope
add_catsum <- add_plpanel[, { f <- sign(y); f[f == 0] <- 1
  .(n = .N, slope = unname(deming(x, y)["slope"]),
    mean_excess = mean(f * y - f * x), frac_supra = mean(f * y > f * x)) }, by = panel][order(panel)]
print(add_catsum)
```

    ##                     panel     n      slope  mean_excess frac_supra
    ##                    <fctr> <int>      <num>        <num>      <num>
    ## 1: CASE significant (all)  2272 0.98847451  0.018333016  0.6637324
    ## 2:              CASE only  1941 1.10061265  0.022741696  0.7042761
    ## 3:              CASE & SE   232 0.90402555 -0.002705815  0.4698276
    ## 4:              CASE & CA    70 0.54847296 -0.014540630  0.3571429
    ## 5:         CASE & SE & CA    29 0.55220385 -0.029083654  0.2413793
    ## 6:      CA & SE (no CASE)    64 0.08678093 -0.013755317  0.4531250

``` r
cat(sprintf("\nGenome-wide orthogonal slope (all loci, reference): %.2f\n",
            unname(deming(add_pl$x, add_pl$y)["slope"])))
```

    ## 
    ## Genome-wide orthogonal slope (all loci, reference): 0.49

``` r
fwrite(add_plpanel, file.path(tab_dir, "Fig3_additivity_perlocus.csv"))
fwrite(add_catsum,  file.path(tab_dir, "Fig3_additivity_by_category.csv"))

rng <- quantile(c(add_pl$x, add_pl$y), c(0.0001, 0.9999), na.rm = TRUE)
#set.seed(1); n_bg <- 100000
#bg_snps <- add_bgpool[sample.int(.N, min(n_bg, .N)), SNP]
bg_snps <- add_bgpool[n_blk == 3, SNP]
bg <- rbindlist(lapply(add_lev, function(L) { d <- add_pl[SNP %in% bg_snps, .(x, y)]; d[, panel := L]; d }))
bg[, panel := factor(panel, levels = add_lev)]
ann <- add_fitdt[add_catsum, on = "panel"][, lab := sprintf("slope %.2f\n%.0f%% above 1:1", slope, 100 * frac_supra)]

fig3_add <- ggplot(mapping = aes(x, y)) +
  geom_bin2d(data = bg, bins = 50) +
  scale_fill_gradient(name = sprintf("non-outlier loci", length(bg_snps)),
                      low = "grey92", high = "grey55", trans = "log10") +
  geom_point(data = add_perblk_panel, aes(colour = block), size = 0.75, alpha = 0.5) +
  scale_colour_manual(name = "block", values = spawn_cols) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
  geom_abline(data = add_fitdt, aes(slope = slope, intercept = intercept), colour = "black", linewidth = 0.6) +
  geom_text(data = ann, aes(label = lab), x = rng[1], y = rng[2], hjust = 0, vjust = 1,
            size = 3, lineheight = 0.9, inherit.aes = FALSE) +
  facet_wrap(~ panel, ncol = 3) +
  coord_equal(xlim = rng, ylim = rng) +
  #scale_x_continuous(expand=c(0,0), limits = rng) +
  #scale_y_continuous(expand=c(0,0), limits = rng) +
  labs(x = "(ΔAF CA − ΔAF CON) + (ΔAF SE − ΔAF CON)   [additive expectation]",
       y = "ΔAF CASE − ΔAF CON   [combined response]",
       title = "Fig 3 — Magnitude synergy: control-anchored additivity",
       subtitle = "dashed = additive 1:1;  solid black = per-panel orthogonal fit;  grey = non-outlier subsample;  colour = block") +
  theme_case(12)
fig3_add
```

![](Final_reproducible_analysis_files/figure-gfm/fig3-additivity-1.png)<!-- -->

``` r
save_case(fig3_add, "Fig3_magnitude_additivity.png", width = 11, height = 7)
```

## Fig S5 — Additivity scatter faceted by replication tier

Robustness of Fig 3 across the replication tiers (Core / Convergent /
Private), reporting OLS alongside the orthogonal fit so the steeper
orthogonal slope is not read as supra-additivity on its own (resolves
the slope ~1.48 caveat). One point per locus (averaged across blocks).

``` r
tier_objs <- list(Core = core.sig, Convergent = convergent.sig, Private = private.sig)
tier_outlier_snps <- function(tier, trt) {
  d <- tier_objs[[tier]]; unique(d$SNP[as.logical(d[[paste0("Sig.", trt)]])])
}
tier_sets_CASE <- lapply(names(tier_objs), function(ti) tier_outlier_snps(ti, "CASE"))
names(tier_sets_CASE) <- names(tier_objs)

case_only_snps <- add_panels[["CASE only"]]

tier_lookup <- rbindlist(lapply(names(tier_sets_CASE), function(ti)
  data.table(snp = tier_sets_CASE[[ti]], tier = ti)))
dup_snps <- tier_lookup[, .N, by = snp][N > 1, snp]
tier_lookup[snp %in% dup_snps, tier := "Multiple tiers"]
tier_lookup <- unique(tier_lookup, by = "snp")

case_only_dt <- merge(data.table(snp = case_only_snps), tier_lookup, by = "snp", all.x = TRUE)
case_only_dt[is.na(tier), tier := "Unclassified"]

cat("CASE-only loci by tier:\n")
```

    ## CASE-only loci by tier:

``` r
print(case_only_dt[, .N, by = tier][order(-N)])
```

    ##          tier     N
    ##        <char> <int>
    ## 1:       Core  1317
    ## 2: Convergent   380
    ## 3:    Private   244

``` r
tier_lev <- c("Core", "Convergent", "Private", "Multiple tiers", "Unclassified")
tier_lev <- tier_lev[tier_lev %in% unique(case_only_dt$tier)]   # drop empty buckets

pl_tier <- merge(add_pl[SNP %in% case_only_snps, .(snp = SNP, x, y, n_blk)], case_only_dt, by = "snp")
pl_tier[, tier := factor(tier, levels = tier_lev)]

per_blk_tier <- merge(add_perblk[SNP %in% case_only_snps, .(snp = SNP, block, x, y)], case_only_dt, by = "snp")
per_blk_tier[, tier := factor(tier, levels = tier_lev)]

fitdt_tier <- pl_tier[, as.list(deming(x, y)), by = tier]
tier_summary <- pl_tier[, {
  f <- sign(y); f[f == 0] <- 1
  .(n = .N, slope = unname(deming(x, y)["slope"]),
    mean_excess = mean(f * y - f * x),
    frac_supra  = mean(f * y > f * x))
}, by = tier][order(tier)]
print(tier_summary)
```

    ##          tier     n     slope mean_excess frac_supra
    ##        <fctr> <int>     <num>       <num>      <num>
    ## 1:       Core  1317 0.9337480  0.01574564  0.6431283
    ## 2: Convergent   380 1.4763767  0.03250475  0.8052632
    ## 3:    Private   244 0.9310882  0.04529844  0.8770492

``` r
fwrite(pl_tier, "./results/CASE_only_additivity_by_tier.csv")
fwrite(tier_summary, "./results/CASE_only_additivity_tier_summary.csv")
```

``` r
add_tier <- rbindlist(list(
  data.table(SNP = core.sig$SNP,       tier = "Core", BLOCK = core.sig$BLOCK),
  data.table(SNP = convergent.sig$SNP, tier = "Convergent", BLOCK = convergent.sig$BLOCK),
  data.table(SNP = private.sig$SNP,    tier = "Private", BLOCK = private.sig$BLOCK)))[!duplicated(SNP)]
add_pl_t <- merge(add_pl, add_tier, by = "SNP")
add_pl_t[, tier := factor(tier, levels = c("Core","Convergent","Private"))]

add_pl_t <- add_pl_t[sCASE & !sCA & !sSE]

add_pl_t$BLOCK <- factor(add_pl_t$BLOCK,levels = c(10, 11, 12),labels = c("B10", "B11", "B12"))


fit_t <- add_pl_t[, {
  ols  <- coef(lm(y ~ x))
  orth <- deming(x, y)
  f    <- sign(y); f[f == 0] <- 1
  .(n = .N,
    ols_int    = ols[1],
    ols_slope  = ols[2],
    orth_int   = orth["intercept"],
    orth_slope = orth["slope"],
    pct_supra  = 100 * mean(f * y > f * x),
    r          = cor(x, y))
}, by = tier][order(tier)]
print(fit_t)
```

    ##          tier     n     ols_int ols_slope     orth_int orth_slope pct_supra
    ##        <fctr> <int>       <num>     <num>        <num>      <num>     <num>
    ## 1:       Core  1317 0.007671365 0.5871317  0.005754033  0.9337480  64.31283
    ## 2: Convergent   380 0.003928169 0.6893050 -0.009998755  1.4763767  80.52632
    ## 3:    Private   244 0.052584857 0.6720712  0.042417267  0.9310882  87.70492
    ##            r
    ##        <num>
    ## 1: 0.6123126
    ## 2: 0.5535231
    ## 3: 0.7068777

``` r
fwrite(fit_t, file.path(tab_dir, "FigS5_additivity_by_tier.csv"))

tier_summary <- pl_tier[, {
  f <- sign(y); f[f == 0] <- 1
  .(n = .N, slope = unname(deming(x, y)["slope"]),
    mean_excess = mean(f * y - f * x),
    frac_supra  = mean(f * y > f * x))
}, by = tier][order(tier)]
print(tier_summary)
```

    ##          tier     n     slope mean_excess frac_supra
    ##        <fctr> <int>     <num>       <num>      <num>
    ## 1:       Core  1317 0.9337480  0.01574564  0.6431283
    ## 2: Convergent   380 1.4763767  0.03250475  0.8052632
    ## 3:    Private   244 0.9310882  0.04529844  0.8770492

``` r
rng_t <- quantile(c(add_pl_t$x, add_pl_t$y), c(0.001, 0.999), na.rm = TRUE)
ann_t  <- fit_t[, .(tier, lab = sprintf("OLS %.2f | orth %.2f\n%.0f%% above | r=%.2f",
                                        ols_slope, orth_slope, pct_supra, r))]
figS5_add <- ggplot(add_pl_t, aes(x, y)) +
  geom_bin2d(data = bg, bins = 200) +
  geom_point(data = add_pl_t, aes(colour = BLOCK), size = 1.5, alpha = 0.5) +
  scale_colour_manual(name = "block", values = spawn_cols) +
  scale_fill_gradient(name = sprintf("non-outlier loci\n(n=%d subsample)", length(bg_snps)),
                      low = "grey92", high = "grey55", trans = "log10") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
  geom_abline(data = fit_t, aes(slope = ols_slope,  intercept = ols_int),  colour = "#D55E00") +
  geom_abline(data = fit_t, aes(slope = orth_slope, intercept = orth_int), colour = "#0072B2", linetype = "dotted") +
  geom_text(data = ann_t, aes(label = lab), x = rng_t[1], y = rng_t[2],
            hjust = 0, vjust = 1, size = 2.8, inherit.aes = FALSE) +
  facet_wrap(~ tier) + coord_equal(xlim = rng_t, ylim = rng_t) +
  labs(x = "(ΔAF CA − ΔAF CON) + (ΔAF SE − ΔAF CON)", y = "ΔAF CASE − ΔAF CON",
       title = "Fig S5 — control-anchored additivity by replication tier",
       subtitle = "dashed = 1:1;  solid = OLS;  dotted = orthogonal") +
  theme_case(12)
save_supp(figS5_add, "FigS5_additivity_by_tier.png", width = 11, height = 4.5)
figS5_add
```

![](Final_reproducible_analysis_files/figure-gfm/figS5-additivity-tier-1.png)<!-- -->

``` r
bg_tier <- rbindlist(lapply(tier_lev, function(L) { d <- add_pl[SNP %in% bg_snps, .(x, y)]; d[, tier := L]; d }))
bg_tier[, tier := factor(tier, levels = tier_lev)]

ann_tier <- fitdt_tier[tier_summary, on = "tier"][
  , lab := sprintf("slope %.2f\n%.0f%% above 1:1", slope, 100 * frac_supra)]

ggplot(mapping = aes(x, y)) +
  geom_bin2d(data = bg_tier, bins = 80) +
  scale_fill_gradient(name = sprintf("non-outlier loci\n(n=%d subsample)", length(bg_snps)),
                       low = "grey92", high = "grey55", trans = "log10") +
  geom_point(data = per_blk_tier, aes(colour = block), size = 0.6, alpha = 0.6) +
  scale_colour_manual(name = "block", values = spawn_cols) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
  geom_abline(data = fitdt_tier, aes(slope = slope, intercept = intercept), colour = "black",
              linewidth = 0.6) +
  geom_text(data = ann_tier, aes(label = lab), x = rng[1], y = rng[2],
            hjust = 0, vjust = 1, size = 3, lineheight = 0.9, inherit.aes = FALSE) +
  facet_wrap(~ tier, nrow = 1) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = "(dCA − dCON) + (dSE − dCON)   [additive expectation]",
       y = "dCASE − dCON   [combined response]",
       title = "CASE-only additivity, faceted by locus replication tier",
       subtitle = "dashed = additive 1:1;  solid black = per-tier orthogonal fit;  grey = same 10,000-locus background;  colour = block") +
  theme_bw()
```

![](Final_reproducible_analysis_files/figure-gfm/hC-tier-plot-1.png)<!-- -->

# Fig 4 — Compositional synergy (novel targets)

## 4.1 Block-12 Manhattan, stacked across CA / SE / CASE

``` r
lighten_color <- function(col, factor = 0.5) {
  v <- col2rgb(col) / 255
  rgb(v[1] + (1 - v[1]) * factor, v[2] + (1 - v[2]) * factor, v[3] + (1 - v[3]) * factor)
}
make_tier_palette <- function(base_color)
  c("Core" = base_color, "Convergent" = lighten_color(base_color, 0.45),
    "Private" = lighten_color(base_color, 0.72))
build_tier_highlights <- function(dfm, block_num, treatment_col) {
  h.core <- core.sig       %>% dplyr::filter(BLOCK == block_num, asL(.data[[treatment_col]]))
  h.conv  <- convergent.sig  %>% dplyr::filter(BLOCK == block_num, asL(.data[[treatment_col]]))
  h.priv   <- private.sig %>% dplyr::filter(asL(.data[[treatment_col]]))
  dplyr::bind_rows(
    dplyr::semi_join(dfm, h.core, by = "SNP") %>% dplyr::mutate(Tier = "Core"),
    dplyr::semi_join(dfm, h.conv,  by = "SNP") %>% dplyr::mutate(Tier = "Convergent"),
    dplyr::semi_join(dfm, h.priv,   by = "SNP") %>% dplyr::mutate(Tier = "Private")
  ) %>% dplyr::mutate(Tier = factor(Tier, levels = c("Core","Convergent","Private")))
}
manhattan_block <- function(pv_df, block_num, treatment_col, qcolumn, base_color,
                            title = NULL, ylim = NULL, alpha_max = 20) {
  tier_fill <- make_tier_palette(base_color)
  con_b <- pv_df$SNP[pv_df$QCON < 0.1]
  pp <- as.data.frame(pv_df)
  in_ac <- pp$SNP %in% all_con$SNP
  in_cf <- pp$SNP %in% con_filter$SNP
  in_sg <- pp$SNP %in% singleton$SNP
  pp$Group <- dplyr::case_when(pp$SNP %in% con_b ~ "CON Filtered",
    in_ac ~ "CON Filtered", in_cf ~ "CON Filtered",
    in_sg ~ "Single-Block", TRUE ~ "Shared SNP")
  man <- ggman(pp, pvalue = qcolumn, chr = "CHR", snp = "SNP", bp = "BP",
               relative.positions = TRUE, sigLine = NA, pointSize = 1.5, title = "", alpha = 0.5)
  dfm <- man[[1]]; dfm$chrom_alt <- factor(dfm$chrom_alt, levels = c(0, 1))
  dfmo <- build_tier_highlights(dfm, block_num, treatment_col)
  dfmsplit <- split(dfm, dfm$chrom)
  xbreaks <- sapply(dfmsplit, function(x) { m <- length(x$index)/2; if (m < 1) m <- 1; x$index[m] })
  ysc <- if (is.null(ylim))
    scale_y_continuous(expand = c(0,0), limits = c(0, max(dfm$marker + 1, na.rm = TRUE)))
  else scale_y_continuous(expand = c(0,0), limits = ylim)
  if (is.null(title)) title <- paste0("Block ", block_num)
  ggplot(dfm, aes(index, marker, colour = as.factor(chrom_alt))) +
    geom_point(aes(alpha = marker, shape = Group, size = marker)) +
    geom_point(data = subset(dfmo, Tier == "Core"),           aes(index, marker, fill = Tier), shape = 21, size = 3, colour = "black", alpha = 0.9) +
    geom_point(data = subset(dfmo, Tier == "Convergent"),     aes(index, marker, fill = Tier), shape = 23, size = 3, colour = "black", alpha = 0.9) +
    geom_point(data = subset(dfmo, Tier == "Private"),        aes(index, marker, fill = Tier), shape = 25, size = 3, colour = "black", alpha = 0.9) +
    scale_fill_manual(values = tier_fill, name = "Locus Tier") +
    scale_x_continuous(breaks = xbreaks, labels = names(xbreaks), expand = c(0,0), limits = c(0, max(dfm$index) + 10)) +
    ysc + scale_size_continuous(range = c(0.25, 3)) +
    scale_alpha_continuous(range = c(0.05, 0.8), limits = c(0, alpha_max)) +
    scale_shape_manual(values = c(15, 21, 16)) +
    scale_color_manual(values = c("grey", "dark grey"), guide = "none") +
    guides(colour = "none", alpha = "none", size = "none",
           fill  = guide_legend(override.aes = list(shape = c(21, 23, 25)), size =4),
           shape = guide_legend(override.aes = list(color = "grey", size = 4, linetype = 0))) +
    # title as a top-left annotation (no plot title) so it lines up with the GO stack's labels
    annotate("text", x = -Inf, y = Inf, label = title, hjust = -0.1, vjust = 1.4,
             fontface = "bold", size = 4) +
    labs(x = "Chromosome", y = expression(-Log10(q-value))) +
    theme_case(13)+
    theme(plot.margin = margin(1, 1, 0, 0), legend.direction = "horizontal", legend.position = "inside", legend.position.inside = c(0.5,0.9), legend.box = "horizontal" , legend.box.margin = margin(t = -0, r = 0, b = 0, l = 0, unit = "pt") , legend.key.size = unit(0.25, "lines"),legend.title = element_text(size = 10),legend.text = element_text(size = 8) )
}
```

``` r
# y capped at 20 for all three; the two CHR5 points above 20 are control-filtered and
# noted in the legend rather than shown.
m12_ca   <- manhattan_block(pv_list$B12, 12, "Sig.CA",   "QCA",   "#E69F00", title = "B12 — CA",   ylim = c(0, 20))
m12_se   <- manhattan_block(pv_list$B12, 12, "Sig.SE",   "QSE",   "#0072B2", title = "B12 — SE",   ylim = c(0, 20))
m12_case <- manhattan_block(pv_list$B12, 12, "Sig.CASE", "QCASE", "#009E73", title = "B12 — CASE", ylim = c(0, 20))
fig3_man <- (m12_ca / m12_se / m12_case) & theme(plot.margin = margin(0, 0, 0, 0))
fig3_man
```

![](Final_reproducible_analysis_files/figure-gfm/fig3-manhattan-1.png)<!-- -->

## 4.1b Per-spawn Manhattan stacks (B10, B11) — supplement

Same stacked layout as Fig 4A (B12), for the other two sequenced spawns.
y capped at 20 for cross-spawn comparability with the main figure.

``` r
mk_manhattan_stack <- function(block_pv, bn, ymax) {
  ca   <- manhattan_block(block_pv, bn, "Sig.CA",   "QCA",   "#E69F00", title = paste0("B", bn, " — CA"),   ylim = c(0, ymax))
  se   <- manhattan_block(block_pv, bn, "Sig.SE",   "QSE",   "#0072B2", title = paste0("B", bn, " — SE"),   ylim = c(0, ymax))
  case <- manhattan_block(block_pv, bn, "Sig.CASE", "QCASE", "#009E73", title = paste0("B", bn, " — CASE"), ylim = c(0, ymax))
  (ca / se / case) & theme(plot.margin = margin(10, 0, 0, 0))
}
figS5 <- mk_manhattan_stack(pv_list$B10, 10,40)
save_supp(figS5, "FigS7_B10_manhattan_3stressor.png", width = 11, height = 10)
figS5
```

![](Final_reproducible_analysis_files/figure-gfm/figS5-S6-manhattan-1.png)<!-- -->

``` r
figS6 <- mk_manhattan_stack(pv_list$B11, 11,12)
save_supp(figS6, "FigS8_B11_manhattan_3stressor.png", width = 11, height = 10)
figS6
```

![](Final_reproducible_analysis_files/figure-gfm/figS5-S6-manhattan-2.png)<!-- -->

## 4.2 GO enrichment — custom plot from the TopGO CSVs

``` r
# Built directly from GO/<stressor>_GO_en_sig_gene.csv (topGO): one point per enriched GO
# term, x grouped by GO category, y = -log10(Fisher p), colour = treatment, shape = category.
# No gprofiler / gostres dependency.
go_files <- c(CA = "CA_GO_en_sig_gene.csv", SE = "SE_GO_en_sig_gene.csv", CASE = "CASE_GO_en_sig_gene.csv")
# find the GO dir regardless of the current working directory (analysis/, project root, etc.)
go_dir <- "./GO"
godat <- rbindlist(lapply(names(go_files), function(tr) {
    d <- fread(file.path(go_dir, go_files[[tr]])); d[, stressor := tr]; d
  }), fill = TRUE)
  godat <- unique(godat[, .(stressor, GO = GeneOntologyIDs, Term, Fisher, ontology)])
  godat[, cat := data.table::fcase(
    grepl("Biolog", ontology), "GO:BP",
    grepl("Molec",  ontology), "GO:MF",
    grepl("Cellul", ontology), "GO:CC", default = NA_character_)]
  godat[, Fisher := suppressWarnings(as.numeric(Fisher))]
  godat <- godat[!is.na(cat) & is.finite(Fisher) & Fisher > 0]
  godat[, neglogp := -log10(Fisher)]
  godat[, stressor := factor(stressor, levels = c("CA","SE","CASE"))]
  godat[, cat := factor(cat, levels = c("GO:BP","GO:CC","GO:MF"))]
  # x position from a GLOBAL term order (by category, then GO id) so the same term sits at the
  # same x in every treatment facet — a Manhattan-style layout sorted by term, not by p.
  # Add a visible gap between GO:BP / GO:CC / GO:MF blocks so the grouping is obvious without
  # needing a shape legend at all.
  terms <- unique(godat[, .(GO, cat)]); setorder(terms, cat, GO)
  terms[, idx := seq_len(.N), by = cat]
  n_per_cat <- terms[, .N, by = cat][order(cat)]
  gap <- max(2, round(max(n_per_cat$N) * 0.08))
  cat_offset <- setNames(c(0, head(cumsum(n_per_cat$N + gap), -1)), n_per_cat$cat)
  terms[, x := idx + cat_offset[as.character(cat)]]
  godat  <- merge(godat, terms[, .(GO, x)], by = "GO")
  catmid <- terms[, .(mid = mean(x)), by = cat]
  # boundaries between groups, for the dashed separator lines
  catbounds <- terms[, .(maxx = max(x)), by = cat][order(cat)]
  sep_x <- head(catbounds$maxx, -1) + gap / 2
  xr     <- range(terms$x)
  ymax_go <- max(godat$neglogp, na.rm = TRUE)
  # single faceted plot (one row per treatment) instead of a nested patchwork stack --
  # facet panels are guaranteed equal-width and fill the plot area, avoiding the
  # wrap_elements() sizing issues we kept hitting with a 3-plot stack.
  label_df <- data.table(stressor = factor(c("CA","SE","CASE"), levels = c("CA","SE","CASE")),
                          x = -Inf, y = Inf, lab = c("CA Enriched GO Terms (all blocks)","SE Enriched GO Terms (all blocks)","CASE Enriched GO Terms (all blocks)"))
  go_panel <- ggplot(godat, aes(x, neglogp)) +
    geom_vline(xintercept = sep_x, linetype = "dashed", colour = "grey80", linewidth = 0.4) +
    geom_point(aes(shape = cat, fill = stressor), size = 1.8, alpha = 0.65, color= "black") +
    geom_text(data = label_df, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
              hjust = -0.2, vjust = 1.3, fontface = "bold", size = 4) +
    facet_wrap(~stressor, ncol = 1, axes = "all") +
    scale_fill_manual(values = unlist(treat_cols[c("CA","SE","CASE")])) +
    scale_shape_manual(values = c(`GO:BP` = 21, `GO:CC` = 24, `GO:MF` = 22), name = "GO category") +
    scale_x_continuous(breaks = catmid$mid, labels = catmid$cat, limits = xr) +
    scale_y_continuous(limits = c(1, ymax_go * 1.05)) +
    coord_cartesian(xlim = xr, expand = FALSE, clip = "off") +
    labs(x = NULL, y = expression(-log[10](Fisher~italic(p)))) +
    theme_case(12) +
    theme(axis.text = element_text(colour = "black"),
          strip.background = element_blank(), strip.text = element_blank(),
          legend.title = element_blank(), plot.margin = margin(18,0, 22, 0),
          axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),panel.spacing.y = unit(2, "lines")) +
    # no legend anywhere -- the gaps + dashed separators make the GO:BP/CC/MF
    # grouping obvious from the x-axis labels alone, and the bold CA/SE/CASE
    # labels (drawn in place of facet strips) identify each row.
    guides(shape = "none", colour = "none", fill = "none")
  cat("GO panel built from:", go_dir, " | terms:", nrow(godat), "\n")
```

    ## GO panel built from: ./GO  | terms: 458

``` r
go_panel
```

![](Final_reproducible_analysis_files/figure-gfm/fig3-go-1.png)<!-- -->

## 4.3 One Venn — SNP and gene counts together (custom ggplot)

``` r
# One triple Venn that carries BOTH counts per region: big number = SNPs, (parenthesised) =
# genes. Counts are computed the same inclusive way for both, so they reconcile. Drawn as a
# native ggplot (3 circles + region labels) so it drops cleanly into patchwork.
gt <- as.data.table(grouped_total.sig)
gt[, `:=`(CA = asL(Sig.CA), SE = asL(Sig.SE), CASE = asL(Sig.CASE))]

disj <- function(d) c(
  CA      = d[ CA & !SE & !CASE, .N], SE   = d[!CA &  SE & !CASE, .N], CASE = d[!CA & !SE & CASE, .N],
  CA_SE   = d[ CA &  SE & !CASE, .N], CA_CASE = d[ CA & !SE & CASE, .N], SE_CASE = d[!CA & SE & CASE, .N],
  ALL     = d[ CA &  SE &  CASE, .N])
snp_d <- disj(gt)

# Gene (LOC) regions: prefer the pipeline's gene lists so the counts match VenLOC.png exactly.
# These files are INCLUSIVE totals/intersections, so convert to disjoint via inclusion-exclusion.
loc_files <- c(CA = "Sig.loci.1.CA.LOC", SE = "Sig.loci.1.SE.LOC", CASE = "Sig.loci.1.CASE.LOC",
               CA_SE = "Sig.loci.1.CA.SE.LOC", CA_CASE = "Sig.loci.1.CA.CASE.LOC",
               SE_CASE = "Sig.loci.1.SE.CASE.LOC", ALL = "Sig.loci.1.SE.CASE.CA.LOC")
nloc <- function(f) nrow(unique(read.table(f, header = FALSE)))
cum <- vapply(loc_files, nloc, integer(1))
gene_d <- c(
  CA      = cum[["CA"]]   - cum[["CA_SE"]]   - cum[["CA_CASE"]] + cum[["ALL"]],
  SE      = cum[["SE"]]   - cum[["CA_SE"]]   - cum[["SE_CASE"]] + cum[["ALL"]],
  CASE    = cum[["CASE"]] - cum[["CA_CASE"]] - cum[["SE_CASE"]] + cum[["ALL"]],
  CA_SE   = cum[["CA_SE"]]   - cum[["ALL"]],
  CA_CASE = cum[["CA_CASE"]] - cum[["ALL"]],
  SE_CASE = cum[["SE_CASE"]] - cum[["ALL"]],
  ALL     = cum[["ALL"]])
cat("Venn SNP regions (CA,SE,CASE,CA&SE,CA&CASE,SE&CASE,all):\n"); print(snp_d)
```

    ## Venn SNP regions (CA,SE,CASE,CA&SE,CA&CASE,SE&CASE,all):

    ##      CA      SE    CASE   CA_SE CA_CASE SE_CASE     ALL 
    ##     375     761    1941      64      70     232      29

``` r
cat("Venn gene regions:\n"); print(gene_d)
```

    ## Venn gene regions:

    ##      CA      SE    CASE   CA_SE CA_CASE SE_CASE     ALL 
    ##     143     291     747      50     102     285     103

``` r
circ <- function(cx, cy, r, lab, n = 160) {
  t <- seq(0, 2*pi, length.out = n); data.frame(x = cx + r*cos(t), y = cy + r*sin(t), grp = lab)
}
# Layout: CA top-left, SE top-right, CASE bottom.
cdf <- rbind(circ(-0.6, 0.4, 0.95, "CA"), circ(0.6, 0.4, 0.95, "SE"), circ(0, -0.6, 0.95, "CASE"))
cdf$grp <- factor(cdf$grp, levels = c("CA","SE","CASE"))
regs <- c("CA","SE","CASE","CA_SE","CA_CASE","SE_CASE","ALL")
lab_pos <- data.frame(
  reg = regs,
  #     CA     SE    CASE  CA_SE CA_CASE SE_CASE ALL
  x = c(-1.00, 1.00, 0,     0,   -0.50,   0.50,   0),
  y = c( 0.62, 0.62, -1.00, 0.72, -0.12, -0.12,   0.08))
lab_pos$lab <- paste0(snp_d[as.character(lab_pos$reg)], "\n(", gene_d[as.character(lab_pos$reg)], ")")
cat_pos <- data.frame(x = c(-1.45, 1.45, 0), y = c(1.12, 1.12, -1.66),
                      lab = c("CA","SE","CASE"), grp = factor(c("CA","SE","CASE"), levels = c("CA","SE","CASE")))
venn_combined <- ggplot() +
  geom_polygon(data = cdf, aes(x, y, fill = grp, group = grp), alpha = 0.35, colour = "grey30") +
  geom_text(data = lab_pos, aes(x, y, label = lab), size = 2.5, lineheight = 0.85) +
  geom_text(data = cat_pos, aes(x, y, label = lab, colour = grp), fontface = "bold", size = 3.4) +
  scale_fill_manual(values = treat_cols[c("CA","SE","CASE")], guide = "none") +
  scale_colour_manual(values = treat_cols[c("CA","SE","CASE")], guide = "none") +
  coord_fixed(clip = "off", xlim = c(-1.8, 1.8), ylim = c(-1.8, 1.5)) + theme_void() +
  labs(title = "Outlier overlap",
       subtitle = "value = SNPs;  (value) = genes") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        plot.margin = margin(2, 2, 2, 2))
venn_combined
```

![](Final_reproducible_analysis_files/figure-gfm/fig3-venn-combined-1.png)<!-- -->

## 4.4b Mechanism panel D — data-driven stress GO clusters

``` r
sm <- melt(go_clusters_stress[, .(cluster, CA, SE, CASE, minF)],
           id.vars = c("cluster","minF"), variable.name = "Treatment", value.name = "n")
sm[, cluster   := factor(cluster, levels = rev(go_clusters_stress[order(minF), cluster]))]
sm[, Treatment := factor(Treatment, levels = c("CA","SE","CASE"))]
fig3d <- ggplot(sm, aes(Treatment, cluster, fill = Treatment, alpha = n)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = n), size = 3.4, colour = "grey15", alpha = 1) +
  scale_fill_manual(values = c(CA = "#E69F00", SE = "#0072B2", CASE = "#009E73"), guide = "none") +
  scale_alpha(range = c(0.15, 1), name = "enriched\nGO terms") +
  labs(x = NULL, y = NULL, title = "Mechanism by stressor — data-driven GO clusters") +
  theme_case(12) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 9))
fig3d
```

![](Final_reproducible_analysis_files/figure-gfm/fig3d-clusters-1.png)<!-- -->

## 4.5 Assemble Fig 4 (patchwork)

``` r
# Manhattan stack is a nested patchwork -> wrap so it stays one area. GO plot, combined
# Venn, and heatmap are native ggplots and drop in directly.
# Both the Manhattan and GO panels are nested 3-plot stacks -> wrap each as one area so
# their three sub-rows line up at equal height.
man_el <- patchwork::wrap_elements(fig3_man)
go_el <- patchwork::wrap_elements(full=go_panel)


# Layout (8 cols): top 4 rows = Manhattan (A, 5/8) | GO stack (B, 3/8)  [taller];
#                  bottom 2 rows = combined Venn (C, 4/8) | heatmap (D, 4/8) -- even split.


design3 <- "
AAB
AAB
AAB
CDD
"
fig3 <- man_el + go_el + venn_combined + fig3d +
  plot_layout(design = design3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.title = element_text(size = 12))
save_case(fig3, "Fig4_compositional_synergy.png", width = 16, height = 15)
#fig3
```

## Fig S6 — Per-spawn outlier overlap (Venns)

``` r
spawn_venn_grob <- function(b) {
  pv <- pv_list[[b]]
  sig <- function(q) pv[get(q) < alpha & QCON > alpha2, SNP]
  A <- sig("QCA"); S <- sig("QSE"); C <- sig("QCASE")
  grid::grid.grabExpr({
    grid::grid.newpage()
    vp <- VennDiagram::draw.triple.venn(
      area1 = length(A), area2 = length(S), area3 = length(C),
      n12 = length(intersect(A, S)), n13 = length(intersect(A, C)),
      n23 = length(intersect(S, C)), n123 = length(Reduce(intersect, list(A, S, C))),
      category = c("CA","SE","CASE"), fill = cbPaletteSmall3, col = cbPaletteSmall3,
      cat.col = cbPaletteSmall3, cex = 1.4, cat.cex = 1.5, margin = 0.08, ind = FALSE)
    grid::grid.draw(vp)
    grid::grid.text(b, y = 0.98, gp = grid::gpar(fontface = "bold", cex = 1.4))
  })
}
figS4 <- patchwork::wrap_elements(spawn_venn_grob("B10")) +
         patchwork::wrap_elements(spawn_venn_grob("B11")) +
         patchwork::wrap_elements(spawn_venn_grob("B12")) +
  plot_annotation(tag_levels = "A",
                  title = "Fig S6 — Per-spawn outlier overlap across stressors (CA / SE / CASE)")
save_supp(figS4, "FigS6_spawn_venns.png", width = 15, height = 5.5)
figS4
```

![](Final_reproducible_analysis_files/figure-gfm/figS4-spawn-venns-1.png)<!-- -->

## Fig S9 — Full GO enrichment by stressor (topGO)

One figure per category (ALL / CA / SE / CASE); the curated subset is in
Fig 4.

``` r
go_dir2 <- "./GO"
go_full_files  <- c(ALL = "ALL_GO_en_sig_gene.csv", CA = "CA_GO_en_sig_gene.csv",
                   SE = "SE_GO_en_sig_gene.csv", CASE = "CASE_GO_en_sig_gene.csv")
go_full_labels <- c(ALL = "All Outlier Loci", CA = "CA Stressor", SE = "SE Stressor",
                    CASE = "CASE Stressor (CA+SE)")
go_full_fignum <- c(ALL = "S9", CA = "S10", SE = "S11", CASE = "S12")
go_cat_cols <- c(`GO:BP` = "#4DAF4A", `GO:CC` = "#377EB8", `GO:MF` = "#984EA3", other = "grey60")
plot_go_full <- function(file, label, fignum) {
  d <- as.data.table(fread(file.path(go_dir2, file)))
  d[, Fisher := suppressWarnings(as.numeric(Fisher))]
  d <- d[is.finite(Fisher) & Fisher > 0]
  ggplot(d, aes(x = Term, y = -log(Fisher), size = prop.sig.genes, fill = -log(Fisher))) +
    expand_limits(y = 1.5) +
    geom_hline(yintercept = -log10(0.01), linetype = "longdash", colour = "black", linewidth = .6) +
    geom_hline(yintercept = -log10(0.001), linetype = "solid", colour = "black", linewidth = .6) +
    geom_point(shape = 21) +
    scale_size(range = c(0.5, 10)) +
    scale_fill_continuous(low = "#1AD3D1FF", high = "#4686FBFF", name = "-log(Fisher)") +
    xlab('') +
    ylab('Enrichment score') +
    ggtitle(paste0("Fig ", fignum, " — ", label, " enriched GO terms")) +
    theme_case(11) +
    theme(axis.text.y = element_text(size = 6)) +
    facet_grid(vars(ontology), scales = "free", space = "free_y") +
    coord_flip()
}
for (tr in names(go_full_files)) {
  p <- plot_go_full(go_full_files[[tr]], go_full_labels[[tr]], go_full_fignum[[tr]])
  save_supp(p, paste0("Fig", go_full_fignum[[tr]], "_GO_", tr, ".png"), width = 10, height = 12)
  print(p)
}
```

![](Final_reproducible_analysis_files/figure-gfm/figS7-go-full-1.png)<!-- -->![](Final_reproducible_analysis_files/figure-gfm/figS7-go-full-2.png)<!-- -->![](Final_reproducible_analysis_files/figure-gfm/figS7-go-full-3.png)<!-- -->![](Final_reproducible_analysis_files/figure-gfm/figS7-go-full-4.png)<!-- -->

## Fig S13 — GO-term overlap and biological-theme enrichment

``` r
# Collapse gene x term -> unique GO terms (best Fisher per term).
dedupe_go <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
    arrange(Fisher) %>%
    distinct(GeneOntologyIDs, .keep_all = TRUE)
}

# Vector of unique GO IDs for a treatment's GO data frame.
go_ids <- function(df) {
  d <- dedupe_go(df)
  if (is.null(d)) character(0) else unique(d$GeneOntologyIDs)
}

# Canonical ordering used across the summary figures.
onto_levels  <- c("Biological Processes", "Cellular Components", "Molecular Functions")
treat_levels <- c("CA", "SE", "CASE")
```

``` r
themes <- list(
  "Cilium / microtubule motility" = c("cilium","cilia","dynein","microtubule",
      "flagell","axonem","motile","myosin","actin","motor activity","motility"),
  "Ion transport / acid-base"     = c("proton","sodium","potassium","calcium",
      "ion transport","ion homeostasis","atpase","bicarbonate"),
  "Translation / ribosome"        = c("translat","ribosom","rrna","eif",
      "polysome","trna"),
  "Mitochondria / oxidative"      = c("mitochond","oxidative","respiratory",
      "electron transport","cytochrome"),
  "Ubiquitin-proteasome"          = c("ubiquitin","proteasom","erad",
      "catabolic process"),
  "Protein folding / chaperone"   = c("fold","chaperone","unfolded",
      "heat shock","hsp","prefoldin","chaperonin"),
  "Apoptosis / stress response"   = c("apopto","stress","death","autophag"),
  "DNA repair / replication"      = c("dna repair","dna replicat","damage",
      "recombinat","nucleotide excision")
)
theme_levels <- names(themes)

# Representative mollusc genes per theme (from GO_analysis_report.md, Section 4).
# Used ONLY for axis labels — turns the count heatmap into a mechanism figure.
# These are literature-derived exemplars, not the CASE candidate-gene list itself.
theme_genes <- c(
  "Cilium / microtubule motility" = "tubulin, tektin, dynein",
  "Ion transport / acid-base"     = "Na/K-ATPase, proton pump, amiloride ch.",
  "Translation / ribosome"        = "eIFs, Rps8/Rps4x",
  "Mitochondria / oxidative"      = "ATP synthase, MCU, SOD/GAPDH",
  "Ubiquitin-proteasome"          = "UPR (IRE1/PERK), ERAD, Ub ligases",
  "Protein folding / chaperone"   = "HSP70, HSP90, UPR",
  "Apoptosis / stress response"   = "IAPs, Bcl-2",
  "DNA repair / replication"      = ""
)
# Two-line axis label: theme on top, representative genes beneath (omitted if none).
theme_lab <- setNames(
  ifelse(theme_genes[theme_levels] == "", theme_levels,
         paste0(theme_levels, "\n(", theme_genes[theme_levels], ")")),
  theme_levels)

# Count unique GO terms matching each theme in one GO data frame.
theme_counts <- function(df) {
  d <- dedupe_go(df)
  if (is.null(d)) return(setNames(rep(0L, length(themes)), theme_levels))
  term_l <- tolower(d$Term)
  vapply(themes, function(kws) sum(Reduce(`|`, lapply(kws, grepl, term_l))),
         integer(1))
}

# Treatment palette — matches Ven.png / cbPaletteSmall3: CA orange, SE blue, CASE green.
treat_pal <- c(CA = "#E69F00", SE = "#0072B2", CASE = "#009E73")

# Shared heatmap builder (gene-annotated y axis).
# Tiles are coloured by treatment (`fill_col` hue) and shaded by count (alpha);
# the count is also printed, so exact values stay legible regardless of shade.
theme_heatmap <- function(df, x_col, x_levels, fill_col, palette, title_label) {
  df$theme     <- factor(df$theme, levels = rev(theme_levels),
                         labels = unname(rev(theme_lab)))
  df[[x_col]]  <- factor(df[[x_col]], levels = x_levels)
  ggplot(df, aes(x = .data[[x_col]], y = theme)) +
    geom_tile(aes(fill = .data[[fill_col]], alpha = n),
              colour = "white", linewidth = 1) +
    geom_text(aes(label = n), colour = "grey15", size = 4.5) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_alpha(range = c(0.18, 1), limits = c(0, NA), name = "Enriched\nGO terms") +
    labs(x = NULL, y = NULL, title = title_label) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot",
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(hjust = 1, lineheight = 0.95, size = 11))
}

# ---- 3a: theme x treatment (Total outlier set) ----
mat_treat <- sapply(list(CA = GO_CA, SE = GO_SE, CASE = GO_CASE), theme_counts)
df_treat <- as.data.frame(mat_treat, check.names = FALSE) %>%
  mutate(theme = theme_levels) %>%
  pivot_longer(-theme, names_to = "treatment", values_to = "n")
p_theme_treat <- theme_heatmap(df_treat, "treatment", treat_levels,
  fill_col = "treatment", palette = treat_pal,
  title_label = "Enriched GO terms by biological theme and stressor")
print(p_theme_treat)
```

![](Final_reproducible_analysis_files/figure-gfm/theme-heatmap-1.png)<!-- -->

``` r
ggsave("./GO/GO_theme_heatmap_treatment.png", p_theme_treat,
       width = 11, height = 6.5, dpi = 300, bg = "transparent")

# ---- 3b: theme x tier, for the CASE treatment ----
mat_tier <- sapply(list(Core = GO_Core_CASE, Convergent = GO_Conv_CASE,
                        "Private" = GO_Priv_CASE), theme_counts)
df_tier <- as.data.frame(mat_tier, check.names = FALSE) %>%
  mutate(theme = theme_levels) %>%
  pivot_longer(-theme, names_to = "tier", values_to = "n")
df_tier$fillgrp <- "CASE"   # tier panel is all CASE -> single CASE-green hue
p_theme_tier <- theme_heatmap(df_tier, "tier",
  c("Core", "Convergent", "Private"),
  fill_col = "fillgrp", palette = treat_pal,
  title_label = "CASE: enriched GO terms by theme and locus tier")
print(p_theme_tier)
```

![](Final_reproducible_analysis_files/figure-gfm/theme-heatmap-2.png)<!-- -->

``` r
ggsave("./GO/GO_theme_heatmap_CASE_by_tier.png", p_theme_tier,
       width = 10, height = 6.5, dpi = 300, bg = "transparent")
```

``` r
# Panel A: GO-term Venn (CA/SE/CASE); Panels B/C: theme x treatment and theme x tier heatmaps.
sA <- go_ids(GO_CA); sB <- go_ids(GO_SE); sC <- go_ids(GO_CASE)
venn_grob <- grid::grid.grabExpr({
  grid::grid.newpage()
  vp <- VennDiagram::draw.triple.venn(
    area1 = length(sA), area2 = length(sB), area3 = length(sC),
    n12 = length(intersect(sA, sB)), n13 = length(intersect(sA, sC)),
    n23 = length(intersect(sB, sC)), n123 = length(Reduce(intersect, list(sA, sB, sC))),
    category = c("CA","SE","CASE"), fill = cbPaletteSmall3, col = cbPaletteSmall3,
    cat.col = cbPaletteSmall3, cex = 1.4, cat.cex = 1.5, margin = 0.08, ind = FALSE)
  grid::grid.draw(vp)
})
# Panels B/C: data-driven cluster enrichment, all clusters, by stressor and by tier
clev <- rev(go_clusters_trt[order(minF), cluster])
bt <- melt(go_clusters_trt[, .(cluster, CA, SE, CASE)], id.vars = "cluster",
           variable.name = "Treatment", value.name = "n")
bt[, cluster := factor(cluster, levels = clev)]
bt[, Treatment := factor(Treatment, levels = c("CA","SE","CASE"))]
p_clust_trt <- ggplot(bt, aes(Treatment, cluster, fill = Treatment, alpha = n)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = n), size = 2.8, colour = "grey15", alpha = 1) +
  scale_fill_manual(values = c(CA="#E69F00", SE="#0072B2", CASE="#009E73"), guide = "none") +
  scale_alpha(range = c(0.12, 1), name = "GO terms") +
  labs(x = NULL, y = NULL, title = "By stressor") +
  theme_case(11) + theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8))
bk <- melt(go_clusters_tier[, .(cluster, Core, Convergent, Private)], id.vars = "cluster",
           variable.name = "Tier", value.name = "n")
bk[, cluster := factor(cluster, levels = clev)]
bk[, Tier := factor(Tier, levels = c("Core","Convergent","Private"))]
p_clust_tier <- ggplot(bk, aes(Tier, cluster, fill = Tier, alpha = n)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = n), size = 2.8, colour = "grey15", alpha = 1) +
  scale_fill_manual(values = tier_cols[c("Core","Convergent","Private")], guide = "none") +
  scale_alpha(range = c(0.12, 1), name = "GO terms") +
  labs(x = NULL, y = NULL, title = "By tier") +
  theme_case(11) + theme(panel.grid = element_blank(), axis.text.y = element_blank())
figS11 <- patchwork::wrap_elements(venn_grob) + p_clust_trt + p_clust_tier +
  plot_layout(widths = c(1, 1.1, 0.9)) +
  plot_annotation(tag_levels = "A",
                  title = "Fig S13 — GO-term overlap and data-driven cluster enrichment (by stressor and by tier)")
save_supp(figS11, "FigS13_GO_synthesis.png", width = 16, height = 7)
figS11
```

![](Final_reproducible_analysis_files/figure-gfm/figS11-assemble-1.png)<!-- -->

# Fig 5 — Diversity, architecture, and convergent function (CA / SE / CASE)

``` r
treat_levels <- c("CON","CA","SE","CASE")   # restored (the GO helpers above redefine it locally)
```

``` r
block_design <- list(b10 = list(gen = c(0,1,0,1,0,1),       repl = c(1,1,2,2,3,3)),
                     b11 = list(gen = c(0,1,0,1,0,1,0,1),   repl = c(1,1,2,2,3,3,4,4)),
                     b12 = list(gen = c(0,1,0,1,0,1,0,1),   repl = c(1,1,2,2,3,3,4,4)))
delta_af_table <- function(af_obj, block) {
  des <- block_design[[block]]; gen <- des$gen; repl <- des$repl
  m <- as.matrix(af_obj); reps <- unique(repl)
  t0 <- vapply(reps, function(r) m[, which(repl == r & gen == 0)], numeric(nrow(m)))
  dd <- vapply(reps, function(r) m[, which(repl == r & gen == 1)] - m[, which(repl == r & gen == 0)], numeric(nrow(m)))
  data.table(SNP = sub("\\.(?=[^.]+$)", "_", rownames(m), perl = TRUE),
             T0 = rowMeans(t0, na.rm = TRUE), dAF = rowMeans(dd, na.rm = TRUE))
}
grid_tb <- CJ(trt = c("CA","CASE","SE"), blk = c("b10","b11","b12"), sorted = FALSE)  # CON dropped
master <- rbindlist(lapply(seq_len(nrow(grid_tb)), function(i) {
  tab <- delta_af_table(get(paste0(tolower(grid_tb$trt[i]), ".", grid_tb$blk[i], ".af")), grid_tb$blk[i])
  tab[, `:=`(treatment = grid_tb$trt[i], block = toupper(grid_tb$blk[i]))]; tab }))
master[, `:=`(MAF0 = pmin(T0, 1 - T0), abs_dAF = abs(dAF),
              treatment = factor(treatment, levels = c("CA","SE","CASE")))]
tier_of <- rbindlist(list(data.table(SNP = core.sig$SNP, tier = "Core"),
                          data.table(SNP = convergent.sig$SNP, tier = "Convergent"),
                          data.table(SNP = private.sig$SNP, tier = "Private")))[!duplicated(SNP)]
master <- merge(master, tier_of, by = "SNP", all.x = TRUE)
master[is.na(tier), tier := "Background"][, tier := factor(tier, levels = tier_levels)]
```

``` r
# Diversity at selected loci, by tier (CON-anchored). For each treatment's outlier loci in a tier,
# the excess change in expected heterozygosity = (treatment post-initial He) - (CON post-initial
# He) at the SAME loci. He_site = (c/(c-1)) * 2pq (read-count corrected). CON-anchoring removes the
# initial-vs-post offset and the ascertainment shared with CON, so a positive excess means the
# treatment raised diversity at its targets (rare variants rising toward intermediate frequency).
# Built from the raw per-treatment af/cov objects (master drops CON, so it cannot be used here).
he_min_cov <- 2; he_n0 <- 40000; he_n1 <- 10000
he_site_corr <- function(p, c, n) (c / (c - 1)) * 2 * p * (1 - p)
he_to_us <- function(x) sub("\\.(?=[^.]+$)", "_", x, perl = TRUE)   # CHROM.POS -> CHROM_POS
he_asL   <- function(x) as.logical(x)
he_tier_objs <- list(Core = core.sig, Convergent = convergent.sig, Private = private.sig)
he_tier_snps <- function(tier, trt) {
  d <- he_tier_objs[[tier]]; unique(d$SNP[he_asL(d[[paste0("Sig.", trt)]])])
}
he_pair <- function(af, cov, i0, i1, idx) {
  p0 <- af[idx, i0]; c0 <- cov[idx, i0]; p1 <- af[idx, i1]; c1 <- cov[idx, i1]
  ok <- is.finite(p0) & is.finite(p1) & is.finite(c0) & is.finite(c1) &
        c0 >= he_min_cov & c1 >= he_min_cov
  if (!any(ok)) return(c(NA_real_, NA_real_))
  c(mean(he_site_corr(p0[ok], c0[ok], he_n0)), mean(he_site_corr(p1[ok], c1[ok], he_n1)))
}
he_rows <- list()
for (blk in c("b10", "b11", "b12")) {
  des <- block_design[[blk]]; gen <- des$gen; repl <- des$repl
  afC <- as.matrix(get(paste0("con.", blk, ".af")));  covC <- as.matrix(get(paste0("con.", blk, ".cov")))
  rnC <- he_to_us(rownames(afC))
  for (tr in c("CA", "SE", "CASE")) {
    afT  <- as.matrix(get(paste0(tolower(tr), ".", blk, ".af")))
    covT <- as.matrix(get(paste0(tolower(tr), ".", blk, ".cov")))
    rnT  <- he_to_us(rownames(afT))
    for (tier in names(he_tier_objs)) {
      common <- intersect(intersect(rnT, rnC), he_tier_snps(tier, tr))
      if (!length(common)) next
      idxT <- match(common, rnT); idxC <- match(common, rnC)
      dT <- dC <- numeric(0)
      for (r in unique(repl)) {
        i0 <- which(repl == r & gen == 0); i1 <- which(repl == r & gen == 1)
        if (length(i0) != 1 || length(i1) != 1) next
        mT <- he_pair(afT, covT, i0, i1, idxT); dT <- c(dT, mT[2] - mT[1])
        mC <- he_pair(afC, covC, i0, i1, idxC); dC <- c(dC, mC[2] - mC[1])
      }
      he_rows[[length(he_rows) + 1]] <- data.table(
        block = toupper(blk), treatment = tr, tier = tier,
        n_out = length(common), excess = mean(dT, na.rm = TRUE) - mean(dC, na.rm = TRUE))
    }
  }
}
he_tier_dt <- rbindlist(he_rows)
he_tier_dt[, treatment := factor(treatment, levels = c("CA", "SE", "CASE"))]
he_tier_dt[, tier := factor(tier, levels = c("Core", "Convergent", "Private"))]
fwrite(he_tier_dt, file.path(tab_dir, "Fig5_He_by_tier.csv"))

y.expression <- expression(H[E])

fig4_he_panel <- ggplot(he_tier_dt, aes(treatment, excess, fill = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(alpha = 0.55, outlier.shape = NA, colour = "black", width = 0.6) +
  geom_jitter(aes(shape = block), width = 0.12, height = 0, size = 2.2, alpha = 0.85) +
  facet_wrap(~ tier, axes = "all") +
  scale_fill_manual(values = unlist(treat_cols[c("CA","SE","CASE")]), name = "Treatment") +
  scale_shape_manual(values = spawn_shapes, name = "Block") +
  labs(x = NULL, y = expression("Excess" ~ italic(H[E]) ~ "compared to CON")) +
  theme_case(12) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.direction = "horizontal", legend.position = "bottom",panel.spacing = unit(1, "lines"), axis.title.y = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 0)), plot.margin = margin(t = -0, r = 0, b = 0, l = -0, unit = "pt"))
fig4_he_panel
```

![](Final_reproducible_analysis_files/figure-gfm/fig4-he-tier-1.png)<!-- -->

``` r
# Frequency/shift-matched neutral null for the CON-anchored He excess, per tier x treatment,
# using the SAME read-count-corrected He and coverage floor as the observed panel above. Overlays
# the null band (2.5-97.5%) and z on the He panel; matching diagnostics -> Fig S14.
he_B <- 1000; he_maf_breaks <- seq(0, 0.5, by = 0.05); he_cal <- 0.15; he_adp_floor <- 0.01
set.seed(1)
he_all_out <- unique(unlist(lapply(names(he_tier_objs), function(ti)
  lapply(c("CA","SE","CASE"), function(tr) he_tier_snps(ti, tr)))))

he_precompute <- function(tr, blk) {
  des <- block_design[[blk]]; gen <- des$gen; repl <- des$repl
  afT <- as.matrix(get(paste0(tolower(tr), ".", blk, ".af"))); covT <- as.matrix(get(paste0(tolower(tr), ".", blk, ".cov")))
  afC <- as.matrix(get(paste0("con.", blk, ".af")));           covC <- as.matrix(get(paste0("con.", blk, ".cov")))
  rnT <- he_to_us(rownames(afT)); rnC <- he_to_us(rownames(afC))
  common <- intersect(rnT, rnC); idxT <- match(common, rnT); idxC <- match(common, rnC); reps <- sort(unique(repl))
  buildHe <- function(af, cov, idx) {
    H0 <- H1 <- matrix(NA_real_, length(idx), length(reps))
    for (ri in seq_along(reps)) { r <- reps[ri]
      i0 <- which(repl == r & gen == 0); i1 <- which(repl == r & gen == 1)
      if (length(i0) != 1 || length(i1) != 1) next
      p0 <- af[idx, i0]; c0 <- cov[idx, i0]; p1 <- af[idx, i1]; c1 <- cov[idx, i1]
      ok <- is.finite(p0) & is.finite(p1) & is.finite(c0) & is.finite(c1) & c0 >= he_min_cov & c1 >= he_min_cov
      h0 <- h1 <- rep(NA_real_, length(idx))
      h0[ok] <- he_site_corr(p0[ok], c0[ok], he_n0); h1[ok] <- he_site_corr(p1[ok], c1[ok], he_n1)
      H0[, ri] <- h0; H1[, ri] <- h1 }
    list(H0 = H0, H1 = H1) }
  Tm <- buildHe(afT, covT, idxT); Cm <- buildHe(afC, covC, idxC)
  i0s <- which(gen == 0); i1s <- which(gen == 1)
  p0m <- rowMeans(afT[idxT, i0s, drop = FALSE], na.rm = TRUE); p1m <- rowMeans(afT[idxT, i1s, drop = FALSE], na.rm = TRUE)
  list(common = common, HT0 = Tm$H0, HT1 = Tm$H1, HC0 = Cm$H0, HC1 = Cm$H1, maf0 = pmin(p0m, 1 - p0m), adp = abs(p1m - p0m))
}
he_excess_S <- function(pc, S) {
  if (length(S) < 1) return(NA_real_)
  mT <- colMeans(pc$HT1[S, , drop = FALSE], na.rm = TRUE) - colMeans(pc$HT0[S, , drop = FALSE], na.rm = TRUE)
  mC <- colMeans(pc$HC1[S, , drop = FALSE], na.rm = TRUE) - colMeans(pc$HC0[S, , drop = FALSE], na.rm = TRUE)
  mean(mT, na.rm = TRUE) - mean(mC, na.rm = TRUE)
}
hePC <- list()
for (tr in c("CA","SE","CASE")) { hePC[[tr]] <- list(); for (blk in c("b10","b11","b12")) hePC[[tr]][[blk]] <- he_precompute(tr, blk) }

he_null_rows <- list(); he_diag_rows <- list(); he_bk <- c("b10","b11","b12")
for (tier in names(he_tier_objs)) for (tr in c("CA","SE","CASE")) {
  obs_b <- rep(NA_real_, 3); null_b <- matrix(NA_real_, he_B, 3)
  for (bi in seq_along(he_bk)) {
    blk <- he_bk[bi]; pc <- hePC[[tr]][[blk]]; common <- pc$common
    out_snps <- intersect(common, he_tier_snps(tier, tr)); neu_snps <- setdiff(common, he_all_out)
    if (length(out_snps) < 3 || length(neu_snps) < 1) next
    out_idx <- match(out_snps, common); neu_idx <- match(neu_snps, common)
    mafb <- cut(pc$maf0, he_maf_breaks, include.lowest = TRUE)
    pools <- lapply(out_idx, function(o) { tol <- max(he_cal * pc$adp[o], he_adp_floor)
      neu_idx[ mafb[neu_idx] == mafb[o] & abs(pc$adp[neu_idx] - pc$adp[o]) <= tol ] })
    matched <- which(lengths(pools) > 0)
    he_diag_rows[[length(he_diag_rows) + 1]] <- data.table(tier = tier, treatment = tr, block = toupper(blk),
      n_out = length(out_idx), n_matched = length(matched),
      pct_matched = round(100 * length(matched) / length(out_idx), 1))
    if (length(matched) < 3) next
    obs_b[bi] <- he_excess_S(pc, out_idx[matched])
    DRAW <- matrix(NA_integer_, length(matched), he_B)
    for (kk in seq_along(matched)) { pool <- pools[[matched[kk]]]; DRAW[kk, ] <- pool[sample.int(length(pool), he_B, replace = TRUE)] }
    for (b in seq_len(he_B)) { s <- DRAW[, b]; s <- s[!is.na(s)]; if (length(s) >= 3) null_b[b, bi] <- he_excess_S(pc, s) }
  }
  if (all(is.na(obs_b))) next
  obs <- mean(obs_b, na.rm = TRUE); nulld <- rowMeans(null_b, na.rm = TRUE); nulld <- nulld[is.finite(nulld)]
  if (!length(nulld)) next
  he_null_rows[[length(he_null_rows) + 1]] <- data.table(tier = tier, treatment = tr, obs = obs,
    null_mean = mean(nulld), null_sd = sd(nulld), z = (obs - mean(nulld)) / sd(nulld),
    band_lo = unname(quantile(nulld, 0.025)), band_hi = unname(quantile(nulld, 0.975)))
}
he_null_dt <- rbindlist(he_null_rows)
he_null_dt[, treatment := factor(treatment, levels = c("CA","SE","CASE"))]
he_null_dt[, tier := factor(tier, levels = c("Core","Convergent","Private"))]
he_diag_dt <- rbindlist(he_diag_rows)
he_diag_dt[, treatment := factor(treatment, levels = c("CA","SE","CASE"))]
he_diag_dt[, tier := factor(tier, levels = c("Core","Convergent","Private"))]
fwrite(he_null_dt, file.path(tab_dir, "Fig5_He_matched_null.csv"))
print(he_null_dt)
```

    ##          tier treatment        obs  null_mean     null_sd        z      band_lo
    ##        <fctr>    <fctr>      <num>      <num>       <num>    <num>        <num>
    ## 1:       Core        CA 0.01906261 0.01108094 0.002275491 3.507666  0.006683559
    ## 2:       Core        SE 0.01726881 0.01054215 0.001527152 4.404707  0.007448281
    ## 3:       Core      CASE 0.01593089 0.01067984 0.001123017 4.675847  0.008472558
    ## 4: Convergent        CA 0.02826545 0.01080095 0.007177404 2.433261 -0.003894264
    ## 5: Convergent        SE 0.04556717 0.02365117 0.004316883 5.076812  0.014812858
    ## 6: Convergent      CASE 0.02727507 0.01513129 0.001822567 6.663012  0.011596394
    ## 7:    Private        CA 0.09432384 0.05869653 0.006210339 5.736773  0.046551763
    ## 8:    Private        SE 0.10975845 0.06871789 0.005345037 7.678256  0.058025409
    ## 9:    Private      CASE 0.10997779 0.06891336 0.004420605 9.289324  0.060265732
    ##       band_hi
    ##         <num>
    ## 1: 0.01548018
    ## 2: 0.01365109
    ## 3: 0.01288355
    ## 4: 0.02531708
    ## 5: 0.03148379
    ## 6: 0.01880278
    ## 7: 0.07062250
    ## 8: 0.07963855
    ## 9: 0.07742953

``` r
# overlay the matched-null band (2.5-97.5%) and z on the observed He panel (rebuilds fig4_he_panel)
fig4_he_panel <- fig4_he_panel +
  geom_crossbar(data = he_null_dt, aes(x = treatment, y = null_mean, ymin = band_lo, ymax = band_hi),
                inherit.aes = FALSE, width = 0.5, fill = "grey80", colour = "grey55", alpha = 0.45) +
  geom_text(data = he_null_dt, aes(x = treatment, y = band_hi, label = sprintf("z=%.1f", z)),
            inherit.aes = FALSE, vjust = -0.5, size = 2.6)
fig4_he_panel
```

![](Final_reproducible_analysis_files/figure-gfm/fig5-he-null-1.png)<!-- -->

``` r
figS13 <- (ggplot(he_null_dt, aes(treatment, z, fill = treatment)) +
             geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
             geom_col(alpha = 0.85) + facet_wrap(~ tier) +
             scale_fill_manual(values = unlist(treat_cols[c("CA","SE","CASE")])) +
             labs(title = "Fig S14A — He excess relative to the frequency/shift-matched null (z)",
                  x = NULL, y = "z") + theme_case(11) + theme(legend.position = "none")) /
           (ggplot(he_diag_dt, aes(treatment, pct_matched, fill = block)) +
             geom_col(position = "dodge") + facet_wrap(~ tier) +
             labs(title = "Fig S14B — matched-null retention (% of outliers with a shift-matched neutral)",
                  x = NULL, y = "% matched") + theme_case(11))
save_supp(figS13, "FigS14_He_null_diagnostics.png", width = 9, height = 8)
figS13
```

![](Final_reproducible_analysis_files/figure-gfm/figS13-he-null-diag-1.png)<!-- -->

``` r
# Per-tier composite: 2D KDE (points coloured by stressor) with MAF0 (top) and |dAF|
# (right) marginal densities. CON dropped. Faceted by Tier (Core/Convergent/Private)
# instead of by stressor; points are coloured by stressor (CA/SE/CASE).
arch <- master[tier != "Background"]; arch[, tier := droplevels(tier)]
stressor_cols <- unlist(treat_cols[c("CA","SE","CASE")])
arch_panel_tier <- function(tr, leg = FALSE) {
  d <- arch[tier == tr]
  center <- ggplot(d, aes(MAF0, abs_dAF, colour = treatment, group = treatment)) +
    #geom_point(alpha = 0.15, size = 1, stroke = 0.5) +
    geom_point(aes(shape=treatment),alpha = 0.1, size = 0.25, stroke = 1) +
    geom_point(aes(shape=treatment),alpha = 0.2, size = 0.33, stroke = 1, col="gray", fill="gray") +
    geom_density_2d() +
    # tier label as an annotation anchored to the scatter's top-left (no plot title,
    # which would otherwise sit awkwardly below the top marginal)
    annotate("text", x = -Inf, y = Inf, label = tr, hjust = -0.25, vjust = 1.4,
             fontface = "bold", size = 4.2) +
    scale_colour_manual(values = stressor_cols, drop = FALSE, guide = "none") +
    scale_y_continuous(limits=c(0,0.4))+
    scale_shape_manual(values=c(23,24,25)) +
    labs(colour = "Treatment", shape = "Treatment")+
    labs(x = "Starting MAF", y = "|ΔAF|", title = NULL) +
    # enlarge the legend keys only (plotted points keep their mapped alpha)
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme_case(12) + theme(legend.position = if (leg) "bottom" else "none")
  # bounds = ... clips the kernel density estimate to the valid range of each quantity
  # (MAF in [0, 0.5]; |dAF| in [0, 1]) so the marginal densities don't bleed past 0 or
  # imply density beyond what's possible.
  # MAF0 is the starting allele frequency of each SNP, measured once before
  # treatment -- it is a property of the SNP/block, not of the stressor. CA/SE/CASE
  # rows for a given SNP+block share (essentially) the same MAF0, so splitting this
  # marginal by treatment just draws three near-identical curves. Use the
  # deduplicated per-SNP/block MAF0 values and a single neutral density instead.
  d_maf0 <- unique(d, by = c("SNP", "block"))
  top <- ggplot(d_maf0, aes(x = MAF0)) +
    geom_density(aes(), colour = cbPalette[2], fill = cbPalette[2], alpha = 0.4,
                  bounds = c(0, 0.5), na.rm = TRUE) +
    theme_void()
  # Ridgeline marginal for |dAF|: each stressor gets its own row (y = treatment), so
  # curves never overlap or hide one another regardless of relative peak height or
  # shape. CASE is plotted on the bottom row (closest to `center`'s axis), CA above
  # it, SE on top.
  right <- ggplot(d, aes(x = abs_dAF, y = treatment, fill = treatment, colour = treatment)) +
    geom_density_ridges(alpha = 0.6, scale = 1.2, 
                          rel_min_height = 0.01, na.rm = TRUE) +
    scale_fill_manual(values = stressor_cols, drop = FALSE, guide = "none") +
    scale_colour_manual(values = stressor_cols, drop = FALSE, guide = "none") +
    scale_y_discrete(expand = expansion(add = c(0.2, 1.2))) +
    coord_flip() + theme_void()
  patch <- top + plot_spacer() + center + right +
    plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4)) 
  patch[[3]] <- patch[[3]] + plot_annotation(tag_levels = 'A')
  patch
}
fig4_arch <- (arch_panel_tier("Core") | arch_panel_tier("Convergent", leg = TRUE) | arch_panel_tier("Private"))
# Stack the tier-architecture row (top) over the diversity-at-selected-loci panel (bottom).
# wrap_elements keeps the nested architecture row as one element so the bottom panel does not get
# pulled into its 2x2 layout. Heights/tags may want a visual tweak after knitting.
fig4 <- patchwork::wrap_elements(plot = fig4_arch) / (fig4_he_panel)+
  plot_layout(heights = c(2, 1.5)) + plot_annotation(tag_levels = 'A')
save_case(fig4, "Fig5_architecture_density2d.png", width = 14, height = 9)
```

    ## Picking joint bandwidth of 0.00804

    ## Picking joint bandwidth of 0.0107

    ## Picking joint bandwidth of 0.0116

``` r
fig4
```

    ## Picking joint bandwidth of 0.00804

    ## Picking joint bandwidth of 0.0107

    ## Picking joint bandwidth of 0.0116

![](Final_reproducible_analysis_files/figure-gfm/fig4-density2d-1.png)<!-- -->

## Fig S15 — Top enriched GO terms per tier (topGO)

``` r
# Top enriched GO terms per tier (Core/Convergent/Private), all stressors pooled,
# lollipop style. Stacked vertically (one tier per row) so term labels have room.
# Reuses go_dir2 / go_cat_cols from Fig S9. Point size = proportion of significant genes
# annotated to that term; colour = GO category (BP/CC/MF). Not split by stressor.
go_tier_files <- c(Core = "Core_ALL_GO_en_sig_gene.csv",
                   Convergent = "Conv_ALL_GO_en_sig_gene.csv",
                   `Private` = "Priv_ALL_GO_en_sig_gene.csv")
plot_go_lollipop <- function(file, label, top_n = 10, leg = FALSE) {
  d <- as.data.table(fread(file.path(go_dir2, file)))
  d <- d[, .(GO = GeneOntologyIDs, Term, Fisher, ontology, prop.sig.genes)]
  d[, Fisher := suppressWarnings(as.numeric(Fisher))]
  d <- d[is.finite(Fisher) & Fisher > 0]
  d <- unique(d, by = "Term")
  d[, cat := data.table::fcase(grepl("Biolog", ontology), "GO:BP",
        grepl("Cellul", ontology), "GO:CC", grepl("Molec", ontology), "GO:MF", default = "other")]
  d[, neglogp := -log10(Fisher)]
  setorder(d, -neglogp)
  d <- head(d, top_n)
  setorder(d, neglogp)
  d[, Term := stringr::str_wrap(Term, width = 35)]
  d[, Term := factor(Term, levels = unique(Term))]
  ggplot(d, aes(x = neglogp, y = Term)) +
    geom_segment(aes(x = 0, xend = neglogp, yend = Term), colour = "grey70") +
    geom_point(aes(size = prop.sig.genes, colour = cat)) +
    scale_colour_manual(values = go_cat_cols, name = "GO category", drop = FALSE) +
    scale_size_continuous(range = c(1.5, 5), limits = c(0, 1), name = "Prop. sig. genes") +
    labs(x = expression(-log[10](Fisher~italic(p))), y = NULL, title = label) +
    guides(colour = guide_legend(nrow = 4), size = guide_legend(nrow = 3)) +
    theme_case(11) + theme(axis.text.y = element_text(size = 8, angle = 15),
                           plot.title = element_text(size = 12, face = "bold"),
                           legend.position = if (leg) "right" else "none", legend.direction = if (leg) "vertical" else "none", theme(legend.box = "horizontal"))
}


fig5 <- (free(plot_go_lollipop(go_tier_files[["Core"]], "Core"), type= "label") /
  free(plot_go_lollipop(go_tier_files[["Convergent"]], "Convergent") , type = "label")/
  free(plot_go_lollipop(go_tier_files[["Private"]], "Private", leg = TRUE), type="label")) / guide_area() + plot_layout(guides = "collect", heights = c(3,3,3,1))
fig5 <- fig5 + plot_annotation(title = "Fig S15 — top enriched GO terms by tier (all stressors pooled)") & theme(legend.box = "horizontal", legend.box.margin = margin(0, 0, -20, -100), legend.key.spacing.x = unit(3.0, "cm"))
save_supp(fig5, "FigS15_GO_by_tier.png", width = 7, height = 11)
fig5
```

![](Final_reproducible_analysis_files/figure-gfm/fig5-go-lollipop-1.png)<!-- -->

## Fig S16 — Tier architecture by spawn (3x3)

``` r
# Same composite cell (2D scatter + MAF0/|dAF| marginals) as fig4, but faceted into a
# 3x3 grid: rows = spawn (B10/B11/B12), cols = treatment (CA/SE/CASE). Density marginals
# kept on every cell per request (no outer-row/col-only trimming).
arch_panel_spawn <- function(trt, blk, leg = FALSE) {
  d <- arch[treatment == trt & block == blk]
  center <- ggplot(d, aes(MAF0, abs_dAF, colour = tier, group = tier)) +
    geom_point(aes(alpha = tier), size = 1, stroke = 0.5) +
    annotate("text", x = -Inf, y = Inf, label = paste0(blk, " — ", trt), hjust = -0.15, vjust = 1.4,
             fontface = "bold", size = 4.0) +
    scale_colour_manual(values = tier_cols, drop = FALSE, name = "Tier") +
    scale_alpha_manual(values = c(0.1, 0.2, 0.2), guide = "none") +
    scale_fill_manual(values = tier_cols, drop = FALSE, name = "Tier") +
    scale_y_continuous(limits = c(0, 0.4)) +
    labs(x = "Starting MAF", y = "|ΔAF|", title = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)),
           fill   = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme_case(11) + theme(legend.position = if (leg) "bottom" else "none")
  top <- ggplot(d, aes(MAF0, colour = tier, fill = tier)) +
    geom_density(alpha = 0.3, linewidth = 0.4, bounds = c(0, 0.5)) +
    scale_colour_manual(values = tier_cols, guide = "none") +
    scale_fill_manual(values = tier_cols, guide = "none") + theme_void()
  right <- ggplot(d, aes(abs_dAF, colour = tier, fill = tier)) +
    geom_density(alpha = 0.3, linewidth = 0.4, bounds = c(0, 1)) +
    scale_colour_manual(values = tier_cols, guide = "none") +
    scale_fill_manual(values = tier_cols, guide = "none") + coord_flip() + theme_void()
  top + plot_spacer() + center + right +
    plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))
}

spawn_rows <- c("B10", "B11", "B12")
treat_cols_fig4 <- c("CA", "SE", "CASE")

# single shared legend: turn it on for the B12/SE cell (bottom row, middle column)
fig4_grid <- lapply(spawn_rows, function(blk) {
  cells <- lapply(treat_cols_fig4, function(trt) {
    arch_panel_spawn(trt, blk, leg = (blk == "B12" && trt == "SE"))
  })
  wrap_plots(cells, nrow = 1)
})

fig4_byspawn <- wrap_plots(fig4_grid, ncol = 1) +
  plot_annotation(title = "Fig S16 — tier architecture by spawn: MAF × |ΔAF| (2D density + marginals)")
save_supp(fig4_byspawn, "FigS16_architecture_byspawn.png", width = 14, height = 14.5)
fig4_byspawn
```

![](Final_reproducible_analysis_files/figure-gfm/fig4-density2d-byspawn-1.png)<!-- -->

# Supplementary tables

``` r
# Table S1 — locus counts by tier x treatment
tier_count <- function(df, tier) data.frame(Tier = tier,
  CA   = sum(asL(df$Sig.CA),   na.rm = TRUE),
  SE   = sum(asL(df$Sig.SE),   na.rm = TRUE),
  CASE = sum(asL(df$Sig.CASE), na.rm = TRUE),
  Total = nrow(df))
tabS1 <- rbind(tier_count(core.sig, "Core"),
               tier_count(convergent.sig, "Convergent"),
               tier_count(private.sig, "Private"))
write.csv(tabS1, file.path(tab_dir, "TableS1_tier_locus_counts.csv"), row.names = FALSE)
print(tabS1)
```

    ##         Tier  CA   SE CASE Total
    ## 1       Core 662 1661 2989  5312
    ## 2 Convergent  62  166  812  1040
    ## 3    Private 233  249  639   972

``` r
# Table S2 — full significant-loci table (tier + treatment flags)
write.csv(as.data.frame(grouped_total.sig),
          file.path(tab_dir, "TableS2_significant_loci.csv"), row.names = FALSE)

# Table S3 — GO enrichment (per-stressor topGO results, combined)
go_all <- rbindlist(lapply(c("ALL","CA","SE","CASE"), function(tr) {
  d <- fread(file.path("./GO", paste0(tr, "_GO_en_sig_gene.csv"))); d$stressor <- tr; d
}), fill = TRUE)
fwrite(go_all, file.path(tab_dir, "TableS3_GO_enrichment.csv"))

# Table S4 — candidate gene lists per stressor
loc_src <- c(CA = "Sig.loci.1.CA.LOC", SE = "Sig.loci.1.SE.LOC", CASE = "Sig.loci.1.CASE.LOC")
for (tr in names(loc_src))
  file.copy(loc_src[[tr]], file.path(tab_dir, paste0("TableS4_candidate_genes_", tr, ".txt")), overwrite = TRUE)
```

## Table S3 — GO enrichment terms by stressor (one row per term)

The Fig S8-S10 dot plots are dense; the tables below give the same topGO
results (deduplicated to one row per GO term) in a readable form, sorted
by Fisher p-value.

``` r
go_term_table <- function(df, stressor_label) {
  d <- df[stressor == stressor_label,
          .(Term, ontology, Annotated, Significant, Expected, Fisher, prop.sig.genes)]
  d <- unique(d, by = "Term")
  d[, Fisher := signif(Fisher, 3)]
  d[, prop.sig.genes := round(prop.sig.genes, 2)]
  setorder(d, Fisher)
  d
}

go_table_labels <- c(ALL = "All Outlier Loci", CA = "CA Stressor",
                      SE = "SE Stressor", CASE = "CASE Stressor (CA+SE)")
go_table_letters <- c(ALL = "a", CA = "b", SE = "c", CASE = "d")

if (nrow(go_all)) {
  for (tr in c("ALL", "CA", "SE", "CASE")) {
    d <- go_term_table(go_all, tr)
    if (nrow(d) > 0) {
      cat(paste0("\n### Table S3", go_table_letters[[tr]], " — ",
                 go_table_labels[[tr]], "\n\n"))
      print(kable(d, row.names = FALSE,
                   caption = paste0("GO enrichment terms — ", go_table_labels[[tr]])) %>%
              kable_styling(full_width = FALSE, font_size = 9, latex_options = "scale_down") %>%
              scroll_box(height = "400px"))
    }
  }
}
```

### Table S3a — All Outlier Loci

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table" style="font-size: 9px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">
GO enrichment terms — All Outlier Loci
</caption>
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
Term
</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
ontology
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Annotated
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Significant
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Expected
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Fisher
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
prop.sig.genes
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
translation repressor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
233
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
12.96
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
250
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
16.21
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin protein ligase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
601
</td>
<td style="text-align:right;">
77
</td>
<td style="text-align:right;">
33.42
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
proteasome-mediated ubiquitin-dependent protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
655
</td>
<td style="text-align:right;">
87
</td>
<td style="text-align:right;">
42.47
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
protein polyubiquitination
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
72
</td>
<td style="text-align:right;">
35.28
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
structural constituent of ribosome
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
174
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
9.68
</td>
<td style="text-align:right;">
1.00e-07
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
unfolded protein binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
3.67
</td>
<td style="text-align:right;">
5.00e-07
</td>
<td style="text-align:right;">
0.24
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosolic small ribosomal subunit
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:right;">
1.20e-06
</td>
<td style="text-align:right;">
0.34
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
160
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
8.90
</td>
<td style="text-align:right;">
1.80e-06
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosolic large ribosomal subunit
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3.21
</td>
<td style="text-align:right;">
1.90e-06
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
actin filament binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
6.51
</td>
<td style="text-align:right;">
6.60e-06
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
protein folding
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
130
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
8.43
</td>
<td style="text-align:right;">
1.10e-05
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
833
</td>
<td style="text-align:right;">
102
</td>
<td style="text-align:right;">
46.32
</td>
<td style="text-align:right;">
3.10e-05
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
motile cilium
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
1.70
</td>
<td style="text-align:right;">
8.30e-05
</td>
<td style="text-align:right;">
0.38
</td>
</tr>
<tr>
<td style="text-align:left;">
GTPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
294
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
16.35
</td>
<td style="text-align:right;">
1.00e-04
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1286
</td>
<td style="text-align:right;">
102
</td>
<td style="text-align:right;">
71.51
</td>
<td style="text-align:right;">
1.60e-04
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoskeleton
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
705
</td>
<td style="text-align:right;">
105
</td>
<td style="text-align:right;">
46.13
</td>
<td style="text-align:right;">
1.70e-04
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
axonemal dynein complex assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1.49
</td>
<td style="text-align:right;">
2.30e-04
</td>
<td style="text-align:right;">
0.35
</td>
</tr>
<tr>
<td style="text-align:left;">
chaperonin-containing T-complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
2.30e-04
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule cytoskeleton
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
453
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
29.64
</td>
<td style="text-align:right;">
2.50e-04
</td>
<td style="text-align:right;">
0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
translation initiation factor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
2.67
</td>
<td style="text-align:right;">
2.60e-04
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
cell projection
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
625
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
40.89
</td>
<td style="text-align:right;">
5.10e-04
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
359
</td>
<td style="text-align:right;">
43
</td>
<td style="text-align:right;">
19.96
</td>
<td style="text-align:right;">
5.40e-04
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
158
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
8.79
</td>
<td style="text-align:right;">
6.60e-04
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
FMN binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
6.90e-04
</td>
<td style="text-align:right;">
0.36
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
165
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
10.70
</td>
<td style="text-align:right;">
7.60e-04
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein light intermediate chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2.06
</td>
<td style="text-align:right;">
8.10e-04
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
cortical actin cytoskeleton
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1.83
</td>
<td style="text-align:right;">
9.90e-04
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium-dependent cell motility
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:right;">
1.02e-03
</td>
<td style="text-align:right;">
0.45
</td>
</tr>
<tr>
<td style="text-align:left;">
calcium import into the mitochondrion
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
1.03e-03
</td>
<td style="text-align:right;">
0.75
</td>
</tr>
<tr>
<td style="text-align:left;">
uniplex complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
1.06e-03
</td>
<td style="text-align:right;">
0.75
</td>
</tr>
<tr>
<td style="text-align:left;">
proton-transporting ATPase activity, rotational mechanism
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1.28
</td>
<td style="text-align:right;">
1.30e-03
</td>
<td style="text-align:right;">
0.26
</td>
</tr>
<tr>
<td style="text-align:left;">
translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
636
</td>
<td style="text-align:right;">
93
</td>
<td style="text-align:right;">
41.24
</td>
<td style="text-align:right;">
1.96e-03
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium movement involved in cell motility
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
1.96e-03
</td>
<td style="text-align:right;">
0.38
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular sodium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
2.46e-03
</td>
<td style="text-align:right;">
0.60
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular potassium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
2.46e-03
</td>
<td style="text-align:right;">
0.60
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium ion export across plasma membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
2.46e-03
</td>
<td style="text-align:right;">
0.60
</td>
</tr>
<tr>
<td style="text-align:left;">
motor neuron axon guidance
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
2.46e-03
</td>
<td style="text-align:right;">
0.60
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial calcium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
2.46e-03
</td>
<td style="text-align:right;">
0.60
</td>
</tr>
<tr>
<td style="text-align:left;">
amino acid transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
3.57
</td>
<td style="text-align:right;">
2.69e-03
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
centrosome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
136
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
8.90
</td>
<td style="text-align:right;">
2.90e-03
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
poly(A) binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.02e-03
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoplasmic translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
3.37
</td>
<td style="text-align:right;">
3.04e-03
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA polymerase II activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
3.09e-03
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inorganic phosphate transmembrane transporter activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
3.09e-03
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
calcium ion binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
846
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
47.05
</td>
<td style="text-align:right;">
3.40e-03
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium movement
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
3.31
</td>
<td style="text-align:right;">
3.45e-03
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
nucleus
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3481
</td>
<td style="text-align:right;">
273
</td>
<td style="text-align:right;">
227.76
</td>
<td style="text-align:right;">
3.55e-03
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
5.23
</td>
<td style="text-align:right;">
3.94e-03
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
protein stabilization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
4.01e-03
</td>
<td style="text-align:right;">
0.36
</td>
</tr>
<tr>
<td style="text-align:left;">
cortical cytoskeleton organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.10
</td>
<td style="text-align:right;">
4.19e-03
</td>
<td style="text-align:right;">
0.24
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular response to cAMP
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
4.20e-03
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of synaptic vesicle exocytosis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
4.20e-03
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
protein N-linked glycosylation via asparagine
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
4.69e-03
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
retrograde protein transport, ER to cytosol
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
4.69e-03
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of canonical Wnt signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
5.71e-03
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
ciliary basal body
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
3.53
</td>
<td style="text-align:right;">
6.05e-03
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP-dependent microtubule motor activity, minus-end-directed
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1.72
</td>
<td style="text-align:right;">
6.48e-03
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
zinc ion binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2113
</td>
<td style="text-align:right;">
143
</td>
<td style="text-align:right;">
117.50
</td>
<td style="text-align:right;">
7.03e-03
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
cell motility
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
122
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
7.91
</td>
<td style="text-align:right;">
7.55e-03
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule-based movement
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
204
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
13.23
</td>
<td style="text-align:right;">
7.61e-03
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosol
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
810
</td>
<td style="text-align:right;">
91
</td>
<td style="text-align:right;">
53.00
</td>
<td style="text-align:right;">
7.83e-03
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
poly(U) RNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
8.93e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
polynucleotide 3’-phosphatase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
8.93e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP-dependent polydeoxyribonucleotide 5’-hydroxyl-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
8.93e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
transmembrane-ephrin receptor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
8.93e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
myosin complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3.79
</td>
<td style="text-align:right;">
9.14e-03
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
motor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
6.12
</td>
<td style="text-align:right;">
9.94e-03
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
chromatin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
141
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
7.84
</td>
<td style="text-align:right;">
9.97e-03
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
protein tyrosine kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
90
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
5.00
</td>
<td style="text-align:right;">
1.00e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
axon guidance
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
86
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
5.58
</td>
<td style="text-align:right;">
1.05e-02
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
139
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
9.09
</td>
<td style="text-align:right;">
1.05e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
histone binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3.95
</td>
<td style="text-align:right;">
1.06e-02
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
protein tag
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
1.12e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin-dependent ERAD pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3.11
</td>
<td style="text-align:right;">
1.13e-02
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
SRP-dependent cotranslational protein targeting to membrane,
translocation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
1.21e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
Rap protein signal transduction
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
1.21e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
translation reinitiation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
1.21e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium:potassium-exchanging ATPase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.23e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
CAF-1 complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.23e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
acrosomal vesicle
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.23e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
nascent polypeptide-associated complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.23e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
nuclear outer membrane
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.23e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
calcium-dependent cysteine-type endopeptidase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.45
</td>
<td style="text-align:right;">
1.30e-02
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
adenylate kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
1.53e-02
</td>
<td style="text-align:right;">
0.30
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA 3’-UTR binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.56e-02
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
potassium ion import across plasma membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
1.71e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
2 iron, 2 sulfur cluster binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.56
</td>
<td style="text-align:right;">
1.78e-02
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
actin cytoskeleton organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
204
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
13.23
</td>
<td style="text-align:right;">
1.88e-02
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
GTP binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
407
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
22.63
</td>
<td style="text-align:right;">
1.96e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
syntaxin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2.17
</td>
<td style="text-align:right;">
1.97e-02
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
cell-cell adhesion mediator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:right;">
2.05e-02
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
recycling endosome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.31
</td>
<td style="text-align:right;">
2.18e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of mRNA processing
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2.08
</td>
<td style="text-align:right;">
2.28e-02
</td>
<td style="text-align:right;">
0.28
</td>
</tr>
<tr>
<td style="text-align:left;">
nucleosome organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
4.09
</td>
<td style="text-align:right;">
2.28e-02
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein intermediate chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3.50
</td>
<td style="text-align:right;">
2.29e-02
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
small GTPase binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
4.17
</td>
<td style="text-align:right;">
2.29e-02
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of mRNA splicing, via spliceosome
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1.69
</td>
<td style="text-align:right;">
2.30e-02
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
actin filament depolymerization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
2.30e-02
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA splice site selection
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
2.31e-02
</td>
<td style="text-align:right;">
0.31
</td>
</tr>
<tr>
<td style="text-align:left;">
endosome organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
2.31e-02
</td>
<td style="text-align:right;">
0.30
</td>
</tr>
<tr>
<td style="text-align:left;">
protein localization to ciliary transition zone
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
2.31e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
formation of cytoplasmic translation initiation complex
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
2.31e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoskeleton-dependent intracellular transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
4.60
</td>
<td style="text-align:right;">
2.32e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
meiotic nuclear membrane microtubule tethering complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
2.35e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
endoplasmic reticulum lumen
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
2.37e-02
</td>
<td style="text-align:right;">
0.30
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule motor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
76
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4.23
</td>
<td style="text-align:right;">
2.44e-02
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
GDP binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
2.59e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
nuclear receptor transcription coactivator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
2.59e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
small ribosomal subunit rRNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
2.76e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
glucosamine 6-phosphate N-acetyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
2.76e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA 5’-UTR binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
2.76e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
protein-glutamine gamma-glutamyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
2.76e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin-dependent protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
750
</td>
<td style="text-align:right;">
100
</td>
<td style="text-align:right;">
48.63
</td>
<td style="text-align:right;">
2.94e-02
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule organizing center
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
200
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
13.09
</td>
<td style="text-align:right;">
2.97e-02
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoplasmic microtubule
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.24
</td>
<td style="text-align:right;">
3.10e-02
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
core mediator complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
3.10e-02
</td>
<td style="text-align:right;">
0.27
</td>
</tr>
<tr>
<td style="text-align:left;">
Set1C/COMPASS complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
3.10e-02
</td>
<td style="text-align:right;">
0.27
</td>
</tr>
<tr>
<td style="text-align:left;">
endoplasmic reticulum exit site
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
3.10e-02
</td>
<td style="text-align:right;">
0.27
</td>
</tr>
<tr>
<td style="text-align:left;">
sleep
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
3.12e-02
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of synaptic transmission, cholinergic
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
3.12e-02
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of voltage-gated potassium channel activity
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
3.12e-02
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
DNA replication-dependent nucleosome assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.68e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
Golgi to endosome transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.68e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of store-operated calcium entry
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.68e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of cilium movement
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
3.68e-02
</td>
<td style="text-align:right;">
0.43
</td>
</tr>
<tr>
<td style="text-align:left;">
synaptic vesicle priming
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.68e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
synaptic vesicle fusion to presynaptic active zone membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.68e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
peptide cross-linking
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.68e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
axonemal dynein complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
3.73e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial proton-transporting ATP synthase complex, coupling factor
F(o)
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.74e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
oligosaccharyltransferase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.74e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
proton-transporting V-type ATPase, V0 domain
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.74e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
vacuolar proton-transporting V-type ATPase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.74e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
nuclear envelope
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
91
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
5.95
</td>
<td style="text-align:right;">
3.80e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
Rho protein signal transduction
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.17
</td>
<td style="text-align:right;">
3.84e-02
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2.98
</td>
<td style="text-align:right;">
3.84e-02
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
phosphatidylinositol phospholipase C activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.99e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
thioredoxin peroxidase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.99e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
steroid 17-alpha-monooxygenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.99e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
17-alpha-hydroxyprogesterone aldolase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
3.99e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
nucleosome assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
59
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3.83
</td>
<td style="text-align:right;">
4.09e-02
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
ribosome binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.33
</td>
<td style="text-align:right;">
4.15e-02
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin protein ligase binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.95
</td>
<td style="text-align:right;">
4.28e-02
</td>
<td style="text-align:right;">
0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
transcription coregulator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
181
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
10.07
</td>
<td style="text-align:right;">
4.40e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
nucleolus
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
181
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
11.84
</td>
<td style="text-align:right;">
4.43e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein heavy chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
4.73e-02
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
integrin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.39
</td>
<td style="text-align:right;">
4.74e-02
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
ribosomal small subunit assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
4.77e-02
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
U2-type prespliceosome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
4.88e-02
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
Z disc
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
4.88e-02
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
kinesin complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
4.88e-02
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
enzyme regulator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
566
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
31.47
</td>
<td style="text-align:right;">
4.92e-02
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
</tbody>
</table>

</div>

### Table S3b — CA Stressor

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table" style="font-size: 9px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">
GO enrichment terms — CA Stressor
</caption>
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
Term
</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
ontology
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Annotated
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Significant
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Expected
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Fisher
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
prop.sig.genes
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
translation repressor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
233
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
2.94
</td>
<td style="text-align:right;">
0.000000
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
250
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
3.47
</td>
<td style="text-align:right;">
0.000000
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
protein polyubiquitination
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
7.54
</td>
<td style="text-align:right;">
0.000000
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
proteasome-mediated ubiquitin-dependent protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
655
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
9.08
</td>
<td style="text-align:right;">
0.000000
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin protein ligase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
601
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
7.57
</td>
<td style="text-align:right;">
0.000000
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosolic small ribosomal subunit
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.000088
</td>
<td style="text-align:right;">
0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
acrosomal vesicle
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.000510
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA splice site selection
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.000680
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoplasmic stress granule
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.000740
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA 3’-UTR binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.001400
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular sodium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.001860
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular potassium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.001860
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium ion export across plasma membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.001860
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
zinc ion binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2113
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
26.63
</td>
<td style="text-align:right;">
0.001900
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
proton-transporting ATPase activity, rotational mechanism
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.002900
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
636
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
8.82
</td>
<td style="text-align:right;">
0.002910
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
215
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2.83
</td>
<td style="text-align:right;">
0.004870
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
ribonuclease MRP complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.005850
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
structural constituent of ribosome
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
174
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2.19
</td>
<td style="text-align:right;">
0.006600
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP-dependent microtubule motor activity, minus-end-directed
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.006800
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
proton transmembrane transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
89
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
0.007900
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein light intermediate chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.011200
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
tRNA 5’-leader removal
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.011500
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
6-phosphofructo-2-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
L-galactose dehydrogenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ornithine decarboxylase inhibitor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ornithine decarboxylase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
calcium sensitive guanylate cyclase activator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
sphingolipid delta-4 desaturase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ephrin receptor binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol tetrakisphosphate 1-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol-1,3,4-trisphosphate 6-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol-1,3,4-trisphosphate 5-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
cohesin ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.012600
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial pyruvate dehydrogenase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013200
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA cleavage stimulating factor complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013200
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
keratin filament
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013200
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ribosomal small subunit assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.013500
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
mesoderm formation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013900
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
heart morphogenesis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013900
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
adiponectin-activated signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013900
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
sensory neuron axon guidance
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013900
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
putrescine biosynthetic process from ornithine
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013900
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
intracellular mRNA localization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013900
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ephrin receptor signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013900
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol trisphosphate metabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.013900
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
potassium ion import across plasma membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.020200
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of alternative mRNA splicing, via spliceosome
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.020200
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
motile cilium
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.022800
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
GDP-L-fucose synthase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
pyruvate dehydrogenase (acetyl-transferring) activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
protein-glycine ligase activity, initiating
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
glucose-6-phosphate dehydrogenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
3-hydroxyisobutyryl-CoA hydrolase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
Atg8-specific protease activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
tachykinin receptor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
methionine adenosyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium:potassium-exchanging ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
FFAT motif binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.025000
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
centriole
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.025200
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
protein sumoylation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.025300
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
integrin complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.025500
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
VCP-NPL4-UFD1 AAA ATPase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.026200
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
protein polyglycylation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.027500
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular response to cAMP
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.027500
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of synaptic vesicle exocytosis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.027500
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
oligopeptide transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.027500
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
mitotic spindle disassembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.027500
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
anterior/posterior axis specification
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.027500
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
endoplasmic reticulum membrane organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.027500
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
cell adhesion mediated by integrin
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.028000
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium movement
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.033600
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin-like protein ligase binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.037200
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
myosin II binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
tRNA (guanine-N7-)-methyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
fructose-2,6-bisphosphate 2-phosphatase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
translation activator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
glucose 6-phosphate:inorganic phosphate antiporter activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
dolichyl-phosphate-mannose-protein mannosyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
L-malate dehydrogenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
ribonuclease P activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
poly(U) RNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
purine-nucleoside phosphorylase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.037300
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium:potassium-exchanging ATPase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.039000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
CAF-1 complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.039000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
fibrillar center
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.039000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
integrin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.039200
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
integrin-mediated signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.040000
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
fructose 2,6-bisphosphate metabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.041000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
Rap protein signal transduction
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.041000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
pentose-phosphate shunt, oxidative branch
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.041000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
bile acid biosynthetic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.041000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
S-adenosylmethionine biosynthetic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.041000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
glucose-6-phosphate transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.041000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of protein ubiquitination
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.041000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
endoplasmic reticulum-plasma membrane tethering
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.041000
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.043200
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein intermediate chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.045200
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
cAMP binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.049500
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
protein kinase C activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.049500
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
SUMO activating enzyme activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.049500
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
ribonuclease P RNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.049500
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
</tbody>
</table>

</div>

### Table S3c — SE Stressor

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table" style="font-size: 9px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">
GO enrichment terms — SE Stressor
</caption>
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
Term
</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
ontology
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Annotated
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Significant
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Expected
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Fisher
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
prop.sig.genes
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
translation repressor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
233
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
5.71
</td>
<td style="text-align:right;">
1.00e-07
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
250
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
6.91
</td>
<td style="text-align:right;">
1.00e-07
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin protein ligase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
601
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
14.72
</td>
<td style="text-align:right;">
3.00e-07
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
proteasome-mediated ubiquitin-dependent protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
655
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
18.11
</td>
<td style="text-align:right;">
3.00e-07
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
protein polyubiquitination
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
15.04
</td>
<td style="text-align:right;">
4.00e-07
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
unfolded protein binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
1.62
</td>
<td style="text-align:right;">
4.40e-06
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
structural constituent of ribosome
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
174
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
4.26
</td>
<td style="text-align:right;">
1.00e-04
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
chaperonin-containing T-complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
1.20e-04
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
actin filament binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
2.87
</td>
<td style="text-align:right;">
1.40e-04
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
proton-transporting ATPase activity, rotational mechanism
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
2.00e-04
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
833
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
20.40
</td>
<td style="text-align:right;">
2.60e-04
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosolic small ribosomal subunit
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
4.10e-04
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosolic large ribosomal subunit
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
1.40
</td>
<td style="text-align:right;">
4.40e-04
</td>
<td style="text-align:right;">
0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1286
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
31.49
</td>
<td style="text-align:right;">
4.70e-04
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
359
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
8.79
</td>
<td style="text-align:right;">
5.50e-04
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
inorganic phosphate transmembrane transporter activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
6.00e-04
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ciliary basal body
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1.54
</td>
<td style="text-align:right;">
6.30e-04
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
cell projection
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
625
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
17.80
</td>
<td style="text-align:right;">
7.00e-04
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
mitotic cell cycle
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
258
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
7.14
</td>
<td style="text-align:right;">
1.10e-03
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
motile cilium
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
1.43e-03
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
cortical actin cytoskeleton
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
1.76e-03
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
160
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3.92
</td>
<td style="text-align:right;">
1.96e-03
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
nascent polypeptide-associated complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
2.38e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
protein folding
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
130
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
3.60
</td>
<td style="text-align:right;">
2.50e-03
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
protein stabilization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
2.90e-03
</td>
<td style="text-align:right;">
0.27
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.27
</td>
<td style="text-align:right;">
3.80e-03
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoskeleton
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
705
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
20.08
</td>
<td style="text-align:right;">
3.98e-03
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
protein localization to ciliary transition zone
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
4.40e-03
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
proton transmembrane transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
89
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2.46
</td>
<td style="text-align:right;">
5.70e-03
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
small ribosomal subunit rRNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
5.70e-03
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
actin cytoskeleton organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
204
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
5.64
</td>
<td style="text-align:right;">
6.00e-03
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
phosphatidylinositol-mediated signaling
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
6.10e-03
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule cytoskeleton
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
453
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
12.90
</td>
<td style="text-align:right;">
6.58e-03
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular sodium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
7.20e-03
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular potassium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
7.20e-03
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium ion export across plasma membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
7.20e-03
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
proton-transporting V-type ATPase, V0 domain
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
7.64e-03
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
vacuolar proton-transporting V-type ATPase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
7.64e-03
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
zinc ion binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2113
</td>
<td style="text-align:right;">
69
</td>
<td style="text-align:right;">
51.74
</td>
<td style="text-align:right;">
7.92e-03
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
chromatin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
141
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3.45
</td>
<td style="text-align:right;">
8.23e-03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
phosphatidylinositol phospholipase C activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
8.41e-03
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of mRNA processing
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.05e-02
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
retrograde protein transport, ER to cytosol
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.06e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
motile cilium assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.06e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
histone binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1.74
</td>
<td style="text-align:right;">
1.09e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
box C/D snoRNP complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.12e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
aconitate hydratase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.16e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein light intermediate chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.24e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
intracellular signal transduction
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
688
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
19.03
</td>
<td style="text-align:right;">
1.34e-02
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
GTPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
294
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
7.20
</td>
<td style="text-align:right;">
1.40e-02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
enzyme regulator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
566
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
13.86
</td>
<td style="text-align:right;">
1.41e-02
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin-dependent protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
750
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
20.74
</td>
<td style="text-align:right;">
1.46e-02
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
actin filament severing
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
1.46e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
syntaxin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.48e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule plus-end binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.52e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of protein secretion
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
1.91e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoplasmic translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.44
</td>
<td style="text-align:right;">
2.20e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
165
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4.56
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
adenylate kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.36e-02
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
creatine kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol pentakisphosphate 2-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
saccharopine dehydrogenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ornithine decarboxylase inhibitor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ornithine decarboxylase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
copper-transporting ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium channel inhibitor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ephrin receptor binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
disordered domain specific binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol tetrakisphosphate 1-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol-1,3,4-trisphosphate 6-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol-1,3,4-trisphosphate 5-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
dihydrolipoyllysine-residue succinyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
choline-phosphate cytidylyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA polymerase III core binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
cohesin ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
very-long-chain-acyl-CoA dehydrogenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
glycerol-3-phosphate O-acyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
2.45e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of transcription by RNA polymerase II
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2.24
</td>
<td style="text-align:right;">
2.46e-02
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
axonemal dynein complex assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
2.47e-02
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
Arp2/3 protein complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
2.55e-02
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
actin cortical patch
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
2.55e-02
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
proteasome regulatory particle, base subcomplex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
2.55e-02
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of gene expression, epigenetic
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
2.76e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
snRNA modification
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
2.76e-02
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
phosphocreatine biosynthetic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquinone-6 biosynthetic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mesoderm formation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
heart morphogenesis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial phosphate ion transmembrane transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
putrescine biosynthetic process from ornithine
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
T cell receptor signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ephrin receptor signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol trisphosphate metabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA destabilization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of transcription by RNA polymerase III
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
glyoxylate catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
L-alanine catabolic process, by transamination
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
box H/ACA snoRNA 3’-end processing
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
membrane fusion
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1.55
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of cilium beat frequency involved in ciliary motility
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of dendritic spine morphogenesis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
postsynaptic actin cytoskeleton organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.77e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
MAP kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
2.84e-02
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
dense core granule
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial pyruvate dehydrogenase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Derlin-1-VIMP complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Derlin-1 retrotranslocation complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
amphisome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
npBAF complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Noc1p-Noc2p complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Noc2p-Noc3p complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
plasma membrane proton-transporting V-type ATPase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
RSF complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
2.85e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of mRNA splicing, via spliceosome
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
2.94e-02
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
removal of superoxide radicals
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
2.96e-02
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
amino acid transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1.52
</td>
<td style="text-align:right;">
2.97e-02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
158
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
3.87
</td>
<td style="text-align:right;">
3.10e-02
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
endoplasmic reticulum lumen
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
3.13e-02
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule organizing center
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
200
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
5.70
</td>
<td style="text-align:right;">
3.60e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoplasmic microtubule
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
3.74e-02
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule plus-end
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
3.76e-02
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP-dependent microtubule motor activity, minus-end-directed
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
3.95e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
peptide receptor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
292
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
7.15
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
centrosome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
136
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3.87
</td>
<td style="text-align:right;">
4.07e-02
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
barbed-end actin filament capping
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
4.19e-02
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of canonical Wnt signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
4.19e-02
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
pyruvate dehydrogenase (acetyl-transferring) activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
glucose-6-phosphate dehydrogenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
gamma-catenin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
ER retention sequence binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
Rho-dependent protein serine/threonine kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
Atg8-specific protease activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
methionine adenosyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium:potassium-exchanging ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
carnosine synthase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
cysteine synthase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
G protein-coupled glutamate receptor binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
phosphotyrosine residue binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
alanine-glyoxylate transaminase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA polymerase II sequence-specific DNA-binding transcription factor
binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
FFAT motif binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
pantothenate kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
pseudouridylate synthase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.84e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA splice site selection
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
4.87e-02
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
</tbody>
</table>

</div>

### Table S3d — CASE Stressor (CA+SE)

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:400px; ">

<table class="table" style="font-size: 9px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">
GO enrichment terms — CASE Stressor (CA+SE)
</caption>
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
Term
</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">
ontology
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Annotated
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Significant
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Expected
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
Fisher
</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">
prop.sig.genes
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
structural constituent of ribosome
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
174
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
7.06
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
160
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
6.49
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosolic large ribosomal subunit
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
2.42
</td>
<td style="text-align:right;">
1.00e-07
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
translation repressor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
233
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
9.45
</td>
<td style="text-align:right;">
1.90e-06
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosolic small ribosomal subunit
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
4.90e-06
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
250
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
11.99
</td>
<td style="text-align:right;">
8.80e-06
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin protein ligase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
601
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
24.38
</td>
<td style="text-align:right;">
1.10e-05
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
translation initiation factor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
1.95
</td>
<td style="text-align:right;">
1.80e-05
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
proteasome-mediated ubiquitin-dependent protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
655
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
31.40
</td>
<td style="text-align:right;">
1.00e-04
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
833
</td>
<td style="text-align:right;">
85
</td>
<td style="text-align:right;">
33.79
</td>
<td style="text-align:right;">
1.00e-04
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
protein folding
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
130
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
6.23
</td>
<td style="text-align:right;">
1.40e-04
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
unfolded protein binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
2.68
</td>
<td style="text-align:right;">
3.10e-04
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
calcium import into the mitochondrion
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
4.20e-04
</td>
<td style="text-align:right;">
0.75
</td>
</tr>
<tr>
<td style="text-align:left;">
uniplex complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
4.60e-04
</td>
<td style="text-align:right;">
0.75
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium movement involved in cell motility
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
7.00e-04
</td>
<td style="text-align:right;">
0.31
</td>
</tr>
<tr>
<td style="text-align:left;">
actin filament binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
4.75
</td>
<td style="text-align:right;">
9.30e-04
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial calcium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
1.02e-03
</td>
<td style="text-align:right;">
0.60
</td>
</tr>
<tr>
<td style="text-align:left;">
poly(A) binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
1.21e-03
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
359
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
14.56
</td>
<td style="text-align:right;">
1.28e-03
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
636
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
30.49
</td>
<td style="text-align:right;">
1.29e-03
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoskeleton
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
705
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
34.87
</td>
<td style="text-align:right;">
1.41e-03
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
motile cilium
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1.29
</td>
<td style="text-align:right;">
1.43e-03
</td>
<td style="text-align:right;">
0.31
</td>
</tr>
<tr>
<td style="text-align:left;">
RNA polymerase II activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
1.64e-03
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inorganic phosphate transmembrane transporter activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
1.64e-03
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
FMN binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.94e-03
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
proton-transporting ATPase activity, rotational mechanism
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
1.98e-03
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
protein polyubiquitination
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
26.08
</td>
<td style="text-align:right;">
2.06e-03
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
GTPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
294
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
11.92
</td>
<td style="text-align:right;">
2.10e-03
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1286
</td>
<td style="text-align:right;">
73
</td>
<td style="text-align:right;">
52.16
</td>
<td style="text-align:right;">
2.22e-03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
cell motility
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
122
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
5.85
</td>
<td style="text-align:right;">
3.23e-03
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein light intermediate chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1.50
</td>
<td style="text-align:right;">
3.46e-03
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
calcium-dependent cysteine-type endopeptidase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.05
</td>
<td style="text-align:right;">
3.51e-03
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium movement
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
2.44
</td>
<td style="text-align:right;">
4.79e-03
</td>
<td style="text-align:right;">
0.24
</td>
</tr>
<tr>
<td style="text-align:left;">
poly(U) RNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
4.80e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
polynucleotide 3’-phosphatase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
4.80e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP-dependent polydeoxyribonucleotide 5’-hydroxyl-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
4.80e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
2 iron, 2 sulfur cluster binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.14
</td>
<td style="text-align:right;">
4.90e-03
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
cell-cell adhesion mediator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.18
</td>
<td style="text-align:right;">
5.73e-03
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
adenylate kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
6.44e-03
</td>
<td style="text-align:right;">
0.30
</td>
</tr>
<tr>
<td style="text-align:left;">
SRP-dependent cotranslational protein targeting to membrane,
translocation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
6.67e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
translation reinitiation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
6.67e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
acrosomal vesicle
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
7.09e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
nuclear outer membrane
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
7.09e-03
</td>
<td style="text-align:right;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin-dependent ERAD pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2.30
</td>
<td style="text-align:right;">
7.50e-03
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
nucleus
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
3481
</td>
<td style="text-align:right;">
212
</td>
<td style="text-align:right;">
172.17
</td>
<td style="text-align:right;">
8.55e-03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
myosin complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2.87
</td>
<td style="text-align:right;">
8.91e-03
</td>
<td style="text-align:right;">
0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
cytoplasmic translation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2.49
</td>
<td style="text-align:right;">
9.88e-03
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule cytoskeleton
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
453
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
22.40
</td>
<td style="text-align:right;">
1.00e-02
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of mRNA splicing, via spliceosome
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.25
</td>
<td style="text-align:right;">
1.02e-02
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
small GTPase binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3.04
</td>
<td style="text-align:right;">
1.09e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
motor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
4.46
</td>
<td style="text-align:right;">
1.11e-02
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
nuclear receptor transcription coactivator activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
1.11e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
158
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
6.41
</td>
<td style="text-align:right;">
1.23e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium-dependent cell motility
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.28e-02
</td>
<td style="text-align:right;">
0.35
</td>
</tr>
<tr>
<td style="text-align:left;">
nucleosome organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
3.02
</td>
<td style="text-align:right;">
1.28e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquitin protein ligase binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.42
</td>
<td style="text-align:right;">
1.28e-02
</td>
<td style="text-align:right;">
0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
actin filament depolymerization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.29e-02
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
formation of cytoplasmic translation initiation complex
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
1.29e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
meiotic nuclear membrane microtubule tethering complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.37e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
core mediator complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
1.48e-02
</td>
<td style="text-align:right;">
0.27
</td>
</tr>
<tr>
<td style="text-align:left;">
endoplasmic reticulum exit site
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
1.48e-02
</td>
<td style="text-align:right;">
0.27
</td>
</tr>
<tr>
<td style="text-align:left;">
ribosome binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
1.49e-02
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
axon guidance
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
86
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
4.12
</td>
<td style="text-align:right;">
1.50e-02
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
small ribosomal subunit rRNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.51e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
glucosamine 6-phosphate N-acetyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.51e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA 5’-UTR binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.51e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
protein-glutamine gamma-glutamyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.51e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
nuclear envelope
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
91
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
4.50
</td>
<td style="text-align:right;">
1.53e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
NAD binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2.43
</td>
<td style="text-align:right;">
1.61e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
cytosol
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
810
</td>
<td style="text-align:right;">
72
</td>
<td style="text-align:right;">
40.06
</td>
<td style="text-align:right;">
1.68e-02
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
centrosome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
136
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
6.73
</td>
<td style="text-align:right;">
1.72e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of canonical Wnt signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.74e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3.96
</td>
<td style="text-align:right;">
1.92e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
syntaxin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.58
</td>
<td style="text-align:right;">
1.99e-02
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular sodium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular potassium ion homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
sodium ion export across plasma membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
DNA replication-dependent nucleosome assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
motor neuron axon guidance
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of store-operated calcium entry
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
synaptic vesicle priming
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
synaptic vesicle fusion to presynaptic active zone membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
peptide cross-linking
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.08e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein heavy chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
2.10e-02
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial proton-transporting ATP synthase complex, coupling factor
F(o)
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
proton-transporting V-type ATPase, V0 domain
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
thioredoxin peroxidase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
steroid 17-alpha-monooxygenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
17-alpha-hydroxyprogesterone aldolase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
2.21e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
cilium assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
165
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
7.91
</td>
<td style="text-align:right;">
2.31e-02
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
U2-type prespliceosome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
2.38e-02
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
protein tyrosine kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
90
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
3.65
</td>
<td style="text-align:right;">
2.50e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
axoneme assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1.77
</td>
<td style="text-align:right;">
2.65e-02
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
calcium ion binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
846
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
34.31
</td>
<td style="text-align:right;">
2.68e-02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
maturation of LSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S
rRNA, LSU-rRNA)
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
2.69e-02
</td>
<td style="text-align:right;">
0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
microtubule-based movement
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
204
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
9.78
</td>
<td style="text-align:right;">
2.84e-02
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
cyclin-dependent protein serine/threonine kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
2.96e-02
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
aconitate hydratase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
3.01e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of actin cytoskeleton organization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3.60
</td>
<td style="text-align:right;">
3.02e-02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
phosphate ion transmembrane transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
3.02e-02
</td>
<td style="text-align:right;">
0.43
</td>
</tr>
<tr>
<td style="text-align:left;">
retrograde protein transport, ER to cytosol
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
3.03e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
progesterone metabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
3.03e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
sperm flagellum
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
3.20e-02
</td>
<td style="text-align:right;">
0.38
</td>
</tr>
<tr>
<td style="text-align:left;">
box C/D snoRNP complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
3.21e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
proton transmembrane transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
89
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
4.27
</td>
<td style="text-align:right;">
3.26e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP-dependent microtubule motor activity, minus-end-directed
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.26
</td>
<td style="text-align:right;">
3.54e-02
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
chromatin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
141
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
5.72
</td>
<td style="text-align:right;">
3.64e-02
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
axoneme
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2.72
</td>
<td style="text-align:right;">
3.67e-02
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
modification-dependent protein catabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
757
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
36.29
</td>
<td style="text-align:right;">
3.78e-02
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
dendrite self-avoidance
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.29
</td>
<td style="text-align:right;">
3.82e-02
</td>
<td style="text-align:right;">
0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
potassium ion import across plasma membrane
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
3.85e-02
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
fatty-acyl-CoA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.91e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
kinesin binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.91e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
nucleoside diphosphate kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
3.91e-02
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
poly-pyrimidine tract binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
4.05e-02
</td>
<td style="text-align:right;">
0.75
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol-polyphosphate 5-phosphatase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
creatine kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol pentakisphosphate 2-kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
saccharopine dehydrogenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inositol-1,4-bisphosphate 1-phosphatase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
pyridoxamine-phosphate oxidase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ornithine decarboxylase inhibitor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
glutamine-fructose-6-phosphate transaminase (isomerizing) activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
xylose isomerase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Rho GDP-dissociation inhibitor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ornithine decarboxylase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
DNA-dependent protein kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
blue light photoreceptor activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ethanolamine-phosphate cytidylyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
copper-transporting ATPase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
MAP-kinase scaffold activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
JUN kinase binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
carnitine O-acetyltransferase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
dolichol kinase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
sphingolipid delta-4 desaturase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ephrin receptor binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
cysteine-type carboxypeptidase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
disordered domain specific binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
eukaryotic initiation factor eIF2 binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
siRNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
macrolide binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
N-acetylglucosaminylphosphatidylinositol deacetylase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mitogen-activated protein kinase kinase binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial promoter sequence-specific DNA binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Atg8 ligase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
3-keto sterol reductase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
UDP-glucose transmembrane transporter activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
very-long-chain-acyl-CoA dehydrogenase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
tripeptidyl-peptidase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquinol-cytochrome-c reductase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
GMP reductase activity
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
4.06e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
actin filament severing
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
4.11e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
formation of translation preinitiation complex
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
4.11e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of proteolysis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
124
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
5.94
</td>
<td style="text-align:right;">
4.13e-02
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
focal adhesion
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
4.17e-02
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
cell projection
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
625
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
30.91
</td>
<td style="text-align:right;">
4.19e-02
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
dynein intermediate chain binding
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2.56
</td>
<td style="text-align:right;">
4.20e-02
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
protein phosphatase type 2A complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
4.35e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial respiratory chain complex III
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
4.35e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
eukaryotic translation initiation factor 4F complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
4.35e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
endocytic recycling
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
4.52e-02
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
cation homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
180
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
8.63
</td>
<td style="text-align:right;">
4.75e-02
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
UDP-N-acetylglucosamine metabolic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
4.78e-02
</td>
<td style="text-align:right;">
0.33
</td>
</tr>
<tr>
<td style="text-align:left;">
membrane fusion
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2.68
</td>
<td style="text-align:right;">
4.78e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
phosphocreatine biosynthetic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ubiquinone-6 biosynthetic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of stress-activated MAPK cascade
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial phosphate ion transmembrane transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of ATPase activity
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
signal complex assembly
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
putrescine biosynthetic process from ornithine
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mitotic cohesin loading
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of chromosome condensation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
regulation of cohesin loading
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
intracellular mRNA localization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
phagosome maturation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
dolichyl monophosphate biosynthetic process
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
muscle cell development
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of sequestering of triglyceride
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
cargo loading into vesicle
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
0.29
</td>
</tr>
<tr>
<td style="text-align:left;">
ephrin receptor signaling pathway
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
protein localization involved in establishment of planar polarity
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
skeletal muscle fiber development
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mRNA destabilization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
B cell homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
T cell homeostasis
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of B cell differentiation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
negative regulation of mRNA polyadenylation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
cellular heat acclimation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
plus-end-directed vesicle transport along microtubule
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of autophagosome maturation
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial outer membrane permeabilization
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
positive regulation of tyrosine phosphorylation of STAT protein
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
maintenance of protein location in nucleus
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
UDP-galactose transmembrane transport
</td>
<td style="text-align:left;">
Biological Processes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.79e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
recycling endosome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.99
</td>
<td style="text-align:right;">
4.86e-02
</td>
<td style="text-align:right;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
protein tag
</td>
<td style="text-align:left;">
Molecular Functions
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
4.89e-02
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
ciliary basal body
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2.67
</td>
<td style="text-align:right;">
4.94e-02
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
dense core granule
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Scc2-Scc4 cohesin loading complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
inner dynein arm
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
GABA-A receptor complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
WICH complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Derlin-1-VIMP complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
Derlin-1 retrotranslocation complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
amphisome
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
npBAF complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
mitochondrial DNA-directed RNA polymerase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
plasma membrane proton-transporting V-type ATPase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
keratin filament
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
GMP reductase complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
RSF complex
</td>
<td style="text-align:left;">
Cellular Components
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
4.95e-02
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
</tbody>
</table>

</div>

## Table 1 — All data-driven GO clusters (main text)

``` r
knitr::kable(go_clusters_full,
             caption = "Table 1. Data-driven GO clusters (all terms; stress flag and per-stressor enrichment). Stress-relevant clusters are shown in Fig 4D.")
```

| cluster                                         | stress | GO           | Term                                                              |    Fisher |  CA |  SE | CASE |
|:------------------------------------------------|:-------|:-------------|:------------------------------------------------------------------|----------:|----:|----:|-----:|
| Ubiquitin-proteasome / translational repression | TRUE   | <GO:0030371> | translation repressor activity                                    | 0.0000000 |   1 |   1 |    1 |
| Ubiquitin-proteasome / translational repression | TRUE   | <GO:0017148> | negative regulation of translation                                | 0.0000000 |   1 |   1 |    1 |
| Ubiquitin-proteasome / translational repression | TRUE   | <GO:0000209> | protein polyubiquitination                                        | 0.0000000 |   1 |   1 |    1 |
| Ubiquitin-proteasome / translational repression | TRUE   | <GO:0043161> | proteasome-mediated ubiquitin-dependent protein catabolic process | 0.0000000 |   1 |   1 |    1 |
| Ubiquitin-proteasome / translational repression | TRUE   | <GO:0061630> | ubiquitin protein ligase activity                                 | 0.0000000 |   1 |   1 |    1 |
| Ubiquitin-proteasome / translational repression | TRUE   | <GO:0008270> | zinc ion binding                                                  | 0.0019000 |   1 |   1 |    0 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0003735> | structural constituent of ribosome                                | 0.0000000 |   1 |   1 |    1 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0022625> | cytosolic large ribosomal subunit                                 | 0.0000001 |   0 |   1 |    1 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0022627> | cytosolic small ribosomal subunit                                 | 0.0000049 |   1 |   1 |    1 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0003723> | RNA binding                                                       | 0.0001000 |   0 |   1 |    1 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0006412> | translation                                                       | 0.0012900 |   1 |   0 |    1 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0070181> | small ribosomal subunit rRNA binding                              | 0.0057000 |   0 |   1 |    1 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0002181> | cytoplasmic translation                                           | 0.0098800 |   0 |   1 |    1 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0000028> | ribosomal small subunit assembly                                  | 0.0134900 |   1 |   0 |    0 |
| Other: structural constituent of ribosome       | FALSE  | <GO:0048027> | mRNA 5’-UTR binding                                               | 0.0151400 |   0 |   0 |    1 |
| Translation initiation                          | TRUE   | <GO:0003729> | mRNA binding                                                      | 0.0000000 |   0 |   1 |    1 |
| Translation initiation                          | TRUE   | <GO:0003743> | translation initiation factor activity                            | 0.0000180 |   0 |   0 |    1 |
| Translation initiation                          | TRUE   | <GO:0002188> | translation reinitiation                                          | 0.0066700 |   0 |   0 |    1 |
| Translation initiation                          | TRUE   | <GO:0048024> | regulation of mRNA splicing, via spliceosome                      | 0.0101700 |   0 |   1 |    1 |
| Translation initiation                          | TRUE   | <GO:0001732> | formation of cytoplasmic translation initiation complex           | 0.0129100 |   0 |   0 |    1 |
| Translation initiation                          | TRUE   | <GO:0001731> | formation of translation preinitiation complex                    | 0.0410500 |   0 |   0 |    1 |
| Translation initiation                          | TRUE   | <GO:0016281> | eukaryotic translation initiation factor 4F complex               | 0.0434600 |   0 |   0 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0051082> | unfolded protein binding                                          | 0.0000044 |   0 |   1 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0005832> | chaperonin-containing T-complex                                   | 0.0001200 |   0 |   1 |    0 |
| Protein folding / chaperone                     | TRUE   | <GO:0006457> | protein folding                                                   | 0.0001400 |   0 |   1 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0005524> | ATP binding                                                       | 0.0004700 |   0 |   1 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0016887> | ATPase activity                                                   | 0.0005500 |   0 |   1 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0050821> | protein stabilization                                             | 0.0029000 |   0 |   1 |    0 |
| Protein folding / chaperone                     | TRUE   | <GO:0051959> | dynein light intermediate chain binding                           | 0.0034600 |   1 |   1 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0008569> | ATP-dependent microtubule motor activity, minus-end-directed      | 0.0068000 |   1 |   1 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0030433> | ubiquitin-dependent ERAD pathway                                  | 0.0075000 |   0 |   0 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0030970> | retrograde protein transport, ER to cytosol                       | 0.0106000 |   0 |   1 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0008017> | microtubule binding                                               | 0.0122700 |   0 |   1 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0030286> | dynein complex                                                    | 0.0191900 |   0 |   0 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0007018> | microtubule-based movement                                        | 0.0283900 |   0 |   0 |    1 |
| Protein folding / chaperone                     | TRUE   | <GO:0005788> | endoplasmic reticulum lumen                                       | 0.0313000 |   0 |   1 |    0 |
| Protein folding / chaperone                     | TRUE   | <GO:0045505> | dynein intermediate chain binding                                 | 0.0420000 |   1 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0051015> | actin filament binding                                            | 0.0001400 |   0 |   1 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0042995> | cell projection                                                   | 0.0007000 |   0 |   1 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0030864> | cortical actin cytoskeleton                                       | 0.0017600 |   0 |   1 |    0 |
| Other: actin filament binding                   | FALSE  | <GO:0048870> | cell motility                                                     | 0.0032300 |   0 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0030036> | actin cytoskeleton organization                                   | 0.0060000 |   0 |   1 |    0 |
| Other: actin filament binding                   | FALSE  | <GO:0005640> | nuclear outer membrane                                            | 0.0070900 |   0 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0031267> | small GTPase binding                                              | 0.0108800 |   0 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0030042> | actin filament depolymerization                                   | 0.0129000 |   0 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0034993> | meiotic nuclear membrane microtubule tethering complex            | 0.0137100 |   0 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0051014> | actin filament severing                                           | 0.0146000 |   0 |   1 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0005635> | nuclear envelope                                                  | 0.0152900 |   0 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0008045> | motor neuron axon guidance                                        | 0.0208400 |   0 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0005885> | Arp2/3 protein complex                                            | 0.0255200 |   0 |   1 |    0 |
| Other: actin filament binding                   | FALSE  | <GO:0030479> | actin cortical patch                                              | 0.0255200 |   0 |   1 |    0 |
| Other: actin filament binding                   | FALSE  | <GO:0032956> | regulation of actin cytoskeleton organization                     | 0.0302100 |   0 |   0 |    1 |
| Other: actin filament binding                   | FALSE  | <GO:0051016> | barbed-end actin filament capping                                 | 0.0419000 |   0 |   1 |    0 |
| Acid-base / V-type ATPase                       | TRUE   | <GO:0046961> | proton-transporting ATPase activity, rotational mechanism         | 0.0002000 |   1 |   1 |    1 |
| Acid-base / V-type ATPase                       | TRUE   | <GO:1902600> | proton transmembrane transport                                    | 0.0057000 |   1 |   1 |    1 |
| Acid-base / V-type ATPase                       | TRUE   | <GO:0016471> | vacuolar proton-transporting V-type ATPase complex                | 0.0076400 |   0 |   1 |    0 |
| Acid-base / V-type ATPase                       | TRUE   | <GO:0033179> | proton-transporting V-type ATPase, V0 domain                      | 0.0076400 |   0 |   1 |    1 |
| Mitochondrial calcium                           | TRUE   | <GO:0036444> | calcium import into the mitochondrion                             | 0.0004200 |   0 |   0 |    1 |
| Mitochondrial calcium                           | TRUE   | <GO:1990246> | uniplex complex                                                   | 0.0004600 |   0 |   0 |    1 |
| Mitochondrial calcium                           | TRUE   | <GO:0051560> | mitochondrial calcium ion homeostasis                             | 0.0010200 |   0 |   0 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0001669> | acrosomal vesicle                                                 | 0.0005100 |   1 |   0 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0036064> | ciliary basal body                                                | 0.0006300 |   0 |   1 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0060294> | cilium movement involved in cell motility                         | 0.0007000 |   0 |   0 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0031514> | motile cilium                                                     | 0.0014300 |   1 |   1 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0003341> | cilium movement                                                   | 0.0047900 |   1 |   0 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0005929> | cilium                                                            | 0.0048700 |   1 |   0 |    0 |
| Cilium / ciliary motility                       | TRUE   | <GO:0015630> | microtubule cytoskeleton                                          | 0.0065800 |   0 |   1 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0044458> | motile cilium assembly                                            | 0.0106000 |   0 |   1 |    0 |
| Cilium / ciliary motility                       | TRUE   | <GO:0060285> | cilium-dependent cell motility                                    | 0.0127600 |   0 |   0 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0045504> | dynein heavy chain binding                                        | 0.0210000 |   0 |   0 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0060271> | cilium assembly                                                   | 0.0221000 |   0 |   1 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0070286> | axonemal dynein complex assembly                                  | 0.0247000 |   0 |   1 |    0 |
| Cilium / ciliary motility                       | TRUE   | <GO:0005814> | centriole                                                         | 0.0252200 |   1 |   0 |    0 |
| Cilium / ciliary motility                       | TRUE   | <GO:0035082> | axoneme assembly                                                  | 0.0265000 |   0 |   0 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0036126> | sperm flagellum                                                   | 0.0320000 |   0 |   0 |    1 |
| Cilium / ciliary motility                       | TRUE   | <GO:0005930> | axoneme                                                           | 0.0366700 |   0 |   0 |    1 |
| Stress granule / mRNA processing                | TRUE   | <GO:0006376> | mRNA splice site selection                                        | 0.0006800 |   1 |   1 |    0 |
| Stress granule / mRNA processing                | TRUE   | <GO:0010494> | cytoplasmic stress granule                                        | 0.0007400 |   1 |   0 |    0 |
| Stress granule / mRNA processing                | TRUE   | <GO:0008143> | poly(A) binding                                                   | 0.0012100 |   0 |   0 |    1 |
| Stress granule / mRNA processing                | TRUE   | <GO:0003730> | mRNA 3’-UTR binding                                               | 0.0014000 |   1 |   0 |    0 |
| Stress granule / mRNA processing                | TRUE   | <GO:0008266> | poly(U) RNA binding                                               | 0.0048000 |   1 |   0 |    1 |
| Stress granule / mRNA processing                | TRUE   | <GO:0000381> | regulation of alternative mRNA splicing, via spliceosome          | 0.0201900 |   1 |   0 |    0 |
| Stress granule / mRNA processing                | TRUE   | <GO:0071004> | U2-type prespliceosome                                            | 0.0237500 |   0 |   0 |    1 |
| Stress granule / mRNA processing                | TRUE   | <GO:0045727> | positive regulation of translation                                | 0.0432200 |   1 |   0 |    0 |
| Other: cytoskeleton                             | FALSE  | <GO:0005856> | cytoskeleton                                                      | 0.0014100 |   0 |   1 |    1 |
| Other: cytoskeleton                             | FALSE  | <GO:0045732> | positive regulation of protein catabolic process                  | 0.0038000 |   0 |   1 |    0 |
| Other: cytoskeleton                             | FALSE  | <GO:1904491> | protein localization to ciliary transition zone                   | 0.0044000 |   0 |   1 |    0 |
| Other: cytoskeleton                             | FALSE  | <GO:0016459> | myosin complex                                                    | 0.0089100 |   0 |   0 |    1 |
| Other: cytoskeleton                             | FALSE  | <GO:0003774> | motor activity                                                    | 0.0110800 |   0 |   0 |    1 |
| Other: cytoskeleton                             | FALSE  | <GO:0051010> | microtubule plus-end binding                                      | 0.0152000 |   0 |   1 |    0 |
| Other: cytoskeleton                             | FALSE  | <GO:0090090> | negative regulation of canonical Wnt signaling pathway            | 0.0174500 |   0 |   1 |    1 |
| Other: cytoskeleton                             | FALSE  | <GO:0005815> | microtubule organizing center                                     | 0.0360100 |   0 |   1 |    0 |
| Other: cytoskeleton                             | FALSE  | <GO:0005881> | cytoplasmic microtubule                                           | 0.0374100 |   0 |   1 |    0 |
| Other: cytoskeleton                             | FALSE  | <GO:0035371> | microtubule plus-end                                              | 0.0375500 |   0 |   1 |    0 |
| Na/K ionoregulation                             | TRUE   | <GO:0006883> | cellular sodium ion homeostasis                                   | 0.0018600 |   1 |   1 |    1 |
| Na/K ionoregulation                             | TRUE   | <GO:0030007> | cellular potassium ion homeostasis                                | 0.0018600 |   1 |   1 |    1 |
| Na/K ionoregulation                             | TRUE   | <GO:0036376> | sodium ion export across plasma membrane                          | 0.0018600 |   1 |   1 |    1 |
| Na/K ionoregulation                             | TRUE   | <GO:1990573> | potassium ion import across plasma membrane                       | 0.0201900 |   1 |   0 |    1 |
| Other: cell-cell adhesion mediator activity     | FALSE  | <GO:0098632> | cell-cell adhesion mediator activity                              | 0.0057300 |   0 |   0 |    1 |
| Other: cell-cell adhesion mediator activity     | FALSE  | <GO:0007411> | axon guidance                                                     | 0.0150100 |   0 |   0 |    1 |
| Other: cell-cell adhesion mediator activity     | FALSE  | <GO:0070593> | dendrite self-avoidance                                           | 0.0382400 |   0 |   0 |    1 |
| Other: chromatin binding                        | FALSE  | <GO:0003682> | chromatin binding                                                 | 0.0082300 |   0 |   1 |    1 |
| Other: chromatin binding                        | FALSE  | <GO:0042393> | histone binding                                                   | 0.0109300 |   0 |   1 |    0 |
| Other: chromatin binding                        | FALSE  | <GO:0000122> | negative regulation of transcription by RNA polymerase II         | 0.0246000 |   0 |   1 |    0 |
| Ubiquitin ligase / proteolysis                  | TRUE   | <GO:0031625> | ubiquitin protein ligase binding                                  | 0.0128200 |   0 |   0 |    1 |
| Ubiquitin ligase / proteolysis                  | TRUE   | <GO:0019941> | modification-dependent protein catabolic process                  | 0.0378400 |   0 |   0 |    1 |
| Ubiquitin ligase / proteolysis                  | TRUE   | <GO:0030162> | regulation of proteolysis                                         | 0.0413300 |   0 |   0 |    1 |
| Ubiquitin ligase / proteolysis                  | TRUE   | <GO:0031386> | protein tag                                                       | 0.0489400 |   0 |   0 |    1 |
| Other: intracellular signal transduction        | FALSE  | <GO:0035556> | intracellular signal transduction                                 | 0.0134000 |   0 |   1 |    0 |
| Other: intracellular signal transduction        | FALSE  | <GO:0004707> | MAP kinase activity                                               | 0.0284400 |   0 |   1 |    0 |
| Other: intracellular signal transduction        | FALSE  | <GO:0001653> | peptide receptor activity                                         | 0.0406400 |   0 |   1 |    0 |
| Other: syntaxin binding                         | FALSE  | <GO:0019905> | syntaxin binding                                                  | 0.0148500 |   0 |   1 |    1 |
| Other: syntaxin binding                         | FALSE  | <GO:0050708> | regulation of protein secretion                                   | 0.0191000 |   0 |   1 |    0 |
| Other: syntaxin binding                         | FALSE  | <GO:0016082> | synaptic vesicle priming                                          | 0.0208400 |   0 |   0 |    1 |
| Other: syntaxin binding                         | FALSE  | <GO:0031629> | synaptic vesicle fusion to presynaptic active zone membrane       | 0.0208400 |   0 |   0 |    1 |
| Steroid / progesterone metabolism               | TRUE   | <GO:0004508> | steroid 17-alpha-monooxygenase activity                           | 0.0221100 |   0 |   0 |    1 |
| Steroid / progesterone metabolism               | TRUE   | <GO:0047442> | 17-alpha-hydroxyprogesterone aldolase activity                    | 0.0221100 |   0 |   0 |    1 |
| Steroid / progesterone metabolism               | TRUE   | <GO:0042448> | progesterone metabolic process                                    | 0.0302700 |   0 |   0 |    1 |
| Other: protein tyrosine kinase activity         | FALSE  | <GO:0004713> | protein tyrosine kinase activity                                  | 0.0250000 |   0 |   0 |    1 |
| Other: protein tyrosine kinase activity         | FALSE  | <GO:0008305> | integrin complex                                                  | 0.0254800 |   1 |   0 |    0 |
| Other: protein tyrosine kinase activity         | FALSE  | <GO:0033627> | cell adhesion mediated by integrin                                | 0.0280100 |   1 |   0 |    0 |
| Other: protein tyrosine kinase activity         | FALSE  | <GO:0005178> | integrin binding                                                  | 0.0392000 |   1 |   0 |    0 |
| Other: protein tyrosine kinase activity         | FALSE  | <GO:0007229> | integrin-mediated signaling pathway                               | 0.0399700 |   1 |   0 |    0 |
| Other: protein tyrosine kinase activity         | FALSE  | <GO:0005925> | focal adhesion                                                    | 0.0416800 |   0 |   0 |    1 |

Table 1. Data-driven GO clusters (all terms; stress flag and
per-stressor enrichment). Stress-relevant clusters are shown in Fig 4D.

``` r
sessionInfo()
```

    ## R version 3.6.0 (2019-04-26)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Scientific Linux 7.4 (Nitrogen)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/local/lib64/R/lib/libRblas.so
    ## LAPACK: /usr/local/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ##  [1] stats4    parallel  grid      stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] igraph_1.2.8         glmmTMB_1.1.2.3      kableExtra_1.4.0    
    ##  [4] topGO_2.38.1         SparseM_1.84-2       GO.db_3.10.0        
    ##  [7] AnnotationDbi_1.48.0 IRanges_2.20.2       S4Vectors_0.24.4    
    ## [10] Biobase_2.46.0       graph_1.64.0         BiocGenerics_0.32.0 
    ## [13] VennDiagram_1.7.3    futile.logger_1.4.3  ggman_0.99.0        
    ## [16] ggrepel_0.9.6        ggrain_0.0.3         ggridges_0.5.7      
    ## [19] patchwork_1.3.0      scales_1.3.0         ACER_1.0            
    ## [22] poolSeq_0.3.5        Rcpp_1.1.0           matrixStats_1.5.0   
    ## [25] stringi_1.8.7        foreach_1.5.2        data.table_1.17.8   
    ## [28] pcadapt_4.4.1        stringr_1.4.0        qvalue_2.18.0       
    ## [31] dplyr_1.1.4          plyr_1.8.9           tidyr_1.3.1         
    ## [34] ggplot2_3.5.2       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-168         gghalves_0.1.4       bit64_4.6.0-1       
    ##  [4] numDeriv_2016.8-1.1  tools_3.6.0          TMB_1.7.22          
    ##  [7] R6_2.6.1             DBI_1.2.3            colorspace_2.1-2    
    ## [10] withr_3.0.2          tidyselect_1.2.1     bit_4.6.0           
    ## [13] compiler_3.6.0       textshaping_0.3.6    cli_3.6.5           
    ## [16] formatR_1.14         xml2_1.5.1           isoband_0.3.0       
    ## [19] labeling_0.4.3       systemfonts_1.0.4    digest_0.6.39       
    ## [22] minqa_1.2.8          rmarkdown_2.30       svglite_2.0.0       
    ## [25] pkgconfig_2.0.3      htmltools_0.5.9      lme4_1.1-25         
    ## [28] fastmap_1.2.0        rlang_1.1.6          rstudioapi_0.17.1   
    ## [31] RSQLite_2.4.5        farver_2.1.2         generics_0.1.4      
    ## [34] gtools_3.9.5         magrittr_2.0.4       Matrix_1.2-17       
    ## [37] munsell_0.5.1        lifecycle_1.0.4      yaml_2.3.12         
    ## [40] MASS_7.3-51.4        blob_1.2.4           lattice_0.20-38     
    ## [43] splines_3.6.0        knitr_1.50           pillar_1.11.1       
    ## [46] boot_1.3-32          reshape2_1.4.5       codetools_0.2-20    
    ## [49] futile.options_1.0.1 glue_1.8.0           evaluate_1.0.5      
    ## [52] lambda.r_1.2.4       vctrs_0.6.5          png_0.1-8           
    ## [55] nloptr_2.2.1         gtable_0.3.5         purrr_1.0.2         
    ## [58] ggpp_0.4.4           cachem_1.1.0         xfun_0.54           
    ## [61] ragg_1.2.5           viridisLite_0.4.2    tibble_3.3.0        
    ## [64] iterators_1.0.14     memoise_2.0.1        statmod_1.5.1
